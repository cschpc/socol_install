/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Fabian Wachsmann

*/

#include <cdi.h>

#include "process_int.h"
#include "cdo_options.h"

void *
EstFreq(void *process)
{

  cdo_initialize(process);

  operator_check_argc(0);

  auto dataIsUnchanged = data_is_unchanged();

  auto streamID1 = cdo_open_read(0);

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto vlistID2 = vlistDuplicate(vlistID1);

  auto taxisID1 = vlistInqTaxis(vlistID1);
  auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  auto streamID2 = cdo_open_write(1);

  cdo_def_vlist(streamID2, vlistID2);

  auto ntsteps = vlistNtsteps(vlistID1);

  Field field;

  VarList varList1;
  varListInit(varList1, vlistID1);

  int tsID = 0;
  int fyear = 0, lyear, fmonth = 0, lmonth, dummy;
  int step_per_year = 0, currentyear, currentmon;

  while (true)
    {
      auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
      if (nrecs == 0) break;

      cdiDate_decode(taxisInqVdatetime(taxisID1).date, &currentyear, &currentmon, &dummy);

      if (tsID == 0)
        {
          fyear = currentyear;
          fmonth = currentmon;
        }

      if (currentyear == fyear)
        {
          // lymonth = currentmon;
          step_per_year++;
        }

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_def_timestep(streamID2, tsID);

      for (int recID = 0; recID < nrecs; ++recID)
        {
          int varID, levelID;
          cdo_inq_record(streamID1, &varID, &levelID);
          cdo_def_record(streamID2, varID, levelID);

          if (dataIsUnchanged) { cdo_copy_record(streamID2, streamID1); }
          else
            {
              field.init(varList1[varID]);
              cdo_read_record(streamID1, field);
              cdo_write_record(streamID2, field);
            }
        }

      tsID++;
    }

  const char *frequency = "";

  if (ntsteps > 2)
    {
      cdo_stream_inq_timestep(streamID2, ntsteps);
      cdiDate_decode(taxisInqVdatetime(taxisID2).date, &lyear, &lmonth, &dummy);
      /* First, estimation by maximal number of time steps divided by covered years between last and first time step */
      if (Options::cdoVerbose)
        printf("Frequency is calculated by dividing the number of time steps '%d' included in the time axis by the covered years"
               " of the time axis\ncomputed by the difference of the year of the last time stamp '%d' and the year of the first"
               " time stamp '%d'.\n",
               ntsteps, lyear, fyear);
      double covered_years = lyear - fyear + 1.0;
      double freq = ntsteps / covered_years;
      if (DBL_IS_EQUAL(freq, 1.))
        frequency = "yr";
      else if (DBL_IS_EQUAL(freq, 12.))
        frequency = "mon";
      else if (IS_EQUAL(freq, 365.) || IS_EQUAL(freq, 365.25) || IS_EQUAL(freq, 366.))
        frequency = "day";
      else if (IS_EQUAL(freq, 365. * 4) || IS_EQUAL(freq, 365.25 * 4) || IS_EQUAL(freq, 366. * 4))
        frequency = "6hr";
      else if (IS_EQUAL(freq, 365. * 8) || IS_EQUAL(freq, 365.25 * 8) || IS_EQUAL(freq, 366. * 8))
        frequency = "3hr";
      else
        {
          int covered_months = lmonth - fmonth + 1;
          if (Options::cdoVerbose)
            printf("The fraction ntsteps / covered_years = '%f' is neither 1, 12, 365, 365.25, 366 nor a multiple of 365 which "
                   "would correspond to frequencies yearly, monthly, daily or subdaily respectively.\n Next try:\n\nFrequency is "
                   "calculated by dividing the number of time steps '%d' in year '%d' by the covered months in that year '%d'.\n",
                   ntsteps / covered_years, step_per_year, fyear, covered_months);
          if (step_per_year > 366 * 8)
            cdo_abort("Step per year '%d' in year '%d' is bigger than 366*8 which corresponds to a frequency of sub-3hourly!"
                      " This is not yet enabled.",
                      step_per_year, fyear);
          else
            {
              freq = (double) step_per_year / (double) covered_months;
              // clang-format off
              if      (freq > 31 * 8) cdo_abort("Frequency is sub-3hourly! Not yet enabled.");
              else if (freq > 31 * 4) frequency = "3hr";
              else if (freq > 31)     frequency = "6hr";
              else if (freq > 1)      frequency = "day";
              else                    frequency = "mon";
              // clang-format on
            }
        }
    }
  else
    cdo_abort("For %d found timesteps no frequency can be computed - at least 3 timesteps are required.", ntsteps);

  if (Options::cdoVerbose) printf("Your file indicates a frequency of '%s'.\n", frequency);
  cdiDefAttTxt(vlistID2, CDI_GLOBAL, "frequency", 3, frequency);

  cdo_stream_close(streamID1);
  cdo_stream_close(streamID2);

  vlistDestroy(vlistID2);

  cdo_finish();

  return nullptr;
}
