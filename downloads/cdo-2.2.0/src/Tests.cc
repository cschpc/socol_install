/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include "process_int.h"
#include "param_conversion.h"
#include "statistic.h"

void *
Tests(void *process)
{
  double degreeOfFreedom = 0, p = 0, q = 0, n = 0, d = 0;

  cdo_initialize(process);

  // clang-format off
  auto NORMAL    = cdo_operator_add("normal",    0, 0, nullptr);
  auto STUDENTT  = cdo_operator_add("studentt",  0, 0, "degree of freedom");
  auto CHISQUARE = cdo_operator_add("chisquare", 0, 0, "degree of freedom");
  auto BETA      = cdo_operator_add("beta",      0, 0, "p and q");
  auto FISHER    = cdo_operator_add("fisher",    0, 0, "degree of freedom of nominator and of denominator");
  // clang-format on

  auto operatorID = cdo_operator_id();

  if (operatorID == STUDENTT || operatorID == CHISQUARE)
    {
      operator_input_arg(cdo_operator_enter(operatorID));

      operator_check_argc(1);

      degreeOfFreedom = parameter_to_double(cdo_operator_argv(0));
      if (degreeOfFreedom <= 0) cdo_abort("degree of freedom must be positive!");
    }
  else if (operatorID == BETA)
    {
      operator_input_arg(cdo_operator_enter(operatorID));

      operator_check_argc(2);

      p = parameter_to_double(cdo_operator_argv(0));
      q = parameter_to_double(cdo_operator_argv(1));

      if (p <= 0 || q <= 0) cdo_abort("p and q must be positive!");
    }
  else if (operatorID == FISHER)
    {
      operator_input_arg(cdo_operator_enter(operatorID));

      operator_check_argc(2);

      n = parameter_to_double(cdo_operator_argv(0));
      d = parameter_to_double(cdo_operator_argv(1));
      if (n <= 0 || d <= 0) cdo_abort("both degrees must be positive!");
    }

  auto streamID1 = cdo_open_read(0);

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto vlistID2 = vlistDuplicate(vlistID1);

  auto taxisID1 = vlistInqTaxis(vlistID1);
  auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  auto streamID2 = cdo_open_write(1);
  cdo_def_vlist(streamID2, vlistID2);

  auto gridsizemax = vlistGridsizeMax(vlistID1);
  Varray<double> array1(gridsizemax), array2(gridsizemax);

  VarList varList1;
  varListInit(varList1, vlistID1);

  int tsID = 0;
  while (true)
    {
      auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
      if (nrecs == 0) break;

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_def_timestep(streamID2, tsID);

      for (int recID = 0; recID < nrecs; ++recID)
        {
          int varID, levelID;
          size_t nmiss;
          cdo_inq_record(streamID1, &varID, &levelID);
          cdo_read_record(streamID1, array1.data(), &nmiss);

          auto gridsize = varList1[varID].gridsize;
          auto missval = varList1[varID].missval;

          if (operatorID == NORMAL)
            {
              for (size_t i = 0; i < gridsize; ++i) array2[i] = dbl_is_equal(array1[i], missval) ? missval : cdo::normal(array1[i]);
            }
          else if (operatorID == STUDENTT)
            {
              for (size_t i = 0; i < gridsize; ++i)
                array2[i] = dbl_is_equal(array1[i], missval) ? missval : cdo::student_t(degreeOfFreedom, array1[i]);
            }
          else if (operatorID == CHISQUARE)
            {
              for (size_t i = 0; i < gridsize; ++i)
                array2[i] = dbl_is_equal(array1[i], missval) ? missval : cdo::chi_square(degreeOfFreedom, array1[i]);
            }
          else if (operatorID == BETA)
            {
              for (size_t i = 0; i < gridsize; ++i)
                {
                  if (array1[i] < 0 || array1[i] > 1) cdo_abort("Value out of range (0-1)!");

                  array2[i] = dbl_is_equal(array1[i], missval) ? missval : cdo::beta_distr(p, q, array1[i]);
                }
            }
          else if (operatorID == FISHER)
            {
              for (size_t i = 0; i < gridsize; ++i)
                array2[i] = dbl_is_equal(array1[i], missval) ? missval : cdo::fisher(n, d, array1[i]);
            }
          else { cdo_abort("Internal problem, operator not implemented!"); }

          cdo_def_record(streamID2, varID, levelID);
          cdo_write_record(streamID2, array2.data(), nmiss);
        }

      tsID++;
    }

  cdo_stream_close(streamID1);
  cdo_stream_close(streamID2);

  vlistDestroy(vlistID2);

  cdo_finish();

  return nullptr;
}
