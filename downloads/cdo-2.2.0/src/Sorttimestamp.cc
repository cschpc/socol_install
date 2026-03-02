/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

     Sorttimestamp    sorttimestamp         Sort all timesteps
*/

#include <algorithm>  // sort

#include <cdi.h>
#include "julian_date.h"

#include "cdo_options.h"
#include "process_int.h"
#include "cdo_vlist.h"
#include "printinfo.h"
#include "field_functions.h"

bool getenv_skip_same_time();

struct TimeInfo
{
  int index;
  double datetime;
};

static bool
cmpdatetime(const TimeInfo &a, const TimeInfo &b)
{
  return a.datetime < b.datetime;
}

void *
Sorttimestamp(void *process)
{
  int varID, levelID;
  int lasttsID = -1;
  int nalloc = 0;
  int vlistID2 = -1, taxisID2 = -1;
  int nvars = 0;

  cdo_initialize(process);

  const auto skipSameTime = getenv_skip_same_time();
  auto unique = false;

  if (cdo_operator_argc() == 1)
    {
      if (cdo_operator_argv(0) == "unique")
        unique = true;
      else
        cdo_abort("Unexpected parameter %s!", cdo_operator_argv(0));
    }
  if (cdo_operator_argc() > 1) operator_check_argc(1);

  FieldVector3D vars;
  std::vector<CdiDateTime> vDateTimes;

  const auto nfiles = cdo_stream_cnt() - 1;

  int xtsID = 0;
  for (int fileID = 0; fileID < nfiles; ++fileID)
    {
      const auto streamID1 = cdo_open_read(fileID);
      const auto vlistID1 = cdo_stream_inq_vlist(streamID1);
      const auto taxisID1 = vlistInqTaxis(vlistID1);

      VarList varList1;
      varListInit(varList1, vlistID1);

      if (fileID == 0)
        {
          vlistID2 = vlistDuplicate(vlistID1);
          taxisID2 = taxisDuplicate(taxisID1);
          if (taxisHasBounds(taxisID2))
            {
              cdo_warning("Time bounds unsupported by this operator, removed!");
              taxisDeleteBounds(taxisID2);
            }
        }
      else { vlist_compare(vlistID2, vlistID1, CMP_ALL); }

      nvars = vlistNvars(vlistID1);

      int tsID = 0;
      while (true)
        {
          const auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
          if (nrecs == 0) break;

          if (xtsID >= nalloc)
            {
              constexpr int NALLOC_INC = 1024;
              nalloc += NALLOC_INC;
              vDateTimes.resize(nalloc);
              vars.resize(nalloc);
            }

          vDateTimes[xtsID] = taxisInqVdatetime(taxisID1);

          fields_from_vlist(vlistID1, vars[xtsID]);

          for (int recID = 0; recID < nrecs; ++recID)
            {
              cdo_inq_record(streamID1, &varID, &levelID);
              auto &field = vars[xtsID][varID][levelID];
              field.init(varList1[varID]);
              cdo_read_record(streamID1, field);
            }

          tsID++;
          xtsID++;
        }

      cdo_stream_close(streamID1);
    }

  int nts = xtsID;

  std::vector<TimeInfo> timeinfo(nts);
  const auto calendar = taxisInqCalendar(taxisID2);

  for (int tsID = 0; tsID < nts; ++tsID)
    {
      const auto julday = date_to_julday(calendar, cdiDate_get(vDateTimes[tsID].date));
      const auto secofday = time_to_sec(cdiTime_get(vDateTimes[tsID].time));
      const double jdatetime = julday + secofday / 86400.0;
      timeinfo[tsID].index = tsID;
      timeinfo[tsID].datetime = jdatetime;
    }

  std::stable_sort(timeinfo.begin(), timeinfo.end(), cmpdatetime);

  vlistDefTaxis(vlistID2, taxisID2);

  const auto streamID2 = cdo_open_write(nfiles);
  cdo_def_vlist(streamID2, vlistID2);

  int tsID2 = 0;
  for (int tsID = 0; tsID < nts; ++tsID)
    {
      xtsID = timeinfo[tsID].index;

      if (tsID > 0 && IS_EQUAL(timeinfo[tsID].datetime, timeinfo[lasttsID].datetime))
        {
          if (skipSameTime)
            {
              if (Options::cdoVerbose)
                cdo_print("Timestep %4d %s already exists, skipped!", xtsID + 1, datetime_to_string(vDateTimes[xtsID]));
              continue;
            }

          if (unique)
            {
              auto lskip = false;
              auto xtsID2 = timeinfo[lasttsID].index;
              const auto &field1 = vars[xtsID][0][0];
              const auto &field2 = vars[xtsID2][0][0];
              if (field1.memType == MemType::Float)
                {
                  if (field1.vec_f == field2.vec_f) lskip = true;
                }
              else
                {
                  if (field1.vec_d == field2.vec_d) lskip = true;
                }

              if (lskip)
                {
                  if (Options::cdoVerbose)
                    cdo_print("Timestep %4d %s already exists with the same data, skipped!", xtsID + 1,
                              datetime_to_string(vDateTimes[xtsID]));
                  continue;
                }
            }
        }

      lasttsID = tsID;

      taxisDefVdatetime(taxisID2, vDateTimes[xtsID]);
      cdo_def_timestep(streamID2, tsID2++);

      for (varID = 0; varID < nvars; ++varID)
        {
          const auto nlevels = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID));
          for (levelID = 0; levelID < nlevels; ++levelID)
            {
              auto &field = vars[xtsID][varID][levelID];
              if (field.hasData())
                {
                  cdo_def_record(streamID2, varID, levelID);
                  cdo_write_record(streamID2, field);
                }
            }
        }
    }

  cdo_stream_close(streamID2);

  cdo_finish();

  return nullptr;
}
