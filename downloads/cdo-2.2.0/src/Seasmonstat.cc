/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Seasmonstat   seasmonmean        Seasonal mean from monthly data
      Seasmonstat   seasmonavg         Seasonal average from monthly data
*/

#include <cdi.h>
#include "calendar.h"

#include "cdo_options.h"
#include "process_int.h"
#include "datetime.h"
#include "printinfo.h"
#include "cdo_season.h"
#include "field_functions.h"

void *
Seasmonstat(void *process)
{
  TimeStat timestat_date = TimeStat::MEAN;
  CdiDateTime vDateTime0{};
  CdiDateTime vDateTime1{};
  int nrecs;
  int seas0 = 0;
  int oldmon = 0;
  int nseason = 0;
  int month0 = 0;
  int year, month, day;

  cdo_initialize(process);

  // clang-format off
  cdo_operator_add("seasmonmean",  FieldFunc_Mean, 0, nullptr);
  cdo_operator_add("seasmonavg",   FieldFunc_Avg,  0, nullptr);
  // clang-format on

  auto operatorID = cdo_operator_id();
  auto operfunc = cdo_operator_f1(operatorID);

  operator_check_argc(0);

  auto seasonStart = get_season_start();
  auto seasonNames = get_season_name();

  auto streamID1 = cdo_open_read(0);

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto vlistID2 = vlistDuplicate(vlistID1);

  auto taxisID1 = vlistInqTaxis(vlistID1);
  auto taxisID2 = taxisDuplicate(taxisID1);
  taxisWithBounds(taxisID2);
  if (taxisInqType(taxisID2) == TAXIS_FORECAST) taxisDefType(taxisID2, TAXIS_RELATIVE);
  vlistDefTaxis(vlistID2, taxisID2);

  auto streamID2 = cdo_open_write(1);
  cdo_def_vlist(streamID2, vlistID2);

  auto nvars = vlistNvars(vlistID1);

  auto maxrecs = vlistNrecs(vlistID1);
  std::vector<RecordInfo> recList(maxrecs);

  auto calendar = taxisInqCalendar(taxisID1);
  DateTimeList dtlist;
  dtlist.set_stat(timestat_date);
  dtlist.set_calendar(calendar);

  VarList varList;
  varListInit(varList, vlistID1);

  Field field;

  FieldVector2D samp1, vars1;
  fields_from_vlist(vlistID1, samp1);
  fields_from_vlist(vlistID1, vars1, FIELD_VEC);

  int tsID = 0;
  int otsID = 0;
  while (true)
    {
      long numSets = 0;
      double dsets = 0.0;
      auto newseas = false;
      while (true)
        {
          nrecs = cdo_stream_inq_timestep(streamID1, tsID);
          if (nrecs <= 0) break;

          dtlist.taxis_inq_timestep(taxisID1, numSets);
          auto vDateTime = dtlist.get_vDateTime(numSets);
          cdiDate_decode(vDateTime.date, &year, &month, &day);

          auto newmon = month;
          if (seasonStart == SeasonStart::DEC && newmon == 12) newmon = 0;

          auto seas = month_to_season(month);

          if (numSets > 0 && month == month0)
            {
              cdo_warning("   last timestep: %s", datetime_to_string(vDateTime0));
              cdo_warning("current timestep: %s", datetime_to_string(vDateTime));
              cdo_abort("Month does not change!");
            }

          auto dpm = days_per_month(calendar, year, month);

          if (numSets == 0)
            {
              nseason++;
              vDateTime0 = vDateTime;
              seas0 = seas;
              oldmon = newmon;
            }

          if (newmon < oldmon) newseas = true;

          if ((seas != seas0) || newseas)
            {
              cdo_add_steps(-1);
              break;
            }

          oldmon = newmon;

          for (int recID = 0; recID < nrecs; ++recID)
            {
              int varID, levelID;
              cdo_inq_record(streamID1, &varID, &levelID);
              const auto &var = varList[varID];

              if (tsID == 0) recList[recID].set(varID, levelID);

              auto &rsamp1 = samp1[varID][levelID];
              auto &rvars1 = vars1[varID][levelID];

              auto fieldsize = rvars1.size;

              if (numSets == 0)
                {
                  cdo_read_record(streamID1, rvars1);

                  fieldc_mul(rvars1, dpm);

                  if (rvars1.nmiss || !rsamp1.empty())
                    {
                      if (rsamp1.empty()) rsamp1.resize(fieldsize);
                      field2_vinit(rsamp1, rvars1, dpm);
                    }
                }
              else
                {
                  field.init(var);
                  cdo_read_record(streamID1, field);

                  fieldc_mul(field, dpm);

                  if (field.nmiss || !rsamp1.empty())
                    {
                      if (rsamp1.empty()) rsamp1.resize(fieldsize, dsets);
                      field2_vincr(rsamp1, field, dpm);
                    }

                  field2_function(rvars1, field, operfunc);
                }
            }

          month0 = month;
          vDateTime1 = vDateTime;
          numSets++;
          dsets += dpm;
          tsID++;
        }

      if (nrecs == 0 && numSets == 0) break;

      for (int varID = 0; varID < nvars; ++varID)
        {
          const auto &var = varList[varID];
          if (var.isConstant) continue;
          for (int levelID = 0; levelID < var.nlevels; ++levelID)
            {
              const auto &rsamp1 = samp1[varID][levelID];
              auto &rvars1 = vars1[varID][levelID];
              if (!rsamp1.empty())
                field2_div(rvars1, rsamp1);
              else
                fieldc_div(rvars1, dsets);
            }
        }

      if (Options::cdoVerbose)
        cdo_print("season: %3d %3s  start: %s  end: %s ntimesteps: %ld", nseason, seasonNames[seas0],
                  datetime_to_string(vDateTime0), datetime_to_string(vDateTime1), numSets);

      dtlist.stat_taxis_def_timestep(taxisID2, numSets);
      cdo_def_timestep(streamID2, otsID);

      if (numSets < 3)
        cdo_warning("Season %3d (%s) has only %d input time step%s!", otsID + 1, date_to_string(vDateTime0.date), numSets,
                    numSets == 1 ? "" : "s");

      for (int recID = 0; recID < maxrecs; ++recID)
        {
          auto [varID, levelID] = recList[recID].get();
          if (otsID && varList[varID].isConstant) continue;

          auto &rvars1 = vars1[varID][levelID];
          cdo_def_record(streamID2, varID, levelID);
          cdo_write_record(streamID2, rvars1);
        }

      if (nrecs == 0) break;
      otsID++;
    }

  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  cdo_finish();

  return nullptr;
}
