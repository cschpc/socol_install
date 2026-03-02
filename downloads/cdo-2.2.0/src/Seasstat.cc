/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Seasstat   seasrange       Seasonal range
      Seasstat   seasmin         Seasonal minimum
      Seasstat   seasmax         Seasonal maximum
      Seasstat   seassum         Seasonal sum
      Seasstat   seasmean        Seasonal mean
      Seasstat   seasavg         Seasonal average
      Seasstat   seasvar         Seasonal variance
      Seasstat   seasvar1        Seasonal variance [Normalize by (n-1)]
      Seasstat   seasstd         Seasonal standard deviation
      Seasstat   seasstd1        Seasonal standard deviation [Normalize by (n-1)]
*/

#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "datetime.h"
#include "printinfo.h"
#include "cdo_season.h"
#include "field_functions.h"

static void
addOperators(void)
{
  // clang-format off
  cdo_operator_add("seasrange", FieldFunc_Range, 0, nullptr);
  cdo_operator_add("seasmin",   FieldFunc_Min,   0, nullptr);
  cdo_operator_add("seasmax",   FieldFunc_Max,   0, nullptr);
  cdo_operator_add("seassum",   FieldFunc_Sum,   0, nullptr);
  cdo_operator_add("seasmean",  FieldFunc_Mean,  0, nullptr);
  cdo_operator_add("seasavg",   FieldFunc_Avg,   0, nullptr);
  cdo_operator_add("seasvar",   FieldFunc_Var,   0, nullptr);
  cdo_operator_add("seasvar1",  FieldFunc_Var1,  0, nullptr);
  cdo_operator_add("seasstd",   FieldFunc_Std,   0, nullptr);
  cdo_operator_add("seasstd1",  FieldFunc_Std1,  0, nullptr);
  // clang-format on
}

void *
Seasstat(void *process)
{
  TimeStat timestat_date = TimeStat::MEAN;
  CdiDateTime vDateTime0{};
  CdiDateTime vDateTime1{};
  int seas0 = 0;
  int oldmon = 0;
  int nseason = 0;

  cdo_initialize(process);

  addOperators();

  auto operatorID = cdo_operator_id();
  auto operfunc = cdo_operator_f1(operatorID);

  operator_check_argc(0);

  auto lrange = (operfunc == FieldFunc_Range);
  auto lmean = (operfunc == FieldFunc_Mean || operfunc == FieldFunc_Avg);
  auto lstd = (operfunc == FieldFunc_Std || operfunc == FieldFunc_Std1);
  auto lvarstd = (lstd || operfunc == FieldFunc_Var || operfunc == FieldFunc_Var1);
  auto lvars2 = (lvarstd || lrange);
  const int divisor = (operfunc == FieldFunc_Std1 || operfunc == FieldFunc_Var1);

  auto field2_stdvar_func = lstd ? field2_std : field2_var;
  auto fieldc_stdvar_func = lstd ? fieldc_std : fieldc_var;

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

  auto maxrecs = vlistNrecs(vlistID1);
  std::vector<RecordInfo> recList(maxrecs);

  DateTimeList dtlist;
  dtlist.set_stat(timestat_date);
  dtlist.set_calendar(taxisInqCalendar(taxisID1));

  VarList varList;
  varListInit(varList, vlistID1);

  int VARS_MEMTYPE = 0;
  if ((operfunc == FieldFunc_Min) || (operfunc == FieldFunc_Max)) VARS_MEMTYPE = FIELD_NAT;

  Field field;

  FieldVector2D samp1, vars1, vars2;
  fields_from_vlist(vlistID1, samp1);
  fields_from_vlist(vlistID1, vars1, FIELD_VEC | VARS_MEMTYPE);
  if (lvars2) fields_from_vlist(vlistID1, vars2, FIELD_VEC);

  int tsID = 0;
  int otsID = 0;
  while (true)
    {
      int nrecs = 0;
      long numSets = 0;
      bool newseas = false;
      while (true)
        {
          nrecs = cdo_stream_inq_timestep(streamID1, tsID);
          if (nrecs == 0) break;

          dtlist.taxis_inq_timestep(taxisID1, numSets);
          auto vDateTime = dtlist.get_vDateTime(numSets);

          auto month = decode_month(vDateTime.date);
          auto newmon = month;
          if (seasonStart == SeasonStart::DEC && newmon == 12) newmon = 0;

          auto seas = month_to_season(month);

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

              if (numSets == 0)
                {
                  cdo_read_record(streamID1, rvars1);
                  if (lrange)
                    {
                      vars2[varID][levelID].nmiss = rvars1.nmiss;
                      vars2[varID][levelID].vec_d = rvars1.vec_d;
                    }

                  if (rvars1.nmiss || !rsamp1.empty())
                    {
                      if (rsamp1.empty()) rsamp1.resize(rvars1.size);
                      field2_vinit(rsamp1, rvars1);
                    }
                }
              else
                {
                  field.init(var);
                  cdo_read_record(streamID1, field);

                  if (field.nmiss || !rsamp1.empty())
                    {
                      if (rsamp1.empty()) rsamp1.resize(rvars1.size, numSets);
                      field2_vincr(rsamp1, field);
                    }

                  // clang-format off
                  if      (lvarstd) field2_sumsumq(rvars1, vars2[varID][levelID], field);
                  else if (lrange)  field2_maxmin(rvars1, vars2[varID][levelID], field);
                  else              field2_function(rvars1, field, operfunc);
                  // clang-format on
                }
            }

          if (numSets == 0 && lvarstd)
            for (int recID = 0; recID < maxrecs; ++recID)
              {
                auto [varID, levelID] = recList[recID].get();
                if (varList[varID].isConstant) continue;

                field2_moq(vars2[varID][levelID], vars1[varID][levelID]);
              }

          vDateTime1 = vDateTime;
          numSets++;
          tsID++;
        }

      if (nrecs == 0 && numSets == 0) break;

      for (int recID = 0; recID < maxrecs; ++recID)
        {
          auto [varID, levelID] = recList[recID].get();
          if (varList[varID].isConstant) continue;

          const auto &rsamp1 = samp1[varID][levelID];
          auto &rvars1 = vars1[varID][levelID];

          if (lmean)
            {
              if (!rsamp1.empty())
                field2_div(rvars1, rsamp1);
              else
                fieldc_div(rvars1, (double) numSets);
            }
          else if (lvarstd)
            {
              if (!rsamp1.empty())
                field2_stdvar_func(rvars1, vars2[varID][levelID], rsamp1, divisor);
              else
                fieldc_stdvar_func(rvars1, vars2[varID][levelID], numSets, divisor);
            }
          else if (lrange) { field2_sub(rvars1, vars2[varID][levelID]); }
        }

      if (Options::cdoVerbose)
        cdo_print("season: %3d %3s  start: %s  end: %s ntimesteps: %ld", nseason, seasonNames[seas0],
                  datetime_to_string(vDateTime0), datetime_to_string(vDateTime1), numSets);

      dtlist.stat_taxis_def_timestep(taxisID2, numSets);
      cdo_def_timestep(streamID2, otsID);

      if (numSets < 3)
        cdo_warning("Season %3d (%s) has only %d input time step%s!", otsID + 1, date_to_string(vDateTime0.date), numSets,
                    (numSets == 1) ? "" : "s");

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
