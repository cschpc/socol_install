/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

   Yhourstat     yhourmin        Multi-year hourly minimum
   Yhourstat     yhourmax        Multi-year hourly maximum
   Yhourstat     yhourrange      Multi-year hourly range
   Yhourstat     yhoursum        Multi-year hourly sum
   Yhourstat     yhourmean       Multi-year hourly mean
   Yhourstat     yhouravg        Multi-year hourly average
   Yhourstat     yhourstd        Multi-year hourly standard deviation
   Yhourstat     yhourstd1       Multi-year hourly standard deviation (n-1)
   Yhourstat     yhourvar        Multi-year hourly variance
   Yhourstat     yhourvar1       Multi-year hourly variance (n-1)

   Dhourstat     dhourmin        Multi-day hourly minimum
   Dhourstat     dhourmax        Multi-day hourly maximum
   Dhourstat     dhourrange      Multi-day hourly range
   Dhourstat     dhoursum        Multi-day hourly sum
   Dhourstat     dhourmean       Multi-day hourly mean
   Dhourstat     dhouravg        Multi-day hourly average
   Dhourstat     dhourstd        Multi-day hourly standard deviation
   Dhourstat     dhourstd1       Multi-day hourly standard deviation (n-1)
   Dhourstat     dhourvar        Multi-day hourly variance
   Dhourstat     dhourvar1       Multi-day hourly variance (n-1)
*/

#include <cdi.h>

#include "cdo_options.h"
#include "datetime.h"
#include "process_int.h"
#include "printinfo.h"
#include "field_functions.h"

static void
addOperators(void)
{
  // clang-format off
  cdo_operator_add("yhourrange", FieldFunc_Range, 0, nullptr);
  cdo_operator_add("yhourmin",   FieldFunc_Min,   0, nullptr);
  cdo_operator_add("yhourmax",   FieldFunc_Max,   0, nullptr);
  cdo_operator_add("yhoursum",   FieldFunc_Sum,   0, nullptr);
  cdo_operator_add("yhourmean",  FieldFunc_Mean,  0, nullptr);
  cdo_operator_add("yhouravg",   FieldFunc_Avg,   0, nullptr);
  cdo_operator_add("yhourvar",   FieldFunc_Var,   0, nullptr);
  cdo_operator_add("yhourvar1",  FieldFunc_Var1,  0, nullptr);
  cdo_operator_add("yhourstd",   FieldFunc_Std,   0, nullptr);
  cdo_operator_add("yhourstd1",  FieldFunc_Std1,  0, nullptr);

  cdo_operator_add("dhourrange", FieldFunc_Range, 0, nullptr);
  cdo_operator_add("dhourmin",   FieldFunc_Min,   1, nullptr);
  cdo_operator_add("dhourmax",   FieldFunc_Max,   1, nullptr);
  cdo_operator_add("dhoursum",   FieldFunc_Sum,   1, nullptr);
  cdo_operator_add("dhourmean",  FieldFunc_Mean,  1, nullptr);
  cdo_operator_add("dhouravg",   FieldFunc_Avg,   1, nullptr);
  cdo_operator_add("dhourvar",   FieldFunc_Var,   1, nullptr);
  cdo_operator_add("dhourvar1",  FieldFunc_Var1,  1, nullptr);
  cdo_operator_add("dhourstd",   FieldFunc_Std,   1, nullptr);
  cdo_operator_add("dhourstd1",  FieldFunc_Std1,  1, nullptr);
  // clang-format on
}

class ModuleYhourStat
{
private:
  bool ldaily;
  bool lrange;
  bool lmean;
  bool lstd;
  bool lvarstd;
  bool lvars2;
  int divisor;
  int MaxHours;
  CdoStreamID streamID1;
  CdoStreamID streamID2;

  int vlistID1;
  int taxisID1;
  int taxisID2;
  std::vector<int> hourot_numSets;  // hour of time
  FieldVector3D vars1;
  FieldVector3D vars2;
  FieldVector3D samp1;

  std::vector<DateTimeList> dtlist;
  VarList varList;

  std::vector<RecordInfo> recList;

  int VARS_MEMTYPE = 0;
  int maxrecs;
  int operfunc;

  Field field;

public:
  void
  init(void *process)
  {
    TimeStat timestat_date = TimeStat::LAST;

    cdo_initialize(process);

    addOperators();

    auto operatorID = cdo_operator_id();
    operfunc = cdo_operator_f1(operatorID);

    ldaily = (cdo_operator_f2(operatorID) == 1);

    lrange = (operfunc == FieldFunc_Range);
    lmean = (operfunc == FieldFunc_Mean || operfunc == FieldFunc_Avg);
    lstd = (operfunc == FieldFunc_Std || operfunc == FieldFunc_Std1);
    lvarstd = (lstd || operfunc == FieldFunc_Var || operfunc == FieldFunc_Var1);
    lvars2 = (lvarstd || lrange);
    divisor = (operfunc == FieldFunc_Std1 || operfunc == FieldFunc_Var1);

    operator_check_argc(0);

    streamID1 = cdo_open_read(0);

    vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = vlistDuplicate(vlistID1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    taxisWithBounds(taxisID2);
    if (taxisInqType(taxisID2) == TAXIS_FORECAST) taxisDefType(taxisID2, TAXIS_RELATIVE);

    MaxHours = ldaily ? 25 : 9301;                   // year: 31*12*25 + 1
    hourot_numSets = std::vector<int>(MaxHours, 0);  // hour of time
    vars1 = FieldVector3D(MaxHours);
    vars2 = FieldVector3D(MaxHours);
    samp1 = FieldVector3D(MaxHours);

    dtlist = std::vector<DateTimeList>(MaxHours);
    for (int hourot = 0; hourot < MaxHours; ++hourot)
      {
        dtlist[hourot].set_stat(timestat_date);
        dtlist[hourot].set_calendar(taxisInqCalendar(taxisID1));
      }

    vlistDefTaxis(vlistID2, taxisID2);

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);

    maxrecs = vlistNrecs(vlistID1);
    recList = std::vector<RecordInfo>(maxrecs);

    varListInit(varList, vlistID1);

    if ((operfunc == FieldFunc_Min) || (operfunc == FieldFunc_Max)) VARS_MEMTYPE = FIELD_NAT;
  }
  void
  run()
  {

    auto field2_stdvar_func = lstd ? field2_std : field2_var;
    auto fieldc_stdvar_func = lstd ? fieldc_std : fieldc_var;
    int tsID = 0;
    int otsID = 0;
    while (true)
      {
        auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
        if (nrecs == 0) break;

        auto vDateTime = taxisInqVdatetime(taxisID1);

        if (Options::cdoVerbose) cdo_print("process timestep: %d %s", tsID + 1, datetime_to_string(vDateTime));

        auto hourot = ldaily ? decode_hour_of_day(vDateTime) : decode_hour_of_year(vDateTime);
        if (hourot < 0 || hourot >= MaxHours)
          cdo_abort("Hour of year %d out of range (%s)!", hourot, datetime_to_string(vDateTime));

        dtlist[hourot].taxis_inq_timestep(taxisID1, hourot_numSets[hourot]);

        if (!vars1[hourot].size())
          {
            fields_from_vlist(vlistID1, samp1[hourot]);
            fields_from_vlist(vlistID1, vars1[hourot], FIELD_VEC | VARS_MEMTYPE);
            if (lvars2) fields_from_vlist(vlistID1, vars2[hourot], FIELD_VEC);
          }

        for (int recID = 0; recID < nrecs; ++recID)
          {
            int varID, levelID;
            cdo_inq_record(streamID1, &varID, &levelID);
            const auto &var = varList[varID];

            if (tsID == 0) recList[recID].set(varID, levelID);

            auto &rsamp1 = samp1[hourot][varID][levelID];
            auto &rvars1 = vars1[hourot][varID][levelID];

            auto numSets = hourot_numSets[hourot];

            if (numSets == 0)
              {
                cdo_read_record(streamID1, rvars1);
                if (lrange)
                  {
                    vars2[hourot][varID][levelID].nmiss = rvars1.nmiss;
                    vars2[hourot][varID][levelID].vec_d = rvars1.vec_d;
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
                if      (lvarstd) field2_sumsumq(rvars1, vars2[hourot][varID][levelID], field);
                else if (lrange)  field2_maxmin(rvars1, vars2[hourot][varID][levelID], field);
                else              field2_function(rvars1, field, operfunc);
                // clang-format on
              }
          }

        if (hourot_numSets[hourot] == 0 && lvarstd)
          for (int recID = 0; recID < maxrecs; ++recID)
            {
              auto [varID, levelID] = recList[recID].get();
              if (varList[varID].isConstant) continue;

              field2_moq(vars2[hourot][varID][levelID], vars1[hourot][varID][levelID]);
            }

        hourot_numSets[hourot]++;
        tsID++;
      }

    for (int hourot = 0; hourot < MaxHours; ++hourot)
      if (hourot_numSets[hourot])
        {
          auto numSets = hourot_numSets[hourot];
          for (int recID = 0; recID < maxrecs; ++recID)
            {
              auto [varID, levelID] = recList[recID].get();
              if (varList[varID].isConstant) continue;

              const auto &rsamp1 = samp1[hourot][varID][levelID];
              auto &rvars1 = vars1[hourot][varID][levelID];

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
                    field2_stdvar_func(rvars1, vars2[hourot][varID][levelID], rsamp1, divisor);
                  else
                    fieldc_stdvar_func(rvars1, vars2[hourot][varID][levelID], numSets, divisor);
                }
              else if (lrange)
                {
                  field2_sub(rvars1, vars2[hourot][varID][levelID]);
                }
            }

          dtlist[hourot].stat_taxis_def_timestep(taxisID2, hourot_numSets[hourot]);
          cdo_def_timestep(streamID2, otsID);

          for (int recID = 0; recID < maxrecs; ++recID)
            {
              auto [varID, levelID] = recList[recID].get();
              if (otsID && varList[varID].isConstant) continue;

              auto &rvars1 = vars1[hourot][varID][levelID];

              cdo_def_record(streamID2, varID, levelID);
              cdo_write_record(streamID2, rvars1);
            }

          otsID++;
        }
  }
  void
  close()
  {
    cdo_stream_close(streamID2);
    cdo_stream_close(streamID1);

    cdo_finish();
  }
};

void *
Yhourstat(void *process)
{
  ModuleYhourStat yhourstat;
  yhourstat.init(process);
  yhourstat.run();
  yhourstat.close();

  return nullptr;
}
