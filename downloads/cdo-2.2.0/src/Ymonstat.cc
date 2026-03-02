/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Ymonstat   ymonrange       Multi-year monthly range
      Ymonstat   ymonmin         Multi-year monthly minimum
      Ymonstat   ymonmax         Multi-year monthly maximum
      Ymonstat   ymonsum         Multi-year monthly sum
      Ymonstat   ymonmean        Multi-year monthly mean
      Ymonstat   ymonavg         Multi-year monthly average
      Ymonstat   ymonvar         Multi-year monthly variance
      Ymonstat   ymonvar1        Multi-year monthly variance [Normalize by (n-1)]
      Ymonstat   ymonstd         Multi-year monthly standard deviation
      Ymonstat   ymonstd1        Multi-year monthly standard deviation [Normalize by (n-1)]
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
  cdo_operator_add("ymonrange", FieldFunc_Range, 0, nullptr);
  cdo_operator_add("ymonmin",   FieldFunc_Min,   0, nullptr);
  cdo_operator_add("ymonmax",   FieldFunc_Max,   0, nullptr);
  cdo_operator_add("ymonmin",   FieldFunc_Min,   0, nullptr);
  cdo_operator_add("ymonmax",   FieldFunc_Max,   0, nullptr);
  cdo_operator_add("ymonsum",   FieldFunc_Sum,   0, nullptr);
  cdo_operator_add("ymonmean",  FieldFunc_Mean,  0, nullptr);
  cdo_operator_add("ymonavg",   FieldFunc_Avg,   0, nullptr);
  cdo_operator_add("ymonvar",   FieldFunc_Var,   0, nullptr);
  cdo_operator_add("ymonvar1",  FieldFunc_Var1,  0, nullptr);
  cdo_operator_add("ymonstd",   FieldFunc_Std,   0, nullptr);
  cdo_operator_add("ymonstd1",  FieldFunc_Std1,  0, nullptr);
  // clang-format on
}

class ModuleYmonstat
{
  const TimeStat timestat_date = TimeStat::LAST;
  static const int MaxMonths = 17;
  int month_numSets[MaxMonths] = { 0 };
  int mon[MaxMonths] = { 0 };
  int nmon = 0;
  FieldVector2D vars1[MaxMonths], vars2[MaxMonths], samp1[MaxMonths];
  int VARS_MEMTYPE = 0;
  DateTimeList dtlists[MaxMonths];

  CdoStreamID streamID1;
  CdoStreamID streamID2;
  int vlistID1;
  int taxisID1;
  int taxisID2;

  std::vector<RecordInfo> recList;
  VarList varList;

  int operfunc;
  int maxrecs;
  bool lrange;
  bool lmean;
  bool lstd;
  bool lvarstd;
  bool lvars2;
  int divisor;

  Field field;

public:
  void
  init(void *process)
  {

    cdo_initialize(process);

    addOperators();

    auto operatorID = cdo_operator_id();
    operfunc = cdo_operator_f1(operatorID);

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
    vlistDefTaxis(vlistID2, taxisID2);

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);

    maxrecs = vlistNrecs(vlistID1);
    recList = std::vector<RecordInfo>(maxrecs);

    for (int month = 0; month < MaxMonths; ++month)
      {
        dtlists[month].set_stat(timestat_date);
        dtlists[month].set_calendar(taxisInqCalendar(taxisID1));
      }

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

        int month = vDateTime.date.month;
        if (month < 0 || month >= MaxMonths) cdo_abort("Month %d out of range!", month);

        dtlists[month].taxis_set_next_timestep(taxisID1);

        if (!vars1[month].size())
          {
            mon[nmon++] = month;
            fields_from_vlist(vlistID1, samp1[month]);
            fields_from_vlist(vlistID1, vars1[month], FIELD_VEC | VARS_MEMTYPE);
            if (lvars2) fields_from_vlist(vlistID1, vars2[month], FIELD_VEC);
          }

        for (int recID = 0; recID < nrecs; ++recID)
          {
            int varID, levelID;
            cdo_inq_record(streamID1, &varID, &levelID);
            const auto &var = varList[varID];

            if (tsID == 0) recList[recID].set(varID, levelID);

            auto &rsamp1 = samp1[month][varID][levelID];
            auto &rvars1 = vars1[month][varID][levelID];

            auto numSets = month_numSets[month];

            if (numSets == 0)
              {
                cdo_read_record(streamID1, rvars1);
                if (lrange)
                  {
                    vars2[month][varID][levelID].nmiss = rvars1.nmiss;
                    vars2[month][varID][levelID].vec_d = rvars1.vec_d;
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
                if      (lvarstd) field2_sumsumq(rvars1, vars2[month][varID][levelID], field);
                else if (lrange)  field2_maxmin(rvars1, vars2[month][varID][levelID], field);
                else              field2_function(rvars1, field, operfunc);
                // clang-format on
              }
          }

        if (month_numSets[month] == 0 && lvarstd)
          for (int recID = 0; recID < maxrecs; ++recID)
            {
              auto [varID, levelID] = recList[recID].get();
              if (varList[varID].isConstant) continue;

              field2_moq(vars2[month][varID][levelID], vars1[month][varID][levelID]);
            }

        month_numSets[month]++;
        tsID++;
      }

    if (nmon == 12)
      {
        int smon = 0;
        for (int month = 1; month <= 12; ++month)
          if (month_numSets[month]) smon++;
        if (smon == 12)
          for (int month = 1; month <= 12; ++month) mon[month - 1] = month;
      }

    for (int i = 0; i < nmon; ++i)
      {
        auto month = mon[i];
        auto numSets = month_numSets[month];
        if (numSets == 0) cdo_abort("Internal problem, numSets[%d] not defined!", month);

        for (int recID = 0; recID < maxrecs; ++recID)
          {
            auto [varID, levelID] = recList[recID].get();
            if (varList[varID].isConstant) continue;

            const auto &rsamp1 = samp1[month][varID][levelID];
            auto &rvars1 = vars1[month][varID][levelID];

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
                  field2_stdvar_func(rvars1, vars2[month][varID][levelID], rsamp1, divisor);
                else
                  fieldc_stdvar_func(rvars1, vars2[month][varID][levelID], numSets, divisor);
              }
            else if (lrange) { field2_sub(rvars1, vars2[month][varID][levelID]); }
          }

        dtlists[month].stat_taxis_def_timestep(taxisID2);
        cdo_def_timestep(streamID2, otsID);

        for (int recID = 0; recID < maxrecs; ++recID)
          {
            auto [varID, levelID] = recList[recID].get();
            if (otsID && varList[varID].isConstant) continue;

            auto &rvars1 = vars1[month][varID][levelID];

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
Ymonstat(void *process)
{
  ModuleYmonstat ymonstat;
  ymonstat.init(process);
  ymonstat.run();
  ymonstat.close();
  return nullptr;
}
