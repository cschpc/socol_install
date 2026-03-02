/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Yseasstat  yseasrange      Multi-year seasonal range
      Yseasstat  yseasmin        Multi-year seasonal minimum
      Yseasstat  yseasmax        Multi-year seasonal maximum
      Yseasstat  yseassum        Multi-year seasonal sum
      Yseasstat  yseasmean       Multi-year seasonal mean
      Yseasstat  yseasavg        Multi-year seasonal average
      Yseasstat  yseasvar        Multi-year seasonal variance
      Yseasstat  yseasvar1       Multi-year seasonal variance [Normalize by (n-1)]
      Yseasstat  yseasstd        Multi-year seasonal standard deviation
      Yseasstat  yseasstd1       Multi-year seasonal standard deviation [Normalize by (n-1)]
*/

#include <cdi.h>

#include "cdo_season.h"
#include "datetime.h"
#include "process_int.h"
#include "field_functions.h"

static void
addOperators(void)
{
  // clang-format off
  cdo_operator_add("yseasrange", FieldFunc_Range, 0, nullptr);
  cdo_operator_add("yseasmin",   FieldFunc_Min,   0, nullptr);
  cdo_operator_add("yseasmax",   FieldFunc_Max,   0, nullptr);
  cdo_operator_add("yseassum",   FieldFunc_Sum,   0, nullptr);
  cdo_operator_add("yseasmean",  FieldFunc_Mean,  0, nullptr);
  cdo_operator_add("yseasavg",   FieldFunc_Avg,   0, nullptr);
  cdo_operator_add("yseasvar",   FieldFunc_Var,   0, nullptr);
  cdo_operator_add("yseasvar1",  FieldFunc_Var1,  0, nullptr);
  cdo_operator_add("yseasstd",   FieldFunc_Std,   0, nullptr);
  cdo_operator_add("yseasstd1",  FieldFunc_Std1,  0, nullptr);
  // clang-format on
}

class ModuleYseasstat
{

  const static int MaxSeasons = 4;
  int seas_numSets[MaxSeasons] = { 0 };
  CdiDateTime vDateTimes[MaxSeasons]{};
  FieldVector2D vars1[MaxSeasons], vars2[MaxSeasons], samp1[MaxSeasons];
  CdoStreamID streamID1;
  CdoStreamID streamID2;
  int taxisID1;
  int taxisID2;
  int vlistID1;

  VarList varList;
  Field field;

  std::vector<RecordInfo> recList;

  int VARS_MEMTYPE = 0;

  int operfunc;
  int maxrecs;

  bool lrange;
  bool lmean;
  bool lstd;
  bool lvarstd;
  bool lvars2;
  int divisor;

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
    if (taxisHasBounds(taxisID2)) taxisDeleteBounds(taxisID2);
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

        auto seas = month_to_season(decode_month(vDateTime.date));

        set_date_time(vDateTimes[seas], vDateTime);

        if (!vars1[seas].size())
          {
            fields_from_vlist(vlistID1, samp1[seas]);
            fields_from_vlist(vlistID1, vars1[seas], FIELD_VEC | VARS_MEMTYPE);
            if (lvars2) fields_from_vlist(vlistID1, vars2[seas], FIELD_VEC);
          }

        for (int recID = 0; recID < nrecs; ++recID)
          {
            int varID, levelID;
            cdo_inq_record(streamID1, &varID, &levelID);
            const auto &var = varList[varID];

            if (tsID == 0) recList[recID].set(varID, levelID);

            auto &rsamp1 = samp1[seas][varID][levelID];
            auto &rvars1 = vars1[seas][varID][levelID];

            auto numSets = seas_numSets[seas];

            if (numSets == 0)
              {
                cdo_read_record(streamID1, rvars1);
                if (lrange) field_copy(rvars1, vars2[seas][varID][levelID]);

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
                if      (lvarstd) field2_sumsumq(rvars1, vars2[seas][varID][levelID], field);
                else if (lrange)  field2_maxmin(rvars1, vars2[seas][varID][levelID], field);
                else              field2_function(rvars1, field, operfunc);
                // clang-format on
              }
          }

        if (seas_numSets[seas] == 0 && lvarstd)
          for (int recID = 0; recID < maxrecs; ++recID)
            {
              auto [varID, levelID] = recList[recID].get();
              if (varList[varID].isConstant) continue;

              field2_moq(vars2[seas][varID][levelID], vars1[seas][varID][levelID]);
            }

        seas_numSets[seas]++;
        tsID++;
      }

    for (int seas = 0; seas < MaxSeasons; ++seas)
      if (seas_numSets[seas])
        {
          auto numSets = seas_numSets[seas];
          for (int recID = 0; recID < maxrecs; ++recID)
            {
              auto [varID, levelID] = recList[recID].get();
              if (varList[varID].isConstant) continue;

              const auto &rsamp1 = samp1[seas][varID][levelID];
              auto &rvars1 = vars1[seas][varID][levelID];

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
                    field2_stdvar_func(rvars1, vars2[seas][varID][levelID], rsamp1, divisor);
                  else
                    fieldc_stdvar_func(rvars1, vars2[seas][varID][levelID], numSets, divisor);
                }
              else if (lrange) { field2_sub(rvars1, vars2[seas][varID][levelID]); }
            }

          taxisDefVdatetime(taxisID2, vDateTimes[seas]);
          cdo_def_timestep(streamID2, otsID);

          for (int recID = 0; recID < maxrecs; ++recID)
            {
              auto [varID, levelID] = recList[recID].get();
              if (otsID && varList[varID].isConstant) continue;

              auto &rvars1 = vars1[seas][varID][levelID];

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
Yseasstat(void *process)
{
  ModuleYseasstat yseasstat;
  yseasstat.init(process);
  yseasstat.run();
  yseasstat.close();

  return nullptr;
}
