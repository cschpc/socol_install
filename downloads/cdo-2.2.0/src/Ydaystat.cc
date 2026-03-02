/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Ydaystat   ydayrange       Multi-year daily range
      Ydaystat   ydaymin         Multi-year daily minimum
      Ydaystat   ydaymax         Multi-year daily maximum
      Ydaystat   ydaysum         Multi-year daily sum
      Ydaystat   ydaymean        Multi-year daily mean
      Ydaystat   ydayavg         Multi-year daily average
      Ydaystat   ydayvar         Multi-year daily variance
      Ydaystat   ydayvar1        Multi-year daily variance [Normalize by (n-1)]
      Ydaystat   ydaystd         Multi-year daily standard deviation
      Ydaystat   ydaystd1        Multi-year daily standard deviation [Normalize by (n-1)]
*/

#include <cdi.h>

#include "cdo_options.h"
#include "cdo_vlist.h"
#include "datetime.h"
#include "process_int.h"
#include "param_conversion.h"
#include "pmlist.h"
#include "printinfo.h"
#include "field_functions.h"

static int yearMode = 0;

static void
setParameter(void)
{
  const auto pargc = cdo_operator_argc();
  if (pargc)
    {
      const auto &pargv = cdo_get_oper_argv();

      KVList kvlist;
      kvlist.name = cdo_module_name();
      if (kvlist.parse_arguments(pargc, pargv) != 0) cdo_abort("Parse error!");
      if (Options::cdoVerbose) kvlist.print();

      for (const auto &kv : kvlist)
        {
          const auto &key = kv.key;
          if (kv.nvalues > 1) cdo_abort("Too many values for parameter key >%s<!", key);
          if (kv.nvalues < 1) cdo_abort("Missing value for parameter key >%s<!", key);
          const auto &value = kv.values[0];

          if (key == "yearMode")
            yearMode = parameter_to_int(value);
          else
            cdo_abort("Invalid parameter key >%s<!", key);
        }
    }
}

static void
addOperators(void)
{
  // clang-format off
  cdo_operator_add("ydayrange", FieldFunc_Range, 0, nullptr);
  cdo_operator_add("ydaymin",   FieldFunc_Min,   0, nullptr);
  cdo_operator_add("ydaymax",   FieldFunc_Max,   0, nullptr);
  cdo_operator_add("ydaysum",   FieldFunc_Sum,   0, nullptr);
  cdo_operator_add("ydaymean",  FieldFunc_Mean,  0, nullptr);
  cdo_operator_add("ydayavg",   FieldFunc_Avg,   0, nullptr);
  cdo_operator_add("ydayvar",   FieldFunc_Var,   0, nullptr);
  cdo_operator_add("ydayvar1",  FieldFunc_Var1,  0, nullptr);
  cdo_operator_add("ydaystd",   FieldFunc_Std,   0, nullptr);
  cdo_operator_add("ydaystd1",  FieldFunc_Std1,  0, nullptr);
  // clang-format on
}
class ModuleYdaystat
{
private:
  static const int MaxDays = 373;

  CdoStreamID streamID1;
  CdoStreamID streamID2;

  int taxisID1;
  int taxisID2;
  int vlistID1;

  int dayOfYear_numSets[MaxDays]{ 0 };
  int VARS_MEMTYPE = 0;

  CdiDateTime vDateTimes[MaxDays]{};
  FieldVector2D vars1[MaxDays], vars2[MaxDays], samp1[MaxDays];

  bool lmean;
  bool lrange;
  bool lstd;
  bool lvars2;
  bool lvarstd;

  int divisor;
  int maxrecs;
  int operfunc;

  std::vector<RecordInfo> recList;
  VarList varList;
  Field field;

public:
  void
  init(void *process)
  {

    cdo_initialize(process);

    setParameter();

    addOperators();

    auto operatorID = cdo_operator_id();
    operfunc = cdo_operator_f1(operatorID);

    lrange = (operfunc == FieldFunc_Range);
    lmean = (operfunc == FieldFunc_Mean || operfunc == FieldFunc_Avg);
    lstd = (operfunc == FieldFunc_Std || operfunc == FieldFunc_Std1);
    lvarstd = (lstd || operfunc == FieldFunc_Var || operfunc == FieldFunc_Var1);
    lvars2 = (lvarstd || lrange);
    divisor = (operfunc == FieldFunc_Std1 || operfunc == FieldFunc_Var1);

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

        if (Options::cdoVerbose) cdo_print("process timestep: %d %s", tsID + 1, datetime_to_string(vDateTime));

        auto dayOfYear = decode_day_of_year(vDateTime.date);
        if (dayOfYear < 0 || dayOfYear >= MaxDays)
          cdo_abort("Day of year %d out of range (%s)!", dayOfYear, datetime_to_string(vDateTime));

        vDateTimes[dayOfYear] = vDateTime;

        if (!vars1[dayOfYear].size())
          {
            fields_from_vlist(vlistID1, samp1[dayOfYear]);
            fields_from_vlist(vlistID1, vars1[dayOfYear], FIELD_VEC | VARS_MEMTYPE);
            if (lvars2) fields_from_vlist(vlistID1, vars2[dayOfYear], FIELD_VEC);
          }

        for (int recID = 0; recID < nrecs; ++recID)
          {
            int varID, levelID;
            cdo_inq_record(streamID1, &varID, &levelID);
            const auto &var = varList[varID];

            if (tsID == 0) recList[recID].set(varID, levelID);

            auto &rsamp1 = samp1[dayOfYear][varID][levelID];
            auto &rvars1 = vars1[dayOfYear][varID][levelID];

            auto numSets = dayOfYear_numSets[dayOfYear];

            if (numSets == 0)
              {
                cdo_read_record(streamID1, rvars1);
                if (lrange)
                  {
                    vars2[dayOfYear][varID][levelID].nmiss = rvars1.nmiss;
                    vars2[dayOfYear][varID][levelID].vec_d = rvars1.vec_d;
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
                if      (lvarstd) field2_sumsumq(rvars1, vars2[dayOfYear][varID][levelID], field);
                else if (lrange)  field2_maxmin(rvars1, vars2[dayOfYear][varID][levelID], field);
                else              field2_function(rvars1, field, operfunc);
                // clang-format on
              }
          }

        if (dayOfYear_numSets[dayOfYear] == 0 && lvarstd)
          for (int recID = 0; recID < maxrecs; ++recID)
            {
              auto [varID, levelID] = recList[recID].get();
              if (varList[varID].isConstant) continue;

              field2_moq(vars2[dayOfYear][varID][levelID], vars1[dayOfYear][varID][levelID]);
            }

        dayOfYear_numSets[dayOfYear]++;
        tsID++;
      }

    // set the year to the minimum of years found on output timestep
    if (yearMode)
      {
        int outyear = 1e9;
        for (int dayOfYear = 0; dayOfYear < MaxDays; dayOfYear++)
          if (dayOfYear_numSets[dayOfYear])
            {
              auto year = vDateTimes[dayOfYear].date.year;
              if (year < outyear) outyear = year;
            }
        for (int dayOfYear = 0; dayOfYear < MaxDays; dayOfYear++)
          if (dayOfYear_numSets[dayOfYear])
            {
              auto year = vDateTimes[dayOfYear].date.year;
              if (year > outyear) vDateTimes[dayOfYear].date.year = outyear;
            }
      }

    for (int dayOfYear = 0; dayOfYear < MaxDays; dayOfYear++)
      if (dayOfYear_numSets[dayOfYear])
        {
          auto numSets = dayOfYear_numSets[dayOfYear];
          for (int recID = 0; recID < maxrecs; ++recID)
            {
              auto [varID, levelID] = recList[recID].get();
              if (varList[varID].isConstant) continue;

              const auto &rsamp1 = samp1[dayOfYear][varID][levelID];
              auto &rvars1 = vars1[dayOfYear][varID][levelID];

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
                    field2_stdvar_func(rvars1, vars2[dayOfYear][varID][levelID], rsamp1, divisor);
                  else
                    fieldc_stdvar_func(rvars1, vars2[dayOfYear][varID][levelID], numSets, divisor);
                }
              else if (lrange) { field2_sub(rvars1, vars2[dayOfYear][varID][levelID]); }
            }

          taxisDefVdatetime(taxisID2, vDateTimes[dayOfYear]);
          cdo_def_timestep(streamID2, otsID);

          for (int recID = 0; recID < maxrecs; ++recID)
            {
              auto [varID, levelID] = recList[recID].get();
              if (otsID && varList[varID].isConstant) continue;

              auto &rvars1 = vars1[dayOfYear][varID][levelID];

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
Ydaystat(void *process)
{
  ModuleYdaystat ydaystat;
  ydaystat.init(process);
  ydaystat.run();
  ydaystat.close();
  return nullptr;
}
