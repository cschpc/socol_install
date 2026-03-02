/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Timstat    timrange        Time range
      Timstat    timmin          Time minimum
      Timstat    timmax          Time maximum
      Timstat    timsum          Time sum
      Timstat    timmean         Time mean
      Timstat    timavg          Time average
      Timstat    timvar          Time variance
      Timstat    timvar1         Time variance [Normalize by (n-1)]
      Timstat    timstd          Time standard deviation
      Timstat    timstd1         Time standard deviation [Normalize by (n-1)]
      Hourstat   hourrange       Hourly range
      Hourstat   hourmin         Hourly minimum
      Hourstat   hourmax         Hourly maximum
      Hourstat   hoursum         Hourly sum
      Hourstat   hourmean        Hourly mean
      Hourstat   houravg         Hourly average
      Hourstat   hourvar         Hourly variance
      Hourstat   hourvar1        Hourly variance [Normalize by (n-1)]
      Hourstat   hourstd         Hourly standard deviation
      Hourstat   hourstd1        Hourly standard deviation [Normalize by (n-1)]
      Daystat    dayrange        Daily range
      Daystat    daymin          Daily minimum
      Daystat    daymax          Daily maximum
      Daystat    daysum          Daily sum
      Daystat    daymean         Daily mean
      Daystat    dayavg          Daily average
      Daystat    dayvar          Daily variance
      Daystat    dayvar1         Daily variance [Normalize by (n-1)]
      Daystat    daystd          Daily standard deviation
      Daystat    daystd1         Daily standard deviation [Normalize by (n-1)]
      Monstat    monrange        Monthly range
      Monstat    monmin          Monthly minimum
      Monstat    monmax          Monthly maximum
      Monstat    monsum          Monthly sum
      Monstat    monmean         Monthly mean
      Monstat    monavg          Monthly average
      Monstat    monvar          Monthly variance
      Monstat    monvar1         Monthly variance [Normalize by (n-1)]
      Monstat    monstd          Monthly standard deviation
      Monstat    monstd1         Monthly standard deviation [Normalize by (n-1)]
      Yearstat   yearrange       Yearly range
      Yearstat   yearmin         Yearly minimum
      Yearstat   yearmax         Yearly maximum
      Yearstat   yearsum         Yearly sum
      Yearstat   yearmean        Yearly mean
      Yearstat   yearavg         Yearly average
      Yearstat   yearvar         Yearly variance
      Yearstat   yearvar1        Yearly variance [Normalize by (n-1)]
      Yearstat   yearstd         Yearly standard deviation
      Yearstat   yearstd1        Yearly standard deviation [Normalize by (n-1)]
*/

#include <cdi.h>

#include "cdo_options.h"
#include "cdo_vlist.h"
#include "process_int.h"
#include "datetime.h"
#include "printinfo.h"
#include "util_date.h"
#include "param_conversion.h"
#include "progress.h"
#include "field_functions.h"

static void
vlist_set_frequency(const int vlistID, const int compareDate)
{
  const char *freq = nullptr;
  // clang-format off
  if      (compareDate == CMP_DAY)   freq = "day";
  else if (compareDate == CMP_MONTH) freq = "mon";
  else if (compareDate == CMP_YEAR)  freq = "year";
  // clang-format on
  if (freq) cdiDefAttTxt(vlistID, CDI_GLOBAL, "frequency", (int) strlen(freq), freq);
}

static void
set_missval(Field &field, const Field &samp, int numSets, double vfrac)
{
  auto fieldsize = field.size;
  auto missval = field.missval;

  size_t irun = 0;
  for (size_t i = 0; i < fieldsize; ++i)
    {
      if ((samp.vec_d[i] / numSets) < vfrac)
        {
          field.vec_d[i] = missval;
          irun++;
        }
    }

  if (irun) field.nmiss = field_num_miss(field);
}

static void
addOperators(void)
{
  // clang-format off
  cdo_operator_add("timrange",   FieldFunc_Range,  CMP_DATE, nullptr);
  cdo_operator_add("timmin",     FieldFunc_Min,    CMP_DATE, nullptr);
  cdo_operator_add("timmax",     FieldFunc_Max,    CMP_DATE, nullptr);
  cdo_operator_add("timminidx",  FieldFunc_Minidx, CMP_DATE, nullptr);
  cdo_operator_add("timmaxidx",  FieldFunc_Maxidx, CMP_DATE, nullptr);
  cdo_operator_add("timsum",     FieldFunc_Sum,    CMP_DATE, nullptr);
  cdo_operator_add("timmean",    FieldFunc_Mean,   CMP_DATE, nullptr);
  cdo_operator_add("timavg",     FieldFunc_Avg,    CMP_DATE, nullptr);
  cdo_operator_add("timvar",     FieldFunc_Var,    CMP_DATE, nullptr);
  cdo_operator_add("timvar1",    FieldFunc_Var1,   CMP_DATE, nullptr);
  cdo_operator_add("timstd",     FieldFunc_Std,    CMP_DATE, nullptr);
  cdo_operator_add("timstd1",    FieldFunc_Std1,   CMP_DATE, nullptr);
  cdo_operator_add("yearrange",  FieldFunc_Range,  CMP_YEAR, nullptr);
  cdo_operator_add("yearmin",    FieldFunc_Min,    CMP_YEAR, nullptr);
  cdo_operator_add("yearmax",    FieldFunc_Max,    CMP_YEAR, nullptr);
  cdo_operator_add("yearminidx", FieldFunc_Minidx, CMP_YEAR, nullptr);
  cdo_operator_add("yearmaxidx", FieldFunc_Maxidx, CMP_YEAR, nullptr);
  cdo_operator_add("yearsum",    FieldFunc_Sum,    CMP_YEAR, nullptr);
  cdo_operator_add("yearmean",   FieldFunc_Mean,   CMP_YEAR, nullptr);
  cdo_operator_add("yearavg",    FieldFunc_Avg,    CMP_YEAR, nullptr);
  cdo_operator_add("yearvar",    FieldFunc_Var,    CMP_YEAR, nullptr);
  cdo_operator_add("yearvar1",   FieldFunc_Var1,   CMP_YEAR, nullptr);
  cdo_operator_add("yearstd",    FieldFunc_Std,    CMP_YEAR, nullptr);
  cdo_operator_add("yearstd1",   FieldFunc_Std1,   CMP_YEAR, nullptr);
  cdo_operator_add("monrange",   FieldFunc_Range,  CMP_MONTH, nullptr);
  cdo_operator_add("monmin",     FieldFunc_Min,    CMP_MONTH, nullptr);
  cdo_operator_add("monmax",     FieldFunc_Max,    CMP_MONTH, nullptr);
  cdo_operator_add("monsum",     FieldFunc_Sum,    CMP_MONTH, nullptr);
  cdo_operator_add("monmean",    FieldFunc_Mean,   CMP_MONTH, nullptr);
  cdo_operator_add("monavg",     FieldFunc_Avg,    CMP_MONTH, nullptr);
  cdo_operator_add("monvar",     FieldFunc_Var,    CMP_MONTH, nullptr);
  cdo_operator_add("monvar1",    FieldFunc_Var1,   CMP_MONTH, nullptr);
  cdo_operator_add("monstd",     FieldFunc_Std,    CMP_MONTH, nullptr);
  cdo_operator_add("monstd1",    FieldFunc_Std1,   CMP_MONTH, nullptr);
  cdo_operator_add("dayrange",   FieldFunc_Range,  CMP_DAY, nullptr);
  cdo_operator_add("daymin",     FieldFunc_Min,    CMP_DAY, nullptr);
  cdo_operator_add("daymax",     FieldFunc_Max,    CMP_DAY, nullptr);
  cdo_operator_add("daysum",     FieldFunc_Sum,    CMP_DAY, nullptr);
  cdo_operator_add("daymean",    FieldFunc_Mean,   CMP_DAY, nullptr);
  cdo_operator_add("dayavg",     FieldFunc_Avg,    CMP_DAY, nullptr);
  cdo_operator_add("dayvar",     FieldFunc_Var,    CMP_DAY, nullptr);
  cdo_operator_add("dayvar1",    FieldFunc_Var1,   CMP_DAY, nullptr);
  cdo_operator_add("daystd",     FieldFunc_Std,    CMP_DAY, nullptr);
  cdo_operator_add("daystd1",    FieldFunc_Std1,   CMP_DAY, nullptr);
  cdo_operator_add("hourrange",  FieldFunc_Range,  CMP_HOUR, nullptr);
  cdo_operator_add("hourmin",    FieldFunc_Min,    CMP_HOUR, nullptr);
  cdo_operator_add("hourmax",    FieldFunc_Max,    CMP_HOUR, nullptr);
  cdo_operator_add("hoursum",    FieldFunc_Sum,    CMP_HOUR, nullptr);
  cdo_operator_add("hourmean",   FieldFunc_Mean,   CMP_HOUR, nullptr);
  cdo_operator_add("houravg",    FieldFunc_Avg,    CMP_HOUR, nullptr);
  cdo_operator_add("hourvar",    FieldFunc_Var,    CMP_HOUR, nullptr);
  cdo_operator_add("hourvar1",   FieldFunc_Var1,   CMP_HOUR, nullptr);
  cdo_operator_add("hourstd",    FieldFunc_Std,    CMP_HOUR, nullptr);
  cdo_operator_add("hourstd1",   FieldFunc_Std1,   CMP_HOUR, nullptr);
  // clang-format on
}

void *
Timstat(void *argument)
{
  const TimeStat timestat_date = TimeStat::MEAN;
  CdiDateTime vDateTime0{};
  CdiDateTime vDateTimeN{};
  CdoStreamID streamID3;
  int vlistID3, taxisID3 = -1;
  bool lvfrac = false;
  double vfrac = 1;

  cdo_initialize(argument);

  addOperators();

  auto operatorID = cdo_operator_id();
  auto operfunc = cdo_operator_f1(operatorID);
  auto compareDate = cdo_operator_f2(operatorID);

  auto lminidx = (operfunc == FieldFunc_Minidx);
  auto lmaxidx = (operfunc == FieldFunc_Maxidx);
  auto lrange = (operfunc == FieldFunc_Range);
  auto lmean = (operfunc == FieldFunc_Mean || operfunc == FieldFunc_Avg);
  auto lstd = (operfunc == FieldFunc_Std || operfunc == FieldFunc_Std1);
  auto lvarstd = (lstd || operfunc == FieldFunc_Var || operfunc == FieldFunc_Var1);
  auto lvars2 = (lvarstd || lrange || lminidx || lmaxidx);
  const int divisor = (operfunc == FieldFunc_Std1 || operfunc == FieldFunc_Var1);

  auto field2_stdvar_func = lstd ? field2_std : field2_var;
  auto fieldc_stdvar_func = lstd ? fieldc_std : fieldc_var;

  if (operfunc == FieldFunc_Mean)
    {
      auto oargc = cdo_operator_argc();
      if (oargc == 1)
        {
          lvfrac = true;
          vfrac = parameter_to_double(cdo_operator_argv(0));
          if (Options::cdoVerbose) cdo_print("Set vfrac to %g", vfrac);
          if (vfrac < 0 || vfrac > 1) cdo_abort("vfrac out of range!");
        }
      else if (oargc > 1)
        cdo_abort("Too many arguments!");
    }
  else { operator_check_argc(0); }

  auto streamID1 = cdo_open_read(0);

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto vlistID2 = vlistDuplicate(vlistID1);

  vlist_define_timestep_type(vlistID2, operfunc);

  vlistDefNtsteps(vlistID2, (compareDate == CMP_DATE) ? 1 : -1);

  auto taxisID1 = vlistInqTaxis(vlistID1);
  auto taxisID2 = taxisDuplicate(taxisID1);
  taxisWithBounds(taxisID2);
  if (taxisInqType(taxisID2) == TAXIS_FORECAST) taxisDefType(taxisID2, TAXIS_RELATIVE);
  vlistDefTaxis(vlistID2, taxisID2);

  auto nvars = vlistNvars(vlistID1);

  if (lminidx || lmaxidx)
    for (int varID = 0; varID < nvars; ++varID) vlistDefVarDatatype(vlistID2, varID, CDI_DATATYPE_INT32);

  vlist_set_frequency(vlistID2, compareDate);

  auto streamID2 = cdo_open_write(1);
  cdo_def_vlist(streamID2, vlistID2);

  if (Options::cdoDiag)
    {
      char filename[8192];
      strcpy(filename, cdo_operator_name(operatorID));
      strcat(filename, "_");
      strcat(filename, cdo_get_stream_name(1));
      streamID3 = cdo_open_write(filename);

      vlistID3 = vlistDuplicate(vlistID1);

      for (int varID = 0; varID < nvars; ++varID)
        {
          vlistDefVarDatatype(vlistID3, varID, CDI_DATATYPE_INT32);
          vlistDefVarMissval(vlistID3, varID, -1);
          cdiDefKeyString(vlistID3, varID, CDI_KEY_UNITS, "");
          cdiDeleteKey(vlistID3, varID, CDI_KEY_ADDOFFSET);
          cdiDeleteKey(vlistID3, varID, CDI_KEY_SCALEFACTOR);
        }

      taxisID3 = taxisDuplicate(taxisID1);
      taxisWithBounds(taxisID3);
      vlistDefTaxis(vlistID3, taxisID3);

      cdo_def_vlist(streamID3, vlistID3);
    }

  auto maxrecs = vlistNrecs(vlistID1);
  std::vector<RecordInfo> recList(maxrecs);

  DateTimeList dtlist;
  dtlist.set_stat(timestat_date);
  dtlist.set_calendar(taxisInqCalendar(taxisID1));

  VarList varList;
  varListInit(varList, vlistID1);

  int VARS_MEMTYPE = 0;
  if ((operfunc == FieldFunc_Min) || (operfunc == FieldFunc_Max)) VARS_MEMTYPE = FIELD_NAT;
  // if ((Options::CDO_???(--single) == MemType::Float) && (operfunc == FieldFunc_Mean)) VARS_MEMTYPE = FIELD_NAT;
  // if (Options::CDO_Memtype == MemType::Float) VARS_MEMTYPE = FIELD_FLT;
  if (Options::cdoDiag || (lvfrac && operfunc == FieldFunc_Mean)) VARS_MEMTYPE = FIELD_DBL;

  Field field;
  Varray<double> samp;

  FieldVector2D samp1, vars1, vars2;
  fields_from_vlist(vlistID1, samp1);
  fields_from_vlist(vlistID1, vars1, FIELD_VEC | VARS_MEMTYPE);
  if (lvars2) fields_from_vlist(vlistID1, vars2, FIELD_VEC);

  auto ntsteps1 = vlistNtsteps(vlistID1);

  if (!Options::cdoVerbose && ntsteps1 > 1) progress::init();

  int tsID = 0;
  int otsID = 0;
  while (true)
    {
      int numSets = 0;
      int nrecs = 0;
      while (true)
        {
          nrecs = cdo_stream_inq_timestep(streamID1, tsID);
          if (nrecs == 0) break;

          if (!Options::cdoVerbose && ntsteps1 > 1) progress::update(0, 1, (tsID + 1.0) / ntsteps1);

          dtlist.taxis_inq_timestep(taxisID1, numSets);
          auto vDateTime = dtlist.get_vDateTime(numSets);

          if (numSets == 0) vDateTime0 = vDateTime;

          if (date_is_neq(vDateTime, vDateTime0, compareDate))
            {
              cdo_add_steps(-1);
              break;
            }

          for (int recID = 0; recID < nrecs; ++recID)
            {
              int varID, levelID;
              cdo_inq_record(streamID1, &varID, &levelID);

              if (tsID == 0) recList[recID].set(varID, levelID);

              auto &rsamp1 = samp1[varID][levelID];
              auto &rvars1 = vars1[varID][levelID];

              if (numSets == 0)
                {
                  cdo_read_record(streamID1, rvars1);
                  if (lrange || lminidx || lmaxidx)
                    {
                      vars2[varID][levelID].nmiss = rvars1.nmiss;
                      vars2[varID][levelID].vec_d = rvars1.vec_d;
                    }

                  if (lminidx || lmaxidx) field_fill(rvars1, 0.0);

                  if (rvars1.nmiss || !rsamp1.empty())
                    {
                      if (rsamp1.empty()) rsamp1.resize(rvars1.size);
                      field2_vinit(rsamp1, rvars1);
                    }
                }
              else
                {
                  field.init(varList[varID]);
                  cdo_read_record(streamID1, field);

                  if (field.nmiss || !rsamp1.empty())
                    {
                      if (rsamp1.empty()) rsamp1.resize(rvars1.size, numSets);
                      field2_vincr(rsamp1, field);
                    }

                  // clang-format off
                  if      (lvarstd) field2_sumsumq(rvars1, vars2[varID][levelID], field);
                  else if (lrange)  field2_maxmin(rvars1, vars2[varID][levelID], field);
                  else if (lminidx) field2_minidx(rvars1, vars2[varID][levelID], field, numSets);
                  else if (lmaxidx) field2_maxidx(rvars1, vars2[varID][levelID], field, numSets);
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

          vDateTimeN = vDateTime;
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

      if (Options::cdoVerbose) cdo_print("%s  vfrac = %g, numSets = %d", datetime_to_string(vDateTimeN), vfrac, numSets);

      if (lvfrac && operfunc == FieldFunc_Mean)
        for (int recID = 0; recID < maxrecs; ++recID)
          {
            auto [varID, levelID] = recList[recID].get();
            if (varList[varID].isConstant) continue;

            const auto &rsamp1 = samp1[varID][levelID];

            if (!rsamp1.empty()) set_missval(vars1[varID][levelID], rsamp1, numSets, vfrac);
          }

      dtlist.stat_taxis_def_timestep(taxisID2, numSets);
      cdo_def_timestep(streamID2, otsID);

      if (Options::cdoDiag)
        {
          dtlist.stat_taxis_def_timestep(taxisID3, numSets);
          cdo_def_timestep(streamID3, otsID);
        }

      for (int recID = 0; recID < maxrecs; ++recID)
        {
          auto [varID, levelID] = recList[recID].get();
          if (otsID && varList[varID].isConstant) continue;

          auto &rvars1 = vars1[varID][levelID];

          cdo_def_record(streamID2, varID, levelID);
          cdo_write_record(streamID2, rvars1);

          if (Options::cdoDiag)
            {
              auto &rsamp1 = samp1[varID][levelID];
              samp.resize(field.size);
              if (!rsamp1.empty())
                samp = rsamp1.vec_d;
              else
                varray_fill(samp, (double) numSets);

              cdo_def_record(streamID3, varID, levelID);
              cdo_write_record(streamID3, samp.data(), 0);
            }
        }

      if (nrecs == 0) break;
      otsID++;
    }

  if (!Options::cdoVerbose && ntsteps1 > 1) progress::update(0, 1, 1);

  if (Options::cdoDiag) cdo_stream_close(streamID3);
  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  cdo_finish();

  return nullptr;
}
