/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Timstat    timmin          Time minimum
      Timstat    timmax          Time maximum
      Timstat    timsum          Time sum
      Timstat    timmean         Time mean
      Timstat    timavg          Time average
      Timstat    timvar          Time variance
      Timstat    timvar1         Time variance [Normalize by (n-1)]
      Timstat    timstd          Time standard deviation
      Timstat    timstd1         Time standard deviation [Normalize by (n-1)]
      Hourstat   hourmin         Hourly minimum
      Hourstat   hourmax         Hourly maximum
      Hourstat   hoursum         Hourly sum
      Hourstat   hourmean        Hourly mean
      Hourstat   houravg         Hourly average
      Hourstat   hourvar         Hourly variance
      Hourstat   hourvar1        Hourly variance [Normalize by (n-1)]
      Hourstat   hourstd         Hourly standard deviation
      Hourstat   hourstd1        Hourly standard deviation [Normalize by (n-1)]
      Daystat    daymin          Daily minimum
      Daystat    daymax          Daily maximum
      Daystat    daysum          Daily sum
      Daystat    daymean         Daily mean
      Daystat    dayavg          Daily average
      Daystat    dayvar          Daily variance
      Daystat    dayvar1         Daily variance [Normalize by (n-1)]
      Daystat    daystd          Daily standard deviation
      Daystat    daystd1         Daily standard deviation [Normalize by (n-1)]
      Monstat    monmin          Monthly minimum
      Monstat    monmax          Monthly maximum
      Monstat    monsum          Monthly sum
      Monstat    monmean         Monthly mean
      Monstat    monavg          Monthly average
      Monstat    monvar          Monthly variance
      Monstat    monvar1         Monthly variance [Normalize by (n-1)]
      Monstat    monstd          Monthly standard deviation
      Monstat    monstd1         Monthly standard deviation [Normalize by (n-1)]
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

#include <utility>

#include "cdo_options.h"
#include "cdo_vlist.h"
#include "process_int.h"
#include "cdo_task.h"
#include "timer.h"
#include "datetime.h"
#include "printinfo.h"
#include "util_date.h"
#include "param_conversion.h"
#include "field_functions.h"

#define USE_CDI_STREAM 1

struct ReadArguments
{
  int tsIDnext;
#ifdef USE_CDI_STREAM
  int streamID;
#else
  CdoStreamID streamID;
#endif
  int nrecs;
  RecordInfo *recList;
  FieldVector2D &vars;

  ReadArguments(FieldVector2D &input_vars) : tsIDnext(0), nrecs(0), recList(nullptr), vars(input_vars){};
};

static int num_recs = 0;

static void *
cdoReadTimestep(void *rarg)
{
  ReadArguments *readarg = (ReadArguments *) rarg;
  auto &input_vars = readarg->vars;
  RecordInfo *recList = readarg->recList;
  const auto streamID = readarg->streamID;
  int tsIDnext = readarg->tsIDnext;
  int nrecs = readarg->nrecs;

  // timer_start(timer_read);

  for (int recID = 0; recID < nrecs; ++recID)
    {
      int varID, levelID;
#ifdef USE_CDI_STREAM
      streamInqRecord(streamID, &varID, &levelID);
#else
      cdo_inq_record(streamID, &varID, &levelID);
#endif

      if (tsIDnext == 1 && recList)
        {
          recList[recID].varID = varID;
          recList[recID].levelID = levelID;
        }

      size_t nmiss;

#ifdef USE_CDI_STREAM
      if (Options::CDO_Memtype == MemType::Float)
        streamReadRecordF(streamID, input_vars[varID][levelID].vec_f.data(), &nmiss);
      else
        streamReadRecord(streamID, input_vars[varID][levelID].vec_d.data(), &nmiss);
#else
      if (Options::CDO_Memtype == MemType::Float)
        cdo_read_record_f(streamID, input_vars[varID][levelID].vec_f.data(), &nmiss);
      else
        cdo_read_record(streamID, input_vars[varID][levelID].vec_d.data(), &nmiss);
#endif
      input_vars[varID][levelID].nmiss = nmiss;
    }

    // timer_stop(timer_read);

#ifdef USE_CDI_STREAM
  num_recs = streamInqTimestep(streamID, tsIDnext);
#else
  num_recs = cdo_stream_inq_timestep(streamID, tsIDnext);
#endif

  return ((void *) &num_recs);
}

static void
vlistSetFrequency(int vlistID, int compareDate)
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
addOperators(void)
{
  // clang-format off
  cdo_operator_add("xtimmin",    FieldFunc_Min,   CMP_DATE, nullptr);
  cdo_operator_add("xtimmax",    FieldFunc_Max,   CMP_DATE, nullptr);
  cdo_operator_add("xtimsum",    FieldFunc_Sum,   CMP_DATE, nullptr);
  cdo_operator_add("xtimmean",   FieldFunc_Mean,  CMP_DATE, nullptr);
  cdo_operator_add("xtimavg",    FieldFunc_Avg,   CMP_DATE, nullptr);
  cdo_operator_add("xtimvar",    FieldFunc_Var,   CMP_DATE, nullptr);
  cdo_operator_add("xtimvar1",   FieldFunc_Var1,  CMP_DATE, nullptr);
  cdo_operator_add("xtimstd",    FieldFunc_Std,   CMP_DATE, nullptr);
  cdo_operator_add("xtimstd1",   FieldFunc_Std1,  CMP_DATE, nullptr);
  cdo_operator_add("xyearmin",   FieldFunc_Min,   CMP_YEAR, nullptr);
  cdo_operator_add("xyearmax",   FieldFunc_Max,   CMP_YEAR, nullptr);
  cdo_operator_add("xyearsum",   FieldFunc_Sum,   CMP_YEAR, nullptr);
  cdo_operator_add("xyearmean",  FieldFunc_Mean,  CMP_YEAR, nullptr);
  cdo_operator_add("xyearavg",   FieldFunc_Avg,   CMP_YEAR, nullptr);
  cdo_operator_add("xyearvar",   FieldFunc_Var,   CMP_YEAR, nullptr);
  cdo_operator_add("xyearvar1",  FieldFunc_Var1,  CMP_YEAR, nullptr);
  cdo_operator_add("xyearstd",   FieldFunc_Std,   CMP_YEAR, nullptr);
  cdo_operator_add("xyearstd1",  FieldFunc_Std1,  CMP_YEAR, nullptr);
  cdo_operator_add("xmonmin",    FieldFunc_Min,   CMP_MONTH, nullptr);
  cdo_operator_add("xmonmax",    FieldFunc_Max,   CMP_MONTH, nullptr);
  cdo_operator_add("xmonsum",    FieldFunc_Sum,   CMP_MONTH, nullptr);
  cdo_operator_add("xmonmean",   FieldFunc_Mean,  CMP_MONTH, nullptr);
  cdo_operator_add("xmonavg",    FieldFunc_Avg,   CMP_MONTH, nullptr);
  cdo_operator_add("xmonvar",    FieldFunc_Var,   CMP_MONTH, nullptr);
  cdo_operator_add("xmonvar1",   FieldFunc_Var1,  CMP_MONTH, nullptr);
  cdo_operator_add("xmonstd",    FieldFunc_Std,   CMP_MONTH, nullptr);
  cdo_operator_add("xmonstd1",   FieldFunc_Std1,  CMP_MONTH, nullptr);
  // clang-format on
}

void *
XTimstat(void *process)
{
  TimeStat timestat_date = TimeStat::MEAN;
  CdiDateTime vDateTime0{};
  CdiDateTime vDateTimeN{};
  CdoStreamID streamID3 = CDO_STREAM_UNDEF;
  int vlistID3, taxisID3 = -1;
  bool lvfrac = false;
  double vfrac = 1;

  cdo_initialize(process);

  addOperators();

  const auto operatorID = cdo_operator_id();
  auto operfunc = cdo_operator_f1(operatorID);
  const auto compareDate = cdo_operator_f2(operatorID);

  const auto lmean = (operfunc == FieldFunc_Mean || operfunc == FieldFunc_Avg);
  const auto lstd = (operfunc == FieldFunc_Std || operfunc == FieldFunc_Std1);
  auto lvarstd = (lstd || operfunc == FieldFunc_Var || operfunc == FieldFunc_Var1);
  const auto lvars2 = lvarstd;
  int divisor = (operfunc == FieldFunc_Std1 || operfunc == FieldFunc_Var1);

  auto field2_stdvar_func = lstd ? field2_std : field2_var;
  auto fieldc_stdvar_func = lstd ? fieldc_std : fieldc_var;

  if (operfunc == FieldFunc_Mean)
    {
      const auto oargc = cdo_operator_argc();
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

#ifdef USE_CDI_STREAM
  const auto streamID1 = streamOpenRead(cdo_get_stream_name(0));
  auto vlistID1 = streamInqVlist(streamID1);
#else
  const auto streamID1 = cdo_open_read(0);
  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
#endif

  const auto vlistID2 = vlistDuplicate(vlistID1);
  vlist_define_timestep_type(vlistID2, operfunc);

  vlistDefNtsteps(vlistID2, (compareDate == CMP_DATE) ? 1 : -1);

  const auto taxisID1 = vlistInqTaxis(vlistID1);
  const auto taxisID2 = taxisDuplicate(taxisID1);
  taxisWithBounds(taxisID2);
  if (taxisInqType(taxisID2) == TAXIS_FORECAST) taxisDefType(taxisID2, TAXIS_RELATIVE);
  vlistDefTaxis(vlistID2, taxisID2);

  const auto nvars = vlistNvars(vlistID1);

  vlistSetFrequency(vlistID2, compareDate);

  const auto streamID2 = cdo_open_write(1);
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

  DateTimeList dtlist;
  dtlist.set_stat(timestat_date);
  dtlist.set_calendar(taxisInqCalendar(taxisID1));

  int FIELD_MEMTYPE = 0;
  if (operfunc == FieldFunc_Mean && Options::CDO_Memtype == MemType::Float) FIELD_MEMTYPE = FIELD_FLT;
  FieldVector3D input_vars(2);
  fields_from_vlist(vlistID1, input_vars[0], FIELD_VEC | FIELD_MEMTYPE);
  fields_from_vlist(vlistID1, input_vars[1], FIELD_VEC | FIELD_MEMTYPE);
  FieldVector2D samp1, vars1, vars2;
  fields_from_vlist(vlistID1, samp1);
  fields_from_vlist(vlistID1, vars1, FIELD_VEC);
  if (lvars2) fields_from_vlist(vlistID1, vars2, FIELD_VEC);

  int curFirst = 0;
  int curSecond = 1;
  (void) curSecond;

  ReadArguments readarg(input_vars[curFirst]);
  readarg.streamID = streamID1;

  bool lparallelread = Options::CDO_Parallel_Read > 0;
  bool ltsfirst = true;
  cdo::Task *read_task = nullptr;
  void *readresult = nullptr;

  if (lparallelread) read_task = new cdo::Task;

  int tsID = 0;
  int otsID = 0;
#ifdef USE_CDI_STREAM
  auto nrecs = streamInqTimestep(streamID1, tsID);
#else
  auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
#endif
  int maxrecs = nrecs;
  std::vector<RecordInfo> recList(maxrecs);

  tsID++;
  while (true)
    {
      int numSets = 0;
      while (nrecs > 0)
        {
          dtlist.taxis_inq_timestep(taxisID1, numSets);
          const auto vDateTime = dtlist.get_vDateTime(numSets);

          if (numSets == 0) vDateTime0 = vDateTime;

          if (date_is_neq(vDateTime, vDateTime0, compareDate))
            {
              cdo_add_steps(-1);
              break;
            }

          readarg.tsIDnext = tsID;
          readarg.nrecs = nrecs;
          readarg.recList = recList.data();
          readarg.vars = input_vars[curFirst];

          if (ltsfirst || !lparallelread)
            {
              ltsfirst = false;
              readresult = cdoReadTimestep(&readarg);
            }
          else { readresult = read_task->wait(); }

          nrecs = *(int *) readresult;

          // std::swap(curFirst, curSecond);

          if (nrecs && lparallelread)
            {
              readarg.vars = input_vars[curFirst];
              readarg.tsIDnext = tsID + 1;
              read_task->start(cdoReadTimestep, &readarg);
            }

          if (numSets == 0)
            {
#ifdef _OPENMP
#pragma omp parallel for default(shared) if (maxrecs > 1)
#endif
              for (int recID = 0; recID < maxrecs; ++recID)
                {
                  const auto varID = recList[recID].varID;
                  const auto levelID = recList[recID].levelID;

                  auto &rsamp1 = samp1[varID][levelID];
                  auto &rvars1 = vars1[varID][levelID];
                  auto &rinput_var = input_vars[curFirst][varID][levelID];

                  const auto nmiss = rinput_var.nmiss;

                  field_copy(rinput_var, rvars1);
                  rvars1.nmiss = nmiss;
                  if (nmiss || !rsamp1.empty())
                    {
                      const auto fieldsize = rvars1.size;
                      if (rsamp1.empty()) rsamp1.resize(fieldsize);

                      for (size_t i = 0; i < fieldsize; ++i) rsamp1.vec_d[i] = !DBL_IS_EQUAL(rvars1.vec_d[i], rvars1.missval);
                    }
                }
            }
          else
            {
#ifdef _OPENMP
#pragma omp parallel for default(shared) if (maxrecs > 1)
#endif
              for (int recID = 0; recID < maxrecs; ++recID)
                {
                  const auto varID = recList[recID].varID;
                  const auto levelID = recList[recID].levelID;

                  auto &rsamp1 = samp1[varID][levelID];
                  auto &rvars1 = vars1[varID][levelID];
                  auto &rinput_var = input_vars[curFirst][varID][levelID];

                  const auto nmiss = rinput_var.nmiss;

                  if (nmiss || !rsamp1.empty())
                    {
                      const auto fieldsize = rvars1.size;
                      if (rsamp1.empty()) rsamp1.resize(fieldsize, numSets);

                      for (size_t i = 0; i < fieldsize; ++i)
                        if (!DBL_IS_EQUAL(rinput_var.vec_d[i], rvars1.missval)) rsamp1.vec_d[i]++;
                    }

                  // clang-format off
                  if (lvarstd) field2_sumsumq(rvars1, vars2[varID][levelID], rinput_var);
                  else         field2_function(rvars1, rinput_var, operfunc);
                  // clang-format on
                }
            }

          if (numSets == 0 && lvarstd)
            for (int recID = 0; recID < maxrecs; ++recID)
              {
                const auto varID = recList[recID].varID;
                const auto levelID = recList[recID].levelID;

                if (vlistInqVarTimetype(vlistID1, varID) == TIME_CONSTANT) continue;

                field2_moq(vars2[varID][levelID], vars1[varID][levelID]);
              }

          vDateTimeN = vDateTime;
          numSets++;
          tsID++;
        }

      if (nrecs == 0 && numSets == 0) break;

      if (lmean)
        {
#ifdef _OPENMP
#pragma omp parallel for default(shared) if (maxrecs > 1)
#endif
          for (int recID = 0; recID < maxrecs; ++recID)
            {
              const auto varID = recList[recID].varID;
              const auto levelID = recList[recID].levelID;

              if (vlistInqVarTimetype(vlistID1, varID) == TIME_CONSTANT) continue;

              if (!samp1[varID][levelID].empty())
                field2_div(vars1[varID][levelID], samp1[varID][levelID]);
              else
                fieldc_div(vars1[varID][levelID], (double) numSets);
            }
        }
      else if (lvarstd)
        {
          for (int recID = 0; recID < maxrecs; ++recID)
            {
              const auto varID = recList[recID].varID;
              const auto levelID = recList[recID].levelID;

              if (vlistInqVarTimetype(vlistID1, varID) == TIME_CONSTANT) continue;

              const auto &rsamp1 = samp1[varID][levelID];
              auto &rvars1 = vars1[varID][levelID];
              const auto &rvars2 = vars2[varID][levelID];

              if (!rsamp1.empty())
                field2_stdvar_func(rvars1, rvars2, rsamp1, divisor);
              else
                fieldc_stdvar_func(rvars1, rvars2, numSets, divisor);
            }
        }

      if (Options::cdoVerbose) cdo_print("%s  vfrac = %g, numSets = %d", datetime_to_string(vDateTimeN), vfrac, numSets);

      if (lvfrac && operfunc == FieldFunc_Mean)
        for (int recID = 0; recID < maxrecs; ++recID)
          {
            const auto varID = recList[recID].varID;
            const auto levelID = recList[recID].levelID;
            auto &rvars1 = vars1[varID][levelID];

            if (vlistInqVarTimetype(vlistID1, varID) == TIME_CONSTANT) continue;

            const auto missval = rvars1.missval;
            if (!samp1[varID][levelID].empty())
              {
                const auto fieldsize = rvars1.size;
                size_t irun = 0;
                for (size_t i = 0; i < fieldsize; ++i)
                  {
                    if ((samp1[varID][levelID].vec_d[i] / numSets) < vfrac)
                      {
                        rvars1.vec_d[i] = missval;
                        irun++;
                      }
                  }

                if (irun) rvars1.nmiss = field_num_miss(rvars1);
              }
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
          const auto varID = recList[recID].varID;
          const auto levelID = recList[recID].levelID;
          auto &rvars1 = vars1[varID][levelID];

          if (otsID && vlistInqVarTimetype(vlistID1, varID) == TIME_CONSTANT) continue;

          cdo_def_record(streamID2, varID, levelID);
          cdo_write_record(streamID2, rvars1.vec_d.data(), rvars1.nmiss);

          if (Options::cdoDiag)
            {
              if (!samp1[varID][levelID].empty())
                {
                  cdo_def_record(streamID3, varID, levelID);
                  cdo_write_record(streamID3, samp1[varID][levelID].vec_d.data(), 0);
                }
            }
        }

      if (nrecs == 0) break;
      otsID++;
    }

  if (Options::cdoDiag) cdo_stream_close(streamID3);
  cdo_stream_close(streamID2);
#ifdef USE_CDI_STREAM
  streamClose(streamID1);
#else
  cdo_stream_close(streamID1);
#endif

  if (read_task) delete read_task;

  cdo_finish();

  return nullptr;
}
