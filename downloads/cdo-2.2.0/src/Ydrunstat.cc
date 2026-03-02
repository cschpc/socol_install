/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Copyright (C) 2006 Brockmann Consult

  Author: Ralf Quast
          Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Ydrunstat    ydrunmin          Multi-year daily running minimum
      Ydrunstat    ydrunmax          Multi-year daily running maximum
      Ydrunstat    ydrunsum          Multi-year daily running sum
      Ydrunstat    ydrunmean         Multi-year daily running mean
      Ydrunstat    ydrunavg          Multi-year daily running average
      Ydrunstat    ydrunvar          Multi-year daily running variance
      Ydrunstat    ydrunvar1         Multi-year daily running variance [Normalize by (n-1)]
      Ydrunstat    ydrunstd          Multi-year daily running standard deviation
      Ydrunstat    ydrunstd1         Multi-year daily running standard deviation [Normalize by (n-1)]
*/

#include "cdi.h"
#include "calendar.h"

#include "process_int.h"
#include "param_conversion.h"
#include "datetime.h"
#include "field_functions.h"

constexpr int MaxDays = 373;

struct YdayStats
{
  CdiDateTime vDateTime[MaxDays]{};
  int numSets[MaxDays]{};
  FieldVector2D vars1[MaxDays];
  FieldVector2D vars2[MaxDays];
  int vlist;

  explicit YdayStats(int vlistID) : vlist(vlistID) {}
};

static void
ydstat_update(YdayStats &stats, CdiDateTime vDateTime, const FieldVector2D &vars1, const FieldVector2D &vars2, int numSets,
              int operfunc)
{
  auto lvarstd = (vars2.size() > 0);

  auto dayOfYear = decode_day_of_year(vDateTime.date);
  if (dayOfYear < 0 || dayOfYear >= MaxDays) cdo_abort("Day %d out of range!", dayOfYear);

  stats.vDateTime[dayOfYear] = vDateTime;

  if (!stats.vars1[dayOfYear].size())
    {
      fields_from_vlist(stats.vlist, stats.vars1[dayOfYear], FIELD_VEC);
      if (lvarstd) fields_from_vlist(stats.vlist, stats.vars2[dayOfYear], FIELD_VEC);
    }

  auto nvars = vlistNvars(stats.vlist);
  for (int varID = 0; varID < nvars; ++varID)
    {
      if (vlistInqVarTimetype(stats.vlist, varID) == TIME_CONSTANT) continue;

      auto nlevels = zaxisInqSize(vlistInqVarZaxis(stats.vlist, varID));
      for (int levelID = 0; levelID < nlevels; ++levelID)
        {
          if (stats.numSets[dayOfYear] == 0)
            {
              field_copy(vars1[varID][levelID], stats.vars1[dayOfYear][varID][levelID]);
              if (lvarstd) field_copy(vars2[varID][levelID], stats.vars2[dayOfYear][varID][levelID]);
            }
          else
            {
              if (lvarstd)
                {
                  field2_sum(stats.vars1[dayOfYear][varID][levelID], vars1[varID][levelID]);
                  field2_sum(stats.vars2[dayOfYear][varID][levelID], vars2[varID][levelID]);
                }
              else { field2_function(stats.vars1[dayOfYear][varID][levelID], vars1[varID][levelID], operfunc); }
            }
        }
    }

  stats.numSets[dayOfYear] += numSets;
}

static void
ydstat_finalize(YdayStats &stats, int operfunc)
{
  auto lmean = (operfunc == FieldFunc_Mean || operfunc == FieldFunc_Avg);
  auto lstd = (operfunc == FieldFunc_Std || operfunc == FieldFunc_Std1);
  auto lvarstd = (lstd || operfunc == FieldFunc_Var || operfunc == FieldFunc_Var1);
  const int divisor = (operfunc == FieldFunc_Std1 || operfunc == FieldFunc_Var1);

  auto fieldc_stdvar_func = lstd ? fieldc_std : fieldc_var;

  for (int dayOfYear = 0; dayOfYear < MaxDays; dayOfYear++)
    if (stats.numSets[dayOfYear])
      {
        auto nvars = vlistNvars(stats.vlist);
        for (int varID = 0; varID < nvars; ++varID)
          {
            if (vlistInqVarTimetype(stats.vlist, varID) == TIME_CONSTANT) continue;

            auto nlevels = zaxisInqSize(vlistInqVarZaxis(stats.vlist, varID));
            for (int levelID = 0; levelID < nlevels; ++levelID)
              {
                auto numSets = stats.numSets[dayOfYear];
                auto &rvars1 = stats.vars1[dayOfYear][varID][levelID];

                if (lmean) { fieldc_div(rvars1, (double) numSets); }
                else if (lvarstd)
                  {
                    const auto &rvars2 = stats.vars2[dayOfYear][varID][levelID];
                    fieldc_stdvar_func(rvars1, rvars2, numSets, divisor);
                  }
              }
          }
      }
}

static void
addOperators(void)
{
  // clang-format off
  cdo_operator_add("ydrunmin",   FieldFunc_Min,   0, nullptr);
  cdo_operator_add("ydrunmax",   FieldFunc_Max,   0, nullptr);
  cdo_operator_add("ydrunsum",   FieldFunc_Sum,   0, nullptr);
  cdo_operator_add("ydrunmean",  FieldFunc_Mean,  0, nullptr);
  cdo_operator_add("ydrunavg",   FieldFunc_Avg,   0, nullptr);
  cdo_operator_add("ydrunvar",   FieldFunc_Var,   0, nullptr);
  cdo_operator_add("ydrunvar1",  FieldFunc_Var1,  0, nullptr);
  cdo_operator_add("ydrunstd",   FieldFunc_Std,   0, nullptr);
  cdo_operator_add("ydrunstd1",  FieldFunc_Std1,  0, nullptr);
  // clang-format on
}

class ModuleYdrunstat
{
private:
  int operfunc;

  CdoStreamID streamID1;
  CdoStreamID streamID2;

  int taxisID1;
  int taxisID2;
  int vlistID1;

  char readMethod;
  bool lvarstd;
  int maxrecs;
  int ndates;
  int dpy;

  VarList varList1;
  FieldVector3D vars1;
  FieldVector3D vars2;
  std::vector<RecordInfo> recList;
  std::vector<CdiDateTime> cdiDateTimes;

public:
  void
  init(void *process)
  {
    cdo_initialize(process);

    addOperators();

    auto operatorID = cdo_operator_id();
    operfunc = cdo_operator_f1(operatorID);

    operator_input_arg("number of timesteps");
    auto nparams = cdo_operator_argc();
    ndates = parameter_to_int(cdo_operator_argv(0));
    readMethod = (nparams == 2) ? cdo_operator_argv(1)[0] : '0';

    lvarstd = (operfunc == FieldFunc_Std || operfunc == FieldFunc_Var || operfunc == FieldFunc_Std1 || operfunc == FieldFunc_Var1);
    streamID1 = cdo_open_read(0);

    vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = vlistDuplicate(vlistID1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    if (taxisHasBounds(taxisID2)) taxisDeleteBounds(taxisID2);
    vlistDefTaxis(vlistID2, taxisID2);

    dpy = calendar_dpy(taxisInqCalendar(taxisID1));

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);

    maxrecs = vlistNrecs(vlistID1);
    recList = std::vector<RecordInfo>(maxrecs);

    varListInit(varList1, vlistID1);

    cdiDateTimes = std::vector<CdiDateTime>(ndates + 1);

    vars1 = FieldVector3D(ndates + 1);
    vars2 = FieldVector3D(ndates + 1);

    for (int its = 0; its < ndates; its++)
      {
        fields_from_vlist(vlistID1, vars1[its], FIELD_VEC);
        if (lvarstd) fields_from_vlist(vlistID1, vars2[its], FIELD_VEC);
      }
  }
  void
  run()
  {
    YdayStats stats = YdayStats(vlistID1);
    int startYear = 0;
    int tsID = 0;

    for (tsID = 0; tsID < ndates; ++tsID)
      {
        auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
        if (nrecs == 0) cdo_abort("File has less then %d timesteps!", ndates);

        cdiDateTimes[tsID] = taxisInqVdatetime(taxisID1);

        if (tsID == 0 && readMethod == 'c') startYear = cdiDateTimes[tsID].date.year;

        for (int recID = 0; recID < nrecs; ++recID)
          {
            int varID, levelID;
            cdo_inq_record(streamID1, &varID, &levelID);

            if (tsID == 0) recList[recID].set(varID, levelID);

            auto &rvars1 = vars1[tsID][varID][levelID];

            cdo_read_record(streamID1, rvars1);

            if (lvarstd)
              {
                field2_moq(vars2[tsID][varID][levelID], rvars1);
                for (int inp = 0; inp < tsID; ++inp) field2_sumsumq(vars1[inp][varID][levelID], vars2[inp][varID][levelID], rvars1);
              }
            else
              {
                for (int inp = 0; inp < tsID; ++inp) field2_function(vars1[inp][varID][levelID], rvars1, operfunc);
              }
          }
      }

    while (true)
      {
        cdiDateTimes[ndates] = datetime_avg(dpy, ndates, cdiDateTimes);

        ydstat_update(stats, cdiDateTimes[ndates], vars1[0], vars2[0], ndates, operfunc);

        cdiDateTimes[ndates] = cdiDateTimes[0];
        vars1[ndates] = vars1[0];
        if (lvarstd) vars2[ndates] = vars2[0];

        for (int inp = 0; inp < ndates; ++inp)
          {
            cdiDateTimes[inp] = cdiDateTimes[inp + 1];
            vars1[inp] = vars1[inp + 1];
            if (lvarstd) vars2[inp] = vars2[inp + 1];
          }

        auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
        if (nrecs == 0) break;

        cdiDateTimes[ndates - 1] = taxisInqVdatetime(taxisID1);

        for (int recID = 0; recID < nrecs; ++recID)
          {
            int varID, levelID;
            cdo_inq_record(streamID1, &varID, &levelID);

            auto &rvars1 = vars1[ndates - 1][varID][levelID];

            cdo_read_record(streamID1, rvars1);

            if (lvarstd)
              {
                field2_moq(vars2[ndates - 1][varID][levelID], rvars1);
                for (int inp = 0; inp < ndates - 1; ++inp)
                  field2_sumsumq(vars1[inp][varID][levelID], vars2[inp][varID][levelID], rvars1);
              }
            else
              {
                for (int inp = 0; inp < ndates - 1; ++inp) field2_function(vars1[inp][varID][levelID], rvars1, operfunc);
              }
          }

        tsID++;
      }

    if (readMethod == 'c' && cdo_assert_files_only())
      {
        auto endYear = cdiDateTimes[ndates - 1].date.year;
        auto cdiStream = streamOpenRead(cdo_get_stream_name(0));
        auto cdiVlistID = streamInqVlist(cdiStream);
        auto cdiTaxisID = vlistInqTaxis(cdiVlistID);
        int missTimes = 0;
        for (missTimes = 0; missTimes < ndates - 1; missTimes++)
          {
            auto nrecs = streamInqTimestep(cdiStream, missTimes);
            if (nrecs == 0) break;

            cdiDateTimes[ndates - 1] = taxisInqVdatetime(cdiTaxisID);
            cdiDateTimes[ndates - 1].date.year = endYear + 1;

            for (int recID = 0; recID < nrecs; ++recID)
              {
                int varID, levelID;
                streamInqRecord(cdiStream, &varID, &levelID);

                auto &rvars1 = vars1[ndates - 1][varID][levelID];

                streamReadRecord(cdiStream, rvars1.vec_d.data(), &rvars1.nmiss);

                if (lvarstd)
                  {
                    field2_moq(vars2[ndates - 1][varID][levelID], rvars1);
                    for (int inp = 0; inp < ndates - 1; ++inp)
                      field2_sumsumq(vars1[inp][varID][levelID], vars2[inp][varID][levelID], rvars1);
                  }
                else
                  {
                    for (int inp = 0; inp < ndates - 1; ++inp) field2_function(vars1[inp][varID][levelID], rvars1, operfunc);
                  }
              }

            cdiDateTimes[ndates] = datetime_avg(dpy, ndates, cdiDateTimes);
            auto vDateTime = cdiDateTimes[ndates];
            if (vDateTime.date.year > endYear) vDateTime.date.year = startYear;

            ydstat_update(stats, vDateTime, vars1[0], vars2[0], ndates, operfunc);

            cdiDateTimes[ndates] = cdiDateTimes[0];
            vars1[ndates] = vars1[0];
            if (lvarstd) vars2[ndates] = vars2[0];

            for (int inp = 0; inp < ndates; ++inp)
              {
                cdiDateTimes[inp] = cdiDateTimes[inp + 1];
                vars1[inp] = vars1[inp + 1];
                if (lvarstd) vars2[inp] = vars2[inp + 1];
              }
          }

        if (missTimes != ndates - 1) cdo_abort("Addding the missing values when using the 'readMethod' method was not possible");

        streamClose(cdiStream);
      }
    else if (readMethod == 'c')
      cdo_warning("Operators cannot be piped in circular mode");

    ydstat_finalize(stats, operfunc);

    int otsID = 0;

    for (int dayOfYear = 0; dayOfYear < MaxDays; dayOfYear++)
      if (stats.numSets[dayOfYear])
        {
          taxisDefVdatetime(taxisID2, stats.vDateTime[dayOfYear]);
          cdo_def_timestep(streamID2, otsID);

          for (int recID = 0; recID < maxrecs; ++recID)
            {
              auto [varID, levelID] = recList[recID].get();
              if (otsID && varList1[varID].isConstant) continue;

              auto &rvars1 = stats.vars1[dayOfYear][varID][levelID];

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
Ydrunstat(void *process)
{
  ModuleYdrunstat ydrunstat;
  ydrunstat.init(process);
  ydrunstat.run();
  ydrunstat.close();

  return nullptr;
}
