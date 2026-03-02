/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Runstat    runrange        Running range
      Runstat    runmin          Running minimum
      Runstat    runmax          Running maximum
      Runstat    runsum          Running sum
      Runstat    runmean         Running mean
      Runstat    runavg          Running average
      Runstat    runvar          Running variance
      Runstat    runvar1         Running variance [Normalize by (n-1)]
      Runstat    runstd          Running standard deviation
      Runstat    runstd1         Running standard deviation [Normalize by (n-1)]
*/

#include <cdi.h>

#include "process_int.h"
#include "param_conversion.h"
#include "datetime.h"
#include "field_functions.h"

static void
addOperators(void)
{
  // clang-format off
  cdo_operator_add("runrange", FieldFunc_Range, 0, nullptr);
  cdo_operator_add("runmin",   FieldFunc_Min,   0, nullptr);
  cdo_operator_add("runmax",   FieldFunc_Max,   0, nullptr);
  cdo_operator_add("runsum",   FieldFunc_Sum,   0, nullptr);
  cdo_operator_add("runmean",  FieldFunc_Mean,  0, nullptr);
  cdo_operator_add("runavg",   FieldFunc_Avg,   0, nullptr);
  cdo_operator_add("runvar",   FieldFunc_Var,   0, nullptr);
  cdo_operator_add("runvar1",  FieldFunc_Var1,  0, nullptr);
  cdo_operator_add("runstd",   FieldFunc_Std,   0, nullptr);
  cdo_operator_add("runstd1",  FieldFunc_Std1,  0, nullptr);
  // clang-format on
}

void *
Runstat(void *process)
{
  const TimeStat timestat_date = TimeStat::MEAN;
  bool runstat_nomiss = false;

  cdo_initialize(process);

  auto envstr = getenv("RUNSTAT_NOMISS");
  if (envstr)
    {
      char *endptr;
      const auto envval = (int) strtol(envstr, &endptr, 10);
      if (envval == 1) runstat_nomiss = true;
    }

  addOperators();

  const auto operatorID = cdo_operator_id();
  auto operfunc = cdo_operator_f1(operatorID);  // used in omp loop

  const auto lrange = (operfunc == FieldFunc_Range);
  const auto lmean = (operfunc == FieldFunc_Mean || operfunc == FieldFunc_Avg);
  const auto lstd = (operfunc == FieldFunc_Std || operfunc == FieldFunc_Std1);
  const auto lvarstd = (lstd || operfunc == FieldFunc_Var || operfunc == FieldFunc_Var1);
  const auto lvars2 = (lvarstd || lrange);
  const int divisor = (operfunc == FieldFunc_Std1 || operfunc == FieldFunc_Var1);

  auto field2_stdvar_func = lstd ? field2_std : field2_var;
  auto fieldc_stdvar_func = lstd ? fieldc_std : fieldc_var;

  operator_input_arg("number of timesteps");
  operator_check_argc(1);
  auto ndates = parameter_to_int(cdo_operator_argv(0));

  const auto streamID1 = cdo_open_read(0);

  const auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  const auto vlistID2 = vlistDuplicate(vlistID1);

  const auto taxisID1 = vlistInqTaxis(vlistID1);
  const auto taxisID2 = taxisDuplicate(taxisID1);
  taxisWithBounds(taxisID2);
  vlistDefTaxis(vlistID2, taxisID2);
  // Number of timestep will be reduced compared to the input error handling in case of not enough timesteps is done per record
  auto nsteps = vlistNtsteps(vlistID1);
  if (nsteps != -1)
    {
      nsteps -= ndates - 1;
      if (nsteps > 0) vlistDefNtsteps(vlistID2, nsteps);
    }

  const auto streamID2 = cdo_open_write(1);
  cdo_def_vlist(streamID2, vlistID2);

  const auto maxrecs = vlistNrecs(vlistID1);
  std::vector<RecordInfo> recList(maxrecs);

  DateTimeList dtlist;
  dtlist.set_stat(timestat_date);
  dtlist.set_calendar(taxisInqCalendar(taxisID1));

  VarList varList;
  varListInit(varList, vlistID1);

  int VARS_MEMTYPE = 0;
  if ((operfunc == FieldFunc_Min) || (operfunc == FieldFunc_Max)) VARS_MEMTYPE = FIELD_NAT;

  FieldVector3D vars1(ndates + 1), vars2, samp1;
  if (!runstat_nomiss) samp1.resize(ndates + 1);
  if (lvars2) vars2.resize(ndates + 1);

  for (int its = 0; its < ndates; its++)
    {
      if (!runstat_nomiss) fields_from_vlist(vlistID1, samp1[its], FIELD_VEC);
      fields_from_vlist(vlistID1, vars1[its], FIELD_VEC | VARS_MEMTYPE);
      if (lvars2) fields_from_vlist(vlistID1, vars2[its], FIELD_VEC);
    }

  const auto gridsizemax = vlistGridsizeMax(vlistID1);
  std::vector<bool> imask(gridsizemax);

  int tsID = 0;
  int otsID = 0;
  int numSteps = 0;
  while (true)
    {
    FILL_FIRST_NDATES:
      const auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
      if (nrecs == 0)
        {
          if (tsID < ndates)
            cdo_abort("File has less then %d timesteps!", ndates);
          else
            break;
        }

      numSteps = (tsID < ndates) ? tsID : ndates - 1;

      dtlist.taxis_inq_timestep(taxisID1, numSteps);

      for (int recID = 0; recID < nrecs; ++recID)
        {
          int varID, levelID;
          cdo_inq_record(streamID1, &varID, &levelID);

          if (tsID == 0) recList[recID].set(varID, levelID);

          auto &rvars1 = vars1[numSteps][varID][levelID];

          auto fieldsize = rvars1.size;  // used in omp loop

          cdo_read_record(streamID1, rvars1);
          if (lrange)
            {
              vars2[numSteps][varID][levelID].nmiss = rvars1.nmiss;
              vars2[numSteps][varID][levelID].vec_d = rvars1.vec_d;
            }

          if (runstat_nomiss && rvars1.nmiss) cdo_abort("Missing values supported was swichted off by env. RUNSTAT_NOMISS!");

          if (!runstat_nomiss)
            {
              const auto missval = rvars1.missval;

              if (rvars1.memType == MemType::Float)
                for (size_t i = 0; i < fieldsize; ++i) imask[i] = !DBL_IS_EQUAL(rvars1.vec_f[i], missval);
              else
                for (size_t i = 0; i < fieldsize; ++i) imask[i] = !DBL_IS_EQUAL(rvars1.vec_d[i], missval);

              for (size_t i = 0; i < fieldsize; ++i) samp1[numSteps][varID][levelID].vec_d[i] = (double) imask[i];

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
              for (int inp = 0; inp < numSteps; ++inp)
                {
                  auto &samp = samp1[inp][varID][levelID].vec_d;
                  for (size_t i = 0; i < fieldsize; ++i)
                    if (imask[i]) samp[i]++;
                }
            }

          if (lvarstd)
            {
              field2_moq(vars2[numSteps][varID][levelID], vars1[numSteps][varID][levelID]);
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
              for (int inp = 0; inp < numSteps; ++inp)
                {
                  field2_sumsumq(vars1[inp][varID][levelID], vars2[inp][varID][levelID], rvars1);
                }
            }
          else if (lrange)
            {
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
              for (int inp = 0; inp < numSteps; ++inp)
                {
                  field2_maxmin(vars1[inp][varID][levelID], vars2[inp][varID][levelID], rvars1);
                }
            }
          else
            {
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
              for (int inp = 0; inp < numSteps; ++inp) { field2_function(vars1[inp][varID][levelID], rvars1, operfunc); }
            }
        }

      tsID++;  // don't move this line

      if (tsID < ndates) goto FILL_FIRST_NDATES;

      for (int recID = 0; recID < maxrecs; ++recID)
        {
          auto [varID, levelID] = recList[recID].get();
          if (varList[varID].isConstant) continue;

          const auto &rsamp1 = samp1[0][varID][levelID];
          auto &rvars1 = vars1[0][varID][levelID];
          const auto numSets = ndates;

          if (lmean)
            {
              if (!runstat_nomiss)
                field2_div(rvars1, rsamp1);
              else
                fieldc_div(rvars1, (double) numSets);
            }
          else if (lvarstd)
            {
              if (!runstat_nomiss)
                field2_stdvar_func(rvars1, vars2[0][varID][levelID], rsamp1, divisor);
              else
                fieldc_stdvar_func(rvars1, vars2[0][varID][levelID], numSets, divisor);
            }
          else if (lrange) { field2_sub(rvars1, vars2[0][varID][levelID]); }
        }

      dtlist.stat_taxis_def_timestep(taxisID2, ndates);
      cdo_def_timestep(streamID2, otsID);

      for (int recID = 0; recID < maxrecs; ++recID)
        {
          auto [varID, levelID] = recList[recID].get();
          if (otsID && varList[varID].isConstant) continue;

          auto &rvars1 = vars1[0][varID][levelID];

          cdo_def_record(streamID2, varID, levelID);
          cdo_write_record(streamID2, rvars1);
        }

      otsID++;

      dtlist.shift();

      vars1[ndates] = vars1[0];
      if (!runstat_nomiss) samp1[ndates] = samp1[0];
      if (lvars2) vars2[ndates] = vars2[0];

      for (int inp = 0; inp < ndates; ++inp)
        {
          vars1[inp] = vars1[inp + 1];
          if (!runstat_nomiss) samp1[inp] = samp1[inp + 1];
          if (lvars2) vars2[inp] = vars2[inp + 1];
        }
    }

  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  cdo_finish();

  return nullptr;
}
