/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Ensstat    ensrange        Ensemble range
      Ensstat    ensmin          Ensemble minimum
      Ensstat    ensmax          Ensemble maximum
      Ensstat    enssum          Ensemble sum
      Ensstat    ensmean         Ensemble mean
      Ensstat    ensavg          Ensemble average
      Ensstat    ensstd          Ensemble standard deviation
      Ensstat    ensstd1         Ensemble standard deviation
      Ensstat    ensvar          Ensemble variance
      Ensstat    ensvar1         Ensemble variance
      Ensstat    enspctl         Ensemble percentiles
*/

#include <atomic>

#include <cdi.h>

#include "cdo_rlimit.h"
#include "process_int.h"
#include "cdo_vlist.h"
#include "param_conversion.h"
#include "cdo_task.h"
#include "cdo_options.h"
#include "util_files.h"
#include "cimdOmp.h"
#include "field_functions.h"

struct EnsstatFile
{
  VarList varList;
  CdoStreamID streamID;
  int vlistID;
  double missval[2];
  Field field[2];
};

struct EnsstatArg
{
  int t;
  int varID[2];
  int levelID[2];
  CdoStreamID streamID2;
  int nfiles;
  EnsstatFile *efData;
  double *array2Data;
  double *count2Data;
  FieldVector fields;
  int operfunc;
  double pn;
  bool lpctl;
  bool withCountData;
  int nvars;
};

static void *
ensstat_func(void *ensarg)
{
  if (Options::CDO_task) cdo_omp_set_num_threads(Threading::ompNumThreads);

  auto arg = (EnsstatArg *) ensarg;
  auto t = arg->t;
  auto nfiles = arg->nfiles;
  auto ef = arg->efData;
  auto &fields = arg->fields;
  auto array2 = arg->array2Data;
  auto count2 = arg->count2Data;
  auto withCountData = arg->withCountData;

  auto hasMissvals = false;
  for (int fileID = 0; fileID < nfiles; ++fileID)
    if (ef[fileID].field[t].nmiss > 0) hasMissvals = true;

  auto varID = arg->varID[t];
  auto gridsize = ef[0].varList[varID].gridsize;
  auto missval = ef[0].varList[varID].missval;
  auto memType = ef[0].field[t].memType;

  std::atomic<size_t> atomicNumMiss{ 0 };
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
  for (size_t i = 0; i < gridsize; ++i)
    {
      auto ompthID = cdo_omp_get_thread_num();

      auto &field = fields[ompthID];
      field.missval = missval;
      field.nmiss = 0;
      if (memType == MemType::Float)
        for (int fileID = 0; fileID < nfiles; ++fileID) field.vec_d[fileID] = ef[fileID].field[t].vec_f[i];
      else
        for (int fileID = 0; fileID < nfiles; ++fileID) field.vec_d[fileID] = ef[fileID].field[t].vec_d[i];

      if (hasMissvals)
        for (int fileID = 0; fileID < nfiles; ++fileID)
          {
            if (dbl_is_equal(field.vec_d[fileID], ef[fileID].missval[t]))
              {
                field.vec_d[fileID] = missval;
                field.nmiss++;
              }
          }

      array2[i] = arg->lpctl ? field_pctl(field, arg->pn) : field_function(field, arg->operfunc);

      if (dbl_is_equal(array2[i], field.missval)) atomicNumMiss++;

      if (withCountData) count2[i] = nfiles - field.nmiss;
    }

  size_t nmiss = atomicNumMiss;

  cdo_def_record(arg->streamID2, arg->varID[t], arg->levelID[t]);
  cdo_write_record(arg->streamID2, array2, nmiss);

  if (withCountData)
    {
      cdo_def_record(arg->streamID2, arg->varID[t] + arg->nvars, arg->levelID[t]);
      cdo_write_record(arg->streamID2, count2, 0);
    }

  return nullptr;
}

static void
addOperators(void)
{
  // clang-format off
  cdo_operator_add("ensrange",  FieldFunc_Range,  0, nullptr);
  cdo_operator_add("ensmin",    FieldFunc_Min,    0, nullptr);
  cdo_operator_add("ensmax",    FieldFunc_Max,    0, nullptr);
  cdo_operator_add("enssum",    FieldFunc_Sum,    0, nullptr);
  cdo_operator_add("ensmean",   FieldFunc_Mean,   0, nullptr);
  cdo_operator_add("ensavg",    FieldFunc_Avg,    0, nullptr);
  cdo_operator_add("ensstd",    FieldFunc_Std,    0, nullptr);
  cdo_operator_add("ensstd1",   FieldFunc_Std1,   0, nullptr);
  cdo_operator_add("ensvar",    FieldFunc_Var,    0, nullptr);
  cdo_operator_add("ensvar1",   FieldFunc_Var1,   0, nullptr);
  cdo_operator_add("ensskew",   FieldFunc_Skew,   0, nullptr);
  cdo_operator_add("enskurt",   FieldFunc_Kurt,   0, nullptr);
  cdo_operator_add("ensmedian", FieldFunc_Median, 0, nullptr);
  cdo_operator_add("enspctl",   FieldFunc_Pctl,   0, nullptr);
  // clang-format on
}

void *
Ensstat(void *process)
{
  cdo::Task *task = Options::CDO_task ? new cdo::Task : nullptr;
  int nrecs0;

  cdo_initialize(process);

  addOperators();

  auto operatorID = cdo_operator_id();
  auto operfunc = cdo_operator_f1(operatorID);

  auto lpctl = (operfunc == FieldFunc_Pctl);

  auto argc = cdo_operator_argc();
  auto nargc = argc;

  double pn = 0.0;
  if (operfunc == FieldFunc_Pctl)
    {
      operator_input_arg("percentile number");
      pn = parameter_to_double(cdo_operator_argv(0));
      argc--;
    }

  auto withCountData = false;
  if (argc == 1)
    {
      if (cdo_operator_argv(nargc - 1) == "count")
        withCountData = true;
      else
        cdo_abort("Unknown parameter: >%s<", cdo_operator_argv(nargc - 1));
    }

  auto nfiles = cdo_stream_cnt() - 1;

  if (Options::cdoVerbose) cdo_print("Ensemble over %d files.", nfiles);

  cdo::set_numfiles(nfiles + 8);

  auto ofilename = cdo_get_stream_name(nfiles);

  if (!Options::cdoOverwriteMode && FileUtils::file_exists(ofilename) && !FileUtils::user_file_overwrite(ofilename))
    cdo_abort("Outputfile %s already exists!", ofilename);

  std::vector<EnsstatFile> ef(nfiles);

  FieldVector fields(Threading::ompNumThreads);
  for (int i = 0; i < Threading::ompNumThreads; ++i) fields[i].resize(nfiles);

  for (int fileID = 0; fileID < nfiles; ++fileID)
    {
      ef[fileID].streamID = cdo_open_read(fileID);
      ef[fileID].vlistID = cdo_stream_inq_vlist(ef[fileID].streamID);
      varListInit(ef[fileID].varList, ef[fileID].vlistID);
    }

  // check that the contents is always the same
  for (int fileID = 1; fileID < nfiles; ++fileID) vlist_compare(ef[0].vlistID, ef[fileID].vlistID, CMP_ALL);

  auto vlistID1 = ef[0].vlistID;
  auto vlistID2 = vlistDuplicate(vlistID1);
  auto taxisID1 = vlistInqTaxis(vlistID1);
  auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  auto gridsizemax = vlistGridsizeMax(vlistID1);
  Varray<double> array2(gridsizemax);

  auto nvars = vlistNvars(vlistID2);
  Varray<double> count2;
  if (withCountData)
    {
      count2.resize(gridsizemax);
      for (int varID = 0; varID < nvars; ++varID)
        {
          auto name = ef[0].varList[varID].name + "_count";
          auto gridID = ef[0].varList[varID].gridID;
          auto zaxisID = ef[0].varList[varID].zaxisID;
          auto timetype = ef[0].varList[varID].timetype;
          auto cvarID = vlistDefVar(vlistID2, gridID, zaxisID, timetype);
          cdiDefKeyString(vlistID2, cvarID, CDI_KEY_NAME, name.c_str());
          vlistDefVarDatatype(vlistID2, cvarID, CDI_DATATYPE_INT16);
          if (cvarID != (varID + nvars)) cdo_abort("Internal error, varIDs do not match!");
        }
    }

  auto streamID2 = cdo_open_write(nfiles);
  cdo_def_vlist(streamID2, vlistID2);

  EnsstatArg ensstatArg;
  ensstatArg.streamID2 = streamID2;
  ensstatArg.nfiles = nfiles;
  ensstatArg.array2Data = array2.data();
  ensstatArg.count2Data = count2.data();
  ensstatArg.fields.resize(Threading::ompNumThreads);
  for (int i = 0; i < Threading::ompNumThreads; ++i) ensstatArg.fields[i].resize(nfiles);
  ensstatArg.operfunc = operfunc;
  ensstatArg.pn = pn;
  ensstatArg.lpctl = lpctl;
  ensstatArg.withCountData = withCountData;
  ensstatArg.nvars = nvars;
  ensstatArg.t = 0;

  auto printWarning = false;
  auto printError = false;
  int t = 0;
  int tsID = 0;
  do {
      nrecs0 = cdo_stream_inq_timestep(ef[0].streamID, tsID);
      for (int fileID = 1; fileID < nfiles; ++fileID)
        {
          auto streamID = ef[fileID].streamID;
          auto nrecs = cdo_stream_inq_timestep(streamID, tsID);
          if (nrecs != nrecs0)
            {
              if (nrecs == 0)
                {
                  printWarning = true;
                  cdo_warning("Inconsistent ensemble file, too few time steps in %s!", cdo_get_stream_name(fileID));
                }
              else if (nrecs0 == 0)
                {
                  printWarning = true;
                  cdo_warning("Inconsistent ensemble file, too few time steps in %s!", cdo_get_stream_name(0));
                }
              else
                {
                  printError = true;
                  cdo_warning("Inconsistent ensemble file, number of records at time step %d of %s and %s differ!", tsID + 1,
                              cdo_get_stream_name(0), cdo_get_stream_name(fileID));
                }
              goto CLEANUP;
            }
        }

      if (nrecs0 > 0)
        {
          cdo_taxis_copy_timestep(taxisID2, taxisID1);
          cdo_def_timestep(streamID2, tsID);
        }

      for (int recID = 0; recID < nrecs0; ++recID)
        {
          int varID = -1, levelID = -1;

          for (int fileID = 0; fileID < nfiles; ++fileID)
            {
              cdo_inq_record(ef[fileID].streamID, &varID, &levelID);
              ef[fileID].missval[t] = ef[fileID].varList[varID].missval;
              ef[fileID].field[t].init(ef[fileID].varList[varID]);
            }
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
          for (int fileID = 0; fileID < nfiles; ++fileID) { cdo_read_record(ef[fileID].streamID, ef[fileID].field[t]); }

          ensstatArg.efData = ef.data();
          ensstatArg.varID[t] = varID;
          ensstatArg.levelID[t] = levelID;
          if (Options::CDO_task)
            {
              task->start(ensstat_func, &ensstatArg);
              task->wait();
              // t = !t;
            }
          else { ensstat_func(&ensstatArg); }
        }

      tsID++;
    }
  while (nrecs0 > 0);

CLEANUP:

  if (printWarning) cdo_warning("Inconsistent ensemble, processed only the first %d timesteps!", tsID);
  if (printError) cdo_abort("Inconsistent ensemble, processed only the first %d timesteps!", tsID);

  for (int fileID = 0; fileID < nfiles; ++fileID) cdo_stream_close(ef[fileID].streamID);

  cdo_stream_close(streamID2);

  if (task) delete task;

  cdo_finish();

  return nullptr;
}
