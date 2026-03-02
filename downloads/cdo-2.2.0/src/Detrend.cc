/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Detrend    detrend         Detrend
*/

#include <cdi.h>

#include "varray.h"
#include "process_int.h"
#include "cdo_vlist.h"
#include "cdo_options.h"
#include "datetime.h"
#include "cimdOmp.h"
#include "pmlist.h"
#include "param_conversion.h"
#include "field_functions.h"

static void
detrend(const long nts, const Varray<double> &deltaTS0, const double missval1, const Varray<double> &array1, Varray<double> &array2)
{
  auto missval2 = missval1;
  double sumj = 0.0, sumjj = 0.0;
  double sumx = 0.0, sumjx = 0.0;
  long n = 0;
  for (long j = 0; j < nts; ++j)
    if (!dbl_is_equal(array1[j], missval1))
      {
        auto zj = deltaTS0[j];
        sumj += zj;
        sumjj += zj * zj;
        sumjx += zj * array1[j];
        sumx += array1[j];
        n++;
      }

  auto work1 = DIVMN(SUBMN(sumjx, DIVMN(MULMN(sumj, sumx), n)), SUBMN(sumjj, DIVMN(MULMN(sumj, sumj), n)));
  auto work2 = SUBMN(DIVMN(sumx, n), MULMN(DIVMN(sumj, n), work1));

  for (long j = 0; j < nts; ++j) array2[j] = SUBMN(array1[j], ADDMN(work2, MULMN(work1, deltaTS0[j])));
}

static void
computeDeltaTS0(bool tstepIsEqual, int nts, int calendar, DateTimeList &dtlist, Varray<double> &deltaTS0)
{
  CheckTimeIncr checkTimeIncr;
  JulianDate julianDate0;
  double deltat1 = 0.0;

  for (int tsID = 0; tsID < nts; ++tsID)
    {
      auto vDateTime = dtlist.get_vDateTime(tsID);
      if (tstepIsEqual) check_time_increment(tsID, calendar, vDateTime, checkTimeIncr);
      deltaTS0[tsID] = tstepIsEqual ? (double) tsID : delta_time_step_0(tsID, calendar, vDateTime, julianDate0, deltat1);
    }
}

static void
detrendGetParameter(bool &tstepIsEqual)
{
  auto pargc = cdo_operator_argc();
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

          // clang-format off
          if      (key == "equal") tstepIsEqual = parameter_to_bool(value);
          else cdo_abort("Invalid parameter key >%s<!", key);
          // clang-format on
        }
    }
}

void *
Detrend(void *process)
{
  int varID, levelID;
  DateTimeList dtlist;

  cdo_initialize(process);

  auto tstepIsEqual = true;
  detrendGetParameter(tstepIsEqual);

  auto streamID1 = cdo_open_read(0);

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto vlistID2 = vlistDuplicate(vlistID1);

  auto taxisID1 = vlistInqTaxis(vlistID1);
  auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  auto streamID2 = cdo_open_write(1);
  cdo_def_vlist(streamID2, vlistID2);

  VarList varList;
  varListInit(varList, vlistID1);

  auto nvars = vlistNvars(vlistID1);
  FieldVector3D vars;

  int tsID = 0;
  while (true)
    {
      auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
      if (nrecs == 0) break;

      constexpr size_t NALLOC_INC = 1024;
      if ((size_t) tsID >= vars.size()) vars.resize(vars.size() + NALLOC_INC);

      dtlist.taxis_inq_timestep(taxisID1, tsID);

      fields_from_vlist(vlistID1, vars[tsID]);

      for (int recID = 0; recID < nrecs; ++recID)
        {
          cdo_inq_record(streamID1, &varID, &levelID);
          auto &field = vars[tsID][varID][levelID];
          field.init(varList[varID]);
          cdo_read_record(streamID1, field);
        }

      tsID++;
    }

  auto nts = tsID;
  Varray<double> deltaTS0(nts);
  Varray2D<double> array1_2D(Threading::ompNumThreads, Varray<double>(nts));
  Varray2D<double> array2_2D(Threading::ompNumThreads, Varray<double>(nts));

  auto calendar = taxisInqCalendar(taxisID1);
  computeDeltaTS0(tstepIsEqual, nts, calendar, dtlist, deltaTS0);

  for (varID = 0; varID < nvars; ++varID)
    {
      const auto &var = varList[varID];
      auto nsteps = var.isConstant ? 1 : nts;
      auto missval = var.missval;
      auto fieldMemType = var.memType;
      auto gridsize = var.gridsize;
      for (levelID = 0; levelID < var.nlevels; ++levelID)
        {
#ifdef _OPENMP
#pragma omp parallel for default(none) schedule(static) \
    shared(fieldMemType, gridsize, nsteps, deltaTS0, missval, array1_2D, array2_2D, vars, varID, levelID)
#endif
          for (size_t i = 0; i < gridsize; ++i)
            {
              auto ompthID = cdo_omp_get_thread_num();
              auto &array1 = array1_2D[ompthID];
              auto &array2 = array2_2D[ompthID];

              if (fieldMemType == MemType::Float)
                for (int k = 0; k < nsteps; ++k) array1[k] = vars[k][varID][levelID].vec_f[i];
              else
                for (int k = 0; k < nsteps; ++k) array1[k] = vars[k][varID][levelID].vec_d[i];

              detrend(nsteps, deltaTS0, missval, array1, array2);

              if (fieldMemType == MemType::Float)
                for (int k = 0; k < nsteps; ++k) vars[k][varID][levelID].vec_f[i] = array2[k];
              else
                for (int k = 0; k < nsteps; ++k) vars[k][varID][levelID].vec_d[i] = array2[k];
            }
        }
    }

  for (tsID = 0; tsID < nts; ++tsID)
    {
      dtlist.taxis_def_timestep(taxisID2, tsID);
      cdo_def_timestep(streamID2, tsID);

      for (varID = 0; varID < nvars; ++varID)
        {
          const auto &var = varList[varID];
          if (tsID && var.isConstant) continue;
          for (levelID = 0; levelID < var.nlevels; ++levelID)
            {
              cdo_def_record(streamID2, varID, levelID);
              auto &field = vars[tsID][varID][levelID];
              cdo_write_record(streamID2, field);
            }
        }
    }

  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  cdo_finish();

  return nullptr;
}
