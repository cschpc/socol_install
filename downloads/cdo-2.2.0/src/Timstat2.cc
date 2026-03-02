/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

        Timstat2        timcor      correlates two data files on the same grid
*/

#include <cdi.h>

#include "process_int.h"
#include "cdo_vlist.h"
#include "cimdOmp.h"
#include "field_functions.h"

// correlation in time
template <typename T1, typename T2>
static void
correlation_init(bool hasMissValues, size_t gridsize, const Varray<T1> &x, const Varray<T2> &y, T1 xmv, T2 ymv,
                 Varray<size_t> &nofvals, Varray<double> &work0, Varray<double> &work1, Varray<double> &work2,
                 Varray<double> &work3, Varray<double> &work4)
{
  if (hasMissValues)
    {
#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
      for (size_t i = 0; i < gridsize; ++i)
        {
          if ((!DBL_IS_EQUAL(x[i], xmv)) && (!DBL_IS_EQUAL(y[i], ymv)))
            {
              const double xx = x[i];
              const double yy = y[i];
              work0[i] += xx;
              work1[i] += yy;
              work2[i] += xx * xx;
              work3[i] += yy * yy;
              work4[i] += xx * yy;
              nofvals[i]++;
            }
        }
    }
  else
    {
#if _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
      for (size_t i = 0; i < gridsize; ++i)
        {
          const double xx = x[i];
          const double yy = y[i];
          work0[i] += xx;
          work1[i] += yy;
          work2[i] += xx * xx;
          work3[i] += yy * yy;
          work4[i] += xx * yy;
          nofvals[i]++;
        }
    }
}

static void
correlation_init(size_t gridsize, const Field &field1, const Field &field2, Varray<size_t> &nofvals, Varray<double> &work0,
                 Varray<double> &work1, Varray<double> &work2, Varray<double> &work3, Varray<double> &work4)
{
  auto hasMissValues = (field1.nmiss > 0 || field2.nmiss > 0);

  if (field1.memType == MemType::Float && field2.memType == MemType::Float)
    correlation_init(hasMissValues, gridsize, field1.vec_f, field2.vec_f, (float) field1.missval, (float) field2.missval, nofvals,
                     work0, work1, work2, work3, work4);
  else if (field1.memType == MemType::Float && field2.memType == MemType::Double)
    correlation_init(hasMissValues, gridsize, field1.vec_f, field2.vec_d, (float) field1.missval, field2.missval, nofvals, work0,
                     work1, work2, work3, work4);
  else if (field1.memType == MemType::Double && field2.memType == MemType::Float)
    correlation_init(hasMissValues, gridsize, field1.vec_d, field2.vec_f, field1.missval, (float) field2.missval, nofvals, work0,
                     work1, work2, work3, work4);
  else
    correlation_init(hasMissValues, gridsize, field1.vec_d, field2.vec_d, field1.missval, field2.missval, nofvals, work0, work1,
                     work2, work3, work4);
}

static size_t
correlation(size_t gridsize, double missval, const Varray<size_t> &nofvals, Varray<double> &work0, const Varray<double> &work1,
            const Varray<double> &work2, const Varray<double> &work3, const Varray<double> &work4)
{
  size_t nmiss = 0;

  for (size_t i = 0; i < gridsize; ++i)
    {
      auto missval1 = missval;
      auto missval2 = missval;
      double cor;
      auto nvals = nofvals[i];
      if (nvals > 0)
        {
          auto temp0 = MULMN(work0[i], work1[i]);
          auto temp1 = SUBMN(work4[i], DIVMN(temp0, nvals));
          auto temp2 = MULMN(work0[i], work0[i]);
          auto temp3 = MULMN(work1[i], work1[i]);
          auto temp4 = SUBMN(work2[i], DIVMN(temp2, nvals));
          auto temp5 = SUBMN(work3[i], DIVMN(temp3, nvals));
          auto temp6 = MULMN(temp4, temp5);

          cor = DIVMN(temp1, SQRTMN(temp6));
          cor = std::min(std::max(cor, -1.0), 1.0);

          if (DBL_IS_EQUAL(cor, missval)) nmiss++;
        }
      else
        {
          nmiss++;
          cor = missval;
        }

      work0[i] = cor;
    }

  return nmiss;
}

// covariance in time
template <typename T1, typename T2>
static void
covariance_init(bool hasMissValues, size_t gridsize, const Varray<T1> &x, const Varray<T2> &y, T1 xmv, T2 ymv,
                Varray<size_t> &nofvals, Varray<double> &work0, Varray<double> &work1, Varray<double> &work2)
{
  if (hasMissValues)
    {
#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
      for (size_t i = 0; i < gridsize; ++i)
        {
          if ((!DBL_IS_EQUAL(x[i], xmv)) && (!DBL_IS_EQUAL(y[i], ymv)))
            {
              const double xx = x[i];
              const double yy = y[i];
              work0[i] += xx;
              work1[i] += yy;
              work2[i] += xx * yy;
              nofvals[i]++;
            }
        }
    }
  else
    {
#if _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
      for (size_t i = 0; i < gridsize; ++i)
        {
          const double xx = x[i];
          const double yy = y[i];
          work0[i] += xx;
          work1[i] += yy;
          work2[i] += xx * yy;
          nofvals[i]++;
        }
    }
}

static void
covariance_init(size_t gridsize, const Field &field1, const Field &field2, Varray<size_t> &nofvals, Varray<double> &work0,
                Varray<double> &work1, Varray<double> &work2)
{
  auto hasMissValues = (field1.nmiss > 0 || field2.nmiss > 0);

  if (field1.memType == MemType::Float && field2.memType == MemType::Float)
    covariance_init(hasMissValues, gridsize, field1.vec_f, field2.vec_f, (float) field1.missval, (float) field2.missval, nofvals,
                    work0, work1, work2);
  else if (field1.memType == MemType::Float && field2.memType == MemType::Double)
    covariance_init(hasMissValues, gridsize, field1.vec_f, field2.vec_d, (float) field1.missval, field2.missval, nofvals, work0,
                    work1, work2);
  else if (field1.memType == MemType::Double && field2.memType == MemType::Float)
    covariance_init(hasMissValues, gridsize, field1.vec_d, field2.vec_f, field1.missval, (float) field2.missval, nofvals, work0,
                    work1, work2);
  else
    covariance_init(hasMissValues, gridsize, field1.vec_d, field2.vec_d, field1.missval, field2.missval, nofvals, work0, work1,
                    work2);
}

static size_t
covariance(size_t gridsize, double missval, const Varray<size_t> &nofvals, Varray<double> &work0, const Varray<double> &work1,
           const Varray<double> &work2)
{
  size_t nmiss = 0;

  for (size_t i = 0; i < gridsize; ++i)
    {
      auto missval1 = missval;
      auto missval2 = missval;
      double covar;
      auto nvals = nofvals[i];
      if (nvals > 0)
        {
          double dnvals = nvals;
          auto temp = DIVMN(MULMN(work0[i], work1[i]), dnvals * dnvals);
          covar = SUBMN(DIVMN(work2[i], dnvals), temp);
          if (DBL_IS_EQUAL(covar, missval)) nmiss++;
        }
      else
        {
          nmiss++;
          covar = missval;
        }

      work0[i] = covar;
    }

  return nmiss;
}

// rms in time
template <typename T1, typename T2>
static void
rmsd_init(size_t gridsize, const Varray<T1> &x, const Varray<T2> &y, T1 xmv, T2 ymv, Varray<size_t> &nofvals, Varray<double> &rmsd)
{
  for (size_t i = 0; i < gridsize; ++i)
    {
      if ((!DBL_IS_EQUAL(x[i], xmv)) && (!DBL_IS_EQUAL(y[i], ymv)))
        {
          const double xx = x[i];
          const double yy = y[i];
          rmsd[i] += ((xx - yy) * (xx - yy));
          nofvals[i]++;
        }
    }
}

static void
rmsd_init(size_t gridsize, const Field &field1, const Field &field2, Varray<size_t> &nofvals, Varray<double> &rmsd)
{
  if (field1.memType == MemType::Float && field2.memType == MemType::Float)
    rmsd_init(gridsize, field1.vec_f, field2.vec_f, (float) field1.missval, (float) field2.missval, nofvals, rmsd);
  else if (field1.memType == MemType::Float && field2.memType == MemType::Double)
    rmsd_init(gridsize, field1.vec_f, field2.vec_d, (float) field1.missval, field2.missval, nofvals, rmsd);
  else if (field1.memType == MemType::Double && field2.memType == MemType::Float)
    rmsd_init(gridsize, field1.vec_d, field2.vec_f, field1.missval, (float) field2.missval, nofvals, rmsd);
  else
    rmsd_init(gridsize, field1.vec_d, field2.vec_d, field1.missval, field2.missval, nofvals, rmsd);
}

static size_t
rmsd_compute(size_t gridsize, double missval, const Varray<size_t> &nofvals, Varray<double> &rmsd)
{
  size_t nmiss = 0;

  for (size_t i = 0; i < gridsize; ++i)
    {
      if (nofvals[i] > 0) { rmsd[i] = std::sqrt(rmsd[i] / (double) nofvals[i]); }
      else
        {
          nmiss++;
          rmsd[i] = missval;
        }
    }

  return nmiss;
}

void *
Timstat2(void *process)
{
  CdiDateTime vDateTime{};

  cdo_initialize(process);

  // clang-format off
  cdo_operator_add("timcor",   FieldFunc_Cor,   5, nullptr);
  cdo_operator_add("timcovar", FieldFunc_Covar, 3, nullptr);
  cdo_operator_add("timrmsd",  FieldFunc_Rmsd,  1, nullptr);
  // clang-format on

  auto operatorID = cdo_operator_id();
  auto operfunc = cdo_operator_f1(operatorID);
  auto nwork = cdo_operator_f2(operatorID);
  auto timeIsConst = (operfunc == FieldFunc_Rmsd);

  operator_check_argc(0);

  auto streamID1 = cdo_open_read(0);
  auto streamID2 = cdo_open_read(1);

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto vlistID2 = cdo_stream_inq_vlist(streamID2);
  auto vlistID3 = vlistDuplicate(vlistID1);

  vlist_compare(vlistID1, vlistID2, CMP_ALL);

  VarList varList1, varList2;
  varListInit(varList1, vlistID1);
  varListInit(varList2, vlistID2);

  auto nvars = vlistNvars(vlistID1);
  auto nrecs1 = vlistNrecs(vlistID1);
  std::vector<int> recVarID(nrecs1), recLevelID(nrecs1);

  auto taxisID1 = vlistInqTaxis(vlistID1);
  // auto taxisID2 = vlistInqTaxis(vlistID2);
  auto taxisID3 = taxisDuplicate(taxisID1);

  if (timeIsConst)
    for (int varID = 0; varID < nvars; ++varID) vlistDefVarTimetype(vlistID3, varID, TIME_CONSTANT);

  Field field1, field2;

  vlistDefTaxis(vlistID3, taxisID3);
  auto streamID3 = cdo_open_write(2);
  cdo_def_vlist(streamID3, vlistID3);

  Varray4D<double> work(nvars);
  Varray3D<size_t> nofvals(nvars);

  for (int varID = 0; varID < nvars; ++varID)
    {
      auto gridsize = varList1[varID].gridsize;
      auto nlevels = varList1[varID].nlevels;

      work[varID].resize(nlevels);
      nofvals[varID].resize(nlevels);

      for (int levelID = 0; levelID < nlevels; ++levelID)
        {
          nofvals[varID][levelID].resize(gridsize, 0);
          work[varID][levelID].resize(nwork);
          for (int i = 0; i < nwork; ++i) work[varID][levelID][i].resize(gridsize, 0.0);
        }
    }

  int tsID = 0;
  while (true)
    {
      auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
      if (nrecs == 0) break;

      vDateTime = taxisInqVdatetime(taxisID1);

      auto nrecs2 = cdo_stream_inq_timestep(streamID2, tsID);
      if (nrecs != nrecs2) cdo_warning("Input streams have different number of records!");

      for (int recID = 0; recID < nrecs; ++recID)
        {
          int varID, levelID;
          cdo_inq_record(streamID1, &varID, &levelID);
          cdo_inq_record(streamID2, &varID, &levelID);

          field1.init(varList1[varID]);
          field2.init(varList2[varID]);

          if (tsID == 0)
            {
              recVarID[recID] = varID;
              recLevelID[recID] = levelID;
            }

          auto gridsize = varList1[varID].gridsize;

          cdo_read_record(streamID1, field1);
          cdo_read_record(streamID2, field2);

          auto &rwork = work[varID][levelID];
          auto &rnofvals = nofvals[varID][levelID];

          if (operfunc == FieldFunc_Cor)
            {
              correlation_init(gridsize, field1, field2, rnofvals, rwork[0], rwork[1], rwork[2], rwork[3], rwork[4]);
            }
          else if (operfunc == FieldFunc_Covar)
            {
              covariance_init(gridsize, field1, field2, rnofvals, rwork[0], rwork[1], rwork[2]);
            }
          else if (operfunc == FieldFunc_Rmsd) { rmsd_init(gridsize, field1, field2, rnofvals, rwork[0]); }
        }

      tsID++;
    }

  tsID = 0;
  taxisDefVdatetime(taxisID3, vDateTime);
  cdo_def_timestep(streamID3, tsID);

  for (int recID = 0; recID < nrecs1; ++recID)
    {
      auto varID = recVarID[recID];
      auto levelID = recLevelID[recID];

      auto gridsize = varList1[varID].gridsize;
      auto missval = varList1[varID].missval;

      auto &rwork = work[varID][levelID];
      const auto &rnofvals = nofvals[varID][levelID];

      size_t nmiss = 0;
      if (operfunc == FieldFunc_Cor)
        {
          nmiss = correlation(gridsize, missval, rnofvals, rwork[0], rwork[1], rwork[2], rwork[3], rwork[4]);
        }
      else if (operfunc == FieldFunc_Covar) { nmiss = covariance(gridsize, missval, rnofvals, rwork[0], rwork[1], rwork[2]); }
      else if (operfunc == FieldFunc_Rmsd) { nmiss = rmsd_compute(gridsize, missval, rnofvals, rwork[0]); }

      cdo_def_record(streamID3, varID, levelID);
      cdo_write_record(streamID3, rwork[0].data(), nmiss);
    }

  cdo_stream_close(streamID3);
  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  cdo_finish();

  return nullptr;
}
