/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida
          Ralf Müller

*/

/*
   This module contains the following operators:

*/

#include <atomic>

#include <cdi.h>

#include "process_int.h"
#include "param_conversion.h"
#include "cdo_wtime.h"
#include <mpim_grid.h>
#include "grid_point_search.h"
#include "cdo_options.h"
#include "progress.h"
#include "cimdOmp.h"
#include "matrix_view.h"

template <typename T, typename CMP_FUNC>
T
fillmiss_kernel(int nfill, bool globgrid, long nx, long ny, long i, long j, T missval, MatrixView<T> &matrix1, CMP_FUNC is_EQ)
{
  if (!is_EQ(matrix1[j][i], missval)) return matrix1[j][i];

  T rval = missval;
  long ir, iu, il, io;
  long k1, k2;
  double s1, s2;

  long kr = 0, ku = 0, kl = 0, ko = 0;
  double xr = 0.0, xu = 0.0, xl = 0.0, xo = 0.0;

  for (ir = i + 1; ir < nx; ir++)
    if (!is_EQ(matrix1[j][ir], missval))
      {
        kr = ir - i;
        xr = matrix1[j][ir];
        break;
      }

  if (globgrid && ir == nx)
    {
      for (ir = 0; ir < i; ir++)
        if (!is_EQ(matrix1[j][ir], missval))
          {
            kr = nx + ir - i;
            xr = matrix1[j][ir];
            break;
          }
    }

  for (il = i - 1; il >= 0; il--)
    if (!is_EQ(matrix1[j][il], missval))
      {
        kl = i - il;
        xl = matrix1[j][il];
        break;
      }

  if (globgrid && il == -1)
    {
      for (il = nx - 1; il > i; il--)
        if (!is_EQ(matrix1[j][il], missval))
          {
            kl = nx + i - il;
            xl = matrix1[j][il];
            break;
          }
    }

  for (iu = j + 1; iu < ny; iu++)
    if (!is_EQ(matrix1[iu][i], missval))
      {
        ku = iu - j;
        xu = matrix1[iu][i];
        break;
      }

  for (io = j - 1; io >= 0; io--)
    if (!is_EQ(matrix1[io][i], missval))
      {
        ko = j - io;
        xo = matrix1[io][i];
        break;
      }

  // printf("%d %d %d %d %d %d %g %g %g %g\n", j,i,kr,kl,ku,ko,xr,xl,xu,xo);

  auto kh = kl + kr;
  auto kv = ko + ku;
  // clang-format off
  if      (kh == 0) { k1 = 0; s1 = 0.0; }
  else if (kl == 0) { k1 = 1; s1 = xr; }
  else if (kr == 0) { k1 = 1; s1 = xl; }
  else              { k1 = 2; s1 = xr * kl / kh + xl * kr / kh; }

  if      (kv == 0) { k2 = 0; s2 = 0.0; }
  else if (ku == 0) { k2 = 1; s2 = xo; }
  else if (ko == 0) { k2 = 1; s2 = xu; }
  else              { k2 = 2; s2 = xu * ko / kv + xo * ku / kv; }

  auto kk = k1 + k2;
  if (kk >= nfill)
    {
      if      (kk == 0) cdo_abort("no point found!");
      else if (k1 == 0) rval = s2;
      else if (k2 == 0) rval = s1;
      else              rval = s1 * k2 / kk + s2 * k1 / kk;
    }
  else
    rval = matrix1[j][i];
  // clang-format on

  return rval;
}

template <typename T, typename CMP_FUNC>
void
fillmiss_x(int gridID, Varray<T> &vIn, Varray<T> &vOut, T missval, int nfill, CMP_FUNC is_EQ)
{
  long nx = gridInqXsize(gridID);
  long ny = gridInqYsize(gridID);
  auto globgrid = (bool) gridIsCircular(gridID);

  auto gridtype = gridInqType(gridID);
  if (!(gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN)) cdo_abort("Unsupported grid type: %s!", gridNamePtr(gridtype));

  MatrixView<T> matrix1(vIn.data(), ny, nx);
  MatrixView<T> matrix2(vOut.data(), ny, nx);

#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
  for (long j = 0; j < ny; ++j)
    for (long i = 0; i < nx; ++i) { matrix2[j][i] = fillmiss_kernel(nfill, globgrid, nx, ny, i, j, missval, matrix1, is_EQ); }
}

static void
fillmiss(Field &field1, Field &field2, int nfill)
{
  if (field1.memType != field2.memType) cdo_abort("Interal error, memType of field1 and field2 differ!");

  if (std::isnan(field1.missval))
    {
      if (field1.memType == MemType::Float)
        fillmiss_x(field1.grid, field1.vec_f, field2.vec_f, (float) field1.missval, nfill, dbl_is_equal);
      else
        fillmiss_x(field1.grid, field1.vec_d, field2.vec_d, field1.missval, nfill, dbl_is_equal);
    }
  else
    {
      if (field1.memType == MemType::Float)
        fillmiss_x(field1.grid, field1.vec_f, field2.vec_f, (float) field1.missval, nfill, is_equal);
      else
        fillmiss_x(field1.grid, field1.vec_d, field2.vec_d, field1.missval, nfill, is_equal);
    }
}

template <typename T, typename CMP_FUNC>
T
fillmiss_one_step_kernel(long nx, long ny, long i, long j, T missval, MatrixView<T> &matrix1, CMP_FUNC is_EQ)
{
  if (!is_EQ(matrix1[j][i], missval)) return matrix1[j][i];

  T rval = missval;
  long ir, iu, il, io;
  long k1, k2;
  T s1, s2;

  long kr = 0, ku = 0, kl = 0, ko = 0;
  T xr = 0.0, xu = 0.0, xl = 0.0, xo = 0.0;

  for (ir = i + 1; ir < nx; ir++)
    if (!is_EQ(matrix1[j][ir], missval))
      {
        kr = ir - i;
        xr = matrix1[j][ir];
        break;
      }

  for (il = i - 1; il >= 0; il--)
    if (!is_EQ(matrix1[j][il], missval))
      {
        kl = i - il;
        xl = matrix1[j][il];
        break;
      }

  for (iu = j + 1; iu < ny; iu++)
    if (!is_EQ(matrix1[iu][i], missval))
      {
        ku = iu - j;
        xu = matrix1[iu][i];
        break;
      }

  for (io = j - 1; io >= 0; io--)
    if (!is_EQ(matrix1[io][i], missval))
      {
        ko = j - io;
        xo = matrix1[io][i];
        break;
      }

  auto kh = kl + kr;
  auto kv = ko + ku;
  // clang-format off
  if      (kh == 0) { s1 = 0.0; k1 = 0; }
  else if (kl == 0) { s1 = xr;  k1 = kr; }
  else if (kr == 0) { s1 = xl;  k1 = kl; }
  else              { s1 = (kl < kr) ? xl : xr;  k1 = (kl < kr) ? kl : kr; }

  if      (kv == 0) { s2 = 0.0; k2 = 0; }
  else if (ku == 0) { s2 = xo;  k2 = ko; }
  else if (ko == 0) { s2 = xu;  k2 = ku; }
  else              { s2 = (ku < ko) ? xu : xo;  k2 = (ku < ko) ? ku : ko; }

  auto kk = k1 + k2;
  if      (kk == 0) rval = matrix1[j][i];
  else if (k1 == 0) rval = s2;
  else if (k2 == 0) rval = s1;
  else              rval = (k1 <= k2) ? s1 : s2;
  // clang-format on

  return rval;
}

template <typename T, typename CMP_FUNC>
void
fillmiss_one_step_x(int gridID, Varray<T> &vIn, Varray<T> &vOut, T missval, int maxfill, CMP_FUNC is_EQ)
{
  long nx = gridInqXsize(gridID);
  long ny = gridInqYsize(gridID);

  MatrixView<T> matrix1(vIn.data(), ny, nx);
  MatrixView<T> matrix2(vOut.data(), ny, nx);

  for (int fill_iterations = 0; fill_iterations < maxfill; fill_iterations++)
    {
#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
      for (long j = 0; j < ny; ++j)
        for (long i = 0; i < nx; ++i) matrix2[j][i] = fillmiss_one_step_kernel(nx, ny, i, j, missval, matrix1, is_EQ);

      if ((fill_iterations + 1) < maxfill)
        for (long j = 0; j < ny; ++j)
          for (long i = 0; i < nx; ++i) matrix1[j][i] = matrix2[j][i];
    }
}

static void
fillmiss_one_step(Field &field1, Field &field2, int maxfill)
{
  if (field1.memType != field2.memType) cdo_abort("Interal error, memType of field1 and field2 differ!");

  if (std::isnan(field1.missval))
    {
      if (field1.memType == MemType::Float)
        fillmiss_one_step_x(field1.grid, field1.vec_f, field2.vec_f, (float) field1.missval, maxfill, dbl_is_equal);
      else
        fillmiss_one_step_x(field1.grid, field1.vec_d, field2.vec_d, field1.missval, maxfill, dbl_is_equal);
    }
  else
    {
      if (field1.memType == MemType::Float)
        fillmiss_one_step_x(field1.grid, field1.vec_f, field2.vec_f, (float) field1.missval, maxfill, is_equal);
      else
        fillmiss_one_step_x(field1.grid, field1.vec_d, field2.vec_d, field1.missval, maxfill, is_equal);
    }
}

template <typename T>
void
setmisstodis(size_t nmiss, int gridID, Varray<T> &vIn, Varray<T> &vOut, T missval, int numNeighbors)
{
  auto gridID0 = gridID;

  auto gridsize = gridInqSize(gridID);
  auto nvals = gridsize - nmiss;
  gridID = generate_full_point_grid(gridID);

  if (!gridHasCoordinates(gridID)) cdo_abort("Cell center coordinates missing!");

  Varray<double> xvals(gridsize), yvals(gridsize);
  gridInqXvals(gridID, xvals.data());
  gridInqYvals(gridID, yvals.data());

  // Convert lat/lon units if required
  cdo_grid_to_radian(gridID, CDI_XAXIS, gridsize, xvals.data(), "grid center lon");
  cdo_grid_to_radian(gridID, CDI_YAXIS, gridsize, yvals.data(), "grid center lat");

  std::vector<size_t> mindex(nmiss, 1), vindex(nvals, 1);
  Varray<double> lons(nvals), lats(nvals);

  size_t nv = 0, nm = 0;
  for (size_t i = 0; i < gridsize; ++i)
    {
      vOut[i] = vIn[i];
      if (DBL_IS_EQUAL(vIn[i], missval))
        {
          mindex[nm] = i;
          nm++;
        }
      else
        {
          if (nv < nvals)
            {
              lons[nv] = xvals[i];
              lats[nv] = yvals[i];
              vindex[nv] = i;
            }
          nv++;
        }
    }

  if (nv != nvals) cdo_abort("Internal problem, number of valid values differ!");

  std::vector<knnWeightsType> knnWeights;
  for (int i = 0; i < Threading::ompNumThreads; ++i) knnWeights.push_back(knnWeightsType(numNeighbors));

  auto start = Options::cdoVerbose ? cdo_get_wtime() : 0.0;

  GridPointSearch gps;

  if (nmiss)
    {
      auto xIsCyclic = false;
      size_t dims[2] = { nvals, 0 };
      grid_point_search_create(gps, xIsCyclic, dims, nvals, lons, lats);
      grid_point_search_extrapolate(gps);
    }

  if (Options::cdoVerbose) cdo_print("Point search created: %.2f seconds", cdo_get_wtime() - start);

  progress::init();

  start = Options::cdoVerbose ? cdo_get_wtime() : 0;

  std::atomic<size_t> atomicCount{ 0 };

#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
  for (size_t i = 0; i < nmiss; ++i)
    {
      atomicCount++;
      if (cdo_omp_get_thread_num() == 0) progress::update(0, 1, (double) atomicCount / nmiss);

      auto ompthID = cdo_omp_get_thread_num();

      grid_search_point(gps, xvals[mindex[i]], yvals[mindex[i]], knnWeights[ompthID]);

      // Compute weights based on inverse distance if mask is false, eliminate those points
      auto nadds = knnWeights[ompthID].computeWeights();
      if (nadds)
        {
          double result = 0.0;
          for (size_t n = 0; n < nadds; ++n) result += vIn[vindex[knnWeights[ompthID].m_addr[n]]] * knnWeights[ompthID].m_dist[n];
          vOut[mindex[i]] = result;
        }
    }

  progress::update(0, 1, 1);

  if (Options::cdoVerbose) cdo_print("Point search nearest: %.2f seconds", cdo_get_wtime() - start);

  grid_point_search_delete(gps);

  if (gridID0 != gridID) gridDestroy(gridID);
}

static void
setmisstodis(Field &field1, Field &field2, int numNeighbors)
{
  if (field1.memType != field2.memType) cdo_abort("Interal error, memType of field1 and field2 differ!");

  if (field1.memType == MemType::Float)
    setmisstodis(field1.nmiss, field1.grid, field1.vec_f, field2.vec_f, (float) field1.missval, numNeighbors);
  else
    setmisstodis(field1.nmiss, field1.grid, field1.vec_d, field2.vec_d, field1.missval, numNeighbors);
}

void *
Fillmiss(void *process)
{
  cdo_initialize(process);

  // clang-format off
  auto FILLMISS        = cdo_operator_add("fillmiss"   ,   0, 0, "nfill");
  auto FILLMISSONESTEP = cdo_operator_add("fillmiss2"  ,   0, 0, "nfill");
  auto SETMISSTONN     = cdo_operator_add("setmisstonn" ,  0, 0, "");
  auto SETMISSTODIS    = cdo_operator_add("setmisstodis" , 0, 0, "number of neighbors");

  auto operatorID = cdo_operator_id();

  void (*fill_method)(Field &, Field &, int) = &setmisstodis;
  if      (operatorID == FILLMISS)        fill_method = &fillmiss;
  else if (operatorID == FILLMISSONESTEP) fill_method = &fillmiss_one_step;
  else if (operatorID == SETMISSTONN)     fill_method = &setmisstodis;
  else if (operatorID == SETMISSTODIS)    fill_method = &setmisstodis;
  // clang-format on

  auto nfill = (operatorID == SETMISSTODIS) ? 4 : 1;

  // Argument handling
  auto oargc = cdo_operator_argc();
  if (oargc == 1)
    {
      nfill = parameter_to_int(cdo_operator_argv(0));
      if (operatorID == FILLMISS && (nfill < 1 || nfill > 4)) cdo_abort("nfill out of range!");
    }
  else if (oargc > 1)
    cdo_abort("Too many arguments!");

  auto streamID1 = cdo_open_read(0);

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto vlistID2 = vlistDuplicate(vlistID1);

  auto taxisID1 = vlistInqTaxis(vlistID1);
  auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  auto streamID2 = cdo_open_write(1);
  cdo_def_vlist(streamID2, vlistID2);

  VarList varList1;
  varListInit(varList1, vlistID1);

  Field field1, field2;

  int tsID = 0;
  while (true)
    {
      auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
      if (nrecs == 0) break;

      cdo_taxis_copy_timestep(taxisID2, taxisID1);

      cdo_def_timestep(streamID2, tsID);

      for (int recID = 0; recID < nrecs; ++recID)
        {
          int varID, levelID;
          cdo_inq_record(streamID1, &varID, &levelID);
          field1.init(varList1[varID]);
          cdo_read_record(streamID1, field1);

          cdo_def_record(streamID2, varID, levelID);

          if (field1.nmiss == 0) { cdo_write_record(streamID2, field1); }
          else
            {
              auto gridtype = gridInqType(varList1[varID].gridID);
              if ((operatorID == FILLMISS || operatorID == FILLMISSONESTEP)
                  && (gridtype == GRID_GME || gridtype == GRID_UNSTRUCTURED))
                cdo_abort("%s data unsupported!", gridNamePtr(gridtype));

              field2.init(varList1[varID]);

              fill_method(field1, field2, nfill);

              field2.nmiss = field_num_mv(field2);

              cdo_write_record(streamID2, field2);
            }
        }

      tsID++;
    }

  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  cdo_finish();

  return nullptr;
}
