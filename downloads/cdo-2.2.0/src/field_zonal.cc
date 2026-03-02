/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include <vector>
#include <algorithm>

#include "process_int.h"
#include "percentiles.h"
#include "field_functions.h"

template <typename T>
static void
varray_copy_zonal(size_t offset, size_t nx, const Varray<T> &v1, Varray<double> &v2)
{
  std::copy(v1.begin() + offset, v1.begin() + offset + nx, v2.begin());
}

static void
copy_latitude_row(size_t offset, size_t nx, const Field &field1, Varray<double> &v)
{
  if (field1.memType == MemType::Float)
    varray_copy_zonal(offset, nx, field1.vec_f, v);
  else
    varray_copy_zonal(offset, nx, field1.vec_d, v);
}

static size_t
fill_reduced_points(int gridID, size_t ny, std::vector<int> &reducedPoints, std::vector<int> &cumreducedPoints)
{
  reducedPoints.resize(ny);
  cumreducedPoints.resize(ny);
  gridInqReducedPoints(gridID, reducedPoints.data());
  cumreducedPoints[0] = 0;
  for (size_t j = 1; j < ny; ++j) cumreducedPoints[j] = cumreducedPoints[j - 1] + reducedPoints[j - 1];
  size_t nx = reducedPoints[0];
  for (size_t j = 1; j < ny; ++j)
    if (reducedPoints[j] > (int) nx) nx = reducedPoints[j];
  return nx;
}

using funcType = double(size_t, const Varray<double> &);
using funcTypeMV = double(size_t, const Varray<double> &, double);
using funcTypeNmissMV = double(size_t, const Varray<double> &, size_t, double);

static void
zonal_kernel_1(const Field &field1, Field &field2, funcTypeNmissMV funcNmissMV)
{
  size_t rnmiss = 0;
  auto nmiss = field1.nmiss;
  auto missval = field1.missval;

  auto ny = gridInqYsize(field1.grid);
  auto isReducedGrid = (gridInqType(field1.grid) == GRID_GAUSSIAN_REDUCED);
  std::vector<int> reducedPoints, cumreducedPoints;
  auto nx = isReducedGrid ? fill_reduced_points(field1.grid, ny, reducedPoints, cumreducedPoints) : gridInqXsize(field1.grid);

  Varray<double> v(nx);

  for (size_t j = 0; j < ny; ++j)
    {
      if (isReducedGrid) nx = reducedPoints[j];
      size_t offset = isReducedGrid ? cumreducedPoints[j] : j * nx;
      copy_latitude_row(offset, nx, field1, v);

      auto result = funcNmissMV(nx, v, nmiss, missval);
      if (dbl_is_equal(result, missval)) rnmiss++;
      field2.vec_d[j] = result;
    }

  field2.nmiss = rnmiss;
}

static void
zonal_kernel_2(const Field &field1, Field &field2, funcType func, funcTypeMV funcMV)
{
  size_t rnmiss = 0;
  auto nmiss = field1.nmiss;
  auto missval = field1.missval;

  auto ny = gridInqYsize(field1.grid);
  auto isReducedGrid = (gridInqType(field1.grid) == GRID_GAUSSIAN_REDUCED);
  std::vector<int> reducedPoints, cumreducedPoints;
  auto nx = isReducedGrid ? fill_reduced_points(field1.grid, ny, reducedPoints, cumreducedPoints) : gridInqXsize(field1.grid);

  Varray<double> v(nx);

  for (size_t j = 0; j < ny; ++j)
    {
      if (isReducedGrid) nx = reducedPoints[j];
      size_t offset = isReducedGrid ? cumreducedPoints[j] : j * nx;
      copy_latitude_row(offset, nx, field1, v);

      auto result = nmiss ? funcMV(nx, v, missval) : func(nx, v);
      if (dbl_is_equal(result, missval)) rnmiss++;
      field2.vec_d[j] = result;
    }

  field2.nmiss = rnmiss;
}

void
zonal_min(const Field &field1, Field &field2)
{
  zonal_kernel_2(field1, field2, varray_min, varray_min_mv);
}

void
zonal_max(const Field &field1, Field &field2)
{
  zonal_kernel_2(field1, field2, varray_max, varray_max_mv);
}

void
zonal_range(const Field &field1, Field &field2)
{
  zonal_kernel_2(field1, field2, varray_range, varray_range_mv);
}

void
zonal_sum(const Field &field1, Field &field2)
{
  zonal_kernel_2(field1, field2, varray_sum, varray_sum_mv);
}

void
zonal_mean(const Field &field1, Field &field2)
{
  zonal_kernel_2(field1, field2, varray_mean, varray_mean_mv);
}

void
zonal_avg(const Field &field1, Field &field2)
{
  zonal_kernel_2(field1, field2, varray_mean, varray_avg_mv);
}

void
zonal_var(const Field &field1, Field &field2)
{
  zonal_kernel_1(field1, field2, varray_var);
}

void
zonal_var1(const Field &field1, Field &field2)
{
  zonal_kernel_1(field1, field2, varray_var_1);
}

void
zonal_std(const Field &field1, Field &field2)
{
  size_t rnmiss = 0;
  auto missval = field2.missval;
  auto ny = gridInqYsize(field2.grid);

  zonal_var(field1, field2);

  for (size_t j = 0; j < ny; ++j)
    {
      auto rstd = var_to_std(field2.vec_d[j], missval);
      if (dbl_is_equal(rstd, missval)) rnmiss++;
      field2.vec_d[j] = rstd;
    }

  field2.nmiss = rnmiss;
}

void
zonal_std1(const Field &field1, Field &field2)
{
  size_t rnmiss = 0;
  auto missval = field2.missval;
  auto ny = gridInqYsize(field2.grid);

  zonal_var1(field1, field2);

  for (size_t j = 0; j < ny; ++j)
    {
      auto rstd = var_to_std(field2.vec_d[j], missval);
      if (dbl_is_equal(rstd, missval)) rnmiss++;
      field2.vec_d[j] = rstd;
    }

  field2.nmiss = rnmiss;
}

void
zonal_skew(const Field &field1, Field &field2)
{
  zonal_kernel_1(field1, field2, varray_skew);
}

void
zonal_kurt(const Field &field1, Field &field2)
{
  zonal_kernel_1(field1, field2, varray_kurt);
}

void
zonal_median(const Field &field1, Field &field2)
{
  zonal_kernel_1(field1, field2, varray_median);
}

void
zonal_pctl(const Field &field1, Field &field2, double pn)
{
  size_t rnmiss = 0;
  auto missval = field2.missval;

  auto ny = gridInqYsize(field1.grid);
  auto isReducedGrid = (gridInqType(field1.grid) == GRID_GAUSSIAN_REDUCED);
  std::vector<int> reducedPoints, cumreducedPoints;
  auto nx = isReducedGrid ? fill_reduced_points(field1.grid, ny, reducedPoints, cumreducedPoints) : gridInqXsize(field1.grid);

  Varray<double> v(nx);

  if (field1.nmiss)
    {
      for (size_t j = 0; j < ny; ++j)
        {
          if (isReducedGrid) nx = reducedPoints[j];
          size_t offset = isReducedGrid ? cumreducedPoints[j] : j * nx;
          copy_latitude_row(offset, nx, field1, v);

          size_t k = 0;
          for (size_t i = 0; i < nx; ++i)
            if (!dbl_is_equal(v[i], missval)) v[k++] = v[i];

          if (k > 0) { field2.vec_d[j] = percentile(v.data(), k, pn); }
          else
            {
              field2.vec_d[j] = missval;
              rnmiss++;
            }
        }
    }
  else
    {
      for (size_t j = 0; j < ny; ++j)
        {
          if (isReducedGrid) nx = reducedPoints[j];
          size_t offset = isReducedGrid ? cumreducedPoints[j] : j * nx;
          copy_latitude_row(offset, nx, field1, v);

          if (nx > 0) { field2.vec_d[j] = percentile(v.data(), nx, pn); }
          else
            {
              field2.vec_d[j] = missval;
              rnmiss++;
            }
        }
    }

  field2.nmiss = rnmiss;
}

void
zonal_function(const Field &field1, Field &field2, int function)
{
  // clang-format off
  switch (function)
    {
    case FieldFunc_Min:    return zonal_min(field1, field2);
    case FieldFunc_Max:    return zonal_max(field1, field2);
    case FieldFunc_Range:  return zonal_range(field1, field2);
    case FieldFunc_Sum:    return zonal_sum(field1, field2);
    case FieldFunc_Mean:   return zonal_mean(field1, field2);
    case FieldFunc_Avg:    return zonal_avg(field1, field2);
    case FieldFunc_Std:    return zonal_std(field1, field2);
    case FieldFunc_Std1:   return zonal_std1(field1, field2);
    case FieldFunc_Var:    return zonal_var(field1, field2);
    case FieldFunc_Var1:   return zonal_var1(field1, field2);
    case FieldFunc_Skew:   return zonal_skew(field1, field2);
    case FieldFunc_Kurt:   return zonal_kurt(field1, field2);
    case FieldFunc_Median: return zonal_median(field1, field2);
    default: cdo_abort("%s: function %d not implemented!",  __func__, function);
    }
  // clang-format on
}
