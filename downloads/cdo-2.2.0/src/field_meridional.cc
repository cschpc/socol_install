/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include "process_int.h"
#include "percentiles.h"
#include "field_functions.h"

using funcType1 = double(size_t, const Varray<double> &);
using funcTypeMV1 = double(size_t, const Varray<double> &, double);
using funcType2 = double(size_t, const Varray<double> &, const Varray<double> &, double);
using funcTypeMV2 = double(size_t, const Varray<double> &, const Varray<double> &, double);
using funcType3 = double(size_t, const Varray<double> &, const Varray<double> &, size_t, double);
using funcType4 = double(size_t, const Varray<double> &, size_t, double);

template <typename T>
static void
varray_copy_meridional(size_t i, size_t nx, size_t ny, const Varray<T> &v1, Varray<double> &v2)
{
  for (size_t j = 0; j < ny; ++j) v2[j] = v1[j * nx + i];
}

static void
meridional_kernel_1(const Field &field1, Field &field2, funcType1 func, funcTypeMV1 funcMV)
{
  size_t rnmiss = 0;
  const auto nmiss = field1.nmiss;
  const auto missval = field1.missval;
  const auto nx = gridInqXsize(field1.grid);
  const auto ny = gridInqYsize(field1.grid);

  Varray<double> v(ny);

  for (size_t i = 0; i < nx; ++i)
    {
      if (field1.memType == MemType::Float)
        varray_copy_meridional(i, nx, ny, field1.vec_f, v);
      else
        varray_copy_meridional(i, nx, ny, field1.vec_d, v);

      const auto result = nmiss ? funcMV(ny, v, missval) : func(ny, v);
      if (DBL_IS_EQUAL(result, missval)) rnmiss++;
      field2.vec_d[i] = result;
    }

  field2.nmiss = rnmiss;
}

static void
meridional_kernel_2(const Field &field1, Field &field2, funcType2 func, funcTypeMV2 funcMV)
{
  size_t rnmiss = 0;
  const auto nmiss = field1.nmiss;
  const auto missval = field1.missval;
  const auto nx = gridInqXsize(field1.grid);
  const auto ny = gridInqYsize(field1.grid);

  Varray<double> v(ny), w(ny);

  for (size_t i = 0; i < nx; ++i)
    {
      varray_copy_meridional(i, nx, ny, field1.weightv, w);
      if (field1.memType == MemType::Float)
        varray_copy_meridional(i, nx, ny, field1.vec_f, v);
      else
        varray_copy_meridional(i, nx, ny, field1.vec_d, v);

      const auto result = nmiss ? funcMV(ny, v, w, missval) : func(ny, v, w, missval);
      if (DBL_IS_EQUAL(result, missval)) rnmiss++;
      field2.vec_d[i] = result;
    }

  field2.nmiss = rnmiss;
}

static void
meridional_kernel_3(const Field &field1, Field &field2, funcType3 func)
{
  size_t rnmiss = 0;
  const auto nmiss = field1.nmiss;
  const auto missval = field1.missval;
  const auto nx = gridInqXsize(field1.grid);
  const auto ny = gridInqYsize(field1.grid);

  Varray<double> v(ny), w(ny);

  for (size_t i = 0; i < nx; ++i)
    {
      varray_copy_meridional(i, nx, ny, field1.weightv, w);
      if (field1.memType == MemType::Float)
        varray_copy_meridional(i, nx, ny, field1.vec_f, v);
      else
        varray_copy_meridional(i, nx, ny, field1.vec_d, v);

      const auto result = func(ny, v, w, nmiss, missval);
      if (DBL_IS_EQUAL(result, missval)) rnmiss++;
      field2.vec_d[i] = result;
    }

  field2.nmiss = rnmiss;
}

static void
meridional_kernel_4(const Field &field1, Field &field2, funcType4 func)
{
  size_t rnmiss = 0;
  const auto nmiss = field1.nmiss;
  const auto missval = field1.missval;
  const auto nx = gridInqXsize(field1.grid);
  const auto ny = gridInqYsize(field1.grid);

  Varray<double> v(ny);

  for (size_t i = 0; i < nx; ++i)
    {
      if (field1.memType == MemType::Float)
        varray_copy_meridional(i, nx, ny, field1.vec_f, v);
      else
        varray_copy_meridional(i, nx, ny, field1.vec_d, v);

      const auto result = func(ny, v, nmiss, missval);
      if (DBL_IS_EQUAL(result, missval)) rnmiss++;
      field2.vec_d[i] = result;
    }

  field2.nmiss = rnmiss;
}

static void
meridional_min(const Field &field1, Field &field2)
{
  meridional_kernel_1(field1, field2, varray_min, varray_min_mv);
}

static void
meridional_max(const Field &field1, Field &field2)
{
  meridional_kernel_1(field1, field2, varray_max, varray_max_mv);
}

static void
meridional_range(const Field &field1, Field &field2)
{
  meridional_kernel_1(field1, field2, varray_range, varray_range_mv);
}

static void
meridional_sum(const Field &field1, Field &field2)
{
  meridional_kernel_1(field1, field2, varray_sum, varray_sum_mv);
}

static void
meridional_meanw(const Field &field1, Field &field2)
{
  meridional_kernel_2(field1, field2, varray_weighted_mean, varray_weighted_mean_mv);
}

static void
meridional_avgw(const Field &field1, Field &field2)
{
  meridional_kernel_2(field1, field2, varray_weighted_mean, varray_weighted_avg_mv);
}

static void
meridional_varw(const Field &field1, Field &field2)
{
  meridional_kernel_3(field1, field2, varray_weighted_var);
}

static void
meridional_var1w(const Field &field1, Field &field2)
{
  meridional_kernel_3(field1, field2, varray_weighted_var_1);
}

static void
meridional_stdw(const Field &field1, Field &field2)
{
  size_t rnmiss = 0;
  const auto missval = field1.missval;

  const auto nx = gridInqXsize(field1.grid);

  meridional_varw(field1, field2);

  for (size_t i = 0; i < nx; ++i)
    {
      const auto rstd = var_to_std(field2.vec_d[i], missval);
      if (DBL_IS_EQUAL(rstd, missval)) rnmiss++;
      field2.vec_d[i] = rstd;
    }

  field2.nmiss = rnmiss;
}

static void
meridional_std1w(const Field &field1, Field &field2)
{
  size_t rnmiss = 0;
  const auto missval = field1.missval;

  const auto nx = gridInqXsize(field1.grid);

  meridional_var1w(field1, field2);

  for (size_t i = 0; i < nx; ++i)
    {
      const auto rstd = var_to_std(field2.vec_d[i], missval);
      if (DBL_IS_EQUAL(rstd, missval)) rnmiss++;
      field2.vec_d[i] = rstd;
    }

  field2.nmiss = rnmiss;
}

static void
meridional_skew(const Field &field1, Field &field2)
{
  meridional_kernel_4(field1, field2, varray_skew);
}

static void
meridional_kurt(const Field &field1, Field &field2)
{
  meridional_kernel_4(field1, field2, varray_kurt);
}

static void
meridional_median(const Field &field1, Field &field2)
{
  meridional_kernel_4(field1, field2, varray_median);
}

void
meridional_pctl(const Field &field1, Field &field2, const double pn)
{
  size_t rnmiss = 0;
  const auto missval = field1.missval;

  const auto nx = gridInqXsize(field1.grid);
  const auto ny = gridInqYsize(field1.grid);

  Varray<double> v(ny);

  if (field1.nmiss)
    {
      for (size_t i = 0; i < nx; ++i)
        {
          size_t k = 0;
          if (field1.memType == MemType::Float)
            {
              for (size_t j = 0; j < ny; ++j)
                if (!DBL_IS_EQUAL(field1.vec_d[j * nx + i], missval)) v[k++] = field1.vec_f[j * nx + i];
            }
          else
            {
              for (size_t j = 0; j < ny; ++j)
                if (!DBL_IS_EQUAL(field1.vec_d[j * nx + i], missval)) v[k++] = field1.vec_d[j * nx + i];
            }

          if (k > 0) { field2.vec_d[i] = percentile(v.data(), k, pn); }
          else
            {
              field2.vec_d[i] = missval;
              rnmiss++;
            }
        }
    }
  else
    {
      for (size_t i = 0; i < nx; ++i)
        {
          if (ny > 0)
            {
              if (field1.memType == MemType::Float)
                varray_copy_meridional(i, nx, ny, field1.vec_f, v);
              else
                varray_copy_meridional(i, nx, ny, field1.vec_d, v);

              field2.vec_d[i] = percentile(v.data(), ny, pn);
            }
          else
            {
              field2.vec_d[i] = missval;
              rnmiss++;
            }
        }
    }

  field2.nmiss = rnmiss;
}

void
meridional_function(const Field &field1, Field &field2, int function)
{
  // clang-format off
  switch (function)
    {
    case FieldFunc_Min:    return meridional_min(field1, field2);
    case FieldFunc_Max:    return meridional_max(field1, field2);
    case FieldFunc_Range:  return meridional_range(field1, field2);
    case FieldFunc_Sum:    return meridional_sum(field1, field2);
    case FieldFunc_Meanw:  return meridional_meanw(field1, field2);
    case FieldFunc_Avgw:   return meridional_avgw(field1, field2);
    case FieldFunc_Stdw:   return meridional_stdw(field1, field2);
    case FieldFunc_Std1w:  return meridional_std1w(field1, field2);
    case FieldFunc_Varw:   return meridional_varw(field1, field2);
    case FieldFunc_Var1w:  return meridional_var1w(field1, field2);
    case FieldFunc_Skew:   return meridional_skew(field1, field2);
    case FieldFunc_Kurt:   return meridional_kurt(field1, field2);
    case FieldFunc_Median: return meridional_median(field1, field2);
    default: cdo_abort("%s: function %d not implemented!", __func__, function);
    }
  // clang-format on
}
