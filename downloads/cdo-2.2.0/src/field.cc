/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include <cassert>
#include <algorithm>

#include "percentiles.h"
#include "varray.h"
#include "field_functions.h"
#include "cdo_output.h"

void
Field::init(const CdoVar &var)
{
  fpeRaised = 0;
  nwpv = 1;
  grid = var.gridID;
  gridsize = var.gridsize;
  nmiss = 0;
  missval = var.missval;
  memType = var.memType;
  size = var.gridsize * var.nwpv;
  m_count = size;
  if (memType == MemType::Float)
    varrayResize(vec_f, size);
  else
    varrayResize(vec_d, size);
}

void
Field::resize(const size_t count)
{
  memType = MemType::Double;
  m_count = count;
  varrayResize(vec_d, m_count);
  if (!size) size = m_count;
}

void
Field::resize(const size_t count, const double value)
{
  memType = MemType::Double;
  m_count = count;
  varrayResizeInit(vec_d, m_count, value);
  if (!size) size = m_count;
}

void
Field::resizef(const size_t count)
{
  memType = MemType::Float;
  m_count = count;
  varrayResize(vec_f, m_count);
  if (!size) size = m_count;
}

void
Field::resizef(const size_t count, const float value)
{
  memType = MemType::Float;
  m_count = count;
  varrayResizeInit(vec_f, m_count, value);
  if (!size) size = m_count;
}

bool
Field::empty() const
{
  return m_count == 0;
}

void
Field::check_gridsize() const
{
  if (size == 0) fprintf(stderr, "Internal problem, size of field not set!\n");
  if (size > m_count) fprintf(stderr, "Internal problem, size of field is greater than allocated size of field!\n");
}

void
Field3D::init(const CdoVar &var)
{
  nlevels = var.nlevels;
  grid = var.gridID;
  gridsize = var.gridsize;
  missval = var.missval;
  memType = var.memType;
  size = var.nlevels * var.gridsize * var.nwpv;
  if (memType == MemType::Float)
    varrayResize(vec_f, size);
  else
    varrayResize(vec_d, size);
}

void
field_fill(Field &field, double value)
{
  field.check_gridsize();

  if (field.memType == MemType::Float)
    std::fill(field.vec_f.begin(), field.vec_f.begin() + field.size, value);
  else
    std::fill(field.vec_d.begin(), field.vec_d.begin() + field.size, value);
}

void
field_ncopy(size_t n, const Field &fieldIn, Field &fieldOut)
{
  if (n > fieldIn.size) cdo_abort("Source field to small (%s)", __func__);
  if (n > fieldOut.size) cdo_abort("Target field to small (%s)", __func__);

  fieldOut.nmiss = fieldIn.nmiss;

  // clang-format off
  if      (memtype_is_float_float  (fieldIn.memType, fieldOut.memType)) varray_copy(n, fieldIn.vec_f, fieldOut.vec_f);
  else if (memtype_is_float_double (fieldIn.memType, fieldOut.memType)) varray_copy(n, fieldIn.vec_f, fieldOut.vec_d);
  else if (memtype_is_double_float (fieldIn.memType, fieldOut.memType)) varray_copy(n, fieldIn.vec_d, fieldOut.vec_f);
  else if (memtype_is_double_double(fieldIn.memType, fieldOut.memType)) varray_copy(n, fieldIn.vec_d, fieldOut.vec_d);
  else cdo_abort("Internal error, unexpected memType");
  // clang-format on
}

void
field_copy(const Field &fieldIn, Field &fieldOut)
{
  if (fieldIn.size > fieldOut.size) cdo_abort("Target field to small (%s)", __func__);

  fieldOut.nmiss = fieldIn.nmiss;

  if (fieldIn.memType == MemType::Float)
    fieldOut.vec_f = fieldIn.vec_f;
  else
    fieldOut.vec_d = fieldIn.vec_d;
}

void
field_copy(const Field3D &fieldIn, Field3D &fieldOut)
{
  if (fieldIn.size > fieldOut.size) cdo_abort("Target field to small (%s)", __func__);

  fieldOut.nmiss = fieldIn.nmiss;

  if (fieldIn.memType == MemType::Float)
    fieldOut.vec_f = fieldIn.vec_f;
  else
    fieldOut.vec_d = fieldIn.vec_d;
}

void
field_copy(const Field3D &fieldIn, const int levelID, Field &fieldOut)
{
  auto size = fieldIn.gridsize * fieldIn.nwpv;
  auto offset = levelID * size;
  if (fieldIn.memType == MemType::Float)
    std::copy(fieldIn.vec_f.begin() + offset, fieldIn.vec_f.begin() + offset + size, fieldOut.vec_f.begin());
  else
    std::copy(fieldIn.vec_d.begin() + offset, fieldIn.vec_d.begin() + offset + size, fieldOut.vec_d.begin());
}

void
field_add(Field &field1, const Field3D &field2, const int levelID)
{
  auto size = field1.gridsize * field1.nwpv;
  auto offset = levelID * size;
  if (field1.memType == MemType::Float)
    for (size_t i = 0; i < size; ++i) field1.vec_f[i] += field2.vec_f[offset + i];
  else
    for (size_t i = 0; i < size; ++i) field1.vec_d[i] += field2.vec_d[offset + i];
}

// functor that returns true if value is equal to the value of the constructor parameter provided
class valueDblIsEqual
{
  double _missval;

public:
  explicit valueDblIsEqual(double missval) : _missval(missval) {}
  bool
  operator()(const double value) const
  {
    return dbl_is_equal(value, _missval);
  }
};

// functor that returns true if value is equal to the value of the constructor parameter provided
class valueIsEqual
{
  double _missval;

public:
  explicit valueIsEqual(double missval) : _missval(missval) {}
  bool
  operator()(const double value) const
  {
    return is_equal(value, _missval);
  }
};

size_t
field_num_miss(const Field &field)
{
  field.check_gridsize();

  auto missval = field.missval;
  const auto &v = field.vec_d;

  if (std::isnan(missval))
    return std::count_if(v.begin(), v.begin() + field.size, valueDblIsEqual(missval));
  else
    return std::count_if(v.begin(), v.begin() + field.size, valueIsEqual(missval));
}

size_t
field_num_mv(Field &field)
{
  if (field.memType == MemType::Float)
    field.nmiss = varray_num_mv(field.size, field.vec_f, (float) field.missval);
  else
    field.nmiss = varray_num_mv(field.size, field.vec_d, field.missval);

  return field.nmiss;
}

MinMax
field_min_max(const Field &field)
{
  if (field.memType == MemType::Float)
    return field.nmiss ? varray_min_max_mv(field.size, field.vec_f, (float) field.missval) : varray_min_max(field.vec_f);
  else
    return field.nmiss ? varray_min_max_mv(field.size, field.vec_d, field.missval) : varray_min_max(field.vec_d);
}

double
field_min(const Field &field)
{
  if (field.memType == MemType::Float)
    return field.nmiss ? varray_min_mv(field.size, field.vec_f, (float) field.missval) : varray_min(field.size, field.vec_f);
  else
    return field.nmiss ? varray_min_mv(field.size, field.vec_d, field.missval) : varray_min(field.size, field.vec_d);
}

double
field_max(const Field &field)
{
  if (field.memType == MemType::Float)
    return field.nmiss ? varray_max_mv(field.size, field.vec_f, (float) field.missval) : varray_max(field.size, field.vec_f);
  else
    return field.nmiss ? varray_max_mv(field.size, field.vec_d, field.missval) : varray_max(field.size, field.vec_d);
}

double
field_range(const Field &field)
{
  if (field.memType == MemType::Float)
    return field.nmiss ? varray_range_mv(field.size, field.vec_f, (float) field.missval) : varray_range(field.size, field.vec_f);
  else
    return field.nmiss ? varray_range_mv(field.size, field.vec_d, field.missval) : varray_range(field.size, field.vec_d);
}

double
field_sum(const Field &field)
{
  if (field.memType == MemType::Float)
    return field.nmiss ? varray_sum_mv(field.size, field.vec_f, (float) field.missval) : varray_sum(field.size, field.vec_f);
  else
    return field.nmiss ? varray_sum_mv(field.size, field.vec_d, field.missval) : varray_sum(field.size, field.vec_d);
}

double
field_mean(const Field &field)
{
  if (field.memType == MemType::Float)
    return field.nmiss ? varray_mean_mv(field.size, field.vec_f, (float) field.missval) : varray_mean(field.size, field.vec_f);
  else
    return field.nmiss ? varray_mean_mv(field.size, field.vec_d, field.missval) : varray_mean(field.size, field.vec_d);
}

double
field_meanw(const Field &field)
{
  if (field.memType == MemType::Float)
    return field.nmiss ? varray_weighted_mean_mv(field.size, field.vec_f, field.weightv, (float) field.missval)
                       : varray_weighted_mean(field.size, field.vec_f, field.weightv, (float) field.missval);
  else
    return field.nmiss ? varray_weighted_mean_mv(field.size, field.vec_d, field.weightv, field.missval)
                       : varray_weighted_mean(field.size, field.vec_d, field.weightv, field.missval);
}

double
field_avg(const Field &field)
{
  if (field.memType == MemType::Float)
    return field.nmiss ? varray_avg_mv(field.size, field.vec_f, (float) field.missval) : varray_mean(field.size, field.vec_f);
  else
    return field.nmiss ? varray_avg_mv(field.size, field.vec_d, field.missval) : varray_mean(field.size, field.vec_d);
}

double
field_avgw(const Field &field)
{
  if (field.memType == MemType::Float)
    return field.nmiss ? varray_weighted_avg_mv(field.size, field.vec_f, field.weightv, (float) field.missval)
                       : varray_weighted_mean(field.size, field.vec_f, field.weightv, (float) field.missval);
  else
    return field.nmiss ? varray_weighted_avg_mv(field.size, field.vec_d, field.weightv, field.missval)
                       : varray_weighted_mean(field.size, field.vec_d, field.weightv, field.missval);
}

double
field_var(const Field &field)
{
  if (field.memType == MemType::Float)
    return varray_var(field.size, field.vec_f, field.nmiss, (float) field.missval);
  else
    return varray_var(field.size, field.vec_d, field.nmiss, field.missval);
}

double
field_var1(const Field &field)
{
  if (field.memType == MemType::Float)
    return varray_var_1(field.size, field.vec_f, field.nmiss, (float) field.missval);
  else
    return varray_var_1(field.size, field.vec_d, field.nmiss, field.missval);
}

double
field_skew(const Field &field)
{
  if (field.memType == MemType::Float)
    return varray_skew(field.size, field.vec_f, field.nmiss, (float) field.missval);
  else
    return varray_skew(field.size, field.vec_d, field.nmiss, field.missval);
}

double
field_kurt(const Field &field)
{
  if (field.memType == MemType::Float)
    return varray_kurt(field.size, field.vec_f, field.nmiss, (float) field.missval);
  else
    return varray_kurt(field.size, field.vec_d, field.nmiss, field.missval);
}

double
field_median(const Field &field)
{
  if (field.memType == MemType::Float)
    return varray_median(field.size, field.vec_f, field.nmiss, (float) field.missval);
  else
    return varray_median(field.size, field.vec_d, field.nmiss, field.missval);
}

double
field_count(const Field &field)
{
  if (field.memType == MemType::Float)
    return varray_count(field.size, field.vec_f, field.nmiss, (float) field.missval);
  else
    return varray_count(field.size, field.vec_d, field.nmiss, field.missval);
}

double
field_varw(const Field &field)
{
  if (field.memType == MemType::Float)
    return varray_weighted_var(field.size, field.vec_f, field.weightv, field.nmiss, (float) field.missval);
  else
    return varray_weighted_var(field.size, field.vec_d, field.weightv, field.nmiss, field.missval);
}

double
field_var1w(const Field &field)
{
  if (field.memType == MemType::Float)
    return varray_weighted_var_1(field.size, field.vec_f, field.weightv, field.nmiss, (float) field.missval);
  else
    return varray_weighted_var_1(field.size, field.vec_d, field.weightv, field.nmiss, field.missval);
}

double
var_to_std(double rvar, double missval)
{
  if (dbl_is_equal(rvar, missval) || rvar < 0) return missval;

  return IS_NOT_EQUAL(rvar, 0) ? std::sqrt(rvar) : 0;
}

double
field_std(const Field &field)
{
  return var_to_std(field_var(field), field.missval);
}

double
field_std1(const Field &field)
{
  return var_to_std(field_var1(field), field.missval);
}

double
field_stdw(const Field &field)
{
  return var_to_std(field_varw(field), field.missval);
}

double
field_std1w(const Field &field)
{
  return var_to_std(field_var1w(field), field.missval);
}

void
field_rms(const Field &field, const Field &field2, Field &field3)
{
  size_t rnmiss = 0;
  auto grid1 = field.grid;
  //  size_t nmiss1   = field.nmiss;
  const auto array1 = field.vec_d.data();
  auto grid2 = field2.grid;
  //  size_t nmiss2   = field2.nmiss;
  const auto array2 = field2.vec_d.data();
  auto missval1 = field.missval;
  auto missval2 = field2.missval;
  const auto &w = field.weightv;
  auto rsum = 0.0, rsumw = 0.0;

  auto len = gridInqSize(grid1);
  if (len != gridInqSize(grid2)) cdo_abort("fields have different size!");

  // if ( nmiss1 )
  {
    for (size_t i = 0; i < len; ++i)
      if (!dbl_is_equal(w[i], missval1))
        {
          rsum = ADDMN(rsum, MULMN(w[i], MULMN(SUBMN(array2[i], array1[i]), SUBMN(array2[i], array1[i]))));
          rsumw = ADDMN(rsumw, w[i]);
        }
  }
  /*
else
  {
    for ( i = 0; i < len; i++ )
      {
        rsum  += w[i] * array1[i];
        rsumw += w[i];
      }
  }
  */

  auto ravg = SQRTMN(DIVMN(rsum, rsumw));

  if (dbl_is_equal(ravg, missval1)) rnmiss++;

  field3.vec_d[0] = ravg;
  field3.nmiss = rnmiss;
}

template <typename T>
double
array_pctl(size_t len, T *array, size_t nmiss, T missval, const double pn)
{
  double pctl = missval;

  if (len != nmiss)
    {
      if (nmiss)
        {
          Varray<T> v(len);

          size_t j = 0;
          for (size_t i = 0; i < len; ++i)
            if (!dbl_is_equal(array[i], missval)) v[j++] = array[i];

          if (nmiss != len - j)
            cdo_warning("Internal problem, inconsistent number of missing values (nmiss: exprected=%zu found=%zu!)", nmiss,
                        len - j);

          pctl = percentile(v.data(), j, pn);
        }
      else { pctl = percentile(array, len, pn); }
    }

  return pctl;
}

double
field_pctl(Field &field, const double pn)
{
  if (field.memType == MemType::Float)
    return array_pctl(field.size, field.vec_f.data(), field.nmiss, (float) field.missval, pn);
  else
    return array_pctl(field.size, field.vec_d.data(), field.nmiss, field.missval, pn);
}

static int
compare_double(const void *const a, const void *const b)
{
  const double *const x = (const double *) a;
  const double *const y = (const double *) b;
  return ((*x > *y) - (*x < *y)) * 2 + (*x > *y) - (*x < *y);
}

double
field_rank(Field &field)
{
  auto res = 0.0;
  // Using first value as reference (observation)
  auto val = field.vec_d[0];
  const auto array = &field.vec_d[1];
  auto len = field.size - 1;

  if (field.nmiss) return field.missval;

  std::qsort(array, len, sizeof(double), compare_double);

  if (val > array[len - 1])
    res = (double) len;
  else
    for (size_t j = 0; j < len; ++j)
      if (array[j] >= val)
        {
          res = (double) j;
          break;
        }

  return res;
}

double
field_function(const Field &field, int function)
{
  // clang-format off
  switch (function)
    {
    case FieldFunc_Min:    return field_min(field);
    case FieldFunc_Max:    return field_max(field);
    case FieldFunc_Range:  return field_range(field);
    case FieldFunc_Sum:    return field_sum(field);
    case FieldFunc_Mean:   return field_mean(field);
    case FieldFunc_Avg:    return field_avg(field);
    case FieldFunc_Std:    return field_std(field);
    case FieldFunc_Std1:   return field_std1(field);
    case FieldFunc_Var:    return field_var(field);
    case FieldFunc_Var1:   return field_var1(field);
    case FieldFunc_Meanw:  return field_meanw(field);
    case FieldFunc_Avgw:   return field_avgw(field);
    case FieldFunc_Stdw:   return field_stdw(field);
    case FieldFunc_Std1w:  return field_std1w(field);
    case FieldFunc_Varw:   return field_varw(field);
    case FieldFunc_Var1w:  return field_var1w(field);
    case FieldFunc_Skew:   return field_skew(field);
    case FieldFunc_Kurt:   return field_kurt(field);
    case FieldFunc_Median: return field_median(field);
    case FieldFunc_Count:  return field_count(field);
    default: cdo_abort("%s: function %d not implemented!", __func__, function);
    }
  // clang-format on
  return 0.0;
}
