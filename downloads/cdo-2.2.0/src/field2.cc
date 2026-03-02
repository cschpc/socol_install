/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <algorithm>
#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "field_functions.h"

// clang-format off

// arithmetic
auto arith_add   = [](const double a, const double b) { return a + b; };
auto arith_sub   = [](const double a, const double b) { return a - b; };
auto arith_mul   = [](const double a, const double b) { return a * b; };
auto arith_div   = [](const double a, const double b) { return a / b; };
auto arith_min   = [](const auto a, const auto b) { return (b > a) ? a : b; };
auto arith_max   = [](const auto a, const auto b) { return (b < a) ? a : b; };
auto arith_sumq  = [](const double a, const double b) { return a + b * b; };
auto arith_atan2 = [](const double a, const double b) { return std::atan2(a, b); };

// arithmetic with missing values
auto arith_add_mv = [](auto &a, const auto mv_a, const auto b, const auto mv_b, auto isEQ)
                    { a = (isEQ(a, mv_a) || isEQ(b, mv_b)) ? mv_a : arith_add(a, b); };
auto arith_sub_mv = [](auto &a, const auto mv_a, const auto b, const auto mv_b, auto isEQ)
                    { a = (isEQ(a, mv_a) || isEQ(b, mv_b)) ? mv_a : arith_sub(a, b); };
auto arith_mul_mv = [](auto &a, const auto mv_a, const auto b, const auto mv_b, auto isEQ)
                    { a = (isEQ(a, 0) || isEQ(b, 0)) ? 0 : ((isEQ(a, mv_a) || isEQ(b, mv_b)) ? mv_a : arith_mul(a, b)); };
auto arith_div_mv = [](auto &a, const auto mv_a, const auto b, const auto mv_b, auto isEQ)
                    { a = (isEQ(a, mv_a) || isEQ(b, mv_b) || isEQ(b, 0)) ? mv_a : arith_div(a, b); };
auto arith_min_mv = [](auto &a, const auto mv_a, const auto b, const auto mv_b, auto isEQ)
                    { a = isEQ(b, mv_b) ? a : (isEQ(a, mv_a) ? b : arith_min(a, b)); };
auto arith_max_mv = [](auto &a, const auto mv_a, const auto b, const auto mv_b, auto isEQ)
                    { a = isEQ(b, mv_b) ? a : (isEQ(a, mv_a) ? b : arith_max(a, b)); };
auto arith_sum_mv = [](auto &a, const auto mv_a, const auto b, const auto mv_b, auto isEQ)
                    { if (!isEQ(b, mv_b)) a = (isEQ(a, mv_a) ? b : arith_add(a, b)); };
auto arith_sumq_mv = [](auto &a, const auto mv_a, const auto b, const auto mv_b, auto isEQ)
                     { if (!isEQ(b, mv_b)) a = (isEQ(a, mv_a) ? arith_mul(b, b) : arith_sumq(a, b)); };
auto arith_atan2_mv = [](auto &a, const auto mv_a, const auto b, const auto mv_b, auto isEQ)
                      { a = (isEQ(a, mv_a) || isEQ(b, mv_b)) ? mv_a : arith_atan2(a, b); };
auto arith_vinit_mv = [](auto &a, const auto mv_a, const auto b, const auto mv_b, auto isEQ)
                      { a = !isEQ(b, mv_b); (void)mv_a; };
auto arith_vincr_mv = [](auto &a, const auto mv_a, const auto b, const auto mv_b, auto isEQ)
                      { if (!isEQ(b, mv_b)) a = a + 1; (void)mv_a; };
// clang-format on

template <typename T1, typename T2>
static void
varray2_div(const size_t len, Varray<T1> &v1, const Varray<T2> &v2, const double missval1)
{
  assert(len > 0);
  assert(v1.size() > 0);
  assert(v2.size() > 0);
  assert(len <= v1.size());
  assert(len <= v2.size());

#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
  for (size_t i = 0; i < len; ++i) { v1[i] = is_equal(v2[i], 0.0) ? missval1 : arith_div(v1[i], v2[i]); }
}

template <typename T1, typename T2, typename FUNC>
static void
varray2_arith(const size_t len, Varray<T1> &v1, const Varray<T2> &v2, FUNC func)
{
  assert(len > 0);
  assert(v1.size() > 0);
  assert(v2.size() > 0);
  assert(len <= v1.size());
  assert(len <= v2.size());

#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
  for (size_t i = 0; i < len; ++i) { v1[i] = func(v1[i], v2[i]); }
}

template <typename T1, typename T2, typename FUNC>
static void
varray2_arith_mv(const size_t len, Varray<T1> &v1, const Varray<T2> &v2, const double missval1, const double missval2, FUNC func)
{
  assert(len > 0);
  assert(v1.size() > 0);
  assert(v2.size() > 0);
  assert(len <= v1.size());
  assert(len <= v2.size());

  const T1 mv1 = missval1;
  const T2 mv2 = missval2;

  if (std::isnan(missval1) || std::isnan(missval2))
    {
#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
      for (size_t i = 0; i < len; ++i) func(v1[i], mv1, v2[i], mv2, dbl_is_equal);
    }
  else
    {
#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
      for (size_t i = 0; i < len; ++i) func(v1[i], mv1, v2[i], mv2, is_equal);
    }
}

// init valid values
void
field2_vinit(Field &field1, const Field &field2)
{
  if (field1.size != field2.size) cdo_abort("Fields have different size (%s)", __func__);

  if (memtype_is_float_float(field1.memType, field2.memType))
    varray2_arith_mv(field1.size, field1.vec_f, field2.vec_f, field1.missval, field2.missval, arith_vinit_mv);
  else if (memtype_is_float_double(field1.memType, field2.memType))
    varray2_arith_mv(field1.size, field1.vec_f, field2.vec_d, field1.missval, field2.missval, arith_vinit_mv);
  else if (memtype_is_double_float(field1.memType, field2.memType))
    varray2_arith_mv(field1.size, field1.vec_d, field2.vec_f, field1.missval, field2.missval, arith_vinit_mv);
  else
    varray2_arith_mv(field1.size, field1.vec_d, field2.vec_d, field1.missval, field2.missval, arith_vinit_mv);

  field1.nmiss = field2.nmiss;
}

// increment valid values
void
field2_vincr(Field &field1, const Field &field2)
{
  if (field1.size != field2.size) cdo_abort("Fields have different size (%s)", __func__);

  if (memtype_is_float_float(field1.memType, field2.memType))
    varray2_arith_mv(field1.size, field1.vec_f, field2.vec_f, field1.missval, field2.missval, arith_vincr_mv);
  else if (memtype_is_float_double(field1.memType, field2.memType))
    varray2_arith_mv(field1.size, field1.vec_f, field2.vec_d, field1.missval, field2.missval, arith_vincr_mv);
  else if (memtype_is_double_float(field1.memType, field2.memType))
    varray2_arith_mv(field1.size, field1.vec_d, field2.vec_f, field1.missval, field2.missval, arith_vincr_mv);
  else
    varray2_arith_mv(field1.size, field1.vec_d, field2.vec_d, field1.missval, field2.missval, arith_vincr_mv);

  field1.nmiss = field2.nmiss;
}

// init valid values
template <typename T1, typename T2>
static void
field2_vinit(const size_t len, Varray<T1> &v1, const Varray<T2> &v2, const T2 missval, int vinit)
{
  for (size_t i = 0; i < len; ++i) v1[i] = (DBL_IS_EQUAL(v2[i], missval)) ? 0 : vinit;
}

void
field2_vinit(Field &field1, const Field &field2, int vinit)
{
  if (field1.size != field2.size) cdo_abort("Fields have different size (%s)", __func__);

  if (memtype_is_float_float(field1.memType, field2.memType))
    field2_vinit(field1.size, field1.vec_f, field2.vec_f, (float) field1.missval, vinit);
  else if (memtype_is_float_double(field1.memType, field2.memType))
    field2_vinit(field1.size, field1.vec_f, field2.vec_d, field1.missval, vinit);
  else if (memtype_is_double_float(field1.memType, field2.memType))
    field2_vinit(field1.size, field1.vec_d, field2.vec_f, (float) field1.missval, vinit);
  else
    field2_vinit(field1.size, field1.vec_d, field2.vec_d, field1.missval, vinit);

  field1.nmiss = field2.nmiss;
}

// increment valid values
template <typename T1, typename T2>
static void
field2_vincr(const size_t len, Varray<T1> &v1, const Varray<T2> &v2, const T2 missval, int vincr)
{
  for (size_t i = 0; i < len; ++i)
    if (!DBL_IS_EQUAL(v2[i], missval)) v1[i] += vincr;
}

void
field2_vincr(Field &field1, const Field &field2, int vincr)
{
  if (field1.size != field2.size) cdo_abort("Fields have different size (%s)", __func__);

  if (memtype_is_float_float(field1.memType, field2.memType))
    field2_vincr(field1.size, field1.vec_f, field2.vec_f, (float) field1.missval, vincr);
  else if (memtype_is_float_double(field1.memType, field2.memType))
    field2_vincr(field1.size, field1.vec_f, field2.vec_d, field1.missval, vincr);
  else if (memtype_is_double_float(field1.memType, field2.memType))
    field2_vincr(field1.size, field1.vec_d, field2.vec_f, (float) field1.missval, vincr);
  else
    field2_vincr(field1.size, field1.vec_d, field2.vec_d, field1.missval, vincr);

  field1.nmiss = field2.nmiss;
}

void
field2_add(Field &field1, const Field &field2)
{
  if (field1.size != field2.size) cdo_abort("Fields have different size (%s)", __func__);

  if (field1.nmiss || field2.nmiss)
    {
      if (memtype_is_float_float(field1.memType, field2.memType))
        varray2_arith_mv(field1.size, field1.vec_f, field2.vec_f, field1.missval, field2.missval, arith_add_mv);
      else if (memtype_is_float_double(field1.memType, field2.memType))
        varray2_arith_mv(field1.size, field1.vec_f, field2.vec_d, field1.missval, field2.missval, arith_add_mv);
      else if (memtype_is_double_float(field1.memType, field2.memType))
        varray2_arith_mv(field1.size, field1.vec_d, field2.vec_f, field1.missval, field2.missval, arith_add_mv);
      else
        varray2_arith_mv(field1.size, field1.vec_d, field2.vec_d, field1.missval, field2.missval, arith_add_mv);

      field_num_mv(field1);
    }
  else
    {
      if (memtype_is_float_float(field1.memType, field2.memType))
        varray2_arith(field1.size, field1.vec_f, field2.vec_f, arith_add);
      else if (memtype_is_float_double(field1.memType, field2.memType))
        varray2_arith(field1.size, field1.vec_f, field2.vec_d, arith_add);
      else if (memtype_is_double_float(field1.memType, field2.memType))
        varray2_arith(field1.size, field1.vec_d, field2.vec_f, arith_add);
      else
        varray2_arith(field1.size, field1.vec_d, field2.vec_d, arith_add);
    }
}

void
field2_sum(Field &field1, const Field &field2)
{
  if (field1.size != field2.size) cdo_abort("Fields have different size (%s)", __func__);

  if (field1.nmiss || field2.nmiss)
    {
      if (memtype_is_float_float(field1.memType, field2.memType))
        varray2_arith_mv(field1.size, field1.vec_f, field2.vec_f, field1.missval, field2.missval, arith_sum_mv);
      else if (memtype_is_float_double(field1.memType, field2.memType))
        varray2_arith_mv(field1.size, field1.vec_f, field2.vec_d, field1.missval, field2.missval, arith_sum_mv);
      else if (memtype_is_double_float(field1.memType, field2.memType))
        varray2_arith_mv(field1.size, field1.vec_d, field2.vec_f, field1.missval, field2.missval, arith_sum_mv);
      else
        varray2_arith_mv(field1.size, field1.vec_d, field2.vec_d, field1.missval, field2.missval, arith_sum_mv);

      field_num_mv(field1);
    }
  else
    {
      if (memtype_is_float_float(field1.memType, field2.memType))
        varray2_arith(field1.size, field1.vec_f, field2.vec_f, arith_add);
      else if (memtype_is_float_double(field1.memType, field2.memType))
        varray2_arith(field1.size, field1.vec_f, field2.vec_f, arith_add);
      else if (memtype_is_double_float(field1.memType, field2.memType))
        varray2_arith(field1.size, field1.vec_d, field2.vec_f, arith_add);
      else
        varray2_arith(field1.size, field1.vec_d, field2.vec_d, arith_add);
    }
}

void
field2_sumw(Field &field1, const Field &field2, double w)
{
  const auto missval1 = field1.missval;
  const auto missval2 = field2.missval;
  auto &array1 = field1.vec_d;
  const auto &array2 = field2.vec_d;

  auto len = field1.size;
  if (len != field2.size) cdo_abort("Fields have different size (%s)", __func__);

  if (field1.nmiss || field2.nmiss)
    {
#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
      for (size_t i = 0; i < len; ++i)
        if (!DBL_IS_EQUAL(array2[i], missval2))
          {
            if (!DBL_IS_EQUAL(array1[i], missval1))
              array1[i] += w * array2[i];
            else
              array1[i] = w * array2[i];
          }

      field1.nmiss = field_num_miss(field1);
    }
  else
    {
#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
      for (size_t i = 0; i < len; ++i) array1[i] += w * array2[i];
    }
}

/*
 * Compute the occurrence of values in field, if they do not equal refval.
 * This can be used to compute the lengths of multiple periods in a timeseries.
 * Missing field values are handled like refval, i.e. they stop a running period.
 * If there is missing data in the occurence field, missing fields values do not
 * change anything (they do not start a non-period by setting occurrence to zero).
 */
void
field2_sumtr(Field &occur, const Field &field, const double refval)
{
  auto omissval = occur.missval;
  auto fmissval = field.missval;
  auto &oarray = occur.vec_d;
  const auto &farray = field.vec_d;

  auto len = occur.size;
  if (len != field.size) cdo_abort("Fields have different size (%s)", __func__);

  if (occur.nmiss || field.nmiss)
    {
#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
      for (size_t i = 0; i < len; ++i)
        if (!DBL_IS_EQUAL(farray[i], fmissval))
          {
            if (!DBL_IS_EQUAL(oarray[i], omissval))
              oarray[i] = (DBL_IS_EQUAL(farray[i], refval)) ? 0.0 : oarray[i] + 1.0;
            else
              oarray[i] = (DBL_IS_EQUAL(farray[i], refval)) ? 0.0 : 1.0;
          }
        else
          {
            if (!DBL_IS_EQUAL(oarray[i], omissval)) oarray[i] = 0.0;
          }

      occur.nmiss = field_num_miss(occur);
    }
  else
    {
#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
      for (size_t i = 0; i < len; ++i) oarray[i] = (DBL_IS_EQUAL(farray[i], refval)) ? 0.0 : oarray[i] + 1.0;
    }
}

void
field2_sumq(Field &field1, const Field &field2)
{
  if (field1.size != field2.size) cdo_abort("Fields have different size (%s)", __func__);

  if (field1.nmiss || field2.nmiss)
    {
      if (memtype_is_float_float(field1.memType, field2.memType))
        varray2_arith_mv(field1.size, field1.vec_f, field2.vec_f, field1.missval, field2.missval, arith_sumq_mv);
      else if (memtype_is_float_double(field1.memType, field2.memType))
        varray2_arith_mv(field1.size, field1.vec_f, field2.vec_d, field1.missval, field2.missval, arith_sumq_mv);
      else if (memtype_is_double_float(field1.memType, field2.memType))
        varray2_arith_mv(field1.size, field1.vec_d, field2.vec_f, field1.missval, field2.missval, arith_sumq_mv);
      else
        varray2_arith_mv(field1.size, field1.vec_d, field2.vec_d, field1.missval, field2.missval, arith_sumq_mv);

      field_num_mv(field1);
    }
  else
    {
      if (memtype_is_float_float(field1.memType, field2.memType))
        varray2_arith(field1.size, field1.vec_f, field2.vec_f, arith_sumq);
      else if (memtype_is_float_double(field1.memType, field2.memType))
        varray2_arith(field1.size, field1.vec_f, field2.vec_d, arith_sumq);
      else if (memtype_is_double_float(field1.memType, field2.memType))
        varray2_arith(field1.size, field1.vec_d, field2.vec_f, arith_sumq);
      else
        varray2_arith(field1.size, field1.vec_d, field2.vec_d, arith_sumq);
    }
}

void
field2_sumsumq(Field &field1, Field &field2, const Field &field3)
{
  field2_sumq(field2, field3);
  field2_sum(field1, field3);
}

void
field2_sumqw(Field &field1, const Field &field2, const double w)
{
  const auto missval1 = field1.missval;
  const auto missval2 = field2.missval;
  auto &array1 = field1.vec_d;
  const auto &array2 = field2.vec_d;

  const auto len = field1.size;
  if (len != field2.size) cdo_abort("Fields have different size (%s)", __func__);

  if (field1.nmiss || field2.nmiss)
    {
#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
      for (size_t i = 0; i < len; ++i)
        if (!DBL_IS_EQUAL(array2[i], missval2))
          {
            if (!DBL_IS_EQUAL(array1[i], missval1))
              array1[i] += w * array2[i] * array2[i];
            else
              array1[i] = w * array2[i] * array2[i];
          }

      field1.nmiss = field_num_miss(field1);
    }
  else
    {
#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
      for (size_t i = 0; i < len; ++i) array1[i] += w * array2[i] * array2[i];
    }
}

void
field2_sub(Field &field1, const Field &field2)
{
  if (field1.size != field2.size) cdo_abort("Fields have different size (%s)", __func__);

  if (field1.nmiss || field2.nmiss)
    {
      if (memtype_is_float_float(field1.memType, field2.memType))
        varray2_arith_mv(field1.size, field1.vec_f, field2.vec_f, field1.missval, field2.missval, arith_sub_mv);
      else if (memtype_is_float_double(field1.memType, field2.memType))
        varray2_arith_mv(field1.size, field1.vec_f, field2.vec_d, field1.missval, field2.missval, arith_sub_mv);
      else if (memtype_is_double_float(field1.memType, field2.memType))
        varray2_arith_mv(field1.size, field1.vec_d, field2.vec_f, field1.missval, field2.missval, arith_sub_mv);
      else
        varray2_arith_mv(field1.size, field1.vec_d, field2.vec_d, field1.missval, field2.missval, arith_sub_mv);

      field_num_mv(field1);
    }
  else
    {
      if (memtype_is_float_float(field1.memType, field2.memType))
        varray2_arith(field1.size, field1.vec_f, field2.vec_f, arith_sub);
      else if (memtype_is_float_double(field1.memType, field2.memType))
        varray2_arith(field1.size, field1.vec_f, field2.vec_d, arith_sub);
      else if (memtype_is_double_float(field1.memType, field2.memType))
        varray2_arith(field1.size, field1.vec_d, field2.vec_f, arith_sub);
      else
        varray2_arith(field1.size, field1.vec_d, field2.vec_d, arith_sub);
    }
}

void
field2_mul(Field &field1, const Field &field2)
{
  if (field1.size != field2.size) cdo_abort("Fields have different size (%s)", __func__);

  if (field1.nmiss || field2.nmiss)
    {
      if (memtype_is_float_float(field1.memType, field2.memType))
        varray2_arith_mv(field1.size, field1.vec_f, field2.vec_f, field1.missval, field2.missval, arith_mul_mv);
      else if (memtype_is_float_double(field1.memType, field2.memType))
        varray2_arith_mv(field1.size, field1.vec_f, field2.vec_d, field1.missval, field2.missval, arith_mul_mv);
      else if (memtype_is_double_float(field1.memType, field2.memType))
        varray2_arith_mv(field1.size, field1.vec_d, field2.vec_f, field1.missval, field2.missval, arith_mul_mv);
      else
        varray2_arith_mv(field1.size, field1.vec_d, field2.vec_d, field1.missval, field2.missval, arith_mul_mv);

      field_num_mv(field1);
    }
  else
    {
      if (memtype_is_float_float(field1.memType, field2.memType))
        varray2_arith(field1.size, field1.vec_f, field2.vec_f, arith_mul);
      else if (memtype_is_float_double(field1.memType, field2.memType))
        varray2_arith(field1.size, field1.vec_f, field2.vec_d, arith_mul);
      else if (memtype_is_double_float(field1.memType, field2.memType))
        varray2_arith(field1.size, field1.vec_d, field2.vec_f, arith_mul);
      else
        varray2_arith(field1.size, field1.vec_d, field2.vec_d, arith_mul);
    }
}

void
field2_div(Field &field1, const Field &field2)
{
  if (field1.size != field2.size) cdo_abort("Fields have different size (%s)", __func__);

  if (field1.nmiss || field2.nmiss)
    {
      if (memtype_is_float_float(field1.memType, field2.memType))
        varray2_arith_mv(field1.size, field1.vec_f, field2.vec_f, field1.missval, field2.missval, arith_div_mv);
      else if (memtype_is_float_double(field1.memType, field2.memType))
        varray2_arith_mv(field1.size, field1.vec_f, field2.vec_d, field1.missval, field2.missval, arith_div_mv);
      else if (memtype_is_double_float(field1.memType, field2.memType))
        varray2_arith_mv(field1.size, field1.vec_d, field2.vec_f, field1.missval, field2.missval, arith_div_mv);
      else
        varray2_arith_mv(field1.size, field1.vec_d, field2.vec_d, field1.missval, field2.missval, arith_div_mv);

      field_num_mv(field1);
    }
  else
    {
      if (memtype_is_float_float(field1.memType, field2.memType))
        varray2_div(field1.size, field1.vec_f, field2.vec_f, field1.missval);
      else if (memtype_is_float_double(field1.memType, field2.memType))
        varray2_div(field1.size, field1.vec_f, field2.vec_d, field1.missval);
      else if (memtype_is_double_float(field1.memType, field2.memType))
        varray2_div(field1.size, field1.vec_d, field2.vec_f, field1.missval);
      else
        varray2_div(field1.size, field1.vec_d, field2.vec_d, field1.missval);

      field_num_mv(field1);
    }
}

void
field2_atan2(Field &field1, const Field &field2)
{
  if (field1.size != field2.size) cdo_abort("Fields have different size (%s)", __func__);

  if (memtype_is_float_float(field1.memType, field2.memType))
    varray2_arith_mv(field1.size, field1.vec_f, field2.vec_f, field1.missval, field2.missval, arith_atan2_mv);
  else if (memtype_is_float_double(field1.memType, field2.memType))
    varray2_arith_mv(field1.size, field1.vec_f, field2.vec_d, field1.missval, field2.missval, arith_atan2_mv);
  else if (memtype_is_double_float(field1.memType, field2.memType))
    varray2_arith_mv(field1.size, field1.vec_d, field2.vec_f, field1.missval, field2.missval, arith_atan2_mv);
  else
    varray2_arith_mv(field1.size, field1.vec_d, field2.vec_d, field1.missval, field2.missval, arith_atan2_mv);

  field_num_mv(field1);
}

void
field2_set_miss(Field &field1, const Field &field2)
{
  const auto missval1 = field1.missval;
  auto &array1 = field1.vec_d;
  const auto &array2 = field2.vec_d;

  const auto len = field1.size;
  if (len != field2.size) cdo_abort("Fields have different size (%s)", __func__);

  if (field1.nmiss)
    {
      for (size_t i = 0; i < len; ++i) array1[i] = DBL_IS_EQUAL(array1[i], missval1) ? array2[i] : array1[i];

      field1.nmiss = field_num_miss(field1);
    }
}

void
field2_min(Field &field1, const Field &field2)
{
  if (field1.size != field2.size) cdo_abort("Fields have different size (%s)", __func__);

  if (field1.nmiss || field2.nmiss)
    {
      if (memtype_is_float_float(field1.memType, field2.memType))
        varray2_arith_mv(field1.size, field1.vec_f, field2.vec_f, field1.missval, field2.missval, arith_min_mv);
      else if (memtype_is_float_double(field1.memType, field2.memType))
        varray2_arith_mv(field1.size, field1.vec_f, field2.vec_d, field1.missval, field2.missval, arith_min_mv);
      else if (memtype_is_double_float(field1.memType, field2.memType))
        varray2_arith_mv(field1.size, field1.vec_d, field2.vec_f, field1.missval, field2.missval, arith_min_mv);
      else
        varray2_arith_mv(field1.size, field1.vec_d, field2.vec_d, field1.missval, field2.missval, arith_min_mv);

      field_num_mv(field1);
    }
  else
    {
      if (memtype_is_float_float(field1.memType, field2.memType))
        varray2_arith(field1.size, field1.vec_f, field2.vec_f, arith_min);
      else if (memtype_is_float_double(field1.memType, field2.memType))
        varray2_arith(field1.size, field1.vec_f, field2.vec_d, arith_min);
      else if (memtype_is_double_float(field1.memType, field2.memType))
        varray2_arith(field1.size, field1.vec_d, field2.vec_f, arith_min);
      else
        varray2_arith(field1.size, field1.vec_d, field2.vec_d, arith_min);
    }
}

void
field2_max(Field &field1, const Field &field2)
{
  if (field1.size != field2.size) cdo_abort("Fields have different size (%s)", __func__);

  if (field1.nmiss || field2.nmiss)
    {
      if (memtype_is_float_float(field1.memType, field2.memType))
        varray2_arith_mv(field1.size, field1.vec_f, field2.vec_f, field1.missval, field2.missval, arith_max_mv);
      else if (memtype_is_float_double(field1.memType, field2.memType))
        varray2_arith_mv(field1.size, field1.vec_f, field2.vec_d, field1.missval, field2.missval, arith_max_mv);
      else if (memtype_is_double_float(field1.memType, field2.memType))
        varray2_arith_mv(field1.size, field1.vec_d, field2.vec_f, field1.missval, field2.missval, arith_max_mv);
      else
        varray2_arith_mv(field1.size, field1.vec_d, field2.vec_d, field1.missval, field2.missval, arith_max_mv);

      field_num_mv(field1);
    }
  else
    {
      if (memtype_is_float_float(field1.memType, field2.memType))
        varray2_arith(field1.size, field1.vec_f, field2.vec_f, arith_max);
      else if (memtype_is_float_double(field1.memType, field2.memType))
        varray2_arith(field1.size, field1.vec_f, field2.vec_d, arith_max);
      else if (memtype_is_double_float(field1.memType, field2.memType))
        varray2_arith(field1.size, field1.vec_d, field2.vec_f, arith_max);
      else
        varray2_arith(field1.size, field1.vec_d, field2.vec_d, arith_max);
    }
}

void
field2_maxmin(Field &field1, Field &field2, const Field &field3)
{
  field2_min(field2, field3);
  field2_max(field1, field3);
}

auto field2_set_index = [](auto &v1, auto &v2, const auto v3, const auto idx) {
  v2 = v3;
  v1 = idx;
};

template <typename T>
void
field2_minidx(size_t nmiss3, size_t len, double missval3, Field &field1, Field &field2, const Varray<T> &v3, int idx)
{
  const auto missval1 = field1.missval;
  const auto missval2 = field2.missval;
  auto &v1 = field1.vec_d;
  auto &v2 = field2.vec_d;

  if (field2.nmiss || nmiss3)
    {
      for (size_t i = 0; i < len; ++i)
        {
          if (DBL_IS_EQUAL(v3[i], missval3))
            {
              if (DBL_IS_EQUAL(v2[i], missval2)) v1[i] = missval1;
            }
          else if (v3[i] < v2[i] || DBL_IS_EQUAL(v2[i], missval2)) { field2_set_index(v1[i], v2[i], v3[i], idx); }
        }

      field_num_mv(field1);
      field_num_mv(field2);
    }
  else
    {
      for (size_t i = 0; i < len; ++i)
        {
          if (v3[i] < v2[i]) field2_set_index(v1[i], v2[i], v3[i], idx);
        }
    }
}

void
field2_minidx(Field &field1, Field &field2, const Field &field3, int idx)
{
  if (field1.size != field3.size) cdo_abort("Fields have different size (%s)", __func__);
  if (field2.size != field3.size) cdo_abort("Fields have different size (%s)", __func__);

  if (field3.memType == MemType::Float)
    field2_minidx(field3.nmiss, field3.size, field3.missval, field1, field2, field3.vec_f, idx);
  else
    field2_minidx(field3.nmiss, field3.size, field3.missval, field1, field2, field3.vec_d, idx);
}

template <typename T>
void
field2_maxidx(size_t nmiss3, size_t len, double missval3, Field &field1, Field &field2, const Varray<T> &v3, int idx)
{
  const auto missval1 = field1.missval;
  const auto missval2 = field2.missval;
  auto &v1 = field1.vec_d;
  auto &v2 = field2.vec_d;

  if (field2.nmiss || nmiss3)
    {
      for (size_t i = 0; i < len; ++i)
        {
          if (DBL_IS_EQUAL(v3[i], missval3))
            {
              if (DBL_IS_EQUAL(v2[i], missval2)) v1[i] = missval1;
            }
          else if (v3[i] > v2[i] || DBL_IS_EQUAL(v2[i], missval2)) { field2_set_index(v1[i], v2[i], v3[i], idx); }
        }

      field_num_mv(field1);
      field_num_mv(field2);
    }
  else
    {
      for (size_t i = 0; i < len; ++i)
        {
          if (v3[i] > v2[i]) field2_set_index(v1[i], v2[i], v3[i], idx);
        }
    }
}

void
field2_maxidx(Field &field1, Field &field2, const Field &field3, int idx)
{
  if (field1.size != field3.size) cdo_abort("Fields have different size (%s)", __func__);
  if (field2.size != field3.size) cdo_abort("Fields have different size (%s)", __func__);

  if (field3.memType == MemType::Float)
    field2_maxidx(field3.nmiss, field3.size, field3.missval, field1, field2, field3.vec_f, idx);
  else
    field2_maxidx(field3.nmiss, field3.size, field3.missval, field1, field2, field3.vec_d, idx);
}

static size_t
field_set_nmiss(const size_t len, Varray<double> &v, double missval)
{
  size_t nmiss = 0;

  if (std::isnan(missval))
    {
      for (size_t i = 0; i < len; ++i)
        if (DBL_IS_EQUAL(v[i], missval) || v[i] < 0)
          {
            v[i] = missval;
            nmiss++;
          }
    }
  else
    {
      for (size_t i = 0; i < len; ++i)
        if (IS_EQUAL(v[i], missval) || v[i] < 0)
          {
            v[i] = missval;
            nmiss++;
          }
    }

  return nmiss;
}

void
field2_var(Field &field1, const Field &field2, const Field &field3, const int divisor)
{
  const auto missval1 = field1.missval;
  const auto missval2 = field2.missval;
  auto &array1 = field1.vec_d;
  const auto &array2 = field2.vec_d;
  const auto &array3 = field3.vec_d;

  const auto len = field1.size;
  if (len != field2.size) cdo_abort("Fields have different size (%s)", __func__);

  if (field1.nmiss || field2.nmiss)
    {
      for (size_t i = 0; i < len; ++i)
        {
          const auto temp = DIVMN(MULMN(array1[i], array1[i]), array3[i]);
          array1[i] = DIVMN(SUBMN(array2[i], temp), array3[i] - divisor);
          if (array1[i] < 0 && array1[i] > -1.e-5) array1[i] = 0;
        }
    }
  else
    {
      for (size_t i = 0; i < len; ++i)
        {
          const auto temp = DIV(MUL(array1[i], array1[i]), array3[i]);
          array1[i] = DIV(SUB(array2[i], temp), array3[i] - divisor);
          if (array1[i] < 0 && array1[i] > -1.e-5) array1[i] = 0;
        }
    }

  field1.nmiss = field_set_nmiss(len, array1, missval1);
}

void
field2_std(Field &field1, const Field &field2, const Field &field3, const int divisor)
{
  const auto missval1 = field1.missval;
  auto &array1 = field1.vec_d;

  const auto len = field1.size;
  if (len != field2.size) cdo_abort("Fields have different size (%s)", __func__);

  field2_var(field1, field2, field3, divisor);

  size_t nmiss = 0;
  for (size_t i = 0; i < len; ++i)
    if (DBL_IS_EQUAL(array1[i], missval1) || array1[i] < 0)
      {
        array1[i] = missval1;
        nmiss++;
      }
    else { array1[i] = IS_NOT_EQUAL(array1[i], 0) ? std::sqrt(array1[i]) : 0; }
  field1.nmiss = nmiss;
}

void
fieldc_var(Field &field1, const Field &field2, const int numSets, const int divisor)
{
  const auto nsetx = numSets - divisor;
  const auto missval1 = field1.missval;
  const auto missval2 = field2.missval;
  auto &array1 = field1.vec_d;
  const auto &array2 = field2.vec_d;

  const auto len = field1.size;
  if (len != field2.size) cdo_abort("Fields have different size (%s)", __func__);

  if (nsetx == 0)
    {
      for (size_t i = 0; i < len; ++i) array1[i] = missval1;
    }
  else if (field1.nmiss || field2.nmiss)
    {
      for (size_t i = 0; i < len; ++i)
        {
          const auto temp = MULMN(array1[i], array1[i]) / numSets;
          array1[i] = SUBMN(array2[i], temp) / nsetx;
          if (array1[i] < 0 && array1[i] > -1.e-5) array1[i] = 0;
        }
    }
  else
    {
      for (size_t i = 0; i < len; ++i)
        {
          const auto temp = MUL(array1[i], array1[i]) / numSets;
          array1[i] = SUB(array2[i], temp) / nsetx;
          if (array1[i] < 0 && array1[i] > -1.e-5) array1[i] = 0;
        }
    }

  field1.nmiss = field_set_nmiss(len, array1, missval1);
}

void
fieldc_std(Field &field1, const Field &field2, const int numSets, const int divisor)
{
  const auto missval1 = field1.missval;
  auto &array1 = field1.vec_d;

  const auto len = field1.size;
  if (len != field2.size) cdo_abort("Fields have different size (%s)", __func__);

  fieldc_var(field1, field2, numSets, divisor);

  size_t nmiss = 0;
  for (size_t i = 0; i < len; ++i)
    if (DBL_IS_EQUAL(array1[i], missval1) || array1[i] < 0)
      {
        array1[i] = missval1;
        nmiss++;
      }
    else { array1[i] = IS_NOT_EQUAL(array1[i], 0) ? std::sqrt(array1[i]) : 0; }

  field1.nmiss = nmiss;
}

void
field2_moq(Field &field1, const Field &field2)
{
  const auto missval1 = field1.missval;
  const auto missval2 = field2.missval;
  auto &array1 = field1.vec_d;
  const auto &array2 = field2.vec_d;

  const auto len = field1.size;
  if (len != field2.size) cdo_abort("Fields have different size (%s)", __func__);

  if (field2.nmiss)
    {
      for (size_t i = 0; i < len; ++i)
        if (!DBL_IS_EQUAL(array2[i], missval2))
          array1[i] = array2[i] * array2[i];
        else
          array1[i] = missval1;

      field1.nmiss = field_num_miss(field1);
    }
  else
    {
      for (size_t i = 0; i < len; ++i) array1[i] = array2[i] * array2[i];
    }
}

void
field2_moqw(Field &field1, const Field &field2, const double w)
{
  const auto missval1 = field1.missval;
  const auto missval2 = field2.missval;
  auto &array1 = field1.vec_d;
  const auto &array2 = field2.vec_d;

  const auto len = field1.size;
  if (len != field2.size) cdo_abort("Fields have different size (%s)", __func__);

  if (field2.nmiss)
    {
      for (size_t i = 0; i < len; ++i)
        if (!DBL_IS_EQUAL(array2[i], missval2))
          array1[i] = w * array2[i] * array2[i];
        else
          array1[i] = missval1;

      field1.nmiss = field_num_miss(field1);
    }
  else
    {
      for (size_t i = 0; i < len; ++i) array1[i] = w * array2[i] * array2[i];
    }
}

/**
 * Counts the number of nonmissing values. The result of the operation is computed according to the following rules:
 *
 * field1  field2  result
 * a       b       a + 1
 * a       miss    a
 * miss    b       1
 * miss    miss    miss
 *
 * @param field1 the 1st input field, also holds the result
 * @param field2 the 2nd input field
 */
template <typename T1, typename T2>
static void
varray2_count_mv(const size_t len, Varray<T1> &v1, const Varray<T2> &v2, const double mv1, const double mv2)
{
  assert(len > 0);
  assert(v1.size() > 0);
  assert(v2.size() > 0);
  assert(len <= v1.size());
  assert(len <= v2.size());

  for (size_t i = 0; i < len; ++i)
    if (!DBL_IS_EQUAL(v2[i], mv2))
      {
        if (!DBL_IS_EQUAL(v1[i], mv1))
          v1[i] += 1.0;
        else
          v1[i] = 1.0;
      }
}

template <typename T>
static void
varray2_count(const size_t len, Varray<T> &v)
{
  assert(len > 0);

  for (size_t i = 0; i < len; ++i) v[i] += 1.0;
}

void
field2_count(Field &field1, const Field &field2)
{
  if (field1.size != field2.size) cdo_abort("Fields have different size (%s)", __func__);

  if (field1.nmiss || field2.nmiss)
    {
      if (memtype_is_float_float(field1.memType, field2.memType))
        varray2_count_mv(field1.size, field1.vec_f, field2.vec_f, field1.missval, field2.missval);
      else if (memtype_is_float_double(field1.memType, field2.memType))
        varray2_count_mv(field1.size, field1.vec_f, field2.vec_d, field1.missval, field2.missval);
      else if (memtype_is_double_float(field1.memType, field2.memType))
        varray2_count_mv(field1.size, field1.vec_d, field2.vec_f, field1.missval, field2.missval);
      else
        varray2_count_mv(field1.size, field1.vec_d, field2.vec_d, field1.missval, field2.missval);

      field1.nmiss = field_num_miss(field1);
    }
  else
    {
      if (field1.memType == MemType::Float)
        varray2_count(field1.size, field1.vec_f);
      else
        varray2_count(field1.size, field1.vec_d);
    }
}

void
field2_function(Field &field1, const Field &field2, int function)
{
  // clang-format off
  switch (function)
    {
    case FieldFunc_Add:     field2_add(field1, field2);   break;
    case FieldFunc_Min:     field2_min(field1, field2);   break;
    case FieldFunc_Max:     field2_max(field1, field2);   break;
    case FieldFunc_Sum:     field2_sum(field1, field2);   break;
    case FieldFunc_Mean:    field2_sum(field1, field2);   break;
    case FieldFunc_Avg:     field2_add(field1, field2);   break;
    case FieldFunc_Sub:     field2_sub(field1, field2);   break;
    case FieldFunc_Mul:     field2_mul(field1, field2);   break;
    case FieldFunc_Div:     field2_div(field1, field2);   break;
    case FieldFunc_Atan2:   field2_atan2(field1, field2); break;
    case FieldFunc_Setmiss: field2_set_miss(field1, field2); break;
    default: cdo_abort("%s: function %d not implemented!", __func__, function);
    }
  // clang-format on
}
