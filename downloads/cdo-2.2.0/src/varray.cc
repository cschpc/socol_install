/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cfloat>
#include <cfenv>
#include <cassert>
#include <limits>

#include "compare.h"
#include "varray.h"
#include "cimdOmp.h"

//#pragma STDC FENV_ACCESS ON

const char *
fpe_errstr(int fpeRaised)
{
  const char *errstr = nullptr;

  // clang-format off
  if      (fpeRaised & FE_DIVBYZERO) errstr = "division by zero";
  else if (fpeRaised & FE_INEXACT)   errstr = "inexact result";
  else if (fpeRaised & FE_INVALID)   errstr = "invalid result";
  else if (fpeRaised & FE_OVERFLOW)  errstr = "overflow";
  else if (fpeRaised & FE_UNDERFLOW) errstr = "underflow";
  // clang-format on

  return errstr;
}

template <typename T>
inline T
min_value(T v1, T v2)
{
  return (v1 < v2) ? v1 : v2;
}

template <typename T>
inline T
max_value(T v1, T v2)
{
  return (v1 > v2) ? v1 : v2;
}

template <typename T>
MinMax
varray_min_max_mv(size_t len, const T *array, T missval)
{
  auto f_minmax_mv = [](auto a, auto mv_a, auto &vmin, auto &vmax, auto &nvals, auto is_EQ) {
    if (!is_EQ(a, mv_a))
      {
        vmin = min_value(vmin, a);
        vmax = max_value(vmax, a);
        nvals++;
      }
  };

  T vmin = std::numeric_limits<T>::max();
  T vmax = -std::numeric_limits<T>::max();

  size_t nvals = 0;
  if (std::isnan(missval))
    for (size_t i = 0; i < len; ++i) f_minmax_mv(array[i], missval, vmin, vmax, nvals, dbl_is_equal);
  else
    for (size_t i = 0; i < len; ++i) f_minmax_mv(array[i], missval, vmin, vmax, nvals, is_equal);

  return MinMax(vmin, vmax, nvals);
}

// Explicit instantiation
template MinMax varray_min_max_mv(size_t len, const float *array, float missval);
template MinMax varray_min_max_mv(size_t len, const double *array, double missval);

template <typename T>
MinMax
varray_min_max_mv(size_t len, const Varray<T> &v, T missval)
{
  return varray_min_max_mv(len, v.data(), missval);
}

// Explicit instantiation
template MinMax varray_min_max_mv(size_t len, const Varray<float> &v, float missval);
template MinMax varray_min_max_mv(size_t len, const Varray<double> &v, double missval);

template <typename T>
MinMaxSum
varray_min_max_sum(size_t len, const Varray<T> &v, MinMaxSum mms)
{
  auto f_minmaxsum = [](auto val, auto &vmin, auto &vmax, auto &vsum) {
    vmin = min_value(vmin, val);
    vmax = max_value(vmax, val);
    vsum += val;
  };

  auto vmin = mms.min;
  auto vmax = mms.max;
  auto vsum = mms.sum;

#ifndef __ICC  // wrong result with icc19
#ifdef HAVE_OPENMP4
#pragma omp parallel for simd if (len > cdoMinLoopSize) default(shared) schedule(static) reduction(min : vmin) reduction(max : vmax) reduction(+ : vsum)
#endif
#endif
  for (size_t i = 0; i < len; ++i) f_minmaxsum((double) v[i], vmin, vmax, vsum);

  return MinMaxSum(vmin, vmax, vsum, len);
}

// Explicit instantiation
template MinMaxSum varray_min_max_sum(size_t len, const Varray<float> &v, MinMaxSum mms);
template MinMaxSum varray_min_max_sum(size_t len, const Varray<double> &v, MinMaxSum mms);

template <typename T>
MinMaxSum
varray_min_max_sum_mv(size_t len, const Varray<T> &v, T missval, MinMaxSum mms)
{
  auto f_minmaxsum_mv = [](auto val, auto mv, auto &vmin, auto &vmax, auto &vsum, auto &nvals, auto is_EQ) {
    if (!is_EQ(val, mv))
      {
        vmin = min_value(vmin, val);
        vmax = max_value(vmax, val);
        vsum += val;
        nvals++;
      }
  };

  auto vmin = mms.min;
  auto vmax = mms.max;
  auto vsum = mms.sum;

  size_t nvals = 0;
  if (std::isnan(missval))
    {
#ifndef __ICC  // wrong result with icc19
#ifdef HAVE_OPENMP4
#pragma omp parallel for simd if (len > cdoMinLoopSize) default(shared) schedule(static) reduction(min : vmin) reduction(max : vmax) reduction(+ : vsum,nvals)
#endif
#endif
      for (size_t i = 0; i < len; ++i) f_minmaxsum_mv((double) v[i], missval, vmin, vmax, vsum, nvals, dbl_is_equal);
    }
  else
    {
#ifndef __ICC  // wrong result with icc19
#ifdef HAVE_OPENMP4
#pragma omp parallel for simd if (len > cdoMinLoopSize) default(shared) schedule(static) reduction(min : vmin) reduction(max : vmax) reduction(+ : vsum,nvals)
#endif
#endif
      for (size_t i = 0; i < len; ++i) f_minmaxsum_mv((double) v[i], missval, vmin, vmax, vsum, nvals, is_equal);
    }

  if (nvals == 0 && is_equal(vmin, std::numeric_limits<double>::max())) vmin = missval;
  if (nvals == 0 && is_equal(vmax, -std::numeric_limits<double>::max())) vmax = missval;

  return MinMaxSum(vmin, vmax, vsum, nvals);
}

// Explicit instantiation
template MinMaxSum varray_min_max_sum_mv(size_t len, const Varray<float> &v, float missval, MinMaxSum mms);
template MinMaxSum varray_min_max_sum_mv(size_t len, const Varray<double> &v, double missval, MinMaxSum mms);

template <typename T>
MinMaxMean
varray_min_max_mean(size_t len, const Varray<T> &v)
{
  auto mms = varray_min_max_sum(len, v, MinMaxSum());
  auto rmean = (len != 0) ? mms.sum / static_cast<double>(len) : 0.0;
  return MinMaxMean(mms.min, mms.max, rmean, len);
}

// Explicit instantiation
template MinMaxMean varray_min_max_mean(size_t len, const Varray<float> &v);
template MinMaxMean varray_min_max_mean(size_t len, const Varray<double> &v);

template <typename T>
MinMaxMean
varray_min_max_mean_mv(size_t len, const Varray<T> &v, T missval)
{
  auto mms = varray_min_max_sum_mv(len, v, missval, MinMaxSum());
  auto rmean = (mms.n != 0) ? mms.sum / static_cast<double>(mms.n) : missval;
  return MinMaxMean(mms.min, mms.max, rmean, mms.n);
}

// Explicit instantiation
template MinMaxMean varray_min_max_mean_mv(size_t len, const Varray<float> &v, float missval);
template MinMaxMean varray_min_max_mean_mv(size_t len, const Varray<double> &v, double missval);

template <typename T>
MinMax
array_min_max_mask(size_t len, const T *const array, const Varray<int> &mask)
{
  T rmin = std::numeric_limits<T>::max();
  T rmax = -std::numeric_limits<T>::max();

  if (!mask.empty())
    {
      for (size_t i = 0; i < len; ++i)
        {
          if (mask[i] == 0)
            {
              rmin = min_value(rmin, array[i]);
              rmax = max_value(rmax, array[i]);
            }
        }
    }
  else
    {
      for (size_t i = 0; i < len; ++i)
        {
          rmin = min_value(rmin, array[i]);
          rmax = max_value(rmax, array[i]);
        }
    }

  return MinMax(rmin, rmax);
}

// Explicit instantiation
template MinMax array_min_max_mask(size_t len, const float *const array, const Varray<int> &mask);
template MinMax array_min_max_mask(size_t len, const double *const array, const Varray<int> &mask);

void
array_add_array(size_t len, double *array1, const double *array2)
{
#ifdef HAVE_OPENMP4
#pragma omp simd
#endif
  for (size_t i = 0; i < len; ++i) array1[i] += array2[i];
}

void
array_add_array_mv(size_t len, double *array1, const double *array2, double missval)
{
  if (std::isnan(missval))
    {
      for (size_t i = 0; i < len; ++i)
        if (!dbl_is_equal(array2[i], missval)) array1[i] = dbl_is_equal(array1[i], missval) ? array2[i] : array1[i] + array2[i];
    }
  else
    {
      for (size_t i = 0; i < len; ++i)
        if (is_not_equal(array2[i], missval)) array1[i] = is_equal(array1[i], missval) ? array2[i] : array1[i] + array2[i];
    }
}

auto count_mv = [](auto val, auto mv, auto &num, auto is_EQ) {
  if (is_EQ(val, mv)) num++;
};

template <typename T>
size_t
array_num_mv(size_t len, const T *array, T missval)
{
  size_t nmiss = 0;

  if (std::isnan(missval))
    {
      for (size_t i = 0; i < len; ++i) count_mv(array[i], missval, nmiss, dbl_is_equal);
    }
  else
    {
      for (size_t i = 0; i < len; ++i) count_mv(array[i], missval, nmiss, is_equal);
    }

  return nmiss;
}

// Explicit instantiation
template size_t array_num_mv(size_t len, const float *array, float missval);
template size_t array_num_mv(size_t len, const double *array, double missval);

template <typename T>
size_t
varray_num_mv(size_t len, const Varray<T> &v, T missval)
{
  assert(len > 0);
  assert(v.size() > 0);
  assert(len <= v.size());

  size_t nmiss = 0;

  if (std::isnan(missval))
    {
      for (size_t i = 0; i < len; ++i) count_mv(v[i], missval, nmiss, dbl_is_equal);
    }
  else
    {
      for (size_t i = 0; i < len; ++i) count_mv(v[i], missval, nmiss, is_equal);
    }

  return nmiss;
}

// Explicit instantiation
template size_t varray_num_mv(size_t len, const Varray<float> &v, float missval);
template size_t varray_num_mv(size_t len, const Varray<double> &v, double missval);

template <typename T>
MinMax
varray_min_max(size_t len, const T *array)
{
  T vmin = std::numeric_limits<T>::max();
  T vmax = -std::numeric_limits<T>::max();

#ifndef __ICC  // wrong result with icc19
#ifdef HAVE_OPENMP4
#pragma omp parallel for simd if (len > cdoMinLoopSize) default(shared) schedule(static) reduction(min : vmin) reduction(max : vmax)
#endif
#endif
  for (size_t i = 0; i < len; ++i)
    {
      vmin = min_value(vmin, array[i]);
      vmax = max_value(vmax, array[i]);
    }

  return MinMax(vmin, vmax);
}

// Explicit instantiation
template MinMax varray_min_max(size_t len, const float *array);
template MinMax varray_min_max(size_t len, const double *array);

template <typename T>
MinMax
varray_min_max(size_t len, const Varray<T> &v)
{
  return varray_min_max(len, v.data());
}

// Explicit instantiation
template MinMax varray_min_max(size_t len, const Varray<float> &v);
template MinMax varray_min_max(size_t len, const Varray<double> &v);

template <typename T>
MinMax
varray_min_max(const Varray<T> &v)
{
  T vmin = std::numeric_limits<T>::max();
  T vmax = -std::numeric_limits<T>::max();

  auto len = v.size();
#ifndef __ICC  // wrong result with icc19
#ifdef HAVE_OPENMP4
#pragma omp parallel for simd if (len > cdoMinLoopSize) default(shared) schedule(static) reduction(min : vmin) reduction(max : vmax)
#endif
#endif
  for (size_t i = 0; i < len; ++i)
    {
      vmin = min_value(vmin, v[i]);
      vmax = max_value(vmax, v[i]);
    }

  return MinMax(vmin, vmax);
}

// Explicit instantiation
template MinMax varray_min_max(const Varray<float> &v);
template MinMax varray_min_max(const Varray<double> &v);

template <typename T>
T
varray_min(size_t len, const Varray<T> &v)
{
  assert(len > 0);
  assert(v.size() > 0);
  assert(len <= v.size());

  auto vmin = v[0];

#ifdef HAVE_OPENMP4
#pragma omp parallel for simd if (len > cdoMinLoopSize) default(shared) schedule(static) reduction(min : vmin)
#endif
  for (size_t i = 0; i < len; ++i) vmin = min_value(vmin, v[i]);

  return vmin;
}

// Explicit instantiation
template float varray_min(size_t len, const Varray<float> &v);
template double varray_min(size_t len, const Varray<double> &v);

template <typename T>
T
varray_max(size_t len, const Varray<T> &v)
{
  assert(len > 0);
  assert(v.size() > 0);
  assert(len <= v.size());

  auto vmax = v[0];

#ifdef HAVE_OPENMP4
#pragma omp parallel for simd if (len > cdoMinLoopSize) default(shared) schedule(static) reduction(max : vmax)
#endif
  for (size_t i = 0; i < len; ++i) vmax = max_value(vmax, v[i]);

  return vmax;
}

// Explicit instantiation
template float varray_max(size_t len, const Varray<float> &v);
template double varray_max(size_t len, const Varray<double> &v);

template <typename T>
T
varray_range(size_t len, const Varray<T> &v)
{
  assert(len > 0);
  assert(v.size() > 0);
  assert(len <= v.size());

  auto vmin = v[0];
  auto vmax = v[0];

#ifdef HAVE_OPENMP4
#pragma omp parallel for simd if (len > cdoMinLoopSize) default(shared) schedule(static) reduction(min : vmin) reduction(max : vmax)
#endif
  for (size_t i = 0; i < len; ++i)
    {
      vmin = min_value(vmin, v[i]);
      vmax = max_value(vmax, v[i]);
    }

  return (vmax - vmin);
}

// Explicit instantiation
template float varray_range(size_t len, const Varray<float> &v);
template double varray_range(size_t len, const Varray<double> &v);

template <typename T>
T
varray_min_mv(size_t len, const Varray<T> &v, T missval)
{
  auto f_min_mv = [](auto a, auto mv_a, auto &vmin, auto is_EQ) {
    if (!is_EQ(a, mv_a)) vmin = min_value(vmin, a);
  };

  assert(len > 0);
  assert(v.size() > 0);
  assert(len <= v.size());

  T vmin = std::numeric_limits<T>::max();

  if (std::isnan(missval))
    {
#ifndef __ICC  // internal error with icc22: lambda not supported
#ifdef HAVE_OPENMP4
#pragma omp parallel for if (len > cdoMinLoopSize) default(shared) schedule(static) reduction(min : vmin)
#endif
#endif
      for (size_t i = 0; i < len; ++i) f_min_mv(v[i], missval, vmin, dbl_is_equal);
    }
  else
    {
#ifndef __ICC  // internal error with icc22: lambda not supported
#ifdef HAVE_OPENMP4
#pragma omp parallel for if (len > cdoMinLoopSize) default(shared) schedule(static) reduction(min : vmin)
#endif
#endif
      for (size_t i = 0; i < len; ++i) f_min_mv(v[i], missval, vmin, is_equal);
    }

  if (is_equal(vmin, std::numeric_limits<T>::max())) vmin = missval;

  return vmin;
}

// Explicit instantiation
template float varray_min_mv(size_t len, const Varray<float> &v, float missval);
template double varray_min_mv(size_t len, const Varray<double> &v, double missval);

template <typename T>
T
varray_max_mv(size_t len, const Varray<T> &v, T missval)
{
  auto f_max_mv = [](auto a, auto mv_a, auto &vmax, auto is_EQ) {
    if (!is_EQ(a, mv_a)) vmax = max_value(vmax, a);
  };

  assert(len > 0);
  assert(v.size() > 0);
  assert(len <= v.size());

  T vmax = -std::numeric_limits<T>::max();

  if (std::isnan(missval))
    {
#ifndef __ICC  // internal error with icc22: lambda not supported
#ifdef HAVE_OPENMP4
#pragma omp parallel for if (len > cdoMinLoopSize) default(shared) schedule(static) reduction(max : vmax)
#endif
#endif
      for (size_t i = 0; i < len; ++i) f_max_mv(v[i], missval, vmax, dbl_is_equal);
    }
  else
    {
#ifndef __ICC  // internal error with icc22: lambda not supported
#ifdef HAVE_OPENMP4
#pragma omp parallel for if (len > cdoMinLoopSize) default(shared) schedule(static) reduction(max : vmax)
#endif
#endif
      for (size_t i = 0; i < len; ++i) f_max_mv(v[i], missval, vmax, is_equal);
    }

  if (is_equal(vmax, -std::numeric_limits<T>::max())) vmax = missval;

  return vmax;
}

// Explicit instantiation
template float varray_max_mv(size_t len, const Varray<float> &v, float missval);
template double varray_max_mv(size_t len, const Varray<double> &v, double missval);

template <typename T>
T
varray_range_mv(size_t len, const Varray<T> &v, T missval)
{
  auto f_minmax_mv = [](auto a, auto mv_a, auto &vmin, auto &vmax, auto is_EQ) {
    if (!is_EQ(a, mv_a))
      {
        vmin = min_value(vmin, a);
        vmax = max_value(vmax, a);
      }
  };

  assert(len > 0);
  assert(v.size() > 0);
  assert(len <= v.size());

  T vmin = std::numeric_limits<T>::max();
  T vmax = -std::numeric_limits<T>::max();

  if (std::isnan(missval))
    {
#ifndef __ICC  // internal error with icc22: lambda not supported
#ifdef HAVE_OPENMP4
#pragma omp parallel for if (len > cdoMinLoopSize) default(shared) schedule(static) reduction(min : vmin) reduction(max : vmax)
#endif
#endif
      for (size_t i = 0; i < len; ++i) f_minmax_mv(v[i], missval, vmin, vmax, dbl_is_equal);
    }
  else
    {
#ifndef __ICC  // internal error with icc22: lambda not supported
#ifdef HAVE_OPENMP4
#pragma omp parallel for if (len > cdoMinLoopSize) default(shared) schedule(static) reduction(min : vmin) reduction(max : vmax)
#endif
#endif
      for (size_t i = 0; i < len; ++i) f_minmax_mv(v[i], missval, vmin, vmax, is_equal);
    }

  return (is_equal(vmin, std::numeric_limits<T>::max()) && is_equal(vmax, -std::numeric_limits<T>::max())) ? missval : vmax - vmin;
}

// Explicit instantiation
template float varray_range_mv(size_t len, const Varray<float> &v, float missval);
template double varray_range_mv(size_t len, const Varray<double> &v, double missval);

double
array_sum(size_t len, const double *array)
{
  double sum = 0.0;
  for (size_t i = 0; i < len; ++i) sum += array[i];

  return sum;
}

template <typename T>
double
varray_sum(size_t len, const Varray<T> &v)
{
  // assert(len > 0); // failed in remapcon
  assert(v.size() > 0);
  assert(len <= v.size());

  double sum = 0.0;
#ifdef HAVE_OPENMP4
#pragma omp parallel for if (len > cdoMinLoopSize) default(shared) schedule(static) reduction(+ : sum)
#endif
  for (size_t i = 0; i < len; ++i) sum += v[i];

  return sum;
}

// Explicit instantiation
template double varray_sum(size_t len, const Varray<float> &v);
template double varray_sum(size_t len, const Varray<double> &v);

template <typename T>
double
varray_sum_mv(size_t len, const Varray<T> &v, T missval)
{
  auto f_sum_mv = [](auto a, auto mv_a, auto &sum, auto &nvals, auto is_EQ) {
    if (!is_EQ(a, mv_a))
      {
        sum += a;
        nvals++;
      }
  };

  assert(len > 0);
  assert(v.size() > 0);
  assert(len <= v.size());

  double sum = 0.0;
  size_t nvals = 0;

  if (std::isnan(missval))
    {
#ifndef __ICC  // internal error with icc22: lambda not supported
#ifdef HAVE_OPENMP4
#pragma omp parallel for if (len > cdoMinLoopSize) default(shared) schedule(static) reduction(+ : sum, nvals)
#endif
#endif
      for (size_t i = 0; i < len; ++i) f_sum_mv(v[i], missval, sum, nvals, dbl_is_equal);
    }
  else
    {
#ifndef __ICC  // internal error with icc22: lambda not supported
#ifdef HAVE_OPENMP4
#pragma omp parallel for if (len > cdoMinLoopSize) default(shared) schedule(static) reduction(+ : sum, nvals)
#endif
#endif
      for (size_t i = 0; i < len; ++i) f_sum_mv(v[i], missval, sum, nvals, is_equal);
    }

  if (!nvals) sum = missval;

  return sum;
}

// Explicit instantiation
template double varray_sum_mv(size_t len, const Varray<float> &v, float missval);
template double varray_sum_mv(size_t len, const Varray<double> &v, double missval);

template <typename T>
double
varray_mean(size_t len, const Varray<T> &v)
{
  assert(len > 0);
  assert(v.size() > 0);
  assert(len <= v.size());

  auto sum = varray_sum(len, v);

  return sum / len;
}

// Explicit instantiation
template double varray_mean(size_t len, const Varray<float> &v);
template double varray_mean(size_t len, const Varray<double> &v);

template <typename T>
double
varray_mean_mv(size_t len, const Varray<T> &v, T missval)
{
  assert(len > 0);
  assert(v.size() > 0);
  assert(len <= v.size());

  double sum = 0.0, sumw = 0.0;

  for (size_t i = 0; i < len; ++i)
    if (!dbl_is_equal(v[i], missval))
      {
        sum += v[i];
        sumw += 1;
      }

  double missval1 = missval, missval2 = missval;
  return DIVMN(sum, sumw);
}

// Explicit instantiation
template double varray_mean_mv(size_t len, const Varray<float> &v, float missval);
template double varray_mean_mv(size_t len, const Varray<double> &v, double missval);

template <typename T>
double
varray_weighted_mean(size_t len, const Varray<T> &v, const Varray<double> &w, T missval)
{
  auto f_weighted_mean = [](auto aw, auto a, auto &sum, auto &sumw) {
    sum += aw * a;
    sumw += aw;
  };

  assert(len > 0);
  assert(v.size() > 0);
  assert(len <= v.size());
  assert(len <= w.size());

  double sum = 0.0, sumw = 0.0;

#ifndef __ICC  // internal error with icc22: lambda not supported
#ifdef HAVE_OPENMP4
#pragma omp parallel for if (len > cdoMinLoopSize) default(shared) schedule(static) reduction(+ : sum, sumw)
#endif
#endif
  for (size_t i = 0; i < len; ++i) f_weighted_mean(w[i], v[i], sum, sumw);

  return is_equal(sumw, 0.0) ? missval : sum / sumw;
}

// Explicit instantiation
template double varray_weighted_mean(size_t len, const Varray<float> &v, const Varray<double> &w, float missval);
template double varray_weighted_mean(size_t len, const Varray<double> &v, const Varray<double> &w, double missval);

template <typename T>
double
varray_weighted_mean_mv(size_t len, const Varray<T> &v, const Varray<double> &w, T missval)
{
  auto f_weighted_mean_mv = [](auto aw, auto a, auto mv_a, auto &sum, auto &sumw, auto is_EQ) {
    if (!is_EQ(a, mv_a) && !is_EQ(aw, mv_a))
      {
        sum += aw * a;
        sumw += aw;
      }
  };

  assert(len > 0);
  assert(v.size() > 0);
  assert(len <= v.size());
  assert(len <= w.size());

  double missval1 = missval, missval2 = missval;
  double sum = 0.0, sumw = 0.0;

  if (std::isnan(missval))
    {
#ifndef __ICC  // internal error with icc22: lambda not supported
#ifdef HAVE_OPENMP4
#pragma omp parallel for if (len > cdoMinLoopSize) default(shared) schedule(static) reduction(+ : sum, sumw)
#endif
#endif
      for (size_t i = 0; i < len; ++i) f_weighted_mean_mv(w[i], v[i], missval1, sum, sumw, dbl_is_equal);
    }
  else
    {
#ifndef __ICC  // internal error with icc22: lambda not supported
#ifdef HAVE_OPENMP4
#pragma omp parallel for if (len > cdoMinLoopSize) default(shared) schedule(static) reduction(+ : sum, sumw)
#endif
#endif
      for (size_t i = 0; i < len; ++i) f_weighted_mean_mv(w[i], v[i], missval1, sum, sumw, is_equal);
    }

  return DIVMN(sum, sumw);
}

// Explicit instantiation
template double varray_weighted_mean_mv(size_t len, const Varray<float> &v, const Varray<double> &w, float missval);
template double varray_weighted_mean_mv(size_t len, const Varray<double> &v, const Varray<double> &w, double missval);

template <typename T>
double
varray_avg_mv(size_t len, const Varray<T> &v, T missval)
{
  assert(len > 0);
  assert(v.size() > 0);
  assert(len <= v.size());

  double missval1 = missval, missval2 = missval;
  double sum = 0.0, sumw = 0.0;

  for (size_t i = 0; i < len; ++i)
    {
      sum = ADDMN(sum, v[i]);
      sumw += 1;
    }

  return DIVMN(sum, sumw);
}

// Explicit instantiation
template double varray_avg_mv(size_t len, const Varray<float> &v, float missval);
template double varray_avg_mv(size_t len, const Varray<double> &v, double missval);

template <typename T>
double
varray_weighted_avg_mv(size_t len, const Varray<T> &v, const Varray<double> &w, T missval)
{
  assert(len > 0);
  assert(v.size() > 0);
  assert(len <= v.size());
  assert(len <= w.size());

  double missval1 = missval, missval2 = missval;
  double sum = 0.0, sumw = 0.0;

  for (size_t i = 0; i < len; ++i)
    if (!dbl_is_equal(w[i], missval))
      {
        sum = ADDMN(sum, MULMN(w[i], v[i]));
        sumw = ADDMN(sumw, w[i]);
      }

  return DIVMN(sum, sumw);
}

// Explicit instantiation
template double varray_weighted_avg_mv(size_t len, const Varray<float> &v, const Varray<double> &w, float missval);
template double varray_weighted_avg_mv(size_t len, const Varray<double> &v, const Varray<double> &w, double missval);

template <typename T>
static void
varray_prevarsum0(size_t len, const Varray<T> &v, double &rsum, double &rsumw)
{
  rsum = 0.0;
#ifdef HAVE_OPENMP4
#pragma omp parallel for if (len > cdoMinLoopSize) default(shared) schedule(static) reduction(+ : rsum)
#endif
  for (size_t i = 0; i < len; ++i) rsum += v[i];
  rsumw = len;
}

template <typename T>
static void
varray_prevarsum0_mv(size_t len, const Varray<T> &v, double missval, double &rsum, double &rsumw)
{
  auto f_prevarsum0_mv = [](auto a, auto mv_a, auto &sum, auto &sumw) {
    if (!dbl_is_equal(a, mv_a))
      {
        sum += a;
        sumw += 1.0;
      }
  };

  rsum = rsumw = 0.0;

#ifndef __ICC  // internal error with icc22: lambda not supported
#ifdef HAVE_OPENMP4
#pragma omp parallel for if (len > cdoMinLoopSize) default(shared) schedule(static) reduction(+ : rsum, rsumw)
#endif
#endif
  for (size_t i = 0; i < len; ++i) f_prevarsum0_mv(v[i], missval, rsum, rsumw);
}

template <typename T>
static void
varray_prevarsum(size_t len, const Varray<T> &v, double &rsum, double &rsumw, double &rsumq, double &rsumwq)
{
  auto f_prevarsum = [](auto a, auto &sum, auto &sumq) {
    sum += a;
    sumq += a * a;
  };

  rsum = rsumq = 0.0;

#ifndef __ICC  // internal error with icc22: lambda not supported
#ifdef HAVE_OPENMP4
#pragma omp parallel for if (len > cdoMinLoopSize) default(shared) schedule(static) reduction(+ : rsum, rsumq)
#endif
#endif
  for (size_t i = 0; i < len; ++i) f_prevarsum((double) v[i], rsum, rsumq);

  rsumw = len;
  rsumwq = len;
}

template <typename T>
static void
varray_prevarsum_mv(size_t len, const Varray<T> &v, T missval, double &rsum, double &rsumw, double &rsumq, double &rsumwq)
{
  auto f_prevarsum = [](auto a, auto mv_a, auto &sum, auto &sumq, auto &sumw, auto &sumwq) {
    if (!dbl_is_equal(a, mv_a))
      {
        double ad = (double) a;
        sum += ad;
        sumq += ad * ad;
        sumw += 1.0;
        sumwq += 1.0;
      }
  };

  rsum = rsumq = rsumw = rsumwq = 0.0;

#ifndef __ICC  // internal error with icc22: lambda not supported
#ifdef HAVE_OPENMP4
#pragma omp parallel for if (len > cdoMinLoopSize) default(shared) schedule(static) reduction(+ : rsum, rsumq, rsumw, rsumwq)
#endif
#endif
  for (size_t i = 0; i < len; ++i) f_prevarsum(v[i], missval, rsum, rsumq, rsumw, rsumwq);
}

template <typename T>
double
varray_var(size_t len, const Varray<T> &v, size_t nmiss, T missval)
{
  double rsum = 0.0, rsumw = 0.0, rsumq = 0.0, rsumwq = 0.0;
  if (nmiss > 0)
    varray_prevarsum_mv(len, v, missval, rsum, rsumw, rsumq, rsumwq);
  else
    varray_prevarsum(len, v, rsum, rsumw, rsumq, rsumwq);

  auto rvar = is_not_equal(rsumw, 0.0) ? (rsumq * rsumw - rsum * rsum) / (rsumw * rsumw) : missval;
  if (rvar < 0.0 && rvar > -1.e-5) rvar = 0.0;

  return rvar;
}

// Explicit instantiation
template double varray_var(size_t len, const Varray<float> &v, size_t nmiss, float missval);
template double varray_var(size_t len, const Varray<double> &v, size_t nmiss, double missval);

template <typename T>
double
varray_var_1(size_t len, const Varray<T> &v, size_t nmiss, T missval)
{
  double rsum = 0.0, rsumw = 0.0, rsumq = 0.0, rsumwq = 0.0;
  if (nmiss > 0)
    varray_prevarsum_mv(len, v, missval, rsum, rsumw, rsumq, rsumwq);
  else
    varray_prevarsum(len, v, rsum, rsumw, rsumq, rsumwq);

  auto rvar = (rsumw * rsumw > rsumwq) ? (rsumq * rsumw - rsum * rsum) / (rsumw * rsumw - rsumwq) : missval;
  if (rvar < 0.0 && rvar > -1.e-5) rvar = 0.0;

  return rvar;
}

// Explicit instantiation
template double varray_var_1(size_t len, const Varray<float> &v, size_t nmiss, float missval);
template double varray_var_1(size_t len, const Varray<double> &v, size_t nmiss, double missval);

template <typename T>
static void
varray_weighted_prevarsum(size_t len, const Varray<T> &v, const Varray<double> &w, double &rsum, double &rsumw, double &rsumq,
                          double &rsumwq)
{
  auto f_weighted_prevarsum = [](auto aw, auto a, auto &sum, auto &sumq, auto &sumw, auto &sumwq) {
    sum += aw * a;
    sumq += aw * a * a;
    sumw += aw;
    sumwq += aw * aw;
  };

  rsum = rsumq = rsumw = rsumwq = 0.0;

#ifndef __ICC  // internal error with icc22: lambda not supported
#ifdef HAVE_OPENMP4
#pragma omp parallel for if (len > cdoMinLoopSize) default(shared) schedule(static) reduction(+ : rsum, rsumq, rsumw, rsumwq)
#endif
#endif
  for (size_t i = 0; i < len; ++i) f_weighted_prevarsum(w[i], (double) v[i], rsum, rsumq, rsumw, rsumwq);
}

template <typename T>
static void
varray_weighted_prevarsum_mv(size_t len, const Varray<T> &v, const Varray<double> &w, double missval, double &rsum, double &rsumw,
                             double &rsumq, double &rsumwq)
{
  auto f_weighted_prevarsum_mv = [](auto aw, auto a, auto mv_a, auto &sum, auto &sumq, auto &sumw, auto &sumwq) {
    if (!dbl_is_equal(a, mv_a) && !dbl_is_equal(aw, mv_a))
      {
        sum += aw * a;
        sumq += aw * a * a;
        sumw += aw;
        sumwq += aw * aw;
      }
  };

  rsum = rsumq = rsumw = rsumwq = 0.0;

#ifndef __ICC  // internal error with icc22: lambda not supported
#ifdef HAVE_OPENMP4
#pragma omp parallel for if (len > cdoMinLoopSize) default(shared) schedule(static) reduction(+ : rsum, rsumq, rsumw, rsumwq)
#endif
#endif
  for (size_t i = 0; i < len; ++i) f_weighted_prevarsum_mv(w[i], (double) v[i], missval, rsum, rsumq, rsumw, rsumwq);
}

template <typename T>
double
varray_weighted_var(size_t len, const Varray<T> &v, const Varray<double> &w, size_t nmiss, T missval)
{
  double rsum = 0.0, rsumw = 0.0, rsumq = 0.0, rsumwq = 0.0;
  if (nmiss > 0)
    varray_weighted_prevarsum_mv(len, v, w, missval, rsum, rsumw, rsumq, rsumwq);
  else
    varray_weighted_prevarsum(len, v, w, rsum, rsumw, rsumq, rsumwq);

  auto rvar = is_not_equal(rsumw, 0) ? (rsumq * rsumw - rsum * rsum) / (rsumw * rsumw) : missval;
  if (rvar < 0.0 && rvar > -1.e-5) rvar = 0.0;

  return rvar;
}

// Explicit instantiation
template double varray_weighted_var(size_t len, const Varray<float> &v, const Varray<double> &w, size_t nmiss,
                                    float missval);
template double varray_weighted_var(size_t len, const Varray<double> &v, const Varray<double> &w, size_t nmiss,
                                    double missval);

template <typename T>
double
varray_weighted_var_1(size_t len, const Varray<T> &v, const Varray<double> &w, size_t nmiss, T missval)
{
  double rsum = 0.0, rsumw = 0.0, rsumq = 0.0, rsumwq = 0.0;
  if (nmiss > 0)
    varray_weighted_prevarsum_mv(len, v, w, missval, rsum, rsumw, rsumq, rsumwq);
  else
    varray_weighted_prevarsum(len, v, w, rsum, rsumw, rsumq, rsumwq);

  auto rvar = (rsumw * rsumw > rsumwq) ? (rsumq * rsumw - rsum * rsum) / (rsumw * rsumw - rsumwq) : missval;
  if (rvar < 0.0 && rvar > -1.e-5) rvar = 0.0;

  return rvar;
}

// Explicit instantiation
template double varray_weighted_var_1(size_t len, const Varray<float> &v, const Varray<double> &w, size_t nmiss,
                                      float missval);
template double varray_weighted_var_1(size_t len, const Varray<double> &v, const Varray<double> &w, size_t nmiss,
                                      double missval);

template <typename T>
static void
varray_prekurtsum(size_t len, const Varray<T> &v, double mean, double &rsum3w, double &rsum2diff, double &rsum4diff)
{
  auto f_prekurtsum = [](auto vdiff, auto &sum2diff, auto &sum4diff) {
    sum2diff += vdiff * vdiff;
    sum4diff += vdiff * vdiff * vdiff * vdiff;
  };

  rsum2diff = rsum4diff = 0.0;

#ifndef __ICC  // internal error with icc22: lambda not supported
#ifdef HAVE_OPENMP4
#pragma omp parallel for if (len > cdoMinLoopSize) default(shared) schedule(static) reduction(+ : rsum2diff, rsum4diff)
#endif
#endif
  for (size_t i = 0; i < len; ++i) f_prekurtsum(v[i] - mean, rsum2diff, rsum4diff);

  rsum3w = len;
}

template <typename T>
static void
varray_prekurtsum_mv(size_t len, const Varray<T> &v, T missval, double mean, double &rsum3w, double &rsum2diff,
                     double &rsum4diff)
{
  auto f_preskewsum_mv = [](auto a, auto mv_a, auto meanval, auto &sum2diff, auto &sum4diff, auto &sum3w) {
    if (!dbl_is_equal(a, mv_a))
      {
        double vdiff = a - meanval;
        sum2diff += vdiff * vdiff;
        sum4diff += vdiff * vdiff * vdiff * vdiff;
        sum3w += 1;
      }
  };

  rsum3w = rsum2diff = rsum4diff = 0.0;

#ifndef __ICC  // internal error with icc22: lambda not supported
#ifdef HAVE_OPENMP4
#pragma omp parallel for if (len > cdoMinLoopSize) default(shared) schedule(static) reduction(+ : rsum2diff, rsum4diff, rsum3w)
#endif
#endif
  for (size_t i = 0; i < len; ++i) f_preskewsum_mv(v[i], missval, mean, rsum2diff, rsum4diff, rsum3w);
}

template <typename T>
double
varray_kurt(size_t len, const Varray<T> &v, size_t nmiss, T missval)
{
  double rsum3w;  // 3rd moment variables
  double rsum2diff, rsum4diff;
  double rsum, rsumw;

  if (nmiss > 0)
    {
      varray_prevarsum0_mv(len, v, missval, rsum, rsumw);
      varray_prekurtsum_mv(len, v, missval, (rsum / rsumw), rsum3w, rsum2diff, rsum4diff);
    }
  else
    {
      varray_prevarsum0(len, v, rsum, rsumw);
      varray_prekurtsum(len, v, (rsum / rsumw), rsum3w, rsum2diff, rsum4diff);
    }

  if (is_equal(rsum3w, 0.0) || is_equal(rsum2diff, 0.0)) return missval;

  auto rkurt = ((rsum4diff / rsum3w) / std::pow(rsum2diff / rsum3w, 2)) - 3.0;
  if (rkurt < 0.0 && rkurt > -1.e-5) rkurt = 0.0;

  return rkurt;
}

// Explicit instantiation
template double varray_kurt(size_t len, const Varray<float> &v, size_t nmiss, float missval);
template double varray_kurt(size_t len, const Varray<double> &v, size_t nmiss, double missval);

template <typename T>
static void
varray_preskewsum(size_t len, const Varray<T> &v, double mean, double &rsum3w, double &rsum3diff, double &rsum2diff)
{
  auto f_preskewsum = [](auto vdiff, auto &sum3diff, auto &sum2diff) {
    sum3diff += vdiff * vdiff * vdiff;
    sum2diff += vdiff * vdiff;
  };

  rsum2diff = 0.0;
  rsum3diff = 0.0;

#ifndef __ICC  // internal error with icc22: lambda not supported
#ifdef HAVE_OPENMP4
#pragma omp parallel for if (len > cdoMinLoopSize) default(shared) schedule(static) reduction(+ : rsum2diff, rsum3diff)
#endif
#endif
  for (size_t i = 0; i < len; ++i) f_preskewsum(v[i] - mean, rsum3diff, rsum2diff);

  rsum3w = len;
}

template <typename T>
static void
varray_preskewsum_mv(size_t len, const Varray<T> &v, T missval, double mean, double &rsum3w, double &rsum3diff,
                     double &rsum2diff)
{
  auto f_preskewsum_mv = [](auto a, auto mv_a, auto meanval, auto &sum3diff, auto &sum2diff, auto &sum3w) {
    if (!dbl_is_equal(a, mv_a))
      {
        double vdiff = a - meanval;
        sum3diff += vdiff * vdiff * vdiff;
        sum2diff += vdiff * vdiff;
        sum3w += 1;
      }
  };

  rsum3w = rsum3diff = rsum2diff = 0.0;

#ifndef __ICC  // internal error with icc22: lambda not supported
#ifdef HAVE_OPENMP4
#pragma omp parallel for if (len > cdoMinLoopSize) default(shared) schedule(static) reduction(+ : rsum2diff, rsum3diff, rsum3w)
#endif
#endif
  for (size_t i = 0; i < len; ++i) f_preskewsum_mv(v[i], missval, mean, rsum3diff, rsum2diff, rsum3w);
}

template <typename T>
double
varray_skew(size_t len, const Varray<T> &v, size_t nmiss, T missval)
{
  double rsum3w;  // 3rd moment variables
  double rsum3diff, rsum2diff;
  double rsum, rsumw;

  if (nmiss > 0)
    {
      varray_prevarsum0_mv(len, v, missval, rsum, rsumw);
      varray_preskewsum_mv(len, v, missval, (rsum / rsumw), rsum3w, rsum3diff, rsum2diff);
    }
  else
    {
      varray_prevarsum0(len, v, rsum, rsumw);
      varray_preskewsum(len, v, (rsum / rsumw), rsum3w, rsum3diff, rsum2diff);
    }

  if (is_equal(rsum3w, 0.0) || is_equal(rsum3w, 1.0) || is_equal(rsum2diff, 0.0)) return missval;

  auto rskew = (rsum3diff / rsum3w) / std::pow((rsum2diff) / (rsum3w - 1.0), 1.5);
  if (rskew < 0.0 && rskew > -1.e-5) rskew = 0.0;

  return rskew;
}

// Explicit instantiation
template double varray_skew(size_t len, const Varray<float> &v, size_t nmiss, float missval);
template double varray_skew(size_t len, const Varray<double> &v, size_t nmiss, double missval);

#include <algorithm>

template <typename T>
static double
get_nth_element(T *array, size_t length, size_t n)
{
  std::nth_element(array, array + n, array + length);
  return array[n];
}

template <typename T>
double
varray_median(size_t len, const Varray<T> &v, size_t nmiss, T missval)
{
  assert(len > 0);
  assert(v.size() > 0);
  assert(len <= v.size());

  double median = missval;

  if (nmiss == 0)
    {
      Varray<T> v2 = v;
      if (len % 2 == 0)
        {
          auto k = len / 2;
          auto vk1 = get_nth_element(v2.data(), len, k - 1);
          auto vk2 = get_nth_element(v2.data(), len, k);
          median = (vk1 + vk2) * 0.5;
        }
      else
        {
          auto k = (len + 1) / 2;
          median = get_nth_element(v2.data(), len, k - 1);
        }
    }

  return median;
}

// Explicit instantiation
template double varray_median(size_t len, const Varray<float> &v, size_t nmiss, float missval);
template double varray_median(size_t len, const Varray<double> &v, size_t nmiss, double missval);

template <typename T>
double
varray_count(size_t len, const Varray<T> &v, size_t nmiss, T missval)
{
  assert(len > 0);
  assert(v.size() > 0);
  assert(len <= v.size());

  size_t count = len;

  if (nmiss > 0)
    {
      count = 0;
      for (size_t i = 0; i < len; ++i)
        {
          if (!dbl_is_equal(v[i], missval)) count++;
        }
    }

  return count;
}

// Explicit instantiation
template double varray_count(size_t len, const Varray<float> &v, size_t nmiss, float missval);
template double varray_count(size_t len, const Varray<double> &v, size_t nmiss, double missval);
