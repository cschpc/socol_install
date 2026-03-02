/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/
#ifndef VARRAY_H
#define VARRAY_H

#include <iostream>  // cerr
#include <vector>
#include <cstddef>
#include <cfloat>
#include <cassert>
#include "compare.h"

// compare
// clang-format off
const auto is_not_equal = [](auto a, auto b) noexcept { return  (a < b || b < a); };
const auto is_equal     = [](auto a, auto b) noexcept { return !(a < b || b < a); };
const auto dbl_is_equal = [](auto a, auto b) noexcept { return (std::isnan(a) || std::isnan(b) ? (std::isnan(a) && std::isnan(b)) : is_equal(a, b)); };
// clang-format on

// unary operators
// clang-format off
const auto unary_op_abs   = [](double x) noexcept { return std::fabs(x); };
const auto unary_op_int   = [](double x) noexcept { return (int)(x); };
const auto unary_op_nint  = [](double x) noexcept { return std::round(x); };
const auto unary_op_sqr   = [](double x) noexcept { return x * x; };
const auto unary_op_nop   = [](double x) noexcept { return x; };
const auto unary_op_reci  = [](double x) noexcept { return 1.0 / x; };
const auto unary_op_not   = [](double x) noexcept { return is_equal(x, 0.0); };
const auto unary_op_exp   = [](double x) noexcept { return std::exp(x); };
const auto unary_op_log   = [](double x) noexcept { return std::log(x); };
const auto unary_op_log10 = [](double x) noexcept { return std::log10(x); };
const auto unary_op_sin   = [](double x) noexcept { return std::sin(x); };
const auto unary_op_cos   = [](double x) noexcept { return std::cos(x); };
const auto unary_op_tan   = [](double x) noexcept { return std::tan(x); };
const auto unary_op_asin  = [](double x) noexcept { return std::asin(x); };
const auto unary_op_acos  = [](double x) noexcept { return std::acos(x); };
const auto unary_op_atan  = [](double x) noexcept { return std::atan(x); };
// clang-format on

// binary operators
// clang-format off
const auto binary_op_LT  = [](double x, double y) noexcept { return x < y; };
const auto binary_op_GT  = [](double x, double y) noexcept { return x > y; };
const auto binary_op_LE  = [](double x, double y) noexcept { return x <= y; };
const auto binary_op_GE  = [](double x, double y) noexcept { return x >= y; };
const auto binary_op_NE  = [](double x, double y) noexcept { return is_not_equal(x, y); };
const auto binary_op_EQ  = [](double x, double y) noexcept { return is_equal(x, y); };
const auto binary_op_LEG = [](double x, double y) noexcept { return (x < y) ? -1.0 : (x > y); };
const auto binary_op_AND = [](double x, double y) noexcept { return is_not_equal(x, 0.0) && is_not_equal(y, 0.0); };
const auto binary_op_OR  = [](double x, double y) noexcept { return is_not_equal(x, 0.0) || is_not_equal(y, 0.0); };
const auto binary_op_POW = [](double x, double y) noexcept { return std::pow(x, y); };
const auto binary_op_ADD = [](double x, double y) noexcept { return x + y; };
const auto binary_op_SUB = [](double x, double y) noexcept { return x - y; };
const auto binary_op_MUL = [](double x, double y) noexcept { return x * y; };
const auto binary_op_DIV = [](double x, double y) noexcept { return x / y; };
// clang-format on

//#define CHECK_UNUSED_VECTOR 1

#ifdef CHECK_UNUSED_VECTOR
// clang-format off
template <typename T>
class
#ifdef WARN_UNUSED
[[gnu::warn_unused]]
#endif
CheckVector
{
 public:
  T *ptr;
  size_t m_count;

  CheckVector() : ptr(dummy), m_count(0) { }
  explicit CheckVector(size_t count) : ptr(dummy), m_count(count) { }
  CheckVector(size_t count, const T& value) : ptr(dummy), m_count(count) { ptr[0] = value; }

  T * begin() noexcept { return &ptr[0]; }
  T * end() noexcept { return &ptr[0] + 1; }
  const T * begin() const noexcept { return &ptr[0]; }
  const T * end() const noexcept { return &ptr[0] + 1; }
  
  bool empty() const { return true; }
  size_t size() const { return m_count; }
  void clear() { }
  void shrink_to_fit() { }

  void resize(size_t count) { m_count = count; }
  void resize(size_t count, const T& value) { m_count = count; ptr[0] = value; }

  void reserve(size_t new_cap) { (void)new_cap; }
  void push_back(const T& value) { ptr[0] = value; }
  void assign(size_t count, const T& value) { (void)count; ptr[0] = value; }

  T * data() noexcept { return ptr; }
  const T * data() const noexcept { return ptr; }

  T& operator[](size_t pos) { (void)pos; return ptr[0]; }
  const T & operator[](size_t pos) const { (void)pos; return ptr[0]; }
  CheckVector& operator=(const CheckVector& other) { *this = other; ptr = dummy; m_count = 0; return *this; }
  CheckVector& operator=(CheckVector&& other) { *this = other; ptr = dummy; m_count = 0; return *this; }
  CheckVector& operator=(const std::vector<T>& other) { (void)other; ptr = dummy; m_count = 0; return *this; }
  CheckVector& operator=(std::vector<T>&& other) { (void)other; ptr = dummy; m_count = 0; return *this; }
  // Copy constructor
  CheckVector(const CheckVector &v2) { m_count = 0; ptr = dummy; ptr[0] = v2.ptr[0]; }
  explicit CheckVector(const std::vector<T> &v2) { (void)v2; m_count = 0; ptr = dummy; }

  bool operator==(const CheckVector& other) { (void)other; return true; }

 private:
  T dummy[1] {};
};
// clang-format on

template <typename T>
using Varray = CheckVector<T>;

#else

template <typename T>
using Varray = std::vector<T>;

#endif

template <typename T>
using Varray2D = Varray<Varray<T>>;

template <typename T>
using Varray3D = Varray2D<Varray<T>>;

template <typename T>
using Varray4D = Varray3D<Varray<T>>;

// clang-format off
struct MinMax
{
  double min = DBL_MAX;
  double max = -DBL_MAX;
  size_t n = 0;
  MinMax() {};
  MinMax(double rmin, double rmax, size_t rn) : min(rmin), max(rmax), n(rn) {};
  MinMax(double rmin, double rmax) : min(rmin), max(rmax), n(0) {};
};

struct MinMaxSum : MinMax
{
  double sum = 0.0;
  MinMaxSum() {};
  MinMaxSum(double rmin, double rmax, double rsum, size_t rn) : sum(rsum) { min = rmin; max = rmax; n = rn; };
  MinMaxSum(double rmin, double rmax, double rsum) : sum(rsum) { min = rmin; max = rmax; n = 0; };
};

struct MinMaxMean : MinMax
{
  double mean = 0.0;
  MinMaxMean() {};
  MinMaxMean(double rmin, double rmax, double rmean, size_t rn) : mean(rmean) { min = rmin; max = rmax; n = rn; };
  MinMaxMean(double rmin, double rmax, double rmean) : mean(rmean) { min = rmin; max = rmax; n = 0; };
};
// clang-format on

template <typename T>
inline void
varray_free(Varray<T> &v)
{
  v.clear();
  v.shrink_to_fit();
}

#define varrayResize(p, s) varray_resize((p), (s), __FILE__, __LINE__)
#define varrayResizeInit(p, s, v) varray_resize((p), (s), (v), __FILE__, __LINE__)

template <typename T>
inline void
varray_resize(Varray<T> &v, size_t count, const char *file, int line)
{
  if (v.size() != count)
    {
      try
        {
          v.resize(count);
        }
      catch (const std::exception &e)
        {
          std::cerr << "Exception caught when trying to allocate " << count << " vector elements: " << e.what() << " in " << file
                    << ":" << line << '\n';
          throw;
        }
    }
}

template <typename T>
inline void
varray_resize(Varray<T> &v, size_t count, T value, const char *file, int line)
{
  try
    {
      v.resize(count, value);
    }
  catch (const std::exception &e)
    {
      std::cerr << "Exception caught when trying to allocate " << count << " vector elements: " << e.what() << " in " << file << ":"
                << line << '\n';
      throw;
    }
}

/*
template <class T, size_t ROW, size_t COL>
using Matrix = std::array<std::array<T, COL>, ROW>;
*/

// clang-format off
#define MADDMN(x, y) ((DBL_IS_EQUAL((x), missval1) || DBL_IS_EQUAL((y), missval2)) ? missval1 : (x) + (y))
#define MSUBMN(x, y) ((DBL_IS_EQUAL((x), missval1) || DBL_IS_EQUAL((y), missval2)) ? missval1 : (x) - (y))
#define MMULMN(x, y) ((DBL_IS_EQUAL((x), 0) || DBL_IS_EQUAL((y), 0)) ? 0 : (DBL_IS_EQUAL((x), missval1) || DBL_IS_EQUAL((y), missval2)) ? missval1 : (x) * (y))
#define MDIVMN(x, y) ((DBL_IS_EQUAL((x), missval1) || DBL_IS_EQUAL((y), missval2) || DBL_IS_EQUAL((y), 0)) ? missval1 : (x) / (y))
#define MPOWMN(x, y) ((DBL_IS_EQUAL((x), missval1) || DBL_IS_EQUAL((y), missval2)) ? missval1 : std::pow((x), (y)))
#define MSQRTMN(x) ((DBL_IS_EQUAL((x), missval1) || (x) < 0) ? missval1 : std::sqrt(x))

#define ADD(x, y) ((x) + (y))
#define SUB(x, y) ((x) - (y))
#define MUL(x, y) ((x) * (y))
#define DIV(x, y) (IS_EQUAL((y), 0.) ? missval1 : (x) / (y))
#define POW(x, y) std::pow((x), (y))
#define SQRT(x) std::sqrt(x)

#define ADDM(x, y) ((IS_EQUAL((x), missval1) || IS_EQUAL((y), missval2)) ? missval1 : (x) + (y))
#define SUBM(x, y) ((IS_EQUAL((x), missval1) || IS_EQUAL((y), missval2)) ? missval1 : (x) - (y))
#define MULM(x, y) ((IS_EQUAL((x), 0) || IS_EQUAL((y), 0)) ? 0 : (IS_EQUAL((x), missval1) || IS_EQUAL((y), missval2)) ? missval1 : (x) * (y))
#define DIVM(x, y) ((IS_EQUAL((x), missval1) || IS_EQUAL((y), missval2) || IS_EQUAL((y), 0)) ? missval1 : (x) / (y))
#define POWM(x, y) ((IS_EQUAL((x), missval1) || IS_EQUAL((y), missval2)) ? missval1 : std::pow((x), (y)))
#define SQRTM(x) ((IS_EQUAL((x), missval1) || (x) < 0) ? missval1 : std::sqrt(x))

#define ADDMN(x, y) FADDMN(x, y, missval1, missval2)
#define SUBMN(x, y) FSUBMN(x, y, missval1, missval2)
#define MULMN(x, y) FMULMN(x, y, missval1, missval2)
#define DIVMN(x, y) FDIVMN(x, y, missval1, missval2)
#define POWMN(x, y) FPOWMN(x, y, missval1, missval2)
#define SQRTMN(x) FSQRTMN(x, missval1)

static inline double FADDMN(double x, double y, double missval1, double missval2)
{
  return MADDMN(x, y);
}
static inline double FSUBMN(double x, double y, double missval1, double missval2)
{
  return MSUBMN(x, y);
}
static inline double FMULMN(double x, double y, double missval1, double missval2)
{
  return MMULMN(x, y);
}
static inline double FDIVMN(double x, double y, double missval1, double missval2)
{
  return MDIVMN(x, y);
}
static inline double FPOWMN(double x, double y, double missval1, double missval2)
{
  return MPOWMN(x, y);
}
static inline double FSQRTMN(double x, double missval1)
{
  return MSQRTMN(x);
}
// clang-format on

const char *fpe_errstr(int fpeRaised);

template <typename T>
MinMaxSum varray_min_max_sum(size_t len, const Varray<T> &v, MinMaxSum mms);

template <typename T>
MinMaxSum varray_min_max_sum_mv(size_t len, const Varray<T> &v, T missval, MinMaxSum mms);

template <typename T>
MinMaxMean varray_min_max_mean(size_t len, const Varray<T> &v);

template <typename T>
MinMaxMean varray_min_max_mean_mv(size_t len, const Varray<T> &v, T missval);

template <typename T>
MinMax array_min_max_mask(size_t len, const T *array, const Varray<int> &mask);

void array_add_array(size_t len, double *array1, const double *array2);
void array_add_array_mv(size_t len, double *array1, const double *array2, double missval);

template <typename T1, typename T2>
void
varray_fill(size_t len, T1 *v, T2 value)
{
  T1 c = value;
  for (size_t i = 0; i < len; ++i) v[i] = c;
}

template <typename T1, typename T2>
void
varray_fill(Varray<T1> &v, T2 value)
{
  varray_fill(v.size(), v.data(), value);
}

template <typename T1, typename T2>
void
varray_fill(size_t len, Varray<T1> &v, T2 value)
{
  varray_fill(len, v.data(), value);
}

template <typename T1, typename T2>
void
array_copy(size_t len, const T1 *array1, T2 *array2)
{
  for (size_t i = 0; i < len; ++i) array2[i] = array1[i];
}

template <typename T1, typename T2>
void
varray_copy(size_t len, const Varray<T1> &v1, Varray<T2> &v2)
{
  for (size_t i = 0; i < len; ++i) v2[i] = v1[i];
}

template <typename T1, typename T2>
void
varray_divc(size_t len, Varray<T1> &v, T2 value)
{
  T1 c = value;
  for (size_t i = 0; i < len; ++i) v[i] /= c;
}

template <typename T, class UnaryOperation>
void
varray_transform(const Varray<T> &vIn, Varray<T> &vOut, UnaryOperation unary_op)
{
  assert(vIn.size() > 0);
  assert(vOut.size() > 0);
  assert(vOut.size() <= vIn.size());

  const auto len = vIn.size();
  for (size_t i = 0; i < len; ++i) vOut[i] = unary_op(vIn[i]);
}

template <typename T>
size_t array_num_mv(size_t len, const T *array, T missval);

template <typename T>
size_t varray_num_mv(size_t len, const Varray<T> &v, T missval);

template <typename T>
MinMax varray_min_max(size_t len, const Varray<T> &v);

template <typename T>
MinMax varray_min_max(size_t len, const T *array);

template <typename T>
MinMax varray_min_max(const Varray<T> &v);

template <typename T>
MinMax varray_min_max_mv(size_t len, const Varray<T> &v, T missval);

template <typename T>
MinMax varray_min_max_mv(size_t len, const T *array, T missval);

template <typename T>
T varray_min(size_t len, const Varray<T> &v);

template <typename T>
T varray_max(size_t len, const Varray<T> &v);

template <typename T>
T varray_range(size_t len, const Varray<T> &v);

template <typename T>
T varray_min_mv(size_t len, const Varray<T> &v, T missval);

template <typename T>
T varray_max_mv(size_t len, const Varray<T> &v, T missval);

template <typename T>
T varray_range_mv(size_t len, const Varray<T> &v, T missval);

double array_sum(size_t len, const double *array);

template <typename T>
double varray_sum(size_t len, const Varray<T> &v);

template <typename T>
double varray_sum_mv(size_t len, const Varray<T> &v, T missval);

template <typename T>
double varray_mean(size_t len, const Varray<T> &v);

template <typename T>
double varray_mean_mv(size_t len, const Varray<T> &v, T missval);

template <typename T>
double varray_weighted_mean(size_t len, const Varray<T> &v, const Varray<double> &w, T missval);

template <typename T>
double varray_weighted_mean_mv(size_t len, const Varray<T> &v, const Varray<double> &w, T missval);

template <typename T>
double varray_avg_mv(size_t len, const Varray<T> &v, T missval);

template <typename T>
double varray_weighted_avg_mv(size_t len, const Varray<T> &v, const Varray<double> &w, T missval);

template <typename T>
double varray_var(size_t len, const Varray<T> &v, size_t nmiss, T missval);

template <typename T>
double varray_var_1(size_t len, const Varray<T> &v, size_t nmiss, T missval);

template <typename T>
double varray_weighted_var(size_t len, const Varray<T> &v, const Varray<double> &w, size_t nmiss, T missval);

template <typename T>
double varray_weighted_var_1(size_t len, const Varray<T> &v, const Varray<double> &w, size_t nmiss, T missval);

template <typename T>
double varray_skew(size_t len, const Varray<T> &v, size_t nmiss, T missval);

template <typename T>
double varray_kurt(size_t len, const Varray<T> &v, size_t nmiss, T missval);

template <typename T>
double varray_median(size_t len, const Varray<T> &v, size_t nmiss, T missval);

template <typename T>
double varray_count(size_t len, const Varray<T> &v, size_t nmiss, T missval);

#endif  //  VARRAY_H
