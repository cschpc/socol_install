#ifndef COMPARE_H
#define COMPARE_H

#include <cmath>
#include <cstring>
#include <string>

//#define DBL_IS_EQUAL(x,y) (!(x < y || y < x))
//#define DBL_IS_EQUAL(x, y) (std::isnan(x) || std::isnan(y) ? (std::isnan(x) && std::isnan(y) ? 1 : 0) : !(x < y || y < x))
static inline bool
DBL_IS_EQUAL(double x, double y)
{
  return (std::isnan(x) || std::isnan(y) ? (std::isnan(x) && std::isnan(y)) : !(x < y || y < x));
}

//#define IS_NOT_EQUAL(x, y) (x < y || y < x)
//#define IS_EQUAL(x, y) (!IS_NOT_EQUAL(x, y))
constexpr bool
IS_NOT_EQUAL(double x, double y) noexcept
{
  return (x < y || y < x);
}

constexpr bool
IS_EQUAL(double x, double y) noexcept
{
  return !(x < y || y < x);
}

static inline bool
cdo_cmpstr(const char *x, const char *y)
{
  return (strcmp(x, y) == 0);
}

static inline bool
cdo_cmpstrLenRhs(const char *lhs, const char *rhs)
{
  return (strncmp(lhs, rhs, strlen(rhs)) == 0);
}

static inline bool
cdo_cmpstrLenRhs(const char *lhs, const char *rhs, size_t &len)
{
  return (strncmp(lhs, rhs, len = strlen(rhs)) == 0);
}

static inline bool
cdo_cmpstr(const std::string &lhs, const std::string &rhs)
{
  return (lhs.compare(rhs) == 0);
}

static inline bool
cdo_cmpstrLenRhs(const std::string &lhs, const std::string &rhs)
{
  return (lhs.compare(0, rhs.length(), rhs) == 0);
}

static inline bool
cdo_cmpstrLenRhs(const std::string &lhs, const std::string &rhs, size_t &len)
{
  return (lhs.compare(0, len = rhs.length(), rhs) == 0);
}

#endif /* COMPARE_H */
