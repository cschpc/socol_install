/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include "param_conversion.h"

#include <inttypes.h>
#include <climits>
#include <cassert>
#include <cctype>
#include <vector>
#include <cdi.h>

#include "cdo_output.h"
#include "cdo_options.h"
#include "compare.h"
#include "util_string.h"
#include "constants.h"
#include "const.h"

const char *
parameter_to_word(const char *cstring)
{
  auto len = strlen(cstring);

  for (size_t i = 0; i < len; ++i)
    {
      const int c = cstring[i];
      if (iscntrl(c) || isblank(c)) cdo_abort("Word parameter >%s< contains invalid character at position %zu!", cstring, i + 1);
    }

  if (len == 0) cdo_abort("Word parameter >%s< is empty!", cstring);

  return cstring;
}

static inline void
parameter_error(const char *name, const char *cstring, const char *endptr)
{
  cdo_abort("%s parameter >%s< contains invalid character at position %d!", name, cstring, (int) (endptr - cstring + 1));
}

long
parameter_to_bytes(const std::string &string)
{
  char *endptr = nullptr;
  auto numBytes = strtoimax(string.c_str(), &endptr, 10);
  if (*endptr)
    {
      switch (tolower((int) *endptr))
        {
          // clang-format off
        case 'k': numBytes *= 1024;       endptr++; break;
        case 'm': numBytes *= 1048576;    endptr++; break;
        case 'g': numBytes *= 1073741824; endptr++; break;
          // clang-format on
        }
    }
  if (*endptr) parameter_error("Bytes", string.c_str(), endptr);
  int len = string.size();
  if ((string[0] != '-' && len > 19) || len > 20)
    cdo_abort("Integer parameter %s too long (min=%ld max=%ld)!", string, LONG_MIN, LONG_MAX);
  return numBytes;
}

double
parameter_to_double(const char *const cstring)
{
  char *endptr = nullptr;
  auto fval = strtod(cstring, &endptr);
  if (*endptr && *endptr == 'f') endptr++;
  if (*endptr) parameter_error("Float", cstring, endptr);
  return fval;
}

static intmax_t
parameter_to_imax(const char *const cstring)
{
  char *endptr = nullptr;
  auto ival = strtoimax(cstring, &endptr, 10);
  if (*endptr) parameter_error("Integer", cstring, endptr);
  int len = strlen(cstring);
  if ((cstring[0] != '-' && len > 19) || len > 20)
    cdo_abort("Integer parameter %s too long (min=%ld max=%ld)!", cstring, LONG_MIN, LONG_MAX);
  return ival;
}

int
parameter_to_int(const char *const cstring)
{
  auto ival = parameter_to_imax(cstring);
  if (ival < INT_MIN || ival > INT_MAX)
    cdo_abort("Integer parameter >%s< out of range (min=%d max=%d)!", cstring, INT_MIN, INT_MAX);
  return static_cast<int>(ival);
}

long
parameter_to_long(const char *const cstring)
{
  auto ival = parameter_to_imax(cstring);
  if (ival < LONG_MIN || ival > LONG_MAX)
    cdo_abort("Integer parameter >%s< out of range (min=%ld max=%ld)!", cstring, LONG_MIN, LONG_MAX);
  return static_cast<long>(ival);
}

size_t
parameter_to_size_t(const char *const cstring)
{
  auto ival = parameter_to_imax(cstring);
  if (ival < 0 || ival > LONG_MAX)
    cdo_abort("Unsigned integer parameter >%s< out of range (min=%ld max=%ld)!", cstring, 0L, LONG_MAX);
  return static_cast<size_t>(ival);
}

int
parameter_to_intlist(const char *const cstring)
{
  char *endptr = nullptr;
  auto ival = (int) strtol(cstring, &endptr, 10);
  if (*endptr && *endptr != '/' && (endptr - cstring) == 0) parameter_error("Integer", cstring, endptr);
  return ival;
}

const std::string &
parameter_to_word(const std::string &string)
{
  auto len = string.size();

  for (size_t i = 0; i < len; ++i)
    {
      const int c = string[i];
      if (iscntrl(c) || isblank(c)) cdo_abort("Word parameter >%s< contains invalid character at position %zu!", string, i + 1);
    }

  if (len == 0) cdo_abort("Word parameter >%s< is empty!", string);

  return string;
}

bool
parameter_to_bool(const std::string &string)
{
  auto lstring = string_to_lower(string);

  if (lstring == "1" || lstring == "t" || lstring == "true") return true;
  if (lstring == "0" || lstring == "f" || lstring == "false") return false;

  cdo_abort("Boolean parameter >%s< contains invalid characters!", string);

  return false;
}

double
parameter_to_double(const std::string &string)
{
  return parameter_to_double(string.c_str());
}

int
parameter_to_int(const std::string &string)
{
  return parameter_to_int(string.c_str());
}

long
parameter_to_long(const std::string &string)
{
  return parameter_to_long(string.c_str());
}

size_t
parameter_to_size_t(const std::string &string)
{
  return parameter_to_size_t(string.c_str());
}

int
parameter_to_intlist(const std::string &string)
{
  return parameter_to_intlist(string.c_str());
}

double
radius_str_to_deg(const std::string &string)
{
  char *endptr = nullptr;
  auto radius = strtod(string.c_str(), &endptr);

  if (*endptr != 0)
    {
      if (strncmp(endptr, "km", 2) == 0)
        radius = 360 * ((radius * 1000) / (2 * PlanetRadius * M_PI));
      else if (strncmp(endptr, "m", 1) == 0)
        radius = 360 * ((radius) / (2 * PlanetRadius * M_PI));
      else if (strncmp(endptr, "deg", 3) == 0)
        ;
      else if (strncmp(endptr, "rad", 3) == 0)
        radius *= RAD2DEG;
      else
        cdo_abort("Float parameter >%s< contains invalid character at position %d!", string, (int) (endptr - string.c_str() + 1));
    }

  if (radius > 180.0) radius = 180.0;

  return radius;
}

int
string_to_param(const char *const paramstr)
{
  int pnum = -1, pcat = 255, pdis = 255;
  std::sscanf(paramstr, "%d.%d.%d", &pnum, &pcat, &pdis);

  return cdiEncodeParam(pnum, pcat, pdis);
}

int
string_to_param(const std::string &paramstr)
{
  return string_to_param(paramstr.c_str());
}

void
param_to_string(int param, char *paramstr, int maxlen)
{
  int dis, cat, num;
  cdiDecodeParam(param, &num, &cat, &dis);

  size_t umaxlen = (maxlen >= 0) ? (unsigned) maxlen : 0U;
  int len;
  if (dis == 255 && (cat == 255 || cat == 0))
    len = std::snprintf(paramstr, umaxlen, "%03d", num);
  else if (dis == 255)
    len = std::snprintf(paramstr, umaxlen, "%03d.%03d", num, cat);
  else
    len = std::snprintf(paramstr, umaxlen, "%03d.%03d.%03d", num, cat, dis);

  if (len >= maxlen || len < 0) cdo_abort("Internal problem (%s): size of input string is too small!", __func__);
}

/* time/date/season converisons */
/* =================================================================================== */
void
season_to_months(const char *season, int *imonths)
{
  const char *const smons = "JFMAMJJASONDJFMAMJJASOND";
  const int imons[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };
  assert(strlen(smons) == (sizeof(imons) / sizeof(int)));

  auto len = strlen(season);
  if (len == 3 && strcmp(season, "ANN") == 0)
    {
      for (size_t k = 0; k < 12; ++k) imonths[k + 1] = 1;
    }
  else
    {
      if (len > 12) cdo_abort("Too many months %d (limit=12)!", (int) len);
      char *const season_u = strdup(season);
      cstr_to_upper(season_u);
      const char *const sstr = strstr(smons, season_u);
      free(season_u);
      if (sstr != nullptr)
        {
          const size_t ks = (size_t) (sstr - smons);
          const size_t ke = ks + len;
          for (size_t k = ks; k < ke; ++k) imonths[imons[k]]++;
        }
      else
        {
          cdo_abort("Season %s not available!", season);
        }
    }
}

double
date_str_to_double(const char *datestr, const int opt)
{
  int year = 1, month = 1, day = 1, hour = 0, minute = 0, second = 0;
  double fval = 0.0;

  auto len = strlen(datestr);

  for (size_t i = 0; i < len; ++i)
    {
      const int c = datestr[i];
      if (!(isdigit(c) || c == '-' || c == ':' || c == '.' || c == 'T'))
        cdo_abort("Date string >%s< contains invalid character at position %zu!", datestr, i + 1);
    }

  if (opt)
    {
      hour = 23;
      minute = 59;
      second = 59;
    }

  if (strchr(datestr, '-') == nullptr) { fval = parameter_to_double(datestr); }
  else if (strchr(datestr, 'T'))
    {
      auto status = std::sscanf(datestr, "%d-%d-%dT%d:%d:%d", &year, &month, &day, &hour, &minute, &second);
      if (status != 6) cdo_abort("Invalid date string >%s<!", datestr);
      fval = cdiEncodeTime(hour, minute, second);
      if (std::fabs(fval) > 0) fval /= 1000000;
      fval += cdiEncodeDate(year, month, day);
    }
  else
    {
      auto status = std::sscanf(datestr, "%d-%d-%d", &year, &month, &day);
      if (status != 3) cdo_abort("Invalid date string >%s<!", datestr);
      fval = cdiEncodeTime(hour, minute, second);
      if (std::fabs(fval) > 0) fval /= 1000000;
      fval += cdiEncodeDate(year, month, day);
    }

  return fval;
}

// argv conversions ==============================================
std::vector<int>
cdo_argv_to_int(const std::vector<std::string> &argv)
{
  std::vector<int> v;

  for (const auto &argument : argv)
    {
      int first, last, inc;
      split_intstring(argument, first, last, inc);

      if (inc >= 0)
        {
          for (auto ival = first; ival <= last; ival += inc) v.push_back(ival);
        }
      else
        {
          for (auto ival = first; ival >= last; ival += inc) v.push_back(ival);
        }
    }

  return v;
}

std::vector<double>
cdo_argv_to_flt(const std::vector<std::string> &argv)
{
  std::vector<double> v;

  for (const auto &argument : argv)
    {
      auto len = (int) argument.size();
      int i;
      for (i = 0; i < len; ++i)
        if (argument[i] != '/' && argument[i] != '-' && !isdigit(argument[i])) break;

      if (i != len) { v.push_back(parameter_to_double(argument)); }
      else
        {
          int first, last, inc;
          split_intstring(argument, first, last, inc);

          if (inc >= 0)
            {
              for (int ival = first; ival <= last; ival += inc) v.push_back((double) ival);
            }
          else
            {
              for (int ival = first; ival >= last; ival += inc) v.push_back((double) ival);
            }
        }
    }

  return v;
}

static void
split_intstring(const char *const intstr, int &first, int &last, int &inc)
{
  auto startptr = intstr;
  char *endptr = nullptr;
  auto ival = (int) strtol(startptr, &endptr, 10);
  if (*endptr != 0 && *endptr != '/' && (endptr - startptr) == 0) parameter_error("Integer", startptr, endptr);
  first = ival;
  last = ival;
  inc = 1;

  if (*endptr == '/')
    {
      startptr = endptr + 1;
      endptr = nullptr;
      ival = (int) strtol(startptr, &endptr, 10);
      if (*endptr != 0 && *endptr != '/' && (endptr - startptr) == 0) parameter_error("Integer", startptr, endptr);
      last = ival;
      if (first > last) inc = -1;

      if (*endptr == '/')
        {
          startptr = endptr + 1;
          endptr = nullptr;
          ival = (int) strtol(startptr, &endptr, 10);
          if (*endptr != 0 && (endptr - startptr) == 0) parameter_error("Integer", startptr, endptr);
          inc = ival;
        }
    }
}

void
split_intstring(const std::string &intstr, int &first, int &last, int &inc)
{
  split_intstring(intstr.c_str(), first, last, inc);
}
