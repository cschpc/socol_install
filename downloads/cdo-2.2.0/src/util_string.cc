/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida
          Oliver Heidmann

*/

#include <cassert>
#include <cstring>
#include <cctype>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <regex>

#include <cdi.h>
#include "util_string.h"

std::vector<std::string>
split_string(const std::string &str, const std::string &delimiter)
{
  std::regex regex(delimiter);
  return { std::sregex_token_iterator(str.begin(), str.end(), regex, -1), std::sregex_token_iterator() };
}

std::string
string_to_upper(std::string str)
{
  std::transform(str.begin(), str.end(), str.begin(), [](unsigned char c) { return std::toupper(c); });
  return str;
}

std::string
string_to_lower(std::string str)
{
  std::transform(str.begin(), str.end(), str.begin(), [](unsigned char c) { return std::tolower(c); });
  return str;
}

void
cstr_to_lower(char *cstr)
{
  if (cstr)
    for (size_t i = 0; cstr[i]; ++i) cstr[i] = (char) tolower((int) cstr[i]);
}

void
cstr_to_upper(char *cstr)
{
  if (cstr)
    for (size_t i = 0; cstr[i]; ++i) cstr[i] = (char) toupper((int) cstr[i]);
}

static void
trim_flt(char *ss)
{
  char *cp = ss;
  if (*cp == '-') cp++;
  while (isdigit((int) *cp) || *cp == '.') cp++;
  if (*--cp == '.') return;

  char *ep = cp + 1;
  while (*cp == '0') cp--;
  cp++;
  if (cp == ep) return;
  while (*ep) *cp++ = *ep++;
  *cp = '\0';

  return;
}

std::string
get_scientific(double p_float_string)
{
  std::stringstream s;
  s << std::defaultfloat << p_float_string;
  return s.str();
}

char *
double_to_att_str(int digits, char *attstr, size_t len, double value)
{
  const auto ret = std::snprintf(attstr, len, "%#.*g", digits, value);
  assert(ret != -1 && ret < (int) len);
  trim_flt(attstr);
  return attstr;
}

const char *
tunit_to_cstr(int tunits)
{
  // clang-format off
  if      (tunits == TUNIT_YEAR)       return "years";
  else if (tunits == TUNIT_MONTH)      return "months";
  else if (tunits == TUNIT_DAY)        return "days";
  else if (tunits == TUNIT_12HOURS)    return "12hours";
  else if (tunits == TUNIT_6HOURS)     return "6hours";
  else if (tunits == TUNIT_3HOURS)     return "3hours";
  else if (tunits == TUNIT_HOUR)       return "hours";
  else if (tunits == TUNIT_30MINUTES)  return "30minutes";
  else if (tunits == TUNIT_QUARTER)    return "15minutes";
  else if (tunits == TUNIT_MINUTE)     return "minutes";
  else if (tunits == TUNIT_SECOND)     return "seconds";
  else                                 return "unknown";
  // clang-format on
}

const char *
calendar_to_cstr(int calendar)
{
  // clang-format off
  if      (calendar == CALENDAR_STANDARD)  return "standard";
  else if (calendar == CALENDAR_GREGORIAN) return "gregorian";
  else if (calendar == CALENDAR_PROLEPTIC) return "proleptic_gregorian";
  else if (calendar == CALENDAR_360DAYS)   return "360_day";
  else if (calendar == CALENDAR_365DAYS)   return "365_day";
  else if (calendar == CALENDAR_366DAYS)   return "366_day";
  else                                     return "unknown";
  // clang-format on
}

std::vector<std::string>
split_with_seperator(const std::string &sourceString, const char seperator)
{
  std::stringstream ss(sourceString);
  std::string token;
  std::vector<std::string> tokens;
  while (std::getline(ss, token, seperator)) { tokens.push_back(token); }
  return tokens;
}

bool
string_is_float(const std::string &str)
{
  if (str.empty() || isspace(str[0])) return 0;

  char *ptr = nullptr;
  strtod(str.c_str(), &ptr);
  return (*ptr == '\0');
}

bool
string_is_int(const std::string &str)
{
  if (str.empty() || isspace(str[0])) return 0;

  char *ptr = nullptr;
  strtol(str.c_str(), &ptr, 10);
  return (*ptr == '\0');
}

#include <iostream>
/** Tokenizes a string with  delimiter ',' and returns a touple where first
 * signals success and the second the result.
 * If any of the tokens is not a integer the result is negative.
 * With a negative result a empty
 * list is returned with the result boolean.
 **/
std::tuple<bool, std::vector<std::string>>
tokenize_comma_seperated_int_list(const std::string &args)
{
  auto tokens = split_with_seperator(args, ',');
  bool res = true;
  for (auto t : tokens)
    {
      if (!string_is_int(t))
        {
          res = false;
          tokens = {};
          break;
        }
    }
  return { res, tokens };
}

// To replace a single char with another single char in a given string

void
cstr_replace_char(char *str_in, char orig_char, char rep_char)
{
  if (strchr(str_in, orig_char) == nullptr) return;

  while (*str_in != '\0')
    {
      if (*str_in == orig_char) *str_in = rep_char;
      str_in++;
    }

  return;
}
