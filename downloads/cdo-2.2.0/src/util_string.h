/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida
          Oliver Heidmann

*/
#ifndef UTIL_STRING_H
#define UTIL_STRING_H

#include <string>
#include <tuple>
#include <vector>

#define ADD_PLURAL(n) ((n) != 1 ? "s" : "")

std::vector<std::string> split_string(const std::string &str, const std::string &regex_str);

std::string string_to_upper(std::string str);
std::string string_to_lower(std::string str);

void cstr_to_lower(char *cstr);
void cstr_to_upper(char *cstr);
char *double_to_att_str(int digits, char *str, size_t len, double value);

const char *tunit_to_cstr(int tunits);
const char *calendar_to_cstr(int calendar);

std::string get_scientific(double p_float_string);

std::vector<std::string> split_with_seperator(const std::string &sourceString, const char seperator);
bool string_is_float(const std::string &str);
bool string_is_int(const std::string &str);
void cstr_replace_char(char *str_in, char orig_char, char rep_char);

std::tuple<bool, std::vector<std::string>> tokenize_comma_seperated_int_list(const std::string &args);

#endif
