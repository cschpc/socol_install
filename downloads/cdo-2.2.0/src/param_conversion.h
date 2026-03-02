#ifndef PARAM_CONVERSION_H
#define PARAM_CONVERSION_H

#include <string>
#include <vector>

long parameter_to_bytes(const std::string &string);

const char *parameter_to_word(const char *cstring);
double parameter_to_double(const char *cstring);
int parameter_to_int(const char *cstring);
long parameter_to_long(const char *cstring);
size_t parameter_to_size_t(const char *cstring);
int parameter_to_intlist(const char *cstring);

const std::string &parameter_to_word(const std::string &string);
double parameter_to_double(const std::string &string);
bool parameter_to_bool(const std::string &string);
int parameter_to_int(const std::string &string);
long parameter_to_long(const std::string &string);
size_t parameter_to_size_t(const std::string &string);
int parameter_to_intlist(const std::string &string);

double radius_str_to_deg(const std::string &string);

int string_to_param(const std::string &paramstr);
int string_to_param(const char *paramstr);
void param_to_string(int param, char *paramstr, int maxlen);

/* time/date/season converisons */
/* =================================================================================== */
void season_to_months(const char *season, int *imonths);
double date_str_to_double(const char *datestr, int opt);

/* argv conversions */
std::vector<int> cdo_argv_to_int(const std::vector<std::string> &argv);
std::vector<double> cdo_argv_to_flt(const std::vector<std::string> &argv);

void split_intstring(const std::string &intstr, int &first, int &last, int &inc);

#endif
