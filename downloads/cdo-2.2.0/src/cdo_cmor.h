#ifndef CDO_CMOR_H
#define CDO_CMOR_H

#include <string>
#include <vector>
#include <cdi.h>

struct CmorVar
{
  bool convert = false;
  bool remove = false;
  // missing value
  bool changemissval = false;
  double missval_old = 0.0;
  //
  bool lfactor = false;
  double factor = 0.0;
  //
  bool checkvalid = false;
  double valid_min = 0.0;
  double valid_max = 0.0;
  //
  bool check_min_mean_abs = false;
  double ok_min_mean_abs = 0.0;
  //
  bool check_max_mean_abs = false;
  double ok_max_mean_abs = 0.0;
  // units
  bool changeunits = false;
  char units_old[CDI_MAX_NAME] = { 0 };
  char units[CDI_MAX_NAME] = { 0 };
  // varname
  std::string name;
  // converter
  void *ut_converter = nullptr;

  double amean = 0;
  long nvals = 0, n_lower_min = 0, n_greater_max = 0;
};

void cmor_check_init(int nvars, std::vector<CmorVar> &vars);
void cmor_check_eval(int vlistID, int nvars, const std::vector<CmorVar> &vars);
void cmor_check_prep(CmorVar &var, const long gridsize, const double missval, const double *const array);

#endif
