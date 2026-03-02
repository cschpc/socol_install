/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include "process_int.h"
#include "util_string.h"
#include "dcw_reader.h"

static void
print_polygons(DCW_Lists &dcw_lists, const std::string &codeNames)
{
  auto codeList = split_string(codeNames, "\\+");

  dcw_sort_countries(dcw_lists);

  printf("# Digital Chart of the World\n");
  printf("# Region for country:");
  for (auto &code : codeList) printf(" %s", code.c_str());
  printf("\n");

  codeList = dcw_expand_code_list(dcw_lists, codeList);

  Region region;
  if (dcw_get_region(dcw_lists, codeList, region)) cdo_abort("dcw_get_region failed!");

  printf("#   West=%g  East=%g  South=%g  North=%g\n", region.west, region.east, region.south, region.north);
  printf("#\n");

  dcw_print_polygons(dcw_lists, codeList);
}

void *
DCW_util(void *process)
{
  cdo_initialize(process);

  DCW_Lists dcw_lists;
  if (dcw_load_lists(dcw_lists)) cdo_abort("dcw_load_lists() failed!");

  if (cdo_operator_argc() == 0) { cdo_abort("Parameter missing (available keywords are path/countries)"); }
  else if (cdo_operator_argc() > 1) { cdo_abort("Too many parameter, max=1!"); }
  else
    {
      const char *param_country = "country=";
      const char *param_dcw = "dcw:";
      auto &param = cdo_operator_argv(0);
      if (param == "path")
        dcw_print_path();
      else if (param == "dir")
        dcw_print_path();
      else if (param == "countries")
        dcw_print_countries(dcw_lists);
      else if (param.compare(0, strlen(param_country), param_country) == 0)
        print_polygons(dcw_lists, param.substr(strlen(param_country)));
      else if (param.compare(0, strlen(param_dcw), param_dcw) == 0)
        print_polygons(dcw_lists, param.substr(strlen(param_dcw)));
      else
        cdo_abort("Unsupported parameter: %s", param);
    }

  cdo_finish();

  return nullptr;
}
