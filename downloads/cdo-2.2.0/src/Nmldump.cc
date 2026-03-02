/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

*/

#include "cdi.h"

#include "cdo_options.h"
#include "process_int.h"
#include "parse_literals.h"
#include "pmlist.h"
#include "util_string.h"

static void
printValues(const int nvalues, const std::vector<std::string> &values)
{
  char fltstr[128];
  if (nvalues)
    {
      const auto dtype = literals_find_datatype(nvalues, values);
      for (int i = 0; i < nvalues; ++i)
        {
          const auto &value = values[i];
          if (i) printf(", ");
          switch (dtype)
            {
            case CDI_DATATYPE_INT8: printf("%db", literal_to_int(value)); break;
            case CDI_DATATYPE_INT16: printf("%ds", literal_to_int(value)); break;
            case CDI_DATATYPE_INT32: printf("%d", literal_to_int(value)); break;
            case CDI_DATATYPE_FLT32:
              printf("%sf", double_to_att_str(Options::CDO_flt_digits, fltstr, sizeof(fltstr), literal_to_double(value)));
              break;
            case CDI_DATATYPE_FLT64:
              printf("%s", double_to_att_str(Options::CDO_dbl_digits, fltstr, sizeof(fltstr), literal_to_double(value)));
              break;
            default: printf("\"%s\"", value.c_str());
            }
        }
    }
}

void
kvldump(const PMList &pmlist)
{
  for (const auto &kvlist : pmlist)
    {
      const auto &listname = kvlist.name;
      if (!listname.empty()) printf("&%s\n", listname.c_str());
      for (const auto &kv : kvlist)
        {
          const auto &key = kv.key;
          if (!listname.empty()) printf("  ");
          printf("%s = ", key.c_str());
          printValues(kv.nvalues, kv.values);
          printf("\n");
        }
      if (!listname.empty()) printf("/\n");
    }
}

void *
Nmldump(void *process)
{
  cdo_initialize(process);

  const auto NMLDUMP = cdo_operator_add("nmldump", 0, 0, nullptr);
  const auto KVLDUMP = cdo_operator_add("kvldump", 0, 0, nullptr);

  const auto operatorID = cdo_operator_id();

  operator_check_argc(0);

  PMList pmlist;
  pmlist.read_namelist(stdin, "STDIN");

  if (operatorID == NMLDUMP)
    pmlist.print();
  else if (operatorID == KVLDUMP)
    kvldump(pmlist);

  cdo_finish();

  return nullptr;
}
