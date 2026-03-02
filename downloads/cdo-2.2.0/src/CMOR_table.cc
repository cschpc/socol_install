/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

*/

#include "cdo_options.h"
#include "process_int.h"
#include "pmlist.h"

static void
dump_cmor_table(const PMList &pmlist)
{
  printf("# Number of lists: %zu\n", pmlist.size());
  int i = 0;
  for (const auto &kvlist : pmlist)
    {
      printf("# list ID: %d;   Number of elements: %zu\n", i, kvlist.size());
      printf("&%s\n", kvlist.name.c_str());
      for (const auto &kv : kvlist) { printf("  %s = %s\n", kv.key.c_str(), kv.values[0].c_str()); }
      printf("/\n");
      ++i;
    }
}

static void
conv_cmor_table(const PMList &pmlist)
{
  const char *hname = "Header";
  const char *vname = "variable";
  // const char *aname = "axis";

  bool hasmissval = false;
  double missval = 0;

  for (const auto &kvlist : pmlist)
    {
      const auto listname = kvlist.name.c_str();

      if (strncmp(listname, hname, strlen(hname)) == 0)
        {
          for (const auto &kv : kvlist)
            {
              const auto ename = kv.key.c_str();
              const auto evalue = kv.values[0].c_str();
              const size_t len = strlen(ename);

              if (strncmp("missing_value", ename, len) == 0)
                {
                  missval = atof(evalue);
                  hasmissval = true;
                }
            }
        }
      else if (strncmp(listname, vname, strlen(vname)) == 0)
        {
          printf("&%s\n", "parameter");
          for (const auto &kv : kvlist)
            {
              const char *ename = kv.key.c_str();
              const char *evalue = kv.values[0].c_str();
              const int len = strlen(ename);
              int vlen = strlen(evalue);

              if (vlen > 1 && evalue[0] == '"' && evalue[vlen - 1] == '"')
                {
                  vlen -= 2;
                  evalue++;
                }

              char *ovalue = strdup(evalue);
              for (int i = 1; i < vlen; ++i)
                {
                  if (ovalue[i - 1] == '"' && ovalue[i] == '"')
                    {
                      ovalue[i - 1] = '\'';
                      for (int j = i + 1; j < vlen; ++j) ovalue[j - 1] = ovalue[j];
                      vlen -= 1;
                    }
                }

              if (vlen)
                {
                  if (strncmp("name", ename, len) == 0 || strncmp("standard_name", ename, len) == 0
                      || strncmp("out_name", ename, len) == 0 || strncmp("type", ename, len) == 0
                      || strncmp("valid_min", ename, len) == 0 || strncmp("valid_max", ename, len) == 0
                      || strncmp("ok_min_mean_abs", ename, len) == 0 || strncmp("ok_max_mean_abs", ename, len) == 0)
                    printf("  %-15s = %s\n", ename, ovalue);
                  else if (strncmp("long_name", ename, len) == 0 || strncmp("units", ename, len) == 0
                           || strncmp("cell_methods", ename, len) == 0 || strncmp("cell_measures", ename, len) == 0
                           || strncmp("comment", ename, len) == 0)
                    printf("  %-15s = \"%.*s\"\n", ename, vlen, ovalue);
                }

              free(ovalue);
            }
          if (hasmissval) printf("  %-15s = %g\n", "missing_value", missval);
          printf("/\n");
        }
    }
}

void *
CMOR_table(void *process)
{
  cdo_initialize(process);

  const auto DUMP_CMOR_TABLE = cdo_operator_add("dump_cmor_table", 0, 0, nullptr);
  const auto CONV_CMOR_TABLE = cdo_operator_add("conv_cmor_table", 0, 0, nullptr);

  const auto operatorID = cdo_operator_id();

  if (cdo_operator_argc() != 1) cdo_abort("Too few arguments!");
  const char *filename = cdo_operator_argv(0).c_str();

  if (Options::cdoVerbose) cdo_print("Parse file: %s", filename);

  auto fp = std::fopen(filename, "r");
  if (fp == nullptr) cdo_abort("Open failed on: %s\n", filename);

  PMList pmlist;
  pmlist.read_cmor_table(fp, filename);
  std::fclose(fp);

  if (operatorID == DUMP_CMOR_TABLE)
    dump_cmor_table(pmlist);
  else if (operatorID == CONV_CMOR_TABLE)
    conv_cmor_table(pmlist);

  cdo_finish();

  return nullptr;
}
