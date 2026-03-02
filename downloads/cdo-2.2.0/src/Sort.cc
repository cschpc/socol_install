/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Sort sortcode  Sort by code number
*/

#include <algorithm>  // sort

#include <cdi.h>

#include "cdo_options.h"
#include "cdo_vlist.h"
#include "process_int.h"
#include "param_conversion.h"
#include "cdo_zaxis.h"

struct LevInfo
{
  int levelID;
  size_t nmiss;
  double level;
};

struct VarInfo
{
  int varID;
  int nlevs;
  int code;
  char param[CDI_MAX_NAME];
  char name[CDI_MAX_NAME];
  std::vector<LevInfo> levInfo;
};

static bool
cmpvarcode(const VarInfo &a, const VarInfo &b)
{
  return a.code < b.code;
}

static bool
cmpvarparam(const VarInfo &a, const VarInfo &b)
{
  return strcmp(a.param, b.param) < 0;
}

static bool
cmpvarname(const VarInfo &a, const VarInfo &b)
{
  return strcmp(a.name, b.name) < 0;
}

static bool
cmpvarlevel(const LevInfo &a, const LevInfo &b)
{
  return a.level < b.level;
}

static bool
cmpvarlevelrev(const LevInfo &a, const LevInfo &b)
{
  return a.level > b.level;
}

static void
setNmiss(int varID, int levelID, int nvars, std::vector<VarInfo> &varInfo, size_t nmiss)
{
  int vindex, lindex;

  for (vindex = 0; vindex < nvars; vindex++)
    if (varInfo[vindex].varID == varID) break;

  if (vindex == nvars) cdo_abort("Internal problem; varID not found!");

  auto nlevels = varInfo[vindex].nlevs;
  for (lindex = 0; lindex < nlevels; lindex++)
    if (varInfo[vindex].levInfo[lindex].levelID == levelID) break;

  if (lindex == nlevels) cdo_abort("Internal problem; levelID not found!");

  varInfo[vindex].levInfo[lindex].nmiss = nmiss;
}

void *
Sort(void *process)
{
  bool (*cmpvarlev)(const LevInfo &a, const LevInfo &b) = cmpvarlevel;

  cdo_initialize(process);

  // clang-format off
  auto SORTCODE  = cdo_operator_add("sortcode",  0, 0, nullptr);
  auto SORTPARAM = cdo_operator_add("sortparam", 0, 0, nullptr);
  auto SORTNAME  = cdo_operator_add("sortname",  0, 0, nullptr);
  auto SORTLEVEL = cdo_operator_add("sortlevel", 0, 0, nullptr);
  // clang-format on

  auto operatorID = cdo_operator_id();

  if (cdo_operator_argc() > 1) cdo_abort("Too many arguments!");

  if (operatorID == SORTLEVEL && cdo_operator_argc() == 1)
    {
      auto iarg = parameter_to_int(cdo_operator_argv(0));
      if (iarg < 0) cmpvarlev = cmpvarlevelrev;
    }

  auto streamID1 = cdo_open_read(0);

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto vlistID2 = vlistDuplicate(vlistID1);

  auto taxisID1 = vlistInqTaxis(vlistID1);
  auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);
  /*
  if ( operatorID == SORTCODE )
      vlistSortCode(vlistID2);
   else if ( operatorID == SORTNAME )
      ;
   else if ( operatorID == SORTLEVEL )
      ;
  */

  auto streamID2 = cdo_open_write(1);
  cdo_def_vlist(streamID2, vlistID2);

  VarList varList1;
  varListInit(varList1, vlistID1);

  auto nvars = vlistNvars(vlistID1);

  std::vector<VarInfo> varInfo(nvars);
  for (int varID = 0; varID < nvars; ++varID)
    {
      const auto &var = varList1[varID];
      varInfo[varID].nlevs = var.nlevels;
      varInfo[varID].levInfo.resize(var.nlevels);
    }

  Varray2D<double> vardata(nvars);
  for (int varID = 0; varID < nvars; ++varID)
    {
      const auto &var = varList1[varID];
      vardata[varID].resize(var.gridsize * var.nlevels);
    }

  int tsID = 0;
  while (true)
    {
      auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
      if (nrecs == 0) break;

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_def_timestep(streamID2, tsID);

      for (int recID = 0; recID < nrecs; ++recID)
        {
          int varID, levelID;
          cdo_inq_record(streamID1, &varID, &levelID);
          const auto &var = varList1[varID];

          if (tsID == 0)
            {
              varInfo[varID].varID = varID;
              varInfo[varID].code = var.code;
              param_to_string(var.param, varInfo[varID].param, sizeof(varInfo[varID].param));
              strcpy(varInfo[varID].name, var.name.c_str());
              varInfo[varID].levInfo[levelID].levelID = levelID;
              varInfo[varID].levInfo[levelID].level = cdo_zaxis_inq_level(var.zaxisID, levelID);
            }

          auto offset = var.gridsize * levelID;
          auto single = &vardata[varID][offset];

          size_t nmiss;
          cdo_read_record(streamID1, single, &nmiss);

          setNmiss(varID, levelID, nvars, varInfo, nmiss);
          // varInfo[varID].levInfo[levelID].nmiss = nmiss;
        }

      if (tsID == 0)
        {
          if (Options::cdoVerbose)
            for (int vindex = 0; vindex < nvars; vindex++)
              {
                auto nlevels = varInfo[vindex].nlevs;
                for (int lindex = 0; lindex < nlevels; ++lindex)
                  printf("sort in: %d %s %d %d %g\n", vindex, varInfo[vindex].name, varInfo[vindex].code, varInfo[vindex].nlevs,
                         varInfo[vindex].levInfo[lindex].level);
              }

          if (operatorID == SORTCODE)
            std::sort(varInfo.begin(), varInfo.end(), cmpvarcode);
          else if (operatorID == SORTPARAM)
            std::sort(varInfo.begin(), varInfo.end(), cmpvarparam);
          else if (operatorID == SORTNAME)
            std::sort(varInfo.begin(), varInfo.end(), cmpvarname);
          else if (operatorID == SORTLEVEL)
            {
              for (int vindex = 0; vindex < nvars; vindex++)
                std::sort(varInfo[vindex].levInfo.begin(), varInfo[vindex].levInfo.end(), cmpvarlev);
            }

          if (Options::cdoVerbose)
            for (int vindex = 0; vindex < nvars; vindex++)
              {
                auto nlevels = varInfo[vindex].nlevs;
                for (int lindex = 0; lindex < nlevels; ++lindex)
                  printf("sort out: %d %s %d %d %g\n", vindex, varInfo[vindex].name, varInfo[vindex].code, varInfo[vindex].nlevs,
                         varInfo[vindex].levInfo[lindex].level);
              }
        }

      for (int vindex = 0; vindex < nvars; vindex++)
        {
          auto varID = varInfo[vindex].varID;
          const auto &var = varList1[varID];
          for (int lindex = 0; lindex < var.nlevels; ++lindex)
            {
              auto levelID = varInfo[vindex].levInfo[lindex].levelID;
              auto nmiss = varInfo[vindex].levInfo[lindex].nmiss;

              if (tsID == 0 || var.isConstant)
                {
                  auto offset = var.gridsize * levelID;
                  auto single = &vardata[varID][offset];

                  cdo_def_record(streamID2, varID, levelID);
                  cdo_write_record(streamID2, single, nmiss);
                }
            }
        }

      tsID++;
    }

  cdo_stream_close(streamID1);
  cdo_stream_close(streamID2);

  cdo_finish();

  return nullptr;
}
