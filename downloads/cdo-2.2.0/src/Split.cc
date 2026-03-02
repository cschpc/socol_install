/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Split      splitcode       Split codes
      Split      splitparam      Split parameters
      Split      splitname       Split variables
      Split      splitlevel      Split levels
      Split      splitgrid       Split grids
      Split      splitzaxis      Split zaxis
      Split      splittabnum     Split table numbers
*/

#include <cdi.h>

#include "process_int.h"
#include "cdo_history.h"
#include "util_files.h"
#include "cdo_zaxis.h"
#include "cdi_lockedIO.h"

#include <cassert>

static void
gen_filename(char *filename, bool swapObase, const std::string &obase, const char *suffix)
{
  if (swapObase) strcat(filename, obase.c_str());
  if (suffix[0]) strcat(filename, suffix);
}

static int
split_code(bool swapObase, char *filesuffix, char *filename, int vlistID1, const VarList &varList1,
           std::vector<CdoStreamID> &streamIDs, std::vector<int> &vlistIDs)
{
  auto nchars = strlen(filename);
  auto nvars = vlistNvars(vlistID1);
  std::vector<int> codes(nvars);

  int nsplit = 0;
  for (int varID = 0; varID < nvars; ++varID)
    {
      const auto &var = varList1[varID];
      int index;
      for (index = 0; index < varID; ++index)
        if (var.code == varList1[index].code) break;

      if (index == varID) codes[nsplit++] = var.code;
    }

  vlistIDs.resize(nsplit);
  streamIDs.resize(nsplit);

  for (int index = 0; index < nsplit; ++index)
    {
      vlistClearFlag(vlistID1);
      for (int varID = 0; varID < nvars; ++varID)
        {
          const auto &var = varList1[varID];
          if (codes[index] == var.code)
            {
              for (int levelID = 0; levelID < var.nlevels; ++levelID)
                {
                  vlistDefIndex(vlistID1, varID, levelID, index);
                  vlistDefFlag(vlistID1, varID, levelID, true);
                }
            }
        }

      auto vlistID2 = vlistCreate();
      cdo_vlist_copy_flag(vlistID2, vlistID1);
      vlistIDs[index] = vlistID2;

      const char *format = (codes[index] > 9999) ? "%05d" : ((codes[index] > 999) ? "%04d" : "%03d");
      sprintf(filename + nchars, format, codes[index]);
      gen_filename(filename, swapObase, cdo_get_obase(), filesuffix);

      streamIDs[index] = cdo_open_write(filename);
    }

  return nsplit;
}

static int
split_param(bool swapObase, char *filesuffix, char *filename, int vlistID1, const VarList &varList1,
            std::vector<CdoStreamID> &streamIDs, std::vector<int> &vlistIDs)
{
  auto nchars = strlen(filename);
  char paramstr[32];
  auto nvars = vlistNvars(vlistID1);
  std::vector<int> params(nvars);

  int nsplit = 0;
  for (int varID = 0; varID < nvars; ++varID)
    {
      const auto &var = varList1[varID];
      int index;
      for (index = 0; index < varID; ++index)
        if (var.param == varList1[index].param) break;

      if (index == varID) params[nsplit++] = var.param;
    }

  vlistIDs.resize(nsplit);
  streamIDs.resize(nsplit);

  for (int index = 0; index < nsplit; ++index)
    {
      vlistClearFlag(vlistID1);
      for (int varID = 0; varID < nvars; ++varID)
        {
          const auto &var = varList1[varID];
          if (params[index] == var.param)
            {
              for (int levelID = 0; levelID < var.nlevels; ++levelID)
                {
                  vlistDefIndex(vlistID1, varID, levelID, index);
                  vlistDefFlag(vlistID1, varID, levelID, true);
                }
            }
        }

      auto vlistID2 = vlistCreate();
      cdo_vlist_copy_flag(vlistID2, vlistID1);
      vlistIDs[index] = vlistID2;

      cdiParamToString(params[index], paramstr, sizeof(paramstr));

      filename[nchars] = '\0';
      strcat(filename, paramstr);
      gen_filename(filename, swapObase, cdo_get_obase(), filesuffix);

      streamIDs[index] = cdo_open_write(filename);
    }

  return nsplit;
}

static int
split_name(bool swapObase, char *filesuffix, char *filename, int vlistID1, const VarList &varList1,
           std::vector<CdoStreamID> &streamIDs, std::vector<int> &vlistIDs)
{
  auto nchars = strlen(filename);
  auto nvars = vlistNvars(vlistID1);
  auto nsplit = nvars;

  vlistIDs.resize(nsplit);
  streamIDs.resize(nsplit);

  for (int index = 0; index < nsplit; ++index)
    {
      vlistClearFlag(vlistID1);
      int varID = index;
      const auto &var = varList1[varID];
      for (int levelID = 0; levelID < var.nlevels; ++levelID)
        {
          vlistDefIndex(vlistID1, varID, levelID, index);
          vlistDefFlag(vlistID1, varID, levelID, true);
        }

      auto vlistID2 = vlistCreate();
      cdo_vlist_copy_flag(vlistID2, vlistID1);
      vlistIDs[index] = vlistID2;

      filename[nchars] = '\0';
      strcat(filename, var.name.c_str());
      gen_filename(filename, swapObase, cdo_get_obase(), filesuffix);

      streamIDs[index] = cdo_open_write(filename);
    }

  return nsplit;
}

static int
split_level(bool swapObase, char *filesuffix, char *filename, int vlistID1, const VarList &varList1,
            std::vector<CdoStreamID> &streamIDs, std::vector<int> &vlistIDs)
{
  auto nchars = strlen(filename);
  auto nvars = vlistNvars(vlistID1);
  auto nzaxis = vlistNzaxis(vlistID1);
  double ftmp[999];

  int nsplit = 0;
  for (int index = 0; index < nzaxis; ++index)
    {
      auto zaxisID = vlistZaxis(vlistID1, index);
      auto nlevels = zaxisInqSize(zaxisID);
      for (int levelID = 0; levelID < nlevels; ++levelID)
        {
          auto level = cdo_zaxis_inq_level(zaxisID, levelID);
          int i;
          for (i = 0; i < nsplit; ++i)
            if (is_equal(level, ftmp[i])) break;
          if (i == nsplit) ftmp[nsplit++] = level;
        }
    }

  vlistIDs.resize(nsplit);
  streamIDs.resize(nsplit);
  Varray<double> levels(nsplit);
  for (int index = 0; index < nsplit; ++index) levels[index] = ftmp[index];

  for (int index = 0; index < nsplit; ++index)
    {
      vlistClearFlag(vlistID1);
      for (int varID = 0; varID < nvars; ++varID)
        {
          const auto &var = varList1[varID];
          for (int levelID = 0; levelID < var.nlevels; ++levelID)
            {
              auto level = cdo_zaxis_inq_level(var.zaxisID, levelID);
              if (is_equal(levels[index], level))
                {
                  vlistDefIndex(vlistID1, varID, levelID, index);
                  vlistDefFlag(vlistID1, varID, levelID, true);
                }
            }
        }

      auto vlistID2 = vlistCreate();
      cdo_vlist_copy_flag(vlistID2, vlistID1);
      vlistIDs[index] = vlistID2;

      sprintf(filename + nchars, "%06g", levels[index]);
      gen_filename(filename, swapObase, cdo_get_obase(), filesuffix);

      streamIDs[index] = cdo_open_write(filename);
    }

  return nsplit;
}

static int
split_grid(bool swapObase, char *filesuffix, char *filename, int vlistID1, const VarList &varList1,
           std::vector<CdoStreamID> &streamIDs, std::vector<int> &vlistIDs)
{
  auto nchars = strlen(filename);
  auto nsplit = vlistNgrids(vlistID1);
  auto nvars = vlistNvars(vlistID1);

  vlistIDs.resize(nsplit);
  streamIDs.resize(nsplit);
  std::vector<int> gridIDs(nsplit);
  for (int index = 0; index < nsplit; ++index) gridIDs[index] = vlistGrid(vlistID1, index);

  for (int index = 0; index < nsplit; ++index)
    {
      vlistClearFlag(vlistID1);
      for (int varID = 0; varID < nvars; ++varID)
        {
          const auto &var = varList1[varID];
          if (gridIDs[index] == var.gridID)
            {
              for (int levelID = 0; levelID < var.nlevels; ++levelID)
                {
                  vlistDefIndex(vlistID1, varID, levelID, index);
                  vlistDefFlag(vlistID1, varID, levelID, true);
                }
            }
        }

      auto vlistID2 = vlistCreate();
      cdo_vlist_copy_flag(vlistID2, vlistID1);
      vlistIDs[index] = vlistID2;

      sprintf(filename + nchars, "%02d", vlistGridIndex(vlistID1, gridIDs[index]) + 1);
      gen_filename(filename, swapObase, cdo_get_obase(), filesuffix);

      streamIDs[index] = cdo_open_write(filename);
    }

  return nsplit;
}

static int
split_zaxis(bool swapObase, char *filesuffix, char *filename, int vlistID1, const VarList &varList1,
            std::vector<CdoStreamID> &streamIDs, std::vector<int> &vlistIDs)
{
  auto nchars = strlen(filename);
  auto nsplit = vlistNzaxis(vlistID1);
  auto nvars = vlistNvars(vlistID1);

  vlistIDs.resize(nsplit);
  streamIDs.resize(nsplit);
  std::vector<int> zaxisIDs(nsplit);
  for (int index = 0; index < nsplit; ++index) zaxisIDs[index] = vlistZaxis(vlistID1, index);

  for (int index = 0; index < nsplit; ++index)
    {
      vlistClearFlag(vlistID1);
      for (int varID = 0; varID < nvars; ++varID)
        {
          const auto &var = varList1[varID];
          if (zaxisIDs[index] == var.zaxisID)
            {
              for (int levelID = 0; levelID < var.nlevels; ++levelID)
                {
                  vlistDefIndex(vlistID1, varID, levelID, index);
                  vlistDefFlag(vlistID1, varID, levelID, true);
                }
            }
        }

      auto vlistID2 = vlistCreate();
      cdo_vlist_copy_flag(vlistID2, vlistID1);
      vlistIDs[index] = vlistID2;

      sprintf(filename + nchars, "%02d", vlistZaxisIndex(vlistID1, zaxisIDs[index]) + 1);
      gen_filename(filename, swapObase, cdo_get_obase(), filesuffix);

      streamIDs[index] = cdo_open_write(filename);
    }

  return nsplit;
}

static int
split_tabnum(bool swapObase, char *filesuffix, char *filename, int vlistID1, const VarList &varList1,
             std::vector<CdoStreamID> &streamIDs, std::vector<int> &vlistIDs)
{
  auto nchars = strlen(filename);
  auto nvars = vlistNvars(vlistID1);
  std::vector<int> tabnums(nvars);

  int nsplit = 0;
  for (int varID = 0; varID < nvars; ++varID)
    {
      auto tabnum = tableInqNum(vlistInqVarTable(vlistID1, varID));
      int index;
      for (index = 0; index < varID; ++index)
        if (tabnum == tableInqNum(vlistInqVarTable(vlistID1, index))) break;

      if (index == varID) tabnums[nsplit++] = tabnum;
    }

  vlistIDs.resize(nsplit);
  streamIDs.resize(nsplit);

  for (int index = 0; index < nsplit; ++index)
    {
      vlistClearFlag(vlistID1);
      for (int varID = 0; varID < nvars; ++varID)
        {
          auto tabnum = tableInqNum(vlistInqVarTable(vlistID1, varID));
          if (tabnums[index] == tabnum)
            {
              const auto &var = varList1[varID];
              for (int levelID = 0; levelID < var.nlevels; ++levelID)
                {
                  vlistDefIndex(vlistID1, varID, levelID, index);
                  vlistDefFlag(vlistID1, varID, levelID, true);
                }
            }
        }

      auto vlistID2 = vlistCreate();
      cdo_vlist_copy_flag(vlistID2, vlistID1);
      vlistIDs[index] = vlistID2;

      sprintf(filename + nchars, "%03d", tabnums[index]);
      gen_filename(filename, swapObase, cdo_get_obase(), filesuffix);

      streamIDs[index] = cdo_open_write(filename);
    }

  return nsplit;
}

void *
Split(void *process)
{
  cdo_initialize(process);

  if (process_self().m_ID != 0) cdo_abort("This operator can't be combined with other operators!");

  auto dataIsUnchanged = data_is_unchanged();

  // clang-format off
  auto SPLITCODE   = cdo_operator_add("splitcode",   0, 0, nullptr);
  auto SPLITPARAM  = cdo_operator_add("splitparam",  0, 0, nullptr);
  auto SPLITNAME   = cdo_operator_add("splitname",   0, 0, nullptr);
  auto SPLITLEVEL  = cdo_operator_add("splitlevel",  0, 0, nullptr);
  auto SPLITGRID   = cdo_operator_add("splitgrid",   0, 0, nullptr);
  auto SPLITZAXIS  = cdo_operator_add("splitzaxis",  0, 0, nullptr);
  auto SPLITTABNUM = cdo_operator_add("splittabnum", 0, 0, nullptr);
  // clang-format on

  auto operatorID = cdo_operator_id();

  auto swapObase = false;
  const char *uuidAttribute = nullptr;
  for (int i = 0; i < cdo_operator_argc(); ++i)
    {
      if (cdo_operator_argv(i) == "swap")
        swapObase = true;
      else if (cdo_operator_argv(i).find("uuid=") == 0)
        uuidAttribute = &cdo_operator_argv(i)[0] + 5;
      else
        cdo_abort("Unknown parameter: >%s<", cdo_operator_argv(0));
    }

  auto streamID1 = cdo_open_read(0);
  auto vlistID1 = cdo_stream_inq_vlist(streamID1);

  VarList varList1;
  varListInit(varList1, vlistID1);

  char filename[8192] = { 0 };
  if (!swapObase) strcpy(filename, cdo_get_obase().c_str());

  char filesuffix[32] = { 0 };
  FileUtils::gen_suffix(filesuffix, sizeof(filesuffix), cdo_inq_filetype(streamID1), vlistID1, cdo_get_stream_name(0));

  std::vector<int> vlistIDs;
  std::vector<CdoStreamID> streamIDs;

  int nsplit = 0;
  if (operatorID == SPLITCODE) { nsplit = split_code(swapObase, filesuffix, filename, vlistID1, varList1, streamIDs, vlistIDs); }
  else if (operatorID == SPLITPARAM)
    {
      nsplit = split_param(swapObase, filesuffix, filename, vlistID1, varList1, streamIDs, vlistIDs);
    }
  else if (operatorID == SPLITTABNUM)
    {
      nsplit = split_tabnum(swapObase, filesuffix, filename, vlistID1, varList1, streamIDs, vlistIDs);
    }
  else if (operatorID == SPLITNAME)
    {
      nsplit = split_name(swapObase, filesuffix, filename, vlistID1, varList1, streamIDs, vlistIDs);
    }
  else if (operatorID == SPLITLEVEL)
    {
      nsplit = split_level(swapObase, filesuffix, filename, vlistID1, varList1, streamIDs, vlistIDs);
    }
  else if (operatorID == SPLITGRID)
    {
      nsplit = split_grid(swapObase, filesuffix, filename, vlistID1, varList1, streamIDs, vlistIDs);
    }
  else if (operatorID == SPLITZAXIS)
    {
      nsplit = split_zaxis(swapObase, filesuffix, filename, vlistID1, varList1, streamIDs, vlistIDs);
    }
  else { cdo_abort("not implemented!"); }

  assert(nsplit > 0);

  for (int index = 0; index < nsplit; ++index)
    {
      if (uuidAttribute) cdo_def_tracking_id(vlistIDs[index], uuidAttribute);

      cdo_def_vlist(streamIDs[index], vlistIDs[index]);
    }

  Field field;

  int tsID = 0;
  while (true)
    {
      auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
      if (nrecs == 0) break;

      for (int index = 0; index < nsplit; ++index) cdo_def_timestep(streamIDs[index], tsID);

      for (int recID = 0; recID < nrecs; ++recID)
        {
          int varID, levelID;
          cdo_inq_record(streamID1, &varID, &levelID);

          auto index = vlistInqIndex(vlistID1, varID, levelID);
          auto vlistID2 = vlistIDs[index];
          auto varID2 = vlistFindVar(vlistID2, varID);
          auto levelID2 = vlistFindLevel(vlistID2, varID, levelID);
          // printf("%d %d %d %d %d %d\n", index, vlistID2, varID, levelID, varID2, levelID2);

          cdo_def_record(streamIDs[index], varID2, levelID2);
          if (dataIsUnchanged) { cdo_copy_record(streamIDs[index], streamID1); }
          else
            {
              const auto &var = varList1[varID];
              field.init(var);
              cdo_read_record(streamID1, field);
              cdo_write_record(streamIDs[index], field);
            }
        }

      tsID++;
    }

  cdo_stream_close(streamID1);

  for (auto &streamID : streamIDs) cdo_stream_close(streamID);
  for (auto &vlistID : vlistIDs) vlistDestroy(vlistID);

  cdo_finish();

  return nullptr;
}
