/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Merge      merge           Merge datasets with different fields
*/

#include <cdi.h>

#include "cdo_options.h"
#include "cdo_rlimit.h"
#include "process_int.h"
#include "cdo_zaxis.h"
#include "util_files.h"
#include "cdo_default_values.h"

static void
check_dup_entry(int vlistID1, int vlistID2, const std::string &filename)
{
  Varray<double> lev1, lev2;

  auto nvars1 = vlistNvars(vlistID1);
  auto nvars2 = vlistNvars(vlistID2);

  VarList varList1, varList2;
  varListInit(varList1, vlistID1);
  varListInit(varList2, vlistID2);

  for (int varID1 = 0; varID1 < nvars1; ++varID1)
    {
      const auto &var1 = varList1[varID1];
      auto gtype1 = gridInqType(var1.gridID);
      auto gsize1 = var1.gridsize;
      auto ztype1 = zaxisInqType(var1.zaxisID);
      size_t nlev1 = var1.nlevels;
      if (nlev1 > lev1.size()) lev1.resize(nlev1);
      cdo_zaxis_inq_levels(var1.zaxisID, lev1.data());

      for (int varID2 = 0; varID2 < nvars2; ++varID2)
        {
          const auto &var2 = varList2[varID2];
          auto gtype2 = gridInqType(var2.gridID);
          auto gsize2 = var2.gridsize;
          auto ztype2 = zaxisInqType(var2.zaxisID);
          size_t nlev2 = var2.nlevels;
          if (gtype1 == gtype2 && gsize1 == gsize2 && ztype1 == ztype2 && nlev1 == nlev2)
            {
              if (nlev2 > lev2.size()) lev2.resize(nlev2);
              cdo_zaxis_inq_levels(var2.zaxisID, lev2.data());

              if (zaxisInqLevels(var1.zaxisID, nullptr) && zaxisInqLevels(var2.zaxisID, nullptr))
                {
                  for (size_t k = 0; k < nlev2; ++k)
                    if (is_not_equal(lev1[k], lev2[k])) return;
                }

              if ((var1.param < 0 || var2.param < 0) && var1.name == var2.name)
                {
                  cdo_warning("Duplicate entry of parameter name %s in %s!", var2.name, filename);
                }
              else if (var1.param >= 0 && var2.param >= 0 && var1.param == var2.param && var1.name == var2.name)
                {
                  char paramstr[32];
                  cdiParamToString(var2.param, paramstr, sizeof(paramstr));
                  cdo_warning("Duplicate entry of parameter %s in %s!", paramstr, filename);
                }
              else if (var1.param != var2.param && var1.name == var2.name)
                {
                  cdo_warning("Duplicate entry of parameter name %s with different IDs in %s!", var2.name, filename);
                }
            }
        }
    }
}

static int
get_taxis_index(const std::vector<int> &vlistIDs)
{
  if (vlistNtsteps(vlistIDs[0]) == 0)
    {
      const int nmerge = vlistIDs.size();
      for (int im = 1; im < nmerge; ++im)
        {
          if (vlistNtsteps(vlistIDs[im]) != 0) return im;
        }
    }

  return 0;
}

int
cdo_inq_base_filetype(CdoStreamID streamID)
{
  auto filetype = cdo_inq_filetype(streamID);
  switch (filetype)
    {
    case CDI_FILETYPE_NCZARR:
    case CDI_FILETYPE_NC5:
    case CDI_FILETYPE_NC4C:
    case CDI_FILETYPE_NC4:
    case CDI_FILETYPE_NC2:
    case CDI_FILETYPE_NC: return CDI_FILETYPE_NC;
    case CDI_FILETYPE_GRB2:
    case CDI_FILETYPE_GRB: return CDI_FILETYPE_GRB;
    }

  return filetype;
}

void *
Merge(void *process)
{
  cdo_initialize(process);

  operator_check_argc(0);

  auto dataIsUnchanged = data_is_unchanged();

  auto streamCnt = cdo_stream_cnt();
  auto nmerge = streamCnt - 1;

  cdo::set_numfiles(nmerge + 8);

  auto ofilename = cdo_get_stream_name(streamCnt - 1);

  if (!Options::cdoOverwriteMode && FileUtils::file_exists(ofilename) && !FileUtils::user_file_overwrite(ofilename))
    cdo_abort("Outputfile %s already exists!", ofilename);

  std::vector<CdoStreamID> streamIDs(nmerge);
  std::vector<int> vlistIDs(nmerge), numRecords(nmerge), numTimeSteps(nmerge);

  auto setFileType = (CdoDefault::FileType == CDI_UNDEFID);
  auto baseFiletype = -1;
  for (int im = 0; im < nmerge; ++im)
    {
      streamIDs[im] = cdo_open_read(im);
      vlistIDs[im] = cdo_stream_inq_vlist(streamIDs[im]);
      if (im == 0) baseFiletype = cdo_inq_base_filetype(streamIDs[im]);
      if (baseFiletype != cdo_inq_base_filetype(streamIDs[im])) dataIsUnchanged = false;
    }

  if (setFileType) CdoDefault::FileType = cdo_inq_filetype(streamIDs[0]);
  // printf("setFileType=%d  filetype=%d\n", setFileType, CdoDefault::FileType);

  auto taxisindex = get_taxis_index(vlistIDs);

  auto taxisID1 = vlistInqTaxis(vlistIDs[taxisindex]);
  auto taxisID2 = taxisDuplicate(taxisID1);

  auto vlistID2 = vlistCreate();
  vlistCopy(vlistID2, vlistIDs[0]);
  for (int im = 1; im < nmerge; ++im)
    {
      std::string commandString = cdo_get_command_from_in_stream(im);
      check_dup_entry(vlistID2, vlistIDs[im], commandString);
      vlistMerge(vlistID2, vlistIDs[im]);
    }

  int numConstVars = 0;
  for (int im = 0; im < nmerge; ++im)
    {
      numTimeSteps[im] = vlistNtsteps(vlistIDs[im]);
      if (numTimeSteps[im] == 0) numTimeSteps[im] = 1;
      if (numTimeSteps[im] == 1) numConstVars++;
    }

  if (numConstVars > 0 && numConstVars < nmerge)
    for (int im = 0; im < nmerge; ++im)
      {
        if (numTimeSteps[im] == 1)
          {
            auto vlistID1 = vlistIDs[im];
            auto nvars = vlistNvars(vlistID1);
            for (int varID = 0; varID < nvars; ++varID)
              vlistDefVarTimetype(vlistID2, vlistMergedVar(vlistID1, varID), TIME_CONSTANT);
          }
      }

  if (Options::cdoVerbose)
    {
      for (int im = 0; im < nmerge; ++im) vlistPrint(vlistIDs[im]);
      vlistPrint(vlistID2);
    }

  auto streamID2 = cdo_open_write(streamCnt - 1);

  vlistDefTaxis(vlistID2, taxisID2);
  cdo_def_vlist(streamID2, vlistID2);

  VarList varList2;
  varListInit(varList2, vlistID2);

  Field field;

  int tsID = 0;
  while (true)
    {
      for (int im = 0; im < nmerge; ++im)
        {
          if (vlistIDs[im] == -1) continue;
          numRecords[im] = cdo_stream_inq_timestep(streamIDs[im], tsID);
        }

      {
        int im;
        for (im = 0; im < nmerge; ++im)
          if (numRecords[im] != 0) break;
        if (im == nmerge) break;  // EOF on all input streams
      }

      if (tsID == 1)
        {
          for (int im = 0; im < nmerge; ++im)
            if (numRecords[im] == 0 && numTimeSteps[im] == 1) vlistIDs[im] = -1;
        }

      if (numRecords[taxisindex] == 0)
        {
          for (int im = 1; im < nmerge; ++im)
            if (vlistIDs[im] != -1 && numRecords[im] != 0)
              cdo_warning("Input stream %d has %d timestep%s. Stream %d has more timesteps, skipped!", taxisindex + 1, tsID,
                          (tsID == 1) ? "" : "s", im + 1);
          break;
        }
      else
        {
          auto lstop = false;
          for (int im = 1; im < nmerge; ++im)
            if (vlistIDs[im] != -1 && numRecords[im] == 0)
              {
                cdo_warning("Input stream %d has %d timestep%s. Stream %d has more timesteps, skipped!", im + 1, tsID,
                            (tsID == 1) ? "" : "s", taxisindex + 1);
                lstop = true;
                break;
              }
          if (lstop) break;
        }

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_def_timestep(streamID2, tsID);

      for (int im = 0; im < nmerge; ++im)
        {
          auto streamID1 = streamIDs[im];
          auto vlistID1 = vlistIDs[im];
          if (vlistID1 == -1) continue;

          auto nrecs = numRecords[im];
          for (int recID = 0; recID < nrecs; ++recID)
            {
              int varID, levelID;
              cdo_inq_record(streamID1, &varID, &levelID);

              auto varID2 = vlistMergedVar(vlistID1, varID);
              auto levelID2 = vlistMergedLevel(vlistID1, varID, levelID);

              if (Options::cdoVerbose)
                cdo_print("fileID=%d  varID=%d levelID=%d   varID2=%d levelID2=%d", im, varID, levelID, varID2, levelID2);

              cdo_def_record(streamID2, varID2, levelID2);
              if (dataIsUnchanged) { cdo_copy_record(streamID2, streamID1); }
              else
                {
                  field.init(varList2[varID2]);
                  cdo_read_record(streamID1, field);
                  cdo_write_record(streamID2, field);
                }
            }
        }

      tsID++;
    }

  for (auto &streamID : streamIDs) cdo_stream_close(streamID);

  cdo_stream_close(streamID2);

  vlistDestroy(vlistID2);

  cdo_finish();

  return nullptr;
}
