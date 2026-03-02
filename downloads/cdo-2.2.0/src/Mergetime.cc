/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Merge      mergetime       Merge datasets sorted by date and time
*/

#include <cdi.h>

#include "cdo_options.h"
#include "cdo_rlimit.h"
#include "process_int.h"
#include "cdo_vlist.h"
#include "util_files.h"
#include "printinfo.h"

struct StreamInfo
{
  CdoStreamID streamID;
  CdiDateTime vDateTime{};
  int vlistID;
  int taxisID;
  int tsID;
  int nrecs;
  VarList varList;
};

bool
getenv_skip_same_time()
{
  const auto envstr = getenv("SKIP_SAME_TIME");
  if (envstr)
    {
      const auto ival = atoi(envstr);
      if (ival == 1)
        {
          if (Options::cdoVerbose) cdo_print("Set SKIP_SAME_TIME to %d", ival);
          return true;
        }
    }

  return false;
}

static void
open_all_files(int nfiles, std::vector<StreamInfo> &streamInfo)
{
  for (int fileID = 0; fileID < nfiles; ++fileID)
    {
      if (Options::cdoVerbose) cdo_print("process: %s", cdo_get_stream_name(fileID));

      auto &si = streamInfo[fileID];
      si.streamID = cdo_open_read(fileID);
      si.vlistID = cdo_stream_inq_vlist(si.streamID);
      si.taxisID = vlistInqTaxis(si.vlistID);
      varListInit(si.varList, si.vlistID);
    }
}

static void
read_first_timestep(int nfiles, std::vector<StreamInfo> &streamInfo)
{
  for (int fileID = 0; fileID < nfiles; ++fileID)
    {
      auto &si = streamInfo[fileID];
      si.tsID = 0;
      si.nrecs = cdo_stream_inq_timestep(si.streamID, si.tsID);
      if (si.nrecs == 0)
        {
          cdo_stream_close(si.streamID);
          si.streamID = CDO_STREAM_UNDEF;
        }
      else { si.vDateTime = taxisInqVdatetime(si.taxisID); }
    }
}

void *
Mergetime(void *process)
{
  int tsID2 = 0;
  int taxisID2 = CDI_UNDEFID;
  CdiDateTime lastDateTime{};

  cdo_initialize(process);

  operator_check_argc(0);

  const auto skipSameTime = getenv_skip_same_time();

  const auto dataIsUnchanged = data_is_unchanged();

  const auto nfiles = cdo_stream_cnt() - 1;
  std::vector<StreamInfo> streamInfo(nfiles);

  cdo::set_numfiles(nfiles + 8);

  open_all_files(nfiles, streamInfo);

  // check that the contents is always the same
  for (int fileID = 1; fileID < nfiles; ++fileID) vlist_compare(streamInfo[0].vlistID, streamInfo[fileID].vlistID, CMP_ALL);

  // read the first time step
  read_first_timestep(nfiles, streamInfo);

  auto ofilename = cdo_get_stream_name(nfiles);
  if (!Options::cdoOverwriteMode && FileUtils::file_exists(ofilename) && !FileUtils::user_file_overwrite(ofilename))
    cdo_abort("Outputfile %s already exists!", ofilename);

  const auto streamID2 = cdo_open_write(nfiles);

  Field field;

  while (true)
    {
      auto processTimestep = true;

      int nextFileID = -1;
      CdiDateTime vDateTime{};
      for (int fileID = 0; fileID < nfiles; ++fileID)
        {
          if (streamInfo[fileID].streamID != CDO_STREAM_UNDEF)
            {
              const auto vdate = cdiDate_get(streamInfo[fileID].vDateTime.date);
              const auto vtime = cdiTime_get(streamInfo[fileID].vDateTime.time);
              if (nextFileID == -1 || vdate < cdiDate_get(vDateTime.date)
                  || (vdate == cdiDate_get(vDateTime.date) && vtime < cdiTime_get(vDateTime.time)))
                {
                  nextFileID = fileID;
                  vDateTime = streamInfo[fileID].vDateTime;
                }
            }
        }

      const auto fileID = nextFileID;
      if (Options::cdoVerbose) cdo_print("nextstep = %d  vDateTime = %s", fileID, datetime_to_string(vDateTime));
      if (fileID == -1) break;

      auto &si = streamInfo[fileID];

      if (skipSameTime && cdiDateTime_isEQ(vDateTime, lastDateTime))
        {
          cdo_print("Timestep %4d in stream %d (%s) already exists, skipped!", si.tsID + 1, si.streamID->get_id(),
                    datetime_to_string(vDateTime));
          processTimestep = false;
        }

      if (processTimestep)
        {
          if (tsID2 == 0)
            {
              const auto vlistID1 = si.vlistID;
              const auto vlistID2 = vlistDuplicate(vlistID1);
              const auto taxisID1 = vlistInqTaxis(vlistID1);
              taxisID2 = taxisDuplicate(taxisID1);
              vlistDefTaxis(vlistID2, taxisID2);

              cdo_def_vlist(streamID2, vlistID2);
            }

          lastDateTime = vDateTime;

          cdo_taxis_copy_timestep(taxisID2, si.taxisID);
          cdo_def_timestep(streamID2, tsID2);

          for (int recID = 0; recID < si.nrecs; ++recID)
            {
              int varID, levelID;
              cdo_inq_record(si.streamID, &varID, &levelID);

              if (tsID2 > 0 && si.tsID == 0 && si.varList[varID].isConstant) continue;

              cdo_def_record(streamID2, varID, levelID);

              if (dataIsUnchanged) { cdo_copy_record(streamID2, si.streamID); }
              else
                {
                  field.init(si.varList[varID]);
                  cdo_read_record(si.streamID, field);
                  cdo_write_record(streamID2, field);
                }
            }

          tsID2++;
        }

      si.nrecs = cdo_stream_inq_timestep(si.streamID, ++si.tsID);
      if (si.nrecs == 0)
        {
          cdo_stream_close(si.streamID);
          si.streamID = CDO_STREAM_UNDEF;
        }
      else { si.vDateTime = taxisInqVdatetime(si.taxisID); }
    }

  cdo_stream_close(streamID2);

  cdo_finish();

  return nullptr;
}
