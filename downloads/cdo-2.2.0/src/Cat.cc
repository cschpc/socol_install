/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Copy       cat             Concatenate datasets
*/

#include <cdi.h>

#include "process_int.h"
#include "cdo_wtime.h"
#include "cdo_vlist.h"
#include "timer.h"
#include "util_files.h"
#include "progress.h"
#include "cdo_options.h"

void *
Cat(void *process)
{
  enum class StreamMode
  {
    APPEND,
    CREATE
  };
  auto streamMode = StreamMode::APPEND;
  auto hasConstVars = true;
  int tsID2 = 0;
  CdoStreamID streamID2;
  int vlistID2 = CDI_UNDEFID;
  int taxisID2 = CDI_UNDEFID;
  Field field;

  cdo_initialize(process);

  operator_check_argc(0);

  auto dataIsUnchanged = data_is_unchanged();

  auto streamCnt = cdo_stream_cnt();
  auto nfiles = streamCnt - 1;

  progress::init();

  for (int indf = 0; indf < nfiles; ++indf)
    {
      auto start = Options::cdoVerbose ? cdo_get_wtime() : 0.0;

      auto streamID1 = cdo_open_read(indf);
      auto vlistID1 = cdo_stream_inq_vlist(streamID1);
      auto taxisID1 = vlistInqTaxis(vlistID1);

      VarList varList1;
      varListInit(varList1, vlistID1);

      if (indf == 0)
        {
          auto nvars = vlistNvars(vlistID1);

          auto ntsteps = vlistNtsteps(vlistID1);
          if (ntsteps == 1 && varList_numVaryingVars(varList1) == 0) ntsteps = 0;

          bool file_exists = Options::cdoOverwriteMode ? false : FileUtils::file_exists(cdo_get_stream_name(nfiles));
          if (file_exists)
            {
              streamID2 = cdo_open_append(nfiles);

              vlistID2 = cdo_stream_inq_vlist(streamID2);
              taxisID2 = vlistInqTaxis(vlistID2);

              vlist_compare(vlistID1, vlistID2, CMP_ALL);

              tsID2 = vlistNtsteps(vlistID2);
              if (tsID2 == 0) tsID2 = 1;  // bug fix for time constant data only

              if (ntsteps == 0) hasConstVars = false;
            }
          else
            {
              if (Options::cdoVerbose) cdo_print("Output file doesn't exist, creating: %s", cdo_get_stream_name(nfiles));

              streamMode = StreamMode::CREATE;
              streamID2 = cdo_open_write(nfiles);

              vlistID2 = vlistDuplicate(vlistID1);
              taxisID2 = taxisDuplicate(taxisID1);
              vlistDefTaxis(vlistID2, taxisID2);

              if (ntsteps == 0 && nfiles > 1)
                {
                  hasConstVars = false;
                  for (int varID = 0; varID < nvars; ++varID) vlistDefVarTimetype(vlistID2, varID, TIME_VARYING);
                }

              cdo_def_vlist(streamID2, vlistID2);
            }
        }
      else { vlist_compare(vlistID1, vlistID2, CMP_ALL); }

      auto ntsteps = vlistNtsteps(vlistID1);

      int tsID1 = 0;
      while (true)
        {
          auto nrecs = cdo_stream_inq_timestep(streamID1, tsID1);
          if (nrecs == 0) break;

          const double fstatus = (ntsteps > 1) ? indf + (tsID1 + 1.0) / ntsteps : indf + 1.0;
          if (!Options::cdoVerbose) progress::update(0, 1, fstatus / nfiles);

          cdo_taxis_copy_timestep(taxisID2, taxisID1);
          cdo_def_timestep(streamID2, tsID2);

          for (int recID = 0; recID < nrecs; ++recID)
            {
              int varID, levelID;
              cdo_inq_record(streamID1, &varID, &levelID);

              if (hasConstVars && tsID2 > 0 && tsID1 == 0)
                if (varList1[varID].isConstant) continue;

              cdo_def_record(streamID2, varID, levelID);

              if (dataIsUnchanged) { cdo_copy_record(streamID2, streamID1); }
              else
                {
                  field.init(varList1[varID]);
                  cdo_read_record(streamID1, field);
                  cdo_write_record(streamID2, field);
                }
            }

          tsID1++;
          tsID2++;
        }

      cdo_stream_close(streamID1);

      if (Options::cdoVerbose) cdo_print("Processed file: %s   %.2f seconds", cdo_get_stream_name(indf), cdo_get_wtime() - start);
    }

  cdo_stream_close(streamID2);

  if (streamMode == StreamMode::CREATE)
    {
      vlistDestroy(vlistID2);
      taxisDestroy(taxisID2);
    }

  cdo_finish();

  return nullptr;
}
