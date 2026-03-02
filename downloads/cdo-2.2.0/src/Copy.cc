/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Copy       copy            Copy datasets
*/

#include <cdi.h>

#include "process_int.h"
#include "cdo_vlist.h"
#include "progress.h"
#include "cdo_options.h"

bool
is_fdb_stream(const char *filename)
{
  auto pathlen = (int) strlen(filename);
  return (pathlen >= 4 && memcmp(filename, "fdb:", 4) == 0);
}

bool
is_fdb_copy(bool dataIsUnchanged, int nfiles)
{
  bool isFdbCopy = false;

  if (dataIsUnchanged)
    {
      isFdbCopy = is_fdb_stream(cdo_get_stream_name(nfiles));
      if (nfiles == 1 && !isFdbCopy) isFdbCopy = is_fdb_stream(cdo_get_stream_name(0));
    }

  return isFdbCopy;
}

void *
Copy(void *process)
{
  auto hasConstantFields = true;
  CdoStreamID streamID2 = CDO_STREAM_UNDEF;
  int vlistID2 = CDI_UNDEFID;
  int taxisID2 = CDI_UNDEFID;
  Field field;

  cdo_initialize(process);

  auto dataIsUnchanged = data_is_unchanged();

  // clang-format off
                 cdo_operator_add("copy",  0, 0, nullptr);
  auto CLONE   = cdo_operator_add("clone", 0, 0, nullptr);
  auto SZIP    = cdo_operator_add("szip",  0, 0, nullptr);
  // clang-format on

  auto operatorID = cdo_operator_id();
  if (operatorID == SZIP)
    {
      Options::cdoCompType = CDI_COMPRESS_SZIP;
      Options::cdoCompLevel = 0;
    }

  operator_check_argc(0);

  auto streamCnt = cdo_stream_cnt();
  auto nfiles = streamCnt - 1;

  auto isFdbCopy = is_fdb_copy(dataIsUnchanged, nfiles);

  progress::init();

  int tsID2 = 0;
  for (int indf = 0; indf < nfiles; ++indf)
    {
      if (Options::cdoVerbose) cdo_print("Process file: %s", cdo_get_stream_name(indf));

      auto streamID1 = cdo_open_read(indf);
      auto vlistID1 = cdo_stream_inq_vlist(streamID1);
      auto taxisID1 = vlistInqTaxis(vlistID1);

      VarList varList1;
      varListInit(varList1, vlistID1);

      if (indf == 0)
        {
          vlistID2 = vlistDuplicate(vlistID1);
          taxisID2 = taxisDuplicate(taxisID1);
          vlistDefTaxis(vlistID2, taxisID2);

          auto nvars = vlistNvars(vlistID1);

          auto ntsteps = vlistNtsteps(vlistID1);
          if (ntsteps == 1 && varList_numVaryingVars(varList1) == 0) ntsteps = 0;

          if (ntsteps == 0 && nfiles > 1)
            {
              hasConstantFields = false;
              for (int varID = 0; varID < nvars; ++varID) vlistDefVarTimetype(vlistID2, varID, TIME_VARYING);
            }
        }
      else { vlist_compare(vlistID1, vlistID2, CMP_ALL); }

      if (streamID2 == CDO_STREAM_UNDEF)
        {
          streamID2 = cdo_open_write(nfiles);
          cdo_def_vlist(streamID2, vlistID2);
        }

      auto ntsteps = vlistNtsteps(vlistID1);

      int tsID1 = 0;
      while (true)
        {
          auto nrecs = cdo_stream_inq_timestep(streamID1, tsID1);
          if (nrecs == 0) break;

          cdo_taxis_copy_timestep(taxisID2, taxisID1);
          cdo_def_timestep(streamID2, tsID2);

          for (int recID = 0; recID < nrecs; ++recID)
            {
              double fstatus = indf + ((ntsteps > 1) ? (tsID1 + (recID + 1.0) / nrecs) / ntsteps : 1.0);
              if (!Options::cdoVerbose) progress::update(0, 1, fstatus / nfiles);

              int varID, levelID;
              cdo_inq_record(streamID1, &varID, &levelID);

              if (hasConstantFields && tsID2 > 0 && tsID1 == 0)
                if (varList1[varID].isConstant) continue;

              cdo_def_record(streamID2, varID, levelID);

              if (dataIsUnchanged && (isFdbCopy || operatorID == CLONE || operatorID == SZIP))
                {
                  cdo_copy_record(streamID2, streamID1);
                }
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
    }

  cdo_stream_close(streamID2);

  if (vlistID2 != CDI_UNDEFID) vlistDestroy(vlistID2);

  cdo_finish();

  return nullptr;
}
