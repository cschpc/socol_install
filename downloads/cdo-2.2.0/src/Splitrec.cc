/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Split      splitrec        Split records
*/

#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "util_files.h"
#include "cdi_lockedIO.h"

void *
Splitrec(void *process)
{
  cdo_initialize(process);

  if (process_self().m_ID != 0) cdo_abort("This operator can't be combined with other operators!");

  operator_check_argc(0);

  auto dataIsUnchanged = data_is_unchanged();

  auto streamID1 = cdo_open_read(0);
  auto vlistID1 = cdo_stream_inq_vlist(streamID1);

  char filename[8192];
  strcpy(filename, cdo_get_obase().c_str());
  int nchars = strlen(filename);

  auto refname = cdo_get_stream_name(0);
  char filesuffix[32] = { 0 };
  FileUtils::gen_suffix(filesuffix, sizeof(filesuffix), cdo_inq_filetype(streamID1), vlistID1, refname);

  Field field;

  VarList varList1;
  varListInit(varList1, vlistID1);

  int index = 0;
  int tsID = 0;
  while (true)
    {
      auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
      if (nrecs == 0) break;

      for (int recID = 0; recID < nrecs; ++recID)
        {
          int varID, levelID;
          cdo_inq_record(streamID1, &varID, &levelID);

          vlistClearFlag(vlistID1);
          vlistDefFlag(vlistID1, varID, levelID, true);

          auto vlistID2 = vlistCreate();
          cdo_vlist_copy_flag(vlistID2, vlistID1);

          index++;
          sprintf(filename + nchars, "%06d", index);
          if (filesuffix[0]) sprintf(filename + nchars + 6, "%s", filesuffix);

          if (Options::cdoVerbose) cdo_print("create file %s", filename);

          auto streamID2 = cdo_open_write(filename);

          cdo_def_vlist(streamID2, vlistID2);

          auto varID2 = vlistFindVar(vlistID2, varID);
          auto levelID2 = vlistFindLevel(vlistID2, varID, levelID);

          cdo_def_timestep(streamID2, 0);
          cdo_def_record(streamID2, varID2, levelID2);
          if (dataIsUnchanged) { cdo_copy_record(streamID2, streamID1); }
          else
            {
              field.init(varList1[varID]);
              cdo_read_record(streamID1, field);
              cdo_write_record(streamID2, field);
            }

          cdo_stream_close(streamID2);
          vlistDestroy(vlistID2);
        }

      tsID++;
    }

  cdo_stream_close(streamID1);

  cdo_finish();

  return nullptr;
}
