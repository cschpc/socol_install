/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include "process_int.h"

void *
Tee(void *process)
{
  cdo_initialize(process);

  operator_check_argc(1);

  auto streamID1 = cdo_open_read(0);

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto taxisID1 = vlistInqTaxis(vlistID1);

  auto streamID2 = cdo_open_write(1);

  auto streamID3 = streamOpenWrite(cdo_operator_argv(0).c_str(), cdo_filetype());

  auto vlistID2 = vlistDuplicate(vlistID1);
  auto vlistID3 = vlistDuplicate(vlistID1);

  auto taxisID2 = taxisDuplicate(taxisID1);
  auto taxisID3 = taxisDuplicate(taxisID1);

  vlistDefTaxis(vlistID2, taxisID2);
  vlistDefTaxis(vlistID3, taxisID3);

  cdo_def_vlist(streamID2, vlistID2);
  streamDefVlist(streamID3, vlistID3);

  Field field;

  VarList varList1;
  varListInit(varList1, vlistID1);

  int tsID = 0;
  while (true)
    {
      auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
      if (nrecs == 0) break;

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_taxis_copy_timestep(taxisID3, taxisID1);

      cdo_def_timestep(streamID2, tsID);
      streamDefTimestep(streamID3, tsID);

      for (int recID = 0; recID < nrecs; ++recID)
        {
          int varID, levelID;
          cdo_inq_record(streamID1, &varID, &levelID);
          field.init(varList1[varID]);
          cdo_read_record(streamID1, field);

          cdo_def_record(streamID2, varID, levelID);
          cdo_write_record(streamID2, field);

          streamDefRecord(streamID3, varID, levelID);
          if (field.memType == MemType::Float)
            streamWriteRecordF(streamID3, field.vec_f.data(), field.nmiss);
          else
            streamWriteRecord(streamID3, field.vec_d.data(), field.nmiss);
        }

      tsID++;
    }

  cdo_stream_close(streamID1);
  cdo_stream_close(streamID2);
  streamClose(streamID3);

  cdo_finish();

  return nullptr;
}
