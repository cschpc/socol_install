/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include "process_int.h"

extern "C" void cdiPrintKeys(int cdiID, int varID);

void *
Unpack(void *process)
{
  cdo_initialize(process);

  operator_check_argc(0);

  auto streamID1 = cdo_open_read(0);

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto vlistID2 = vlistDuplicate(vlistID1);

  auto nvars = vlistNvars(vlistID1);
  for (int varID = 0; varID < nvars; ++varID)
    {
      auto datatype = vlistInqVarDatatype(vlistID1, varID);
      double addoffset = 0.0, scalefactor = 1.0;
      auto haveAddoffset = (cdiInqKeyFloat(vlistID1, varID, CDI_KEY_ADDOFFSET, &addoffset) == CDI_NOERR);
      auto haveScalefactor = (cdiInqKeyFloat(vlistID1, varID, CDI_KEY_SCALEFACTOR, &scalefactor) == CDI_NOERR);
      if (haveAddoffset || haveScalefactor)
        {
          if (haveAddoffset) cdiDeleteKey(vlistID2, varID, CDI_KEY_ADDOFFSET);
          if (haveScalefactor) cdiDeleteKey(vlistID2, varID, CDI_KEY_SCALEFACTOR);
          if (datatype != CDI_DATATYPE_FLT64) vlistDefVarDatatype(vlistID2, varID, CDI_DATATYPE_FLT32);
        }
    }

  auto taxisID1 = vlistInqTaxis(vlistID1);
  auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  auto streamID2 = cdo_open_write(1);
  cdo_def_vlist(streamID2, vlistID2);

  Field field;

  VarList varList1;
  varListInit(varList1, vlistID1);

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
          field.init(varList1[varID]);
          cdo_read_record(streamID1, field);

          cdo_def_record(streamID2, varID, levelID);
          cdo_write_record(streamID2, field);
        }

      tsID++;
    }

  cdo_stream_close(streamID1);
  cdo_stream_close(streamID2);

  vlistDestroy(vlistID2);

  cdo_finish();

  return nullptr;
}
