/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include "process_int.h"

void *
Recttocomplex(void *process)
{
  cdo_initialize(process);

  operator_check_argc(0);

  auto streamID1 = cdo_open_read(0);
  auto vlistID1 = cdo_stream_inq_vlist(streamID1);

  auto streamID2 = cdo_open_read(1);
  auto vlistID2 = cdo_stream_inq_vlist(streamID2);

  auto vlistID3 = vlistDuplicate(vlistID1);

  auto taxisID1 = vlistInqTaxis(vlistID1);
  auto taxisID3 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID3, taxisID3);

  auto nvars = vlistNvars(vlistID3);
  for (int varID = 0; varID < nvars; ++varID)
    {
      auto datatype = vlistInqVarDatatype(vlistID2, varID);
      datatype = (datatype == CDI_DATATYPE_FLT64) ? CDI_DATATYPE_CPX64 : CDI_DATATYPE_CPX32;
      vlistDefVarDatatype(vlistID3, varID, datatype);
    }

  auto streamID3 = cdo_open_write(2);
  cdo_def_vlist(streamID3, vlistID3);

  auto gridsizemax = vlistGridsizeMax(vlistID1);
  Varray<double> array1(gridsizemax), array2(gridsizemax), array3(2 * gridsizemax);

  VarList varList1;
  varListInit(varList1, vlistID1);

  int tsID = 0;
  while (true)
    {
      auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
      if (nrecs == 0) break;
      auto nrecs2 = cdo_stream_inq_timestep(streamID2, tsID);
      if (nrecs2 == 0) break;

      cdo_taxis_copy_timestep(taxisID3, taxisID1);
      cdo_def_timestep(streamID3, tsID);

      for (int recID = 0; recID < nrecs; ++recID)
        {
          int varID, levelID;
          cdo_inq_record(streamID1, &varID, &levelID);
          cdo_def_record(streamID3, varID, levelID);

          cdo_inq_record(streamID2, &varID, &levelID);

          size_t nmiss;
          cdo_read_record(streamID1, array1.data(), &nmiss);
          cdo_read_record(streamID2, array2.data(), &nmiss);

          auto gridsize = varList1[varID].gridsize;
          for (size_t i = 0; i < gridsize; ++i)
            {
              array3[2 * i] = array1[i];
              array3[2 * i + 1] = array2[i];
            }

          cdo_write_record(streamID3, array3.data(), nmiss);
        }

      tsID++;
    }

  cdo_stream_close(streamID3);
  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  vlistDestroy(vlistID3);

  cdo_finish();

  return nullptr;
}
