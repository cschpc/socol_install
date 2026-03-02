/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

*/

#include <cdi.h>

#include "process_int.h"
#include <mpim_grid.h>
#include "field_functions.h"

void *
Fldrms(void *process)
{
  int lastgrid = -1;

  cdo_initialize(process);

  operator_check_argc(0);

  auto needWeights = true;

  auto streamID1 = cdo_open_read(0);
  auto streamID2 = cdo_open_read(1);

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto vlistID2 = cdo_stream_inq_vlist(streamID2);
  auto vlistID3 = vlistDuplicate(vlistID1);

  auto taxisID1 = vlistInqTaxis(vlistID1);
  auto taxisID3 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID3, taxisID3);

  double slon = 0.0, slat = 0.0;
  auto gridID3 = gridCreate(GRID_LONLAT, 1);
  gridDefXsize(gridID3, 1);
  gridDefYsize(gridID3, 1);
  gridDefXvals(gridID3, &slon);
  gridDefYvals(gridID3, &slat);

  auto ngrids = vlistNgrids(vlistID1);
  int ndiffgrids = 0;
  for (int index = 1; index < ngrids; ++index)
    if (vlistGrid(vlistID1, 0) != vlistGrid(vlistID1, index)) ndiffgrids++;

  auto gridID1 = vlistGrid(vlistID1, 0);
  auto gridID2 = vlistGrid(vlistID2, 0);

  if (gridInqSize(gridID1) != gridInqSize(gridID2)) cdo_abort("Fields have different grid size!");

  if (needWeights && gridInqType(gridID1) != GRID_LONLAT && gridInqType(gridID1) != GRID_GAUSSIAN)
    cdo_abort("Unsupported gridtype: %s", gridNamePtr(gridInqType(gridID1)));

  for (int index = 0; index < ngrids; ++index) vlistChangeGridIndex(vlistID3, index, gridID3);

  if (ndiffgrids > 0) cdo_abort("Too many different grids!");

  auto streamID3 = cdo_open_write(2);
  cdo_def_vlist(streamID3, vlistID3);

  auto gridsizemax = vlistGridsizeMax(vlistID1);

  VarList varList1, varList2;
  varListInit(varList1, vlistID1);
  varListInit(varList2, vlistID2);

  Field field1, field2, field3;
  field1.resize(gridsizemax);
  if (needWeights) field1.weightv.resize(gridsizemax);

  field2.resize(gridsizemax);

  field3.resize(1);
  field3.grid = gridID3;

  int tsID = 0;
  while (true)
    {
      auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
      if (nrecs == 0) break;

      auto nrecs2 = cdo_stream_inq_timestep(streamID2, tsID);
      if (nrecs2 == 0) cdo_abort("Input streams have different number of timesteps!");

      cdo_taxis_copy_timestep(taxisID3, taxisID1);

      cdo_def_timestep(streamID3, tsID);

      for (int recID = 0; recID < nrecs; ++recID)
        {
          int varID, levelID;
          cdo_inq_record(streamID1, &varID, &levelID);
          cdo_read_record(streamID1, field1.vec_d.data(), &field1.nmiss);
          cdo_inq_record(streamID2, &varID, &levelID);
          cdo_read_record(streamID2, field2.vec_d.data(), &field2.nmiss);

          field1.grid = varList1[varID].gridID;
          field2.grid = varList2[varID].gridID;

          if (needWeights && field1.grid != lastgrid)
            {
              lastgrid = field1.grid;
              field1.weightv[0] = 1;
              if (field1.size > 1)
                {
                  auto wstatus = gridcell_weights(field1.grid, field1.weightv);
                  if (wstatus != 0 && tsID == 0 && levelID == 0)
                    cdo_warning("Grid cell bounds not available, using constant grid cell area weights for variable %s!",
                                varList1[varID].name);
                }
            }

          field1.missval = varList1[varID].missval;
          field2.missval = varList1[varID].missval;
          field3.missval = varList1[varID].missval;

          field_rms(field1, field2, field3);

          cdo_def_record(streamID3, varID, levelID);
          cdo_write_record(streamID3, field3.vec_d.data(), field3.nmiss);
        }

      tsID++;
    }

  cdo_stream_close(streamID3);
  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  cdo_finish();

  return nullptr;
}
