/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      writegrid Write grid
*/

#include <cdi.h>

#include "process_int.h"
#include <mpim_grid.h>
#include "griddes.h"

class ModuleWriteGrid
{
  CdoStreamID streamID;
  int gridID;
  std::vector<int> mask;


public:
  void
  init(void *process)
  {
    cdo_initialize(process);

    operator_check_argc(0);

    streamID = cdo_open_read(0);
    const auto vlistID = cdo_stream_inq_vlist(streamID);

    gridID = vlistGrid(vlistID, 0);

    const auto gridsize = gridInqSize(gridID);

    gridID = generate_full_cell_grid(gridID);

    if (!gridHasCoordinates(gridID)) cdo_abort("Cell corner coordinates missing!");

    mask =std::vector<int>(gridsize);

    if (gridInqMask(gridID, nullptr)) { gridInqMask(gridID, mask.data()); }
    else
      {
        for (size_t i = 0; i < gridsize; ++i) mask[i] = 1;
      }
  }
  void
  run()
  {
    write_nc_grid(cdo_get_stream_name(1), gridID, mask.data());
  }
  void
  close()
  {
    cdo_stream_close(streamID);
  }
};

void *
Writegrid(void *process)
{

  ModuleWriteGrid writegrid;
  writegrid.init(process);
  writegrid.run();
  writegrid.close();
  cdo_finish();

  return nullptr;
}
