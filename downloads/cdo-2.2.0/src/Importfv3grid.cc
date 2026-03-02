/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include "process_int.h"

void *
Importfv3grid(void *process)
{
  cdo_initialize(process);

  operator_check_argc(0);

  auto streamID1 = cdo_open_read(0);

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto nvars = vlistNvars(vlistID1);
  if (nvars != 5) cdo_abort("Found %d variables, expected 5 variables!", nvars);

  std::vector<std::string> vars(nvars);
  vars[0] = "grid_lon";
  vars[1] = "grid_lat";
  vars[2] = "grid_lont";
  vars[3] = "grid_latt";
  vars[4] = "area";

  for (int varID = 0; varID < nvars; ++varID)
    {
      auto varname = cdo::inq_var_name(vlistID1, varID);
      if (varname == vars[varID]) cdo_abort("Found variable %s, expected variable %s!", varname, vars[varID]);
    }

  auto ngrids = vlistNgrids(vlistID1);
  if (ngrids != 2) cdo_abort("Found %d grids, expected 2 grids!", nvars);

  auto gridIDi1 = vlistGrid(vlistID1, 0);
  auto gridIDi2 = vlistGrid(vlistID1, 1);

  auto nx = gridInqXsize(gridIDi1);
  auto ny = gridInqYsize(gridIDi1);
  auto gridsize1 = gridInqSize(gridIDi1);
  auto gridsize2 = gridInqSize(gridIDi2);

  cdo_stream_inq_timestep(streamID1, 0);

  auto gridIDo = gridCreate(GRID_UNSTRUCTURED, gridsize2);
  gridDefNvertex(gridIDo, 4);

  Varray<double> buffer(gridsize1);
  int varID, levelID;
  size_t nmiss;

  {
    Varray<double> grid_corner(4 * gridsize2);

    cdo_inq_record(streamID1, &varID, &levelID);  // grid_lon
    cdo_read_record(streamID1, buffer.data(), &nmiss);

    for (size_t j = 1; j < ny; ++j)
      for (size_t i = 1; i < nx; ++i)
        {
          grid_corner[((j - 1) * (nx - 1) + (i - 1)) * 4 + 0] = buffer[j * nx + i];
          grid_corner[((j - 1) * (nx - 1) + (i - 1)) * 4 + 1] = buffer[j * nx + (i - 1)];
          grid_corner[((j - 1) * (nx - 1) + (i - 1)) * 4 + 2] = buffer[(j - 1) * nx + (i - 1)];
          grid_corner[((j - 1) * (nx - 1) + (i - 1)) * 4 + 3] = buffer[(j - 1) * nx + i];
        }

    gridDefXbounds(gridIDo, grid_corner.data());

    cdo_inq_record(streamID1, &varID, &levelID);  // grid_lat
    cdo_read_record(streamID1, buffer.data(), &nmiss);

    for (size_t j = 1; j < ny; ++j)
      for (size_t i = 1; i < nx; ++i)
        {
          grid_corner[((j - 1) * (nx - 1) + (i - 1)) * 4 + 0] = buffer[j * nx + i];
          grid_corner[((j - 1) * (nx - 1) + (i - 1)) * 4 + 1] = buffer[j * nx + (i - 1)];
          grid_corner[((j - 1) * (nx - 1) + (i - 1)) * 4 + 2] = buffer[(j - 1) * nx + (i - 1)];
          grid_corner[((j - 1) * (nx - 1) + (i - 1)) * 4 + 3] = buffer[(j - 1) * nx + i];
        }

    gridDefYbounds(gridIDo, grid_corner.data());
  }

  cdo_inq_record(streamID1, &varID, &levelID);  // grid_lont
  cdo_read_record(streamID1, buffer.data(), &nmiss);
  gridDefXvals(gridIDo, buffer.data());

  cdo_inq_record(streamID1, &varID, &levelID);  // grid_latt
  cdo_read_record(streamID1, buffer.data(), &nmiss);
  gridDefYvals(gridIDo, buffer.data());

  cdo_inq_record(streamID1, &varID, &levelID);  // area
  cdo_read_record(streamID1, buffer.data(), &nmiss);

  double sfclevel = 0;
  auto surfaceID = zaxisCreate(ZAXIS_SURFACE, 1);
  zaxisDefLevels(surfaceID, &sfclevel);

  auto vlistID2 = vlistCreate();
  varID = vlistDefVar(vlistID2, gridIDo, surfaceID, TIME_CONSTANT);
  cdiDefKeyString(vlistID2, varID, CDI_KEY_NAME, "area");
  cdiDefKeyString(vlistID2, varID, CDI_KEY_STDNAME, "area");
  cdiDefKeyString(vlistID2, varID, CDI_KEY_LONGNAME, "cell area");
  cdiDefKeyString(vlistID2, varID, CDI_KEY_UNITS, "m2");
  vlistDefVarDatatype(vlistID2, varID, CDI_DATATYPE_FLT32);

  auto taxisID = cdo_taxis_create(TAXIS_ABSOLUTE);
  vlistDefTaxis(vlistID2, taxisID);

  auto streamID2 = cdo_open_write(1);
  cdo_def_vlist(streamID2, vlistID2);
  cdo_def_timestep(streamID2, 0);
  cdo_def_record(streamID2, 0, 0);
  cdo_write_record(streamID2, buffer.data(), 0);

  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  vlistDestroy(vlistID2);
  gridDestroy(gridIDo);
  zaxisDestroy(surfaceID);

  cdo_finish();

  return nullptr;
}
