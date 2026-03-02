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
#include "griddes.h"

void
genGridIndex(int gridID1, int gridID2, std::vector<long> &index)
{
  size_t i1;

  const auto gridtype1 = gridInqType(gridID1);
  const auto gridtype2 = gridInqType(gridID2);

  const auto gridsize2 = gridInqSize(gridID2);

  if (gridtype1 != gridtype2) cdo_abort("Input streams have different grid types!");

  for (size_t i = 0; i < gridsize2; ++i) index[i] = -1;

  if (gridtype1 == GRID_LONLAT || gridtype1 == GRID_GAUSSIAN)
    {
      const auto nlon1 = gridInqXsize(gridID1);
      const auto nlat1 = gridInqYsize(gridID1);

      const auto nlon2 = gridInqXsize(gridID2);
      const auto nlat2 = gridInqYsize(gridID2);

      if (!gridHasCoordinates(gridID1)) cdo_abort("Grid 1 has no values!");
      if (!gridHasCoordinates(gridID2)) cdo_abort("Grid 2 has no values!");

      Varray<double> xvals1(nlon1), yvals1(nlat1);
      Varray<double> xvals2(nlon2), yvals2(nlat2);
      std::vector<long> xindex(nlon2), yindex(nlat2);

      gridInqXvals(gridID1, xvals1.data());
      gridInqYvals(gridID1, yvals1.data());

      // Convert lat/lon units if required
      cdo_grid_to_degree(gridID1, CDI_XAXIS, nlon1, xvals1.data(), "grid1 center lon");
      cdo_grid_to_degree(gridID1, CDI_YAXIS, nlat1, yvals1.data(), "grid1 center lat");

      gridInqXvals(gridID2, xvals2.data());
      gridInqYvals(gridID2, yvals2.data());

      // Convert lat/lon units if required
      cdo_grid_to_degree(gridID2, CDI_XAXIS, nlon2, xvals2.data(), "grid2 center lon");
      cdo_grid_to_degree(gridID2, CDI_YAXIS, nlat2, yvals2.data(), "grid2 center lat");

      for (size_t i2 = 0; i2 < nlat2; ++i2)
        {
          for (i1 = 0; i1 < nlat1; ++i1)
            if (std::fabs(yvals2[i2] - yvals1[i1]) < 0.001) break;

          yindex[i2] = (i1 == nlat1) ? -1 : (long) i1;
        }

      for (size_t i2 = 0; i2 < nlon2; ++i2)
        {
          for (i1 = 0; i1 < nlon1; ++i1)
            if (std::fabs(xvals2[i2] - xvals1[i1]) < 0.001) break;

          if (i1 == nlon1)
            {
              if (xvals2[i2] < 0)
                {
                  for (i1 = 0; i1 < nlon1; ++i1)
                    if (std::fabs(xvals2[i2] + 360 - xvals1[i1]) < 0.001) break;
                }
              else if (xvals2[i2] > 180)
                {
                  for (i1 = 0; i1 < nlon1; ++i1)
                    if (std::fabs(xvals2[i2] - 360 - xvals1[i1]) < 0.001) break;
                }
            }

          xindex[i2] = (i1 == nlon1) ? -1 : (long) i1;
        }

      // for ( i2 = 0; i2 < nlon2; i2++ ) printf("x %d %d\n", i2, xindex[i2]);
      // for ( i2 = 0; i2 < nlat2; i2++ ) printf("y %d %d\n", i2, yindex[i2]);
      size_t k = 0;
      for (size_t j = 0; j < nlat2; ++j)
        for (size_t i = 0; i < nlon2; ++i)
          {
            if (xindex[i] == -1 || yindex[j] == -1)
              index[k++] = -1;
            else
              index[k++] = yindex[j] * nlon1 + xindex[i];
          }
    }
  else
    cdo_abort("Unsupported grid type: %s", gridNamePtr(gridtype1));
}

void *
Enlargegrid(void *process)
{
  cdo_initialize(process);

  operator_input_arg("grid description file or name");
  if (cdo_operator_argc() < 1) cdo_abort("Too few arguments!");
  if (cdo_operator_argc() > 2) cdo_abort("Too many arguments!");

  const auto gridID2 = cdo_define_grid(cdo_operator_argv(0));

  const auto streamID1 = cdo_open_read(0);

  const auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  const auto taxisID1 = vlistInqTaxis(vlistID1);
  const auto taxisID2 = taxisDuplicate(taxisID1);

  int ndiffgrids = 0;
  for (int index = 1; index < vlistNgrids(vlistID1); ++index)
    if (vlistGrid(vlistID1, 0) != vlistGrid(vlistID1, index)) ndiffgrids++;

  if (ndiffgrids > 0) cdo_abort("Too many different grids in %s!", cdo_get_stream_name(0));

  const auto gridID1 = vlistGrid(vlistID1, 0);

  const auto gridsize1 = gridInqSize(gridID1);
  const auto gridsize2 = gridInqSize(gridID2);

  Varray<double> array1(gridsize1), array2(gridsize2);
  std::vector<long> gindex(gridsize1);

  genGridIndex(gridID2, gridID1, gindex);

  const auto vlistID2 = vlistDuplicate(vlistID1);

  const auto ngrids = vlistNgrids(vlistID1);
  for (int index = 0; index < ngrids; ++index) vlistChangeGridIndex(vlistID2, index, gridID2);

  const auto streamID2 = cdo_open_write(1);

  vlistDefTaxis(vlistID2, taxisID2);
  cdo_def_vlist(streamID2, vlistID2);

  int tsID = 0;
  while (true)
    {
      const auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
      if (nrecs == 0) break;

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_def_timestep(streamID2, tsID);

      for (int recID = 0; recID < nrecs; ++recID)
        {
          int varID, levelID;
          cdo_inq_record(streamID1, &varID, &levelID);
          size_t nmiss;
          cdo_read_record(streamID1, array1.data(), &nmiss);

          const auto missval1 = vlistInqVarMissval(vlistID1, varID);

          for (size_t i = 0; i < gridsize2; ++i) array2[i] = missval1;
          for (size_t i = 0; i < gridsize1; ++i)
            if (gindex[i] >= 0) array2[gindex[i]] = array1[i];

          nmiss = varray_num_mv(gridsize2, array2, missval1);
          cdo_def_record(streamID2, varID, levelID);
          cdo_write_record(streamID2, array2.data(), nmiss);
        }

      tsID++;
    }

  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  cdo_finish();

  return nullptr;
}
