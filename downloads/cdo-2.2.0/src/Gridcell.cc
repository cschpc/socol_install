/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Gridcell   gridarea        Grid cell area in m^2
      Gridcell   gridweights     Grid cell weights
      Gridcell   gridmask        Grid mask
*/

#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "cdo_zaxis.h"
#include <mpim_grid.h>
#include "constants.h"

static inline double
orthodrome(double px1, double py1, double px2, double py2)
{
  return std::acos(std::sin(py1) * std::sin(py2) + std::cos(py1) * std::cos(py2) * std::cos(px2 - px1));
}

double gridGetPlanetRadius(int gridID);

double
getPlanetRadius(int gridID)
{
  static auto radiusFromEnv = true;
  if (radiusFromEnv)
    {
      radiusFromEnv = false;
      auto envstr = getenv("PLANET_RADIUS");
      if (envstr)
        {
          auto fval = atof(envstr);
          if (fval > 0) PlanetRadius = fval;
        }
    }

  auto planetRadius = PlanetRadius;

  auto gridRadius = gridGetPlanetRadius(gridID);
  if (gridRadius > 1) planetRadius = gridRadius;

  return planetRadius;
}

void
gridcell_areas(int gridID, Varray<double> &array)
{
  auto gridtype = gridInqType(gridID);

  if (gridProjIsSupported(gridID) || gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN || gridtype == GRID_GME
      || gridtype == GRID_CURVILINEAR || gridtype == GRID_UNSTRUCTURED || gridtype == GRID_GAUSSIAN_REDUCED)
    {
      if (gridHasArea(gridID))
        {
          if (Options::cdoVerbose) cdo_print("Using existing grid cell area!");
          gridInqArea(gridID, array.data());
        }
      else
        {
          auto status = gridGenArea(gridID, array.data());
          if (status == 1)
            cdo_abort("%s: Cell corner coordinates missing!", __func__);
          else if (status == 2)
            cdo_abort("%s: Can't compute grid cell area for this grid!", __func__);

          auto planetRadius = getPlanetRadius(gridID);
          if (Options::cdoVerbose) cdo_print("Planet radius: %.2f", planetRadius);
          auto ngp = gridInqSize(gridID);
          for (size_t i = 0; i < ngp; ++i) array[i] *= planetRadius * planetRadius;
        }
    }
  else { cdo_abort("%s: Unsupported grid type: %s", __func__, gridNamePtr(gridtype)); }
}

static void
grid_dx(int gridID, Varray<double> &array)
{
  auto gridsize = gridInqSize(gridID);
  auto xsize = gridInqXsize(gridID);
  auto ysize = gridInqYsize(gridID);

  if (!gridHasCoordinates(gridID)) cdo_abort("Cell center coordinates missing!");

  Varray<double> xv(gridsize), yv(gridsize);
  gridInqXvals(gridID, xv.data());
  gridInqYvals(gridID, yv.data());

  // Convert lat/lon units if required
  cdo_grid_to_radian(gridID, CDI_XAXIS, gridsize, &xv[0], "grid longitudes");
  cdo_grid_to_radian(gridID, CDI_YAXIS, gridsize, &yv[0], "grid latitudes");

  auto planetRadius = getPlanetRadius(gridID);
  if (Options::cdoVerbose) cdo_print("Planet radius: %.2f", planetRadius);

  for (size_t j = 0; j < ysize; ++j)
    {
      auto joff = j * xsize;
      for (size_t i = 0; i < xsize; ++i)
        {
          double len1, len2;
          if (i == 0)
            {
              len2 = orthodrome(xv[joff + i], yv[joff + i], xv[joff + i + 1], yv[joff + i + 1]);
              len1 = len2;
            }
          else if (i == (xsize - 1))
            {
              len1 = orthodrome(xv[joff + i - 1], yv[joff + i - 1], xv[joff + i], yv[joff + i]);
              len2 = len1;
            }
          else
            {
              len1 = orthodrome(xv[joff + i - 1], yv[joff + i - 1], xv[joff + i], yv[joff + i]);
              len2 = orthodrome(xv[joff + i], yv[joff + i], xv[joff + i + 1], yv[joff + i + 1]);
            }

          array[joff + i] = 0.5 * (len1 + len2) * planetRadius;
        }
    }
}

static void
grid_dy(int gridID, Varray<double> &array)
{
  auto gridsize = gridInqSize(gridID);
  auto xsize = gridInqXsize(gridID);
  auto ysize = gridInqYsize(gridID);

  if (!gridHasCoordinates(gridID)) cdo_abort("Cell center coordinates missing!");

  Varray<double> xv(gridsize), yv(gridsize);
  gridInqXvals(gridID, xv.data());
  gridInqYvals(gridID, yv.data());

  // Convert lat/lon units if required
  cdo_grid_to_radian(gridID, CDI_XAXIS, gridsize, &xv[0], "grid longitudes");
  cdo_grid_to_radian(gridID, CDI_YAXIS, gridsize, &yv[0], "grid latitudes");

  auto planetRadius = getPlanetRadius(gridID);
  if (Options::cdoVerbose) cdo_print("Planet radius: %.2f", planetRadius);

  for (size_t i = 0; i < xsize; ++i)
    {
      for (size_t j = 0; j < ysize; ++j)
        {
          auto joff = j * xsize;
          auto joffp1 = (j + 1) * xsize;
          auto joffm1 = (j - 1) * xsize;
          double len1, len2;
          if (j == 0)
            {
              len2 = orthodrome(xv[joff + i], yv[joff + i], xv[joffp1 + i], yv[joffp1 + i]);
              len1 = len2;
            }
          else if (j == (ysize - 1))
            {
              len1 = orthodrome(xv[joffm1 + i], yv[joffm1 + i], xv[joff + i], yv[joff + i]);
              len2 = len1;
            }
          else
            {
              len1 = orthodrome(xv[joffm1 + i], yv[joffm1 + i], xv[joff + i], yv[joff + i]);
              len2 = orthodrome(xv[joff + i], yv[joff + i], xv[joffp1 + i], yv[joffp1 + i]);
            }

          array[joff + i] = 0.5 * (len1 + len2) * planetRadius;
        }
    }
}

void *
Gridcell(void *process)
{
  cdo_initialize(process);

  // clang-format off
  auto GRIDAREA    = cdo_operator_add("gridarea",     1,  0, nullptr);
  auto GRIDWGTS    = cdo_operator_add("gridweights",  1,  0, nullptr);
  auto GRIDMASK    = cdo_operator_add("gridmask",     0,  0, nullptr);
  auto GRIDDX      = cdo_operator_add("griddx",       1,  0, nullptr);
  auto GRIDDY      = cdo_operator_add("griddy",       1,  0, nullptr);
  auto GRIDCELLIDX = cdo_operator_add("gridcellidx",  0,  0, nullptr);
  // clang-format on

  auto operatorID = cdo_operator_id();

  operator_check_argc(0);

  // const bool needRadius = cdo_operator_f1(operatorID) > 0;

  auto streamID1 = cdo_open_read(0);
  auto vlistID1 = cdo_stream_inq_vlist(streamID1);

  auto ngrids = vlistNgrids(vlistID1);

  if (ngrids > 1) cdo_warning("Found more than 1 grid, using the first one!");

  auto gridID = vlistGrid(vlistID1, 0);
  auto zaxisID = zaxis_from_name("surface");

  auto vlistID2 = vlistCreate();
  auto varID = vlistDefVar(vlistID2, gridID, zaxisID, TIME_CONSTANT);
  vlistDefNtsteps(vlistID2, 0);

  if (operatorID == GRIDAREA)
    {
      cdiDefKeyString(vlistID2, varID, CDI_KEY_NAME, "cell_area");
      cdiDefKeyString(vlistID2, varID, CDI_KEY_STDNAME, "area");
      cdiDefKeyString(vlistID2, varID, CDI_KEY_LONGNAME, "area of grid cell");
      cdiDefKeyString(vlistID2, varID, CDI_KEY_UNITS, "m2");
      vlistDefVarDatatype(vlistID2, varID, CDI_DATATYPE_FLT64);
    }
  else if (operatorID == GRIDWGTS)
    {
      cdiDefKeyString(vlistID2, varID, CDI_KEY_NAME, "cell_weights");
      vlistDefVarDatatype(vlistID2, varID, CDI_DATATYPE_FLT64);
    }
  else if (operatorID == GRIDMASK)
    {
      cdiDefKeyString(vlistID2, varID, CDI_KEY_NAME, "grid_mask");
      vlistDefVarDatatype(vlistID2, varID, CDI_DATATYPE_UINT8);
    }
  else if (operatorID == GRIDDX)
    {
      cdiDefKeyString(vlistID2, varID, CDI_KEY_NAME, "dx");
      cdiDefKeyString(vlistID2, varID, CDI_KEY_LONGNAME, "delta x");
      cdiDefKeyString(vlistID2, varID, CDI_KEY_UNITS, "m");
    }
  else if (operatorID == GRIDDY)
    {
      cdiDefKeyString(vlistID2, varID, CDI_KEY_NAME, "dy");
      cdiDefKeyString(vlistID2, varID, CDI_KEY_LONGNAME, "delta y");
      cdiDefKeyString(vlistID2, varID, CDI_KEY_UNITS, "m");
    }
  else if (operatorID == GRIDCELLIDX)
    {
      cdiDefKeyString(vlistID2, varID, CDI_KEY_NAME, "gridcellidx");
      cdiDefKeyString(vlistID2, varID, CDI_KEY_LONGNAME, "grid cell index");
      vlistDefVarDatatype(vlistID2, varID, CDI_DATATYPE_INT32);
    }

  auto taxisID = cdo_taxis_create(TAXIS_ABSOLUTE);
  vlistDefTaxis(vlistID2, taxisID);

  auto gridsize = gridInqSize(gridID);
  Varray<double> array(gridsize);

  if (operatorID == GRIDAREA) { gridcell_areas(gridID, array); }
  else if (operatorID == GRIDWGTS)
    {
      auto status = gridcell_weights(gridID, array);
      if (status != 0) cdo_warning("Grid cell bounds not available, using constant grid cell area weights!");
    }
  else if (operatorID == GRIDMASK)
    {
      std::vector<int> mask(gridsize, 1);
      if (gridInqMask(gridID, nullptr)) gridInqMask(gridID, &mask[0]);

      for (size_t i = 0; i < gridsize; ++i) array[i] = mask[i];
    }
  else if (operatorID == GRIDCELLIDX)
    {
      for (size_t i = 0; i < gridsize; ++i) array[i] = i + 1;
    }
  else if (operatorID == GRIDDX || operatorID == GRIDDY)
    {
      auto gridtype = gridInqType(gridID);
      if (gridProjIsSupported(gridID) || gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN || gridtype == GRID_CURVILINEAR)
        {
          if (gridtype != GRID_CURVILINEAR) gridID = gridToCurvilinear(gridID, NeedCorners::Yes);

          (operatorID == GRIDDX) ? grid_dx(gridID, array) : grid_dy(gridID, array);
        }
      else { cdo_abort("Unsupported grid type: %s", gridNamePtr(gridtype)); }
    }

  auto streamID2 = cdo_open_write(1);
  cdo_def_vlist(streamID2, vlistID2);
  cdo_def_timestep(streamID2, 0);
  cdo_def_record(streamID2, 0, 0);
  cdo_write_record(streamID2, &array[0], 0);

  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  cdo_finish();

  return nullptr;
}
