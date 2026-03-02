/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include "process_int.h"
#include "param_conversion.h"
#include <mpim_grid.h>
#include "matrix_view.h"

static void
shiftx(bool fillCyclic, int numberOfShifts, int nx, int ny, Varray<double> &v1, Varray<double> &v2, double missval)
{
  MatrixView<double> mv1(v1.data(), ny, nx);
  MatrixView<double> mv2(v2.data(), ny, nx);

  for (int i = 0; i < nx; ++i)
    {
      auto isCyclic = false;
      auto ins = i + numberOfShifts % nx;
      while (ins >= nx)
        {
          ins -= nx;
          isCyclic = true;
        }
      while (ins < 0)
        {
          ins += nx;
          isCyclic = true;
        }

      if (!fillCyclic && isCyclic)
        {
          for (int j = 0; j < ny; ++j) mv2[j][ins] = missval;
        }
      else
        {
          for (int j = 0; j < ny; ++j) mv2[j][ins] = mv1[j][i];
        }
    }
}

static void
shifty(bool fillCyclic, int numberOfShifts, int nx, int ny, Varray<double> &v1, Varray<double> &v2, double missval)
{
  MatrixView<double> mv1(v1.data(), ny, nx);
  MatrixView<double> mv2(v2.data(), ny, nx);

  for (int j = 0; j < ny; ++j)
    {
      auto isCyclic = false;
      auto jns = j + numberOfShifts % ny;

      while (jns >= ny)
        {
          jns -= ny;
          isCyclic = true;
        }
      while (jns < 0)
        {
          jns += ny;
          isCyclic = true;
        }

      if (!fillCyclic && isCyclic)
        {
          for (int i = 0; i < nx; ++i) mv2[jns][i] = missval;
        }
      else
        {
          for (int i = 0; i < nx; ++i) mv2[jns][i] = mv1[j][i];
        }
    }
}

static int
shiftx_coord(bool fillCyclic, int numberOfShifts, int gridID1)
{
  const auto gridID2 = gridDuplicate(gridID1);

  const auto nx = gridInqXsize(gridID1);
  auto ny = gridInqYsize(gridID1);
  if (gridInqType(gridID1) != GRID_CURVILINEAR) ny = 1;

  Varray<double> v1(nx * ny), v2(nx * ny);
  gridInqXvals(gridID1, v1.data());
  shiftx(fillCyclic, numberOfShifts, nx, ny, v1, v2, 0);
  gridDefXvals(gridID2, v2.data());

  if (gridInqXbounds(gridID1, nullptr))
    {
      const size_t nv = (gridInqType(gridID1) != GRID_CURVILINEAR) ? 2 : 4;

      Varray<double> bounds(nx * ny * nv);
      gridInqXbounds(gridID1, bounds.data());
      for (size_t k = 0; k < nv; ++k)
        {
          for (size_t i = 0; i < nx * ny; ++i) v1[i] = bounds[i * nv + k];
          shiftx(fillCyclic, numberOfShifts, nx, ny, v1, v2, 0);
          for (size_t i = 0; i < nx * ny; ++i) bounds[i * nv + k] = v2[i];
        }
      gridDefXbounds(gridID2, bounds.data());
    }

  return gridID2;
}

static int
shifty_coord(bool fillCyclic, int numberOfShifts, int gridID1)
{
  const auto gridID2 = gridDuplicate(gridID1);

  auto nx = gridInqXsize(gridID1);
  const auto ny = gridInqYsize(gridID1);
  if (gridInqType(gridID1) != GRID_CURVILINEAR) nx = 1;

  Varray<double> v1(nx * ny), v2(nx * ny);
  gridInqYvals(gridID1, v1.data());
  shifty(fillCyclic, numberOfShifts, nx, ny, v1, v2, 0);
  gridDefYvals(gridID2, v2.data());

  if (gridInqYbounds(gridID1, nullptr))
    {
      const size_t nv = (gridInqType(gridID1) != GRID_CURVILINEAR) ? 2 : 4;

      Varray<double> bounds(nx * ny * nv);
      gridInqYbounds(gridID1, bounds.data());
      for (size_t k = 0; k < nv; ++k)
        {
          for (size_t i = 0; i < nx * ny; ++i) v1[i] = bounds[i * nv + k];
          shifty(fillCyclic, numberOfShifts, nx, ny, v1, v2, 0);
          for (size_t i = 0; i < nx * ny; ++i) bounds[i * nv + k] = v2[i];
        }
      gridDefYbounds(gridID2, bounds.data());
    }

  return gridID2;
}

void *
Shiftxy(void *process)
{
  cdo_initialize(process);

  const auto SHIFTX = cdo_operator_add("shiftx", 0, 0, nullptr);
  const auto SHIFTY = cdo_operator_add("shifty", 0, 0, nullptr);

  const auto operatorID = cdo_operator_id();

  int numberOfShifts = 1;
  auto fillCyclic = false;
  auto shiftCoords = false;
  if (cdo_operator_argc() > 0)
    {
      numberOfShifts = parameter_to_int(cdo_operator_argv(0));
      auto pargc = cdo_operator_argc();
      auto pargv = cdo_get_oper_argv();
      for (int ic = 1; ic < pargc; ++ic)
        {
          if (pargv[ic] == "cyclic")
            fillCyclic = true;
          else if (pargv[ic] == "coord")
            shiftCoords = true;
        }
    }

  const auto streamID1 = cdo_open_read(0);

  const auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  const auto vlistID2 = vlistDuplicate(vlistID1);

  const auto taxisID1 = vlistInqTaxis(vlistID1);
  const auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  const auto nvars = vlistNvars(vlistID1);
  std::vector<bool> vars(nvars, false);

  const auto ngrids = vlistNgrids(vlistID1);
  for (int index = 0; index < ngrids; ++index)
    {
      const auto gridID1 = vlistGrid(vlistID1, index);
      const auto gridtype = gridInqType(gridID1);

      if (gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN || gridtype == GRID_CURVILINEAR
          || (gridtype == GRID_PROJECTION && gridInqProjType(gridID1) == CDI_PROJ_RLL)
          || (gridtype == GRID_GENERIC && gridInqXsize(gridID1) > 0 && gridInqYsize(gridID1) > 0))
        {
          if (shiftCoords)
            {
              int gridID2 = -1;
              if (operatorID == SHIFTX)
                gridID2 = shiftx_coord(fillCyclic, numberOfShifts, gridID1);
              else if (operatorID == SHIFTY)
                gridID2 = shifty_coord(fillCyclic, numberOfShifts, gridID1);

              vlistChangeGridIndex(vlistID2, index, gridID2);
            }

          for (int varID = 0; varID < nvars; ++varID)
            if (gridID1 == vlistInqVarGrid(vlistID1, varID)) vars[varID] = true;
        }
      else if (gridtype == GRID_GENERIC && gridInqXsize(gridID1) <= 1 && gridInqYsize(gridID1) <= 1) {}
      else { cdo_abort("Unsupported grid type: %s", gridNamePtr(gridtype)); }
    }

  {
    int varID;
    for (varID = 0; varID < nvars; ++varID)
      if (vars[varID]) break;
    if (varID >= nvars) cdo_warning("No variables selected!");
  }

  const auto streamID2 = cdo_open_write(1);
  cdo_def_vlist(streamID2, vlistID2);

  const auto gridsizemax = vlistGridsizeMax(vlistID1);
  Varray<double> array1(gridsizemax), array2(gridsizemax);

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

          cdo_def_record(streamID2, varID, levelID);

          if (vars[varID])
            {
              const auto gridID1 = vlistInqVarGrid(vlistID1, varID);
              const auto gridsize = gridInqSize(gridID1);
              const auto missval = vlistInqVarMissval(vlistID2, varID);

              const auto nx = gridInqXsize(gridID1);
              const auto ny = gridInqYsize(gridID1);

              if (operatorID == SHIFTX)
                shiftx(fillCyclic, numberOfShifts, nx, ny, array1, array2, missval);
              else if (operatorID == SHIFTY)
                shifty(fillCyclic, numberOfShifts, nx, ny, array1, array2, missval);

              nmiss = varray_num_mv(gridsize, array2, missval);
              cdo_write_record(streamID2, array2.data(), nmiss);
            }
          else { cdo_write_record(streamID2, array1.data(), nmiss); }
        }
      tsID++;
    }

  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  vlistDestroy(vlistID2);

  cdo_finish();

  return nullptr;
}
