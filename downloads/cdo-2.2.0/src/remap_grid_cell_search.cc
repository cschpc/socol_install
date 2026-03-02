/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

//#define WITH_XYZCOORDS 1

#include "remap_grid_cell_search.h"
#include "cdo_output.h"
#include "compare.h"
#include "varray.h"
#include "cdo_options.h"
#include "cimdOmp.h"
#include <mpim_grid.h>

extern "C"
{
#include "lib/yac/sphere_part.h"
}

CellSearchMethod cellSearchMethod(CellSearchMethod::spherepart);

void
set_cell_search_method(const std::string &methodString)
{
  // clang-format off
  if      (methodString == "spherepart") cellSearchMethod = CellSearchMethod::spherepart;
  else if (methodString == "latbins")    cellSearchMethod = CellSearchMethod::latbins;
  else cdo_abort("Grid cell search method %s not available!", methodString);
  // clang-format on
}

static void
gridBoundboxReg2d(size_t nx, size_t ny, const Varray<double> &reg2d_corner_lon, const Varray<double> &reg2d_corner_lat,
                  double *grid_bound_box)
{
  grid_bound_box[0] = reg2d_corner_lat[0];
  grid_bound_box[1] = reg2d_corner_lat[ny];
  if (grid_bound_box[0] > grid_bound_box[1])
    {
      grid_bound_box[0] = reg2d_corner_lat[ny];
      grid_bound_box[1] = reg2d_corner_lat[0];
    }
  grid_bound_box[2] = reg2d_corner_lon[0];
  grid_bound_box[3] = reg2d_corner_lon[nx];
}

void
grid_cell_search_create_reg_2d(GridCellSearch &gcs, size_t dims[2], const Varray<double> &reg2d_corner_lon,
                               const Varray<double> &reg2d_corner_lat)
{
  gcs.is_reg2d = true;
  gcs.dims[0] = dims[0];
  gcs.dims[1] = dims[1];
  const auto nx = dims[0];
  const auto ny = dims[1];
  const auto nxp1 = nx + 1;
  const auto nyp1 = ny + 1;

  gcs.reg2d_corner_lon.resize(nxp1);
  gcs.reg2d_corner_lat.resize(nyp1);

  varray_copy(nxp1, reg2d_corner_lon, gcs.reg2d_corner_lon);
  varray_copy(nyp1, reg2d_corner_lat, gcs.reg2d_corner_lat);

  gridBoundboxReg2d(nx, ny, gcs.reg2d_corner_lon, gcs.reg2d_corner_lat, gcs.gridBoundboxReg2d);

  gcs.in_use = true;
}

void
grid_cell_search_create(GridCellSearch &gcs, size_t numCells, size_t numCellCorners, Varray<double> &cellCornerLon,
                        Varray<double> &cellCornerLat)
{
  gcs.method = cellSearchMethod;

  Varray<enum yac_edge_type> edgeTypes(numCellCorners, GREAT_CIRCLE_EDGE);
  Varray<struct grid_cell> cells(Threading::ompNumThreads);
  for (int i = 0; i < Threading::ompNumThreads; ++i)
    {
      cells[i].coordinates_xyz = new double[numCellCorners][3];
      cells[i].edge_type = edgeTypes.data();
      cells[i].num_corners = numCellCorners;
      cells[i].array_size = numCellCorners;
    }

  auto bndCircles = new struct bounding_circle[numCells];

#ifdef WITH_XYZCOORDS
  gcs.xyzCoords = new double[numCells * numCellCorners][3];
#endif

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
  for (size_t i = 0; i < numCells; ++i)
    {
      const auto ompthID = cdo_omp_get_thread_num();
      auto &cell = cells[ompthID];
      auto xyz = cell.coordinates_xyz;

      for (size_t k = 0; k < numCellCorners; ++k)
        gcLLtoXYZ(cellCornerLon[i * numCellCorners + k], cellCornerLat[i * numCellCorners + k], xyz[k]);

      if (numCellCorners == 3)
        yac_get_cell_bounding_circle_unstruct_triangle(xyz[0], xyz[1], xyz[2], &bndCircles[i]);
      else
        yac_get_cell_bounding_circle(cell, &bndCircles[i]);

#ifdef WITH_XYZCOORDS
      const auto offset = i * numCellCorners;
      for (size_t k = 0; k < numCellCorners; ++k)
        for (size_t l = 0; l < 3; ++l) gcs.xyzCoords[offset + k][l] = xyz[k][l];
#endif
    }

  gcs.yacSearch = yac_bnd_sphere_part_search_new(bndCircles, numCells);
  gcs.yacBndCircles = bndCircles;

  for (int i = 0; i < Threading::ompNumThreads; ++i) delete[] cells[i].coordinates_xyz;

  gcs.in_use = true;
}

void
grid_cell_search_delete(GridCellSearch &gcs)
{
  if (gcs.in_use)
    {
#ifdef WITH_XYZCOORDS
      if (gcs.xyzCoords) delete[] gcs.xyzCoords;
#endif
      varray_free(gcs.reg2d_corner_lon);
      varray_free(gcs.reg2d_corner_lat);

      if (gcs.yacSearch) yac_bnd_sphere_part_search_delete((struct bnd_sphere_part_search *) gcs.yacSearch);
      if (gcs.yacBndCircles) delete[](struct bounding_circle *) gcs.yacBndCircles;

      gcs.in_use = false;
    }
}

static size_t
doGridCellSearchYac(struct bnd_sphere_part_search *yacSearch, struct bounding_circle *yacBndCircles, bool isReg2dCell,
                    const grid_cell &yacGridCell, Varray<size_t> &srchAddr)
{
  size_t numCellCorners = yacGridCell.num_corners;
  struct bounding_circle bndCircle;
  auto xyz = yacGridCell.coordinates_xyz;

  if (numCellCorners == 4 && isReg2dCell)
    yac_get_cell_bounding_circle_reg_quad(xyz[0], xyz[1], xyz[2], &bndCircle);
  else if (numCellCorners == 3)
    yac_get_cell_bounding_circle_unstruct_triangle(xyz[0], xyz[1], xyz[2], &bndCircle);
  else
    yac_get_cell_bounding_circle(yacGridCell, &bndCircle);

  size_t numSearchCells;
  size_t *currNeighs;
  yac_bnd_sphere_part_search_do_bnd_circle_search(yacSearch, &bndCircle, 1, &currNeighs, &numSearchCells);

  if (srchAddr.size() < numSearchCells) srchAddr.resize(numSearchCells);

  size_t k = 0;
  // for (size_t i = 0; i < numSearchCells; ++i) srchAddr[i] = currNeighs[i];
  for (size_t i = 0; i < numSearchCells; ++i)
    {
      if (yac_extents_overlap(&bndCircle, &yacBndCircles[currNeighs[i]])) srchAddr[k++] = currNeighs[i];
    }
  numSearchCells = k;
  free(currNeighs);

  return numSearchCells;
}

size_t
do_grid_cell_search(GridCellSearch &gcs, bool isReg2dCell, GridCell &gridCell, Varray<size_t> &srchAddr)
{
  if (gcs.in_use)
    {
      // clang-format off
      if (gcs.method == CellSearchMethod::spherepart)
        return doGridCellSearchYac((struct bnd_sphere_part_search*)gcs.yacSearch, (struct bounding_circle *)gcs.yacBndCircles,
                                   isReg2dCell, gridCell.yacGridCell, srchAddr);
      else
        cdo_abort("%s::method undefined!", __func__);
      // clang-format on
    }

  return 0;
}
