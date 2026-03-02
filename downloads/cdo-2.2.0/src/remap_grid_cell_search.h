/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/
#ifndef GRID_CELL_SEARCH_H
#define GRID_CELL_SEARCH_H

#include <cstddef>
#include <string>
#include "varray.h"

extern "C"
{
#include "lib/yac/grid_cell.h"
}

struct GridCell
{
  double *coordinates_x;
  double *coordinates_y;
  grid_cell yacGridCell;
};

enum class CellSearchMethod
{
  spherepart,
  latbins
};

struct GridCellSearch
{
  bool in_use = false;

  CellSearchMethod method{ CellSearchMethod::spherepart };

  bool is_reg2d = false;
  size_t dims[2] = { 0 };

  // reg2d search
  double gridBoundboxReg2d[4] = { 0 };
  Varray<double> reg2d_corner_lon, reg2d_corner_lat;

  double (*xyzCoords)[3] = nullptr;
  void *yacBndCircles = nullptr;
  void *yacSearch = nullptr;
};

void grid_cell_search_create_reg_2d(GridCellSearch &gcs, size_t dims[2], const Varray<double> &reg2d_corner_lon,
                                    const Varray<double> &reg2d_corner_lat);
void grid_cell_search_create(GridCellSearch &gcs, size_t numCells, size_t numCellCorners, Varray<double> &cellCornerLon,
                             Varray<double> &cellCornerLat);
void grid_cell_search_delete(GridCellSearch &gcs);
size_t do_grid_cell_search(GridCellSearch &gcs, bool isLonLatCell, GridCell &gridCell, Varray<size_t> &srchAddr);
void set_cell_search_method(const std::string &methodString);

#endif
