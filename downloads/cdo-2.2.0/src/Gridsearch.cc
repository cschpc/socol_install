/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include "cdo_options.h"
#include "dmemory.h"
#include "process_int.h"
#include <mpim_grid.h>
#include "griddes.h"

struct grid_type
{
  int gridID;
  long size;
  long num_cell_corners;
  double *cell_corner_lon;
  double *cell_corner_lat;
};

struct cellsearch_type
{
  grid_type *srcGrid;
  grid_type *tgtGrid;
  float *src_cell_bound_box;
};

static grid_type *
grid_new(int gridID, const char *txt)
{
  bool lgrid_destroy = false;
  const auto gridtype = gridInqType(gridID);

  if (gridtype == GRID_GME)
    {
      lgrid_destroy = true;
      auto gridID_gme = gridToUnstructured(gridID, NeedCorners::Yes);
      gridCompress(gridID_gme);
      gridID = gridID_gme;
    }

  if (gridtype != GRID_UNSTRUCTURED && gridtype != GRID_CURVILINEAR)
    {
      lgrid_destroy = true;
      gridID = gridToCurvilinear(gridID, NeedCorners::Yes);
    }

  if (!gridHasCoordinates(gridID)) cdo_abort("%s grid corner missing!", txt);

  grid_type *grid = (grid_type *) Malloc(sizeof(grid_type));

  grid->gridID = gridID;
  grid->size = gridInqSize(grid->gridID);
  grid->num_cell_corners = (gridInqType(grid->gridID) == GRID_UNSTRUCTURED) ? gridInqNvertex(grid->gridID) : 4;

  // printf("%s grid size %ld nv %ld\n", txt, grid->size, grid->num_cell_corners);
  grid->cell_corner_lon = (double *) Malloc(grid->num_cell_corners * grid->size * sizeof(double));
  grid->cell_corner_lat = (double *) Malloc(grid->num_cell_corners * grid->size * sizeof(double));
  gridInqXbounds(grid->gridID, grid->cell_corner_lon);
  gridInqYbounds(grid->gridID, grid->cell_corner_lat);

  cdo_grid_to_radian(gridID, CDI_XAXIS, grid->num_cell_corners * grid->size, grid->cell_corner_lon, "grid corner lon");
  cdo_grid_to_radian(gridID, CDI_YAXIS, grid->num_cell_corners * grid->size, grid->cell_corner_lat, "grid corner lat");

  if (lgrid_destroy) gridDestroy(gridID);

  return grid;
}

static void
grid_delete(grid_type *grid)
{
  if (grid)
    {
      if (grid->cell_corner_lon) Free(grid->cell_corner_lon);
      if (grid->cell_corner_lat) Free(grid->cell_corner_lat);
      Free(grid);
    }
}

void
boundbox_from_corners1r(long ic, long nc, const double *corner_lon, const double *corner_lat, float *bound_box)
{
  const auto inc = ic * nc;

  auto clat = corner_lat[inc];
  auto clon = corner_lon[inc];

  bound_box[0] = clat;
  bound_box[1] = clat;
  bound_box[2] = clon;
  bound_box[3] = clon;

  for (long j = 1; j < nc; ++j)
    {
      clat = corner_lat[inc + j];
      clon = corner_lon[inc + j];

      if (clat < bound_box[0]) bound_box[0] = clat;
      if (clat > bound_box[1]) bound_box[1] = clat;
      if (clon < bound_box[2]) bound_box[2] = clon;
      if (clon > bound_box[3]) bound_box[3] = clon;
    }

  /*
  if ( std::fabs(bound_box[3] - bound_box[2]) > PI )
    {
      bound_box[2] = 0;
      bound_box[3] = PI2;
    }
  */
}

void
boundbox_from_corners(long size, long nc, const double *corner_lon, const double *corner_lat, float *bound_box)
{
  for (long i = 0; i < size; ++i)
    {
      const auto i4 = i << 2;  // *4
      const auto inc = i * nc;
      auto clat = corner_lat[inc];
      auto clon = corner_lon[inc];
      bound_box[i4] = clat;
      bound_box[i4 + 1] = clat;
      bound_box[i4 + 2] = clon;
      bound_box[i4 + 3] = clon;
      for (long j = 1; j < nc; ++j)
        {
          clat = corner_lat[inc + j];
          clon = corner_lon[inc + j];
          if (clat < bound_box[i4]) bound_box[i4] = clat;
          if (clat > bound_box[i4 + 1]) bound_box[i4 + 1] = clat;
          if (clon < bound_box[i4 + 2]) bound_box[i4 + 2] = clon;
          if (clon > bound_box[i4 + 3]) bound_box[i4 + 3] = clon;
        }
    }
}

static cellsearch_type *
cellsearch_new(grid_type *srcGrid, grid_type *tgtGrid)
{
  cellsearch_type *cellsearch = (cellsearch_type *) Malloc(sizeof(cellsearch_type));

  cellsearch->srcGrid = srcGrid;
  cellsearch->tgtGrid = tgtGrid;

  float *src_cell_bound_box = (float *) Malloc(4 * srcGrid->size * sizeof(double));

  boundbox_from_corners(srcGrid->size, srcGrid->num_cell_corners, srcGrid->cell_corner_lon, srcGrid->cell_corner_lat,
                        src_cell_bound_box);

  cellsearch->src_cell_bound_box = src_cell_bound_box;

  return cellsearch;
}

static void
cellsearch_delete(cellsearch_type *cellsearch)
{
  if (cellsearch) Free(cellsearch);
}

static long
search_cells(cellsearch_type *cellsearch, long tgtCellIndex, long *srch_add)
{
  grid_type *srcGrid = cellsearch->srcGrid;
  grid_type *tgtGrid = cellsearch->tgtGrid;
  float *src_cell_bound_box = cellsearch->src_cell_bound_box;

  float tgt_cell_bound_box[4];
  boundbox_from_corners1r(tgtCellIndex, tgtGrid->num_cell_corners, tgtGrid->cell_corner_lon, tgtGrid->cell_corner_lat,
                          tgt_cell_bound_box);

  const auto bound_box_lat1 = tgt_cell_bound_box[0];
  const auto bound_box_lat2 = tgt_cell_bound_box[1];
  const auto bound_box_lon1 = tgt_cell_bound_box[2];
  const auto bound_box_lon2 = tgt_cell_bound_box[3];

  long numSearchCells = 0;
  for (long srcCellIndex = 0; srcCellIndex < srcGrid->size; ++srcCellIndex)
    {
      const auto srcCellIndexM4 = srcCellIndex << 2;
      if ((src_cell_bound_box[srcCellIndexM4 + 2] <= bound_box_lon2) && (src_cell_bound_box[srcCellIndexM4 + 3] >= bound_box_lon1))
        {
          if ((src_cell_bound_box[srcCellIndexM4] <= bound_box_lat2) && (src_cell_bound_box[srcCellIndexM4 + 1] >= bound_box_lat1))
            {
              srch_add[numSearchCells] = srcCellIndex;
              numSearchCells++;
            }
        }
    }

  return numSearchCells;
}

static void
cell_search(int gridIDsrc, int gridIDtgt)
{
  grid_type *srcGrid = grid_new(gridIDsrc, "source");
  grid_type *tgtGrid = grid_new(gridIDtgt, "target");

  long *srch_add = (long *) Malloc(srcGrid->size * sizeof(long));

  cellsearch_type *cellsearch = cellsearch_new(srcGrid, tgtGrid);

  for (long tgtCellIndex = 0; tgtCellIndex < tgtGrid->size; ++tgtCellIndex)
    {
      long numSearchCells = search_cells(cellsearch, tgtCellIndex, srch_add);

      if (Options::cdoVerbose && numSearchCells > 0)
        {
          printf("tgt cell %ld: found %ld src cells\n", tgtCellIndex, numSearchCells);
          for (long n = 0; n < numSearchCells; ++n) printf("   %ld: %ld\n", n + 1, srch_add[n]);
        }
    }

  cellsearch_delete(cellsearch);
  grid_delete(srcGrid);
  grid_delete(tgtGrid);
  Free(srch_add);
}

void *
Gridsearch(void *process)
{
  cdo_initialize(process);

  // clang-format off
  const auto PSEARCH = cdo_operator_add("testpointsearch",  0,   0, nullptr);
  (void)PSEARCH;
  const auto CSEARCH = cdo_operator_add("testcellsearch",   0,   0, nullptr);
  // clang-format on

  const auto operatorID = cdo_operator_id();

  operator_input_arg("source and target grid description file or name");
  operator_check_argc(2);

  const auto gridID1 = cdo_define_grid(cdo_operator_argv(0));
  const auto gridID2 = cdo_define_grid(cdo_operator_argv(1));

  if (operatorID == CSEARCH) cell_search(gridID1, gridID2);

  cdo_finish();

  return nullptr;
}
