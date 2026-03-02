/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <atomic>
#include <algorithm>  // std::sort

#include "process_int.h"
#include "cdo_wtime.h"
#include "remap.h"
#include "remap_store_link.h"
#include "cdo_options.h"
#include "progress.h"
#include "cimdOmp.h"
#include <mpim_grid.h>

extern "C"
{
#include "lib/yac/clipping.h"
#include "lib/yac/area.h"
#include "lib/yac/geometry.h"
}

struct CellSearch
{
  enum yac_edge_type *edgeType = nullptr;
  size_t numCellCorners = 0;
  size_t maxCells = 0;
  Varray<double> partialAreas;
  Varray<struct grid_cell> gridCells;
  Varray<struct grid_cell> overlapCells;
};

static void
cellsearch_realloc(size_t numCells, CellSearch &search)
{
  if (numCells > search.maxCells)
    {
      search.partialAreas.resize(numCells);
      search.overlapCells.resize(numCells);
      search.gridCells.resize(numCells);

      for (size_t i = search.maxCells; i < numCells; ++i)
        {
          search.overlapCells[i].array_size = 0;
          search.overlapCells[i].num_corners = 0;
          search.overlapCells[i].edge_type = nullptr;
          search.overlapCells[i].coordinates_xyz = nullptr;

          search.gridCells[i].array_size = search.numCellCorners;
          search.gridCells[i].num_corners = search.numCellCorners;
          search.gridCells[i].edge_type = search.edgeType;
          search.gridCells[i].coordinates_xyz = new double[search.numCellCorners][3];
        }

      search.maxCells = numCells;
    }
}

static void
cellsearch_free(CellSearch &search)
{
  for (size_t i = 0; i < search.maxCells; ++i)
    {
      if (search.overlapCells[i].array_size > 0)
        {
          if (search.overlapCells[i].coordinates_xyz) free(search.overlapCells[i].coordinates_xyz);
          if (search.overlapCells[i].edge_type) free(search.overlapCells[i].edge_type);
        }

      delete[] search.gridCells[i].coordinates_xyz;
    }

  varray_free(search.partialAreas);
  varray_free(search.overlapCells);
  varray_free(search.gridCells);
}

static void
boundbox_from_corners1r(size_t ic, size_t nc, const Varray<double> &cornerLon, const Varray<double> &cornerLat, float *boundBox)
{
  auto inc = ic * nc;

  float clat = cornerLat[inc];
  float clon = cornerLon[inc];

  boundBox[0] = clat;
  boundBox[1] = clat;
  boundBox[2] = clon;
  boundBox[3] = clon;

  for (size_t j = 1; j < nc; ++j)
    {
      clat = cornerLat[inc + j];
      clon = cornerLon[inc + j];

      if (clat < boundBox[0]) boundBox[0] = clat;
      if (clat > boundBox[1]) boundBox[1] = clat;
      if (clon < boundBox[2]) boundBox[2] = clon;
      if (clon > boundBox[3]) boundBox[3] = clon;
    }

  if (std::fabs(boundBox[3] - boundBox[2]) > PI_f)
    {
      boundBox[2] = 0;
      boundBox[3] = PI2_f;
    }
  /*
  if (std::fabs(boundBox[3] - boundBox[2]) > PI_f)
    {
      if (boundBox[3] > boundBox[2] && (boundBox[3]-PI2_f) < 0.0f)
        {
          float tmp = boundBox[2];
          boundBox[2] = boundBox[3] - PI2_f;
          boundBox[3] = tmp;
        }
    }
  */
}

static inline double
gridcell_area(const struct grid_cell &cell)
{
  return yac_huiliers_area(cell);
}

static void
cdo_compute_overlap_areas(size_t numSearchCells, CellSearch &search, const GridCell &gridCell)
{
  auto &overlapCells = search.overlapCells;

  // Do the clipping and get the cell for the overlapping area
  yac_cell_clipping(numSearchCells, search.gridCells.data(), gridCell.yacGridCell, overlapCells.data());

  // Get the partial areas for the overlapping regions
  for (size_t i = 0; i < numSearchCells; ++i) search.partialAreas[i] = gridcell_area(overlapCells[i]);

#ifdef VERBOSE
  for (size_t i = 0; i < numSearchCells; ++i) cdo_print("overlap area : %lf", search.partialAreas[i]);
#endif
}

static double
get_edge_direction(double *ref_corner, double *corner_a, double *corner_b)
{
  double edge_norm[3];
  crossproduct_ld(corner_a, corner_b, edge_norm);
  normalise_vector(edge_norm);

  // sine of the angle between the edge and the reference corner
  double angle = edge_norm[0] * ref_corner[0] + edge_norm[1] * ref_corner[1] + edge_norm[2] * ref_corner[2];

  // if the reference corner is directly on the edge
  // (for small angles sin(x)==x)
  if (fabs(angle) < yac_angle_tol) return 0.0;

  return copysign(1.0, angle);
}

static void
cdo_compute_concave_overlap_areas(size_t numSearchCells, CellSearch &search, const GridCell &gridCell)
{
  auto targetCell = gridCell.yacGridCell;
  auto &overlapAreas = search.partialAreas;
  auto &overlapCells = search.overlapCells;
  auto &sourceCells = search.gridCells;

  // common node point to all partial target cells
  double *baseCorner = targetCell.coordinates_xyz[0];

  // the triangulation algorithm only works for cells that only have great circle edges
  enum yac_edge_type edgeTypes[3] = { GREAT_CIRCLE_EDGE, GREAT_CIRCLE_EDGE, GREAT_CIRCLE_EDGE };
  double coordinates_xyz[3][3] = { { -1.0, -1.0, -1.0 }, { -1.0, -1.0, -1.0 }, { -1.0, -1.0, -1.0 } };
  coordinates_xyz[0][0] = baseCorner[0];
  coordinates_xyz[0][1] = baseCorner[1];
  coordinates_xyz[0][2] = baseCorner[2];

  // data structure to hold the triangles of the target cell
  struct grid_cell partialCell;
  partialCell.array_size = 3;
  partialCell.num_corners = 3;
  partialCell.coordinates_xyz = coordinates_xyz;
  partialCell.edge_type = edgeTypes;

  // Do the clipping and get the cell for the overlapping area

  for (size_t i = 0; i < numSearchCells; ++i) overlapAreas[i] = 0.0;

  // for all triangles of the target cell
  // (triangles a formed by first corner of the target cells and each edge of
  //  the cell; the first and last edge of the cell already, contain the
  //  first corner, therefore we can skip them)
  for (size_t cornerIdx = 1; cornerIdx < targetCell.num_corners - 1; ++cornerIdx)
    {
      auto cornerA = targetCell.coordinates_xyz[cornerIdx];
      auto cornerB = targetCell.coordinates_xyz[(cornerIdx + 1)];

      // if the current edge has a length of zero
      if (points_are_identically(cornerA, cornerB)) continue;

      auto edgeDirection = get_edge_direction(baseCorner, cornerA, cornerB);

      partialCell.coordinates_xyz[1][0] = cornerA[0];
      partialCell.coordinates_xyz[1][1] = cornerA[1];
      partialCell.coordinates_xyz[1][2] = cornerA[2];
      partialCell.coordinates_xyz[2][0] = cornerB[0];
      partialCell.coordinates_xyz[2][1] = cornerB[1];
      partialCell.coordinates_xyz[2][2] = cornerB[2];

      // clip the current target cell triangle with all source cells
      yac_cell_clipping(numSearchCells, sourceCells.data(), partialCell, overlapCells.data());

      // Get the partial areas for the overlapping regions as sum over the partial target cells.
      for (size_t i = 0; i < numSearchCells; ++i)
        {
          if (overlapCells[i].num_corners == 0) continue;

          overlapAreas[i] += gridcell_area(overlapCells[i]) * edgeDirection;
        }
    }

  for (size_t i = 0; i < numSearchCells; ++i)
    {
      if (overlapAreas[i] < 0.0) overlapAreas[i] = -overlapAreas[i];
    }

#ifdef VERBOSE
  for (size_t i = 0; i < numSearchCells; ++i) cdo_print("overlap area %zu: %lf", i, partialAreas[i]);
#endif
}

static void
set_cell_coordinates_yac(RemapGridType remapGridType, size_t cellIndex, size_t numCorners, const RemapGrid *remapGrid,
                         const struct grid_cell &yacGridCell, double *x, double *y)
{
  auto storeXY = (x && y);
  auto xyz = yacGridCell.coordinates_xyz;

  if (remapGridType == RemapGridType::Reg2D)
    {
      auto nx = remapGrid->dims[0];
      auto iy = cellIndex / nx;
      auto ix = cellIndex - iy * nx;
      auto reg2d_corner_lon = &remapGrid->reg2d_corner_lon[ix];
      auto reg2d_corner_lat = &remapGrid->reg2d_corner_lat[iy];
      constexpr int xi[4] = { 0, 1, 1, 0 };
      constexpr int yi[4] = { 0, 0, 1, 1 };
      for (int k = 0; k < 4; ++k)
        {
          auto lon = reg2d_corner_lon[xi[k]];
          auto lat = reg2d_corner_lat[yi[k]];
          gcLLtoXYZ(lon, lat, xyz[k]);
          if (storeXY)
            {
              x[k] = lon;
              y[k] = lat;
            }
        }
    }
  else
    {
      auto cell_corner_lon = &remapGrid->cell_corner_lon[cellIndex * numCorners];
      auto cell_corner_lat = &remapGrid->cell_corner_lat[cellIndex * numCorners];
      for (size_t i = 0; i < numCorners; ++i) gcLLtoXYZ(cell_corner_lon[i], cell_corner_lat[i], xyz[i]);

      if (storeXY)
        {
          for (size_t k = 0; k < numCorners; ++k)
            {
              x[k] = cell_corner_lon[k];
              y[k] = cell_corner_lat[k];
            }
        }
    }
}

static void
set_cell_coordinates(RemapGridType remapGridType, size_t cellIndex, size_t numCorners, const RemapGrid *remapGrid,
                     const GridCell &gridCell)
{
  set_cell_coordinates_yac(remapGridType, cellIndex, numCorners, remapGrid, gridCell.yacGridCell, gridCell.coordinates_x,
                           gridCell.coordinates_y);
}

static void
set_coordinates_yac(size_t numCells, RemapGridType remapGridType, const Varray<size_t> &cellIndices, size_t numCorners,
                    const RemapGrid *remapGrid, const Varray<grid_cell> &yacGridCell)
{
  for (size_t i = 0; i < numCells; ++i)
    set_cell_coordinates_yac(remapGridType, cellIndices[i], numCorners, remapGrid, yacGridCell[i], nullptr, nullptr);
}

static void
set_coordinates_yac(const double (*xyzCoords)[3], size_t numCells, const Varray<size_t> &cellIndices, size_t numCorners,
                    const Varray<grid_cell> &yacGridCell)
{
  for (size_t i = 0; i < numCells; ++i)
    {
      auto offset = cellIndices[i] * numCorners;
      auto xyz = yacGridCell[i].coordinates_xyz;
      for (size_t k = 0; k < numCorners; ++k)
        for (size_t l = 0; l < 3; ++l) xyz[k][l] = xyzCoords[offset + k][l];
    }
}

static void
gridcell_init_yac(GridCell &gridCell, size_t numCorners, enum yac_edge_type *edgeType)
{
  gridCell.yacGridCell.array_size = numCorners;
  gridCell.yacGridCell.num_corners = numCorners;
  gridCell.yacGridCell.edge_type = edgeType;
  gridCell.yacGridCell.coordinates_xyz = new double[numCorners][3];
  gridCell.coordinates_x = new double[numCorners];
  gridCell.coordinates_y = new double[numCorners];
}

static void
gridcell_free_yac(const GridCell &gridCell)
{
  delete[] gridCell.yacGridCell.coordinates_xyz;
  delete[] gridCell.coordinates_x;
  delete[] gridCell.coordinates_y;
}

static int
get_lonlat_circle_index(RemapGridType remapGridType, size_t gridsize, size_t numCorners, const Varray<double> &clon,
                        const Varray<double> &clat)
{
  int lonlatCircleIndex = -1;

  if (numCorners == 4)
    {
      if (remapGridType == RemapGridType::Reg2D) { lonlatCircleIndex = 1; }
      else
        {
          size_t iadd = (gridsize < 100) ? 1 : gridsize / 30 - 1;
          size_t num_i = 0, num_eq0 = 0, num_eq1 = 0;

          for (size_t i = 0; i < gridsize; i += iadd)
            {
              auto i4 = i * 4;
              num_i++;
              // clang-format off
              if (is_equal(clon[i4 + 1], clon[i4 + 2]) && is_equal(clon[i4 + 3], clon[i4 + 0]) &&
                  is_equal(clat[i4 + 0], clat[i4 + 1]) && is_equal(clat[i4 + 2], clat[i4 + 3]))
                {
                  num_eq1++;
                }
              else if (is_equal(clon[i4 + 0], clon[i4 + 1]) && is_equal(clon[i4 + 2], clon[i4 + 3]) &&
                       is_equal(clat[i4 + 1], clat[i4 + 2]) && is_equal(clat[i4 + 3], clat[i4 + 0]))
                {
                  num_eq0++;
                }
              // clang-format on
            }

          if (num_i == num_eq1) lonlatCircleIndex = 1;
          if (num_i == num_eq0) lonlatCircleIndex = 0;
        }
    }

  // printf("lonlatCircleIndex %d\n", lonlatCircleIndex);

  return lonlatCircleIndex;
}

static int
get_lonlat_circle_index(RemapGrid *remapGrid)
{
  auto gridsize = remapGrid->size;
  auto numCorners = remapGrid->num_cell_corners;
  const auto &clon = remapGrid->cell_corner_lon;
  const auto &clat = remapGrid->cell_corner_lat;
  return get_lonlat_circle_index(remapGrid->type, gridsize, numCorners, clon, clat);
}

static void
normalize_weights(RemapGrid *tgtGrid, RemapVars &rv)
{
  // Include centroids in weights and normalize using target cell area if requested
  auto numLinks = rv.numLinks;
  auto num_wts = rv.num_wts;

  if (rv.normOpt == NormOpt::DESTAREA)
    {
      const auto &cellArea = tgtGrid->cell_area;
#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
      for (size_t i = 0; i < numLinks; ++i)
        {
          auto index = rv.tgtCellIndices[i];  // current linear address for target grid cell
          auto normFactor = is_not_equal(cellArea[index], 0.0) ? 1.0 / cellArea[index] : 0.0;
          rv.wts[i * num_wts] *= normFactor;
        }
    }
  else if (rv.normOpt == NormOpt::FRACAREA)
    {
      const auto &cellFrac = tgtGrid->cell_frac;
#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
      for (size_t i = 0; i < numLinks; ++i)
        {
          auto index = rv.tgtCellIndices[i];  // current linear address for target grid cell
          auto normFactor = is_not_equal(cellFrac[index], 0.0) ? 1.0 / cellFrac[index] : 0.0;
          rv.wts[i * num_wts] *= normFactor;
        }
    }
  else if (rv.normOpt == NormOpt::NONE)
    {
    }
}

static void
normalize_weights(NormOpt normOpt, double cellArea, double cellFrac, size_t numWeights, Varray<double> &weights)
{
  if (normOpt == NormOpt::DESTAREA)
    {
      auto normFactor = is_not_equal(cellArea, 0.0) ? 1.0 / cellArea : 0.0;
      for (size_t i = 0; i < numWeights; ++i) weights[i] *= normFactor;
    }
  else if (normOpt == NormOpt::FRACAREA)
    {
      auto normFactor = is_not_equal(cellFrac, 0.0) ? 1.0 / cellFrac : 0.0;
      for (size_t i = 0; i < numWeights; ++i) weights[i] *= normFactor;
    }
  else if (normOpt == NormOpt::NONE)
    {
    }
}

static void
correct_weights(double cellArea, size_t numWeights, Varray<double> &weights)
{
  for (size_t i = 0; i < numWeights; ++i) weights[i] /= cellArea;
  yac_correct_weights(numWeights, weights.data());
  for (size_t i = 0; i < numWeights; ++i) weights[i] *= cellArea;
}

static void
sort_weights(size_t numWeights, Varray<size_t> &indices, Varray<double> &weights)
{
  if (is_sorted_list(numWeights, indices.data())) return;

  if (numWeights > 1)
    {
      struct IndexWeightX
      {
        size_t index;
        double weight;
      };

      Varray<IndexWeightX> indexWeights(numWeights);

      for (size_t i = 0; i < numWeights; ++i)
        {
          indexWeights[i].index = indices[i];
          indexWeights[i].weight = weights[i];
        }

      auto compareIndex = [](const auto &a, const auto &b) noexcept { return a.index < b.index; };
      std::sort(indexWeights.begin(), indexWeights.end(), compareIndex);

      for (size_t i = 0; i < numWeights; ++i)
        {
          indices[i] = indexWeights[i].index;
          weights[i] = indexWeights[i].weight;
        }
    }
}
/*
static void
reg2d_bound_box(RemapGrid *remapGrid, double *gridBoundBox)
{
  auto nx = remapGrid->dims[0];
  auto ny = remapGrid->dims[1];
  const auto &reg2d_corner_lon = remapGrid->reg2d_corner_lon;
  const auto &reg2d_corner_lat = remapGrid->reg2d_corner_lat;

  gridBoundBox[0] = reg2d_corner_lat[0];
  gridBoundBox[1] = reg2d_corner_lat[ny];
  if (gridBoundBox[0] > gridBoundBox[1])
    {
      gridBoundBox[0] = reg2d_corner_lat[ny];
      gridBoundBox[1] = reg2d_corner_lat[0];
    }
  gridBoundBox[2] = reg2d_corner_lon[0];
  gridBoundBox[3] = reg2d_corner_lon[nx];
}
*/

static void
scale_cellfrac(size_t numCells, Varray<double> &cellFrac, const Varray<double> &cellArea)
{
  for (size_t i = 0; i < numCells; ++i)
    if (is_not_equal(cellArea[i], 0)) cellFrac[i] /= cellArea[i];
}

static void
vec_add_weights(Varray<double> &vec, size_t numWeights, const Varray<double> &weights, const Varray<size_t> &indices)
{
  for (size_t i = 0; i < numWeights; ++i)
    {
      auto weight = weights[i];
      auto index = indices[i];
#ifndef __PGI
#ifdef _OPENMP
#pragma omp atomic
#endif
      vec[index] += weight;
#endif
    }
}

static size_t
remove_invalid_areas(size_t numSearchCells, Varray<double> &areas, Varray<size_t> &indices)
{
  size_t n = 0;
  for (size_t i = 0; i < numSearchCells; ++i)
    {
      if (areas[i] > 0.0)
        {
          areas[n] = areas[i];
          indices[n] = indices[i];
          n++;
        }
    }

  return n;
}

static size_t
remove_invalid_weights(size_t gridSize, size_t numWeights, Varray<double> &weights, Varray<size_t> &indices)
{
  size_t n = 0;
  for (size_t i = 0; i < numWeights; ++i)
    {
      auto index = (weights[i] > 0.0) ? indices[i] : gridSize;
      if (index != gridSize)
        {
          weights[n] = weights[i];
          indices[n] = index;
          n++;
        }
    }

  return n;
}

static size_t
remove_unmask_weights(const Varray<short> &gridMask, size_t numWeights, Varray<double> &weights, Varray<size_t> &indices)
{
  size_t n = 0;
  for (size_t i = 0; i < numWeights; ++i)
    {
      auto index = indices[i];
      /*
        Store the appropriate addresses and weights.
        Also add contributions to cell areas.
        The source grid mask is the master mask.
      */
      if (gridMask[index])
        {
          weights[n] = weights[i];
          indices[n] = index;
          n++;
        }
    }

  return n;
}

void
remap_conserv_weights(RemapSearch &remapSearch, RemapVars &rv)
{
  auto srcGrid = remapSearch.srcGrid;
  auto tgtGrid = remapSearch.tgtGrid;

  auto doCheck = true;

  auto srcGridType = srcGrid->type;
  auto tgtGridType = tgtGrid->type;

  if (Options::cdoVerbose) cdo_print("Called %s()", __func__);

  progress::init();

  auto start = Options::cdoVerbose ? cdo_get_wtime() : 0.0;

  auto srcGridSize = srcGrid->size;
  auto tgtGridSize = tgtGrid->size;

  auto srcNumCorners = srcGrid->num_cell_corners;
  auto tgtNumCorners = tgtGrid->num_cell_corners;

  enum yac_edge_type lonlatCircleType[] = { LON_CIRCLE_EDGE, LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE, LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE };
  Varray<enum yac_edge_type> greatCircleType(std::max(srcNumCorners, tgtNumCorners), GREAT_CIRCLE_EDGE);

  auto srcEdgeType = greatCircleType.data();
  auto tgtEdgeType = greatCircleType.data();

  enum yac_cell_type tgtCellType = MIXED_CELL;

  if (srcNumCorners == 4)
    {
      auto lonlatCircleIndex = get_lonlat_circle_index(srcGrid);
      if (lonlatCircleIndex >= 0) srcEdgeType = &lonlatCircleType[lonlatCircleIndex];
    }

  if (tgtNumCorners == 4)
    {
      auto lonlatCircleIndex = get_lonlat_circle_index(tgtGrid);
      if (lonlatCircleIndex >= 0)
        {
          tgtCellType = LON_LAT_CELL;
          tgtEdgeType = &lonlatCircleType[lonlatCircleIndex];
        }
    }

  Varray<GridCell> tgtGridCell2(Threading::ompNumThreads);
  for (int i = 0; i < Threading::ompNumThreads; ++i) gridcell_init_yac(tgtGridCell2[i], tgtNumCorners, tgtEdgeType);

  Varray<CellSearch> cellSearch2(Threading::ompNumThreads);
  for (int i = 0; i < Threading::ompNumThreads; ++i)
    {
      cellSearch2[i].numCellCorners = srcNumCorners;
      cellSearch2[i].edgeType = srcEdgeType;
    }

  auto numCellCorners = srcNumCorners;  // num of corners of search cells

  // double srcGridBoundBox[4];
  // if (srcGridType == RemapGridType::Reg2D) reg2d_bound_box(srcGrid, srcGridBoundBox);

  std::vector<WeightLinks> weightLinks(tgtGridSize);

  std::atomic<size_t> atomicCount{ 0 };

  size_t numSearchCellsStat[3] = { 0, 100000, 0 };

  extern CellSearchMethod cellSearchMethod;
  auto useCellsearch = (cellSearchMethod == CellSearchMethod::spherepart || srcGridType == RemapGridType::Reg2D);

  Varray<Varray<size_t>> indices2(Threading::ompNumThreads);
  if (!useCellsearch)
    for (int i = 0; i < Threading::ompNumThreads; ++i) indices2[i].resize(srcGridSize);

  // Loop over target grid

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) default(shared)
#endif
  for (size_t tgtCellIndex = 0; tgtCellIndex < tgtGridSize; ++tgtCellIndex)
    {
      auto ompthID = cdo_omp_get_thread_num();
      auto &cellSearch = cellSearch2[ompthID];
      auto &indices = indices2[ompthID];
      auto &tgtGridCell = tgtGridCell2[ompthID];

      atomicCount++;
      if (ompthID == 0) progress::update(0, 1, (double) atomicCount / tgtGridSize);

      weightLinks[tgtCellIndex].nlinks = 0;

      if (!tgtGrid->mask[tgtCellIndex]) continue;

      set_cell_coordinates(tgtGridType, tgtCellIndex, tgtNumCorners, tgtGrid, tgtGridCell);

      // Get search cells
      size_t numSearchCells;
      if (useCellsearch)
        {
          // numSearchCells = remap_search_cells(remapSearch, tgtCellType == LON_LAT_CELL, tgtGridCell, indices);
          numSearchCells = remap_search_cells(remapSearch, tgtGridType == RemapGridType::Reg2D, tgtGridCell, indices);
        }
      else
        {
          float tgtCellBoundBoxR[4];
          boundbox_from_corners1r(tgtCellIndex, tgtNumCorners, tgtGrid->cell_corner_lon, tgtGrid->cell_corner_lat,
                                  tgtCellBoundBoxR);

          numSearchCells = get_srch_cells(tgtCellIndex, remapSearch.tgtBins, remapSearch.srcBins, tgtCellBoundBoxR, indices);
        }

      if (1 && Options::cdoVerbose)
        {
          numSearchCellsStat[0] += numSearchCells;
          numSearchCellsStat[1] = std::min(numSearchCellsStat[1], numSearchCells);
          numSearchCellsStat[2] = std::max(numSearchCellsStat[2], numSearchCells);
        }

      if (0 && Options::cdoVerbose) cdo_print("tgtCellIndex %zu  numSearchCells %zu", tgtCellIndex, numSearchCells);

      if (numSearchCells == 0) continue;

      // Create search arrays

      cellsearch_realloc(numSearchCells, cellSearch);

      if (remapSearch.gcs.xyzCoords)
        set_coordinates_yac(remapSearch.gcs.xyzCoords, numSearchCells, indices, numCellCorners, cellSearch.gridCells);
      else
        set_coordinates_yac(numSearchCells, srcGridType, indices, numCellCorners, srcGrid, cellSearch.gridCells);

      if (tgtNumCorners < 4 || tgtCellType == LON_LAT_CELL)
        cdo_compute_overlap_areas(numSearchCells, cellSearch, tgtGridCell);
      else
        cdo_compute_concave_overlap_areas(numSearchCells, cellSearch, tgtGridCell);

      auto &partialWeights = cellSearch.partialAreas;

      auto numWeights = remove_invalid_areas(numSearchCells, partialWeights, indices);

      auto tgtCellArea = gridcell_area(tgtGridCell.yacGridCell);
      tgtGrid->cell_area[tgtCellIndex] = tgtCellArea;

      if (rv.normOpt == NormOpt::FRACAREA) correct_weights(tgtCellArea, numWeights, partialWeights);

      numWeights = remove_invalid_weights(srcGridSize, numWeights, partialWeights, indices);

      vec_add_weights(srcGrid->cell_area, numWeights, partialWeights, indices);

      numWeights = remove_unmask_weights(srcGrid->mask, numWeights, partialWeights, indices);

      vec_add_weights(srcGrid->cell_frac, numWeights, partialWeights, indices);

      tgtGrid->cell_frac[tgtCellIndex] = varray_sum(numWeights, partialWeights);

      store_weightlinks(1, numWeights, indices.data(), partialWeights.data(), tgtCellIndex, weightLinks);
    }

  progress::update(0, 1, 1);

  if (1 && Options::cdoVerbose)
    {
      cdo_print("Num search cells min,mean,max :  %zu  %3.1f  %zu", numSearchCellsStat[1],
                numSearchCellsStat[0] / (double) tgtGridSize, numSearchCellsStat[2]);
    }

  // Finished with all cells: deallocate search arrays
  for (auto ompthID = 0; ompthID < Threading::ompNumThreads; ++ompthID)
    {
      cellsearch_free(cellSearch2[ompthID]);
      gridcell_free_yac(tgtGridCell2[ompthID]);
    }

  weight_links_to_remap_links(1, tgtGridSize, weightLinks, rv);

  // Normalize weights using target cell area if requested
  normalize_weights(tgtGrid, rv);

  if (Options::cdoVerbose) cdo_print("Total number of links = %zu", rv.numLinks);

  scale_cellfrac(srcGridSize, srcGrid->cell_frac, srcGrid->cell_area);
  scale_cellfrac(tgtGridSize, tgtGrid->cell_frac, tgtGrid->cell_area);

  // Perform some error checking on final weights
  if (doCheck)
    {
      remap_check_area(srcGridSize, srcGrid->cell_area, "Source");
      remap_check_area(tgtGridSize, tgtGrid->cell_area, "Target");

      remap_vars_check_weights(rv);
    }

  if (Options::cdoVerbose) cdo_print("Cells search: %.2f seconds", cdo_get_wtime() - start);
}  // remap_conserv_weights

template <typename T>
static double
conserv_remap(const Varray<T> &srcArray, size_t numWeights, const Varray<double> &weights, const Varray<size_t> &srcIndices)
{
  double tgtPoint = 0.0;
  for (size_t i = 0; i < numWeights; ++i) tgtPoint += srcArray[srcIndices[i]] * weights[i];

  return tgtPoint;
}

template <typename T>
static void
remap_conserv(NormOpt normOpt, RemapSearch &remapSearch, const Varray<T> &srcArray, Varray<T> &tgtArray, T missval, size_t nmiss)
{
  auto srcGrid = remapSearch.srcGrid;
  auto tgtGrid = remapSearch.tgtGrid;

  auto doCheck = true;

  // Variables necessary if segment manages to hit pole
  auto srcGridType = srcGrid->type;
  auto tgtGridType = tgtGrid->type;

  if (Options::cdoVerbose) cdo_print("Called %s()", __func__);

  progress::init();

  auto start = Options::cdoVerbose ? cdo_get_wtime() : 0.0;

  auto srcGridSize = srcGrid->size;
  auto tgtGridSize = tgtGrid->size;

  Varray<short> srcGridMask;
  if (nmiss)
    {
      srcGridMask.resize(srcGridSize, 1);
      remap_set_mask(srcGridSize, srcArray, missval, srcGridMask);
    }

  auto srcNumCorners = srcGrid->num_cell_corners;
  auto tgtNumCorners = tgtGrid->num_cell_corners;

  enum yac_edge_type lonlatCircleType[] = { LON_CIRCLE_EDGE, LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE, LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE };
  Varray<enum yac_edge_type> greatCircleType(std::max(srcNumCorners, tgtNumCorners), GREAT_CIRCLE_EDGE);

  auto srcEdgeType = greatCircleType.data();
  auto tgtEdgeType = greatCircleType.data();

  enum yac_cell_type tgtCellType = MIXED_CELL;

  if (srcNumCorners == 4)
    {
      auto lonlatCircleIndex = get_lonlat_circle_index(srcGrid);
      if (lonlatCircleIndex >= 0) srcEdgeType = &lonlatCircleType[lonlatCircleIndex];
    }

  if (tgtNumCorners == 4)
    {
      auto lonlatCircleIndex = get_lonlat_circle_index(tgtGrid);
      if (lonlatCircleIndex >= 0)
        {
          tgtCellType = LON_LAT_CELL;
          tgtEdgeType = &lonlatCircleType[lonlatCircleIndex];
        }
    }

  Varray<GridCell> tgtGridCell2(Threading::ompNumThreads);
  for (int i = 0; i < Threading::ompNumThreads; ++i) gridcell_init_yac(tgtGridCell2[i], tgtNumCorners, tgtEdgeType);

  Varray<CellSearch> cellSearch2(Threading::ompNumThreads);
  for (int i = 0; i < Threading::ompNumThreads; ++i)
    {
      cellSearch2[i].numCellCorners = srcNumCorners;
      cellSearch2[i].edgeType = srcEdgeType;
    }

  auto numCellCorners = srcNumCorners;  // num of corners of search cells

  // double srcGridBoundBox[4];
  // if (srcGridType == RemapGridType::Reg2D) reg2d_bound_box(srcGrid, srcGridBoundBox);

  std::atomic<size_t> atomicCount{ 0 };

  size_t numSearchCellsStat[3] = { 0, 100000, 0 };

  extern CellSearchMethod cellSearchMethod;
  auto useCellsearch = (cellSearchMethod == CellSearchMethod::spherepart || srcGridType == RemapGridType::Reg2D);

  Varray<Varray<size_t>> indices2(Threading::ompNumThreads);
  if (!useCellsearch)
    for (int i = 0; i < Threading::ompNumThreads; ++i) indices2[i].resize(srcGridSize);

      // Loop over target grid

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) default(shared)
#endif
  for (size_t tgtCellIndex = 0; tgtCellIndex < tgtGridSize; ++tgtCellIndex)
    {
      auto ompthID = cdo_omp_get_thread_num();
      auto &cellSearch = cellSearch2[ompthID];
      auto &indices = indices2[ompthID];
      auto &tgtGridCell = tgtGridCell2[ompthID];

      atomicCount++;
      if (ompthID == 0) progress::update(0, 1, (double) atomicCount / tgtGridSize);

      tgtArray[tgtCellIndex] = missval;

      if (!tgtGrid->mask[tgtCellIndex]) continue;

      set_cell_coordinates(tgtGridType, tgtCellIndex, tgtNumCorners, tgtGrid, tgtGridCell);

      // Get search cells
      size_t numSearchCells;
      if (useCellsearch)
        {
          // numSearchCells = remap_search_cells(remapSearch, tgtCellType == LON_LAT_CELL, tgtGridCell, indices);
          numSearchCells = remap_search_cells(remapSearch, tgtGridType == RemapGridType::Reg2D, tgtGridCell, indices);
        }
      else
        {
          float tgtCellBoundBoxR[4];
          boundbox_from_corners1r(tgtCellIndex, tgtNumCorners, tgtGrid->cell_corner_lon, tgtGrid->cell_corner_lat,
                                  tgtCellBoundBoxR);

          numSearchCells = get_srch_cells(tgtCellIndex, remapSearch.tgtBins, remapSearch.srcBins, tgtCellBoundBoxR, indices);
        }

      if (1 && Options::cdoVerbose)
        {
          numSearchCellsStat[0] += numSearchCells;
          numSearchCellsStat[1] = std::min(numSearchCellsStat[1], numSearchCells);
          numSearchCellsStat[2] = std::max(numSearchCellsStat[2], numSearchCells);
        }

      if (0 && Options::cdoVerbose) cdo_print("tgtCellIndex %zu  numSearchCells %zu", tgtCellIndex, numSearchCells);

      if (numSearchCells == 0) continue;

      // Create search arrays

      cellsearch_realloc(numSearchCells, cellSearch);

      if (remapSearch.gcs.xyzCoords)
        set_coordinates_yac(remapSearch.gcs.xyzCoords, numSearchCells, indices, numCellCorners, cellSearch.gridCells);
      else
        set_coordinates_yac(numSearchCells, srcGridType, indices, numCellCorners, srcGrid, cellSearch.gridCells);

      if (tgtNumCorners < 4 || tgtCellType == LON_LAT_CELL)
        cdo_compute_overlap_areas(numSearchCells, cellSearch, tgtGridCell);
      else
        cdo_compute_concave_overlap_areas(numSearchCells, cellSearch, tgtGridCell);

      auto &partialWeights = cellSearch.partialAreas;

      auto numWeights = remove_invalid_areas(numSearchCells, partialWeights, indices);

      auto tgtCellArea = gridcell_area(tgtGridCell.yacGridCell);
      tgtGrid->cell_area[tgtCellIndex] = tgtCellArea;

      if (normOpt == NormOpt::FRACAREA) correct_weights(tgtCellArea, numWeights, partialWeights);

      numWeights = remove_invalid_weights(srcGridSize, numWeights, partialWeights, indices);
      if (srcGridMask.size() > 0) numWeights = remove_unmask_weights(srcGridMask, numWeights, partialWeights, indices);

      tgtGrid->cell_frac[tgtCellIndex] = varray_sum(numWeights, partialWeights);

      if (numWeights)
        {
          sort_weights(numWeights, indices, partialWeights);
          // Normalize weights using cell target area if requested
          normalize_weights(normOpt, tgtCellArea, tgtGrid->cell_frac[tgtCellIndex], numWeights, partialWeights);
          tgtArray[tgtCellIndex] = conserv_remap(srcArray, numWeights, partialWeights, indices);
        }
    }

  progress::update(0, 1, 1);

  if (1 && Options::cdoVerbose)
    {
      cdo_print("Num search cells min,mean,max :  %zu  %3.1f  %zu", numSearchCellsStat[1],
                numSearchCellsStat[0] / (double) tgtGridSize, numSearchCellsStat[2]);
    }

  // Finished with all cells: deallocate search arrays

  for (auto ompthID = 0; ompthID < Threading::ompNumThreads; ++ompthID)
    {
      cellsearch_free(cellSearch2[ompthID]);
      gridcell_free_yac(tgtGridCell2[ompthID]);
    }

  scale_cellfrac(tgtGridSize, tgtGrid->cell_frac, tgtGrid->cell_area);

  // Perform some error checking on final weights
  if (doCheck) remap_check_area(tgtGridSize, tgtGrid->cell_area, "Target");

  if (Options::cdoVerbose) cdo_print("Cells search: %.2f seconds", cdo_get_wtime() - start);
}  // remap_conserv

void
remap_conserv(NormOpt normOpt, RemapSearch &remapSearch, const Field &field1, Field &field2)
{
  if (field1.memType != field2.memType) cdo_abort("Interal error, memType of field1 and field2 differ!");

  if (field1.memType == MemType::Float)
    remap_conserv(normOpt, remapSearch, field1.vec_f, field2.vec_f, (float) field1.missval, field1.nmiss);
  else
    remap_conserv(normOpt, remapSearch, field1.vec_d, field2.vec_d, field1.missval, field1.nmiss);
}

template <typename T>
static size_t
remove_missing_weights(const Varray<T> &srcArray, T missval, size_t numWeights, Varray<double> &partialWeights,
                       Varray<size_t> &indices)
{
  size_t n = 0;
  for (size_t i = 0; i < numWeights; ++i)
    {
      auto cellIndex = indices[i];
      if (!dbl_is_equal(srcArray[cellIndex], missval))
        {
          partialWeights[n] = partialWeights[i];
          indices[n] = cellIndex;
          n++;
        }
    }

  return n;
}

static double
sphere_segment_area(double latInRadian)
{
  return 2.0 * M_PI * (1.0 - std::cos(M_PI * 0.5 - latInRadian));
}

static double
latitude_area(double latMin, double latMax)
{
  return sphere_segment_area(latMin) - sphere_segment_area(latMax);
}

void
remap_weights_zonal_mean(int gridID1, int gridID2, Varray2D<size_t> &remapIndices, Varray2D<double> &remapWeights)
{
  constexpr double scaleFactor = 1000000000.0;
  auto gridsize1 = gridInqSize(gridID1);
  size_t nv1 = gridInqNvertex(gridID1);
  Varray<double> xbounds1(gridsize1 * nv1), ybounds1(gridsize1 * nv1);
  gridInqXbounds(gridID1, xbounds1.data());
  gridInqYbounds(gridID1, ybounds1.data());

  // Convert lonlat units if required
  cdo_grid_to_radian(gridID1, CDI_XAXIS, gridsize1 * nv1, xbounds1.data(), "source grid longitude bounds");
  cdo_grid_to_radian(gridID1, CDI_YAXIS, gridsize1 * nv1, ybounds1.data(), "source grid latitude bounds");

  Varray<int> ymin1(gridsize1), ymax1(gridsize1);

#ifdef _OPENMP
#pragma omp parallel for schedule(static) default(shared)
#endif
  for (size_t i = 0; i < gridsize1; ++i)
    {
      auto minval = ybounds1[i * nv1];
      auto maxval = ybounds1[i * nv1];
      for (size_t k = 1; k < nv1; ++k)
        {
          auto val = ybounds1[i * nv1 + k];
          maxval = std::max(maxval, val);
          minval = std::min(minval, val);
        }
      ymin1[i] = scaleFactor * minval;
      ymax1[i] = scaleFactor * maxval;
    }

  auto ysize2 = gridInqYsize(gridID2);
  Varray<double> ybounds2(ysize2 * 2);
  gridInqYbounds(gridID2, ybounds2.data());

  // Convert lat units if required
  cdo_grid_to_radian(gridID2, CDI_YAXIS, ysize2 * 2, ybounds2.data(), "target grid latitude bounds");

  remapIndices.resize(ysize2);
  remapWeights.resize(ysize2);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) default(shared)
#endif
  for (size_t i2 = 0; i2 < ysize2; ++i2)
    {
      int latBounds[2] = { static_cast<int>(scaleFactor * ybounds2[2 * i2]), static_cast<int>(scaleFactor * ybounds2[2 * i2 + 1]) };
      if (latBounds[0] > latBounds[1]) std::swap(latBounds[0], latBounds[1]);

      size_t numSearchCells = 0;
      for (size_t i1 = 0; i1 < gridsize1; ++i1)
        {
          if (ymin1[i1] < latBounds[1] && ymax1[i1] > latBounds[0]) numSearchCells++;
        }

      remapIndices[i2].resize(numSearchCells);
      size_t n = 0;
      for (size_t i1 = 0; i1 < gridsize1; ++i1)
        {
          if (ymin1[i1] < latBounds[1] && ymax1[i1] > latBounds[0]) remapIndices[i2][n++] = i1;
        }

      // printf("lat %zu found %zu of %zu\n", i2 + 1, numSearchCells, gridsize1);
    }

  enum yac_edge_type lonlatCircleType[] = { LON_CIRCLE_EDGE, LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE, LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE };
  Varray<enum yac_edge_type> greatCircleType(nv1, GREAT_CIRCLE_EDGE);

  auto srcEdgeType = greatCircleType.data();

  if (nv1 == 4)
    {
      auto lonlatCircleIndex = get_lonlat_circle_index(RemapGridType::Undefined, gridsize1, nv1, xbounds1, ybounds1);
      if (lonlatCircleIndex >= 0) srcEdgeType = &lonlatCircleType[lonlatCircleIndex];
    }

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) default(shared)
#endif
  for (size_t i2 = 0; i2 < ysize2; ++i2)
    {
      auto normOpt(NormOpt::FRACAREA);
      double latBounds[2] = { ybounds2[2 * i2], ybounds2[2 * i2 + 1] };
      if (latBounds[0] > latBounds[1]) std::swap(latBounds[0], latBounds[1]);

      auto tgtCellArea = latitude_area(latBounds[0], latBounds[1]);
      // printf("tgtCellArea %zu %g\n", i2 + 1, tgtCellArea);

      auto numSearchCells = remapIndices[i2].size();

      CellSearch search;
      search.numCellCorners = nv1;
      search.edgeType = srcEdgeType;

      cellsearch_realloc(numSearchCells, search);

      for (size_t j = 0; j < numSearchCells; ++j)
        {
          auto cellIndex = remapIndices[i2][j];
          auto xyz = search.gridCells[j].coordinates_xyz;
          auto cell_corner_lon = &xbounds1[cellIndex * nv1];
          auto cell_corner_lat = &ybounds1[cellIndex * nv1];
          for (size_t i = 0; i < nv1; ++i) gcLLtoXYZ(cell_corner_lon[i], cell_corner_lat[i], xyz[i]);
        }

      auto &partialWeights = search.partialAreas;
      auto &overlapCells = search.overlapCells;

      // Do the clipping and get the cell for the overlapping area
      yac_cell_lat_clipping(numSearchCells, search.gridCells.data(), latBounds, overlapCells.data());

      // Get the partial areas for the overlapping regions
      for (size_t i = 0; i < numSearchCells; ++i) partialWeights[i] = gridcell_area(overlapCells[i]);

      auto numWeights = remove_invalid_areas(numSearchCells, partialWeights, remapIndices[i2]);
      // printf("numWeights: %zu %zu\n", numSearchCells, numWeights);

      if (normOpt == NormOpt::FRACAREA) correct_weights(tgtCellArea, numWeights, partialWeights);

      numWeights = remove_invalid_weights(gridsize1, numWeights, partialWeights, remapIndices[i2]);

      remapWeights[i2].resize(numWeights);
      for (size_t i = 0; i < numWeights; ++i) remapWeights[i2][i] = partialWeights[i];

      cellsearch_free(search);
    }
}

template <typename T1, typename T2>
static size_t
remap_zonal_mean(const Varray2D<size_t> &remapIndices, const Varray2D<double> &remapWeights, const Varray<T1> &srcArray,
                 Varray<T2> &tgtArray, T1 missval)
{
  auto ysize2 = remapIndices.size();

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) default(shared)
#endif
  for (size_t i2 = 0; i2 < ysize2; ++i2)
    {
      auto numWeights = remapWeights[i2].size();
      Varray<size_t> indices(numWeights);
      for (size_t i = 0; i < numWeights; ++i) indices[i] = remapIndices[i2][i];
      Varray<double> partialWeights(numWeights);
      for (size_t i = 0; i < numWeights; ++i) partialWeights[i] = remapWeights[i2][i];

      numWeights = remove_missing_weights(srcArray, missval, numWeights, partialWeights, indices);

      tgtArray[i2] = missval;

      if (numWeights)
        {
          auto normOpt(NormOpt::FRACAREA);
          auto tgtCellArea = 0.0;  // not needed for NormOpt::FRACAREA
          auto tgtCellFrac = varray_sum(numWeights, partialWeights);
          sort_weights(numWeights, indices, partialWeights);
          // Normalize weights using cell target area if requested
          normalize_weights(normOpt, tgtCellArea, tgtCellFrac, numWeights, partialWeights);
          tgtArray[i2] = conserv_remap(srcArray, numWeights, partialWeights, indices);
        }
    }

  size_t nmiss = 0;
  for (size_t i2 = 0; i2 < ysize2; ++i2)
    if (dbl_is_equal(tgtArray[i2], missval)) nmiss++;

  return nmiss;
}

void
remap_zonal_mean(const Varray2D<size_t> &remapIndices, const Varray2D<double> &remapWeights, const Field &field1, Field &field2)
{
  // clang-format off
  if      (field1.memType == MemType::Float && field2.memType == MemType::Float)
    field2.nmiss = remap_zonal_mean(remapIndices, remapWeights, field1.vec_f, field2.vec_f, (float)field1.missval);
  else if (field1.memType == MemType::Float && field2.memType == MemType::Double)
    field2.nmiss = remap_zonal_mean(remapIndices, remapWeights, field1.vec_f, field2.vec_d, (float)field1.missval);
  else if (field1.memType == MemType::Double && field2.memType == MemType::Float)
    field2.nmiss = remap_zonal_mean(remapIndices, remapWeights, field1.vec_d, field2.vec_f, field1.missval);
  else
    field2.nmiss = remap_zonal_mean(remapIndices, remapWeights, field1.vec_d, field2.vec_d, field1.missval);
  // clang-format on
}
