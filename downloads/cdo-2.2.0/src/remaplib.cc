/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/
/*
  This is a C library of the Fortran SCRIP version 1.4

  ===>>> Please send bug reports to <https://mpimet.mpg.de/cdo> <<<===

  Spherical Coordinate Remapping and Interpolation Package (SCRIP)
  ================================================================

  SCRIP is a software package which computes addresses and weights for
  remapping and interpolating fields between grids in spherical coordinates.
  It was written originally for remapping fields to other grids in a coupled
  climate model, but is sufficiently general that it can be used in other
  applications as well. The package should work for any grid on the surface
  of a sphere. SCRIP currently supports four remapping options:

  Conservative remapping
  ----------------------
  First- and second-order conservative remapping as described in
  Jones (1999, Monthly Weather Review, 127, 2204-2210).

  Bilinear interpolation
  ----------------------
  Slightly generalized to use a local bilinear approximation
  (only logically-rectangular grids).

  Bicubic interpolation
  ----------------------
  Similarly generalized (only logically-rectangular grids).

  Distance-weighted averaging
  ---------------------------
  Distance-weighted average of a user-specified number of nearest neighbor
  values.

  Documentation
  =============

  http://climate.lanl.gov/Software/SCRIP/SCRIPusers.pdf

*/
/*
  2013-11-08 Uwe Schulzweida: split remapgrid class to srcGrid and tgtGrid
  2012-01-16 Uwe Schulzweida: alloc grid2_bound_box only for conservative remapping
  2011-01-07 Uwe Schulzweida: Changed remap weights from 2D to 1D array
  2009-05-25 Uwe Schulzweida: Changed restrict data type from double to int
  2009-01-11 Uwe Schulzweida: OpenMP parallelization
 */

#include <cdi.h>

#include "cdo_options.h"
#include "cimdOmp.h"
#include "process_int.h"
#include "cdo_wtime.h"
#include <mpim_grid.h>
#include "gridreference.h"
#include "remap.h"

static bool remap_gen_weights = true;
static bool remap_write_remap = false;
static int remap_num_srch_bins = 180;

constexpr size_t DEFAULT_MAX_ITER = 100;
size_t remap_max_iter = DEFAULT_MAX_ITER;  // Max iteration count for i, j iteration

int
remap_check_mask_indices(const size_t (&indices)[4], const Varray<short> &mask)
{
  int searchResult = 1;
  if (mask.size() > 0)
    {
      for (int n = 0; n < 4; ++n)
        if (mask[indices[n]] == 0) searchResult = 0;
    }
  return searchResult;
}

void
remap_set_int(int remapvar, int value)
{
  // clang-format off
  if      (remapvar == REMAP_WRITE_REMAP)   remap_write_remap = (value > 0);
  else if (remapvar == REMAP_MAX_ITER)      remap_max_iter = value;
  else if (remapvar == REMAP_NUM_SRCH_BINS) remap_num_srch_bins = value;
  else if (remapvar == REMAP_GENWEIGHTS)    remap_gen_weights = (value > 0);
  else    cdo_abort("Unsupported remap variable (%d)!", remapvar);
  // clang-format on
}

/*****************************************************************************/
static void
boundboxFromCorners(size_t size, size_t nc, const Varray<double> &corner_lon, const Varray<double> &corner_lat, float *bound_box)
{
#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
  for (size_t i = 0; i < size; ++i)
    {
      size_t i4 = (i << 2);  // *4
      size_t inc = i * nc;
      float clat = corner_lat[inc];
      float clon = corner_lon[inc];
      bound_box[i4 + 0] = clat;
      bound_box[i4 + 1] = clat;
      bound_box[i4 + 2] = clon;
      bound_box[i4 + 3] = clon;
      for (size_t j = 1; j < nc; ++j)
        {
          clat = corner_lat[inc + j];
          clon = corner_lon[inc + j];
          if (clat < bound_box[i4 + 0]) bound_box[i4 + 0] = clat;
          if (clat > bound_box[i4 + 1]) bound_box[i4 + 1] = clat;
          if (clon < bound_box[i4 + 2]) bound_box[i4 + 2] = clon;
          if (clon > bound_box[i4 + 3]) bound_box[i4 + 3] = clon;
        }
    }
}

static void
boundboxFromCenter(bool lonIsCyclic, size_t size, size_t nx, size_t ny, const Varray<double> &center_lon,
                   const Varray<double> &center_lat, float *bound_box)
{
#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
  for (size_t n = 0; n < size; ++n)
    {
      size_t n4 = (n << 2);

      // Find N,S and NE points to this grid point

      size_t j = n / nx;
      size_t i = n - j * nx;

      size_t ip1 = (i < (nx - 1)) ? i + 1 : lonIsCyclic ? 0 : i;
      size_t jp1 = (j < (ny - 1)) ? j + 1 : j;

      size_t idx[4];
      idx[0] = n;
      idx[1] = j * nx + ip1;    // east
      idx[2] = jp1 * nx + ip1;  // north-east
      idx[3] = jp1 * nx + i;    // north

      // Find N,S and NE lat/lon coords and check bounding box

      float tmp_lats[4], tmp_lons[4];
      for (int k = 0; k < 4; ++k) tmp_lons[k] = center_lon[idx[k]];
      for (int k = 0; k < 4; ++k) tmp_lats[k] = center_lat[idx[k]];

      bound_box[n4 + 0] = tmp_lats[0];
      bound_box[n4 + 1] = tmp_lats[0];
      bound_box[n4 + 2] = tmp_lons[0];
      bound_box[n4 + 3] = tmp_lons[0];

      for (int k = 1; k < 4; ++k)
        {
          if (tmp_lats[k] < bound_box[n4 + 0]) bound_box[n4 + 0] = tmp_lats[k];
          if (tmp_lats[k] > bound_box[n4 + 1]) bound_box[n4 + 1] = tmp_lats[k];
          if (tmp_lons[k] < bound_box[n4 + 2]) bound_box[n4 + 2] = tmp_lons[k];
          if (tmp_lons[k] > bound_box[n4 + 3]) bound_box[n4 + 3] = tmp_lons[k];
        }
    }
}

LonLatPoint
remapgrid_get_lonlat(RemapGrid *grid, size_t index)
{
  double lon, lat;

  if (grid->type == RemapGridType::Reg2D)
    {
      auto nx = grid->dims[0];
      auto iy = index / nx;
      auto ix = index - iy * nx;
      lat = grid->reg2d_center_lat[iy];
      lon = grid->reg2d_center_lon[ix];
      if (lon < 0) lon += PI2;
    }
  else if (grid->type == RemapGridType::HealPix)
    {
      hp_index_to_lonlat(grid->order, grid->nside, index, &lon, &lat);
      //if (lon < 0) lon += PI2;
    }
  else
    {
      lat = grid->cell_center_lat[index];
      lon = grid->cell_center_lon[index];
    }

  return LonLatPoint(lon, lat);
}

void
check_lon_range(const char *txt, size_t nlons, Varray<double> &lons)
{
  assert(!lons.empty());

  if (txt)
    {
      double minval = 1.e36;
      double maxval = -1.e36;
#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static) reduction(min : minval) reduction(max : maxval)
#endif
      for (size_t i = 0; i < nlons; ++i)
        {
          minval = std::min(minval, lons[i]);
          maxval = std::max(maxval, lons[i]);
        }
      if (minval < -PI2 || maxval > 2 * PI2)
        cdo_warning("%s grid cell center longitudes out of range (min=%.3g/max=%.3g)!", txt, RAD2DEG * minval, RAD2DEG * maxval);
    }

#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
  for (size_t i = 0; i < nlons; ++i)
    {
      // remove missing values
      if (lons[i] < -PI2) lons[i] = 0;
      if (lons[i] > 2 * PI2) lons[i] = PI2;

      if (lons[i] > PI2) lons[i] -= PI2;
      if (lons[i] < 0.0) lons[i] += PI2;
    }
}

void
check_lat_range(const char *txt, size_t nlats, Varray<double> &lats)
{
  assert(!lats.empty());

  if (txt)
    {
      double minval = 1.e36;
      double maxval = -1.e36;
#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static) reduction(min : minval) reduction(max : maxval)
#endif
      for (size_t i = 0; i < nlats; ++i)
        {
          minval = std::min(minval, lats[i]);
          maxval = std::max(maxval, lats[i]);
        }
      if (minval < -(PIH + 0.0001) || maxval > (PIH + 0.0001))
        cdo_warning("%s grid cell center latitudes out of range (min=%.3g/max=%.3g)!", txt, RAD2DEG * minval, RAD2DEG * maxval);
    }

#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
  for (size_t i = 0; i < nlats; ++i)
    {
      if (lats[i] > PIH) lats[i] = PIH;
      if (lats[i] < -PIH) lats[i] = -PIH;
    }
}

static void
check_lon_boundbox_range(size_t nlons, float *bound_box)
{
  assert(bound_box != nullptr);

#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
  for (size_t n = 0; n < nlons; ++n)
    {
      auto n4 = n << 2;
      if (std::fabs(bound_box[n4 + 3] - bound_box[n4 + 2]) > PI_f)
        {
          bound_box[n4 + 2] = 0.0f;
          bound_box[n4 + 3] = PI2_f;
        }
    }
}

static void
check_lat_boundbox_range(size_t nlats, float *bound_box, Varray<double> &lats)
{
  assert(bound_box != nullptr);

#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
  for (size_t n = 0; n < nlats; ++n)
    {
      auto n4 = n << 2;
      if ((float) lats[n] < bound_box[n4]) bound_box[n4] = -PIH_f;
      if ((float) lats[n] > bound_box[n4 + 1]) bound_box[n4 + 1] = PIH_f;
    }
}

/*****************************************************************************/

static void
grid_check_lat_borders_rad(size_t n, Varray<double> &ybounds)
{
  constexpr double YMAX = PIH;
  constexpr double YLIM = 88 * DEG2RAD;
  auto lrev = (ybounds[0] > ybounds[n - 1]);
  if (lrev)
    {
      if (ybounds[0] > ybounds[1])
        {
          if (ybounds[0] > YLIM) ybounds[0] = YMAX;
          if (ybounds[n - 1] < -YLIM) ybounds[n - 1] = -YMAX;
        }
      else
        {
          if (ybounds[1] > YLIM) ybounds[1] = YMAX;
          if (ybounds[n - 2] < -YLIM) ybounds[n - 2] = -YMAX;
        }
    }
  else
    {
      if (ybounds[0] < ybounds[1])
        {
          if (ybounds[0] < -YLIM) ybounds[0] = -YMAX;
          if (ybounds[n - 1] > YLIM) ybounds[n - 1] = YMAX;
        }
      else
        {
          if (ybounds[1] < -YLIM) ybounds[1] = -YMAX;
          if (ybounds[n - 2] > YLIM) ybounds[n - 2] = YMAX;
        }
    }
}

static void
convert_bounds_reg2d(size_t n, const Varray<double> &boundsIn, Varray<double> &boundsOut)
{
  auto lrev = (boundsIn[0] > boundsIn[2 * n - 1]);
  if (boundsIn[0] > boundsIn[1]) lrev = !lrev;
  if (lrev)
    {
      boundsOut[0] = boundsIn[1];
      for (size_t i = 0; i < n; ++i) boundsOut[i + 1] = boundsIn[2 * i];
    }
  else
    {
      boundsOut[0] = boundsIn[0];
      for (size_t i = 0; i < n; ++i) boundsOut[i + 1] = boundsIn[2 * i + 1];
    }
}

static void
remap_define_reg2d(int gridID, RemapGrid &grid, bool conservMapping, const char *txt)
{
  auto nx = grid.dims[0];
  auto ny = grid.dims[1];
  auto nxp1 = nx + 1;
  auto nyp1 = ny + 1;

  auto nxm = nx;
  if (grid.isCyclic) nxm++;

  if (grid.size != nx * ny) cdo_abort("Internal error, wrong dimensions!");

  grid.reg2d_center_lon.resize(nxm);
  grid.reg2d_center_lat.resize(ny);

  grid.reg2d_center_lon[0] = 0.0;
  grid.reg2d_center_lat[0] = 0.0;
  gridInqXvals(gridID, grid.reg2d_center_lon.data());
  gridInqYvals(gridID, grid.reg2d_center_lat.data());

  // Convert lat/lon units if required
  cdo_grid_to_radian(gridID, CDI_XAXIS, nx, grid.reg2d_center_lon.data(), "grid reg2d center lon");
  cdo_grid_to_radian(gridID, CDI_YAXIS, ny, grid.reg2d_center_lat.data(), "grid reg2d center lat");

  if (grid.reg2d_center_lon[nx - 1] < grid.reg2d_center_lon[0])
    for (size_t i = 1; i < nx; ++i)
      if (grid.reg2d_center_lon[i] < grid.reg2d_center_lon[i - 1]) grid.reg2d_center_lon[i] += PI2;

  if (grid.isCyclic) grid.reg2d_center_lon[nx] = grid.reg2d_center_lon[0] + PI2;

  grid.reg2d_corner_lon.resize(nxp1);
  grid.reg2d_corner_lat.resize(nyp1);

  if (gridInqXbounds(gridID, nullptr))
    {
      Varray<double> xbounds(2 * nx);
      gridInqXbounds(gridID, xbounds.data());
      convert_bounds_reg2d(nx, xbounds, grid.reg2d_corner_lon);
      cdo_grid_to_radian(gridID, CDI_XAXIS, nx + 1, grid.reg2d_corner_lon.data(), "grid reg2d corner lon");
    }
  else
    {
      if (conservMapping && nx == 1) cdo_abort("Longitude bounds of %s grid missing!", txt);
      grid_gen_corners(nx, grid.reg2d_center_lon.data(), grid.reg2d_corner_lon.data());
    }

  if (gridInqYbounds(gridID, nullptr))
    {
      Varray<double> ybounds(2 * ny);
      gridInqYbounds(gridID, ybounds.data());
      convert_bounds_reg2d(ny, ybounds, grid.reg2d_corner_lat);
      cdo_grid_to_radian(gridID, CDI_YAXIS, ny + 1, grid.reg2d_corner_lat.data(), "grid reg2d corner lat");
    }
  else
    {
      if (conservMapping && ny == 1) cdo_abort("Latitude bounds of %s grid missing!", txt);
      grid_gen_corners(ny, grid.reg2d_center_lat.data(), grid.reg2d_corner_lat.data());
      grid_check_lat_borders_rad(ny + 1, grid.reg2d_corner_lat);
    }
}

static void
init_mask(int gridID, RemapGrid &grid)
{
#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
  for (size_t i = 0; i < grid.size; ++i) grid.mask[i] = 1;

  if (gridInqMask(gridID, nullptr))
    {
      std::vector<int> mask(grid.size);
      gridInqMask(gridID, &mask[0]);
      for (size_t i = 0; i < grid.size; ++i)
        if (mask[i] == 0) grid.mask[i] = 0;
    }
}

static void
remap_define_grid(RemapMethod mapType, int gridID, RemapGrid &grid, const char *txt)
{
  bool destroyGrid = false;
  int gridID_gme = -1;

  auto gridtype = gridInqType(grid.gridID);
  auto isHealpixGrid = (grid.type == RemapGridType::HealPix);

  if (isHealpixGrid)
    {
      auto [nside, order] = cdo::get_healpix_params(gridID);
      grid.nside = nside;
      grid.order = order;
    }

  if (gridtype != GRID_UNSTRUCTURED && gridtype != GRID_CURVILINEAR && !isHealpixGrid)
    {
      if (gridtype == GRID_GME)
        {
          gridID_gme = gridToUnstructured(grid.gridID, NeedCorners::Yes);
          grid.nvgp = gridInqSize(gridID_gme);
          gridID = gridDuplicate(gridID_gme);
          gridCompress(gridID);
          grid.useCellCorners = true;
        }
      else if (gridtype == GRID_GAUSSIAN_REDUCED || is_healpix_grid(gridID))
        {
          destroyGrid = true;
          gridID = gridToUnstructured(grid.gridID, NeedCorners::Yes);
        }
      else if (remap_write_remap || grid.type != RemapGridType::Reg2D)
        {
          destroyGrid = true;
          gridID = gridToCurvilinear(grid.gridID, NeedCorners::Yes);
        }
    }

  grid.size = gridInqSize(gridID);

  grid.dims[0] = isHealpixGrid ? grid.size : gridInqXsize(gridID);
  grid.dims[1] = gridInqYsize(gridID);
  if (gridtype != GRID_UNSTRUCTURED && gridInqType(gridID) != GRID_UNSTRUCTURED && !isHealpixGrid)
    {
      if (grid.dims[0] == 0) cdo_abort("%s grid without longitude coordinates!", gridNamePtr(gridtype));
      if (grid.dims[1] == 0) cdo_abort("%s grid without latitude coordinates!", gridNamePtr(gridtype));
    }

  grid.isCyclic = (gridIsCircular(gridID) > 0);

  grid.rank = (gridInqType(gridID) == GRID_UNSTRUCTURED || isHealpixGrid) ? 1 : 2;

  grid.num_cell_corners = (gridInqType(gridID) == GRID_UNSTRUCTURED) ? gridInqNvertex(gridID) : 4;

  remap_grid_alloc(mapType, grid);

  // Initialize logical mask
  init_mask(gridID, grid);

  if (!remap_write_remap && grid.type == RemapGridType::Reg2D) return;

  if (!isHealpixGrid && !gridHasCoordinates(gridID))
    cdo_abort("%s grid cell center coordinates missing!", txt);

  gridInqXvals(gridID, grid.cell_center_lon.data());
  gridInqYvals(gridID, grid.cell_center_lat.data());

  if (grid.needCellCorners)
    {
      if (!gridHasBounds(gridID)) cdo_abort("%s grid cell corner coordinates missing!", txt);

      gridInqXbounds(gridID, grid.cell_corner_lon.data());
      gridInqYbounds(gridID, grid.cell_corner_lat.data());
    }

  if (gridInqType(grid.gridID) == GRID_GME) gridInqMaskGME(gridID_gme, &grid.vgpm[0]);

  // Convert lat/lon units if required
  if (!isHealpixGrid)
    {
      cdo_grid_to_radian(gridID, CDI_XAXIS, grid.size, grid.cell_center_lon.data(), "grid center lon");
      cdo_grid_to_radian(gridID, CDI_YAXIS, grid.size, grid.cell_center_lat.data(), "grid center lat");
      if (grid.num_cell_corners && grid.needCellCorners)
        {
          cdo_grid_to_radian(gridID, CDI_XAXIS, grid.num_cell_corners * grid.size, grid.cell_corner_lon.data(), "grid corner lon");
          cdo_grid_to_radian(gridID, CDI_YAXIS, grid.num_cell_corners * grid.size, grid.cell_corner_lat.data(), "grid corner lat");
        }

      // Convert longitudes to 0,2pi interval
      check_lon_range(txt, grid.size, grid.cell_center_lon);
      if (grid.num_cell_corners && grid.needCellCorners)
        check_lon_range(nullptr, grid.num_cell_corners * grid.size, grid.cell_corner_lon);

      // Make sure input latitude range is within the machine values for +/- pi/2
      check_lat_range(txt, grid.size, grid.cell_center_lat);
      if (grid.num_cell_corners && grid.needCellCorners)
        check_lat_range(nullptr, grid.num_cell_corners * grid.size, grid.cell_corner_lat);
    }

  if (destroyGrid) gridDestroy(gridID);
}

// Compute bounding boxes for restricting future grid searches
static void
cell_bounding_boxes(RemapGrid &grid, float *cell_bound_box, int remap_grid_basis)
{
  if (grid.useCellCorners)
    {
      if (grid.needCellCorners)
        {
          if (Options::cdoVerbose) cdo_print("Grid: boundboxFromCorners");
          boundboxFromCorners(grid.size, grid.num_cell_corners, grid.cell_corner_lon, grid.cell_corner_lat, cell_bound_box);
        }
      else  // full grid search
        {
          if (Options::cdoVerbose) cdo_print("Grid: bounds missing -> full grid search!");

          auto gridsize = grid.size;
          for (size_t i = 0; i < gridsize; ++i)
            {
              cell_bound_box[i * 4] = -PIH_f;
              cell_bound_box[i * 4 + 1] = PIH_f;
              cell_bound_box[i * 4 + 2] = 0.0f;
              cell_bound_box[i * 4 + 3] = PI2_f;
            }
        }
    }
  else if (remap_grid_basis == REMAP_GRID_BASIS_SRC)
    {
      if (Options::cdoVerbose) cdo_print("Grid: boundboxFromCenter");
      if (grid.rank != 2) cdo_abort("Internal problem, grid rank = %d!", grid.rank);

      auto nx = grid.dims[0];
      auto ny = grid.dims[1];
      boundboxFromCenter(grid.isCyclic, grid.size, nx, ny, grid.cell_center_lon, grid.cell_center_lat, cell_bound_box);
    }

  if (remap_grid_basis == REMAP_GRID_BASIS_SRC || grid.needCellCorners) check_lon_boundbox_range(grid.size, cell_bound_box);

  // Try to check for cells that overlap poles
  if (remap_grid_basis == REMAP_GRID_BASIS_SRC || grid.needCellCorners)
    check_lat_boundbox_range(grid.size, cell_bound_box, grid.cell_center_lat);
}

void
remap_grid_alloc(RemapMethod mapType, RemapGrid &grid)
{
  if (grid.nvgp) grid.vgpm.resize(grid.nvgp);

  // only needed for srcGrid and remap_gen_weights
  grid.mask.resize(grid.size);

  if (remap_write_remap || (grid.type != RemapGridType::Reg2D && grid.type != RemapGridType::HealPix))
    {
      grid.cell_center_lon.resize(grid.size);
      grid.cell_center_lat.resize(grid.size);
    }

  auto needCellarea = (mapType == RemapMethod::CONSERV_SCRIP || mapType == RemapMethod::CONSERV);
  if (needCellarea) grid.cell_area.resize(grid.size, 0.0);

  if (remap_gen_weights || mapType == RemapMethod::CONSERV) { grid.cell_frac.resize(grid.size, 0.0); }

  if (grid.needCellCorners && grid.num_cell_corners > 0)
    {
      auto nalloc = grid.num_cell_corners * grid.size;
      grid.cell_corner_lon.resize(nalloc, 0);
      grid.cell_corner_lat.resize(nalloc, 0);
    }
}

void
remap_grid_free(RemapGrid &grid)
{
  varray_free(grid.vgpm);
  varray_free(grid.mask);

  varray_free(grid.reg2d_center_lat);
  varray_free(grid.reg2d_center_lon);
  varray_free(grid.reg2d_corner_lat);
  varray_free(grid.reg2d_corner_lon);

  varray_free(grid.cell_center_lat);
  varray_free(grid.cell_center_lon);
  varray_free(grid.cell_corner_lat);
  varray_free(grid.cell_corner_lon);

  varray_free(grid.cell_area);
  varray_free(grid.cell_frac);

  if (grid.tmpgridID != -1) gridDestroy(grid.tmpgridID);
}

static void
check_for_convex_cells(RemapGrid &tgtGrid)
{
  if (tgtGrid.type == RemapGridType::Reg2D) return;

  auto numCellCorners = tgtGrid.num_cell_corners;
  if (numCellCorners <= 4) return;

  auto numCells = tgtGrid.size;
  if (numCells > 1000) numCells = 1000;

  for (size_t i = 0; i < numCells; ++i) {}
}

void
remap_search_init(RemapMethod mapType, RemapSearch &search, RemapGrid &srcGrid, RemapGrid &tgtGrid)
{
  extern PointSearchMethod pointSearchMethod;
  extern CellSearchMethod cellSearchMethod;

  search.srcGrid = &srcGrid;
  search.tgtGrid = &tgtGrid;

  search.srcBins.ncells = srcGrid.size;
  search.tgtBins.ncells = tgtGrid.size;

  search.srcBins.nbins = remap_num_srch_bins;
  search.tgtBins.nbins = remap_num_srch_bins;

  auto usePointsearch = (mapType == RemapMethod::DISTWGT);
  if (srcGrid.type != RemapGridType::Reg2D && pointSearchMethod != PointSearchMethod::latbins)
    {
      usePointsearch |= (mapType == RemapMethod::BILINEAR);
      usePointsearch |= (mapType == RemapMethod::BICUBIC);
    }

  auto useCellsearch = (mapType == RemapMethod::CONSERV)
                       && (cellSearchMethod == CellSearchMethod::spherepart || srcGrid.type == RemapGridType::Reg2D);

  const char *searchMethodStr = NULL;
  auto start = Options::cdoVerbose ? cdo_get_wtime() : 0.0;

  if (usePointsearch)
    {
      searchMethodStr = "Point search";
      auto xIsCyclic = srcGrid.isCyclic;
      if (srcGrid.type == RemapGridType::HealPix)
        grid_point_search_create_healpix(search.gps, srcGrid.size, srcGrid.nside, srcGrid.order);
      else if (srcGrid.type == RemapGridType::Reg2D)
        grid_point_search_create_reg_2d(search.gps, xIsCyclic, srcGrid.dims, srcGrid.reg2d_center_lon, srcGrid.reg2d_center_lat);
      else
        grid_point_search_create(search.gps, xIsCyclic, srcGrid.dims, srcGrid.size, srcGrid.cell_center_lon,
                                 srcGrid.cell_center_lat);

      if (srcGrid.doExtrapolate) grid_point_search_extrapolate(search.gps);
    }
  else if (useCellsearch)
    {
      searchMethodStr = "Cell search";
      if (srcGrid.type == RemapGridType::Reg2D)
        grid_cell_search_create_reg_2d(search.gcs, srcGrid.dims, srcGrid.reg2d_corner_lon, srcGrid.reg2d_corner_lat);
      else
        grid_cell_search_create(search.gcs, srcGrid.size, srcGrid.num_cell_corners, srcGrid.cell_corner_lon,
                                srcGrid.cell_corner_lat);

      // check_for_convex_cells(tgtGrid);
    }
  else if (!(srcGrid.type == RemapGridType::Reg2D || tgtGrid.type == RemapGridType::Reg2D))
    {
      searchMethodStr = "Latitude bins";
      search.srcBins.cell_bound_box.resize(4 * srcGrid.size);
      if (tgtGrid.useCellCorners) search.tgtBins.cell_bound_box.resize(4 * tgtGrid.size);

      cell_bounding_boxes(srcGrid, &search.srcBins.cell_bound_box[0], REMAP_GRID_BASIS_SRC);
      cell_bounding_boxes(tgtGrid, &search.tgtBins.cell_bound_box[0], REMAP_GRID_BASIS_TGT);
      // Set up and assign address ranges to search bins in order to further restrict later searches
      calc_lat_bins(search.srcBins);
      if (mapType == RemapMethod::CONSERV_SCRIP || mapType == RemapMethod::CONSERV)
        {
          calc_lat_bins(search.tgtBins);
          varray_free(search.tgtBins.bin_lats);
          varray_free(search.srcBins.bin_lats);
          if (mapType == RemapMethod::CONSERV) varray_free(search.tgtBins.cell_bound_box);
        }
    }

  if (Options::cdoVerbose && searchMethodStr) cdo_print("%s created: %.2f seconds", searchMethodStr, cdo_get_wtime() - start);
}

void
remap_search_free(RemapSearch &search)
{
  varray_free(search.srcBins.bin_addr);
  varray_free(search.srcBins.bin_lats);
  varray_free(search.srcBins.cell_bound_box);

  varray_free(search.tgtBins.bin_addr);
  varray_free(search.tgtBins.bin_lats);
  varray_free(search.tgtBins.cell_bound_box);

  grid_point_search_delete(search.gps);
  grid_cell_search_delete(search.gcs);
}

void
remap_init_grids(RemapMethod mapType, bool doExtrapolate, int gridID1, RemapGrid &srcGrid, int gridID2, RemapGrid &tgtGrid)
{
  auto reg2d_srcGridID = gridID1;
  auto reg2d_tgtGridID = gridID2;

  if (mapType == RemapMethod::BILINEAR || mapType == RemapMethod::BICUBIC || mapType == RemapMethod::DISTWGT
      || mapType == RemapMethod::CONSERV)
    {
      if (is_reg2d_grid(gridID1))
        srcGrid.type = RemapGridType::Reg2D;
      else if (is_healpix_grid(gridID1) && mapType == RemapMethod::BILINEAR && (!remap_gen_weights || Options::test))
        srcGrid.type = RemapGridType::HealPix;
      else if (is_healpix_grid(gridID1) && mapType == RemapMethod::DISTWGT && !remap_gen_weights)
        srcGrid.type = RemapGridType::HealPix;
    }

  //#ifdef YAC_CELL_SEARCH
  // if (is_reg2d_grid(gridID2) && mapType == RemapMethod::CONSERV) tgtGrid.type = RemapGridType::Reg2D;
  //#else
  if (srcGrid.type == RemapGridType::Reg2D)
    {
      if (is_reg2d_grid(gridID2) && mapType == RemapMethod::CONSERV) tgtGrid.type = RemapGridType::Reg2D;
      // else srcGrid.type = -1;
    }
  //#endif

  if (!remap_gen_weights && is_reg2d_grid(gridID2) && tgtGrid.type != RemapGridType::Reg2D)
    {
      if (mapType == RemapMethod::DISTWGT) tgtGrid.type = RemapGridType::Reg2D;
      if (mapType == RemapMethod::BILINEAR && (srcGrid.type == RemapGridType::Reg2D || srcGrid.type == RemapGridType::HealPix))
        tgtGrid.type = RemapGridType::Reg2D;
    }

  if (!remap_gen_weights && is_healpix_grid(gridID2))
    {
      if (mapType == RemapMethod::BILINEAR || mapType == RemapMethod::DISTWGT)
        tgtGrid.type = RemapGridType::HealPix;
    }

  srcGrid.doExtrapolate = doExtrapolate;

  if (mapType == RemapMethod::CONSERV_SCRIP || mapType == RemapMethod::CONSERV)
    {
      if (srcGrid.type != RemapGridType::Reg2D && srcGrid.type != RemapGridType::HealPix)
        {
          srcGrid.useCellCorners = true;
          srcGrid.needCellCorners = true;
        }

      if (tgtGrid.type != RemapGridType::Reg2D)
        {
          tgtGrid.useCellCorners = true;
          tgtGrid.needCellCorners = true;
        }
    }

  srcGrid.gridID = gridID1;
  tgtGrid.gridID = gridID2;

  if (gridInqType(gridID1) == GRID_UNSTRUCTURED && !gridHasCoordinates(gridID1))
    {
      auto reference = dereferenceGrid(gridID1);
      if (reference.isValid) srcGrid.gridID = gridID1 = reference.gridID;
      if (reference.notFound) cdo_abort("Reference to source grid not found!");
    }

  if (gridInqType(gridID2) == GRID_UNSTRUCTURED && !gridHasCoordinates(gridID2))
    {
      auto reference = dereferenceGrid(gridID2);
      if (reference.isValid) tgtGrid.gridID = gridID2 = reference.gridID;
      if (reference.notFound) cdo_abort("Reference to target grid not found!");
    }

  auto sgridID = srcGrid.gridID;
  if (gridInqSize(sgridID) > 1 && gridProjIsSupported(sgridID) && srcGrid.type != RemapGridType::HealPix)
    {
      auto needCorners = srcGrid.needCellCorners ? NeedCorners::Yes : NeedCorners::No;
      if (is_healpix_grid(sgridID))
        gridID1 = gridToUnstructured(srcGrid.gridID, needCorners);
      else
        gridID1 = gridToCurvilinear(srcGrid.gridID, needCorners);
      srcGrid.gridID = gridID1;
      srcGrid.tmpgridID = srcGrid.gridID;
    }

  // if (srcGrid.type != RemapGridType::Reg2D)
  remap_define_grid(mapType, gridID1, srcGrid, "Source");
  remap_define_grid(mapType, gridID2, tgtGrid, "Target");

  auto conservMapping = (mapType == RemapMethod::CONSERV);
  if (srcGrid.type == RemapGridType::Reg2D) remap_define_reg2d(reg2d_srcGridID, srcGrid, conservMapping, "source");
  if (tgtGrid.type == RemapGridType::Reg2D) remap_define_reg2d(reg2d_tgtGridID, tgtGrid, conservMapping, "target");
}

/*****************************************************************************/

template <typename T1, typename T2>
static void
remap_stat(int remapOrder, RemapGrid &srcGrid, RemapGrid &tgtGrid, RemapVars &rv, const Varray<T1> &array1,
           const Varray<T2> &array2, double missval)
{
  T1 mv1 = missval;
  T2 mv2 = missval;

  cdo_print("%s order mapping from grid1 to grid2:", (remapOrder == 2) ? "Second" : "First");
  cdo_print("----------------------------------------------");

  auto mmm = varray_min_max_mean_mv(srcGrid.size, array1, mv1);
  cdo_print("  Grid1 min,mean,max: %g %g %g", mmm.min, mmm.mean, mmm.max);

  mmm = varray_min_max_mean_mv(tgtGrid.size, array2, mv2);
  cdo_print("  Grid2 min,mean,max: %g %g %g", mmm.min, mmm.mean, mmm.max);

  // Conservation Test

  if (srcGrid.cell_area.size())
    {
      cdo_print("  Conservation:");
      double sum = 0.0;
      for (size_t n = 0; n < srcGrid.size; ++n)
        if (!dbl_is_equal(array1[n], mv1)) sum += array1[n] * srcGrid.cell_area[n] * srcGrid.cell_frac[n];
      cdo_print("  Grid1 Integral = %g", sum);

      sum = 0;
      for (size_t n = 0; n < tgtGrid.size; ++n)
        if (!dbl_is_equal(array2[n], mv2)) sum += array2[n] * tgtGrid.cell_area[n] * tgtGrid.cell_frac[n];
      cdo_print("  Grid2 Integral = %g", sum);
      /*
      for ( n = 0; n < srcGrid.size; n++ )
       fprintf(stderr, "1 %d %g %g %g\n", n, array1[n], srcGrid.cell_area[n],
      srcGrid.cell_frac[n]); for ( n = 0; n < tgtGrid.size; n++ )
        fprintf(stderr, "2 %d %g %g %g\n", n, array2[n], tgtGrid.cell_area[n],
      tgtGrid.cell_frac[n]);
      */
    }

  cdo_print("  Number of weights %zu", rv.num_wts);
  cdo_print("  Number of sparse matrix entries %zu", rv.numLinks);
  cdo_print("  Total number of dest cells %zu", tgtGrid.size);

  std::vector<size_t> tgt_count(tgtGrid.size, 0);

#ifdef HAVE_OPENMP4
#pragma omp simd
#endif
  for (size_t n = 0; n < rv.numLinks; ++n) tgt_count[rv.tgtCellIndices[n]]++;

  size_t imin = SIZE_MAX;
  size_t imax = 0;
  for (size_t n = 0; n < tgtGrid.size; ++n)
    {
      if (tgt_count[n] > 0)
        {
          if (tgt_count[n] < imin) imin = tgt_count[n];
          if (tgt_count[n] > imax) imax = tgt_count[n];
        }
    }

  size_t idiff = (imax - imin) / 10 + 1;
  size_t icount = 0;
  for (size_t i = 0; i < tgtGrid.size; ++i)
    if (tgt_count[i] > 0) icount++;

  cdo_print("  Number of cells participating in remap %zu", icount);

  if (icount)
    {
      cdo_print("  Min no of entries/row = %zu", imin);
      cdo_print("  Max no of entries/row = %zu", imax);

      imax = imin + idiff;
      for (size_t n = 0; n < 10; ++n)
        {
          icount = 0;
          for (size_t i = 0; i < tgtGrid.size; ++i)
            if (tgt_count[i] >= imin && tgt_count[i] < imax) icount++;

          if (icount) cdo_print("  Num of rows with entries between %zu - %zu  %zu", imin, imax - 1, icount);

          imin = imin + idiff;
          imax = imax + idiff;
        }
    }

  if (rv.sort_add) cdo_print("  Sparse matrix entries are explicitly sorted.");
}

void
remap_stat(int remapOrder, RemapGrid &srcGrid, RemapGrid &tgtGrid, RemapVars &rv, const Field &field1, const Field &field2)
{
  if (memtype_is_float_float(field1.memType, field2.memType))
    remap_stat(remapOrder, srcGrid, tgtGrid, rv, field1.vec_f, field2.vec_f, field1.missval);
  else if (memtype_is_float_double(field1.memType, field2.memType))
    remap_stat(remapOrder, srcGrid, tgtGrid, rv, field1.vec_f, field2.vec_d, field1.missval);
  else if (memtype_is_double_float(field1.memType, field2.memType))
    remap_stat(remapOrder, srcGrid, tgtGrid, rv, field1.vec_d, field2.vec_f, field1.missval);
  else
    remap_stat(remapOrder, srcGrid, tgtGrid, rv, field1.vec_d, field2.vec_d, field1.missval);
}

/*****************************************************************************/

template <typename T>
void
remap_gradients(RemapGrid &grid, const Varray<short> &mask, const Varray<T> &array, RemapGradients &gradients)
{
  if (grid.rank != 2) cdo_abort("Internal problem (remap_gradients), grid rank = %d!", grid.rank);

  auto gridSize = grid.size;
  auto nx = grid.dims[0];
  auto ny = grid.dims[1];

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
  for (size_t n = 0; n < gridSize; ++n)
    {
      if (mask[n])
        {
          // clang-format off
          auto delew = 0.5;
          auto delns = 0.5;

          auto j = n / nx + 1;
          auto i = n - (j - 1) * nx + 1;

          auto ip1 = i + 1;
          auto im1 = i - 1;
          auto jp1 = j + 1;
          auto jm1 = j - 1;

          if (ip1 > nx) ip1 = ip1 - nx;
          if (im1 < 1)  im1 = nx;
          if (jp1 > ny) { jp1 = j; delns = 1.0; }
          if (jm1 < 1)  { jm1 = j; delns = 1.0; }

          auto in = (jp1 - 1) * nx + i - 1;
          auto is = (jm1 - 1) * nx + i - 1;
          auto ie = (j - 1) * nx + ip1 - 1;
          auto iw = (j - 1) * nx + im1 - 1;

          auto ine = (jp1 - 1) * nx + ip1 - 1;
          auto inw = (jp1 - 1) * nx + im1 - 1;
          auto ise = (jm1 - 1) * nx + ip1 - 1;
          auto isw = (jm1 - 1) * nx + im1 - 1;

          // Compute i-gradient
          if (!mask[ie]) { ie = n; delew = 1.0; }
          if (!mask[iw]) { iw = n; delew = 1.0; }

          gradients.grad_lat[n] = delew * (array[ie] - array[iw]);

          // Compute j-gradient
          if (!mask[in]) { in = n; delns = 1.0; }
          if (!mask[is]) { is = n; delns = 1.0; }

          gradients.grad_lon[n] = delns * (array[in] - array[is]);
          // clang-format on

          // Compute ij-gradient
          delew = 0.5;
          delns = (jp1 == j || jm1 == j) ? 1.0 : 0.5;

          if (!mask[ine])
            {
              if (in != n)
                {
                  ine = in;
                  delew = 1.0;
                }
              else if (ie != n)
                {
                  ine = ie;
                  inw = iw;
                  if (inw == n) delew = 1.0;
                  delns = 1.0;
                }
              else
                {
                  ine = n;
                  inw = iw;
                  delew = 1.0;
                  delns = 1.0;
                }
            }

          if (!mask[inw])
            {
              if (in != n)
                {
                  inw = in;
                  delew = 1.0;
                }
              else if (iw != n)
                {
                  inw = iw;
                  ine = ie;
                  if (ie == n) delew = 1.0;
                  delns = 1.0;
                }
              else
                {
                  inw = n;
                  ine = ie;
                  delew = 1.0;
                  delns = 1.0;
                }
            }

          auto grad_lat_zero = delew * (array[ine] - array[inw]);

          if (!mask[ise])
            {
              if (is != n)
                {
                  ise = is;
                  delew = 1.0;
                }
              else if (ie != n)
                {
                  ise = ie;
                  isw = iw;
                  if (isw == n) delew = 1.0;
                  delns = 1.0;
                }
              else
                {
                  ise = n;
                  isw = iw;
                  delew = 1.0;
                  delns = 1.0;
                }
            }

          if (!mask[isw])
            {
              if (is != n)
                {
                  isw = is;
                  delew = 1.0;
                }
              else if (iw != n)
                {
                  isw = iw;
                  ise = ie;
                  if (ie == n) delew = 1.0;
                  delns = 1.0;
                }
              else
                {
                  isw = n;
                  ise = ie;
                  delew = 1.0;
                  delns = 1.0;
                }
            }

          auto grad_lon_zero = delew * (array[ise] - array[isw]);
          gradients.grad_latlon[n] = delns * (grad_lat_zero - grad_lon_zero);
        }
      else
        {
          gradients.grad_lat[n] = 0.0;
          gradients.grad_lon[n] = 0.0;
          gradients.grad_latlon[n] = 0.0;
        }
    }
}  // remap_gradients

// Explicit instantiation
template void remap_gradients(RemapGrid &grid, const Varray<short> &mask, const Varray<float> &array, RemapGradients &gradients);
template void remap_gradients(RemapGrid &grid, const Varray<short> &mask, const Varray<double> &array, RemapGradients &gradients);

void
remap_gradients(RemapGrid &grid, const Field &field, RemapGradients &gradients)
{
  auto gridSize = grid.size;
  Varray<short> mask(gridSize);
#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
  for (size_t i = 0; i < gridSize; ++i) mask[i] = (grid.mask[i] > 0);

  if (field.memType == MemType::Float)
    remap_gradients(grid, mask, field.vec_f, gradients);
  else
    remap_gradients(grid, mask, field.vec_d, gradients);
}

/*****************************************************************************/

void
remap_check_area(size_t gridSize, const Varray<double> &cell_area, const char *name)
{
  for (size_t n = 0; n < gridSize; ++n)
    {
      if (cell_area[n] < -0.01) cdo_print("%s grid area error: %zu %g", name, n, cell_area[n]);
    }
}
