/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include "process_int.h"
#include <mpim_grid.h>
#include "remap.h"

int rect_grid_search2(long &imin, long &imax, double xmin, double xmax, long nxm, const Varray<double> &xm);

static long
set_srch_add(size_t numSearchCells, long nx, long imin, long imax, long jmin, long jmax, Varray<size_t> &srch_add)
{
  const size_t maxSize = numSearchCells + (jmax - jmin + 1) * (imax - imin + 1);
  if (srch_add.size() < maxSize) srch_add.resize(maxSize);
  for (long jm = jmin; jm <= jmax; ++jm)
    for (long im = imin; im <= imax; ++im) srch_add[numSearchCells++] = jm * nx + im;

  return numSearchCells;
}

static void
check_lon_bounds(double offset, double src_lon_min, double src_lon_max, double &bound_lon1, double &bound_lon2)
{
  bound_lon1 += offset;
  bound_lon2 += offset;
  if (bound_lon1 < src_lon_min && bound_lon2 > src_lon_min) bound_lon1 = src_lon_min;
  if (bound_lon2 > src_lon_max && bound_lon1 < src_lon_max) bound_lon2 = src_lon_max;
}

static void
debug_message(const char *txt, long imin, long imax, long jmin, long jmax, const Varray<double> &src_corner_lon)
{
  printf("%s:  lonMin=%g lonMax=%g  iMin=%ld iMax=%ld  jMin=%ld jMax %ld  numCells=%ld\n", txt, RAD2DEG * src_corner_lon[imin],
         RAD2DEG * src_corner_lon[imax + 1], imin, imax, jmin, jmax, (jmax - jmin + 1) * (imax - imin + 1));
}

static size_t
getSrchCellsReg2d(const size_t *src_grid_dims, const Varray<double> &src_corner_lat, const Varray<double> &src_corner_lon,
                  const double *tgt_cell_bound_box, Varray<size_t> &srch_add)
{
  auto debug = false;
  long nx = src_grid_dims[0];
  long ny = src_grid_dims[1];
  size_t numSearchCells = 0;  // num cells in restricted search arrays

  const long nxp1 = nx + 1;
  const long nyp1 = ny + 1;

  auto src_lon_min = src_corner_lon[0];
  auto src_lon_max = src_corner_lon[nx];

  long imin = nxp1, imax = -1, jmin = nyp1, jmax = -1;

  int lfound = rect_grid_search2(jmin, jmax, tgt_cell_bound_box[0], tgt_cell_bound_box[1], nyp1, src_corner_lat);
  if (!lfound) return 0;
  // printf("lfound, jmin, jmax %d %ld %ld\n", lfound, jmin, jmax);
  // if (jmin > 0) jmin--;
  // if (jmax < (ny-2)) jmax++;

  auto bound_lon1 = tgt_cell_bound_box[2];
  auto bound_lon2 = tgt_cell_bound_box[3];
  // debug = (bound_lon1 <= 0 && bound_lon2 >= 0);

  if (bound_lon1 <= src_lon_max && bound_lon2 >= src_lon_min)
    {
      check_lon_bounds(0.0, src_lon_min, src_lon_max, bound_lon1, bound_lon2);
      lfound = rect_grid_search2(imin, imax, bound_lon1, bound_lon2, nxp1, src_corner_lon);
      if (lfound)
        {
          if (debug) debug_message("1", imin, imax, jmin, jmax, src_corner_lon);
          numSearchCells = set_srch_add(numSearchCells, nx, imin, imax, jmin, jmax, srch_add);
        }
    }

  bound_lon1 = tgt_cell_bound_box[2];
  bound_lon2 = tgt_cell_bound_box[3];
  if (bound_lon1 < src_lon_min || bound_lon2 > src_lon_max)
    {
      if (bound_lon1 <= src_lon_min && bound_lon2 >= src_lon_min)
        {
          check_lon_bounds(2.0 * M_PI, src_lon_min, src_lon_max, bound_lon1, bound_lon2);
          long imin2 = nxp1, imax2 = -1;
          lfound = rect_grid_search2(imin2, imax2, bound_lon1, bound_lon2, nxp1, src_corner_lon);
          if (lfound)
            {
              if (imax != -1 && imin2 <= imax) imin2 = imax + 1;
              if (imax != -1 && imax2 <= imax) imax2 = imax + 1;
              if (imin2 >= 0 && imax2 < nx)
                {
                  if (debug) debug_message("2", imin2, imax2, jmin, jmax, src_corner_lon);
                  numSearchCells = set_srch_add(numSearchCells, nx, imin2, imax2, jmin, jmax, srch_add);
                }
            }
        }

      bound_lon1 = tgt_cell_bound_box[2];
      bound_lon2 = tgt_cell_bound_box[3];
      if (bound_lon1 <= src_lon_max && bound_lon2 >= src_lon_max)
        {
          check_lon_bounds(-2.0 * M_PI, src_lon_min, src_lon_max, bound_lon1, bound_lon2);
          long imin3 = nxp1, imax3 = -1;
          lfound = rect_grid_search2(imin3, imax3, bound_lon1, bound_lon2, nxp1, src_corner_lon);
          if (lfound)
            {
              if (imin != nxp1 && imin3 >= imin) imin3 = imin - 1;
              if (imax != nxp1 && imax3 >= imin) imax3 = imin - 1;
              if (imin3 >= 0 && imin3 < nx)
                {
                  if (debug) debug_message("3", imin3, imax3, jmin, jmax, src_corner_lon);
                  numSearchCells = set_srch_add(numSearchCells, nx, imin3, imax3, jmin, jmax, srch_add);
                }
            }
        }
    }

  if (debug) printf(" numSearchCells: %zu\n", numSearchCells);

  return numSearchCells;
}

static void
restrict_boundbox(const double *grid_bound_box, double *bound_box)
{
  if (bound_box[0] < grid_bound_box[0] && bound_box[1] > grid_bound_box[0]) bound_box[0] = grid_bound_box[0];
  if (bound_box[1] > grid_bound_box[1] && bound_box[0] < grid_bound_box[1]) bound_box[1] = grid_bound_box[1];

  if (bound_box[2] >= grid_bound_box[3] && (bound_box[3] - 2 * M_PI) > grid_bound_box[2])
    {
      bound_box[2] -= 2 * M_PI;
      bound_box[3] -= 2 * M_PI;
    }
  if (bound_box[3] <= grid_bound_box[2] && (bound_box[2] - 2 * M_PI) < grid_bound_box[3])
    {
      bound_box[2] += 2 * M_PI;
      bound_box[3] += 2 * M_PI;
    }
}

static void
boundboxFromCornersReg2d(const GridCell &gridCell, double *bound_box)
{
  const auto coordinates_x = gridCell.coordinates_x;
  const auto coordinates_y = gridCell.coordinates_y;

  const auto clat1 = coordinates_y[0];
  const auto clat2 = coordinates_y[2];

  bound_box[0] = (clat2 > clat1) ? clat1 : clat2;
  bound_box[1] = (clat2 > clat1) ? clat2 : clat1;
  bound_box[2] = coordinates_x[0];
  bound_box[3] = coordinates_x[1];
}

static void
boundboxFromCorners1(const GridCell &gridCell, double *bound_box)
{
  const auto coordinates_x = gridCell.coordinates_x;
  const auto coordinates_y = gridCell.coordinates_y;

  auto clon = coordinates_x[0];
  auto clat = coordinates_y[0];

  bound_box[0] = clat;
  bound_box[1] = clat;
  bound_box[2] = clon;
  bound_box[3] = clon;

  const auto nc = gridCell.yacGridCell.num_corners;
  for (size_t j = 1; j < nc; ++j)
    {
      clon = coordinates_x[j];
      clat = coordinates_y[j];

      if (clat < bound_box[0]) bound_box[0] = clat;
      if (clat > bound_box[1]) bound_box[1] = clat;
      if (clon < bound_box[2]) bound_box[2] = clon;
      if (clon > bound_box[3]) bound_box[3] = clon;
    }

  if (std::fabs(bound_box[3] - bound_box[2]) > PI)
    {
      bound_box[2] = 0;
      bound_box[3] = PI2;
    }
}

static size_t
gridSearchCellsReg2d(const GridCellSearch &gcs, bool isReg2dCell, GridCell &gridCell, Varray<size_t> &srchAddr)
{
  double tgt_cell_bound_box[4];
  if (isReg2dCell)
    boundboxFromCornersReg2d(gridCell, tgt_cell_bound_box);
  else
    boundboxFromCorners1(gridCell, tgt_cell_bound_box);

  restrict_boundbox(gcs.gridBoundboxReg2d, tgt_cell_bound_box);

  auto numSearchCells = getSrchCellsReg2d(gcs.dims, gcs.reg2d_corner_lat, gcs.reg2d_corner_lon, tgt_cell_bound_box, srchAddr);

  if (numSearchCells == 1 && gcs.dims[0] == 1 && gcs.dims[1] == 1 && IS_EQUAL(gcs.reg2d_corner_lat[0], gcs.reg2d_corner_lat[1])
      && IS_EQUAL(gcs.reg2d_corner_lon[0], gcs.reg2d_corner_lon[1]))
    numSearchCells = 0;

  return numSearchCells;
}

static size_t
gridSearchCells(GridCellSearch &gcs, bool isReg2dCell, GridCell &gridCell, Varray<size_t> &srchAddr)
{
  return do_grid_cell_search(gcs, isReg2dCell, gridCell, srchAddr);
}

size_t
remap_search_cells(RemapSearch &rsearch, bool isReg2dCell, GridCell &gridCell, Varray<size_t> &srchAddr)
{
  if (rsearch.srcGrid->type == RemapGridType::Reg2D)
    return gridSearchCellsReg2d(rsearch.gcs, isReg2dCell, gridCell, srchAddr);
  else
    return gridSearchCells(rsearch.gcs, isReg2dCell, gridCell, srchAddr);
}
