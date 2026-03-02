/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include "remap.h"

static void
nbr_add_and_dist(int &searchResult, size_t jj, size_t ii, size_t nx, double &dist_min, double distance, size_t *nbr_add,
                 double *nbr_dist)
{
  if (distance < dist_min)
    {
      const auto srch_add = jj * nx + ii;
      for (int n = 0; n < 4; ++n)
        {
          if (distance < nbr_dist[n])
            {
              for (int i = 3; i > n; --i) nbr_add[i] = nbr_add[i - 1];
              for (int i = 3; i > n; --i) nbr_dist[i] = nbr_dist[i - 1];
              searchResult = -1;
              nbr_add[n] = srch_add;
              nbr_dist[n] = distance;
              dist_min = nbr_dist[3];
              break;
            }
        }
    }
}

int
grid_search_square_reg_2d_NN(size_t nx, size_t ny, size_t *nbr_add, double *nbr_dist, double plat, double plon,
                             const Varray<double> &src_center_lat, const Varray<double> &src_center_lon)
{
  int searchResult = 0;

  const auto tgtCoslat = std::cos(plat);
  const auto tgtSinlat = std::sin(plat);
  const auto tgtCoslon = std::cos(plon);
  const auto tgtSinlon = std::sin(plon);

  double dist_min = DBL_MAX;
  for (int n = 0; n < 4; ++n) nbr_dist[n] = DBL_MAX;

  size_t jjf = 0, jjl = ny - 1;
  if (plon >= src_center_lon[0] && plon <= src_center_lon[nx - 1])
    {
      if (src_center_lat[0] < src_center_lat[ny - 1])
        {
          if (plat <= src_center_lat[0])
            jjl = (ny == 1) ? 0 : 1;
          else
            jjf = (ny == 1) ? 0 : ny - 2;
        }
      else
        {
          if (plat >= src_center_lat[0])
            jjl = (ny == 1) ? 0 : 1;
          else
            jjf = (ny == 1) ? 0 : ny - 2;
        }
    }

  std::vector<double> sincoslon(nx);

  for (size_t ii = 0; ii < nx; ++ii)
    sincoslon[ii] = tgtCoslon * std::cos(src_center_lon[ii]) + tgtSinlon * std::sin(src_center_lon[ii]);

  for (size_t jj = jjf; jj <= jjl; ++jj)
    {
      const auto coslat = tgtCoslat * std::cos(src_center_lat[jj]);
      const auto sinlat = tgtSinlat * std::sin(src_center_lat[jj]);

      const auto jjskip = (jj > 1 && jj < (ny - 2));

      if (jjskip)
        {
          size_t ii = 0;
          auto distance = std::acos(coslat * sincoslon[ii] + sinlat);
          nbr_add_and_dist(searchResult, jj, ii, nx, dist_min, distance, nbr_add, nbr_dist);
          ii = nx - 1;
          distance = std::acos(coslat * sincoslon[ii] + sinlat);
          nbr_add_and_dist(searchResult, jj, ii, nx, dist_min, distance, nbr_add, nbr_dist);
        }
      else
        {
          for (size_t ii = 0; ii < nx; ++ii)
            {
              const auto distance = std::acos(coslat * sincoslon[ii] + sinlat);
              nbr_add_and_dist(searchResult, jj, ii, nx, dist_min, distance, nbr_add, nbr_dist);
            }
        }
    }

  for (int n = 0; n < 4; ++n) nbr_dist[n] = 1.0 / (nbr_dist[n] + TINY);
  double distance = 0.0;
  for (int n = 0; n < 4; ++n) distance += nbr_dist[n];
  for (int n = 0; n < 4; ++n) nbr_dist[n] /= distance;

  return searchResult;
}

int
grid_search_square_reg_2d(RemapGrid *srcGrid, size_t (&src_add)[4], double (&srcLats)[4], double (&srcLons)[4], double plat,
                          double plon)
{
  /*
    Input variables:

      plat : latitude  of the search point
      plon : longitude of the search point

    Output variables:

      src_add[4]  : address of each corner point enclosing P
      srcLats[4] : latitudes  of the four corner points
      srcLons[4] : longitudes of the four corner points
  */
  int searchResult = 0;
  const auto &src_center_lat = srcGrid->reg2d_center_lat;
  const auto &src_center_lon = srcGrid->reg2d_center_lon;

  for (int n = 0; n < 4; ++n) src_add[n] = 0;

  const auto nx = srcGrid->dims[0];
  const auto ny = srcGrid->dims[1];

  const auto nxm = srcGrid->isCyclic ? nx + 1 : nx;

  if (plon < src_center_lon[0]) plon += PI2;
  if (plon > src_center_lon[nxm - 1]) plon -= PI2;

  size_t ii, jj;
  int lfound = rect_grid_search(ii, jj, plon, plat, nxm, ny, src_center_lon, src_center_lat);
  if (lfound)
    {
      auto iix = ii;
      if (srcGrid->isCyclic && iix == (nxm - 1)) iix = 0;
      src_add[0] = (jj - 1) * nx + (ii - 1);
      src_add[1] = (jj - 1) * nx + (iix);
      src_add[2] = (jj) *nx + (iix);
      src_add[3] = (jj) *nx + (ii - 1);

      srcLons[0] = src_center_lon[ii - 1];
      srcLons[1] = src_center_lon[iix];
      // For consistency, we must make sure all lons are in same 2pi interval
      if (srcLons[0] > PI2) srcLons[0] -= PI2;
      if (srcLons[0] < 0) srcLons[0] += PI2;
      if (srcLons[1] > PI2) srcLons[1] -= PI2;
      if (srcLons[1] < 0) srcLons[1] += PI2;
      srcLons[2] = srcLons[1];
      srcLons[3] = srcLons[0];

      srcLats[0] = src_center_lat[jj - 1];
      srcLats[1] = srcLats[0];
      srcLats[2] = src_center_lat[jj];
      srcLats[3] = srcLats[2];

      searchResult = 1;

      return searchResult;
    }

  /*
    If no cell found, point is likely either in a box that straddles either pole
    or is outside the grid. Fall back to a distance-weighted average of the four
    closest points. Go ahead and compute weights here, but store in srcLats and
    return -add to prevent the parent routine from computing bilinear weights.
  */
  if (!srcGrid->doExtrapolate) return searchResult;

  searchResult = grid_search_square_reg_2d_NN(nx, ny, src_add, srcLats, plat, plon, src_center_lat, src_center_lon);

  return searchResult;
}
