/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include "cdo_math.h"
#include "remap.h"
#include <mpim_grid.h>

static size_t
fill_src_add(bool isCyclic, long nx, long ny, long ii, long jj, long k, size_t *psrc_add)
{
  k /= 2;

  auto j0 = jj - k;
  auto jn = jj + k;
  auto i0 = ii - k;
  auto in = ii + k;
  if (j0 < 0) j0 = 0;
  if (jn >= ny) jn = ny - 1;
  if ((in - i0) > nx)
    {
      i0 = 0;
      in = nx - 1;
    }

  size_t num_add = 0;

  for (long j = j0; j <= jn; ++j)
    for (long i = i0; i <= in; ++i)
      {
        auto ix = i;
        if (isCyclic && ix < 0) ix += nx;
        if (isCyclic && ix >= nx) ix -= nx;
        if (ix >= 0 && ix < nx && j < ny) psrc_add[num_add++] = j * nx + ix;
      }

  return num_add;
}

static void
store_distance_healpix(GridPointSearch &gps, double plon, double plat, knnWeightsType &knnWeights, size_t num_add, size_t *indices, double *srcLons, double *srcLats)
{
  double xyz[3], query_pt[3];
  gcLLtoXYZ(plon, plat, query_pt);
  auto sqrSearchRadius = cdo::sqr(gps.searchRadius);

  for (size_t na = 0; na < num_add; ++na)
    {
      auto nadd = indices[na];
      auto lon = srcLons[na];
      auto lat = srcLats[na];

      gcLLtoXYZ(lon, lat, xyz);

      // Find distance to this point
      double sqrDist = (float) squareDistance(query_pt, xyz);
      if (sqrDist <= sqrSearchRadius)
        {
          // Store the address and distance if this is one of the smallest so far
          knnWeights.storeDistance(nadd, std::sqrt(sqrDist));
        }
    }

  knnWeights.checkDistance();
}

static void
store_distance_reg2d(GridPointSearch &gps, double plon, double plat, knnWeightsType &knnWeights, size_t nx, size_t num_add,
                     size_t *psrc_add)
{
  const auto &coslon = gps.coslon;
  const auto &sinlon = gps.sinlon;
  const auto &coslat = gps.coslat;
  const auto &sinlat = gps.sinlat;

  double xyz[3], query_pt[3];
  gcLLtoXYZ(plon, plat, query_pt);
  auto sqrSearchRadius = cdo::sqr(gps.searchRadius);

  for (size_t na = 0; na < num_add; ++na)
    {
      auto nadd = psrc_add[na];
      auto iy = nadd / nx;
      auto ix = nadd - iy * nx;

      xyz[0] = coslat[iy] * coslon[ix];
      xyz[1] = coslat[iy] * sinlon[ix];
      xyz[2] = sinlat[iy];
      // Find distance to this point
      double sqrDist = (float) squareDistance(query_pt, xyz);
      if (sqrDist <= sqrSearchRadius)
        {
          // Store the address and distance if this is one of the smallest so far
          knnWeights.storeDistance(nadd, std::sqrt(sqrDist));
        }
    }

  knnWeights.checkDistance();
}

static void
gridSearchPointHealpix(GridPointSearch &gps, double plon, double plat, knnWeightsType &knnWeights)
{
  /*
    Input variables:

      plat : latitude  of the search point
      plon : longitude of the search point

    Output variables:

      knnWeights.m_addr[numNeighbors] :  address of each of the closest points
      knnWeights.m_dist[numNeighbors] : distance to each of the closest points
  */
  auto numNeighbors = knnWeights.maxNeighbors();

  // Initialize distance and address arrays
  knnWeights.initAddr();
  knnWeights.initDist();

  auto index = hp_lonlat_to_index(gps.order, gps.nside, plon, plat);

  int64_t neighbours[8];
  hp_get_neighbours(gps.order, gps.nside, index, neighbours);

  size_t indices[9];
  indices[0] = index;
  size_t numWeights = 1;
  for (int i = 0; i < 8; ++i)
    if (neighbours[i] >= 0) indices[numWeights++] = neighbours[i];

  double srcLons[9], srcLats[9];
  for (size_t i = 0; i < numWeights; ++i)
    hp_index_to_lonlat(gps.order, gps.nside, indices[i], &srcLons[i], &srcLats[i]);
    
  store_distance_healpix(gps, plon, plat, knnWeights, numWeights, indices, srcLons, srcLats);
}

// This routine finds the closest numNeighbor points to a search point and computes a distance to each of the neighbors
static void
gridSearchPointReg2d(GridPointSearch &gps, double plon, double plat, knnWeightsType &knnWeights)
{
  /*
    Input variables:

      plat : latitude  of the search point
      plon : longitude of the search point

    Output variables:

      knnWeights.m_addr[numNeighbors] :  address of each of the closest points
      knnWeights.m_dist[numNeighbors] : distance to each of the closest points
  */
  auto numNeighbors = knnWeights.maxNeighbors();
  auto &nbr_add = knnWeights.m_addr;
  auto &nbr_dist = knnWeights.m_dist;

  // Initialize distance and address arrays
  knnWeights.initAddr();
  knnWeights.initDist();

  const auto &src_center_lon = gps.reg2d_center_lon;
  const auto &src_center_lat = gps.reg2d_center_lat;

  long nx = gps.dims[0];
  long ny = gps.dims[1];
  size_t nxm = gps.isCyclic ? nx + 1 : nx;

  if (plon < src_center_lon[0]) plon += PI2;
  if (plon > src_center_lon[nxm - 1]) plon -= PI2;

  size_t ii, jj;
  auto lfound = rect_grid_search(ii, jj, plon, plat, nxm, ny, src_center_lon, src_center_lat);
  if (lfound)
    {
      if (gps.isCyclic && ii == (nxm - 1)) ii = 0;

      constexpr size_t MAX_SEARCH_CELLS = 25;
      size_t src_add[MAX_SEARCH_CELLS];
      size_t *psrc_add = src_add;

      size_t k;
      for (k = 3; k < 10000; k += 2)
        if (numNeighbors <= (size_t) (k - 2) * (k - 2)) break;

      std::vector<size_t> src_add_tmp;
      if ((k * k) > MAX_SEARCH_CELLS)
        {
          src_add_tmp.resize(k * k);
          psrc_add = src_add_tmp.data();
        }

      auto num_add = fill_src_add(gps.isCyclic, nx, ny, ii, jj, k, psrc_add);

      store_distance_reg2d(gps, plon, plat, knnWeights, nx, num_add, psrc_add);
    }
  else if (gps.extrapolate)
    {
      int searchResult = 0;

      if (numNeighbors < 4)
        {
          size_t nbr_add4[4];
          double nbr_dist4[4];
          for (size_t n = 0; n < numNeighbors; ++n) nbr_add4[n] = SIZE_MAX;
          searchResult = grid_search_square_reg_2d_NN(nx, ny, nbr_add4, nbr_dist4, plat, plon, src_center_lat, src_center_lon);
          if (searchResult < 0)
            {
              for (size_t n = 0; n < numNeighbors; ++n) nbr_add[n] = nbr_add4[n];
              for (size_t n = 0; n < numNeighbors; ++n) nbr_dist[n] = nbr_dist4[n];
            }
        }
      else
        {
          searchResult
              = grid_search_square_reg_2d_NN(nx, ny, nbr_add.data(), nbr_dist.data(), plat, plon, src_center_lat, src_center_lon);
        }

      if (searchResult >= 0)
        for (size_t n = 0; n < numNeighbors; ++n) nbr_add[n] = SIZE_MAX;
    }
}

void
grid_search_point(GridPointSearch &gps, double plon, double plat, knnWeightsType &knnWeights)
{
  /*
    Input variables:

      plat : latitude  of the search point
      plon : longitude of the search point

    Output variables:

      knnWeights.m_addr[numNeighbors] :  address of each of the closest points
      knnWeights.m_dist[numNeighbors] : distance to each of the closest points
  */
  auto numNeighbors = knnWeights.maxNeighbors();

  // check some more points if distance is the same use the smaller index (nadd)
  auto ndist = (numNeighbors > 8) ? numNeighbors + 8 : numNeighbors * 2;
  if (ndist > gps.n) ndist = gps.n;

  if (knnWeights.m_tmpaddr.empty()) knnWeights.m_tmpaddr.resize(ndist);
  if (knnWeights.m_tmpdist.empty()) knnWeights.m_tmpdist.resize(ndist);
  auto &adds = knnWeights.m_tmpaddr;
  auto &dist = knnWeights.m_tmpdist;

  size_t nadds = 0;
  if (numNeighbors == 1)
    nadds = grid_point_search_nearest(gps, plon, plat, adds.data(), dist.data());
  else
    nadds = grid_point_search_qnearest(gps, plon, plat, ndist, adds.data(), dist.data());

  ndist = nadds;
  if (ndist < numNeighbors) numNeighbors = ndist;

  // Initialize distance and address arrays
  knnWeights.initAddr();
  knnWeights.initDist();
  for (size_t i = 0; i < ndist; ++i) knnWeights.storeDistance(adds[i], dist[i], numNeighbors);

  knnWeights.checkDistance();
}

void
grid_search_point_smooth(GridPointSearch &gps, double plon, double plat, knnWeightsType &knnWeights)
{
  /*
    Input variables:

      plat : latitude  of the search point
      plon : longitude of the search point

    Output variables:

      knnWeights.m_addr[numNeighbors] :  address of each of the closest points
      knnWeights.m_dist[numNeighbors] : distance to each of the closest points
  */
  auto numNeighbors = knnWeights.maxNeighbors();
  auto checkDistance = (numNeighbors <= 32);

  // check some more points if distance is the same use the smaller index (nadd)
  auto ndist = checkDistance ? ((numNeighbors > 8) ? numNeighbors + 8 : numNeighbors * 2) : numNeighbors;
  if (ndist > gps.n) ndist = gps.n;

  if (knnWeights.m_tmpaddr.empty()) knnWeights.m_tmpaddr.resize(ndist);
  if (knnWeights.m_tmpdist.empty()) knnWeights.m_tmpdist.resize(ndist);
  auto &adds = knnWeights.m_tmpaddr;
  auto &dist = knnWeights.m_tmpdist;

  size_t nadds = 0;
  if (numNeighbors == 1)
    nadds = grid_point_search_nearest(gps, plon, plat, adds.data(), dist.data());
  else
    nadds = grid_point_search_qnearest(gps, plon, plat, ndist, adds.data(), dist.data());

  ndist = nadds;

  if (checkDistance)
    {
      if (ndist < numNeighbors) numNeighbors = ndist;

      // Initialize distance and address arrays
      knnWeights.initAddr(numNeighbors);
      knnWeights.initDist(numNeighbors);
      for (size_t i = 0; i < ndist; ++i) knnWeights.storeDistance(adds[i], dist[i], numNeighbors);
    }
  else
    {
      knnWeights.m_numNeighbors = ndist;
      for (size_t i = 0; i < ndist; ++i) knnWeights.m_addr[i] = adds[i];
      for (size_t i = 0; i < ndist; ++i) knnWeights.m_dist[i] = dist[i];
    }

  knnWeights.checkDistance();
}

void
remap_search_points(RemapSearch &rsearch, const LonLatPoint &llpoint, knnWeightsType &knnWeights)
{
  if (rsearch.srcGrid->type == RemapGridType::HealPix)
    gridSearchPointHealpix(rsearch.gps, llpoint.lon, llpoint.lat, knnWeights);
  else if (rsearch.srcGrid->type == RemapGridType::Reg2D)
    gridSearchPointReg2d(rsearch.gps, llpoint.lon, llpoint.lat, knnWeights);
  else
    grid_search_point(rsearch.gps, llpoint.lon, llpoint.lat, knnWeights);
}

static int
gridSearchSquareCurv2d(GridPointSearch &gps, RemapGrid *rgrid, size_t (&src_add)[4], double (&srcLats)[4], double (&srcLons)[4],
                       double plat, double plon)
{
  /*
    Input variables:

      plat : latitude  of the search point
      plon : longitude of the search point

    Output variables:

      src_add[4] :   address of each corner point enclosing P
      srcLats[4] :  latitudes  of the four corner points
      srcLons[4] :  longitudes of the four corner points
  */
  int searchResult = 0;

  for (int n = 0; n < 4; ++n) src_add[n] = 0;

  double dist = 0.0;
  size_t addr = 0;
  size_t nadds = grid_point_search_nearest(gps, plon, plat, &addr, &dist);
  if (nadds > 0)
    {
      auto nx = rgrid->dims[0];
      auto ny = rgrid->dims[1];

      for (int k = 0; k < 4; ++k)
        {
          // Determine neighbor addresses
          auto j = addr / nx;
          auto i = addr - j * nx;
          if (k == 0 || k == 2) i = (i > 0) ? i - 1 : rgrid->isCyclic ? nx - 1 : 0;
          if (k == 0 || k == 1) j = (j > 0) ? j - 1 : 0;
          if (point_in_quad(rgrid->isCyclic, nx, ny, i, j, src_add, srcLons, srcLats, plon, plat, rgrid->cell_center_lon.data(),
                            rgrid->cell_center_lat.data()))
            {
              searchResult = 1;
              return searchResult;
            }
        }
    }

  /*
    If no cell found, point is likely either in a box that straddles either pole or is outside the grid.
    Fall back to a distance-weighted average of the four closest points. Go ahead and compute weights here,
    but store in srcLats and return -add to prevent the parent routine from computing bilinear weights.
  */
  if (!rgrid->doExtrapolate) return searchResult;

  size_t ndist = 4;
  nadds = grid_point_search_qnearest(gps, plon, plat, ndist, src_add, srcLats);
  if (nadds == 4)
    {
      for (int n = 0; n < 4; ++n) srcLats[n] = 1.0 / (srcLats[n] + TINY);
      double distance = 0.0;
      for (int n = 0; n < 4; ++n) distance += srcLats[n];
      for (int n = 0; n < 4; ++n) srcLats[n] /= distance;
      searchResult = -1;
    }

  return searchResult;
}

int
remap_search_square(RemapSearch &rsearch, const LonLatPoint &llpoint, size_t (&src_add)[4], double (&srcLats)[4],
                    double (&srcLons)[4])
{
  if (rsearch.srcGrid->type == RemapGridType::Reg2D)
    return grid_search_square_reg_2d(rsearch.srcGrid, src_add, srcLats, srcLons, llpoint.lat, llpoint.lon);
  else if (rsearch.gps.in_use)
    return gridSearchSquareCurv2d(rsearch.gps, rsearch.srcGrid, src_add, srcLats, srcLons, llpoint.lat, llpoint.lon);
  else
    return grid_search_square_curv_2d_scrip(rsearch.srcGrid, src_add, srcLats, srcLons, llpoint.lat, llpoint.lon, rsearch.srcBins);
}
