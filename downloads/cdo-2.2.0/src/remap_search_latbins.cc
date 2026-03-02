/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <algorithm>

#include "cdo_output.h"
#include "cdo_options.h"
#include <mpim_grid.h>
#include "remap.h"

static void
calcBinAddr(GridSearchBins &searchBins)
{
  const size_t ncells = searchBins.ncells;
  const size_t nbins = searchBins.nbins;
  size_t *bin_addr = &searchBins.bin_addr[0];
  const float *bin_lats = &searchBins.bin_lats[0];
  const float *cell_bound_box = &searchBins.cell_bound_box[0];

  for (size_t n = 0; n < nbins; ++n)
    {
      bin_addr[n * 2] = ncells;
      bin_addr[n * 2 + 1] = 0;
    }

  for (size_t nele = 0; nele < ncells; ++nele)
    {
      const size_t nele4 = nele << 2;
      const float cell_bound_box_lat1 = cell_bound_box[nele4];
      const float cell_bound_box_lat2 = cell_bound_box[nele4 + 1];
      for (size_t n = 0; n < nbins; ++n)
        {
          const size_t n2 = n << 1;
          if (cell_bound_box_lat1 <= bin_lats[n2 + 1] && cell_bound_box_lat2 >= bin_lats[n2])
            {
              bin_addr[n2] = std::min(nele, bin_addr[n2]);
              bin_addr[n2 + 1] = std::max(nele, bin_addr[n2 + 1]);
            }
        }
    }
}

void
calc_lat_bins(GridSearchBins &searchBins)
{
  const size_t nbins = searchBins.nbins;
  const double dlat = PI / nbins;  // lat interval for search bins

  if (Options::cdoVerbose) cdo_print("Using %zu latitude bins to restrict search.", nbins);

  if (nbins > 0)
    {
      searchBins.bin_lats.resize(2 * nbins);
      for (size_t n = 0; n < nbins; ++n)
        {
          searchBins.bin_lats[n * 2] = (n) *dlat - PIH;
          searchBins.bin_lats[n * 2 + 1] = (n + 1) * dlat - PIH;
        }

      searchBins.bin_addr.resize(2 * nbins);
      calcBinAddr(searchBins);
    }
}

size_t
get_srch_cells(size_t tgtCellIndex, GridSearchBins &tgtBins, GridSearchBins &srcBins, float *tgt_cell_bound_box,
               Varray<size_t> &srch_add)
{
  const size_t nbins = srcBins.nbins;
  const size_t srcGridSize = srcBins.ncells;
  const size_t *bin_addr1 = &tgtBins.bin_addr[0];
  const size_t *bin_addr2 = &srcBins.bin_addr[0];

  // Restrict searches first using search bins

  size_t min_add = srcGridSize - 1;
  size_t max_add = 0;

  for (size_t n = 0; n < nbins; ++n)
    {
      size_t n2 = n << 1;
      if (tgtCellIndex >= bin_addr1[n2] && tgtCellIndex <= bin_addr1[n2 + 1])
        {
          if (bin_addr2[n2] < min_add) min_add = bin_addr2[n2];
          if (bin_addr2[n2 + 1] > max_add) max_add = bin_addr2[n2 + 1];
        }
    }

  // Further restrict searches using bounding boxes

  float bound_box_lat1 = tgt_cell_bound_box[0];
  float bound_box_lat2 = tgt_cell_bound_box[1];
  float bound_box_lon1 = tgt_cell_bound_box[2];
  float bound_box_lon2 = tgt_cell_bound_box[3];

  const float *src_cell_bound_box = &srcBins.cell_bound_box[0];
  size_t srcCellIndexM4;
  size_t numSearchCells = 0;
  for (size_t srcCellIndex = min_add; srcCellIndex <= max_add; ++srcCellIndex)
    {
      srcCellIndexM4 = srcCellIndex << 2;
      if ((src_cell_bound_box[srcCellIndexM4 + 2] <= bound_box_lon2) && (src_cell_bound_box[srcCellIndexM4 + 3] >= bound_box_lon1))
        {
          if ((src_cell_bound_box[srcCellIndexM4] <= bound_box_lat2) && (src_cell_bound_box[srcCellIndexM4 + 1] >= bound_box_lat1))
            {
              srch_add[numSearchCells] = srcCellIndex;
              numSearchCells++;
            }
        }
    }

  if (bound_box_lon1 < 0.0f || bound_box_lon2 > PI2_f)
    {
      if (bound_box_lon1 < 0.0f)
        {
          bound_box_lon1 += PI2_f;
          bound_box_lon2 += PI2_f;
        }
      else
        {
          bound_box_lon1 -= PI2_f;
          bound_box_lon2 -= PI2_f;
        }

      for (size_t srcCellIndex = min_add; srcCellIndex <= max_add; ++srcCellIndex)
        {
          srcCellIndexM4 = srcCellIndex << 2;
          if ((src_cell_bound_box[srcCellIndexM4 + 2] <= bound_box_lon2)
              && (src_cell_bound_box[srcCellIndexM4 + 3] >= bound_box_lon1))
            {
              if ((src_cell_bound_box[srcCellIndexM4] <= bound_box_lat2)
                  && (src_cell_bound_box[srcCellIndexM4 + 1] >= bound_box_lat1))
                {
                  size_t ii;
                  for (ii = 0; ii < numSearchCells; ++ii)
                    if (srch_add[ii] == srcCellIndex) break;

                  if (ii == numSearchCells)
                    {
                      srch_add[numSearchCells] = srcCellIndex;
                      numSearchCells++;
                    }
                }
            }
        }
    }

  return numSearchCells;
}

static int
gridSearchSquareCurv2dNNScrip(size_t min_add, size_t max_add, size_t *nbr_add, double *nbr_dist, double plat, double plon,
                              const double *src_center_lat, const double *src_center_lon)
{
  int searchResult = 0;

  const auto tgtCoslat = std::cos(plat);
  const auto tgtSinlat = std::sin(plat);
  const auto tgtCoslon = std::cos(plon);
  const auto tgtSinlon = std::sin(plon);

  double dist_min = DBL_MAX;
  for (int n = 0; n < 4; ++n) nbr_dist[n] = DBL_MAX;

  double distance;
  for (size_t srch_add = min_add; srch_add <= max_add; ++srch_add)
    {
      distance = std::acos(tgtCoslat * std::cos(src_center_lat[srch_add])
                               * (tgtCoslon * std::cos(src_center_lon[srch_add]) + tgtSinlon * std::sin(src_center_lon[srch_add]))
                           + tgtSinlat * std::sin(src_center_lat[srch_add]));

      if (distance < dist_min)
        {
          for (int n = 0; n < 4; ++n)
            {
              if (distance < nbr_dist[n])
                {
                  for (int i = 3; i > n; --i)
                    {
                      nbr_add[i] = nbr_add[i - 1];
                      nbr_dist[i] = nbr_dist[i - 1];
                    }
                  searchResult = -1;
                  nbr_add[n] = srch_add;
                  nbr_dist[n] = distance;
                  dist_min = nbr_dist[3];
                  break;
                }
            }
        }
    }

  for (int n = 0; n < 4; ++n) nbr_dist[n] = 1.0 / (nbr_dist[n] + TINY);
  distance = 0.0;
  for (int n = 0; n < 4; ++n) distance += nbr_dist[n];
  for (int n = 0; n < 4; ++n) nbr_dist[n] /= distance;

  return searchResult;
}

static int
quadCrossProducts(double plon, double plat, double *lons, double *lats)
{
  int n;
  // Vectors for cross-product check
  double vec1_lat, vec1_lon;
  double vec2_lat, vec2_lon;

  // clang-format off
  // For consistency, we must make sure all lons are in same 2pi interval
  vec1_lon = lons[0] - plon;
  if      (vec1_lon >  PI) lons[0] -= PI2;
  else if (vec1_lon < -PI) lons[0] += PI2;

  for (n = 1; n < 4; ++n)
    {
      vec1_lon = lons[n] - lons[0];
      if      (vec1_lon >  PI) lons[n] -= PI2;
      else if (vec1_lon < -PI) lons[n] += PI2;
    }

  constexpr double crossEps = 1.e-20;
  int scross[4], scrossLast = 0;
  // corner_loop
  for (n = 0; n < 4; ++n)
    {
      const int next_n = (n + 1) % 4;
      /*
        Here we take the cross product of the vector making up each box side
        with the vector formed by the vertex and search point.
        If all the cross products are positive, the point is contained in the box.
      */
      vec1_lat = lats[next_n] - lats[n];
      vec1_lon = lons[next_n] - lons[n];
      vec2_lat = plat - lats[n];
      vec2_lon = plon - lons[n];

      // Check for 0,2pi crossings
      if      (vec1_lon >  3.0 * PIH) vec1_lon -= PI2;
      else if (vec1_lon < -3.0 * PIH) vec1_lon += PI2;

      if      (vec2_lon >  3.0 * PIH) vec2_lon -= PI2;
      else if (vec2_lon < -3.0 * PIH) vec2_lon += PI2;

      const double crossProduct = vec1_lon * vec2_lat - vec2_lon * vec1_lat;

      // If cross product is less than ZERO, this cell doesn't work
      // 2008-10-16 Uwe Schulzweida: bug fix for cross_product eq zero
      // 2022-06-11 Uwe Schulzweida: replace zero by crossEps
      scross[n] = (crossProduct < -crossEps) ? -1 : (crossProduct > crossEps) ? 1 : 0;
      if (n == 0) scrossLast = scross[n];
      if ((scross[n] < 0 && scrossLast > 0) || (scross[n] > 0 && scrossLast < 0)) break;
      scrossLast = scross[n];
    }

  if (n >= 4)
    {
      n = 0;
      if      (scross[0] >= 0 && scross[1] >= 0 && scross[2] >= 0 && scross[3] >= 0) n = 4;
      else if (scross[0] <= 0 && scross[1] <= 0 && scross[2] <= 0 && scross[3] <= 0) n = 4;
    }
  // clang-format on

  return n;
}

bool
point_in_quad(bool isCyclic, size_t nx, size_t ny, size_t i, size_t j, size_t adds[4], double lons[4], double lats[4], double plon,
              double plat, const double *centerLon, const double *centerLat)
{
  bool searchResult = false;
  const size_t ip1 = (i < (nx - 1)) ? i + 1 : isCyclic ? 0 : i;
  const size_t jp1 = (j < (ny - 1)) ? j + 1 : j;

  if (i == ip1 || j == jp1) return searchResult;

  size_t idx[4];
  idx[0] = j * nx + i;
  idx[1] = j * nx + ip1;    // east
  idx[2] = jp1 * nx + ip1;  // north-east
  idx[3] = jp1 * nx + i;    // north

  for (int k = 0; k < 4; ++k) lons[k] = centerLon[idx[k]];
  for (int k = 0; k < 4; ++k) lats[k] = centerLat[idx[k]];

  const int n = quadCrossProducts(plon, plat, lons, lats);

  // If cross products all same sign, we found the location
  if (n >= 4)
    {
      for (int k = 0; k < 4; ++k) adds[k] = idx[k];
      searchResult = true;
    }

  return searchResult;
}

int
grid_search_square_curv_2d_scrip(RemapGrid *srcGrid, size_t (&src_add)[4], double (&srcLats)[4], double (&srcLons)[4],
                                 double plat, double plon, GridSearchBins &srcBins)
{
  /*
    Output variables:

    int    src_add[4]              ! address of each corner point enclosing P
    double srcLats[4]             ! latitudes  of the four corner points
    double srcLons[4]             ! longitudes of the four corner points

    Input variables:

    double plat                    ! latitude  of the search point
    double plon                    ! longitude of the search point

    int src_grid_dims[2]           ! size of each src grid dimension

    double src_center_lat[]        ! latitude  of each src grid center
    double src_center_lon[]        ! longitude of each src grid center

    float srcGridBoundBox[][4]  ! bound box for source grid

    int src_bin_addr[][2]          ! latitude bins for restricting
  */
  int searchResult = 0;

  const size_t nbins = srcBins.nbins;
  const size_t *src_bin_addr = &srcBins.bin_addr[0];
  const float *bin_lats = &srcBins.bin_lats[0];
  const float *srcGridBoundBox = &srcBins.cell_bound_box[0];

  const float rlat = plat;
  const float rlon = plon;

  // restrict search first using bins

  for (int n = 0; n < 4; ++n) src_add[n] = 0;

  // addresses for restricting search
  size_t min_add = srcGrid->size - 1;
  size_t max_add = 0;

  for (size_t n = 0; n < nbins; ++n)
    {
      const size_t n2 = n << 1;
      if (rlat >= bin_lats[n2] && rlat <= bin_lats[n2 + 1])
        {
          if (src_bin_addr[n2] < min_add) min_add = src_bin_addr[n2];
          if (src_bin_addr[n2 + 1] > max_add) max_add = src_bin_addr[n2 + 1];
        }
    }

  /* Now perform a more detailed search */

  const size_t nx = srcGrid->dims[0];
  const size_t ny = srcGrid->dims[1];

  for (size_t srch_add = min_add; srch_add <= max_add; ++srch_add)
    {
      const size_t srch_add4 = srch_add << 2;
      /* First check bounding box */
      if (rlon >= srcGridBoundBox[srch_add4 + 2] && rlon <= srcGridBoundBox[srch_add4 + 3] && rlat >= srcGridBoundBox[srch_add4]
          && rlat <= srcGridBoundBox[srch_add4 + 1])
        {
          /* We are within bounding box so get really serious */

          /* Determine neighbor addresses */
          size_t j = srch_add / nx;
          size_t i = srch_add - j * nx;

          if (point_in_quad(srcGrid->isCyclic, nx, ny, i, j, src_add, srcLons, srcLats, plon, plat,
                            srcGrid->cell_center_lon.data(), srcGrid->cell_center_lat.data()))
            {
              searchResult = 1;
              return searchResult;
            }

          /* Otherwise move on to next cell */

        } /* Bounding box check */
    }

  /*
    If no cell found, point is likely either in a box that straddles either pole or is outside the grid.
    Fall back to a distance-weighted average of the four closest points. Go ahead and compute weights here,
    but store in srcLats and return -add to prevent the parent routine from computing bilinear weights.
  */
  if (!srcGrid->doExtrapolate) return searchResult;

  /*
    printf("Could not find location for %g %g\n", plat*RAD2DEG, plon*RAD2DEG);
    printf("Using nearest-neighbor average for this point\n");
  */
  searchResult = gridSearchSquareCurv2dNNScrip(min_add, max_add, src_add, srcLats, plat, plon, srcGrid->cell_center_lat.data(),
                                                srcGrid->cell_center_lon.data());

  return searchResult;
}
