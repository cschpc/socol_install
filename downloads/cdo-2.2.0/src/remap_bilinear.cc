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
#include "grid_healpix.h"

// bilinear interpolation

static inline void
limit_dphi_bounds(double &dphi)
{
  if (dphi > 3.0 * PIH) dphi -= PI2;
  if (dphi < -3.0 * PIH) dphi += PI2;
}

std::pair<double, double>
remap_find_weights(const LonLatPoint &llpoint, const double (&srcLons)[4], const double (&srcLats)[4])
{
  constexpr double converge = 1.0e-10;  // Convergence criterion
  extern long remap_max_iter;

  // Iterate to find xfrac,yfrac for bilinear approximation

  // some latitude  differences
  auto dth1 = srcLats[1] - srcLats[0];
  auto dth2 = srcLats[3] - srcLats[0];
  auto dth3 = srcLats[2] - srcLats[1] - dth2;

  // some longitude differences
  auto dph1 = srcLons[1] - srcLons[0];
  auto dph2 = srcLons[3] - srcLons[0];
  auto dph3 = srcLons[2] - srcLons[1];

  limit_dphi_bounds(dph1);
  limit_dphi_bounds(dph2);
  limit_dphi_bounds(dph3);

  dph3 -= dph2;

  // current guess for bilinear coordinate
  double xguess = 0.5;
  double yguess = 0.5;

  long iter = 0;  // iteration counters
  for (iter = 0; iter < remap_max_iter; ++iter)
    {
      auto dthp = llpoint.lat - srcLats[0] - dth1 * xguess - dth2 * yguess - dth3 * xguess * yguess;
      auto dphp = llpoint.lon - srcLons[0];

      limit_dphi_bounds(dphp);

      dphp -= dph1 * xguess + dph2 * yguess + dph3 * xguess * yguess;

      auto mat1 = dth1 + dth3 * yguess;
      auto mat2 = dth2 + dth3 * xguess;
      auto mat3 = dph1 + dph3 * yguess;
      auto mat4 = dph2 + dph3 * xguess;

      auto determinant = mat1 * mat4 - mat2 * mat3;

      auto deli = (dthp * mat4 - dphp * mat2) / determinant;
      auto delj = (dphp * mat1 - dthp * mat3) / determinant;

      if (std::fabs(deli) < converge && std::fabs(delj) < converge) break;

      xguess += deli;
      yguess += delj;
    }

  if (iter >= remap_max_iter) xguess = yguess = -1.0;

  return std::make_pair(xguess, yguess);
}

static void
bilinear_set_weights(double xfrac, double yfrac, double (&weights)[4])
{
  // clang-format off
  weights[0] = (1.0 - xfrac) * (1.0 - yfrac);
  weights[1] =        xfrac  * (1.0 - yfrac);
  weights[2] =        xfrac  *        yfrac;
  weights[3] = (1.0 - xfrac) *        yfrac;
  // clang-format on
}

int
num_src_points(const Varray<short> &mask, const size_t (&indices)[4], double (&srcLats)[4])
{
  int icount = 4;

  for (int i = 0; i < 4; ++i)
    {
      if (mask[indices[i]] == 0)
        {
          icount--;
          srcLats[i] = 0.0;
        }
    }

  return icount;
}

static void
renormalize_weights(const double (&srcLats)[4], double (&weights)[4])
{
  // sum of weights for normalization
  auto sumWeights = std::fabs(srcLats[0]) + std::fabs(srcLats[1]) + std::fabs(srcLats[2]) + std::fabs(srcLats[3]);
  for (int i = 0; i < 4; ++i) weights[i] = std::fabs(srcLats[i]) / sumWeights;
}

static void
bilinear_sort_weights(size_t (&indices)[4], double (&weights)[4])
{
  constexpr size_t numWeights = 4;

  if (is_sorted_list(numWeights, indices)) return;

  struct IndexWeightX
  {
    size_t index;
    double weight;
  };

  std::array<IndexWeightX, numWeights> indexWeights;

  for (size_t i = 0; i < numWeights; ++i)
    {
      indexWeights[i].index = indices[i];
      indexWeights[i].weight = weights[i];
    }

  auto comp_index = [](const auto &a, const auto &b) noexcept { return a.index < b.index; };
  std::sort(indexWeights.begin(), indexWeights.end(), comp_index);

  for (size_t i = 0; i < numWeights; ++i)
    {
      indices[i] = indexWeights[i].index;
      weights[i] = indexWeights[i].weight;
    }
}

static void
bilinear_warning()
{
  static auto printWarning = true;
  if (Options::cdoVerbose || printWarning)
    {
      printWarning = false;
      cdo_warning("Bilinear interpolation failed for some grid points - used a distance-weighted average instead!");
    }
}

static void
remap_bilinear_weights_regular(RemapSearch &rsearch, const Varray<short> &srcGridMask, const LonLatPoint &llpoint,
                               double &tgtCellFrac, size_t tgtCellIndex, std::vector<WeightLinks> &weightLinks)
{
  double srcLats[4];  //  latitudes  of four bilinear corners
  double srcLons[4];  //  longitudes of four bilinear corners
  double weights[4];  //  bilinear weights for four corners
  size_t indices[4];  //  address for the four source points

  // Find nearest square of grid points on source grid
  auto searchResult = remap_search_square(rsearch, llpoint, indices, srcLats, srcLons);

  // Check to see if points are mask points
  if (searchResult > 0) searchResult = remap_check_mask_indices(indices, srcGridMask);

  // If point found, find local xfrac, yfrac coordinates for weights
  if (searchResult > 0)
    {
      tgtCellFrac = 1.0;

      auto [xfrac, yfrac] = remap_find_weights(llpoint, srcLons, srcLats);
      if (xfrac >= 0.0 && yfrac >= 0.0)
        {
          // Successfully found xfrac, yfrac - compute weights
          bilinear_set_weights(xfrac, yfrac, weights);
          store_weightlinks(0, 4, indices, weights, tgtCellIndex, weightLinks);
        }
      else
        {
          bilinear_warning();
          searchResult = -1;
        }
    }

  /*
    Search for bilinear failed - use a distance-weighted average instead
    (this is typically near the pole) Distance was stored in srcLats!
  */
  if (searchResult < 0)
    {
      if (num_src_points(srcGridMask, indices, srcLats) > 0)
        {
          tgtCellFrac = 1.0;
          renormalize_weights(srcLats, weights);
          store_weightlinks(0, 4, indices, weights, tgtCellIndex, weightLinks);
        }
    }
}

static void
remap_bilinear_weights_healpix(const RemapSearch &rsearch, const Varray<short> &srcGridMask, const LonLatPoint &llpoint,
                               double &tgtCellFrac, size_t tgtCellIndex, std::vector<WeightLinks> &weightLinks)
{
  auto nside = rsearch.srcGrid->nside;
  auto order = rsearch.srcGrid->order;

  double weights[4];  //  bilinear weights for four corners
  size_t indices[4];  //  address for the four source points

  hp_bilinear_interpolate_weights(llpoint.lon, llpoint.lat, indices, weights, nside, order);

  // Check to see if points are mask points
  int searchResult = remap_check_mask_indices(indices, srcGridMask);

  if (searchResult > 0)
    {
      tgtCellFrac = 1.0;
      bilinear_sort_weights(indices, weights);
      store_weightlinks(0, 4, indices, weights, tgtCellIndex, weightLinks);
    }
}

// This routine computes the weights for a bilinear interpolation.
void
remap_bilinear_weights(RemapSearch &rsearch, RemapVars &rv)
{
  auto srcGrid = rsearch.srcGrid;
  auto tgtGrid = rsearch.tgtGrid;

  auto isHealpixGrid = (srcGrid->type == RemapGridType::HealPix);

  if (Options::cdoVerbose) cdo_print("Called %s()", __func__);

  if (!isHealpixGrid && srcGrid->rank != 2)
    cdo_abort("Can't do bilinear interpolation if the source grid is not a regular 2D grid!");

  auto start = Options::cdoVerbose ? cdo_get_wtime() : 0.0;

  progress::init();

  // Compute mappings from source to target grid

  auto tgtGridSize = tgtGrid->size;

  std::vector<WeightLinks> weightLinks(tgtGridSize);
  weight_links_alloc(4, tgtGridSize, weightLinks);

  std::atomic<size_t> atomicCount{ 0 };

  // Loop over target grid

#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
  for (size_t tgtCellIndex = 0; tgtCellIndex < tgtGridSize; ++tgtCellIndex)
    {
      atomicCount++;
      if (cdo_omp_get_thread_num() == 0) progress::update(0, 1, (double) atomicCount / tgtGridSize);

      weightLinks[tgtCellIndex].nlinks = 0;

      if (!tgtGrid->mask[tgtCellIndex]) continue;

      auto llpoint = remapgrid_get_lonlat(tgtGrid, tgtCellIndex);

      auto &tgtCellFrac = tgtGrid->cell_frac[tgtCellIndex];
      if (isHealpixGrid)
        remap_bilinear_weights_healpix(rsearch, srcGrid->mask, llpoint, tgtCellFrac, tgtCellIndex, weightLinks);
      else
        remap_bilinear_weights_regular(rsearch, srcGrid->mask, llpoint, tgtCellFrac, tgtCellIndex, weightLinks);
    }

  progress::update(0, 1, 1);

  weight_links_to_remap_links(0, tgtGridSize, weightLinks, rv);

  if (Options::cdoVerbose) cdo_print("%s: %.2f seconds", __func__, cdo_get_wtime() - start);
}  // remap_bilinear_weights

template <typename T>
static inline T
bilinear_remap(const Varray<T> &srcArray, const double (&weights)[4], const size_t (&indices)[4])
{
  // *tgtPoint = 0.0;
  // for (int i = 0; i < 4; ++i) *tgtPoint += srcArray[indices[i]]*weights[i];
  return srcArray[indices[0]] * weights[0] + srcArray[indices[1]] * weights[1] + srcArray[indices[2]] * weights[2]
         + srcArray[indices[3]] * weights[3];
}

template <typename T>
void
remap_set_mask(size_t gridsize, const Varray<T> &v, T missval, Varray<short> &mask)
{
  if (std::isnan(missval))
    {
#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
      for (size_t i = 0; i < gridsize; ++i) mask[i] = !dbl_is_equal(v[i], missval);
    }
  else
    {
#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
      for (size_t i = 0; i < gridsize; ++i) mask[i] = !is_equal(v[i], missval);
    }
}

// Explicit instantiation
template void remap_set_mask(size_t gridsize, const Varray<float> &v, float missval, Varray<short> &mask);
template void remap_set_mask(size_t gridsize, const Varray<double> &v, double missval, Varray<short> &mask);

template <typename T>
static void
remap_bilinear_regular(RemapSearch &rsearch, const Varray<T> &srcArray, const Varray<short> &srcGridMask,
                       const LonLatPoint &llpoint, T &tgtValue)
{
  double srcLats[4];  //  latitudes  of four bilinear corners
  double srcLons[4];  //  longitudes of four bilinear corners
  double weights[4];  //  bilinear weights for four corners
  size_t indices[4];  //  address for the four source points

  // Find nearest square of grid points on source grid
  auto searchResult = remap_search_square(rsearch, llpoint, indices, srcLats, srcLons);

  // Check to see if points are mask points
  if (searchResult > 0) searchResult = remap_check_mask_indices(indices, srcGridMask);

  // If point found, find local xfrac, yfrac coordinates for weights
  if (searchResult > 0)
    {
      auto [xfrac, yfrac] = remap_find_weights(llpoint, srcLons, srcLats);
      if (xfrac >= 0.0 && yfrac >= 0.0)
        {
          // Successfully found xfrac, yfrac - compute weights
          bilinear_set_weights(xfrac, yfrac, weights);
          bilinear_sort_weights(indices, weights);
          tgtValue = bilinear_remap(srcArray, weights, indices);
        }
      else
        {
          bilinear_warning();
          searchResult = -1;
        }
    }

  /*
    Search for bilinear failed - use a distance-weighted average instead
    (this is typically near the pole) Distance was stored in srcLats!
  */
  if (searchResult < 0)
    {
      if (srcGridMask.size() == 0 || num_src_points(srcGridMask, indices, srcLats) > 0)
        {
          renormalize_weights(srcLats, weights);
          bilinear_sort_weights(indices, weights);
          tgtValue = bilinear_remap(srcArray, weights, indices);
        }
    }
}

template <typename T>
static void
remap_bilinear_healpix(const RemapSearch &rsearch, const Varray<T> &srcArray, const Varray<short> &srcGridMask,
                       const LonLatPoint &llpoint, T &tgtValue)
{
  auto nside = rsearch.srcGrid->nside;
  auto order = rsearch.srcGrid->order;

  double weights[4];  //  bilinear weights for four corners
  size_t indices[4];  //  address for the four source points

  hp_bilinear_interpolate_weights(llpoint.lon, llpoint.lat, indices, weights, nside, order);

  // Check to see if points are mask points
  int searchResult = remap_check_mask_indices(indices, srcGridMask);

  if (searchResult > 0)
    {
      bilinear_sort_weights(indices, weights);
      tgtValue = bilinear_remap(srcArray, weights, indices);
    }
}

// This routine computes and apply the weights for a bilinear interpolation.
template <typename T>
static void
remap_bilinear(RemapSearch &rsearch, const Varray<T> &srcArray, Varray<T> &tgtArray, T missval, size_t nmiss)
{
  auto srcGrid = rsearch.srcGrid;
  auto tgtGrid = rsearch.tgtGrid;

  auto isHealpixGrid = (srcGrid->type == RemapGridType::HealPix);

  if (Options::cdoVerbose) cdo_print("Called %s()", __func__);

  if (!isHealpixGrid && srcGrid->rank != 2)
    cdo_abort("Can't do bilinear interpolation if the source grid is not a regular 2D grid!");

  auto start = Options::cdoVerbose ? cdo_get_wtime() : 0.0;

  progress::init();

  auto tgtGridSize = tgtGrid->size;
  auto srcGridSize = srcGrid->size;

  Varray<short> srcGridMask;
  if (nmiss)
    {
      srcGridMask.resize(srcGridSize, 1);
      remap_set_mask(srcGridSize, srcArray, missval, srcGridMask);
    }

  // Compute mappings from source to target grid

  std::atomic<size_t> atomicCount{ 0 };

  // Loop over target grid

#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
  for (size_t tgtCellIndex = 0; tgtCellIndex < tgtGridSize; ++tgtCellIndex)
    {
      atomicCount++;
      if (cdo_omp_get_thread_num() == 0) progress::update(0, 1, (double) atomicCount / tgtGridSize);

      auto &tgtValue = tgtArray[tgtCellIndex];
      tgtValue = missval;

      if (!tgtGrid->mask[tgtCellIndex]) continue;

      auto llpoint = remapgrid_get_lonlat(tgtGrid, tgtCellIndex);

      if (isHealpixGrid)
        remap_bilinear_healpix(rsearch, srcArray, srcGridMask, llpoint, tgtValue);
      else
        remap_bilinear_regular(rsearch, srcArray, srcGridMask, llpoint, tgtValue);
    }

  progress::update(0, 1, 1);

  if (Options::cdoVerbose) cdo_print("%s: %.2f seconds", __func__, cdo_get_wtime() - start);
}  // remap_bilinear

void
remap_bilinear(RemapSearch &rsearch, const Field &field1, Field &field2)
{
  if (field1.memType != field2.memType) cdo_abort("Interal error, memType of field1 and field2 differ!");

  if (field1.memType == MemType::Float)
    remap_bilinear(rsearch, field1.vec_f, field2.vec_f, (float) field1.missval, field1.nmiss);
  else
    remap_bilinear(rsearch, field1.vec_d, field2.vec_d, field1.missval, field1.nmiss);
}
