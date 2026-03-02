/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <atomic>
#include <algorithm>  // std::sort

#include "process_int.h"
#include "cdo_wtime.h"
#include <mpim_grid.h>
#include "cdo_options.h"
#include "remap.h"
#include "remap_store_link.h"
#include "progress.h"
#include "cimdOmp.h"

// bicubic interpolation

static void
bicubic_set_weights(double xfrac, double yfrac, double (&weights)[4][4])
{
  auto xfrac1 = xfrac * xfrac * (xfrac - 1.0);
  auto xfrac2 = xfrac * (xfrac - 1.0) * (xfrac - 1.0);
  auto xfrac3 = xfrac * xfrac * (3.0 - 2.0 * xfrac);
  auto yfrac1 = yfrac * yfrac * (yfrac - 1.0);
  auto yfrac2 = yfrac * (yfrac - 1.0) * (yfrac - 1.0);
  auto yfrac3 = yfrac * yfrac * (3.0 - 2.0 * yfrac);
  // clang-format off
  weights[0][0] = (1.0 - yfrac3) * (1.0 - xfrac3);
  weights[1][0] = (1.0 - yfrac3) *        xfrac3;
  weights[2][0] =        yfrac3  *        xfrac3;
  weights[3][0] =        yfrac3  * (1.0 - xfrac3);
  weights[0][1] = (1.0 - yfrac3) *        xfrac2;
  weights[1][1] = (1.0 - yfrac3) *        xfrac1;
  weights[2][1] =        yfrac3  *        xfrac1;
  weights[3][1] =        yfrac3  *        xfrac2;
  weights[0][2] =        yfrac2  * (1.0 - xfrac3);
  weights[1][2] =        yfrac2  *        xfrac3;
  weights[2][2] =        yfrac1  *        xfrac3;
  weights[3][2] =        yfrac1  * (1.0 - xfrac3);
  weights[0][3] =        yfrac2  *        xfrac2;
  weights[1][3] =        yfrac2  *        xfrac1;
  weights[2][3] =        yfrac1  *        xfrac1;
  weights[3][3] =        yfrac1  *        xfrac2;
  // clang-format on
}

int num_src_points(const Varray<short> &mask, const size_t (&indices)[4], double (&srcLats)[4]);

static void
renormalize_weights(const double (&srcLats)[4], double (&weights)[4][4])
{
  // sum of weights for normalization
  auto sumWeights = std::fabs(srcLats[0]) + std::fabs(srcLats[1]) + std::fabs(srcLats[2]) + std::fabs(srcLats[3]);
  for (int i = 0; i < 4; ++i) weights[i][0] = std::fabs(srcLats[i]) / sumWeights;
  for (int i = 0; i < 4; ++i) weights[i][1] = 0.0;
  for (int i = 0; i < 4; ++i) weights[i][2] = 0.0;
  for (int i = 0; i < 4; ++i) weights[i][3] = 0.0;
}

static void
bicubic_sort_weights(size_t (&indices)[4], double (&weights)[4][4])
{
  constexpr size_t numWeights = 4;

  if (is_sorted_list(numWeights, indices)) return;

  struct IndexWeightX
  {
    size_t index;
    double weight[4];
  };

  std::array<IndexWeightX, numWeights> indexWeights;

  for (size_t i = 0; i < numWeights; ++i)
    {
      indexWeights[i].index = indices[i];
      for (int k = 0; k < 4; ++k) indexWeights[i].weight[k] = weights[i][k];
    }

  auto comp_index = [](const auto &a, const auto &b) noexcept { return a.index < b.index; };
  std::sort(indexWeights.begin(), indexWeights.end(), comp_index);

  for (size_t i = 0; i < numWeights; ++i)
    {
      indices[i] = indexWeights[i].index;
      for (int k = 0; k < 4; ++k) weights[i][k] = indexWeights[i].weight[k];
    }
}

static void
bicubic_warning()
{
  static auto printWarning = true;

  if (Options::cdoVerbose || printWarning)
    {
      printWarning = false;
      cdo_warning("Bicubic interpolation failed for some grid points - used a distance-weighted average instead!");
    }
}

/*
  -----------------------------------------------------------------------

  This routine computes the weights for a bicubic interpolation.

  -----------------------------------------------------------------------
*/
void
remap_bicubic_weights(RemapSearch &rsearch, RemapVars &rv)
{
  auto srcGrid = rsearch.srcGrid;
  auto tgtGrid = rsearch.tgtGrid;

  if (Options::cdoVerbose) cdo_print("Called %s()", __func__);

  if (srcGrid->rank != 2) cdo_abort("Can't do bicubic interpolation if the source grid is not a regular 2D grid!");

  auto start = Options::cdoVerbose ? cdo_get_wtime() : 0.0;

  progress::init();

  // Compute mappings from source to target grid

  auto tgtGridSize = tgtGrid->size;

  std::vector<WeightLinks4> weightLinks(tgtGridSize);
  weight_links_4_alloc(tgtGridSize, weightLinks);

  std::atomic<size_t> atomicCount{ 0 };

  // Loop over target grid

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
  for (size_t tgtCellIndex = 0; tgtCellIndex < tgtGridSize; ++tgtCellIndex)
    {
      atomicCount++;
      if (cdo_omp_get_thread_num() == 0) progress::update(0, 1, (double) atomicCount / tgtGridSize);

      weightLinks[tgtCellIndex].nlinks = 0;

      if (!tgtGrid->mask[tgtCellIndex]) continue;

      auto llpoint = remapgrid_get_lonlat(tgtGrid, tgtCellIndex);

      double srcLats[4];     //  latitudes  of four bilinear corners
      double srcLons[4];     //  longitudes of four bilinear corners
      double weights[4][4];  //  bicubic weights for four corners
      size_t indices[4];     //  address for the four source points

      // Find nearest square of grid points on source grid
      auto searchResult = remap_search_square(rsearch, llpoint, indices, srcLats, srcLons);

      // Check to see if points are mask points
      if (searchResult > 0) searchResult = remap_check_mask_indices(indices, srcGrid->mask);

      // If point found, find local xfrac, yfrac coordinates for weights
      if (searchResult > 0)
        {
          tgtGrid->cell_frac[tgtCellIndex] = 1.0;

          auto [xfrac, yfrac] = remap_find_weights(llpoint, srcLons, srcLats);
          if (xfrac >= 0.0 && yfrac >= 0.0)
            {
              // Successfully found xfrac, yfrac - compute weights
              bicubic_set_weights(xfrac, yfrac, weights);
              store_weightlinks_bicubic(indices, weights, tgtCellIndex, weightLinks);
            }
          else
            {
              bicubic_warning();
              searchResult = -1;
            }
        }

      // Search for bicubic failed - use a distance-weighted average instead
      // (this is typically near the pole) Distance was stored in srcLats!
      if (searchResult < 0)
        {
          if (num_src_points(srcGrid->mask, indices, srcLats) > 0)
            {
              tgtGrid->cell_frac[tgtCellIndex] = 1.0;
              renormalize_weights(srcLats, weights);
              store_weightlinks_bicubic(indices, weights, tgtCellIndex, weightLinks);
            }
        }
    }

  progress::update(0, 1, 1);

  weight_links_4_to_remap_links(tgtGridSize, weightLinks, rv);

  if (Options::cdoVerbose) cdo_print("%s: %.2f seconds", __func__, cdo_get_wtime() - start);
}  // remap_bicubic_weights

/*
  -----------------------------------------------------------------------

  This routine computes and apply the weights for a bicubic interpolation.

  -----------------------------------------------------------------------
*/

template <typename T>
static T
bicubic_remap(const Varray<T> &srcArray, const double (&weights)[4][4], const size_t (&indices)[4], const RemapGradients &gradients)
{
  const auto &glat = gradients.grad_lat;
  const auto &glon = gradients.grad_lon;
  const auto &glatlon = gradients.grad_latlon;

  double tgtPoint = 0.0;
  for (int i = 0; i < 4; ++i)
    tgtPoint += srcArray[indices[i]] * weights[i][0] + glat[indices[i]] * weights[i][1] + glon[indices[i]] * weights[i][2]
                + glatlon[indices[i]] * weights[i][3];

  return tgtPoint;
}

template <typename T>
static void
remap_bicubic(RemapSearch &rsearch, const Varray<T> &srcArray, Varray<T> &tgtArray, T missval, size_t nmiss)
{
  auto srcGrid = rsearch.srcGrid;
  auto tgtGrid = rsearch.tgtGrid;

  if (Options::cdoVerbose) cdo_print("Called %s()", __func__);

  if (srcGrid->rank != 2) cdo_abort("Can't do bicubic interpolation if the source grid is not a regular 2D grid!");

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

  RemapGradients gradients(srcGrid->size);
  remap_gradients(*srcGrid, srcGridMask, srcArray, gradients);

  std::atomic<size_t> atomicCount{ 0 };

  // Loop over target grid

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
  for (size_t tgtCellIndex = 0; tgtCellIndex < tgtGridSize; ++tgtCellIndex)
    {
      atomicCount++;
      if (cdo_omp_get_thread_num() == 0) progress::update(0, 1, (double) atomicCount / tgtGridSize);

      tgtArray[tgtCellIndex] = missval;

      if (!tgtGrid->mask[tgtCellIndex]) continue;

      auto llpoint = remapgrid_get_lonlat(tgtGrid, tgtCellIndex);

      double srcLats[4];     //  latitudes  of four bilinear corners
      double srcLons[4];     //  longitudes of four bilinear corners
      double weights[4][4];  //  bicubic weights for four corners
      size_t indices[4];     //  address for the four source points

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
              bicubic_set_weights(xfrac, yfrac, weights);
              bicubic_sort_weights(indices, weights);
              tgtArray[tgtCellIndex] = bicubic_remap(srcArray, weights, indices, gradients);
            }
          else
            {
              bicubic_warning();
              searchResult = -1;
            }
        }

      // Search for bicubic failed - use a distance-weighted average instead
      // (this is typically near the pole) Distance was stored in srcLats!
      if (searchResult < 0)
        {
          if (srcGridMask.size() == 0 || num_src_points(srcGridMask, indices, srcLats) > 0)
            {
              renormalize_weights(srcLats, weights);
              bicubic_sort_weights(indices, weights);
              tgtArray[tgtCellIndex] = bicubic_remap(srcArray, weights, indices, gradients);
            }
        }
    }

  progress::update(0, 1, 1);

  if (Options::cdoVerbose) cdo_print("%s: %.2f seconds", __func__, cdo_get_wtime() - start);
}  // remap_bicubic

void
remap_bicubic(RemapSearch &rsearch, const Field &field1, Field &field2)
{
  if (field1.memType != field2.memType) cdo_abort("Interal error, memType of field1 and field2 differ!");

  if (field1.memType == MemType::Float)
    remap_bicubic(rsearch, field1.vec_f, field2.vec_f, (float) field1.missval, field1.nmiss);
  else
    remap_bicubic(rsearch, field1.vec_d, field2.vec_d, field1.missval, field1.nmiss);
}
