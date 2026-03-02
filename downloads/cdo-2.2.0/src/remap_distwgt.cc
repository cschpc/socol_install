/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <atomic>

#include "process_int.h"
#include "cdo_wtime.h"
#include "remap.h"
#include "remap_store_link.h"
#include "cdo_options.h"
#include "progress.h"
#include "cimdOmp.h"

// Interpolation using a distance-weighted average

// This routine computes the inverse-distance weights for a nearest-neighbor interpolation
void
remap_distwgt_weights(size_t numNeighbors, RemapSearch &rsearch, RemapVars &rv)
{
  auto srcGrid = rsearch.srcGrid;
  auto tgtGrid = rsearch.tgtGrid;

  if (Options::cdoVerbose) cdo_print("Called %s()", __func__);

  progress::init();

  // Compute mappings from source to target grid

  auto tgtGridSize = tgtGrid->size;

  std::vector<WeightLinks> weightLinks(tgtGridSize);
  weight_links_alloc(numNeighbors, tgtGridSize, weightLinks);

  std::vector<knnWeightsType> knnWeights;
  knnWeights.reserve(Threading::ompNumThreads);
  for (int i = 0; i < Threading::ompNumThreads; ++i) knnWeights.push_back(knnWeightsType(numNeighbors));

  auto start = Options::cdoVerbose ? cdo_get_wtime() : 0.0;

  // Loop over target grid

  std::atomic<size_t> atomicCount{ 0 };

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
  for (size_t tgtCellIndex = 0; tgtCellIndex < tgtGridSize; ++tgtCellIndex)
    {
      auto ompthID = cdo_omp_get_thread_num();

      atomicCount++;
      if (ompthID == 0) progress::update(0, 1, (double) atomicCount / tgtGridSize);

      weightLinks[tgtCellIndex].nlinks = 0;

      if (!tgtGrid->mask[tgtCellIndex]) continue;

      auto &knnWgt = knnWeights[ompthID];
      auto llpoint = remapgrid_get_lonlat(tgtGrid, tgtCellIndex);

      // Find nearest grid points on source grid and distances to each point
      remap_search_points(rsearch, llpoint, knnWgt);

      // Compute weights based on inverse distance if mask is false, eliminate those points
      auto nadds = knnWgt.computeWeights(srcGrid->mask);

      for (size_t i = 0; i < nadds; ++i)
        if (knnWgt.m_mask[i]) tgtGrid->cell_frac[tgtCellIndex] = 1.0;

      // Store the link
      store_weightlinks(0, nadds, knnWgt.m_addr.data(), knnWgt.m_dist.data(), tgtCellIndex, weightLinks);
    }

  progress::update(0, 1, 1);

  grid_point_search_delete(rsearch.gps);

  weight_links_to_remap_links(0, tgtGridSize, weightLinks, rv);

  if (numNeighbors == 1) rv.linksPerValue = numNeighbors;

  if (Options::cdoVerbose) cdo_print("Point search nearest: %.2f seconds", cdo_get_wtime() - start);
}  // remap_distwgt_weights

template <typename T>
static void
remap_dist_wgt(size_t numNeighbors, RemapSearch &rsearch, const Varray<T> &srcArray, Varray<T> &tgtArray, T missval, size_t nmiss)
{
  auto srcGrid = rsearch.srcGrid;
  auto tgtGrid = rsearch.tgtGrid;

  if (Options::cdoVerbose) cdo_print("Called %s()", __func__);

  progress::init();

  // Compute mappings from source to target grid

  auto tgtGridSize = tgtGrid->size;
  auto srcGridSize = srcGrid->size;

  Varray<short> srcGridMask;
  if (nmiss)
    {
      srcGridMask.resize(srcGridSize, 1);
      remap_set_mask(srcGridSize, srcArray, missval, srcGridMask);
    }

  std::vector<knnWeightsType> knnWeights;
  knnWeights.reserve(Threading::ompNumThreads);
  for (int i = 0; i < Threading::ompNumThreads; ++i) knnWeights.push_back(knnWeightsType(numNeighbors));

  auto start = Options::cdoVerbose ? cdo_get_wtime() : 0.0;

  // Loop over target grid

  std::atomic<size_t> atomicCount{ 0 };

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
  for (size_t tgtCellIndex = 0; tgtCellIndex < tgtGridSize; ++tgtCellIndex)
    {
      auto ompthID = cdo_omp_get_thread_num();

      atomicCount++;
      if (ompthID == 0) progress::update(0, 1, (double) atomicCount / tgtGridSize);

      auto &tgtValue = tgtArray[tgtCellIndex];
      tgtValue = missval;

      if (!tgtGrid->mask[tgtCellIndex]) continue;

      auto &knnWgt = knnWeights[ompthID];
      auto llpoint = remapgrid_get_lonlat(tgtGrid, tgtCellIndex);

      // Find nearest grid points on source grid and distances to each point
      remap_search_points(rsearch, llpoint, knnWgt);

      // Compute weights based on inverse distance if mask is false, eliminate those points
      auto nadds
          = (srcGridMask.size() > 0) ? knnWgt.computeWeights(srcGridMask) : knnWgt.computeWeights();
      if (nadds) tgtValue = knnWgt.arrayWeightsSum(srcArray);
    }

  progress::update(0, 1, 1);

  if (Options::cdoVerbose) cdo_print("Point search nearest: %.2f seconds", cdo_get_wtime() - start);
}  // remap_dist_wgt

void
remap_dist_wgt(size_t numNeighbors, RemapSearch &rsearch, const Field &field1, Field &field2)
{
  if (field1.memType != field2.memType) cdo_abort("Interal error, memType of field1 and field2 differ!");

  if (field1.memType == MemType::Float)
    remap_dist_wgt(numNeighbors, rsearch, field1.vec_f, field2.vec_f, (float) field1.missval, field1.nmiss);
  else
    remap_dist_wgt(numNeighbors, rsearch, field1.vec_d, field2.vec_d, field1.missval, field1.nmiss);
}

void
intgriddis(Field &field1, Field &field2, size_t numNeighbors)
{
  auto mapType = RemapMethod::DISTWGT;
  auto gridID1 = field1.grid;
  auto gridID2 = field2.grid;
  auto srcMissval = field1.missval;
  auto tgtMissval = field2.missval;
  const auto &srcArray = field1.vec_d;
  auto &tgtArray = field2.vec_d;

  if (Options::cdoVerbose) cdo_print("Called %s()", __func__);

  progress::init();

  // Interpolate from source to target grid

  RemapType remap;

  constexpr bool remap_extrapolate{ false };
  remap_set_int(REMAP_GENWEIGHTS, 0);
  remap_init_grids(mapType, remap_extrapolate, gridID1, remap.srcGrid, gridID2, remap.tgtGrid);

  auto srcGridSize = remap.srcGrid.size;
  auto tgtGridSize = remap.tgtGrid.size;

  Varray<short> srcGridMask;
  if (field1.nmiss)
    {
      srcGridMask.resize(srcGridSize, 1);
      remap_set_mask(srcGridSize, srcArray, srcMissval, srcGridMask);
    }

  std::vector<knnWeightsType> knnWeights;
  knnWeights.reserve(Threading::ompNumThreads);
  for (int i = 0; i < Threading::ompNumThreads; ++i) knnWeights.push_back(knnWeightsType(numNeighbors));

  remap_search_init(mapType, remap.search, remap.srcGrid, remap.tgtGrid);

  auto start = Options::cdoVerbose ? cdo_get_wtime() : 0.0;

  // Loop over target grid

  std::atomic<size_t> nmiss{ 0 };
  std::atomic<size_t> atomicCount{ 0 };

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
  for (size_t tgtCellIndex = 0; tgtCellIndex < tgtGridSize; ++tgtCellIndex)
    {
      auto ompthID = cdo_omp_get_thread_num();

      atomicCount++;
      if (ompthID == 0) progress::update(0, 1, (double) atomicCount / tgtGridSize);

      auto &tgtValue = tgtArray[tgtCellIndex];
      tgtValue = tgtMissval;

      // if (!tgt_mask[tgtCellIndex]) continue;

      auto &knnWgt = knnWeights[ompthID];
      auto llpoint = remapgrid_get_lonlat(&remap.tgtGrid, tgtCellIndex);

      // Find nearest grid points on source grid and distances to each point
      remap_search_points(remap.search, llpoint, knnWgt);

      // Compute weights based on inverse distance if mask is false, eliminate those points
      auto nadds
          = (srcGridMask.size() > 0) ? knnWgt.computeWeights(srcGridMask) : knnWgt.computeWeights();
      if (nadds)
        tgtMissval = knnWgt.arrayWeightsSum(srcArray);
      else
        nmiss++;
    }

  progress::update(0, 1, 1);

  field2.nmiss = nmiss;

  remap_grid_free(remap.srcGrid);
  remap_grid_free(remap.tgtGrid);
  remap_search_free(remap.search);

  if (Options::cdoVerbose) cdo_print("Point search nearest: %.2f seconds", cdo_get_wtime() - start);
}  // intgriddis

void
intgridnn(Field &field1, Field &field2)
{
  intgriddis(field1, field2, 1);
}
