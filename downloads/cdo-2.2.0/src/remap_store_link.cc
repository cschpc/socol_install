/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <algorithm>
#include <array>

#include "dmemory.h"
#include "cdo_options.h"
#include "remap.h"
#include "remap_store_link.h"

bool
is_sorted_list(size_t n, const size_t *list)
{
  for (size_t i = 1; i < n; ++i)
    if (list[i] < list[i - 1]) return false;
  return true;
}

static int
qcompare_index(const void *a, const void *b)
{
  return (((const IndexWeight *) a)->index < ((const IndexWeight *) b)->index);
}

static int
qcompare_index4(const void *a, const void *b)
{
  return (((const IndexWeight4 *) a)->index < ((const IndexWeight4 *) b)->index);
}

static void
sort_indexWeights(size_t numWeights, IndexWeight *indexWeights)
{
  std::qsort(indexWeights, numWeights, sizeof(IndexWeight), qcompare_index);
}

static void
sort_indexWeights4(IndexWeight4 *indexWeights)
{
  std::qsort(indexWeights, 4, sizeof(IndexWeight4), qcompare_index4);
}

void
store_weightlinks(int doAlloc, size_t numWeights, size_t *indices, double *weights, size_t cellIndex,
                  std::vector<WeightLinks> &weightLinks)
{
  weightLinks[cellIndex].nlinks = 0;
  weightLinks[cellIndex].offset = 0;

  if (numWeights)
    {
      auto indexWeights = doAlloc ? (IndexWeight *) Malloc(numWeights * sizeof(IndexWeight)) : weightLinks[cellIndex].indexWeights;

      for (size_t i = 0; i < numWeights; ++i)
        {
          indexWeights[i].index = indices[i];
          indexWeights[i].weight = weights[i];
        }

      if (numWeights > 1 && !is_sorted_list(numWeights, indices)) sort_indexWeights(numWeights, indexWeights);

      weightLinks[cellIndex].nlinks = numWeights;

      if (doAlloc) weightLinks[cellIndex].indexWeights = indexWeights;
    }
}

void
store_weightlinks_bicubic(size_t *indices, double (&weights)[4][4], size_t cellIndex, std::vector<WeightLinks4> &weightLinks)
{
  weightLinks[cellIndex].nlinks = 0;
  weightLinks[cellIndex].offset = 0;

  auto indexWeights = weightLinks[cellIndex].indexWeights;

  for (int i = 0; i < 4; ++i)
    {
      indexWeights[i].index = indices[i];
      for (int k = 0; k < 4; ++k) indexWeights[i].weight[k] = weights[i][k];
    }

  if (!is_sorted_list(4, indices)) sort_indexWeights4(indexWeights);

  weightLinks[cellIndex].nlinks = 4;
}

void
weight_links_to_remap_links(int doAlloc, size_t gridSize, std::vector<WeightLinks> &weightLinks, RemapVars &rv)
{
  size_t nlinks = 0;
  for (size_t i = 0; i < gridSize; ++i)
    {
      if (weightLinks[i].nlinks)
        {
          weightLinks[i].offset = nlinks;
          nlinks += weightLinks[i].nlinks;
        }
    }

  rv.maxLinks = nlinks;
  rv.numLinks = nlinks;

  if (nlinks)
    {
      auto num_wts = rv.num_wts;
      rv.srcCellIndices.resize(nlinks);
      rv.tgtCellIndices.resize(nlinks);
      rv.wts.resize(nlinks * num_wts);
      auto &srcCellIndices = rv.srcCellIndices;
      auto &tgtCellIndices = rv.tgtCellIndices;
      auto &wts = rv.wts;

#ifdef _OPENMP
#pragma omp parallel for schedule(static) default(shared)
#endif
      for (size_t i = 0; i < gridSize; ++i)
        {
          const auto numLinks = weightLinks[i].nlinks;
          if (numLinks)
            {
              const auto offset = weightLinks[i].offset;
              IndexWeight *indexWeights = weightLinks[i].indexWeights;
              for (size_t ilink = 0; ilink < numLinks; ++ilink)
                {
                  srcCellIndices[offset + ilink] = indexWeights[ilink].index;
                  tgtCellIndices[offset + ilink] = i;
                  wts[(offset + ilink) * num_wts] = indexWeights[ilink].weight;
                }
            }
        }

      if (doAlloc)
        {
          for (size_t i = 0; i < gridSize; ++i)
            {
              const auto numLinks = weightLinks[i].nlinks;
              if (numLinks) Free(weightLinks[i].indexWeights);
            }
        }
      else
        {
          Free(weightLinks[0].indexWeights);
        }
    }
}

void
weight_links_4_to_remap_links(size_t gridSize, std::vector<WeightLinks4> &weightLinks, RemapVars &rv)
{
  size_t nlinks = 0;
  for (size_t i = 0; i < gridSize; ++i)
    {
      if (weightLinks[i].nlinks)
        {
          weightLinks[i].offset = nlinks;
          nlinks += weightLinks[i].nlinks;
        }
    }

  rv.maxLinks = nlinks;
  rv.numLinks = nlinks;
  if (nlinks)
    {
      rv.srcCellIndices.resize(nlinks);
      rv.tgtCellIndices.resize(nlinks);
      rv.wts.resize(4 * nlinks);
      auto &srcCellIndices = rv.srcCellIndices;
      auto &tgtCellIndices = rv.tgtCellIndices;
      auto &wts = rv.wts;

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
      for (size_t i = 0; i < gridSize; ++i)
        {
          const auto numLinks = weightLinks[i].nlinks;
          if (numLinks)
            {
              const auto offset = weightLinks[i].offset;
              const auto indexWeights = weightLinks[i].indexWeights;
              for (size_t ilink = 0; ilink < numLinks; ++ilink)
                {
                  srcCellIndices[offset + ilink] = indexWeights[ilink].index;
                  tgtCellIndices[offset + ilink] = i;
                  for (size_t k = 0; k < 4; ++k) wts[(offset + ilink) * 4 + k] = indexWeights[ilink].weight[k];
                }
            }
        }

      Free(weightLinks[0].indexWeights);
    }
}

void
weight_links_alloc(size_t numNeighbors, size_t gridSize, std::vector<WeightLinks> &weightLinks)
{
  weightLinks[0].indexWeights = (IndexWeight *) Malloc(numNeighbors * gridSize * sizeof(IndexWeight));
  for (size_t i = 1; i < gridSize; ++i) weightLinks[i].indexWeights = weightLinks[0].indexWeights + numNeighbors * i;
}

void
weight_links_4_alloc(size_t gridSize, std::vector<WeightLinks4> &weightLinks)
{
  weightLinks[0].indexWeights = (IndexWeight4 *) Malloc(4 * gridSize * sizeof(IndexWeight4));
  for (size_t i = 1; i < gridSize; ++i) weightLinks[i].indexWeights = weightLinks[0].indexWeights + 4 * i;
}
