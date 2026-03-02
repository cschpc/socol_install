/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include "field.h"
#include "process_int.h"
#include "cdo_wtime.h"
#include "cdo_options.h"
#include "remap_vars.h"
#include "timer.h"
#include "cimdOmp.h"

/*
  -----------------------------------------------------------------------

  Performs the remapping based on weights computed elsewhere

  -----------------------------------------------------------------------
*/

template <typename T1, typename T2>
static void
remapFirstOrder(Varray<T2> &tgtArray, const RemapVars &rv, const Varray<T1> &srcArray)
{
  auto numLinks = rv.numLinks;
  auto num_wts = rv.num_wts;
  const auto &map_wts = rv.wts;
  const auto &tgt_add = rv.tgtCellIndices;
  const auto &src_add = rv.srcCellIndices;
  const auto &links = rv.links;
  auto linksPerValue = rv.linksPerValue;

  if (links.option)
    {
#ifdef HAVE_OPENMP4
#pragma omp parallel for simd default(shared) schedule(static)
#endif
      for (size_t n = 0; n < numLinks; ++n) tgtArray[tgt_add[n]] = static_cast<T2>(0.0);

      for (size_t j = 0; j < links.numBlocks; ++j)
        {
          const auto &tgt_addx = links.tgt_add[j];
          const auto &src_addx = links.src_add[j];
          const auto &windex = links.w_index[j];
          auto nlinks = links.numLinks[j];

#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
          for (size_t n = 0; n < nlinks; ++n) { tgtArray[tgt_addx[n]] += srcArray[src_addx[n]] * map_wts[num_wts * windex[n]]; }
        }
    }
  else
    {
      auto lpv = linksPerValue;
      if (lpv > 0)
        {
          size_t nlinks = numLinks / lpv;

          if (lpv == 1)
            {
#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
              for (size_t n = 0; n < nlinks; ++n) { tgtArray[tgt_add[n]] = srcArray[src_add[n]] * map_wts[n]; }
            }
          else if (lpv == 4)
            {
#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
              for (size_t n = 0; n < nlinks; ++n)
                {
                  auto noff = n * 4;
                  tgtArray[tgt_add[noff]] = srcArray[src_add[noff]] * map_wts[num_wts * (noff)]
                                            + srcArray[src_add[noff + 1]] * map_wts[num_wts * (noff + 1)]
                                            + srcArray[src_add[noff + 2]] * map_wts[num_wts * (noff + 2)]
                                            + srcArray[src_add[noff + 3]] * map_wts[num_wts * (noff + 3)];
                }
            }
          else
            {
#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
              for (size_t n = 0; n < nlinks; ++n)
                {
                  auto noff = n * lpv;
                  tgtArray[tgt_add[noff]] = srcArray[src_add[noff]] * map_wts[num_wts * noff];
                  for (size_t k = 1; k < (size_t) lpv; ++k)
                    tgtArray[tgt_add[noff]] += srcArray[src_add[noff + k]] * map_wts[num_wts * (noff + k)];
                }
            }
        }
      else
        {
#ifdef SX
#pragma cdir nodep
#endif
          for (size_t n = 0; n < numLinks; ++n) tgtArray[tgt_add[n]] = static_cast<T2>(0.0);

          for (size_t n = 0; n < numLinks; ++n)
            {
              // printf("%5zu %5zu %5zu %g # tgt_add src_add n\n", tgt_add[n], src_add[n], n, map_wts[num_wts*n]);
              tgtArray[tgt_add[n]] += srcArray[src_add[n]] * map_wts[num_wts * n];
            }
        }
    }
}

template <typename T1, typename T2>
static void
remapSecondOrder(Varray<T2> &tgtArray, const RemapVars &rv, const Varray<T1> &srcArray, RemapGradients &gradients)
{
  const auto &grad1 = gradients.grad_lat;
  const auto &grad2 = gradients.grad_lon;
  const auto &grad3 = gradients.grad_latlon;

  auto numLinks = rv.numLinks;
  auto num_wts = rv.num_wts;
  const auto &map_wts = rv.wts;
  const auto &tgt_add = rv.tgtCellIndices;
  const auto &src_add = rv.srcCellIndices;

#ifdef SX
#pragma cdir nodep
#endif
  for (size_t n = 0; n < numLinks; ++n) tgtArray[tgt_add[n]] = static_cast<T2>(0.0);

  if (num_wts == 3)
    {
      for (size_t n = 0; n < numLinks; ++n)
        {
          auto i = src_add[n];
          auto w = &map_wts[3 * n];
          tgtArray[tgt_add[n]] += srcArray[i] * w[0] + grad1[i] * w[1] + grad2[i] * w[2];
          // printf("%zu %zu %.5f %.5f %.5f %.5f %.5f\n", n, src_add[n], grad1[i], grad2[i], w[0], w[1], w[2]);
        }
    }
  else if (num_wts == 4)
    {
      for (size_t n = 0; n < numLinks; ++n)
        {
          auto i = src_add[n];
          auto w = &map_wts[4 * n];
          tgtArray[tgt_add[n]] += srcArray[i] * w[0] + grad1[i] * w[1] + grad2[i] * w[2] + grad3[i] * w[3];
        }
    }
}

template <typename T1, typename T2>
static void
remap(Varray<T2> &tgtArray, T2 missval, size_t tgtSize, const RemapVars &rv, const Varray<T1> &srcArray, RemapGradients &gradients)
{
  /*
    Input arrays:

      tgt_add    target address for each link
      src_add    source address for each link
      num_wts    num of weights used in remapping
      map_wts    remapping weights for each link
      srcArray  array with source field to be remapped

    Optional:

      gradients  gradient arrays on source grid necessary for higher-order remappings

    Output variables:

      tgtArray  array for remapped field on target grid
  */
  extern int timer_remap;

  auto start = Options::cdoVerbose ? cdo_get_wtime() : 0.0;

  // Check the order of the interpolation

  auto firstOrder = (gradients.grad_lat.size() == 0);

#ifdef HAVE_OPENMP4
#pragma omp parallel for simd default(shared) schedule(static)
#endif
  for (size_t n = 0; n < tgtSize; ++n) tgtArray[n] = missval;

  if (Options::Timer) timer_start(timer_remap);

  if (firstOrder)  // First order remapping
    {
      remapFirstOrder(tgtArray, rv, srcArray);
    }
  else  // Second order remapping
    {
      remapSecondOrder(tgtArray, rv, srcArray, gradients);
    }

  if (Options::cdoVerbose) cdo_print("Remap: %.2f seconds", cdo_get_wtime() - start);

  if (Options::Timer) timer_stop(timer_remap);
}

void
remap_field(Field &field2, double missval, size_t gridsize2, const RemapVars &rv, const Field &field1, RemapGradients &gradients)
{
  if (memtype_is_float_float(field1.memType, field2.memType))
    remap(field2.vec_f, (float) missval, gridsize2, rv, field1.vec_f, gradients);
  else if (memtype_is_float_double(field1.memType, field2.memType))
    remap(field2.vec_d, missval, gridsize2, rv, field1.vec_f, gradients);
  else if (memtype_is_double_float(field1.memType, field2.memType))
    remap(field2.vec_f, (float) missval, gridsize2, rv, field1.vec_d, gradients);
  else
    remap(field2.vec_d, missval, gridsize2, rv, field1.vec_d, gradients);
}

static size_t
get_max_add(size_t numLinks, size_t size, const Varray<size_t> &add)
{
  std::vector<size_t> isum(size, 0);

  for (size_t n = 0; n < numLinks; ++n) isum[add[n]]++;

  size_t max_add = 0;
  for (size_t i = 0; i < size; ++i)
    if (isum[i] > max_add) max_add = isum[i];

  return max_add;
}

static size_t
binary_search_int(const Varray<size_t> &array, size_t len, size_t value)
{
  int64_t low = 0, high = len - 1;

  while (low <= high)
    {
      auto midpoint = low + (high - low) / 2;

      // check to see if value is equal to item in array
      if (value == array[midpoint]) return midpoint;

      if (value < array[midpoint])
        high = midpoint - 1;
      else
        low = midpoint + 1;
    }

  // item was not found
  return len;
}

/*
  -----------------------------------------------------------------------

  Performs the remapping based on weights computed elsewhere

  -----------------------------------------------------------------------
*/
template <typename T1, typename T2>
static void
remap_laf(Varray<T2> &tgtArray, T2 missval, size_t tgtSize, const RemapVars &rv, const Varray<T1> &srcArray)
{
  /*
    Input:
      srcArray : array with source field to be remapped

    Output:
      tgtArray : array for remapped field on target grid
  */
  auto numLinks = rv.numLinks;
  auto num_wts = rv.num_wts;                // num of weights used in remapping
  const auto &map_wts = rv.wts;             // remapping weights for each link
  const auto &tgt_add = rv.tgtCellIndices;  // target address for each link
  const auto &src_add = rv.srcCellIndices;  // source address for each link

  varray_fill(tgtSize, tgtArray, missval);

  if (numLinks == 0) return;

  auto max_cls = get_max_add(numLinks, tgtSize, tgt_add);

#ifdef _OPENMP
  Varray2D<T1> src_cls2(Threading::ompNumThreads, Varray<T1>(max_cls));
  Varray2D<double> src_wts2(Threading::ompNumThreads, Varray<double>(max_cls));
#else
  Varray<T1> src_cls(max_cls);
  Varray<double> src_wts(max_cls);
#endif

  for (size_t n = 0; n < numLinks; ++n)
    if (dbl_is_equal(tgtArray[tgt_add[n]], missval)) tgtArray[tgt_add[n]] = 0.0;

#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(dynamic, 1)
#endif
  for (size_t i = 0; i < tgtSize; ++i)
    {
      size_t k, ncls;
#ifdef _OPENMP
      auto ompthID = cdo_omp_get_thread_num();
      auto &src_cls = src_cls2[ompthID];
      auto &src_wts = src_wts2[ompthID];
#endif
      varray_fill(src_cls, static_cast<T1>(0.0));
      varray_fill(src_wts, 0.0);
      /*
      ncls = 0;
      for ( n = 0; n < numLinks; n++ )
        {
          if ( i == tgt_add[n] )
            {
              for ( k = 0; k < ncls; k++ )
                if ( is_equal(srcArray[src_add[n]], src_cls[k]) ) break;

              if ( k == ncls )
                {
                  src_cls[k] = srcArray[src_add[n]];
                  ncls++;
                }

              src_wts[k] += map_wts[num_wts*n];
            }
        }
      */
      // only for sorted tgt_add!
      {
        size_t min_add = 1, max_add = 0;

        auto n = binary_search_int(tgt_add, numLinks, i);

        if (n < numLinks)
          {
            min_add = n;

            for (n = min_add + 1; n < numLinks; ++n)
              if (i != tgt_add[n]) break;

            max_add = n;

            for (n = min_add; n > 0; --n)
              if (i != tgt_add[n - 1]) break;

            min_add = n;
          }

        ncls = 0;
        for (n = min_add; n < max_add; ++n)
          {
            auto value = srcArray[src_add[n]];

            for (k = 0; k < ncls; ++k)
              if (is_equal(value, src_cls[k])) break;

            if (k == ncls)
              {
                src_cls[k] = value;
                ncls++;
              }

            src_wts[k] += map_wts[num_wts * n];
          }
        // printf("i, min_add, max_add, ncls %zu %zu %zu %zu\n", i, min_add, max_add, ncls);
      }

      if (ncls)
        {
          size_t imax = 0;
          auto wts = src_wts[0];
          for (k = 1; k < ncls; ++k)
            {
              if (src_wts[k] > wts)
                {
                  wts = src_wts[k];
                  imax = k;
                }
            }

          // for (k = 0; k < ncls; ++k) printf(" i  k, src_wts[k],  src_cls[k] %zu %zu %g %g\n", i, k, src_wts[k],  src_cls[k]);
          // printf("imax, src_wts[imax],  src_cls[imax] %zu %zu %g %g\n", i , imax, src_wts[imax],  src_cls[imax]);
          tgtArray[i] = src_cls[imax];
        }
    }
}

void
remap_laf(Field &field2, double missval, size_t gridsize2, const RemapVars &rv, const Field &field1)
{
  if (memtype_is_float_float(field1.memType, field2.memType))
    remap_laf(field2.vec_f, (float) missval, gridsize2, rv, field1.vec_f);
  else if (memtype_is_float_double(field1.memType, field2.memType))
    remap_laf(field2.vec_d, missval, gridsize2, rv, field1.vec_f);
  else if (memtype_is_double_float(field1.memType, field2.memType))
    remap_laf(field2.vec_f, (float) missval, gridsize2, rv, field1.vec_d);
  else
    remap_laf(field2.vec_d, missval, gridsize2, rv, field1.vec_d);
}

template <typename T1, typename T2>
static void
remap_avg(Varray<T2> &tgtArray, T2 missval, size_t tgtSize, const RemapVars &rv, const Varray<T1> &srcArray)
{
  auto numLinks = rv.numLinks;
  auto num_wts = rv.num_wts;                // num of weights used in remapping
  const auto &map_wts = rv.wts;             // remapping weights for each link
  const auto &tgt_add = rv.tgtCellIndices;  // target address for each link
  const auto &src_add = rv.srcCellIndices;  // source address for each link

  /*
  for (size_t n = 0; n < tgtSize; ++n) tgtArray[n] = missval;

  std::vector<int> count(tgtSize, 0);

#ifdef SX
#pragma cdir nodep
#endif
  for (size_t n = 0; n < numLinks; ++n)
    if (dbl_is_equal(tgtArray[tgt_add[n]], missval)) tgtArray[tgt_add[n]] = 0.0;

  for (size_t n = 0; n < numLinks; ++n)
    {
      // printf("%5d %5d %5d %g # tgt_add src_add n\n", tgt_add[n], src_add[n], n, map_wts[num_wts*n]);
      // tgtArray[tgt_add[n]] += srcArray[src_add[n]]*map_wts[num_wts*n];
      tgtArray[tgt_add[n]] += srcArray[src_add[n]];
      count[tgt_add[n]] += 1;
      if (src_cell_frac[src_add[n]] < 1.0)
        printf("%zu %zu %zu %g %g %g %g\n", n, tgt_add[n], src_add[n], srcArray[src_add[n]], map_wts[num_wts * n],
               tgtArray[tgt_add[n]], src_cell_frac[src_add[n]]);
    }

  for (size_t i = 0; i < tgtSize; ++i)
    {
      if (count[i] > 0) tgtArray[i] /= count[i];
    }
  */
  /*
#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(dynamic, 1)
#endif
  */
  // size_t max_values = 0;
  for (size_t i = 0; i < tgtSize; ++i)
    {
      size_t nvalues = 0;
      double sum = 0.0;
      // double wts = 0.0;
      // only for sorted tgt_add!
      {
        size_t min_add = 1, max_add = 0;

        auto n = binary_search_int(tgt_add, numLinks, i);

        if (n < numLinks)
          {
            min_add = n;

            for (n = min_add + 1; n < numLinks; ++n)
              if (i != tgt_add[n]) break;

            max_add = n;

            for (n = min_add; n > 0; --n)
              if (i != tgt_add[n - 1]) break;

            min_add = n;
          }

        // auto nadds = (max_add - min_add) + 1;
        // double lim = 0.1 / nadds;
        constexpr double lim = 0.0;
        for (n = min_add; n < max_add; ++n)
          {
            auto value = srcArray[src_add[n]];
            if (map_wts[num_wts * n] > lim && !dbl_is_equal(value, missval))
              {
                sum += value;
                // wts += map_wts[num_wts * n];
                nvalues++;
              }
          }
      }

      tgtArray[i] = (nvalues > 0) ? sum / nvalues : missval;
      // printf("%zu %zu %g %g\n", i+1, nvalues, wts, sum);
      // max_values += nvalues;
    }
  // printf("max_values = %zu  numLinks = %zu\n", max_values, numLinks);
}

void
remap_avg(Field &field2, double missval, size_t gridsize2, const RemapVars &rv, const Field &field1)
{
  if (memtype_is_float_float(field1.memType, field2.memType))
    remap_avg(field2.vec_f, (float) missval, gridsize2, rv, field1.vec_f);
  else if (memtype_is_float_double(field1.memType, field2.memType))
    remap_avg(field2.vec_d, missval, gridsize2, rv, field1.vec_f);
  else if (memtype_is_double_float(field1.memType, field2.memType))
    remap_avg(field2.vec_f, (float) missval, gridsize2, rv, field1.vec_d);
  else
    remap_avg(field2.vec_d, missval, gridsize2, rv, field1.vec_d);
}

void
remap_vars_init(RemapMethod mapType, int remapOrder, RemapVars &rv)
{
  // Initialize all pointer
  if (!rv.pinit) rv.pinit = true;

  rv.sort_add = (mapType == RemapMethod::CONSERV_SCRIP);

  // Determine the number of weights
  rv.num_wts = (mapType == RemapMethod::CONSERV_SCRIP) ? 3 : ((mapType == RemapMethod::BICUBIC) ? 4 : 1);
  if (mapType == RemapMethod::CONSERV && remapOrder == 2) rv.num_wts = 3;

  rv.linksPerValue = -1;
  rv.numLinks = 0;
  rv.maxLinks = 0;
  rv.resize_increment = 1024;

  rv.links.option = false;
  rv.links.maxLinks = 0;
  rv.links.numBlocks = 0;
}

void
remap_vars_ensure_size(RemapVars &rv, size_t size)
{
  if (size >= rv.maxLinks)
    {
      while (size >= rv.maxLinks) rv.maxLinks += rv.resize_increment;

      rv.srcCellIndices.resize(rv.maxLinks);
      rv.tgtCellIndices.resize(rv.maxLinks);
      rv.wts.resize(rv.num_wts * rv.maxLinks);
    }
}

void
remap_vars_resize(RemapVars &rv, size_t size)
{
  rv.maxLinks = size;

  rv.srcCellIndices.resize(rv.maxLinks);
  rv.tgtCellIndices.resize(rv.maxLinks);
  rv.wts.resize(rv.num_wts * rv.maxLinks);
}

void
remap_vars_reorder(RemapVars &rv)
{
  size_t nval = 0, numBlocks = 0;

  auto numLinks = rv.numLinks;

  printf("remap_vars_reorder\n");
  printf("  numLinks %zu\n", numLinks);
  rv.links.option = true;

  size_t lastval = -1;
  size_t maxLinks = 0;
  for (size_t n = 0; n < numLinks; ++n)
    {
      if (rv.tgtCellIndices[n] == lastval)
        nval++;
      else
        {
          if (nval > numBlocks) numBlocks = nval;
          nval = 1;
          maxLinks++;
          lastval = rv.tgtCellIndices[n];
        }
    }

  if (numBlocks)
    {
      rv.links.maxLinks = maxLinks;
      rv.links.numBlocks = numBlocks;

      printf("numLinks %zu  maxLinks %zu  numBlocks %zu\n", rv.numLinks, maxLinks, numBlocks);

      rv.links.numLinks.resize(numBlocks);
      rv.links.tgt_add.resize(numBlocks);
      rv.links.src_add.resize(numBlocks);
      rv.links.w_index.resize(numBlocks);
    }

  for (size_t j = 0; j < numBlocks; ++j)
    {
      rv.links.tgt_add[j].resize(maxLinks);
      rv.links.src_add[j].resize(maxLinks);
      rv.links.w_index[j].resize(maxLinks);
    }

  for (size_t j = 0; j < numBlocks; ++j)
    {
      nval = 0;
      lastval = -1;
      size_t nlinks = 0;

      for (size_t n = 0; n < numLinks; ++n)
        {
          if (rv.tgtCellIndices[n] == lastval)
            nval++;
          else
            {
              nval = 1;
              lastval = rv.tgtCellIndices[n];
            }

          if (nval == j + 1)
            {
              rv.links.tgt_add[j][nlinks] = rv.tgtCellIndices[n];
              rv.links.src_add[j][nlinks] = rv.srcCellIndices[n];
              rv.links.w_index[j][nlinks] = n;
              nlinks++;
            }
        }

      rv.links.numLinks[j] = nlinks;
      printf("loop %zu  nlinks %zu\n", j + 1, nlinks);
    }
}

void
remap_vars_free(RemapVars &rv)
{
  if (rv.pinit)
    {
      rv.pinit = false;
      rv.sort_add = false;

      varray_free(rv.srcCellIndices);
      varray_free(rv.tgtCellIndices);
      varray_free(rv.wts);

      if (rv.links.option)
        {
          rv.links.option = false;

          if (rv.links.numBlocks)
            {
              varray_free(rv.links.numLinks);
              auto numBlocks = rv.links.numBlocks;
              for (size_t i = 0; i < numBlocks; ++i)
                {
                  varray_free(rv.links.src_add[i]);
                  varray_free(rv.links.tgt_add[i]);
                  varray_free(rv.links.w_index[i]);
                }
              varray_free(rv.links.src_add);
              varray_free(rv.links.tgt_add);
              varray_free(rv.links.w_index);
            }
        }
    }
  else
    fprintf(stderr, "%s Warning: vars not initialized!\n", __func__);

}  // remap_vars_free

void
remap_vars_check_weights(const RemapVars &rv)
{
  auto numLinks = rv.numLinks;
  auto num_wts = rv.num_wts;
  auto normOpt = rv.normOpt;
  const auto &src_add = rv.srcCellIndices;
  const auto &tgt_add = rv.tgtCellIndices;
  const auto &wts = rv.wts;

  for (size_t n = 0; n < numLinks; ++n)
    {
      if (wts[n * num_wts] < -0.01)
        cdo_print("Map weight < 0! grid1idx=%zu grid2idx=%zu nlink=%zu wts=%g", src_add[n], tgt_add[n], n, wts[n * num_wts]);

      if (normOpt != NormOpt::NONE && wts[n * num_wts] > 1.01)
        cdo_print("Map weight > 1! grid1idx=%zu grid2idx=%zu nlink=%zu wts=%g", src_add[n], tgt_add[n], n, wts[n * num_wts]);
    }
}
