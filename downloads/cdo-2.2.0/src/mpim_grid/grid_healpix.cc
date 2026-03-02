#ifdef __cplusplus
extern "C"
{
#endif
#include "lib/healpix/healpix.h"
#include "lib/healpix/interpolation.h"
#ifdef __cplusplus
}
#endif

#include "grid_healpix.h"
#include <cstdio>
#include <cstring>
#include <climits>

HpOrder
hp_get_order(const std::string &orderName)
{
  // clang-format off
  if      (orderName == "xy")     return HpOrder::XY;
  else if (orderName == "ring")   return HpOrder::Ring;
  else if (orderName == "nest")   return HpOrder::Nested;
  else if (orderName == "nested") return HpOrder::Nested;
  // clang-format on
  return HpOrder::Undef;
}

static int64_t
hp_xy_to_xy(int64_t index, int nside)
{
  (void) nside;
  return index;
}

int64_t
hp_lonlat_to_index(HpOrder order, int nside, double xval, double yval)
{
  int64_t (*xy_to_order)(int64_t index, int nside);
  // clang-format off
  xy_to_order = (order == HpOrder::Ring)   ? &healpixl_xy_to_ring :
                (order == HpOrder::Nested) ? &healpixl_xy_to_nested : &hp_xy_to_xy;
  // clang-format on

  return xy_to_order(radec_to_healpixl(xval, yval, nside), nside);
}

void
hp_index_to_lonlat(HpOrder order, int nside, int64_t index, double *xval, double *yval)
{
  int64_t (*order_to_xy)(int64_t index, int nside);
  // clang-format off
  order_to_xy = (order == HpOrder::Ring)   ? &healpixl_ring_to_xy :
                (order == HpOrder::Nested) ? &healpixl_nested_to_xy : &hp_xy_to_xy;
  // clang-format on

  healpixl_to_radec(order_to_xy(index, nside), nside, 0.5, 0.5, xval, yval);
}

void
hp_get_neighbours(HpOrder order, int nside, int64_t index, int64_t *neighbours)
{
  int64_t (*order_to_xy)(int64_t index, int nside);
  // clang-format off
  order_to_xy = (order == HpOrder::Ring)   ? &healpixl_ring_to_xy :
                (order == HpOrder::Nested) ? &healpixl_nested_to_xy : &hp_xy_to_xy;
  // clang-format on

  healpixl_get_neighbours(order_to_xy(index, nside), neighbours, nside);

  int64_t (*xy_to_order)(int64_t index, int nside);
  // clang-format off
  xy_to_order = (order == HpOrder::Ring)   ? &healpixl_xy_to_ring :
                (order == HpOrder::Nested) ? &healpixl_xy_to_nested : &hp_xy_to_xy;
  // clang-format on

  for (int i = 0; i < 8; ++i)
    if (neighbours[i] >= 0) neighbours[i] = xy_to_order(neighbours[i], nside);
}

void
hp_bilinear_interpolate_weights(double lon, double lat, size_t *indices, double *weights, int nside, HpOrder order)
{
  int64_t (*xy_to_order)(int64_t index, int nside);
  // clang-format off
  xy_to_order = (order == HpOrder::Ring)   ? &healpixl_xy_to_ring :
                (order == HpOrder::Nested) ? &healpixl_xy_to_nested : &hp_xy_to_xy;
  // clang-format on

  int64_t ringIndices[4];
  interpolate_weights(lon, lat, ringIndices, weights, nside);

  if (order == HpOrder::Ring)
    for (int i = 0; i < 4; ++i) indices[i] = ringIndices[i];
  else
    for (int i = 0; i < 4; ++i) indices[i] = xy_to_order(healpixl_ring_to_xy(ringIndices[i], nside), nside);
}

void
hp_generate_coords(HpOrder order, int nside, int64_t nvals, double *xvals, double *yvals, bool withBounds, double *xbounds,
                   double *ybounds)
{
  int64_t (*order_to_xy)(int64_t index, int nside);
  // clang-format off
  order_to_xy = (order == HpOrder::Ring)   ? &healpixl_ring_to_xy :
                (order == HpOrder::Nested) ? &healpixl_nested_to_xy : &hp_xy_to_xy;
  // clang-format on

#ifdef _OPENMP
#pragma omp parallel for if (nvals > 99999) default(shared) schedule(static)
#endif
  for (int64_t index = 0; index < nvals; index++)
    {
      auto xyIndex = order_to_xy(index, nside);
      healpixl_to_radec(xyIndex, nside, 0.5, 0.5, &xvals[index], &yvals[index]);

      if (withBounds)
        {
          healpixl_to_radec(xyIndex, nside, 0.0, 0.0, &xbounds[index * 4 + 0], &ybounds[index * 4 + 0]);
          healpixl_to_radec(xyIndex, nside, 1.0, 0.0, &xbounds[index * 4 + 1], &ybounds[index * 4 + 1]);
          healpixl_to_radec(xyIndex, nside, 1.0, 1.0, &xbounds[index * 4 + 2], &ybounds[index * 4 + 2]);
          healpixl_to_radec(xyIndex, nside, 0.0, 1.0, &xbounds[index * 4 + 3], &ybounds[index * 4 + 3]);
        }
    }
}

void
hp_generate_latitudes(int nside, std::vector<double> &latitudes)
{
  for (int ringNumber = 1; ringNumber < 4 * nside; ringNumber++)
    {
      // Find the longitude/latitude and ring index of this pixel
      auto ringIndex = healpixl_compose_ring(ringNumber, 0, nside);
      auto xyIndex = healpixl_ring_to_xy(ringIndex, nside);
      double lon, lat;
      healpixl_to_radec(xyIndex, nside, 0.5, 0.5, &lon, &lat);
      latitudes[ringNumber - 1] = lat;
    }
}

static int
num_in_ring(int nside, int ringNumber)
{
  int numInRing = 0;
  // Now figure out again how many pixels are in the ring
  if (ringNumber < nside) { numInRing = 4 * ringNumber; }
  else if (ringNumber < 3 * nside) { numInRing = 4 * nside; }
  else { numInRing = (int) (4 * (4 * (int64_t) nside - (int64_t) ringNumber)); }
  return numInRing;
}

void
hp_generate_ring_indices(HpOrder order, int nside, size_t gridsize, std::vector<int> &ringIndices, std::vector<int> &ringRows)
{
  assert(gridsize <= INT_MAX && "Large grid size unsupported!");

  int64_t (*xy_to_order)(int64_t index, int nside);
  // clang-format off
  xy_to_order = (order == HpOrder::Ring)   ? &healpixl_xy_to_ring :
                (order == HpOrder::Nested) ? &healpixl_xy_to_nested : &hp_xy_to_xy;
  // clang-format on

  int numRows = 4 * nside - 1;
  if (ringRows.size() < (size_t)numRows) ringRows.resize(numRows);
  if (order != HpOrder::Ring && ringIndices.size() < gridsize) ringIndices.resize(gridsize);

  int index = 0;
  for (int ringNumber = 1; ringNumber < 4 * nside; ringNumber++)
    {
      auto numInRing = num_in_ring(nside, ringNumber);
      ringRows[ringNumber - 1] = numInRing;

      if (order != HpOrder::Ring)
        {
          for (int i = 0; i < numInRing; i++)
            {
              // Find the ring index of this pixel
              int ringIndex = healpixl_compose_ring(ringNumber, i, nside);
              if (order != HpOrder::Ring) ringIndex = xy_to_order(healpixl_ring_to_xy(ringIndex, nside), nside);
              ringIndices[index++] = ringIndex;
            }
        }
    }
}

template <typename T>
void hp_ring_to_nested(int nside, size_t gridsize, T *arrayIn, T *arrayOut)
{
  assert(gridsize <= INT_MAX && "Large grid size unsupported!");

  for (size_t i = 0; i < gridsize; ++i)
    {
      auto nestedIndex = healpixl_xy_to_nested(healpixl_ring_to_xy(i, nside), nside);
      arrayOut[nestedIndex] = arrayIn[i];
    }
}

// Explicit instantiation
template void hp_ring_to_nested(int nside, size_t gridsize, float *arrayIn, float *arrayOut);
template void hp_ring_to_nested(int nside, size_t gridsize, double *arrayIn, double *arrayOut);


template <typename T>
void hp_nested_to_ring(int nside, size_t gridsize, T *arrayIn, T *arrayOut)
{
  assert(gridsize <= INT_MAX && "Large grid size unsupported!");

  for (size_t i = 0; i < gridsize; ++i)
    {
      auto ringIndex = healpixl_xy_to_ring(healpixl_nested_to_xy(i, nside), nside);
      arrayOut[ringIndex] = arrayIn[i];
    }
}

// Explicit instantiation
template void hp_nested_to_ring(int nside, size_t gridsize, float *arrayIn, float *arrayOut);
template void hp_nested_to_ring(int nside, size_t gridsize, double *arrayIn, double *arrayOut);
