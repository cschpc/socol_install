#ifndef GRID_HEALPIX_H
#define GRID_HEALPIX_H

#include <cassert>
#include <cinttypes>
#include <string>
#include <vector>

enum class HpOrder
{
  Undef,
  XY,
  Ring,
  Nested
};

HpOrder hp_get_order(const std::string &orderName);
int64_t hp_lonlat_to_index(HpOrder order, int nside, double xval, double yval);
void hp_index_to_lonlat(HpOrder order, int nside, int64_t index, double *xval, double *yval);
void hp_get_neighbours(HpOrder order, int nside, int64_t index, int64_t *neighbours);
void hp_bilinear_interpolate_weights(double lon, double lat, size_t *indices, double *weights, int nside, HpOrder order);
void hp_generate_coords(HpOrder order, int nside, int64_t nvals, double *xvals, double *yvals, bool withBounds, double *xbounds,
                        double *ybounds);
void hp_generate_latitudes(int nside, std::vector<double> &latitudes);
void hp_generate_ring_indices(HpOrder order, int nside, size_t gridsize, std::vector<int> &ringIndices, std::vector<int> &ringRows);

template <typename T>
void hp_ring_to_nested(int nside, size_t gridsize, T *arrayIn, T *arrayOut);

template <typename T>
void hp_nested_to_ring(int nside, size_t gridsize, T *arrayIn, T *arrayOut);

#endif
