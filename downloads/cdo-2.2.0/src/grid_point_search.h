/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/
#ifndef GRID_POINT_SEARCH_H
#define GRID_POINT_SEARCH_H

#include <cstddef>
#include "knn_weights.h"
#include "varray.h"
#include "mpim_grid/grid_healpix.h"

#define GPS_NOT_FOUND SIZE_MAX

constexpr double
square(const double x) noexcept
{
  return x * x;
}

constexpr double
squareDistance(const double *a, const double *b) noexcept
{
  return square(a[0] - b[0]) + square(a[1] - b[1]) + square(a[2] - b[2]);
}

enum class PointSearchMethod
{
  undefined,
  full,
  nanoflann,
  kdtree,
  spherepart,
  latbins
};

struct GridPointSearch
{
  bool in_use = false;
  bool extrapolate = false;
  bool isCyclic = false;
  bool is_reg2d = false;
  bool is_curve = false;
  PointSearchMethod method = PointSearchMethod::nanoflann;
  size_t n = 0;
  size_t dims[2] = { 0 };

  void *search_container = nullptr;
  double searchRadius = 0;
  double searchArcRadius = 0;

  // healpix
  int nside = 0;
  HpOrder order = HpOrder::Undef;

  // reg2d search
  Varray<double> reg2d_center_lon, reg2d_center_lat;
  Varray<double> coslat, sinlat;  // cosine, sine of grid lats (for distance)
  Varray<double> coslon, sinlon;  // cosine, sine of grid lons (for distance)

  const double *plons = nullptr, *plats = nullptr;

  double lonmin = 0.0, lonmax = 0.0, latmin = 0.0, latmax = 0.0;
  float min[3] = { 0 }, max[3] = { 0 };
  void *pointCloud = nullptr;

  double (*coordinates_xyz)[3];
  ~GridPointSearch();
};

void grid_point_search_set_arc_radius(GridPointSearch &gps, double arcRadius);
void grid_point_search_set_chord_radius(GridPointSearch &gps, double chordRadius);

void grid_search_point(GridPointSearch &gps, double plon, double plat, knnWeightsType &knnWeights);
void grid_search_point_smooth(GridPointSearch &gps, double plon, double plat, knnWeightsType &knnWeights);

void grid_point_search_create_healpix(GridPointSearch &gps, size_t size, int nside, HpOrder order);
void grid_point_search_create_reg_2d(GridPointSearch &gps, bool xIsCyclic, size_t dims[2], const Varray<double> &lons,
                                     const Varray<double> &lats);
void grid_point_search_create(GridPointSearch &gps, bool xIsCyclic, size_t dims[2], size_t n, const Varray<double> &lons,
                              const Varray<double> &lats);
void grid_point_search_create(GridPointSearch &gps, const Varray<double> &lons, const Varray<double> &lats,
                              PointSearchMethod method = PointSearchMethod::nanoflann);
void grid_point_search_delete(GridPointSearch &gps);
size_t grid_point_search_nearest(const GridPointSearch &gps, double lon, double lat, size_t *addr, double *dist);
size_t grid_point_search_qnearest(const GridPointSearch &gps, double lon, double lat, size_t nnn, size_t *adds, double *dist);
size_t grid_point_search_distance_qnearest(const GridPointSearch &gps, double searchRadius, double lon, double lat, size_t nnn,
                                           size_t *adds, double *dist);
void grid_point_search_extrapolate(GridPointSearch &gps);

#endif
