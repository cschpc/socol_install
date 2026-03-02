/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

// std includes
#include <cstdio>

// sub lib includes
#include <mpim_grid.h>

// local includes
#include "cdo_output.h"
#include "cdo_math.h"
#include "cdo_options.h"
#include "cimdOmp.h"
#include "grid_point_search.h"
#include "kdtreelib/kdtree.h"
#include "nanoflann.hpp"
extern "C"
{
#include "lib/yac/sphere_part.h"
}

constexpr double PI = M_PI;
constexpr double PI2 = 2.0 * PI;

PointSearchMethod pointSearchMethod(PointSearchMethod::undefined);

struct gpsFull
{
  size_t n = 0;
  const double *plons = nullptr;
  const double *plats = nullptr;
  double (*pts)[3] = nullptr;
};

template <typename T>
struct PointCloud
{
  struct Point
  {
    T x, y, z;
  };
  std::vector<Point> pts;
  T min[3], max[3];

  // Must return the number of data points
  inline size_t
  kdtree_get_point_count() const
  {
    return pts.size();
  }

  // Returns the dim'th component of the idx'th point in the class:
  // Since this is inlined and the "dim" argument is typically an immediate
  // value, the
  //  "if/else's" are actually solved at compile time.
  inline T
  kdtree_get_pt(const size_t idx, int dim) const
  {
    // clang-format off
    if      (dim == 0) return pts[idx].x;
    else if (dim == 1) return pts[idx].y;
    else               return pts[idx].z;
    // clang-format on
  }

  // Optional bounding-box computation: return false to default to a standard
  // bbox computation loop.
  //   Return true if the BBOX was already computed by the class and returned in
  //   "bb" so it can be avoided to redo it again. Look at bb.size() to find out
  //   the expected dimensionality (e.g. 2 or 3 for point clouds)
  template <class BBOX>
  bool
  kdtree_get_bbox(BBOX &bb) const
  {
    for (int i = 0; i < 3; ++i) bb[i].low = min[i];
    for (int i = 0; i < 3; ++i) bb[i].high = max[i];
    return true;
  }
  // bool kdtree_get_bbox(BBOX& /* bb */) const { return false; }
};

using nfTree_t
    = nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, PointCloud<double>>, PointCloud<double>, 3>;

static double
cdoDefaultRadius(void)
{
  extern double pointSearchRadius;

  auto searchRadius = pointSearchRadius;
  searchRadius = std::min(std::max(searchRadius, 0.), 180.);
  searchRadius *= DEG2RAD;

  return searchRadius;
}

void
set_point_search_method(const std::string &methodstr)
{
  // clang-format off
  if      (methodstr == "kdtree")     pointSearchMethod = PointSearchMethod::kdtree;
  else if (methodstr == "nanoflann")  pointSearchMethod = PointSearchMethod::nanoflann;
  else if (methodstr == "spherepart") pointSearchMethod = PointSearchMethod::spherepart;
  else if (methodstr == "full")       pointSearchMethod = PointSearchMethod::full;
  else if (methodstr == "latbins")    pointSearchMethod = PointSearchMethod::latbins;
  else cdo_abort("Grid point search method %s not available!", methodstr);
  // clang-format on
}

static double
arc_to_chord_length(double arcLength)
{
  return 2.0 * std::sin(arcLength / 2.0);
}

static double
chord_to_arc_length(double chordLength)
{
  return 2.0 * std::asin(std::min(std::max(chordLength / 2.0, -1.0), 1.0));
}

void
grid_point_search_set_arc_radius(GridPointSearch &gps, double arcRadius)
{
  gps.searchArcRadius = arcRadius;
  gps.searchRadius = arc_to_chord_length(arcRadius);
}

void
grid_point_search_set_chord_radius(GridPointSearch &gps, double chordRadius)
{
  gps.searchRadius = chordRadius;
  gps.searchArcRadius = chord_to_arc_length(chordRadius);
}

void
grid_point_search_extrapolate(GridPointSearch &gps)
{
  gps.extrapolate = true;
}

void
grid_point_search_create_healpix(GridPointSearch &gps, size_t size, int nside, HpOrder order)
{
  gps.isCyclic = true;
  gps.dims[0] = size;
  gps.nside = nside;
  gps.order = order;

  gps.searchRadius = cdoDefaultRadius();

  gps.in_use = true;
}

void
grid_point_search_create_reg_2d(GridPointSearch &gps, bool xIsCyclic, size_t dims[2], const Varray<double> &lons,
                                const Varray<double> &lats)
{
  gps.isCyclic = xIsCyclic;
  gps.is_reg2d = true;
  gps.dims[0] = dims[0];
  gps.dims[1] = dims[1];
  auto nx = dims[0];
  auto ny = dims[1];

  auto nxm = xIsCyclic ? nx + 1 : nx;

  gps.reg2d_center_lon.resize(nxm);
  gps.reg2d_center_lat.resize(ny);

  varray_copy(nxm, lons, gps.reg2d_center_lon);
  varray_copy(ny, lats, gps.reg2d_center_lat);

  gps.coslon.resize(nx);
  gps.sinlon.resize(nx);
  gps.coslat.resize(ny);
  gps.sinlat.resize(ny);

  auto &coslon = gps.coslon;
  auto &sinlon = gps.sinlon;
  auto &coslat = gps.coslat;
  auto &sinlat = gps.sinlat;

  for (size_t n = 0; n < nx; ++n)
    {
      auto rlon = lons[n];
      if (rlon > PI2) rlon -= PI2;
      if (rlon < 0) rlon += PI2;
      coslon[n] = std::cos(rlon);
      sinlon[n] = std::sin(rlon);
    }

  for (size_t n = 0; n < ny; ++n)
    {
      coslat[n] = std::cos(lats[n]);
      sinlat[n] = std::sin(lats[n]);
    }

  gps.searchRadius = cdoDefaultRadius();

  gps.in_use = true;
}

static inline void
min_point(double *min, double *point)
{
  for (int i = 0; i < 3; ++i) min[i] = (point[i] < min[i]) ? point[i] : min[i];
}

static inline void
max_point(double *max, double *point)
{
  for (int i = 0; i < 3; ++i) max[i] = (point[i] > max[i]) ? point[i] : max[i];
}

template <typename T>
static void
adjust_bbox_min(T *min)
{
  for (int i = 0; i < 3; ++i) min[i] = (min[i] < 0) ? min[i] * 1.001 : min[i] * 0.999;
}

template <typename T>
static void
adjust_bbox_max(T *max)
{
  for (int i = 0; i < 3; ++i) max[i] = (max[i] < 0) ? max[i] * 0.999 : max[i] * 1.001;
}

static void *
gps_create_kdtree(size_t n, const Varray<double> &lons, const Varray<double> &lats, GridPointSearch &gps, bool useBoundBox = false)
{
  std::vector<kd_point> pointlist(n);
  // see  example_cartesian.c

  kdata_t min[3] = { 1.e9, 1.e9, 1.e9 };
  kdata_t max[3] = { -1.e9, -1.e9, -1.e9 };

#ifdef HAVE_OPENMP45
#pragma omp parallel for schedule(static) reduction(min : min[:3]) reduction(max : max[:3])
#endif
  for (size_t i = 0; i < n; ++i)
    {
      auto &point = pointlist[i].point;
      gcLLtoXYZ(lons[i], lats[i], point);
      min_point(min, point);
      max_point(max, point);
      pointlist[i].index = i;
    }

  if (!useBoundBox) min[0] = min[1] = min[2] = -1;
  if (!useBoundBox) max[0] = max[1] = max[2] = 1;

  adjust_bbox_min(min);
  adjust_bbox_max(max);
  for (int i = 0; i < 3; ++i) gps.min[i] = min[i];
  for (int i = 0; i < 3; ++i) gps.max[i] = max[i];

  if (Options::cdoVerbose) cdo_print("BBOX: min=%g/%g/%g  max=%g/%g/%g", min[0], min[1], min[2], max[0], max[1], max[2]);

  auto kdt = kd_buildTree(pointlist.data(), n, min, max, 3, Threading::ompNumThreads);
  if (kdt == nullptr) cdo_abort("kd_buildTree failed!");

  return (void *) kdt;
}

static void *
gps_create_nanoflann(size_t n, const Varray<double> &lons, const Varray<double> &lats, GridPointSearch &gps,
                     bool useBoundBox = false)
{
  PointCloud<double> *pointCloud = new PointCloud<double>();

  double min[3] = { 1.e9, 1.e9, 1.e9 };
  double max[3] = { -1.e9, -1.e9, -1.e9 };

  // Generating Point Cloud
  pointCloud->pts.resize(n);
#ifdef HAVE_OPENMP45
#pragma omp parallel for schedule(static) reduction(min : min[:3]) reduction(max : max[:3])
#endif
  for (size_t i = 0; i < n; ++i)
    {
      double point[3];
      gcLLtoXYZ(lons[i], lats[i], point);
      pointCloud->pts[i].x = point[0];
      pointCloud->pts[i].y = point[1];
      pointCloud->pts[i].z = point[2];
      min_point(min, point);
      max_point(max, point);
    }

  gps.pointCloud = (void *) pointCloud;

  if (!useBoundBox) min[0] = min[1] = min[2] = -1;
  if (!useBoundBox) max[0] = max[1] = max[2] = 1;

  adjust_bbox_min(min);
  adjust_bbox_max(max);
  for (int i = 0; i < 3; ++i) gps.min[i] = min[i];
  for (int i = 0; i < 3; ++i) gps.max[i] = max[i];

  if (Options::cdoVerbose) cdo_print("BBOX: min=%g/%g/%g  max=%g/%g/%g", min[0], min[1], min[2], max[0], max[1], max[2]);

  for (int i = 0; i < 3; ++i) pointCloud->min[i] = min[i];
  for (int i = 0; i < 3; ++i) pointCloud->max[i] = max[i];

  // construct a kd-tree index:
  nfTree_t *nft = new nfTree_t(3 /*dim*/, *pointCloud, nanoflann::KDTreeSingleIndexAdaptorParams(50 /* max leaf */));
  nft->buildIndex();

  return (void *) nft;
}

static void *
gps_create_spherepart(size_t n, const Varray<double> &lons, const Varray<double> &lats, GridPointSearch &gps,
                      bool useBoundBox = false)
{
  gps.coordinates_xyz = new double[n][3];

  double min[3] = { 1.e9, 1.e9, 1.e9 };
  double max[3] = { -1.e9, -1.e9, -1.e9 };

#ifdef HAVE_OPENMP45
#pragma omp parallel for schedule(static) reduction(min : min[:3]) reduction(max : max[:3])
#endif
  for (size_t i = 0; i < n; ++i)
    {
      double *point = gps.coordinates_xyz[i];
      gcLLtoXYZ(lons[i], lats[i], point);
      min_point(min, point);
      max_point(max, point);
    }

  if (!useBoundBox) min[0] = min[1] = min[2] = -1;
  if (!useBoundBox) max[0] = max[1] = max[2] = 1;

  adjust_bbox_min(min);
  adjust_bbox_max(max);
  for (int i = 0; i < 3; ++i) gps.min[i] = min[i];
  for (int i = 0; i < 3; ++i) gps.max[i] = max[i];

  if (Options::cdoVerbose) cdo_print("BBOX: min=%g/%g/%g  max=%g/%g/%g", min[0], min[1], min[2], max[0], max[1], max[2]);

  size_t *global_ids = new size_t[n];
  for (size_t i = 0; i < n; ++i) global_ids[i] = i;
  auto yacPointSearch = yac_point_sphere_part_search_new(n, gps.coordinates_xyz, global_ids);
  delete[] global_ids;
  return (void *) yacPointSearch;
}

static void
gps_destroy_kdtree(void *search_container)
{
  auto kdt = (kdTree_t *) search_container;
  if (kdt) kd_destroyTree(kdt);
}

static void
gps_destroy_full(void *search_container)
{
  auto full = (gpsFull *) search_container;
  if (full)
    {
      if (full->pts) delete[] full->pts;
      delete full;
    }
}

static void
gps_destroy_spherepart(void *search_container)
{
  yac_delete_point_sphere_part_search((point_sphere_part_search *) search_container);
}

static void *
gps_create_full(size_t n, const Varray<double> &lons, const Varray<double> &lats)
{
  if (Options::cdoVerbose) cdo_print("Init full grid search: n=%zu", n);

  auto full = new gpsFull;
  full->pts = new double[n][3];

#ifdef HAVE_OPENMP4
#pragma omp simd
#endif
  for (size_t i = 0; i < n; ++i) { gcLLtoXYZ(lons[i], lats[i], full->pts[i]); }

  full->n = n;
  full->plons = lons.data();
  full->plats = lats.data();

  return (void *) full;
}

static void
print_method(PointSearchMethod method)
{
  // clang-format off
  if      (method == PointSearchMethod::kdtree)     cdo_print("Point search method: kdtree");
  else if (method == PointSearchMethod::full)       cdo_print("Point search method: full");
  else if (method == PointSearchMethod::nanoflann)  cdo_print("Point search method: nanoflann");
  else if (method == PointSearchMethod::spherepart) cdo_print("Point search method: spherepart");
  // clang-format on
}

void
grid_point_search_create(GridPointSearch &gps, const Varray<double> &lons, const Varray<double> &lats, PointSearchMethod method)
{
  auto n = lons.size();

  gps.isCyclic = false;
  gps.is_curve = false;
  gps.dims[0] = n;
  gps.dims[1] = 0;

  gps.method = method;
  if (pointSearchMethod != PointSearchMethod::undefined) gps.method = pointSearchMethod;
  if (gps.method == PointSearchMethod::latbins) gps.method = PointSearchMethod::spherepart;

  gps.n = n;
  if (n == 0) return;

  gps.plons = lons.data();
  gps.plats = lats.data();

  if (Options::cdoVerbose) print_method(gps.method);

  // clang-format off
  if      (gps.method == PointSearchMethod::kdtree)     gps.search_container = gps_create_kdtree(n, lons, lats, gps);
  else if (gps.method == PointSearchMethod::full)       gps.search_container = gps_create_full(n, lons, lats);
  else if (gps.method == PointSearchMethod::nanoflann)  gps.search_container = gps_create_nanoflann(n, lons, lats, gps);
  else if (gps.method == PointSearchMethod::spherepart) gps.search_container = gps_create_spherepart(n, lons, lats, gps);
  else cdo_abort("%s::method undefined!", __func__);
  // clang-format on

  gps.searchRadius = cdoDefaultRadius();

  gps.in_use = true;
}

void
grid_point_search_create(GridPointSearch &gps, bool xIsCyclic, size_t dims[2], size_t n, const Varray<double> &lons,
                         const Varray<double> &lats)
{
  gps.isCyclic = xIsCyclic;
  gps.is_curve = (n != 1 && n == dims[0] * dims[1]);
  gps.dims[0] = dims[0];
  gps.dims[1] = dims[1];

  if (pointSearchMethod != PointSearchMethod::undefined) gps.method = pointSearchMethod;
  if (gps.method == PointSearchMethod::latbins) gps.method = PointSearchMethod::nanoflann;

  gps.n = n;
  if (n == 0) return;

  gps.plons = lons.data();
  gps.plats = lats.data();

  if (Options::cdoVerbose) print_method(gps.method);

  // clang-format off
  constexpr auto useBoundBox = true;
  if      (gps.method == PointSearchMethod::kdtree)     gps.search_container = gps_create_kdtree(n, lons, lats, gps, useBoundBox);
  else if (gps.method == PointSearchMethod::full)       gps.search_container = gps_create_full(n, lons, lats);
  else if (gps.method == PointSearchMethod::nanoflann)  gps.search_container = gps_create_nanoflann(n, lons, lats, gps, useBoundBox);
  else if (gps.method == PointSearchMethod::spherepart) gps.search_container = gps_create_spherepart(n, lons, lats, gps, useBoundBox);
  else cdo_abort("%s::method undefined!", __func__);
  // clang-format on

  gps.searchRadius = cdoDefaultRadius();

  gps.in_use = true;
}

void
grid_point_search_delete(GridPointSearch &gps)
{
  if (gps.in_use)
    {
      varray_free(gps.reg2d_center_lon);
      varray_free(gps.reg2d_center_lat);

      varray_free(gps.coslat);
      varray_free(gps.coslon);
      varray_free(gps.sinlat);
      varray_free(gps.sinlon);

      // clang-format off
      if      (gps.method == PointSearchMethod::kdtree)     gps_destroy_kdtree(gps.search_container);
      else if (gps.method == PointSearchMethod::nanoflann)
        {
          delete ((PointCloud<double> *) gps.pointCloud);
          delete ((nfTree_t *) gps.search_container);
        }
      else if (gps.method == PointSearchMethod::spherepart)
        {
          delete[] gps.coordinates_xyz;
          gps_destroy_spherepart(gps.search_container);
        }
      else if (gps.method == PointSearchMethod::full)       gps_destroy_full(gps.search_container);
      // clang-format on

      gps.search_container = nullptr;
      gps.in_use = false;
    }
}

static size_t
gps_nearest_kdtree(void *search_container, double lon, double lat, double searchRadius, size_t *addr, double *dist,
                   const GridPointSearch &gps)
{
  auto kdt = (kdTree_t *) search_container;
  if (kdt == nullptr) return 0;

  auto sqrDistMax = cdo::sqr(searchRadius);
  auto sqrDist = sqrDistMax;

  double query_pt[3];
  gcLLtoXYZ(lon, lat, query_pt);

  if (!gps.extrapolate)
    for (int i = 0; i < 3; ++i)
      if (query_pt[i] < gps.min[i] || query_pt[i] > gps.max[i]) return 0;

  auto node = kd_nearest(kdt->node, query_pt, &sqrDist, 3);

  if (node && sqrDist < sqrDistMax)
    {
      *addr = node->index;
      *dist = std::sqrt(sqrDist);
      return 1;
    }

  return 0;
}

static size_t
gps_nearest_nanoflann(void *search_container, double lon, double lat, double searchRadius, size_t *addr, double *dist,
                      const GridPointSearch &gps)
{
  auto nft = (nfTree_t *) search_container;
  if (nft == nullptr) return 0;

  auto sqrDistMax = cdo::sqr(searchRadius);

  double query_pt[3];
  gcLLtoXYZ(lon, lat, query_pt);

  if (!gps.extrapolate)
    for (int i = 0; i < 3; ++i)
      if (query_pt[i] < gps.min[i] || query_pt[i] > gps.max[i]) return 0;

  const size_t num_results = 1;
  size_t retIndex;
  double sqrDist;
  nanoflann::KNNResultSet<double> resultSet(sqrDistMax, num_results);
  resultSet.init(&retIndex, &sqrDist);
  nft->findNeighbors(resultSet, query_pt, nanoflann::SearchParams(10));

  if (retIndex != GPS_NOT_FOUND)
    {
      *addr = retIndex;
      *dist = std::sqrt(sqrDist);
      return 1;
    }

  return 0;
}

static size_t
gps_nearest_spherepart(void *search_container, double lon, double lat, double searchArcRadius, size_t *addr, double *dist,
                       const GridPointSearch &gps)
{
  double query_pt[1][3];
  gcLLtoXYZ(lon, lat, query_pt[0]);

  if (!gps.extrapolate)
    for (int i = 0; i < 3; ++i)
      if (query_pt[0][i] < gps.min[i] || query_pt[0][i] > gps.max[i]) return 0;

  size_t local_point_ids_array_size = 0;
  size_t num_local_point_ids;
  size_t *local_point_ids = nullptr;
  double cos_angle;

  yac_point_sphere_part_search_NN((point_sphere_part_search *) search_container, 1, query_pt, &cos_angle, nullptr, nullptr,
                                  &local_point_ids, &local_point_ids_array_size, &num_local_point_ids);

  size_t nadd = 0;
  if (num_local_point_ids > 0)
    {
      *dist = std::acos(cos_angle);
      if (*dist <= searchArcRadius)
        {
          nadd = 1;
          *addr = local_point_ids[0];
          for (size_t i = 1; i < num_local_point_ids; ++i)
            if (local_point_ids[i] < *addr) *addr = local_point_ids[i];
        }
    }

  if (local_point_ids) free(local_point_ids);

  return nadd;
}

static size_t
gps_nearest_full(void *search_container, double lon, double lat, double searchRadius, size_t *addr, double *dist)
{
  auto full = (gpsFull *) search_container;
  if (full == nullptr) return 0;

  auto sqrDistMax = cdo::sqr(searchRadius);

  double query_pt[3];
  gcLLtoXYZ(lon, lat, query_pt);

  size_t n = full->n;
  size_t closestpt = n;
  auto pts = full->pts;
  double sqrDist = FLT_MAX;
  for (size_t i = 0; i < n; ++i)
    {
      double d = (float) squareDistance(query_pt, pts[i]);
      if (closestpt >= n || d < sqrDist || (d <= sqrDist && i < closestpt))
        {
          sqrDist = d;
          closestpt = i;
        }
    }

  if (closestpt < n && sqrDist < sqrDistMax)
    {
      *addr = closestpt;
      *dist = std::sqrt(sqrDist);
      return 1;
    }

  return 0;
}

bool point_in_quad(bool isCyclic, size_t nx, size_t ny, size_t i, size_t j, size_t adds[4], double lons[4], double lats[4],
                   double plon, double plat, const double *centerLon, const double *centerLat);

static size_t
llindex_in_quad(const GridPointSearch &gps, size_t index, double lon, double lat)
{
  if (index != GPS_NOT_FOUND)
    {
      auto nx = gps.dims[0];
      auto ny = gps.dims[1];
      size_t adds[4];
      double lons[4], lats[4];
      auto isCyclic = gps.isCyclic;
      for (int k = 0; k < 4; ++k)
        {
          // Determine neighbor addresses
          auto j = index / nx;
          auto i = index - j * nx;
          if (k == 1 || k == 3) i = (i > 0) ? i - 1 : (isCyclic) ? nx - 1 : 0;
          if (k == 2 || k == 3) j = (j > 0) ? j - 1 : 0;

          if (point_in_quad(isCyclic, nx, ny, i, j, adds, lons, lats, lon, lat, gps.plons, gps.plats)) return index;
        }
    }

  return GPS_NOT_FOUND;
}

size_t
grid_point_search_nearest(const GridPointSearch &gps, double lon, double lat, size_t *addr, double *dist)
{
  if (gps.in_use)
    {
      size_t nadds = 0;
      auto searchRadius = gps.searchRadius;
      void *sc = gps.search_container;
      // clang-format off
      if      (gps.method == PointSearchMethod::kdtree)     nadds = gps_nearest_kdtree(sc, lon, lat, searchRadius, addr, dist, gps);
      else if (gps.method == PointSearchMethod::nanoflann)  nadds = gps_nearest_nanoflann(sc, lon, lat, searchRadius, addr, dist, gps);
      else if (gps.method == PointSearchMethod::spherepart) nadds = gps_nearest_spherepart(sc, lon, lat, chord_to_arc_length(searchRadius), addr, dist, gps);
      else if (gps.method == PointSearchMethod::full)       nadds = gps_nearest_full(sc, lon, lat, searchRadius, addr, dist);
      else cdo_abort("%s::method undefined!", __func__);
      // clang-format on

      if (nadds > 0)
        {
          auto index = *addr;
          if (!gps.extrapolate && gps.is_curve) index = llindex_in_quad(gps, *addr, lon, lat);
          if (index != GPS_NOT_FOUND) return 1;
        }
    }

  return 0;
}

static size_t
gps_qnearest_kdtree(const GridPointSearch &gps, double lon, double lat, double searchRadius, size_t nnn, size_t *adds, double *dist)
{
  size_t nadds = 0;

  auto kdt = (kdTree_t *) gps.search_container;
  if (kdt == nullptr) return nadds;

  auto sqrDistMax = cdo::sqr(searchRadius);

  kdata_t query_pt[3];
  gcLLtoXYZ(lon, lat, query_pt);

  if (!gps.extrapolate)
    for (int i = 0; i < 3; ++i)
      if (query_pt[i] < gps.min[i] || query_pt[i] > gps.max[i]) return nadds;

  if (gps.in_use)
    {
      kdata_t sqrDist = sqrDistMax;
      auto result = kd_qnearest(kdt->node, query_pt, &sqrDist, nnn, 3);
      if (result)
        {
          resItem *p;
          while (pqremove_min(result, &p))
            {
              if (p->dist_sq < sqrDistMax)
                {
                  adds[nadds] = p->node->index;
                  dist[nadds] = std::sqrt(p->dist_sq);
                  nadds++;
                }

              free(p);  // Free the result node taken from the heap
            }
          free(result->d);  // free the heap
          free(result);     // and free the heap information structure
        }
    }

  return nadds;
}

static size_t
gps_qnearest_nanoflann(const GridPointSearch &gps, double lon, double lat, double searchRadius, size_t nnn, size_t *adds,
                       double *dist)
{
  size_t nadds = 0;

  auto nft = (nfTree_t *) gps.search_container;
  if (nft == nullptr) return nadds;

  auto sqrDistMax = cdo::sqr(searchRadius);

  double query_pt[3];
  gcLLtoXYZ(lon, lat, query_pt);

  if (!gps.extrapolate)
    for (int i = 0; i < 3; ++i)
      if (query_pt[i] < gps.min[i] || query_pt[i] > gps.max[i]) return nadds;

  nadds = nft->knnRangeSearch(&query_pt[0], sqrDistMax, nnn, &adds[0], &dist[0]);
  for (size_t i = 0; i < nadds; ++i) dist[i] = std::sqrt(dist[i]);

  return nadds;
}

static size_t
gps_qnearest_spherepart(const GridPointSearch &gps, double lon, double lat, double searchArcRadius, size_t nnn, size_t *adds,
                        double *dist)
{
  size_t nadds = 0;

  if (gps.in_use)
    {
      double query_pt[3];
      gcLLtoXYZ(lon, lat, query_pt);

      if (!gps.extrapolate)
        for (int i = 0; i < 3; ++i)
          if (query_pt[i] < gps.min[i] || query_pt[i] > gps.max[i]) return nadds;

      size_t local_point_ids_array_size = 0;
      size_t num_local_point_ids;
      size_t *local_point_ids = nullptr;

      size_t cos_angles_array_size = 0;
      double *cos_angles = nullptr;

      yac_point_sphere_part_search_NNN((point_sphere_part_search *) gps.search_container, 1, &query_pt, nnn, &cos_angles,
                                       &cos_angles_array_size, nullptr, nullptr, &local_point_ids, &local_point_ids_array_size,
                                       &num_local_point_ids);

      if (num_local_point_ids > 0)
        {
          auto maxadds = (num_local_point_ids < nnn) ? num_local_point_ids : nnn;
          nadds = 0;
          for (size_t i = 0; i < maxadds; ++i)
            {
              auto angle = std::acos(cos_angles[i]);
              if (angle < searchArcRadius)
                {
                  adds[nadds] = local_point_ids[i];
                  dist[nadds] = angle;
                  nadds++;
                }
            }
        }

      if (cos_angles) free(cos_angles);
      if (local_point_ids) free(local_point_ids);
    }

  return nadds;
}

size_t
grid_point_search_qnearest(const GridPointSearch &gps, double lon, double lat, size_t nnn, size_t *adds, double *dist)
{
  size_t nadds = 0;

  if (gps.in_use)
    {
      auto searchRadius = gps.searchRadius;
      // clang-format off
      if      (gps.method == PointSearchMethod::kdtree)     nadds = gps_qnearest_kdtree(gps, lon, lat, searchRadius, nnn, adds, dist);
      else if (gps.method == PointSearchMethod::nanoflann)  nadds = gps_qnearest_nanoflann(gps, lon, lat, searchRadius, nnn, adds, dist);
      else if (gps.method == PointSearchMethod::spherepart) nadds = gps_qnearest_spherepart(gps, lon, lat, chord_to_arc_length(searchRadius), nnn, adds, dist);
      else cdo_abort("%s::method undefined!", __func__);
      // clang-format on

      if (!gps.extrapolate && gps.is_curve)
        {
          auto naddsmax = nadds;
          nadds = 0;
          for (size_t i = 0; i < naddsmax; ++i)
            {
              auto index = llindex_in_quad(gps, adds[i], lon, lat);
              if (index != GPS_NOT_FOUND)
                {
                  adds[nadds] = adds[i];
                  dist[nadds] = dist[i];
                  nadds++;
                }
            }
        }
    }

  return nadds;
}

size_t
grid_point_search_distance_qnearest(const GridPointSearch &gps, double searchRadius, double lon, double lat, size_t nnn,
                                    size_t *adds, double *dist)
{
  size_t nadds = 0;

  if (gps.in_use)
    {
      // clang-format off
      if      (gps.method == PointSearchMethod::kdtree)     nadds = gps_qnearest_kdtree(gps, lon, lat, searchRadius, nnn, adds, dist);
      else if (gps.method == PointSearchMethod::nanoflann)  nadds = gps_qnearest_nanoflann(gps, lon, lat, searchRadius, nnn, adds, dist);
      else if (gps.method == PointSearchMethod::spherepart) nadds = gps_qnearest_spherepart(gps, lon, lat, searchRadius, nnn, adds, dist);
      else cdo_abort("%s::method undefined!", __func__);
      // clang-format on

      if (!gps.extrapolate && gps.is_curve)
        {
          auto naddsmax = nadds;
          nadds = 0;
          for (size_t i = 0; i < naddsmax; ++i)
            {
              auto index = llindex_in_quad(gps, adds[i], lon, lat);
              if (index != GPS_NOT_FOUND)
                {
                  adds[nadds] = adds[i];
                  dist[nadds] = dist[i];
                  nadds++;
                }
            }
        }
    }

  return nadds;
}

GridPointSearch::~GridPointSearch()
{
  if (search_container != nullptr) { cdo_warning("search_container was not freed or not set to nullptr after free!"); }
}
