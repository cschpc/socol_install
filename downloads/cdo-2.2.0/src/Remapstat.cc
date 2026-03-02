/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

*/

#include <cdi.h>

#include "cdo_wtime.h"
#include "process_int.h"
#include "griddes.h"
#include "grid_point_search.h"
#include <mpim_grid.h>
#include "gridreference.h"
#include "cimdOmp.h"
#include "verifygrid.h"
#include "field_functions.h"

// #define TESTIMPL 1
// #define USE_YAC  1
#ifdef TESTIMPL
#ifdef USE_YAC
extern "C"
{
#include "lib/yac/geometry.h"
}
#endif
#endif

constexpr double PI2 = M_PI * 2.0;

double radiusDegToKm(const double radiusInDeg);

struct StatInfo
{
  size_t n = 0;
  size_t min_svals = SIZE_MAX, max_svals = 0, sum_svals = 0;
  size_t min_nvals = SIZE_MAX, max_nvals = 0, sum_nvals = 0;
  double min_radius = 1.e33, max_radius = 0.0, sum_radius = 0.0;

  void
  add(size_t svals, size_t nvals, double radius)
  {
    radius *= RAD2DEG;
    n++;
    sum_svals += svals;
    sum_nvals += nvals;
    sum_radius += radius;
    min_svals = std::min(min_svals, svals);
    min_nvals = std::min(min_nvals, nvals);
    min_radius = std::min(min_radius, radius);
    max_svals = std::max(max_svals, svals);
    max_nvals = std::max(max_nvals, nvals);
    max_radius = std::max(max_radius, radius);
  };

  void
  print()
  {
    cdo_print("N=%zu", n);
    cdo_print("Min:   svals=%3zu  nvals=%2zu  radius=%.3gdeg(%.3gkm)", min_svals, min_nvals, min_radius, radiusDegToKm(min_radius));
    cdo_print("Mean:  svals=%3zu  nvals=%2zu  radius=%.3gdeg(%.3gkm)", sum_svals / n, sum_nvals / n, sum_radius / n,
              radiusDegToKm(sum_radius / n));
    cdo_print("Max:   svals=%3zu  nvals=%2zu  radius=%.3gdeg(%.3gkm)", max_svals, max_nvals, max_radius, radiusDegToKm(max_radius));
  };
};

static size_t
read_target_cell_bounds(int gridID, Varray<double> &xbounds, Varray<double> &ybounds)
{
  if (!gridHasBounds(gridID)) cdo_abort("Target cell corner coordinates missing!");

  auto nv = gridInqNvertex(gridID);
  auto gridsize = gridInqSize(gridID);
  xbounds.resize(nv * gridsize);
  ybounds.resize(nv * gridsize);

  gridInqXbounds(gridID, xbounds.data());
  gridInqYbounds(gridID, ybounds.data());

  // Convert lat/lon units if required
  cdo_grid_to_radian(gridID, CDI_XAXIS, nv * gridsize, xbounds.data(), "grid corner lon");
  cdo_grid_to_radian(gridID, CDI_YAXIS, nv * gridsize, ybounds.data(), "grid corner lat");

  return nv;
}

static void
read_coordinates(int gridID, Varray<double> &xvals, Varray<double> &yvals)
{
  auto gridID0 = gridID;

  gridID = generate_full_point_grid(gridID);
  if (!gridHasCoordinates(gridID)) cdo_abort("Cell center coordinates missing!");

  gridInqXvals(gridID, xvals.data());
  gridInqYvals(gridID, yvals.data());

  auto gridsize = gridInqSize(gridID);
  // Convert lat/lon units if required
  cdo_grid_to_radian(gridID, CDI_XAXIS, gridsize, xvals.data(), "grid center lon");
  cdo_grid_to_radian(gridID, CDI_YAXIS, gridsize, yvals.data(), "grid center lat");

  for (size_t i = 0; i < gridsize; ++i)
    {
      if (xvals[i] > PI2) xvals[i] -= PI2;
      if (xvals[i] < 0.0) xvals[i] += PI2;
    }

  if (gridID0 != gridID) gridDestroy(gridID);
}

static void
read_xbounds_reg2d(int gridID, std::vector<double> &xbounds2d)
{
  auto nlon = gridInqXsize(gridID);
  xbounds2d.resize(nlon * 2);

  if (gridInqXbounds(gridID, nullptr)) { gridInqXbounds(gridID, xbounds2d.data()); }
  else
    {
      std::vector<double> xvals(nlon);
      gridInqXvals(gridID, xvals.data());
      grid_gen_bounds(nlon, xvals, xbounds2d);
    }

  for (size_t i = 0; i < 2 * nlon; ++i) xbounds2d[i] *= DEG2RAD;

  for (size_t i = 0; i < nlon; ++i)
    {
      if (xbounds2d[2 * i + 1] > PI2)
        {
          xbounds2d[2 * i] -= PI2;
          xbounds2d[2 * i + 1] -= PI2;
        }
      if (xbounds2d[2 * i + 1] < 0.0)
        {
          xbounds2d[2 * i] += PI2;
          xbounds2d[2 * i + 1] += PI2;
        }
    }
}

static void
read_ybounds_reg2d(int gridID, std::vector<double> &ybounds2d)
{
  auto nlat = gridInqYsize(gridID);
  ybounds2d.resize(nlat * 2);

  if (gridInqYbounds(gridID, nullptr)) { gridInqYbounds(gridID, ybounds2d.data()); }
  else
    {
      std::vector<double> yvals(nlat);
      gridInqYvals(gridID, yvals.data());
      grid_gen_bounds(nlat, yvals, ybounds2d);
      grid_check_lat_borders(2 * nlat, ybounds2d.data());
    }

  for (size_t i = 0; i < 2 * nlat; ++i) ybounds2d[i] *= DEG2RAD;
}

static double
calc_maxdist(size_t i, size_t nv, double plon, double plat, const Varray<double> &xbounds, const Varray<double> &ybounds)
{
  double p1[3], p2[3];
  gcLLtoXYZ(plon, plat, p1);

  double maxdist = 0.0;
  for (size_t k = 0; k < nv; ++k)
    {
      auto lons = &xbounds[nv * i];
      auto lats = &ybounds[nv * i];
      gcLLtoXYZ(lons[k], lats[k], p2);
      auto sdist = squareDistance(p1, p2);
      maxdist = std::max(maxdist, sdist);
    }

  maxdist = 1.01 * std::sqrt(maxdist);

  return maxdist;
}

static double
calc_maxdist_rec2d(size_t i, size_t nlon, double plon, double plat, const std::vector<double> &xbounds,
                   const std::vector<double> &ybounds)
{
  double p1[3], p2[3];
  gcLLtoXYZ(plon, plat, p1);

  constexpr int nv = 4;
  double lons[nv], lats[nv];
  auto iy = i / nlon;
  auto ix = i - iy * nlon;
  lons[0] = xbounds[2 * ix];
  lons[1] = xbounds[2 * ix];
  lons[2] = xbounds[2 * ix + 1];
  lons[3] = xbounds[2 * ix + 1];
  lats[0] = ybounds[2 * iy + 1];
  lats[1] = ybounds[2 * iy];
  lats[2] = ybounds[2 * iy];
  lats[3] = ybounds[2 * iy + 1];

  double maxdist = 0.0;
  for (size_t k = 0; k < nv; ++k)
    {
      gcLLtoXYZ(lons[k], lats[k], p2);
      auto sdist = squareDistance(p1, p2);
      maxdist = std::max(maxdist, sdist);
    }
  maxdist = 1.01 * std::sqrt(maxdist);

  return maxdist;
}

#ifdef USE_YAC
static size_t
find_points_yac(std::vector<char> &vmask, size_t cell_no, size_t ncorner, size_t nadds, Varray<size_t> &adds,
                const Varray<double> &xvals1, const Varray<double> &yvals1, const Varray<double> &xbounds,
                const Varray<double> &ybounds)
{
#ifndef TESTIMPL
  cdo_abort("Internal error: find_points() not implemented!");
#endif

  struct grid_cell cell;
  cell.array_size = ncorner;
  cell.num_corners = ncorner;
  cell.edge_type = new enum yac_edge_type[ncorner];
  cell.coordinates_xyz = new double[ncorner][3];

  // enum yac_edge_type lonlat_circle_type[] = { LON_CIRCLE_EDGE, LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE, LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE
  // };

  for (size_t k = 0; k < ncorner; ++k)
    {
      auto lon = xbounds[cell_no * ncorner + k];
      auto lat = ybounds[cell_no * ncorner + k];
      cell.edge_type[k] = GREAT_CIRCLE_EDGE;
      // cell.edge_type[k] = lonlat_circle_type[k + 1];
      gcLLtoXYZ(lon, lat, cell.coordinates_xyz[k]);
    }

  double centerCoordinates[3];
  size_t nvalues = 0;
  for (size_t k = 0; k < nadds; ++k)
    {
      auto index1 = adds[k];
      auto lon = xvals1[index1];
      auto lat = yvals1[index1];

      gcLLtoXYZ(lon, lat, centerCoordinates);

      if (yac_point_in_cell(centerCoordinates, cell))
        {
          if (vmask[index1] == 0)
            {
              adds[nvalues] = adds[k];
              nvalues++;
            }
          if (vmask[index1] < 127) vmask[index1]++;
        }
    }

  delete[] cell.edge_type;
  delete[] cell.coordinates_xyz;

  return nvalues;
}
#endif

static size_t
find_points(std::vector<char> &vmask, size_t cell_no, size_t ncorner, size_t nadds, Varray<size_t> &adds,
            const Varray<double> &xvals1, const Varray<double> &yvals1, const Varray<double> &xbounds,
            const Varray<double> &ybounds)
{
#ifndef TESTIMPL
  cdo_abort("Internal error: find_points() not implemented!");
#endif

  Point3D centerPoint3D;
  Varray<Point> cellCornersPlaneProjection(ncorner + 1);
  Varray<Point3D> cellCorners3D(ncorner + 1);
  Varray<Point3D> cell_corners_xyz_open_cell(ncorner);
  std::vector<bool> marked_duplicate_indices(ncorner);

  set_cell_corners_3D(ncorner, &xbounds[cell_no * ncorner], &ybounds[cell_no * ncorner], cell_corners_xyz_open_cell);

  auto actualNumberOfCorners = get_actual_number_of_corners(ncorner, cell_corners_xyz_open_cell);

  auto no_duplicates = get_no_duplicates(actualNumberOfCorners, cell_corners_xyz_open_cell, marked_duplicate_indices);

  copy_unique_corners(actualNumberOfCorners, cell_corners_xyz_open_cell, marked_duplicate_indices, cellCorners3D);

  actualNumberOfCorners -= no_duplicates;

  cellCorners3D[actualNumberOfCorners] = cellCorners3D[0];

  if (actualNumberOfCorners < 3) return 0;

  auto coordinateToIgnore = find_coordinate_to_ignore(cellCorners3D);

  auto cval
      = (coordinateToIgnore == 1) ? cellCorners3D[0].X : ((coordinateToIgnore == 2) ? cellCorners3D[0].Y : cellCorners3D[0].Z);
  auto invertResult = (cval < 0.0);

  set_cell_corners_plane_projection(coordinateToIgnore, actualNumberOfCorners, cellCorners3D, cellCornersPlaneProjection);

  auto isClockwise = are_polygon_vertices_arranged_in_clockwise_order(cellCornersPlaneProjection, actualNumberOfCorners + 1);

  if (invertResult) isClockwise = !isClockwise;
  if (isClockwise) return 0;

  double centerCoordinates[3];
  size_t nvalues = 0;
  for (size_t k = 0; k < nadds; ++k)
    {
      auto index1 = adds[k];
      auto lon = xvals1[index1];
      auto lat = yvals1[index1];

      gcLLtoXYZ(lon, lat, centerCoordinates);
      centerPoint3D.X = centerCoordinates[0];
      centerPoint3D.Y = centerCoordinates[1];
      centerPoint3D.Z = centerCoordinates[2];

      auto centerPoint2D = set_center_point_plane_projection(coordinateToIgnore, centerPoint3D);

      auto windingNumber = winding_numbers_algorithm(cellCornersPlaneProjection, actualNumberOfCorners + 1, centerPoint2D);

      if (windingNumber != 0)
        {
          if (vmask[index1] == 0)
            {
              adds[nvalues] = adds[k];
              nvalues++;
            }
          if (vmask[index1] < 127) vmask[index1]++;
        }
    }

  return nvalues;
}

static size_t
find_points_rec2d(std::vector<char> &vmask, size_t i, size_t nlon2, size_t nadds, Varray<size_t> &adds,
                  const Varray<double> &xvals1, const Varray<double> &yvals1, const std::vector<double> &xbounds2d,
                  const std::vector<double> &ybounds2d)
{
  auto iy = i / nlon2;
  auto ix = i - iy * nlon2;

  size_t nvalues = 0;
  for (size_t k = 0; k < nadds; ++k)
    {
      auto index1 = adds[k];
      auto x = xvals1[index1];
      auto y = yvals1[index1];
      if (y >= ybounds2d[2 * iy] && y < ybounds2d[2 * iy + 1])
        if ((x >= xbounds2d[2 * ix] && x < xbounds2d[2 * ix + 1])
            || ((x - PI2) >= xbounds2d[2 * ix] && (x - PI2) < xbounds2d[2 * ix + 1]))
          {
            if (vmask[index1] == 0)
              {
                adds[nvalues] = index1;
                nvalues++;
              }
            if (vmask[index1] < 127) vmask[index1]++;
          }
    }

  return nvalues;
}

static void
check_vmask(std::vector<char> vmask)
{
  constexpr int max_vals = 128;
  size_t vm[max_vals] = { 0 };
  auto size = vmask.size();
  for (size_t i = 0; i < size; ++i)
    if (vmask[i] >= 0) vm[(int) vmask[i]]++;
  for (int i = 0; i < max_vals; ++i)
    if (vm[i]) cdo_print("Number of source values: %d --> %zu/%zu", i, vm[i], size);
  size_t sum = 0;
  for (int i = 1; i < max_vals; ++i) sum += vm[i];
  cdo_print("Sum of used source values:     %zu/%zu", sum, size);
}

static Varray2D<size_t>
gen_mapdata(int gridID1, int gridID2)
{
  auto gridsize1 = gridInqSize(gridID1);
  auto gridsize2 = gridInqSize(gridID2);

  Varray2D<size_t> mapdata(gridsize2);

  Varray<double> xvals1(gridsize1), yvals1(gridsize1);
  read_coordinates(gridID1, xvals1, yvals1);

  Varray<double> xvals2(gridsize2), yvals2(gridsize2);
  read_coordinates(gridID2, xvals2, yvals2);

  auto gridtype2 = gridInqType(gridID2);
  auto grid2_is_reg2d = (gridtype2 == GRID_GAUSSIAN || gridtype2 == GRID_LONLAT);

  auto nlon2 = grid2_is_reg2d ? gridInqXsize(gridID2) : 0;

  int gridID2x = -1;

  Varray<double> xbounds, ybounds;
  std::vector<double> xbounds2d, ybounds2d;
  size_t nv = 4;
  if (grid2_is_reg2d)
    {
      read_xbounds_reg2d(gridID2, xbounds2d);
      read_ybounds_reg2d(gridID2, ybounds2d);
    }
  else
    {
      gridID2 = gridID2x = generate_full_cell_grid(gridID2);
      nv = read_target_cell_bounds(gridID2, xbounds, ybounds);
    }

  // for (size_t i = 0; i < nlon2; ++i) printf("%zu %g %g\n", i+1, xbounds2d[2*i], xbounds2d[2*i+1]);
  // for (size_t i = 0; i < nlat2; ++i) printf("%zu %g %g\n", i+1, ybounds2d[2*i], ybounds2d[2*i+1]);

  auto start = Options::cdoVerbose ? cdo_get_wtime() : 0.0;

  auto xIsCyclic = false;
  size_t dims[2] = { gridsize1, 0 };
  GridPointSearch gps;
  grid_point_search_create(gps, xIsCyclic, dims, gridsize1, xvals1, yvals1);

  if (Options::cdoVerbose) cdo_print("Point search created: %.2f seconds", cdo_get_wtime() - start);

  start = Options::cdoVerbose ? cdo_get_wtime() : 0;

  auto ndist_max = gridsize1;
  if (gridsize1 > 1000000) ndist_max /= 4;
  std::vector<char> vmask(gridsize1, 0);
  Varray2D<size_t> adds_2D(Threading::ompNumThreads, Varray<size_t>(ndist_max));
  Varray2D<double> dist_2D(Threading::ompNumThreads, Varray<double>(ndist_max));

  StatInfo statInfo;

#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(dynamic)
#endif
  for (size_t i = 0; i < gridsize2; ++i)
    {
      auto ompthID = cdo_omp_get_thread_num();
      auto &adds = adds_2D[ompthID];
      auto &dist = dist_2D[ompthID];

      auto lon2 = xvals2[i];
      auto lat2 = yvals2[i];

      auto maxdist = grid2_is_reg2d ? calc_maxdist_rec2d(i, nlon2, lon2, lat2, xbounds2d, ybounds2d)
                                    : calc_maxdist(i, nv, lon2, lat2, xbounds, ybounds);

      auto nadds = grid_point_search_distance_qnearest(gps, maxdist, lon2, lat2, ndist_max, adds.data(), dist.data());
      // printf("%zu nadds %zu\n", i+1, nadds);

      auto nvalues = grid2_is_reg2d ? find_points_rec2d(vmask, i, nlon2, nadds, adds, xvals1, yvals1, xbounds2d, ybounds2d)
#ifdef USE_YAC
                                    : find_points_yac(vmask, i, nv, nadds, adds, xvals1, yvals1, xbounds, ybounds);
#else
                                    : find_points(vmask, i, nv, nadds, adds, xvals1, yvals1, xbounds, ybounds);
#endif

      if (nvalues)
        {
          mapdata[i].resize(nvalues);
          for (size_t k = 0; k < nvalues; ++k) mapdata[i][k] = adds[k];
        }

      // if (Options::cdoVerbose) printf("%zu nadds %zu nvalues %zu  maxdist %g\n", i+1, nadds, nvalues, maxdist);
      if (Options::cdoVerbose) statInfo.add(nadds, nvalues, maxdist);
    }

  if (Options::cdoVerbose) statInfo.print();

  if (Options::cdoVerbose) check_vmask(vmask);

  if (Options::cdoVerbose) cdo_print("Point search qnearest: %.2f seconds", cdo_get_wtime() - start);

  grid_point_search_delete(gps);

  if (gridID2x != -1) gridDestroy(gridID2x);

  return mapdata;
}

template <typename T>
static T
remap_kernel(int operfunc, const Varray<size_t> &adds, size_t &nmiss2, Field &field, Varray<T> &fieldvec, const Varray<T> &vec1,
             const T missval)
{
  T value;
  auto nvalues = adds.size();
  if (nvalues)
    {
      field.nmiss = 0;
      for (size_t k = 0; k < nvalues; ++k)
        {
          auto v1 = vec1[adds[k]];
          fieldvec[k] = v1;
          if (DBL_IS_EQUAL(v1, missval)) field.nmiss++;
        }

      field.size = nvalues;
      field.missval = missval;
      value = field_function(field, operfunc);
      if (DBL_IS_EQUAL(value, missval)) nmiss2++;
    }
  else
    {
      value = missval;
      nmiss2++;
    }

  return value;
}

static void
remap_field(const Varray2D<size_t> &mapdata, const Field &field1, Field &field2, int operfunc)
{
  std::vector<Field> fields(Threading::ompNumThreads);

  auto gridsize2 = gridInqSize(field2.grid);
  auto start = Options::cdoVerbose ? cdo_get_wtime() : 0.0;

  size_t nmiss2 = 0;
#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static) reduction(+ : nmiss2)
#endif
  for (size_t i = 0; i < gridsize2; ++i)
    {
      const auto &adds = mapdata[i];
      auto nvalues = adds.size();
      auto missval = field1.missval;

      auto ompthID = cdo_omp_get_thread_num();
      auto &field = fields[ompthID];
      field.memType = field1.memType;
      if (field1.memType == MemType::Float)
        field.vec_f.resize(nvalues);
      else
        field.vec_d.resize(nvalues);

      double rvalue = 0.0;
      if (field1.memType == MemType::Float)
        rvalue = remap_kernel(operfunc, adds, nmiss2, field, field.vec_f, field1.vec_f, (float) missval);
      else
        rvalue = remap_kernel(operfunc, adds, nmiss2, field, field.vec_d, field1.vec_d, missval);

      if (field2.memType == MemType::Float)
        field2.vec_f[i] = rvalue;
      else
        field2.vec_d[i] = rvalue;
    }

  field2.nmiss = nmiss2;

  if (Options::cdoVerbose) cdo_print("Remap: %.3f seconds", cdo_get_wtime() - start);
}

static void
add_operators(void)
{
  // clang-format off
  cdo_operator_add("remaprange",  FieldFunc_Range,  0, nullptr);
  cdo_operator_add("remapmin",    FieldFunc_Min,    0, nullptr);
  cdo_operator_add("remapmax",    FieldFunc_Max,    0, nullptr);
  cdo_operator_add("remapsum",    FieldFunc_Sum,    0, nullptr);
  cdo_operator_add("remapmean",   FieldFunc_Mean,   0, nullptr);
  cdo_operator_add("remapavg",    FieldFunc_Avg,    0, nullptr);
  cdo_operator_add("remapstd",    FieldFunc_Std,    0, nullptr);
  cdo_operator_add("remapstd1",   FieldFunc_Std1,   0, nullptr);
  cdo_operator_add("remapvar",    FieldFunc_Var,    0, nullptr);
  cdo_operator_add("remapvar1",   FieldFunc_Var1,   0, nullptr);
  cdo_operator_add("remapskew",   FieldFunc_Skew,   0, nullptr);
  cdo_operator_add("remapkurt",   FieldFunc_Kurt,   0, nullptr);
  cdo_operator_add("remapmedian", FieldFunc_Median, 0, nullptr);
  // clang-format on
}

void *
Remapstat(void *process)
{
  cdo_initialize(process);

  add_operators();

  auto operatorID = cdo_operator_id();
  auto operfunc = cdo_operator_f1(operatorID);

  operator_input_arg("grid description file or name");
  auto gridID2 = cdo_define_grid(cdo_operator_argv(0));
  auto gridtype2 = gridInqType(gridID2);

  {
#ifdef TESTIMPL
    auto hasProjParams = ((gridtype2 == GRID_PROJECTION) && grid_has_proj_params(gridID2));
    if (!gridProjIsSupported(gridID2) && !hasProjParams && gridtype2 != GRID_LONLAT && gridtype2 != GRID_GAUSSIAN
        && gridtype2 != GRID_CURVILINEAR && gridtype2 != GRID_UNSTRUCTURED)
#else
    if (gridtype2 != GRID_LONLAT && gridtype2 != GRID_GAUSSIAN)
#endif
      cdo_abort("Remapping to %s data unsupported!", gridNamePtr(gridtype2));
  }

  auto streamID1 = cdo_open_read(0);

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto vlistID2 = vlistDuplicate(vlistID1);

  auto taxisID1 = vlistInqTaxis(vlistID1);
  auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  auto ngrids = vlistNgrids(vlistID1);
  auto gridID1 = vlistGrid(vlistID1, 0);
  for (int index = 0; index < ngrids; ++index)
    {
      if (index > 0) cdo_abort("Too many different grids!");

      // auto gridID1 = vlistGrid(vlistID1, index);
      auto gridtype = gridInqType(gridID1);
      auto hasProjParams = ((gridtype == GRID_PROJECTION) && grid_has_proj_params(gridID1));
      if (!gridProjIsSupported(gridID1) && !hasProjParams && gridtype != GRID_LONLAT && gridtype != GRID_GAUSSIAN
          && gridtype != GRID_GME && gridtype != GRID_CURVILINEAR && gridtype != GRID_UNSTRUCTURED)
        cdo_abort("Interpolation of %s data unsupported!", gridNamePtr(gridtype));

      vlistChangeGridIndex(vlistID2, index, gridID2);
    }

  auto streamID2 = cdo_open_write(1);
  cdo_def_vlist(streamID2, vlistID2);

  Varray2D<size_t> mapdata = gen_mapdata(gridID1, gridID2);

  VarList varList1, varList2;
  varListInit(varList1, vlistID1);
  varListInit(varList2, vlistID2);

  Field field1, field2;

  int tsID = 0;
  while (true)
    {
      auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
      if (nrecs == 0) break;

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_def_timestep(streamID2, tsID);

      for (int recID = 0; recID < nrecs; ++recID)
        {
          int varID, levelID;
          cdo_inq_record(streamID1, &varID, &levelID);
          field1.init(varList1[varID]);
          cdo_read_record(streamID1, field1);

          field2.init(varList2[varID]);

          remap_field(mapdata, field1, field2, operfunc);

          cdo_def_record(streamID2, varID, levelID);
          cdo_write_record(streamID2, field2);
        }
      tsID++;
    }

  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  cdo_finish();

  return nullptr;
}
