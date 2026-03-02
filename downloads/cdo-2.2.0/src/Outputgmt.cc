/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

   Output field with grid cell center or cell bounds for plotting with GMT

    - outputcenter
    - outputbounds
    - outputboundscpt
    - outputvector
*/

#ifdef HAVE_CONFIG_H
#include "config.h" /* VERSION */
#endif

#include <cdi.h>

#include "varray.h"
#include "cdo_options.h"
#include "process_int.h"
#include "param_conversion.h"
#include <mpim_grid.h>
#include "color.h"
#include "printinfo.h"
#include "interpol.h"
#include "cdo_zaxis.h"

static int
check_ncorner(int ncorner, const double *lon_bounds, const double *lat_bounds)
{
  auto ncorner_new = ncorner;

  int k;
  for (k = ncorner - 1; k > 0; --k)
    if (IS_NOT_EQUAL(lon_bounds[k], lon_bounds[k - 1]) || IS_NOT_EQUAL(lat_bounds[k], lat_bounds[k - 1])) break;

  if (k < ncorner - 1) ncorner_new = k + 1;

  return ncorner_new;
}

static void
check_lonbounds(int ncorner, double *lon_bounds)
{
  auto isLtM90 = false;
  auto isGtP90 = false;
  for (int k = 0; k < ncorner; ++k)
    {
      if (lon_bounds[k] < -90.0) isLtM90 = true;
      if (lon_bounds[k] > 90.0) isGtP90 = true;
    }

  if (isLtM90 && isGtP90)
    {
      for (int k = 0; k < ncorner; ++k)
        if (lon_bounds[k] < -90.0) lon_bounds[k] += 360.0;

      // printf("XXXX %g %g %g\n", lon_bounds[0], lon_bounds[1], lon_bounds[2]);
    }
}

static void
make_cyclic(double *array1, double *array2, long nlon, long nlat)
{
  for (long j = 0; j < nlat; ++j)
    for (long i = 0; i < nlon; ++i)
      {
        auto ij1 = j * nlon + i;
        auto ij2 = j * (nlon + 1) + i;
        array2[ij2] = array1[ij1];
      }

  for (long j = 0; j < nlat; ++j)
    {
      auto ij2 = j * (nlon + 1);
      array2[ij2 + nlon] = array2[ij2];
    }
}

static void
output_zon(double levmin, double levmax, const double *cell_corner_lat)
{
  auto latmin = cell_corner_lat[0];
  auto latmax = cell_corner_lat[0];
  for (int ic = 1; ic < 4; ++ic) latmin = std::min(latmin, cell_corner_lat[ic]);
  for (int ic = 1; ic < 4; ++ic) latmax = std::max(latmax, cell_corner_lat[ic]);
  const double xlev[4] = { levmin, levmax, levmax, levmin };
  const double xlat[4] = { latmin, latmin, latmax, latmax };
  for (int ic = 0; ic < 4; ++ic) fprintf(stdout, "   %g  %g\n", xlat[ic], xlev[ic]);
  fprintf(stdout, "   %g  %g\n", xlat[0], xlev[0]);
}

static void
output_mer(double levmin, double levmax, const double *cell_corner_lon)
{
  auto lonmin = cell_corner_lon[0];
  auto lonmax = cell_corner_lon[0];
  for (int ic = 1; ic < 4; ++ic) lonmin = std::min(lonmin, cell_corner_lon[ic]);
  for (int ic = 1; ic < 4; ++ic) lonmax = std::max(lonmax, cell_corner_lon[ic]);
  const double xlev[4] = { levmin, levmin, levmax, levmax };
  const double xlon[4] = { lonmin, lonmax, lonmax, lonmin };
  for (int ic = 0; ic < 4; ++ic) fprintf(stdout, "   %g  %g\n", xlon[ic], xlev[ic]);
  fprintf(stdout, "   %g  %g\n", xlon[0], xlev[0]);
}

static const int *
get_rgb(double value, double missval, const CPT &cpt)
{
  if (!DBL_IS_EQUAL(value, missval))
    {
      int n;
      for (n = 0; n < cpt.ncolors; ++n)
        if (value > cpt.lut[n].z_low && value <= cpt.lut[n].z_high) break;

      return (n == cpt.ncolors) ? cpt.bfn[0].rgb : cpt.lut[n].rgb_high;
    }
  else { return cpt.bfn[2].rgb; }
}

static void
output_vrml(long nlon, long nlat, long ngp, const Varray<double> &array, double missval, const CPT &cpt)
{
  auto mm = varray_min_max_mv(ngp, array, missval);
  auto dx = 10.0 / nlon;

  printf("Viewpoint {\n");
  printf("  description \"viewpoint1\"\n");
  printf("  orientation 0 0 1 0\n");
  printf("  position 0.0 0.0 10.0\n");
  printf("}\n");
  printf("\n");
  printf("Background {\n");
  printf("  skyColor [\n");
  printf("    0.0 0.1 0.8,\n");
  printf("    0.0 0.5 1.0,\n");
  printf("    1.0 1.0 1.0\n");
  printf("  ]\n");
  printf("  skyAngle [0.785, 1.571]\n");
  printf("\n");
  printf("  groundColor [\n");
  printf("    0.0 0.0 0.0,\n");
  printf("    0.3 0.3 0.3,\n");
  printf("    0.5 0.5 0.5\n");
  printf("  ]\n");
  printf("  groundAngle [0.785, 1.571]\n");
  printf("}\n");
  printf("\n");
  printf("Transform {\n");
  printf("  children [\n");
  printf("    Shape {\n");
  printf("      appearance Appearance {\n");
  printf("        material Material {}\n");
  printf("      }\n");
  printf("      geometry ElevationGrid {\n");
  printf("        colorPerVertex true\n");
  printf("        solid false\n");
  printf("        xDimension %ld\n", nlon);
  printf("        zDimension %ld\n", nlat);
  printf("        xSpacing %g\n", dx);
  printf("        zSpacing %g\n", dx);
  printf("        color Color {\n");
  printf("          color [\n");
  for (long j = nlat - 1; j >= 0; --j)
    for (long i = 0; i < nlon; ++i)
      {
        int r = 0, g = 0, b = 0;
        auto val = array[j * nlon + i];

        if (!DBL_IS_EQUAL(val, missval))
          {
            int n;
            for (n = 0; n < cpt.ncolors; ++n)
              if (val > cpt.lut[n].z_low && val <= cpt.lut[n].z_high) break;

            if (n == cpt.ncolors)
              {
                r = cpt.bfn[0].rgb[0];
                g = cpt.bfn[0].rgb[1];
                b = cpt.bfn[0].rgb[2];
              }
            else
              {
                //  r = cpt.lut[n].rgb_high[0];  g = cpt.lut[n].rgb_high[1];  b = cpt.lut[n].rgb_high[2];
                r = intlin(val, cpt.lut[n].rgb_low[0], cpt.lut[n].z_low, cpt.lut[n].rgb_high[0], cpt.lut[n].z_high);
                g = intlin(val, cpt.lut[n].rgb_low[1], cpt.lut[n].z_low, cpt.lut[n].rgb_high[1], cpt.lut[n].z_high);
                b = intlin(val, cpt.lut[n].rgb_low[2], cpt.lut[n].z_low, cpt.lut[n].rgb_high[2], cpt.lut[n].z_high);
              }
          }
        else
          {
            r = cpt.bfn[2].rgb[0];
            g = cpt.bfn[2].rgb[1];
            b = cpt.bfn[2].rgb[2];
          }
        printf(" %.3g %.3g %.3g,\n", r / 255., g / 255., b / 255.);
      }
  printf("          ]\n");
  printf("        }\n");
  printf("        height [\n");

  for (long j = nlat - 1; j >= 0; --j)
    for (long i = 0; i < nlon; ++i) printf("%g,\n", array[j * nlon + i]);

  printf("        ]\n");
  printf("      }\n");
  printf("    }\n");
  printf("  ]\n");
  printf("  translation -5 0 %g\n", -5. * nlat / nlon);
  printf("  rotation 0.0 0.0 0.0 0.0\n");
  printf("  scale 1.0 %g 1.0\n", 0.5 / (mm.max - mm.min));
  printf("}\n");
}

static void
output_kml(size_t ngp, const Varray<double> &array, int ncorner, const Varray<double> &grid_corner_lat,
           Varray<double> &grid_corner_lon, const std::vector<int> &grid_mask, double missval, const CPT &cpt)
{
  fprintf(stdout, "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n");
  fprintf(stdout, "<kml xmlns=\"http://www.opengis.net/kml/2.2\" \n");
  fprintf(stdout, "     xmlns:atom=\"http://www.w3.org/2005/Atom\">\n");
  fprintf(stdout, "<Document>\n");
  fprintf(stdout, " <name>CDO plot</name>\n");
  fprintf(stdout, " <open>1</open>\n");
#ifdef VERSION
  fprintf(stdout, " <atom:generator>CDO %s</atom:generator>\n", VERSION);
#endif
  fprintf(stdout, " <description>\n");
  fprintf(stdout, " <![CDATA[Generated by CDO]]>\n");
  fprintf(stdout, " </description>\n");
  fprintf(stdout, " <LookAt>\n");
  fprintf(stdout, "	<longitude>0</longitude>\n");
  fprintf(stdout, "	<latitude>40</latitude>\n");
  fprintf(stdout, "	<range>6e+06</range>\n");
  fprintf(stdout, "	<tilt>0</tilt>\n");
  fprintf(stdout, "	<heading>0</heading>\n");
  fprintf(stdout, "	<altitudeMode>absolute</altitudeMode>\n");
  fprintf(stdout, " </LookAt>\n");
  fprintf(stdout, " <Style id=\"check-hide-children\">\n");
  fprintf(stdout, "  <ListStyle>\n");
  fprintf(stdout, "   <listItemType>checkHideChildren</listItemType>\n");
  fprintf(stdout, "  </ListStyle>\n");
  fprintf(stdout, " </Style>\n");
  fprintf(stdout, "<Folder>\n");
  fprintf(stdout, "<name>Layer:Page</name>\n");
  fprintf(stdout, "<open>0</open>\n");
  fprintf(stdout, " <styleUrl>#check-hide-children</styleUrl>\n");
  fprintf(stdout, "<TimeSpan>\n");
  fprintf(stdout, " <begin></begin>\n");
  fprintf(stdout, " <end></end>\n");
  fprintf(stdout, "</TimeSpan>\n");
  fprintf(stdout, "<description><![CDATA[Layer:Page]]></description>\n");
  fprintf(stdout, "<Folder>\n");
  fprintf(stdout, "<name>Layer:no_name</name>\n");
  fprintf(stdout, "<open>0</open>\n");
  fprintf(stdout, " <styleUrl>#check-hide-children</styleUrl>\n");
  fprintf(stdout, "<TimeStamp>\n");
  fprintf(stdout, " <when>2021-11-04T13:19:00Z</when>\n");
  fprintf(stdout, "</TimeStamp>\n");
  fprintf(stdout, "<styleUrl>#hiker-icon</styleUrl>\n");
  fprintf(stdout, "<description><![CDATA[Layer:no_name]]></description>\n");

  int height = 5000;
  for (size_t i = 0; i < ngp; ++i)
    {
      if (grid_mask.size() && grid_mask[i] == 0) continue;

      auto lonBounds = &grid_corner_lon[i * ncorner];
      auto latBounds = &grid_corner_lat[i * ncorner];
      auto ncornerNew = check_ncorner(ncorner, lonBounds, latBounds);
      check_lonbounds(ncornerNew, lonBounds);

      auto rgb = get_rgb(array[i], missval, cpt);

      fprintf(stdout, "<Placemark>\n");
      fprintf(stdout, "<visibility>1</visibility>\n");
      fprintf(stdout, "<open>0</open>\n");
      fprintf(stdout, "<Style>\n");
      fprintf(stdout, "<PolyStyle>\n");
      fprintf(stdout, "	<!-- r:%d g:%d b:%d -->\n", rgb[0], rgb[1], rgb[2]);
      fprintf(stdout, "	<color>fe%x%x%x</color>\n", (unsigned char) rgb[0], (unsigned char) rgb[1], (unsigned char) rgb[2]);
      fprintf(stdout, "	<fill>1</fill>\n");
      fprintf(stdout, "</PolyStyle>\n");
      fprintf(stdout, "<LineStyle>\n");
      fprintf(stdout, "	<width>2</width>\n");
      fprintf(stdout, "	<!-- r:%d g:%d b:%d -->\n", rgb[0], rgb[1], rgb[2]);
      fprintf(stdout, "	<color>fe%x%x%x</color>\n", (unsigned char) rgb[0], (unsigned char) rgb[1], (unsigned char) rgb[2]);
      fprintf(stdout, "</LineStyle>\n");
      fprintf(stdout, "</Style>\n");
      fprintf(stdout, "<MultiGeometry>\n");

      fprintf(stdout, "<Polygon>\n");
      fprintf(stdout, "<extrude>1</extrude>\n");
      fprintf(stdout, "<altitudeMode>clampToGround</altitudeMode>\n");
      fprintf(stdout, "<tessellate>0</tessellate>\n");
      fprintf(stdout, " <outerBoundaryIs>\n");
      fprintf(stdout, "  <LinearRing>\n");
      fprintf(stdout, "   <coordinates>\n");

      for (int ic = 0; ic < ncornerNew; ++ic) fprintf(stdout, "     %g,%g,%d\n", lonBounds[ic], latBounds[ic], height);
      fprintf(stdout, "     %g,%g,%d\n", lonBounds[0], latBounds[0], height);

      fprintf(stdout, "   </coordinates>\n");
      fprintf(stdout, "  </LinearRing>\n");
      fprintf(stdout, " </outerBoundaryIs>\n");
      fprintf(stdout, "</Polygon>\n");
      fprintf(stdout, "</MultiGeometry>\n");
      fprintf(stdout, "</Placemark>\n");
    }

  fprintf(stdout, "</Folder>\n");
  fprintf(stdout, "<ScreenOverlay id=\"legend\">\n");
  fprintf(stdout, "<name>Legend</name>\n");
  fprintf(stdout, "<Icon>\n");
  fprintf(stdout, " <href>legend.png</href>\n");
  fprintf(stdout, "</Icon>\n");
  fprintf(stdout, "<overlayXY x=\"0\" y=\"0\" xunits=\"fraction\" yunits=\"fraction\"/>\n");
  fprintf(stdout, "<screenXY x=\"0\" y=\"0\" xunits=\"fraction\" yunits=\"fraction\"/>\n");
  fprintf(stdout, "<size x=\"-1\" y=\"0.1\" xunits=\"fraction\" yunits=\"fraction\"/>\n");
  fprintf(stdout, "</ScreenOverlay>\n");
  fprintf(stdout, "</Folder>\n");
  fprintf(stdout, "</Document>\n");
  fprintf(stdout, "</kml>\n");
}

static void
output_vector(long nlon, long nlat, int ninc, const Varray<double> &lon, const Varray<double> &lat, const Varray<double> &uf,
              const Varray<double> &vf)
{
  for (long j = 0; j < nlat; j += ninc)
    for (long i = 0; i < nlon; i += ninc)
      {
        auto idx = j * nlon + i;

        // compute length of velocity vector
        auto auv = std::sqrt(uf[idx] * uf[idx] + vf[idx] * vf[idx]);

        auto alpha = std::atan2(vf[idx], uf[idx]);
        alpha = 90. - alpha * RAD2DEG;

        if (alpha < 0) alpha += 360;
        if (alpha > 360) alpha -= 360;

        if (std::fabs(auv) > 0) fprintf(stdout, " %g  %g  %g  %g\n", lon[idx], lat[idx], alpha, auv);
      }

  fprintf(stdout, "#\n");
}

void *
Outputgmt(void *process)
{
  int varID0;
  size_t gridsize2 = 0;
  size_t nmiss;
  int ninc = 1;
  bool lzon = false, lmer = false, lhov = false;
  Varray<double> grid_center_lat2, grid_center_lon2;
  Varray<double> grid_corner_lat, grid_corner_lon;
  std::vector<int> grid_mask;
  CPT cpt;

  cdo_initialize(process);

  // clang-format off
  auto OUTPUTCENTER    = cdo_operator_add("gmtxyz",          0, 0, nullptr);
  auto OUTPUTCENTER2   = cdo_operator_add("outputcenter2",   0, 0, nullptr);
  auto OUTPUTCENTERCPT = cdo_operator_add("outputcentercpt", 0, 0, nullptr);
  auto OUTPUTBOUNDS    = cdo_operator_add("gmtcells",        0, 0, nullptr);
  auto OUTPUTBOUNDSCPT = cdo_operator_add("outputboundscpt", 0, 0, nullptr);
  auto OUTPUTVECTOR    = cdo_operator_add("outputvector",    0, 0, nullptr);
  auto OUTPUTTRI       = cdo_operator_add("outputtri",       0, 0, nullptr);
  auto OUTPUTVRML      = cdo_operator_add("outputvrml",      0, 0, nullptr);
  auto OUTPUTKML       = cdo_operator_add("outputkml",       0, 0, nullptr);
  // clang-format on

  auto operatorID = cdo_operator_id();

  auto printHeader = (operatorID != OUTPUTTRI && operatorID != OUTPUTKML);

  if (operatorID == OUTPUTVECTOR)
    {
      operator_input_arg("increment");
      operator_check_argc(1);
      ninc = parameter_to_int(cdo_operator_argv(0));
      if (ninc < 1) cdo_abort("Increment must be greater than 0!");
    }

  auto needCellCorners = (operatorID == OUTPUTBOUNDS || operatorID == OUTPUTBOUNDSCPT || operatorID == OUTPUTKML);

  if (operatorID == OUTPUTCENTERCPT || operatorID == OUTPUTBOUNDSCPT || operatorID == OUTPUTKML || operatorID == OUTPUTVRML)
    {
      operator_check_argc(1);
      auto cpt_file = cdo_operator_argv(0).c_str();

      auto cpt_fp = std::fopen(cpt_file, "r");
      if (cpt_fp == nullptr) cdo_abort("Open failed on color palette table %s", cpt_file);

      auto status = cpt_read(cpt_fp, &cpt);
      if (status != 0) cdo_abort("Error during read of color palette table %s", cpt_file);

      if (Options::cdoVerbose) cpt_write(stderr, cpt);
    }

  auto streamID = cdo_open_read(0);

  auto vlistID = cdo_stream_inq_vlist(streamID);
  auto taxisID = vlistInqTaxis(vlistID);

  int varID = 0;
  VarList varList;
  varListInit(varList, vlistID);
  const auto &var = varList[varID];

  auto gridID = var.gridID;
  auto zaxisID = var.zaxisID;
  auto missval = var.missval;

  gridID = generate_full_cell_grid(gridID);

  if (!gridHasCoordinates(gridID)) cdo_abort("Cell center coordinates missing!");

  auto gridsize = gridInqSize(gridID);

  long nlon = gridInqXsize(gridID);
  long nlat = gridInqYsize(gridID);
  long nlev = zaxisInqSize(zaxisID);

  if (gridInqMaskGME(gridID, nullptr))
    {
      grid_mask.resize(gridsize);
      gridInqMaskGME(gridID, grid_mask.data());
    }

  if (gridInqType(gridID) != GRID_UNSTRUCTURED)
    {
      if (nlon == 1 && nlat > 1 && nlev == 1) lhov = true;
      if (nlon == 1 && nlat > 1 && nlev > 1) lzon = true;
      if (nlon > 1 && nlat == 1 && nlev > 1) lmer = true;
    }
  else { nlat = 1; }

  if (Options::cdoVerbose && lhov) cdo_print("Process hovmoeller data");
  if (Options::cdoVerbose && lzon) cdo_print("Process zonal data");
  if (Options::cdoVerbose && lmer) cdo_print("Process meridional data");
  /*
  if ( lzon || lmer )
    {
      if ( operatorID == OUTPUTBOUNDS || operatorID == OUTPUTBOUNDSCPT )
        cdo_abort("Bounds not available for zonal/meridional data!");
    }
  */
  if (lhov)
    {
      if (operatorID == OUTPUTBOUNDS || operatorID == OUTPUTBOUNDSCPT) cdo_abort("Bounds not available hovmoeller data!");
    }

  const int ncorner = (gridInqType(gridID) == GRID_UNSTRUCTURED) ? gridInqNvertex(gridID) : 4;

  const bool grid_is_circular = gridIsCircular(gridID);

  Varray<double> grid_center_lat(gridsize), grid_center_lon(gridsize);
  gridInqYvals(gridID, grid_center_lat.data());
  gridInqXvals(gridID, grid_center_lon.data());

  // Convert lat/lon units if required
  cdo_grid_to_degree(gridID, CDI_XAXIS, gridsize, grid_center_lon.data(), "grid center lon");
  cdo_grid_to_degree(gridID, CDI_YAXIS, gridsize, grid_center_lat.data(), "grid center lat");

  long nvals = gridsize;
  double *plon = grid_center_lon.data();
  double *plat = grid_center_lat.data();

  if (operatorID == OUTPUTCENTER2 && grid_is_circular)
    {
      gridsize2 = nlat * (nlon + 1);

      grid_center_lat2.resize(gridsize2);
      grid_center_lon2.resize(gridsize2);

      make_cyclic(grid_center_lat.data(), grid_center_lat2.data(), nlon, nlat);
      make_cyclic(grid_center_lon.data(), grid_center_lon2.data(), nlon, nlat);

      for (long j = 0; j < nlat; ++j)
        {
          const long ij2 = j * (nlon + 1);
          grid_center_lon2[ij2 + nlon] += 360;
        }

      nvals = gridsize2;
      plon = grid_center_lon2.data();
      plat = grid_center_lat2.data();
    }

  Varray<double> zaxis_center_lev(nlev), zaxis_lower_lev(nlev), zaxis_upper_lev(nlev);

  cdo_zaxis_inq_levels(zaxisID, zaxis_center_lev.data());

  if (needCellCorners)
    {
      if (ncorner == 0) cdo_abort("Number of cell corners undefined!");
      const size_t nalloc = ncorner * gridsize;
      grid_corner_lat.resize(nalloc);
      grid_corner_lon.resize(nalloc);

      if (!gridHasBounds(gridID)) cdo_abort("Cell corner coordinates missing!");

      gridInqYbounds(gridID, grid_corner_lat.data());
      gridInqXbounds(gridID, grid_corner_lon.data());

      cdo_grid_to_degree(gridID, CDI_XAXIS, ncorner * gridsize, grid_corner_lon.data(), "grid corner lon");
      cdo_grid_to_degree(gridID, CDI_YAXIS, ncorner * gridsize, grid_corner_lat.data(), "grid corner lat");

      if (zaxisInqLbounds(zaxisID, nullptr) && zaxisInqUbounds(zaxisID, nullptr))
        {
          zaxisInqLbounds(zaxisID, zaxis_lower_lev.data());
          zaxisInqUbounds(zaxisID, zaxis_upper_lev.data());
        }
      else
        {
          zaxis_lower_lev[0] = zaxis_center_lev[0];
          for (int i = 1; i < nlev; ++i) zaxis_lower_lev[i] = 0.5 * (zaxis_center_lev[i] + zaxis_center_lev[i - 1]);

          zaxis_upper_lev[nlev - 1] = zaxis_center_lev[nlev - 1];
          for (int i = 0; i < nlev - 1; ++i) zaxis_upper_lev[i] = zaxis_lower_lev[i + 1];

          if (Options::cdoVerbose)
            for (int i = 0; i < nlev; ++i)
              fprintf(stderr, "level: %d %g %g %g\n", i + 1, zaxis_lower_lev[i], zaxis_center_lev[i], zaxis_upper_lev[i]);
        }
    }

  Varray<double> array(gridsize);
  double *parray = array.data();

  Varray<double> array2;
  if (operatorID == OUTPUTCENTER2 && grid_is_circular)
    {
      array2.resize(nlat * (nlon + 1));
      parray = array2.data();
    }

  Varray<double> uf, vf;
  if (operatorID == OUTPUTVECTOR) uf.resize(gridsize);
  if (operatorID == OUTPUTVECTOR) vf.resize(gridsize);

  int tsID = 0;
  while (true)
    {
      auto nrecs = cdo_stream_inq_timestep(streamID, tsID);
      if (nrecs == 0) break;

      auto vDateTime = taxisInqVdatetime(taxisID);
      auto vdateString = date_to_string(vDateTime.date);
      auto vtimeString = time_to_string(vDateTime.time);

      if (tsID == 0 && printHeader)
        {
          if (operatorID == OUTPUTVRML) printf("#VRML V2.0 utf8\n\n");
#ifdef VERSION
          fprintf(stdout, "# Generated by CDO version %s\n", VERSION);
          fprintf(stdout, "#\n");
#endif
          fprintf(stdout, "# Operator = %s\n", cdo_operator_name(operatorID));
          // clang-format off
          if      (lhov) fprintf(stdout, "# Mode     = hovmoeller\n");
          else if (lzon) fprintf(stdout, "# Mode     = zonal\n");
          else if (lmer) fprintf(stdout, "# Mode     = meridional\n");
          else           fprintf(stdout, "# Mode     = horizonal\n");
          // clang-format on

          if (operatorID == OUTPUTVECTOR) fprintf(stdout, "# Increment = %d\n", ninc);
          fprintf(stdout, "#\n");
          fprintf(stdout, "# Stream = %s\n", cdo_get_stream_name(0));
          fprintf(stdout, "# Date   = %s\n", vdateString.c_str());
          fprintf(stdout, "# Time   = %s\n", vtimeString.c_str());
          fprintf(stdout, "# Name   = %s\n", var.name.c_str());
          fprintf(stdout, "# Code   = %d\n", var.code);
        }

      varID0 = varID;

      for (int recID = 0; recID < nrecs; ++recID)
        {
          int levelID;
          cdo_inq_record(streamID, &varID, &levelID);

          if (varID != varID0) continue;
          if (recID > 0 && !lzon && !lmer) continue;

          cdo_read_record(streamID, array.data(), &nmiss);

          if (operatorID == OUTPUTCENTER2 && grid_is_circular) make_cyclic(array.data(), array2.data(), nlon, nlat);

          auto level = zaxis_center_lev[levelID];

          if ((tsID == 0 || lzon || lmer) && printHeader) fprintf(stdout, "# Level  = %g\n", level);
          if (lhov) fprintf(stdout, "# Timestep = %d\n", tsID + 1);

          if (printHeader) fprintf(stdout, "#\n");

          if (operatorID == OUTPUTCENTER || operatorID == OUTPUTCENTER2 || operatorID == OUTPUTCENTERCPT)
            {
              if (Options::cdoVerbose)
                {
                  auto mm = varray_min_max_mv(gridsize, array, missval);
                  auto range = mm.max - mm.min;
                  fprintf(stderr, "makecpt -T%g/%g/%g -Crainbow > gmt.cpt\n", mm.min, mm.max, range / 20);
                  fprintf(stderr, "pscontour -K -JQ0/10i -Rd -I -Cgmt.cpt data.gmt > gmtplot.ps\n");
                  fprintf(stderr, "pscoast -O -J -R -Dc -W -B40g20 >> gmtplot.ps\n");
                }

              for (long i = 0; i < nvals; ++i)
                {
                  if (grid_mask.size() && grid_mask[i] == 0) continue;

                  if (operatorID == OUTPUTCENTER)
                    {
                      if (lzon)
                        fprintf(stdout, " %g  %g  %g\n", grid_center_lat[i], level, array[i]);
                      else if (lmer)
                        fprintf(stdout, " %g  %g  %g\n", grid_center_lon[i], level, array[i]);
                      else if (lhov)
                        fprintf(stdout, " %d  %g  %g\n", tsID + 1, grid_center_lat[i], array[i]);
                      else
                        fprintf(stdout, " %g  %g  %g\n", grid_center_lon[i], grid_center_lat[i], array[i]);
                    }
                  else if (operatorID == OUTPUTCENTER2) { fprintf(stdout, " %g  %g  %g\n", plon[i], plat[i], parray[i]); }
                  else
                    {
                      if (lzon)
                        fprintf(stdout, " %g  %g  %g\n", grid_center_lat[i], level, array[i]);
                      else if (lmer)
                        fprintf(stdout, " %g  %g  %g\n", grid_center_lon[i], level, array[i]);
                      else
                        fprintf(stdout, " %g  %g  %g  %g\n", grid_center_lon[i], grid_center_lat[i], array[i], array[i]);
                    }
                }
              fprintf(stdout, "#\n");
            }
          else if (operatorID == OUTPUTTRI)
            {
              if (gridInqType(gridID) != GRID_CURVILINEAR) cdo_abort("Unsupported grid!");

              const long mlon = nlon - 1;
              // if ( gridIsCircular(gridID) ) mlon = nlon;
              for (long j = 0; j < nlat - 1; ++j)
                for (long i = 0; i < mlon; ++i)
                  {
                    int ip1 = i + 1;
                    if (i == nlon - 1) ip1 = 0;
                    int c1 = (j) *nlon + ip1;
                    int c2 = (j) *nlon + i;
                    int c3 = (j + 1) * nlon + i;
                    fprintf(stdout, "%d   %d   %d\n", c1, c2, c3);
                    c1 = (j) *nlon + i + 1;
                    c2 = (j + 1) * nlon + i;
                    c3 = (j + 1) * nlon + ip1;
                    fprintf(stdout, "%d   %d   %d\n", c1, c2, c3);
                  }
            }
          else if (operatorID == OUTPUTKML)
            {
              output_kml(gridsize, array, ncorner, grid_corner_lat, grid_corner_lon, grid_mask, missval, cpt);
            }
          else if (operatorID == OUTPUTVECTOR)
            {
              if (nrecs < 2) cdo_abort("Too few fields!");

              varray_copy(gridsize, array, uf);
              cdo_inq_record(streamID, &varID, &levelID);
              cdo_read_record(streamID, vf.data(), &nmiss);

              output_vector(nlon, nlat, ninc, grid_center_lon, grid_center_lat, uf, vf);

              break;
            }
          else if (operatorID == OUTPUTVRML) { output_vrml(nlon, nlat, gridsize, array, missval, cpt); }
          else if (operatorID == OUTPUTBOUNDS || operatorID == OUTPUTBOUNDSCPT)
            {
              if (Options::cdoVerbose)
                {
                  auto mm = varray_min_max_mv(gridsize, array, missval);
                  auto range = mm.max - mm.min;
                  fprintf(stderr, "makecpt -T%g/%g/%g -Crainbow > gmt.cpt\n", mm.min, mm.max, range / 20);
                  fprintf(stderr, "psxy -K -JQ0/10i -Rd -L -Cgmt.cpt -m data.gmt > gmtplot.ps\n");
                  // fprintf(stderr, "psxy -K -Jx0.028id -Rd -L -Cgmt.cpt -m
                  // data.gmt > gmtplot.ps\n"); fprintf(stderr, "psxy -K
                  // -JN0/10i -Rd -L -Cgmt.cpt -m data.gmt > gmtplot.ps\n");
                  fprintf(stderr, "pscoast -O -J -R -Dc -W -B40g20 >> gmtplot.ps\n");
                  fprintf(stderr, "ps2pdf gmtplot.ps\n");
                }

              for (size_t i = 0; i < gridsize; ++i)
                {
                  if (grid_mask.size() && grid_mask[i] == 0) continue;

                  if (!DBL_IS_EQUAL(array[i], missval))
                    fprintf(stdout, "> -Z%g", array[i]);
                  else
                    fprintf(stdout, "> -ZNaN");

                  if (operatorID == OUTPUTBOUNDSCPT)
                    {
                      auto rgb = get_rgb(array[i], missval, cpt);
                      fprintf(stdout, " -G%d/%d/%d", rgb[0], rgb[1], rgb[2]);
                    }

                  fprintf(stdout, "\n");

                  if (lzon) { output_zon(zaxis_lower_lev[levelID], zaxis_upper_lev[levelID], &grid_corner_lat[i * 4]); }
                  else if (lmer) { output_mer(zaxis_lower_lev[levelID], zaxis_upper_lev[levelID], &grid_corner_lon[i * 4]); }
                  else if (lhov) { cdo_abort("Implementation for hovmoeller data missing!"); }
                  else
                    {
                      const double *lon_bounds = grid_corner_lon.data() + i * ncorner;
                      const double *lat_bounds = grid_corner_lat.data() + i * ncorner;
                      const int ncorner_new = check_ncorner(ncorner, lon_bounds, lat_bounds);

                      for (int ic = 0; ic < ncorner_new; ++ic) fprintf(stdout, "   %g  %g\n", lon_bounds[ic], lat_bounds[ic]);
                      fprintf(stdout, "   %g  %g\n", lon_bounds[0], lat_bounds[0]);
                    }
                }
              fprintf(stdout, "\n");
            }
        }

      if (!lhov) break;

      tsID++;
    }

  cdo_stream_close(streamID);

  cdo_finish();

  return nullptr;
}
