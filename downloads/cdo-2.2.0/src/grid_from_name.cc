/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include "process_int.h"
#include "mpim_grid.h"
#include "grid_healpix.h"
#include "gaussian_latitudes.h"
#include "griddes.h"
#include "util_string.h"
#include "dcw_reader.h"

size_t genIcosphereCoords(int subdivisions, bool withBounds, std::vector<double> &xvals, std::vector<double> &yvals,
                          std::vector<double> &xbounds, std::vector<double> &ybounds);

static void
gen_grid_icosphere(GridDesciption &grid, const char *pline)
{
  int gridtype = GRID_UNSTRUCTURED;
  bool withBounds = true;
  long b = 0;

  if (*pline != 0)
    {
      if (*pline == 'r')
        pline++;
      else
        return;

      if (*pline == 0) return;
      if (!isdigit((int) *pline)) return;

      char *endptr = (char *) pline;
      auto r = strtol(pline, &endptr, 10);
      if (*endptr == 0 || r != 2) return;
      pline = endptr;

      if (*pline == 'b')
        pline++;
      else
        return;

      if (*pline == 0) return;
      if (!isdigit((int) *pline)) return;

      endptr = (char *) pline;
      b = strtol(pline, &endptr, 10);

      if (*endptr != 0)
        {
          pline = endptr;
          if (*pline != '_') return;
          pline++;
          if (*pline == 0) return;
          if (*pline == '0')
            {
              withBounds = false;
              pline++;
            }
          if (*pline != 0) return;
        }
    }

  grid.type = gridtype;
  if (withBounds) grid.nvertex = 3;

  auto ncells = genIcosphereCoords(b + 1, withBounds, grid.xvals, grid.yvals, grid.xbounds, grid.ybounds);
  grid.xsize = ncells;
  grid.ysize = ncells;
  grid.xname = "clon";
  grid.yname = "clat";
  grid.xunits = "radian";
  grid.yunits = "radian";
}

static void
gen_grid_zonal(GridDesciption &grid, const char *pline, double inc, double lon1, double lon2, double lat1, double lat2)
{
  int gridtype = GRID_LONLAT;
  bool withBounds = true;

  if (*pline != 0)
    {
      if (*pline == '_')
        pline++;
      else
        return;

      if (*pline == 0) return;

      if (!isdigit((int) *pline) && !ispunct((int) *pline)) return;

      char *endptr = (char *) pline;
      inc = strtod(pline, &endptr);
      if (*endptr != 0)
        {
          pline = endptr;
          if (*pline == '_')
            pline++;
          else
            return;

          if (*pline == 0) return;
          if (*pline != 0) return;
        }

      if (inc < 1.0e-9) inc = 1.0;
      if (inc > 180.0) cdo_abort("Zonal latitude increment out if range (max=180)!");
    }

  grid.type = gridtype;

  if (lon1 >= lon2 || lat1 >= lat2) cdo_abort("Invalid grid box: lon1=%g lon2=%g lat1=%g lat2=%g", lon1, lon2, lat1, lat2);

  auto nlon = 1;
  auto nlat = (size_t) ((lat2 - lat1) / inc + 0.5);

  grid.xvals.resize(nlon, 0.0);
  grid.yvals.resize(nlat);

  for (size_t i = 0; i < nlat; ++i) grid.yvals[i] = lat1 + inc / 2.0 + i * inc;

  grid.xsize = nlon;
  grid.ysize = nlat;

  if (withBounds)
    {
      grid.xbounds.resize(2);
      grid.xbounds[0] = lon1;
      grid.xbounds[1] = lon2;
      grid.ybounds.resize(2 * nlat);
      grid.ybounds[0] = lat1;
      grid.ybounds[1] = lat2;
      if (nlat > 1) grid_gen_bounds(nlat, grid.yvals, grid.ybounds);
    }
}

static void
gen_grid_lonlat(GridDesciption &grid, const char *pline, double inc, double lon1, double lon2, double lat1, double lat2)
{
  int gridtype = GRID_LONLAT;
  bool withBounds = true;

  if (*pline != 0 && (*pline == '+' || *pline == '-') && (isdigit((int) *(pline + 1)) || ispunct((int) *(pline + 1))))
    {
      char *endptr = (char *) pline;
      auto off = strtod(pline, &endptr);
      pline = endptr;

      lon1 -= off;
      lon2 += off;
      lat1 -= off;
      lat2 += off;
      lat1 = std::max(lat1, -90.0);
      lat2 = std::min(lat2, 90.0);
    }

  if (*pline != 0)
    {
      if (*pline == '_')
        pline++;
      else
        return;

      if (*pline == 0) return;

      if (!isdigit((int) *pline) && !ispunct((int) *pline)) return;

      char *endptr = (char *) pline;
      inc = strtod(pline, &endptr);
      if (*endptr != 0)
        {
          pline = endptr;
          if (*pline == '_')
            pline++;
          else
            return;

          if (*pline == 0) return;
          if (*pline == 'c')
            {
              gridtype = GRID_CURVILINEAR;
              pline++;
              if (*pline == '0')
                {
                  withBounds = false;
                  pline++;
                }
            }
          else if (*pline == 'u')
            {
              gridtype = GRID_UNSTRUCTURED;
              pline++;
              if (*pline == '0')
                {
                  withBounds = false;
                  pline++;
                }
            }
          if (*pline != 0) return;
        }

      if (inc < 1.0e-9) inc = 1.0;
      if (inc > 180.0) cdo_abort("LonLat increment out if range (max=180)!");
    }

  grid.type = gridtype;

  if (lon1 >= lon2 || lat1 >= lat2) cdo_abort("Invalid grid box: lon1=%g lon2=%g lat1=%g lat2=%g", lon1, lon2, lat1, lat2);

  auto nlon = (size_t) ((lon2 - lon1) / inc + 0.5);
  auto nlat = (size_t) ((lat2 - lat1) / inc + 0.5);

  grid.xvals.resize(nlon);
  grid.yvals.resize(nlat);

  for (size_t i = 0; i < nlon; ++i) grid.xvals[i] = lon1 + inc * 0.5 + i * inc;
  for (size_t i = 0; i < nlat; ++i) grid.yvals[i] = lat1 + inc * 0.5 + i * inc;

  if (gridtype == GRID_LONLAT)
    {
      grid.xsize = nlon;
      grid.ysize = nlat;
    }
  else
    {
      std::vector<double> yvals(nlat);
      for (size_t j = 0; j < nlat; ++j) yvals[j] = grid.yvals[j];
      auto gridsize = nlon * nlat;
      grid.xvals.resize(gridsize);
      grid.yvals.resize(gridsize);
      for (size_t j = 0; j < nlat; ++j)
        for (size_t i = 0; i < nlon; ++i)
          {
            grid.xvals[j * nlon + i] = grid.xvals[i];
            grid.yvals[j * nlon + i] = yvals[j];
          }

      if (gridtype == GRID_CURVILINEAR)
        {
          grid.xsize = nlon;
          grid.ysize = nlat;
        }
      else
        {
          grid.xsize = gridsize;
          grid.ysize = gridsize;
          if (withBounds) grid.nvertex = 4;
        }

      if (withBounds && nlon > 1 && nlat > 1)
        {
          std::vector<double> xbounds(2 * nlon), ybounds(2 * nlat);

          grid_gen_bounds(nlon, grid.xvals, xbounds);
          grid_gen_bounds(nlat, yvals, ybounds);
          grid_check_lat_borders(2 * nlat, ybounds.data());

          grid.xbounds.resize(4 * gridsize);
          grid.ybounds.resize(4 * gridsize);
          grid_gen_xbounds2D(nlon, nlat, xbounds, grid.xbounds);
          grid_gen_ybounds2D(nlon, nlat, ybounds, grid.ybounds);
        }
    }
}

static void
gen_grid_dcw(GridDesciption &grid, const char *pline, double inc)
{
  auto param1 = pline;
  auto param2 = strstr(pline, "_");
  auto param1len = param2 ? param2 - pline : strlen(pline);

  if (param2)
    {
      pline = param2 + 1;

      if (*pline == 0) return;
      if (!isdigit((int) *pline) && !ispunct((int) *pline)) return;

      char *endptr = (char *) pline;
      inc = strtod(pline, &endptr);
    }

  const std::string codeNames(string_to_upper({ param1, param1len }));

  DCW_Lists dcw_lists;
  if (dcw_load_lists(dcw_lists)) cdo_abort("dcw_load_lists failed!");

  auto codeList = split_string(codeNames, "\\+");

  dcw_sort_countries(dcw_lists);

  codeList = dcw_expand_code_list(dcw_lists, codeList);

  Region region;
  if (dcw_get_region(dcw_lists, codeList, region)) cdo_abort("dcw_get_region failed!");

  // printf("lon1, lon2, lat1, lat2 %g %g %g %g\n", region.west, region.east, region.south, region.north);
  auto lon1 = std::round(region.west / inc - 0.5) * inc;
  auto lon2 = std::round(region.east / inc + 0.5) * inc;
  auto lat1 = std::round(region.south / inc - 0.5) * inc;
  auto lat2 = std::round(region.north / inc + 0.5) * inc;
  // printf("lon1, lon2, lat1, lat2 %g %g %g %g\n", lon1, lon2, lat1, lat2);

  const char *param = param2 ? param2 : "";
  gen_grid_lonlat(grid, param, inc, lon1 - inc * 0.5, lon2 + inc * 0.5, lat1 - inc * 0.5, lat2 + inc * 0.5);
}

static void
gen_grid_gme(GridDesciption &grid, const char *pline)
{
  if (isdigit((int) *pline))
    {
      char *endptr = (char *) pline;
      auto ni = strtol(pline, &endptr, 10);
      if (*endptr == 0)
        {
          grid.type = GRID_GME;
          grid.ni = ni;
          grid.nd = 10;
          gme_factorni(grid.ni, &grid.ni2, &grid.ni3);
          grid.size = (grid.ni + 1) * (grid.ni + 1) * 10;
        }
    }
}

static void
gen_grid_healpix(GridDesciption &grid, const char *pline)
{
  if (isdigit((int) *pline))
    {
      char *endptr = (char *) pline;
      size_t nside = strtol(pline, &endptr, 10);
      if (*endptr == 0 || endptr != pline)
        {
          HpOrder hpOrder(HpOrder::Nested);

          pline = endptr;
          /*
          bool withBounds = false;
          if (*pline == 'b')
            {
              pline++;
              withBounds = true;
            }
          */
          if (*pline == '_')
            {
              pline++;
              hpOrder = hp_get_order(pline);
              if (hpOrder == HpOrder::Undef || hpOrder == HpOrder::XY) return;
            }

          size_t ncells = 12 * nside * nside;
          grid.type = GRID_PROJECTION;
          grid.size = ncells;
          grid.projection = "healpix";
          grid.healpixNside = nside;
          grid.healpixOrder = (hpOrder == HpOrder::Ring) ? "ring" : "nested";

          /*
          grid.type = GRID_UNSTRUCTURED;
          grid.size = ncells;
          grid.xunits = "radian";
          grid.yunits = "radian";
          grid.xvals.resize(ncells);
          grid.yvals.resize(ncells);

          if (withBounds)
            {
              grid.nvertex = 4;
              grid.xbounds.resize(4 * ncells);
              grid.ybounds.resize(4 * ncells);
            }

          hp_generate_coords(hpOrder, nside, ncells, grid.xvals.data(), grid.yvals.data(), withBounds, grid.xbounds.data(),
                             grid.ybounds.data());
          */
        }
    }
}

void
gaussian_latitudes_in_degrees(std::vector<double> &lats, std::vector<double> &lat_bounds, size_t nlat)
{
  // lats(nlat)
  // lat_bounds(nlat+1)
  std::vector<double> latw(nlat), latw_cumsum(nlat);

  gaussian_latitudes(nlat, lats.data(), latw.data());

  for (size_t j = 0; j < nlat; ++j) lats[j] = RAD2DEG * std::asin(lats[j]);

  latw_cumsum[0] = latw[0];
  for (size_t j = 1; j < nlat; ++j) latw_cumsum[j] = latw_cumsum[j - 1] + latw[j];

  lat_bounds[0] = 1.0;
  for (size_t j = 1; j < nlat; ++j) lat_bounds[j] = 1.0 - latw_cumsum[j - 1];
  lat_bounds[nlat] = -1.0;

  for (size_t j = 0; j < nlat + 1; ++j) lat_bounds[j] = RAD2DEG * std::asin(lat_bounds[j]);
}

static void
gen_grid_gea(GridDesciption &grid, const char *pline)
{
  if (isdigit((int) *pline))
    {
      auto endptr = (char *) pline;
      auto dx = strtod(pline, &endptr);
      if (*endptr != 0) return;

      auto dy = dx;
      constexpr auto re = 6378.137;
      constexpr auto f = 1.0 / 298.257223563;
      constexpr auto rp = re * (1.0 - f);
      constexpr auto polar_circumference = 2.0 * M_PI * rp;
      constexpr auto equator_circumference = 2.0 * M_PI * re;

      size_t nlat = 0.5 * polar_circumference / dy;
      if (nlat % 2) nlat++;

      std::vector<double> lats(nlat), lat_bounds(nlat + 1);
      gaussian_latitudes_in_degrees(lats, lat_bounds, nlat);

      std::vector<double> cell_height(nlat);
      for (size_t j = 0; j < nlat; ++j) cell_height[j] = 0.25 * polar_circumference * (lat_bounds[j] - lat_bounds[j + 1]) / 90.0;

      size_t nlone = equator_circumference / dx;
      if (nlone % 2) nlone++;

      std::vector<int> reducedPoints(nlat);
      size_t ncells = 0;
      for (size_t j = 0; j < nlat; ++j)
        {
          auto rlat = re * std::cos(DEG2RAD * lats[j]);
          auto circumference = 2.0 * M_PI * rlat;
          auto dx_to_use = dx * dy / cell_height[j];

          size_t nlon = std::max((int) std::lround(circumference / dx_to_use), 1);
          if (nlon % 2) nlon++;

          reducedPoints[j] = nlon;
          ncells += nlon;
        }

      // printf("%zu %zu %zu %zu %g\n", ncells, nlone, nlat, nlone*nlat, 100.0*ncells/(nlone*nlat));

      std::vector<double> lons(ncells);
      size_t ij = 0;
      for (size_t j = 0; j < nlat; ++j)
        {
          size_t nlon = reducedPoints[j];
          for (size_t i = 0; i < nlon; ++i) lons[ij++] = i * 360. / nlon;
        }

      grid.type = GRID_GAUSSIAN_REDUCED;
      grid.size = ncells;
      grid.xsize = ncells;
      grid.ysize = nlat;
      grid.numLPE = nlat / 2;
      grid.xvals.resize(ncells);
      grid.yvals.resize(nlat);
      grid.ybounds.resize(nlat * 2);
      grid.reducedPoints.resize(nlat);
      for (size_t i = 0; i < ncells; ++i) grid.xvals[i] = lons[i];
      for (size_t j = 0; j < nlat; ++j) grid.yvals[j] = lats[j];
      for (size_t j = 0; j < nlat; ++j) grid.ybounds[j * 2 + 1] = lat_bounds[j];
      for (size_t j = 0; j < nlat; ++j) grid.ybounds[j * 2] = lat_bounds[j + 1];
      for (size_t j = 0; j < nlat; ++j) grid.reducedPoints[j] = reducedPoints[j];
    }
}

static void
gen_grid_zonal(GridDesciption &grid, const char *pline)
{
  if (isdigit((int) *pline))
    {
      constexpr size_t nextra = 1;
      grid.type = GRID_UNSTRUCTURED;
      auto nlats = (size_t) atol(pline);
      auto nlons = nlats * 2;
      auto gridsize = nlats;
      grid.size = gridsize;
      grid.xsize = gridsize;
      grid.ysize = gridsize;
      grid.xvals.resize(gridsize);
      grid.yvals.resize(gridsize);
      for (size_t i = 0; i < nlats; ++i) grid.xvals[i] = 180.0;
      auto dlat = 180.0 / nlats;
      auto dlon = 360.0 / nlons;
      // printf("dlat %g dlon %g\n", dlat, dlon);
      for (size_t i = 0; i < nlats; ++i) grid.yvals[i] = -90.0 + i * dlat + dlat / 2.0;
      auto nv = (nlons + 1) * 2;
      grid.nvertex = nv;
      grid.xbounds.resize(nv * gridsize);
      grid.ybounds.resize(nv * gridsize);
      std::vector<double> xbounds(nlons + 1), ybounds(nlats + 1);
      for (size_t i = 0; i <= nlons; ++i) xbounds[i] = 0.0 + i * dlon;
      for (size_t i = 0; i <= nlats; ++i) ybounds[i] = -90.0 + i * dlat;
      // for (size_t i = 0; i <= nlons; ++i)  printf("lon %zu %g\n", i, xbounds[i]);
      // for (size_t i = 0; i <= nlats; ++i)  printf("lat %zu %g\n", i, ybounds[i]);
      size_t k = 0;
      for (size_t j = 0; j < nlats; ++j)
        {
          for (size_t i = nlons; i > 0; i--)
            {
              grid.xbounds[k] = xbounds[i];
              grid.ybounds[k] = ybounds[j + 1];
              k++;
            }
          for (size_t i = 0; i < nextra; ++i)
            {
              grid.xbounds[k] = xbounds[0];
              grid.ybounds[k] = ybounds[j + 1];
              k++;
            }
          for (size_t i = 0; i < nlons; ++i)
            {
              grid.xbounds[k] = xbounds[i];
              grid.ybounds[k] = ybounds[j];
              k++;
            }
          for (size_t i = 0; i < nextra; ++i)
            {
              grid.xbounds[k] = xbounds[nlons];
              grid.ybounds[k] = ybounds[j];
              k++;
            }
        }
    }
}

static void
gen_grid_reg2d(GridDesciption &grid, const char *pline)
{
  if (isdigit((int) *pline))
    {
      grid.type = GRID_LONLAT;
      grid.xsize = atol(pline);
      while (isdigit((int) *pline)) pline++;
      if (*pline == 'x' || *pline == '/' || *pline == '_')
        pline++;
      else
        {
          grid.type = CDI_UNDEFID;
          return;
        }
      grid.ysize = atol(pline);
      while (isdigit((int) *pline)) pline++;

      grid.xfirst = 0.0;
      grid.yfirst = 0.0;
    }
}

static void
gen_grid_point(GridDesciption &grid, const char *pline)
{
  if (isdigit((int) *pline) || ispunct((int) *pline) || *pline == '-')
    {
      grid.type = GRID_LONLAT;
      grid.xsize = 1;
      grid.ysize = 1;
      grid.xvals.resize(1);
      grid.yvals.resize(1);
      grid.xvals[0] = atof(pline);
      while (isdigit((int) *pline) || ispunct((int) *pline) || *pline == '-') pline++;
      if (*pline == '_') pline++;
      if (strncmp(pline, "lat=", 4) != 0)
        {
          grid.type = CDI_UNDEFID;
          return;
        }
      pline += 4;
      if (isdigit((int) *pline) || ispunct((int) *pline) || *pline == '-')
        grid.yvals[0] = atof(pline);
      else
        grid.type = CDI_UNDEFID;
    }
}

int
grid_from_name(const char *gridnameptr)
{
  const char *pline;
  int gridID = CDI_UNDEFID;
  GridDesciption grid;
  size_t len;
  char *endptr;

  char *gridname = strdup(gridnameptr);
  cstr_to_lower(gridname);

  if (gridname[0] == 't')  // t<RES>grid or t<RES>spec
    {
      int off = 0;
      int type = 'q';
      if (gridname[1] == 'l')
        {
          type = 'l';
          off = 1;
        }
      else if (gridname[1] == 'c')
        {
          type = 'c';
          off = 1;
        }

      pline = &gridname[off + 1];
      if (isdigit((int) *pline))
        {
          grid.ntr = atol(pline);
          while (isdigit((int) *pline)) pline++;
          // clang-format off
          if      (cdo_cmpstrLenRhs(pline, "grid", len)) grid.type = GRID_GAUSSIAN;
          else if (cdo_cmpstrLenRhs(pline, "zon", len))  grid.type = GRID_GAUSSIAN;
          else if (cdo_cmpstrLenRhs(pline, "spec", len)) grid.type = GRID_SPECTRAL;
          else if (cdo_cmpstrLenRhs(pline, "", len))     grid.type = GRID_SPECTRAL;
          // clang-format on

          if (pline[len] != 0) return gridID;

          if (grid.type == GRID_GAUSSIAN)
            {
              if (type == 'l')
                grid.ysize = ntr_to_nlat_linear(grid.ntr);
              else if (type == 'c')
                grid.ysize = ntr_to_nlat_cubic(grid.ntr);
              else
                grid.ysize = ntr_to_nlat(grid.ntr);

              grid.numLPE = grid.ysize / 2;
              grid.xsize = (cdo_cmpstrLenRhs(pline, "zon"))
                               ? 1
                               : ((type == 'c') ? nlat_to_nlon_cubic(grid.ysize) : nlat_to_nlon(grid.ysize));

              grid.xfirst = 0.0;
              grid.yfirst = 0.0;
              grid.yvals.resize(grid.ysize);
              grid.ybounds.resize(grid.ysize * 2);

              auto nlat = grid.ysize;
              std::vector<double> lats(nlat), lat_bounds(nlat + 1);
              gaussian_latitudes_in_degrees(lats, lat_bounds, nlat);

              for (size_t j = 0; j < nlat; ++j) grid.yvals[j] = lats[j];
              for (size_t j = 0; j < nlat; ++j) grid.ybounds[j * 2 + 1] = lat_bounds[j];
              for (size_t j = 0; j < nlat; ++j) grid.ybounds[j * 2] = lat_bounds[j + 1];
            }
        }
    }
  else if (gridname[0] == 'r')  // r<LON>x<LAT>; regular 2D grid
    {
      gen_grid_reg2d(grid, &gridname[1]);
    }
  else if (cdo_cmpstrLenRhs(gridname, "lon=", len))  // lon=<LON>_lat=<LAT>; one gridpoint
    {
      gen_grid_point(grid, &gridname[len]);
    }
  else if (gridname[0] == 'g' && gridname[1] == 'm' && gridname[2] == 'e')  // gme<NI>
    {
      gen_grid_gme(grid, &gridname[3]);
    }
  else if (gridname[0] == 'n' && gridname[1] == 'i')  // ni<NI>
    {
      gen_grid_gme(grid, &gridname[2]);
    }
  else if (gridname[0] == 'h' && gridname[1] == 'p')  // healpix  hp<nside>[_order] (order=[nested|ring|xy])
    {
      gen_grid_healpix(grid, &gridname[2]);
    }
  else if (gridname[0] == 'g' && gridname[1] == 'e' && gridname[2] == 'a')
    {
      gen_grid_gea(grid, &gridname[3]);  // gea<DX>: gaussian reduced equal area; DX in km
    }
  else if ((gridname[0] == 'f' || gridname[0] == 'n') && isdigit((int) gridname[1]))
    {
      // FXXX - full (regular) Gaussian grid with XXX latitude lines between the pole and equator
      pline = &gridname[1];
      auto numLPE = strtol(pline, &endptr, 10);
      pline = endptr;

      if (*pline == 'b')
        {
          grid.genBounds = true;
          pline++;
        }

      if (*pline == '_') pline++;

      if (cdo_cmpstrLenRhs(pline, "zon", len))
        {
          grid.xsize = 1;
          pline += len;
        }

      if (*pline == 0)
        {
          grid.type = GRID_GAUSSIAN;
          grid.numLPE = numLPE;
          grid.ysize = numLPE * 2;
          if (!grid.xsize) grid.xsize = nlat_to_nlon(grid.ysize);

          grid.xfirst = 0.0;
          grid.yfirst = 0.0;
          /* this will change the result of remapcon
          grid.yvals.resize(grid.ysize);
          grid.ybounds.resize(grid.ysize * 2);

          size_t nlat = grid.ysize;
          std::vector<double> lats(nlat), lat_bounds(nlat + 1);
          gaussian_latitudes_in_degrees(lats, lat_bounds, nlat);

          for (size_t j = 0; j < nlat; ++j) grid.yvals[j] = lats[j];
          for (size_t j = 0; j < nlat; ++j) grid.ybounds[j * 2 + 1] = lat_bounds[j];
          for (size_t j = 0; j < nlat; ++j) grid.ybounds[j * 2] = lat_bounds[j + 1];
          */
        }
    }
  else if (gridname[0] == 'o' && isdigit((int) gridname[1]))  // O<xxx>
    {
      pline = &gridname[1];
      auto numLPE = strtol(pline, &endptr, 10);
      pline = endptr;

      if (*pline == 'b')
        {
          grid.genBounds = true;
          pline++;
        }

      if (*pline == '_') pline++;

      if (cdo_cmpstrLenRhs(pline, "zon", len))
        {
          grid.xsize = 1;
          pline += len;
        }

      if (*pline == 0)
        {
          grid.type = GRID_GAUSSIAN;
          grid.numLPE = numLPE;
          grid.ysize = numLPE * 2;
          if (!grid.xsize) grid.xsize = nlat_to_nlon(grid.ysize) + 16;

          grid.xfirst = 0.0;
          grid.yfirst = 0.0;
        }
    }
  else if (gridname[0] == 'g' && isdigit(gridname[1]))  // g<LON>x<LAT> or g<SIZE>
    {
      pline = &gridname[1];
      if (isdigit((int) *pline))
        {
          grid.type = GRID_GENERIC;
          grid.xsize = atol(pline);
          while (isdigit((int) *pline)) pline++;
          if (*pline)
            {
              pline++;
              grid.ysize = atol(pline);
              while (isdigit((int) *pline)) pline++;
            }
          else if (grid.xsize == 1)
            {
              grid.size = 1;
              grid.xsize = 0;
            }
        }
    }
  else if (cdo_cmpstrLenRhs(gridname, "dcw:", len))  // dcw:code_Xdeg
    {
      gen_grid_dcw(grid, &gridname[len], 0.1);
    }
  else if (cdo_cmpstrLenRhs(gridname, "germany", len))  // germany_Xdeg
    {
      gen_grid_lonlat(grid, &gridname[len], 0.1, 5.6, 15.2, 47.1, 55.1);
    }
  else if (cdo_cmpstrLenRhs(gridname, "europe", len))  // europe_Xdeg
    {
      gen_grid_lonlat(grid, &gridname[len], 1, -30, 60, 30, 80);
    }
  else if (cdo_cmpstrLenRhs(gridname, "africa", len))  // africa_Xdeg
    {
      gen_grid_lonlat(grid, &gridname[len], 1, -20, 60, -40, 40);
    }
  else if (cdo_cmpstrLenRhs(gridname, "global", len))  // global_Xdeg
    {
      gen_grid_lonlat(grid, &gridname[len], 1, -180, 180, -90, 90);
    }
  else if (cdo_cmpstrLenRhs(gridname, "zonal", len))  // zonal_Xdeg
    {
      gen_grid_zonal(grid, &gridname[len], 1, -180, 180, -90, 90);
    }
  else if (gridname[0] == 'z')  // z<LAT>; zonal unstructured grid with <LAT> latitudes
    {
      gen_grid_zonal(grid, &gridname[1]);
    }
  else if (cdo_cmpstrLenRhs(gridname, "ico", len))  // icoR02BXX
    {
      gen_grid_icosphere(grid, &gridname[len]);
    }

  if (grid.type != -1) gridID = grid_define(grid);

  free(gridname);

  return gridID;
}
