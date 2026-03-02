/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <string.h>
#include <assert.h>
#include "stdnametable.h"

struct stdnametable_t
{
  int varid;
  int echamcode;
  const char *name;
  const char *stdname;  // Standard name
  const char *units;    // Units
};

static const stdnametable_t stdnametable[] = {
  // clang-format off
  // varid                                  code    name           standard name                             units
  { air_pressure,                              1,  "apres",       "air_pressure",                           "Pa" },
  { pressure_thickness,                        2,  "dpress",      "pressure_thickness",                     "Pa" },
  { surface_geopotential,                    129,  "geosp",       "surface_geopotential",                   "m2 s-2" },
  { geopotential,                            129,  "z",           "geopotential",                           "m2 s-2" },
  { air_temperature,                         130,  "ta",          "air_temperature",                        "K" },
  { specific_humidity,                       133,  "hus",         "specific_humidity",                      "1" },
  { surface_air_pressure,                    134,  "aps",         "surface_air_pressure",                   "Pa" },
  { air_pressure_at_sea_level,               151,  "psl",         "air_pressure_at_sea_level",              "Pa" },
  { geopotential_height,                     156,  "zh",          "geopotential_height",                    "m" },
  { geometric_height_at_full_level_center,    21,  "zg",          "geometric_height_at_full_level_center",  "m" },
  { geometric_height_at_half_level_center,    22,  "zghalf",      "geometric_height_at_half_level_center",  "m" },
  // clang-format on
};

static int
stdnametable_idx(int varid)
{
  int num_entries = (int) (sizeof(stdnametable) / sizeof(stdnametable_t));

  int idx;
  for (idx = 0; idx < num_entries; ++idx)
    if (stdnametable[idx].varid == varid) break;

  assert(idx < num_entries);

  return idx;
}

int
var_echamcode(int varid)
{
  return stdnametable[stdnametable_idx(varid)].echamcode;
}

const char *
var_name(int varid)
{
  return stdnametable[stdnametable_idx(varid)].name;
}

const char *
var_stdname(int varid)
{
  return stdnametable[stdnametable_idx(varid)].stdname;
}

const char *
var_units(int varid)
{
  return stdnametable[stdnametable_idx(varid)].units;
}

int
stdname_to_echamcode(const std::string &stdname)
{
  int code = -1;

  // clang-format off
  if      (stdname == var_stdname(surface_geopotential))      code = 129;
  else if (stdname == var_stdname(geopotential))              code = 129;
  else if (stdname == var_stdname(air_temperature))           code = 130;
  else if (stdname == var_stdname(specific_humidity))         code = 133;
  else if (stdname == var_stdname(surface_air_pressure))      code = 134;
  else if (stdname == var_stdname(air_pressure_at_sea_level)) code = 151;
  else if (stdname == var_stdname(geopotential_height))       code = 156;
  // clang-format on

  return code;
}

void
echam_gribcodes(gribcode_t *gribcodes)
{
  // clang-format off
  gribcodes->geopot  =  129;
  gribcodes->temp    =  130;
  gribcodes->hum     =  133;
  gribcodes->ps      =  134;
  gribcodes->lsp     =  152;
  gribcodes->gheight =  156;
  gribcodes->wind    =    0;  
  gribcodes->uwind   =  131;
  gribcodes->vwind   =  132;
  // clang-format on
}

void
wmo_gribcodes(gribcode_t *gribcodes)
{
  // clang-format off
  gribcodes->geopot  =   6;
  gribcodes->temp    =  11;
  gribcodes->hum     =   0;
  gribcodes->ps      =   1;
  gribcodes->lsp     =   0;
  gribcodes->gheight =   7;
  gribcodes->wind    =  10;  
  gribcodes->uwind   = 131;
  gribcodes->vwind   = 132;
  // clang-format on
  /*  ECMWF (IFS) GLOBAL Model
   *
   *
   http://old.ecmwf.int/publications/manuals/d/gribapi/param/filter=grib1/order=paramId/order_type=asc/p=1/search=wind/table=128/
   Wind speed	ws	m s**-1	10	grib1	grib2
   U component of wind	u	m s**-1	131	grib1	grib2	netcdf
   V component of wind	v	m s**-1	132	grib1	grib2	netcdf
   10 metre U wind component	10u	m s**-1	165	grib1	grib2
   10 metre V wind component	10v	m s**-1	166	grib1	grib2
   10 metre wind speed	10si	m s**-1	207	grib1	grib2
  */
}

void
hirlam_harmonie_gribcodes(gribcode_t *gribcodes)
{
  // clang-format off
  gribcodes->geopot  =   6; // Geopotential [m2/s2]
  gribcodes->temp    =  11; // Temperature [K]
  gribcodes->hum     =  51; // Specific humidity [kg/kg]
  gribcodes->ps      =   1; // Pressure [Pa]
  gribcodes->lsp     =   0; // - not available -
  gribcodes->gheight =   7; // 	007 Geopotential height [Gpm]; 008 Geometric height [m]
  gribcodes->wind    =  32;  
  gribcodes->uwind   =  33;
  gribcodes->vwind   =  34;
  // clang-format on
  /*
   1/103/0 	"mean sea level pressure" pressure_msl   [Pa]
   1/105/0 	"surface pressure"  pressure_surf  [Pa]
   032	Wind speed          	m/s 	WIND
   033	u-component of wind 	m/s 	UGRD
   034	v-component of wind 	m/s 	VGRD
   051	Specific humidity
   052	Relative humidity
   NOT AVAILABLE:
   - Natural log of surface pressure 	ln(kPa) 	NLGSP
  */
}
