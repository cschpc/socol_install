/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/
#ifndef STDNAMETABLE_H
#define STDNAMETABLE_H

#include <string>

enum stdnameid
{
  air_pressure,
  pressure_thickness,
  surface_geopotential,
  geopotential,
  air_temperature,
  specific_humidity,
  surface_air_pressure,
  air_pressure_at_sea_level,
  geopotential_height,
  geometric_height_at_full_level_center,
  geometric_height_at_half_level_center
};

int var_echamcode(int varid);
const char *var_name(int varid);
const char *var_stdname(int varid);
const char *var_units(int varid);

int stdname_to_echamcode(const std::string &stdname);

struct gribcode_t
{
  int geopot = 0;
  int temp = 0;
  int hum = 0;
  int ps = 0;
  int lsp = 0;
  int gheight = 0;
  int wind = 0;
  int uwind = 0;
  int vwind = 0;
};

void echam_gribcodes(gribcode_t *gribcodes);
void wmo_gribcodes(gribcode_t *gribcodes);
void hirlam_harmonie_gribcodes(gribcode_t *gribcodes);

#endif
