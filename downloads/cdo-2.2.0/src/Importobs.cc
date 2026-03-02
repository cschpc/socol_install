/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "varray.h"
#include "readline.h"
#include <mpim_grid.h>
#include "griddes.h"

static void
init_vars(int vlistID, int gridID, int zaxisID, int nvars)
{
  const int code[] = { 11, 17, 33, 34, 1, 2 /*, 3*/ };
  const char *name[] = { "temp", "depoint", "u", "v", "height", "pressure" /*, "station"*/ };
  const char *units[] = { "Celsius", "", "m/s", "m/s", "m", "hPa" /*, ""*/ };

  for (int i = 0; i < nvars; ++i)
    {
      const auto varID = vlistDefVar(vlistID, gridID, zaxisID, TIME_VARYING);
      vlistDefVarParam(vlistID, varID, cdiEncodeParam(code[i], 255, 255));
      cdiDefKeyString(vlistID, varID, CDI_KEY_NAME, name[i]);
      cdiDefKeyString(vlistID, varID, CDI_KEY_UNITS, units[i]);
      vlistDefVarDatatype(vlistID, varID, CDI_DATATYPE_FLT32);
    }
}

static void
init_data(int vlistID, int nvars, Varray2D<double> &data)
{
  for (int varID = 0; varID < nvars; ++varID)
    {
      const auto gridsize = gridInqSize(vlistInqVarGrid(vlistID, varID));
      const auto missval = vlistInqVarMissval(vlistID, varID);
      for (size_t i = 0; i < gridsize; ++i) data[varID][i] = missval;
    }
}

static void
write_data(CdoStreamID streamID, int vlistID, int nvars, Varray2D<double> &data)
{
  for (int varID = 0; varID < nvars; ++varID)
    {
      const auto gridsize = gridInqSize(vlistInqVarGrid(vlistID, varID));
      const auto missval = vlistInqVarMissval(vlistID, varID);
      const auto nmiss = varray_num_mv(gridsize, data[varID], missval);
      cdo_def_record(streamID, varID, 0);
      cdo_write_record(streamID, data[varID].data(), nmiss);
    }
}

static int
get_date(const char *name)
{
  const char *pname = strchr(name, '_');
  const int date = pname ? atoi(pname + 1) : 0;
  return date;
}

#define MAX_LINE_LEN 4096

void *
Importobs(void *process)
{
  char line[MAX_LINE_LEN];
  size_t i, j;
  long index;
  constexpr int nvars = 6;
  int vtime = 0;
  char dummy[32], station[32], datetime[32];
  float lat, lon, height1, pressure, height2, value;
  double latmin = 90, latmax = -90, lonmin = 360, lonmax = -360;
  int code;

  cdo_initialize(process);

  cdo_operator_add("import_obs", 0, 0, "grid description file or name");

  const auto operatorID = cdo_operator_id();

  operator_input_arg(cdo_operator_enter(operatorID));

  const auto gridID = cdo_define_grid(cdo_operator_argv(0));

  if (gridInqType(gridID) != GRID_LONLAT) cdo_abort("Unsupported grid type: %s", gridNamePtr(gridInqType(gridID)));

  const auto gridsize = gridInqSize(gridID);
  const auto xsize = gridInqXsize(gridID);
  const auto ysize = gridInqYsize(gridID);

  Varray<double> xvals(gridsize), yvals(gridsize);
  gridInqXvals(gridID, xvals.data());
  gridInqYvals(gridID, yvals.data());

  // Convert lat/lon units if required
  cdo_grid_to_degree(gridID, CDI_XAXIS, gridsize, xvals.data(), "grid center lon");
  cdo_grid_to_degree(gridID, CDI_YAXIS, gridsize, yvals.data(), "grid center lat");

  auto fp = std::fopen(cdo_get_stream_name(0), "r");
  if (fp == nullptr)
    {
      perror(cdo_get_stream_name(0));
      exit(EXIT_FAILURE);
    }

  auto vdate = get_date(cdo_get_stream_name(0));
  if (vdate <= 999999) vdate = vdate * 100 + 1;

  const auto streamID = cdo_open_write(1);

  const auto zaxisID = zaxisCreate(ZAXIS_SURFACE, 1);

  const auto vlistID = vlistCreate();

  const auto taxisID = cdo_taxis_create(TAXIS_ABSOLUTE);
  vlistDefTaxis(vlistID, taxisID);

  Varray2D<double> data(nvars);
  for (i = 0; i < nvars; ++i) data[i].resize(gridsize);

  init_vars(vlistID, gridID, zaxisID, nvars);

  cdo_def_vlist(streamID, vlistID);

  int vdate0 = 0;
  int vtime0 = 0;
  // ntime = 0;
  int tsID = 0;
  while (cdo::readline(fp, line, MAX_LINE_LEN))
    {
      std::sscanf(line, "%s %s %s %g %g %g %d %g %g %g", dummy, station, datetime, &lat, &lon, &height1, &code, &pressure, &height2,
             &value);
      long vdate_l;
      std::sscanf(datetime, "%ld_%d", &vdate_l, &vtime);
      vdate = vdate_l;

      if (vdate != vdate0 || vtime != vtime0)
        {
          if (tsID > 0) write_data(streamID, vlistID, nvars, data);

          vdate0 = vdate;
          vtime0 = vtime;
          // printf("%s %d %d %g %g %g %d %g %g %g\n", station, vdate, vtime, lat, lon, height1, code, pressure, height2, value);
          CdiDateTime vDateTime{};
          vDateTime.date = cdiDate_set(vdate);
          vDateTime.time = cdiTime_set(vtime);
          taxisDefVdatetime(taxisID, vDateTime);
          cdo_def_timestep(streamID, tsID);

          init_data(vlistID, nvars, data);

          tsID++;
        }

      if (lon < lonmin) lonmin = lon;
      if (lon > lonmax) lonmax = lon;
      if (lat < latmin) latmin = lat;
      if (lat > latmax) latmax = lat;

      const auto dy = yvals[1] - yvals[0];
      for (j = 0; j < ysize; ++j)
        if (lat >= (yvals[j] - dy / 2) && lat < (yvals[j] + dy / 2)) break;

      const auto dx = xvals[1] - xvals[0];
      if (lon < (xvals[0] - dx / 2) && lon < 0) lon += 360;
      for (i = 0; i < xsize; ++i)
        if (lon >= (xvals[i] - dx / 2) && lon < (xvals[i] + dx / 2)) break;

      index = -1;
      if (code == 11) index = 0;
      if (code == 17) index = 1;
      if (code == 33) index = 2;
      if (code == 34) index = 3;

      // printf("%d %d %d %g %g %g %g\n", i, j, index, dx, dy, lon, lat);
      if (i < xsize && j < ysize && index >= 0)
        {
          char *pstation = station;
          while (isalpha(*pstation)) pstation++;
          // printf("station %s %d\n", pstation, atoi(pstation));
          data[index][j * xsize + i] = value;
          data[4][j * xsize + i] = height1;
          data[5][j * xsize + i] = pressure;
          // data[    6][j*xsize+i] = atoi(pstation);
        }

      /*
        printf("%s %d %d %g %g %g %d %g %g %g\n", station, vdate, vtime, lat, lon, height1, code, pressure, height2, value);
      */
    }

  write_data(streamID, vlistID, nvars, data);

  if (Options::cdoVerbose) printf("lonmin=%g, lonmax=%g, latmin=%g, latmax=%g\n", lonmin, lonmax, latmin, latmax);

  process_def_var_num(vlistNvars(vlistID));

  cdo_stream_close(streamID);

  std::fclose(fp);

  vlistDestroy(vlistID);
  gridDestroy(gridID);
  zaxisDestroy(zaxisID);
  taxisDestroy(taxisID);

  cdo_finish();

  return nullptr;
}
