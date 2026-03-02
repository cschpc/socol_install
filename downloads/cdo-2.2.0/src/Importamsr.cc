/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include "process_int.h"

#define NLON 1440
#define NLAT 720
#define MAX_VARS 6

static void
init_amsr_day(int vlistID, int gridID, int zaxisID, int nvars)
{
  /*
    Version-5 RSS AMSR-E or AMSR-J daily files

    filename  with path in form satname_yyyymmdd_v5.gz
    where satname  = name of satellite (amsre or amsr)
             yyyy	= year
               mm	= month
               dd	= day of month

    1:time	time of measurement in fractional hours GMT
    2:sst     	sea surface temperature in deg Celcius
    3:wind	10m surface wind in meters/second
    4:vapor	columnar water vapor in millimeters
    5:cloud	cloud liquid water in millimeters
    6:rain   	rain rate in millimeters/hour
  */
  const char *name[] = { "hours", "sst", "wind", "vapor", "cloud", "rain" };
  const char *units[] = { "h", "deg Celcius", "m/s", "mm", "mm", "mm/h" };
  const double xscale[] = { 0.1, 0.15, 0.2, 0.3, 0.01, 0.1 };
  const double xminval[] = { 0., -3., 0., 0., 0., 0. };

  for (int i = 0; i < nvars; ++i)
    {
      const int varID = vlistDefVar(vlistID, gridID, zaxisID, TIME_VARYING);
      cdiDefKeyString(vlistID, varID, CDI_KEY_NAME, name[i]);
      cdiDefKeyString(vlistID, varID, CDI_KEY_UNITS, units[i]);
      vlistDefVarDatatype(vlistID, varID, CDI_DATATYPE_INT16);
      vlistDefVarMissval(vlistID, varID, 254);
      cdiDefKeyFloat(vlistID, varID, CDI_KEY_SCALEFACTOR, xscale[i]);
      cdiDefKeyFloat(vlistID, varID, CDI_KEY_ADDOFFSET, xminval[i]);
    }
}

static void
init_amsr_averaged(int vlistID, int gridID, int zaxisID, int nvars)
{
  /*
    Version-5 AMSR-E or AMSR-J time-averaged files including:
          3-day		(average of 3 days ending on file date)
          weekly	(average of 7 days ending on Saturday of file date)
          monthly	(average of all days in month)


     filename
           format of file names are:
                3-day	satname_yyyymmddv5_d3d.gz
                weekly	satname_yyyymmddv5.gz
                monthly	satname_yyyymmv5.gz

        where	satname	=name of satellite (amsre or amsr)
                        yyyy	=year
                        mm	=month
                        dd	=day of month

    1:sst       sea surface temperature in deg Celcius
    2:wind      10m surface wind in meters/second
    3:vapor	columnar water vapor in millimeters
    4:cloud     cloud liquid water in millimeters
    5:rain	rain rate in millimeters/hour
  */
  const char *name[] = { "sst", "wind", "vapor", "cloud", "rain" };
  const char *units[] = { "deg Celcius", "m/s", "mm", "mm", "mm/h" };
  const double xscale[] = { 0.15, 0.2, 0.3, 0.01, 0.1 };
  const double xminval[] = { -3., 0., 0., 0., 0. };

  for (int i = 0; i < nvars; ++i)
    {
      // varID = vlistDefVar(vlistID, gridID, zaxisID, TIME_CONSTANT);
      const int varID = vlistDefVar(vlistID, gridID, zaxisID, TIME_VARYING);
      cdiDefKeyString(vlistID, varID, CDI_KEY_NAME, name[i]);
      cdiDefKeyString(vlistID, varID, CDI_KEY_UNITS, units[i]);
      vlistDefVarDatatype(vlistID, varID, CDI_DATATYPE_INT16);
      vlistDefVarMissval(vlistID, varID, 254);
      cdiDefKeyFloat(vlistID, varID, CDI_KEY_SCALEFACTOR, xscale[i]);
      cdiDefKeyFloat(vlistID, varID, CDI_KEY_ADDOFFSET, xminval[i]);
    }
}

static void
read_amsr(FILE *fp, int vlistID, int nvars, Varray2D<double> &data, size_t *nmiss)
{
  std::vector<unsigned char> amsr_data;

  for (int varID = 0; varID < nvars; ++varID)
    {
      const size_t gridsize = gridInqSize(vlistInqVarGrid(vlistID, varID));
      amsr_data.resize(gridsize);
      const size_t size = fread(amsr_data.data(), 1, gridsize, fp);
      if (size != gridsize) cdo_abort("Read error!");

      const double missval = vlistInqVarMissval(vlistID, varID);
      double xminval = 0.0, xscale = 1.0;
      cdiInqKeyFloat(vlistID, varID, CDI_KEY_SCALEFACTOR, &xscale);
      cdiInqKeyFloat(vlistID, varID, CDI_KEY_ADDOFFSET, &xminval);

      nmiss[varID] = 0;
      for (size_t i = 0; i < gridsize; ++i)
        {
          if (amsr_data[i] <= 250) { data[varID][i] = amsr_data[i] * xscale + xminval; }
          else
            {
              data[varID][i] = missval;
              nmiss[varID]++;
            }
        }
    }
}

static void
write_data(CdoStreamID streamID, int nvars, Varray2D<double> &data, size_t *nmiss)
{
  for (int varID = 0; varID < nvars; ++varID)
    {
      cdo_def_record(streamID, varID, 0);
      cdo_write_record(streamID, data[varID].data(), nmiss[varID]);
    }
}

static int
get_date(const char *name)
{
  const char *pname = strchr(name, '_');
  const int date = pname ? atoi(pname + 1) : 0;
  return date;
}

void *
Importamsr(void *process)
{
  int tsID;
  int nvars;
  int vtime = 0;
  double xvals[NLON], yvals[NLAT];
  Varray2D<double> data(MAX_VARS);
  size_t nmiss[MAX_VARS];

  cdo_initialize(process);

  operator_check_argc(0);

  auto fp = std::fopen(cdo_get_stream_name(0), "r");
  if (fp == nullptr)
    {
      perror(cdo_get_stream_name(0));
      exit(EXIT_FAILURE);
    }

  fseek(fp, 0L, SEEK_END);
  const size_t fsize = (size_t) ftell(fp);
  fseek(fp, 0L, SEEK_SET);

  auto vdate = get_date(cdo_get_stream_name(0));
  if (vdate <= 999999) vdate = vdate * 100 + 1;

  const auto streamID = cdo_open_write(1);

  /*
    Longitude  is 0.25*xdim-0.125    degrees east
    Latitude   is 0.25*ydim-90.125
  */
  const size_t gridsize = NLON * NLAT;
  const auto gridID = gridCreate(GRID_LONLAT, gridsize);
  gridDefXsize(gridID, NLON);
  gridDefYsize(gridID, NLAT);
  for (int i = 0; i < NLON; ++i) xvals[i] = 0.25 * (i + 1) - 0.125;
  for (int i = 0; i < NLAT; ++i) yvals[i] = 0.25 * (i + 1) - 90.125;
  gridDefXvals(gridID, xvals);
  gridDefYvals(gridID, yvals);

  const auto zaxisID = zaxisCreate(ZAXIS_SURFACE, 1);

  const auto vlistID = vlistCreate();
  const auto taxisID = cdo_taxis_create(TAXIS_ABSOLUTE);
  vlistDefTaxis(vlistID, taxisID);

  if (fsize == 12441600)
    {
      nvars = 6;
      for (int i = 0; i < nvars; ++i) data[i].resize(gridsize);

      init_amsr_day(vlistID, gridID, zaxisID, nvars);

      cdo_def_vlist(streamID, vlistID);

      vtime = 13000;  // 1:30:00
      for (tsID = 0; tsID < 2; ++tsID)
        {
          CdiDateTime vDateTime{};
          vDateTime.date = cdiDate_set(vdate);
          vDateTime.time = cdiTime_set(vtime);
          taxisDefVdatetime(taxisID, vDateTime);
          vtime += 120000;  // 13:30:00
          cdo_def_timestep(streamID, tsID);

          read_amsr(fp, vlistID, nvars, data, nmiss);

          write_data(streamID, nvars, data, nmiss);
        }
    }
  else if (fsize == 5184000)
    {
      nvars = 5;
      for (int i = 0; i < nvars; ++i) data[i].resize(gridsize);

      init_amsr_averaged(vlistID, gridID, zaxisID, nvars);

      // vlistDefNtsteps(vlistID, 0);
      cdo_def_vlist(streamID, vlistID);

      CdiDateTime vDateTime{};
      vDateTime.date = cdiDate_set(vdate);
      vDateTime.time = cdiTime_set(vtime);
      taxisDefVdatetime(taxisID, vDateTime);
      tsID = 0;
      cdo_def_timestep(streamID, tsID);

      read_amsr(fp, vlistID, nvars, data, nmiss);

      write_data(streamID, nvars, data, nmiss);
    }
  else
    cdo_abort("Unexpected file size for AMSR data!");

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
