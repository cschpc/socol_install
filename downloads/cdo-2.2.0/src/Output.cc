/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Output     output          ASCII output
      Output     outputf         Formatted output
      Output     outputint       Integer output
      Output     outputsrv       SERVICE output
      Output     outputext       EXTRA output
      Output     outputtab       Table output
*/

#include <cdi.h>

#include <unordered_map>

#include "cdo_options.h"
#include "process_int.h"
#include "param_conversion.h"
#include <mpim_grid.h>
#include "printinfo.h"
#include "cdo_zaxis.h"

static void
outputarr(int dig, size_t gridsize, const Varray<double> &array)
{
  for (size_t i = 0; i < gridsize; ++i) { fprintf(stdout, "  arr[%zu] = %.*g;\n", i, dig, array[i]); }
}

static void
outputsp(size_t gridsize, const Varray<double> &array, long ntr)
{
  auto mm = varray_min_max(gridsize, array);
  if (/* T11 */ mm.min >= -1 && mm.max <= 12)
    {
      auto spc = array.data();
      for (long m = 0; m <= ntr; ++m)
        {
          for (long n = m; n <= ntr; ++n)
            {
              fprintf(stdout, "%3d", (int) *spc++);
              fprintf(stdout, "%3d", (int) *spc++);
            }
          fprintf(stdout, "\n");
        }
    }
}

static void
output(size_t gridsize, const Varray<double> &array)
{
  int nout = 0;
  for (size_t i = 0; i < gridsize; ++i)
    {
      if (nout == 6)
        {
          nout = 0;
          fprintf(stdout, "\n");
        }
      fprintf(stdout, " %12.6g", array[i]);
      nout++;
    }
  fprintf(stdout, "\n");
}

static void
outputxyz(size_t gridsize, const Varray<double> &array, double missval, size_t nlon, size_t nlat, const Varray<double> &lon,
          const Varray<double> &lat)
{
  double fmin = 0.0;
  double x, y, z;
  for (size_t i = 0; i < gridsize; ++i)
    if (!dbl_is_equal(array[i], missval))
      {
        if (array[i] < fmin) fmin = array[i];
        fprintf(stdout, "%g\t%g\t%g\t%g\n", lon[i], lat[i], array[i], array[i]);
      }
  const char *fname = "frontplane.xyz";
  auto fp = std::fopen(fname, "w");
  if (fp == nullptr) cdo_abort("Open failed on %s", fname);
  // first front plane
  auto dx = (lon[1] - lon[0]);
  auto x0 = lon[0] - dx / 2;
  auto y0 = lat[0] - dx / 2;
  auto z0 = fmin;
  fprintf(fp, ">\n");
  for (size_t i = 0; i < nlon; ++i)
    {
      x = x0;
      y = y0;
      z = z0;
      fprintf(fp, "%g %g %g\n", x, y, z);
      x = x0;
      y = y0;
      z = array[i];
      fprintf(fp, "%g %g %g\n", x, y, z);
      x = x0 + dx;
      y = y0;
      fprintf(fp, "%g %g %g\n", x, y, z);
      x0 = x; /*y0 = y0;*/
      z0 = z;
    }
  x = x0;
  y = y0;
  z = fmin;
  fprintf(fp, "%g %g %g\n", x, y, z);
  x = lon[0] - dx / 2;
  fprintf(fp, "%g %g %g\n", x, y, z);

  // second front plane
  x0 = lon[0] - dx / 2;
  y0 = lat[0] - dx / 2;
  z0 = fmin;
  fprintf(fp, ">\n");
  for (size_t i = 0; i < nlat; ++i)
    {
      x = x0;
      y = y0;
      z = z0;
      fprintf(fp, "%g %g %g\n", x, y, z);
      x = x0;
      y = y0;
      z = array[i * nlon];
      fprintf(fp, "%g %g %g\n", x, y, z);
      x = x0;
      y = y0 + dx;
      fprintf(fp, "%g %g %g\n", x, y, z);
      /*x0 = x0;*/ y0 = y;
      z0 = z;
    }
  x = x0;
  y = y0;
  z = fmin;
  fprintf(fp, "%g %g %g\n", x, y, z);
  y = lat[0] - dx / 2;
  fprintf(fp, "%g %g %g\n", x, y, z);

  std::fclose(fp);
}

static void
read_xy_coordinates(bool hasRegxyCoordinates, int gridID0, Varray<double> &grid_xvals, Varray<double> &grid_yvals)
{
  auto gridsize = gridInqSize(gridID0);
  auto xsize = gridInqXsize(gridID0);
  auto ysize = gridInqYsize(gridID0);

  if (hasRegxyCoordinates)
    {
      grid_xvals.resize(xsize);
      grid_yvals.resize(ysize);
      gridInqXvals(gridID0, grid_xvals.data());
      gridInqYvals(gridID0, grid_yvals.data());
    }
  else
    {
      auto gridIDx = generate_full_point_grid(gridID0);
      if (!gridHasCoordinates(gridIDx)) cdo_abort("Cell center coordinates missing!");

      grid_xvals.resize(gridsize);
      grid_yvals.resize(gridsize);
      gridInqXvals(gridIDx, grid_xvals.data());
      gridInqYvals(gridIDx, grid_yvals.data());

      if (gridIDx != gridID0) gridDestroy(gridIDx);
    }
}

static void
read_lonlat_coordinates(int gridID0, Varray<double> &grid_center_lon, Varray<double> &grid_center_lat)
{
  auto gridsize = gridInqSize(gridID0);
  auto gridIDx = generate_full_point_grid(gridID0);
  if (!gridHasCoordinates(gridIDx)) cdo_abort("Cell center coordinates missing!");

  grid_center_lon.resize(gridsize);
  grid_center_lat.resize(gridsize);
  gridInqXvals(gridIDx, grid_center_lon.data());
  gridInqYvals(gridIDx, grid_center_lat.data());

  // Convert lat/lon units if required
  cdo_grid_to_degree(gridIDx, CDI_XAXIS, gridsize, grid_center_lon.data(), "grid center lon");
  cdo_grid_to_degree(gridIDx, CDI_YAXIS, gridsize, grid_center_lat.data(), "grid center lat");

  if (gridIDx != gridID0) gridDestroy(gridIDx);
}

static void
outputint(size_t gridsize, const Varray<double> &array)
{
  int nout = 0;
  for (size_t i = 0; i < gridsize; ++i)
    {
      if (nout == 8)
        {
          nout = 0;
          fprintf(stdout, "\n");
        }
      fprintf(stdout, " %8d", (int) array[i]);
      nout++;
    }
  fprintf(stdout, "\n");
}

static void
outputf(int nelem, const std::string &format, size_t gridsize, const Varray<double> &array)
{
  int nout = 0;
  for (size_t i = 0; i < gridsize; ++i)
    {
      if (nout == nelem)
        {
          nout = 0;
          fprintf(stdout, "\n");
        }
      fprintf(stdout, format.c_str(), array[i]);
      nout++;
    }
  fprintf(stdout, "\n");
}

void *
Output(void *process)
{
  size_t nmiss;
  int nelem = 1;
  int index;
  std::string format;
  char paramstr[32];
  int year, month, day;
  std::vector<int> keyIndices;

  // clang-format off
  enum {knohead,   kvalue,  kparam,  kcode,  kname,  kx,  ky,  klon,  klat,  klev,  kbin,  kxind,  kyind,  ktimestep,  kdate,  ktime,  kyear,  kmonth,  kday};
  struct KeyLenEntry {
    std::string key;
    int         idx;
    int         len;
  };
  std::vector<KeyLenEntry> keyMap = {
    { "nohead",   knohead,    0 },
    { "value",    kvalue,     8 },
    { "param",    kparam,    11 },
    { "code",     kcode,      4 },
    { "name",     kname,      8 },
    { "x",        kx,         6 },
    { "y",        ky,         6 },
    { "lon",      klon,       6 },
    { "lat",      klat,       6 },
    { "lev",      klev,       6 },
    { "bin",      kbin,       6 },
    { "xind",     kxind,      4 },
    { "yind",     kyind,      4 },
    { "timestep", ktimestep,  6 },
    { "date",     kdate,     10 },
    { "time",     ktime,      8 },
    { "year",     kyear,      5 },
    { "month",    kmonth,     2 },
    { "day",      kday,       2 },
  };

  cdo_initialize(process);

  auto OUTPUT    = cdo_operator_add("output",    0, 1, nullptr);
  auto OUTPUTINT = cdo_operator_add("outputint", 0, 0, nullptr);
  auto OUTPUTSRV = cdo_operator_add("outputsrv", 0, 0, nullptr);
  auto OUTPUTEXT = cdo_operator_add("outputext", 0, 0, nullptr);
  auto OUTPUTF   = cdo_operator_add("outputf",   0, 0, nullptr);
  auto OUTPUTTS  = cdo_operator_add("outputts",  0, 0, nullptr);
  auto OUTPUTFLD = cdo_operator_add("outputfld", 0, 0, nullptr);
  auto OUTPUTARR = cdo_operator_add("outputarr", 0, 0, nullptr);
  auto OUTPUTXYZ = cdo_operator_add("outputxyz", 0, 0, nullptr);
  auto OUTPUTTAB = cdo_operator_add("outputtab", 0, 0, nullptr);
  // clang-format on

  (void) (OUTPUT);  // unused

  auto operatorID = cdo_operator_id();
  auto opercplx = (cdo_operator_f2(operatorID) == 1);

  if (operatorID == OUTPUTF)
    {
      operator_input_arg("format and number of elements [optional]");

      if (cdo_operator_argc() < 1) cdo_abort("Too few arguments!");

      format = cdo_operator_argv(0);
      if (cdo_operator_argc() == 2) nelem = parameter_to_int(cdo_operator_argv(1));
    }
  else if (operatorID == OUTPUTTAB)
    {
      auto lhead = true;

      operator_input_arg("keys to print");

      auto npar = cdo_operator_argc();
      auto parnames = cdo_get_oper_argv();

      if (Options::cdoVerbose)
        for (int i = 0; i < npar; ++i) cdo_print("key%d=%s", i + 1, parnames[i]);

      keyIndices.reserve(npar);
      for (int i = 0; i < npar; ++i)
        {
          const char *currentName = parnames[i].c_str();
          size_t k;
          for (k = 0; k < keyMap.size(); ++k)
            {
              // int len = strlen(parnames[i]);
              int len = keyMap[k].key.size();
              if (len < 3) len = 3;
              if (strncmp(currentName, keyMap[k].key.c_str(), len) == 0)
                {
                  int len2 = strlen(currentName);
                  if (len2 > len && currentName[len] != ':')
                    cdo_abort("Key parameter >%s< contains invalid character at position %d!", currentName, len + 1);

                  if (keyMap[k].idx == knohead)
                    lhead = false;
                  else
                    {
                      keyIndices.push_back(k);
                      if (len2 > len && currentName[len] == ':' && isdigit(currentName[len + 1]))
                        keyMap[k].len = atoi(&currentName[len + 1]);
                    }
                  break;
                }
            }

          if (k == keyMap.size()) cdo_abort("Key >%s< unsupported!", currentName);
        }

      if (Options::cdoVerbose)
        for (auto ki : keyIndices) cdo_print("idx=%d/%d  len=%d  name=%s", ki, keyMap[ki].idx, keyMap[ki].len, keyMap[ki].key);

      if (lhead)
        {
          fprintf(stdout, "#");
          for (auto ki : keyIndices) fprintf(stdout, "%*s ", keyMap[ki].len, keyMap[ki].key.c_str());
          fprintf(stdout, "\n");
        }
    }
  else { operator_check_argc(0); }

  for (int indf = 0; indf < cdo_stream_cnt(); indf++)
    {
      auto streamID = cdo_open_read(indf);
      auto vlistID = cdo_stream_inq_vlist(streamID);

      VarList varList;
      varListInit(varList, vlistID);

      auto ngrids = vlistNgrids(vlistID);
      int ndiffgrids = 0;
      for (index = 1; index < ngrids; ++index)
        if (vlistGrid(vlistID, 0) != vlistGrid(vlistID, index)) ndiffgrids++;

      if (ndiffgrids > 0) cdo_abort("Too many different grids!");

      auto gridID0 = vlistGrid(vlistID, 0);
      auto gridtype = gridInqType(gridID0);
      auto hasRegxyCoordinates = (gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN || gridtype == GRID_PROJECTION);
      auto gridsize = gridInqSize(gridID0);
      auto xsize = gridInqXsize(gridID0);
      size_t nwpv = (vlistNumber(vlistID) == CDI_COMP) ? 2 : 1;
      if (nwpv == 2 && !opercplx) cdo_abort("Fields with complex numbers are not supported by this operator!");
      auto gridsizemax = nwpv * gridsize;
      Varray<double> array(gridsizemax);
      Varray<double> grid_center_lon, grid_center_lat;
      Varray<double> grid_xvals, grid_yvals;

      if (operatorID == OUTPUTTAB) read_xy_coordinates(hasRegxyCoordinates, gridID0, grid_xvals, grid_yvals);

      if (operatorID == OUTPUTFLD || operatorID == OUTPUTXYZ || operatorID == OUTPUTTAB)
        read_lonlat_coordinates(gridID0, grid_center_lon, grid_center_lat);

      int tsID = 0;
      auto taxisID = vlistInqTaxis(vlistID);
      while (true)
        {
          auto nrecs = cdo_stream_inq_timestep(streamID, tsID);
          if (nrecs == 0) break;

          auto vDateTime = taxisInqVdatetime(taxisID);
          auto vDateStr = date_to_string(vDateTime.date);
          auto vTimeStr = time_to_string(vDateTime.time);

          cdiDate_decode(vDateTime.date, &year, &month, &day);

          for (int recID = 0; recID < nrecs; ++recID)
            {
              int varID, levelID;
              cdo_inq_record(streamID, &varID, &levelID);
              const auto &var = varList[varID];

              auto code = var.code;
              auto gridID = var.gridID;
              auto dig = (var.datatype == CDI_DATATYPE_FLT64) ? Options::CDO_dbl_digits : Options::CDO_flt_digits;
              gridsize = var.nwpv * var.gridsize;
              auto nlon = gridInqXsize(gridID);
              auto nlat = gridInqYsize(gridID);
              auto level = cdo_zaxis_inq_level(var.zaxisID, levelID);
              auto missval = var.missval;

              cdiParamToString(var.param, paramstr, sizeof(paramstr));

              if (nlon * nlat != gridsize)
                {
                  nlon = gridsize;
                  nlat = 1;
                }

              cdo_read_record(streamID, array.data(), &nmiss);

              auto vdate = cdiDate_get(vDateTime.date);
              auto vtime = cdiTime_get(vDateTime.time);

              if (operatorID == OUTPUTSRV)
                fprintf(stdout, "%4d %8g %8ld %4d %8zu %8zu %d %d\n", code, level, (long) vdate, vtime, nlon, nlat, 0, 0);

              if (operatorID == OUTPUTEXT) fprintf(stdout, "%8ld %4d %8g %8zu\n", (long) vdate, code, level, gridsize);

              if (operatorID == OUTPUTINT) { outputint(gridsize, array); }
              else if (operatorID == OUTPUTF) { outputf(nelem, format, gridsize, array); }
              else if (operatorID == OUTPUTTS)
                {
                  if (gridsize > 1) cdo_abort("operator works only with one gridpoint!");

                  fprintf(stdout, "%s %s %.*g\n", vDateStr.c_str(), vTimeStr.c_str(), dig, array[0]);
                }
              else if (operatorID == OUTPUTFLD)
                {
                  int hour, minute, second, ms;
                  cdiTime_decode(vDateTime.time, &hour, &minute, &second, &ms);
                  double xdate = vdate - (vdate / 100) * 100 + (hour * 3600 + minute * 60 + second) / 86400.0;
                  for (size_t i = 0; i < gridsize; ++i)
                    if (!dbl_is_equal(array[i], missval))
                      fprintf(stdout, "%g\t%g\t%g\t%.*g\n", xdate, grid_center_lat[i], grid_center_lon[i], dig, array[i]);
                }
              else if (operatorID == OUTPUTTAB)
                {
                  auto is2dGrid = (gridtype == GRID_CURVILINEAR || gridtype == GRID_LONLAT || gridtype == GRID_PROJECTION);
                  for (size_t i = 0; i < gridsize; ++i)
                    {
                      auto yind = i;
                      auto xind = i;
                      if (is2dGrid)
                        {
                          yind /= xsize;
                          xind -= yind * xsize;
                        }
                      auto lon = grid_center_lon[i];
                      auto lat = grid_center_lat[i];

                      auto xval = hasRegxyCoordinates ? grid_xvals[xind] : grid_xvals[i];
                      auto yval = hasRegxyCoordinates ? grid_yvals[yind] : grid_yvals[i];

                      for (auto ki : keyIndices)
                        {
                          auto len = keyMap[ki].len;
                          // clang-format off
                          switch (keyMap[ki].idx)
                            {
                            case kvalue:    fprintf(stdout, "%*.*g ", len, dig, array[i]); break;
                            case kx:        fprintf(stdout, "%*.*g ", len, dig, xval); break;
                            case ky:        fprintf(stdout, "%*.*g ", len, dig, yval); break;
                            case klon:      fprintf(stdout, "%*.*g ", len, dig, lon); break;
                            case klat:      fprintf(stdout, "%*.*g ", len, dig, lat); break;
                            case klev:      fprintf(stdout, "%*.*g ", len, dig, level); break;
                            case kbin:      fprintf(stdout, "%*.*g ", len, dig, level); break;
                            case kparam:    fprintf(stdout, "%*s ", len, paramstr); break;
                            case kcode:     fprintf(stdout, "%*d ", len, code); break;
                            case kname:     fprintf(stdout, "%*s ", len, var.name.c_str()); break;
                            case kxind:     fprintf(stdout, "%*zu ", len, xind + 1); break;
                            case kyind:     fprintf(stdout, "%*zu ", len, yind + 1); break;
                            case ktimestep: fprintf(stdout, "%*d ", len, tsID + 1); break;
                            case kdate:     fprintf(stdout, "%*s ", len, vDateStr.c_str()); break;
                            case ktime:     fprintf(stdout, "%*s ", len, vTimeStr.c_str()); break;
                            case kyear:     fprintf(stdout, "%*d ", len, year); break;
                            case kmonth:    fprintf(stdout, "%*d ", len, month); break;
                            case kday:      fprintf(stdout, "%*d ", len, day); break;
                            }
                          // clang-format on
                        }
                      fprintf(stdout, "\n");
                    }
                }
              else if (operatorID == OUTPUTXYZ)
                {
                  if (tsID == 0 && recID == 0) outputxyz(gridsize, array, missval, nlon, nlat, grid_center_lon, grid_center_lat);
                }
              else if (operatorID == OUTPUTARR) { outputarr(dig, gridsize, array); }
              else
                {
                  if (gridInqType(gridID) == GRID_SPECTRAL && gridsize <= 156)
                    outputsp(gridsize, array, gridInqTrunc(gridID));
                  else
                    output(gridsize, array);
                }
            }

          tsID++;
        }

      cdo_stream_close(streamID);
    }

  cdo_finish();

  return nullptr;
}
