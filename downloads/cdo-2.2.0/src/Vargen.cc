/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida
          Ralf Müller

*/

/*
   This module contains the following operators:

      Vargen     const           Create a constant field
      Vargen     random          Field with random values
      Vargen     stdatm          Field values for pressure and temperature for the standard atmosphere
*/

#include <cstdlib>
#include <cassert>
#include <cdi.h>
#include "julian_date.h"

#include "cdo_options.h"
#include "process_int.h"
#include "cdo_zaxis.h"
#include "param_conversion.h"
#include <mpim_grid.h>
#include "grid_healpix.h"
#include "griddes.h"
#include "constants.h"
#include "stdnametable.h"
#include "param_conversion.h"

static constexpr double etopo_scale = 3.0;
static constexpr double etopo_offset = 11000.0;
static constexpr unsigned short etopo[] = {
#include "etopo.dat"
};

static constexpr double temp_scale = 500.0;
static constexpr double temp_offset = -220.0;
static constexpr unsigned short temp[] = {
#include "temp.dat"
};

static constexpr double mask_scale = 1.0;
static constexpr double mask_offset = 0.0;
static constexpr unsigned short mask[] = {
#include "mask.dat"
};

// Some Constants for creating temperatur and pressure for the standard atmosphere
constexpr double T_ZERO = 213.0;
constexpr double T_DELTA = 75.0;
constexpr double SCALEHEIGHT = 10000.0;  // [m]

static double
std_atm_temperatur(double height)
{
  // Compute the temperatur for the given height (in meters) according to the solution of the hydrostatic atmosphere
  return (T_ZERO + T_DELTA * std::exp((-1) * (height / SCALEHEIGHT)));
}

static double
std_atm_pressure(double height)
{
  constexpr double P_ZERO = 1013.25;  // surface pressure [hPa]
  constexpr double CC_R = 287.05;     // specific gas constant for air
  constexpr double TMP4PRESSURE = (C_EARTH_GRAV * SCALEHEIGHT) / (CC_R * T_ZERO);

  // Compute the pressure for the given height (in meters) according to the solution of the hydrostatic atmosphere
  return (P_ZERO
          * std::exp((-1) * TMP4PRESSURE * std::log((std::exp(height / SCALEHEIGHT) * T_ZERO + T_DELTA) / (T_ZERO + T_DELTA))));
}

static void
conv_generic_grid(int gridID, size_t gridsize, Varray<double> &xvals2D, Varray<double> &yvals2D)
{
  auto xsize = gridInqXsize(gridID);
  auto ysize = gridInqYsize(gridID);

  assert(gridsize == xsize * ysize);

  Varray<double> xcoord(xsize), ycoord(ysize);
  gridInqXvals(gridID, &xcoord[0]);
  gridInqYvals(gridID, &ycoord[0]);

  auto xrange = varray_range(xsize, xcoord);
  auto yrange = varray_range(ysize, ycoord);

  for (size_t j = 0; j < ysize; ++j)
    for (size_t i = 0; i < xsize; ++i)
      {
        xvals2D[j * xsize + i] = xcoord[i] * M_PI / xrange;
        yvals2D[j * xsize + i] = ycoord[j] * M_PI / yrange;
      }
}

static size_t
calc_index_ii(size_t nx, double xval)
{
  if (xval >= 180.0) xval -= 360.0;
  if (xval < -180.0) xval += 360.0;
  size_t ii = (xval + 180.0) * 2.0;
  if (ii >= nx) ii = nx - 1;
  return ii;
}

static size_t
calc_index_jj(size_t ny, double yval)
{
  size_t jj = (yval + 90.0) * 2.0;
  if (jj >= ny) jj = ny - 1;
  return jj;
}

static void
remap_nn_reg2d_to_reg2d(size_t nx, size_t ny, const Varray<float> &data, int gridID, Varray<float> &array)
{
  auto gridtype = gridInqType(gridID);
  if (gridtype != GRID_LONLAT && gridtype != GRID_GAUSSIAN) cdo_abort("Internal error, wrong grid type!");

  auto nxvals = gridInqXsize(gridID);
  auto nyvals = gridInqYsize(gridID);
  Varray<double> xvals(nxvals), yvals(nyvals);

  gridInqXvals(gridID, xvals.data());
  gridInqYvals(gridID, yvals.data());

  // Convert lat/lon units if required
  cdo_grid_to_degree(gridID, CDI_XAXIS, nxvals, xvals.data(), "grid center lon");
  cdo_grid_to_degree(gridID, CDI_YAXIS, nyvals, yvals.data(), "grid center lat");

#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
  for (size_t j = 0; j < nyvals; ++j)
    {
      auto jj = calc_index_jj(ny, yvals[j]);
      for (size_t i = 0; i < nxvals; ++i)
        {
          auto ii = calc_index_ii(nx, xvals[i]);
          array[j * nxvals + i] = data[jj * nx + ii];
        }
    }
}

static void
remap_nn_reg2d_to_nonreg2d(size_t nx, size_t ny, const Varray<float> &data, int gridID, Varray<float> &array)
{
  auto gridsize = gridInqSize(gridID);
  Varray<double> xvals(gridsize), yvals(gridsize);

  auto gridID2 = generate_full_point_grid(gridID);
  if (!gridHasCoordinates(gridID2)) cdo_abort("Target cell center coordinates missing!");

  gridInqXvals(gridID2, xvals.data());
  gridInqYvals(gridID2, yvals.data());

  // Convert lat/lon units if required
  cdo_grid_to_degree(gridID2, CDI_XAXIS, gridsize, xvals.data(), "grid center lon");
  cdo_grid_to_degree(gridID2, CDI_YAXIS, gridsize, yvals.data(), "grid center lat");

#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
  for (size_t i = 0; i < gridsize; ++i)
    {
      auto jj = calc_index_jj(ny, yvals[i]);
      auto ii = calc_index_ii(nx, xvals[i]);
      array[i] = data[jj * nx + ii];
    }

  if (gridID != gridID2) gridDestroy(gridID2);
}

static void
remap_nn_reg2d_to_healpix(size_t nx, size_t ny, const Varray<float> &data, int gridID, Varray<float> &array)
{
  auto gridID2 = gridID;
  auto gridsize = gridInqSize(gridID2);
  auto [r_nside, r_order] = cdo::get_healpix_params(gridID2);
  int nside = r_nside;
  HpOrder order = r_order;

#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
  for (size_t i = 0; i < gridsize; ++i)
    {
      double xval, yval;
      hp_index_to_lonlat(order, nside, i, &xval, &yval);
      auto jj = calc_index_jj(ny, yval * RAD2DEG);
      auto ii = calc_index_ii(nx, xval * RAD2DEG);
      array[i] = data[jj * nx + ii];
    }
}

static void
remap_nn_reg2d(size_t nx, size_t ny, const Varray<float> &data, int gridID, Varray<float> &array)
{
  auto gridtype = gridInqType(gridID);
  if (gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN)
    remap_nn_reg2d_to_reg2d(nx, ny, data, gridID, array);
  else if (is_healpix_grid(gridID))
    remap_nn_reg2d_to_healpix(nx, ny, data, gridID, array);
  else
    remap_nn_reg2d_to_nonreg2d(nx, ny, data, gridID, array);
}

static int
random_init(int operatorID)
{
  unsigned int seed = Options::Random_Seed;
  operator_input_arg(cdo_operator_enter(operatorID));
  if (cdo_operator_argc() < 1) cdo_abort("Too few arguments!");
  if (cdo_operator_argc() > 2) cdo_abort("Too many arguments!");
  auto gridID = cdo_define_grid(cdo_operator_argv(0));
  if (cdo_operator_argc() == 2)
    {
      auto idum = parameter_to_int(cdo_operator_argv(1));
      if (idum >= 0 && idum < 0x7FFFFFFF) seed = idum;
    }
  std::srand(seed);
  return gridID;
}

static void
random_compute(size_t gridsize, Varray<float> &array)
{
  for (size_t i = 0; i < gridsize; ++i) array[i] = ((double) std::rand()) / ((double) RAND_MAX);
}

static void
sincos_compute(size_t gridsize, Varray<float> &array, const Varray<double> &xvals, const Varray<double> &yvals)
{
  for (size_t i = 0; i < gridsize; ++i) array[i] = std::cos(1.0 * xvals[i]) * std::sin(2.0 * yvals[i]);
}

static void
coshill_compute(size_t gridsize, Varray<float> &array, const Varray<double> &xvals, const Varray<double> &yvals)
{
  for (size_t i = 0; i < gridsize; ++i) array[i] = 2.0 - std::cos(std::acos(std::cos(xvals[i]) * std::cos(yvals[i])) / 1.2);
}

static void
testfield_compute(size_t gridsize, Varray<float> &array, const Varray<double> &xvals, const Varray<double> &yvals)
{
  double xyz[3];
  for (size_t i = 0; i < gridsize; ++i)
    {
      gcLLtoXYZ(xvals[i], yvals[i], xyz);
      auto x = xyz[0];
      auto y = xyz[1];
      auto z = xyz[2];
      array[i] = 1.0 + std::pow(x, 8.0) + std::exp(2.0 * y * y * y) + std::exp(2.0 * x * x) + 10.0 * x * y * z;
    }
}

static void
unpack_data(size_t datasize, Varray<float> &data, double scale, double offset, const unsigned short *zdata)
{
  for (size_t i = 0; i < datasize; ++i) data[i] = zdata[i] / scale - offset;
}

static int
define_point_grid()
{
  auto gridID = gridCreate(GRID_LONLAT, 1);
  gridDefXsize(gridID, 1);
  gridDefYsize(gridID, 1);
  double value = 0.0;
  gridDefXvals(gridID, &value);
  gridDefYvals(gridID, &value);

  return gridID;
}

static int
define_zaxis(bool lstdatm, int nlevels, double *levels)
{
  int zaxisID = -1;

  if (lstdatm)
    {
      zaxisID = zaxisCreate(ZAXIS_HEIGHT, nlevels);
      zaxisDefLevels(zaxisID, levels);
      cdiDefKeyString(zaxisID, CDI_GLOBAL, CDI_KEY_NAME, "level");
      cdiDefKeyString(zaxisID, CDI_GLOBAL, CDI_KEY_LONGNAME, "Level");
      cdiDefKeyString(zaxisID, CDI_GLOBAL, CDI_KEY_UNITS, "m");
    }
  else { zaxisID = zaxis_from_name("surface"); }

  return zaxisID;
}

static void
define_pressure_attributes(int vlistID, int varID)
{
  vlistDefVarParam(vlistID, varID, cdiEncodeParam(1, 255, 255));
  cdiDefKeyString(vlistID, varID, CDI_KEY_NAME, "P");
  cdiDefKeyString(vlistID, varID, CDI_KEY_STDNAME, "air_pressure");
  cdiDefKeyString(vlistID, varID, CDI_KEY_LONGNAME, "pressure");
  cdiDefKeyString(vlistID, varID, CDI_KEY_UNITS, "hPa");
}

static void
define_temperature_attributes(int vlistID, int varID)
{
  vlistDefVarParam(vlistID, varID, cdiEncodeParam(130, 128, 255));
  cdiDefKeyString(vlistID, varID, CDI_KEY_NAME, "T");
  cdiDefKeyString(vlistID, varID, CDI_KEY_STDNAME, var_stdname(air_temperature));
  cdiDefKeyString(vlistID, varID, CDI_KEY_LONGNAME, "temperature");
  cdiDefKeyString(vlistID, varID, CDI_KEY_UNITS, "K");
}

class ModuleVargen
{
private:
  int RANDOM;
  int SINCOS;
  int COSHILL;
  int TESTFIELD;
  int CONST;
  int SEQ;
  int TOPO;
  int TEMP;
  int MASK;
  int STDATM;

  static constexpr size_t nlat = 360, nlon = 720;
  double lon[nlon], lat[nlat];
  int nlevels = 1;
  int gridID = -1, gridIDdata = -1;
  double rstart = 0.0, rstop = 0.0, rinc = 0.0;
  double rconst = 0.0;
  std::vector<double> levels;
  int ntimesteps;
  int julday;

  int taxisID;
  CdoStreamID streamID;
  int varID;
  int varID2;
  int vlistID;
  int operatorID;
  size_t gridsize;
  Varray<float> array;
  int nvars;

public:
  void
  init(void *process)
  {
    cdo_initialize(process);

    // clang-format off
    RANDOM    = cdo_operator_add("random",    0, 0, "grid description file or name, <seed>");
    SINCOS    = cdo_operator_add("sincos",    0, 0, "grid description file or name");
    COSHILL   = cdo_operator_add("coshill",   0, 0, "grid description file or name");
    //not usesd todo: make unavailable for non developers
    TESTFIELD = cdo_operator_add("testfield", 0, 0, "grid description file or name");
    CONST     = cdo_operator_add("const",     0, 0, "constant value, grid description file or name");
    SEQ       = cdo_operator_add("seq",       0, 0, "start, end, <increment>");
    TOPO      = cdo_operator_add("topo",      0, 0, nullptr);
    TEMP      = cdo_operator_add("temp",      0, 0, nullptr);
    MASK      = cdo_operator_add("mask",      0, 0, nullptr);
    STDATM    = cdo_operator_add("stdatm",    0, 0, "height levels [m]");
    // clang-format on

    operatorID = cdo_operator_id();

    if (operatorID == RANDOM) { gridID = random_init(operatorID); }
    else if (operatorID == SINCOS || operatorID == COSHILL || operatorID == TESTFIELD)
      {
        operator_input_arg(cdo_operator_enter(operatorID));
        operator_check_argc(1);
        gridID = cdo_define_grid(cdo_operator_argv(0));
      }
    else if (operatorID == CONST)
      {
        operator_input_arg(cdo_operator_enter(operatorID));
        operator_check_argc(2);
        rconst = parameter_to_double(cdo_operator_argv(0));
        gridID = cdo_define_grid(cdo_operator_argv(1));
      }
    else if (operatorID == TOPO || operatorID == TEMP || operatorID == MASK)
      {
        gridIDdata = gridCreate(GRID_LONLAT, nlon * nlat);
        gridDefXsize(gridIDdata, nlon);
        gridDefYsize(gridIDdata, nlat);

        for (size_t i = 0; i < nlon; ++i) lon[i] = -179.75 + i * 0.5;
        for (size_t i = 0; i < nlat; ++i) lat[i] = -89.75 + i * 0.5;

        gridDefXvals(gridIDdata, lon);
        gridDefYvals(gridIDdata, lat);

        gridID = gridIDdata;

        if (cdo_operator_argc() == 1) gridID = cdo_define_grid(cdo_operator_argv(0));
        if (cdo_operator_argc() > 1) cdo_abort("Too many arguments!");
      }
    else if (operatorID == SEQ)
      {
        operator_input_arg(cdo_operator_enter(operatorID));
        if (cdo_operator_argc() < 2) cdo_abort("Too few arguments!");
        if (cdo_operator_argc() > 3) cdo_abort("Too many arguments!");

        rstart = parameter_to_double(cdo_operator_argv(0));
        rstop = parameter_to_double(cdo_operator_argv(1));
        rinc = (cdo_operator_argc() == 3) ? parameter_to_double(cdo_operator_argv(2)) : 1;
        if (DBL_IS_EQUAL(rinc, 0.0)) cdo_abort("Increment is zero!");

        gridID = define_point_grid();
      }
    else if (operatorID == STDATM)
      {
        operator_input_arg(cdo_operator_enter(operatorID));
        levels = cdo_argv_to_flt(cdo_get_oper_argv());
        nlevels = levels.size();

        if (Options::cdoVerbose)
          for (int i = 0; i < nlevels; ++i) printf("levels %d: %g\n", i, levels[i]);

        gridID = define_point_grid();
      }

    auto zaxisID = define_zaxis(operatorID == STDATM, nlevels, levels.data());

    vlistID = vlistCreate();

    auto timetype = (operatorID == SEQ) ? TIME_VARYING : TIME_CONSTANT;

    varID = vlistDefVar(vlistID, gridID, zaxisID, timetype);
    /*
       For the standard atmosphere two output variables are generated: pressure and temperature.
       The first (varID) is pressure, second (varID2) is temperature. Add an additional variable for the standard atmosphere.
    */
    varID2 = (operatorID == STDATM) ? vlistDefVar(vlistID, gridID, zaxisID, TIME_CONSTANT) : -1;

    if (operatorID == MASK) vlistDefVarDatatype(vlistID, varID, CDI_DATATYPE_INT8);

    if (operatorID == STDATM)
      {
        define_pressure_attributes(vlistID, varID);
        define_temperature_attributes(vlistID, varID2);
      }
    else
      {
        cdiDefKeyString(vlistID, varID, CDI_KEY_NAME, cdo_operator_name(operatorID));
        if (operatorID == TOPO) cdiDefKeyString(vlistID, varID, CDI_KEY_UNITS, "m");
        if (operatorID == TEMP) cdiDefKeyString(vlistID, varID, CDI_KEY_UNITS, "K");
      }

    taxisID = cdo_taxis_create(TAXIS_RELATIVE);
    vlistDefTaxis(vlistID, taxisID);

    if (operatorID != SEQ) vlistDefNtsteps(vlistID, 1);

    streamID = cdo_open_write(0);

    cdo_def_vlist(streamID, vlistID);

    gridsize = gridInqSize(gridID);
    array = Varray<float>(gridsize);

    ntimesteps = (operatorID == SEQ) ? 1.001 + ((rstop - rstart) / rinc) : 1;
    if (operatorID != SEQ) vlistDefNtsteps(vlistID, 0);

    julday = date_to_julday(CALENDAR_PROLEPTIC, 10101);

    nvars = vlistNvars(vlistID);
  }

  void
  run()
  {
    for (int tsID = 0; tsID < ntimesteps; ++tsID)
      {
        auto rval = rstart + rinc * tsID;
        CdiDateTime vDateTime{};
        vDateTime.date = cdiDate_set(julday_to_date(CALENDAR_PROLEPTIC, julday + tsID));
        taxisDefVdatetime(taxisID, vDateTime);
        cdo_def_timestep(streamID, tsID);

        // this should either be 1 or 2, two for atmosphere
        for (varID = 0; varID < nvars; ++varID)
          {
            nlevels = zaxisInqSize(vlistInqVarZaxis(vlistID, varID));
            for (int levelID = 0; levelID < nlevels; ++levelID)
              {
                cdo_def_record(streamID, varID, levelID);

                if (operatorID == RANDOM) { random_compute(gridsize, array); }
                else if (operatorID == SINCOS || operatorID == COSHILL || operatorID == TESTFIELD)
                  {
                    Varray<double> xvals(gridsize), yvals(gridsize);

                    if (grid_is_distance_generic(gridID)) { conv_generic_grid(gridID, gridsize, xvals, yvals); }
                    else
                      {
                        gridID = generate_full_point_grid(gridID);
                        if (!gridHasCoordinates(gridID)) cdo_abort("Target cell center coordinates missing!");

                        gridInqXvals(gridID, xvals.data());
                        gridInqYvals(gridID, yvals.data());

                        // Convert lat/lon units if required
                        cdo_grid_to_radian(gridID, CDI_XAXIS, gridsize, xvals.data(), "grid center lon");
                        cdo_grid_to_radian(gridID, CDI_YAXIS, gridsize, yvals.data(), "grid center lat");
                      }

                    // clang-format off
                    if      (operatorID == SINCOS)    sincos_compute(gridsize, array, xvals, yvals);
                    else if (operatorID == COSHILL)   coshill_compute(gridsize, array, xvals, yvals);
                    else if (operatorID == TESTFIELD) testfield_compute(gridsize, array, xvals, yvals);
                    // clang-format on
                  }
                else if (operatorID == CONST)
                  {
                    for (size_t i = 0; i < gridsize; ++i) array[i] = rconst;
                  }
                else if (operatorID == TOPO || operatorID == TEMP || operatorID == MASK)
                  {
                    auto datasize = gridInqSize(gridIDdata);
                    Varray<float> data(datasize);

                    // clang-format off
                    if      (operatorID == TOPO) unpack_data(datasize, data, etopo_scale, etopo_offset, etopo);
                    else if (operatorID == TEMP) unpack_data(datasize, data, temp_scale, temp_offset, temp);
                    else if (operatorID == MASK) unpack_data(datasize, data, mask_scale, mask_offset, mask);
                    // clang-format on

                    if (gridID != gridIDdata && gridIDdata != -1) { remap_nn_reg2d(nlon, nlat, data, gridID, array); }
                    else
                      {
                        for (size_t i = 0; i < gridsize; ++i) array[i] = data[i];
                      }
                  }
                else if (operatorID == SEQ) { array[0] = rval; }
                else if (operatorID == STDATM)
                  {
                    array[0] = (varID == varID2) ? std_atm_temperatur(levels[levelID]) : std_atm_pressure(levels[levelID]);
                  }

                cdo_write_record_f(streamID, array.data(), 0);
              }
          }
      }
  }

  void
  close()
  {
    cdo_stream_close(streamID);

    vlistDestroy(vlistID);

    cdo_finish();
  }
};

void *
Vargen(void *process)
{
  ModuleVargen vargen;
  vargen.init(process);
  vargen.run();
  vargen.close();
  return nullptr;
}
