/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

*/

#include <cdi.h>

#include "process_int.h"
#include "cdo_math.h"
#include "cdo_vlist.h"
#include "progress.h"
#include "cdo_options.h"
#include "grid_healpix.h"
#include "param_conversion.h"
#include "pmlist.h"
#include "mpim_grid.h"

enum class Stat
{
  Mean = 1,
  Avg = 2
};

struct HealpixParams
{
  int fact = 1;
  int nsideIn = 0;
  int nsideOut = 0;
  HpOrder orderIn = HpOrder::Undef;
  HpOrder orderOut = HpOrder::Undef;
  Stat stat = Stat::Mean;
  double power = 0.0;
  bool doDegrade = true;
};

template <typename T>
static T
stat_avg_mv(const T *v, size_t n, T missval, double scale)
{
  double sum = 0.0;
  size_t nOut = 0;
  for (size_t i = 0; i < n; ++i)
    if (!dbl_is_equal(v[i], missval))
      {
        sum += v[i];
        nOut++;
      }

  return (nOut == n) ? (sum / nOut) * scale : missval;
}

template <typename T>
static T
stat_mean_mv(const T *v, size_t n, T missval, double scale)
{
  double sum = 0.0;
  size_t nOut = 0;
  for (size_t i = 0; i < n; ++i)
    if (!dbl_is_equal(v[i], missval))
      {
        sum += v[i];
        nOut++;
      }

  return (nOut > 0) ? (sum / nOut) * scale : missval;
}

template <typename T>
static T
stat_mean(const T *v, size_t n)
{
  double sum = 0.0;
  for (size_t i = 0; i < n; ++i) sum += v[i];
  return sum / n;
}

static double
get_scalefactor(const HealpixParams &params)
{
  double nsideQuot = static_cast<double>(params.nsideIn) / params.nsideOut;
  return (std::fabs(params.power) > 0.0) ? std::pow(nsideQuot, -params.power) : 1.0;
}

template <typename T>
static void
degrade(const Varray<T> &v1, size_t gridsize2, Varray<T> &v2, bool hasMissvals, T missval, const HealpixParams &params)
{
  auto scale = get_scalefactor(params);
  size_t nvals = params.fact * params.fact;
  if (hasMissvals)
    {
      if (params.stat == Stat::Mean)
        for (size_t i = 0; i < gridsize2; ++i) v2[i] = stat_mean_mv(&v1[i * nvals], nvals, missval, scale);
      else
        for (size_t i = 0; i < gridsize2; ++i) v2[i] = stat_avg_mv(&v1[i * nvals], nvals, missval, scale);
    }
  else
    {
      for (size_t i = 0; i < gridsize2; ++i) v2[i] = stat_mean(&v1[i * nvals], nvals) * scale;
    }
}

static void
hp_degrade(const Field &field1, Field &field2, const HealpixParams &params)
{
  auto hasMissvals = (field1.nmiss > 0);
  if (field1.memType == MemType::Float)
    degrade(field1.vec_f, field2.gridsize, field2.vec_f, hasMissvals, (float) field1.missval, params);
  else
    degrade(field1.vec_d, field2.gridsize, field2.vec_d, hasMissvals, field1.missval, params);

  if (hasMissvals) field_num_mv(field2);
}

template <typename T>
static void
upgrade(size_t gridsize1, const Varray<T> &v1, Varray<T> &v2, bool hasMissvals, T missval, const HealpixParams &params)
{
  auto scale = get_scalefactor(params);
  size_t nvals = params.fact * params.fact;
  if (hasMissvals)
    {
      for (size_t i = 0; i < gridsize1; ++i)
        for (size_t k = 0; k < nvals; ++k) v2[i * nvals + k] = dbl_is_equal(v1[i], missval) ? missval : v1[i] * scale;
    }
  else
    {
      for (size_t i = 0; i < gridsize1; ++i)
        for (size_t k = 0; k < nvals; ++k) v2[i * nvals + k] = v1[i] * scale;
    }
}

static void
hp_upgrade(const Field &field1, Field &field2, const HealpixParams &params)
{
  auto hasMissvals = (field1.nmiss > 0);
  if (field1.memType == MemType::Float)
    upgrade(field1.gridsize, field1.vec_f, field2.vec_f, hasMissvals, (float) field1.missval, params);
  else
    upgrade(field1.gridsize, field1.vec_d, field2.vec_d, hasMissvals, field1.missval, params);

  if (hasMissvals) field_num_mv(field2);
}

template <typename T>
static void
ring_to_nested(int nside, size_t gridsize, Varray<T> &v)
{
  Varray<T> vtmp = v;
  hp_ring_to_nested(nside, gridsize, vtmp.data(), v.data());
}

static void
ring_to_nested(Field &field, int nside)
{
  if (field.memType == MemType::Float)
    ring_to_nested(nside, field.gridsize, field.vec_f);
  else
    ring_to_nested(nside, field.gridsize, field.vec_d);
}

template <typename T>
static void
nested_to_ring(int nside, size_t gridsize, Varray<T> &v)
{
  Varray<T> vtmp = v;
  hp_nested_to_ring(nside, gridsize, vtmp.data(), v.data());
}

static void
nested_to_ring(Field &field, int nside)
{
  if (field.memType == MemType::Float)
    nested_to_ring(nside, field.gridsize, field.vec_f);
  else
    nested_to_ring(nside, field.gridsize, field.vec_d);
}

static Stat
set_stat(const std::string &statString)
{
  if (statString == "mean")
    return Stat::Mean;
  else if (statString == "avg")
    return Stat::Avg;
  else
    cdo_abort("Parameter value stat=%s unsupported!", statString);

  return Stat::Mean;
}

static HealpixParams
get_parameter(void)
{
  HealpixParams params;

  auto pargc = cdo_operator_argc();
  if (pargc)
    {
      auto &pargv = cdo_get_oper_argv();

      KVList kvlist;
      kvlist.name = cdo_module_name();
      if (kvlist.parse_arguments(pargc, pargv) != 0) cdo_abort("Parse error!");
      if (Options::cdoVerbose) kvlist.print();

      for (const auto &kv : kvlist)
        {
          const auto &key = kv.key;
          if (kv.nvalues > 1) cdo_abort("Too many values for parameter key >%s<!", key);
          if (kv.nvalues < 1) cdo_abort("Missing value for parameter key >%s<!", key);
          const auto &value = kv.values[0];

          // clang-format off
          if      (key == "nside") params.nsideOut = parameter_to_int(value);
          else if (key == "order") params.orderOut = hp_get_order(parameter_to_word(value));
          else if (key == "fact")  params.fact     = parameter_to_int(value);
          else if (key == "stat")  params.stat     = set_stat(parameter_to_word(value));
          else if (key == "power") params.power    = parameter_to_double(value);
          else cdo_abort("Invalid parameter key >%s<!", key);
          // clang-format on
        }
    }

  return params;
}

static void
verify_parameter(const HealpixParams &params)
{
  if (params.fact > 1 && params.nsideOut > 0) cdo_abort("Parameter 'fact' can't be combined with 'nside'!");
}

static int
define_healpix_grid(size_t gridsize, int nside, HpOrder order)
{
  auto orderString = (order == HpOrder::Ring) ? "ring" : "nested";
  auto projection = "healpix";
  auto gridID = gridCreate(GRID_PROJECTION, gridsize);
  cdiDefKeyString(gridID, CDI_GLOBAL, CDI_KEY_DIMNAME, "cells");
  cdiDefKeyString(gridID, CDI_GLOBAL, CDI_KEY_GRIDMAP_VARNAME, projection);
  cdiDefKeyString(gridID, CDI_GLOBAL, CDI_KEY_GRIDMAP_NAME, projection);
  cdiDefAttTxt(gridID, CDI_GLOBAL, "grid_mapping_name", (int) strlen(projection), projection);
  cdiDefAttInt(gridID, CDI_GLOBAL, "healpix_nside", CDI_DATATYPE_INT32, 1, &nside);
  cdiDefAttTxt(gridID, CDI_GLOBAL, "healpix_order", (int) strlen(orderString), orderString);

  return gridID;
}

static int
hp_define_grid(int gridID1, HealpixParams &params)
{
  int gridIDout = -1;

  auto [nside, order] = cdo::get_healpix_params(gridID1);
  params.nsideIn = nside;
  params.orderIn = order;

  if (!cdo::is_power_of_two(params.nsideIn)) cdo_abort("Input healpix: nside must be a power of two!");

  if (params.nsideOut == 0)
    {
      auto fact = params.fact;
      params.nsideOut = (fact > 1) ? (params.doDegrade ? nside / fact : nside * fact) : nside;
    }
  else
    {
      if (params.doDegrade)
        {
          if (params.nsideOut > params.nsideIn) cdo_abort("Parameter nside must be less than input nside=%d!", params.nsideIn);
          params.fact = params.nsideIn / params.nsideOut;
        }
      else
        {
          if (params.nsideOut < params.nsideIn) cdo_abort("Parameter nside must be greater than input nside=%d!", params.nsideIn);
          params.fact = params.nsideOut / params.nsideIn;
        }
    }

  if (!cdo::is_power_of_two(params.nsideOut)) cdo_abort("Parameter nside must be a power of two!");

  if (params.orderOut == HpOrder::Undef) params.orderOut = params.orderIn;

  size_t gridsize = 12 * params.nsideOut * params.nsideOut;
  gridIDout = define_healpix_grid(gridsize, params.nsideOut, params.orderOut);

  return gridIDout;
}

void *
Healpix(void *process)
{
  cdo_initialize(process);

  // clang-format off
                   cdo_operator_add("hpupgrade", 0, 0, nullptr);
  auto HPDEGRADE = cdo_operator_add("hpdegrade", 0, 0, nullptr);
  // clang-format on

  auto operatorID = cdo_operator_id();
  auto doDegrade = (operatorID == HPDEGRADE);

  auto params = get_parameter();
  params.doDegrade = doDegrade;
  verify_parameter(params);

  auto streamID1 = cdo_open_read(0);
  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto taxisID1 = vlistInqTaxis(vlistID1);

  auto ngrids = vlistNgrids(vlistID1);
  if (ngrids > 1) cdo_abort("Too many different grids!");

  auto gridID = vlistGrid(vlistID1, 0);
  if (!is_healpix_grid(gridID)) cdo_abort("Input grid is not healpix!");

  auto vlistID2 = vlistDuplicate(vlistID1);
  auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  auto gridID2 = hp_define_grid(gridID, params);
  for (int index = 0; index < ngrids; ++index) vlistChangeGridIndex(vlistID2, index, gridID2);

  VarList varList1, varList2;
  varListInit(varList1, vlistID1);
  varListInit(varList2, vlistID2);

  auto streamID2 = cdo_open_write(1);
  cdo_def_vlist(streamID2, vlistID2);

  Field field1, field2;

  int tsID1 = 0;
  while (true)
    {
      auto nrecs = cdo_stream_inq_timestep(streamID1, tsID1);
      if (nrecs == 0) break;

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_def_timestep(streamID2, tsID1);

      for (int recID = 0; recID < nrecs; ++recID)
        {
          int varID, levelID;
          cdo_inq_record(streamID1, &varID, &levelID);
          field1.init(varList1[varID]);
          cdo_read_record(streamID1, field1);

          field2.init(varList2[varID]);

          if (params.orderIn == HpOrder::Ring) ring_to_nested(field1, params.nsideIn);

          doDegrade ? hp_degrade(field1, field2, params) : hp_upgrade(field1, field2, params);

          if (params.orderOut == HpOrder::Ring) nested_to_ring(field2, params.nsideOut);

          cdo_def_record(streamID2, varID, levelID);
          cdo_write_record(streamID2, field2);
        }

      tsID1++;
    }

  cdo_stream_close(streamID1);
  cdo_stream_close(streamID2);

  vlistDestroy(vlistID2);

  cdo_finish();

  return nullptr;
}
