/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Setmiss    setmissval      Set a new missing value
      Setmiss    setctomiss      Set constant to missing value
      Setmiss    setmisstoc      Set missing value to constant
      Setmiss    setrtomiss      Set range to missing value
      Setmiss    setvrange       Set range of valid value
*/

#include <cmath>
#include <cdi.h>

#include "process_int.h"
#include "param_conversion.h"

template <typename T>
static size_t
set_missval(size_t gridsize, Varray<T> &array, T missval, T new_missval)
{
  size_t nmiss = 0;
  for (size_t i = 0; i < gridsize; ++i)
    if (DBL_IS_EQUAL(array[i], missval) || DBL_IS_EQUAL(array[i], (float) missval) || DBL_IS_EQUAL(array[i], new_missval)
        || DBL_IS_EQUAL(array[i], (float) new_missval))
      {
        array[i] = new_missval;
        nmiss++;
      }

  return nmiss;
}

static void
set_missval(Field &field, double new_missval)
{
  if (field.memType == MemType::Float)
    field.nmiss = set_missval(field.size, field.vec_f, (float) field.missval, (float) new_missval);
  else
    field.nmiss = set_missval(field.size, field.vec_d, field.missval, new_missval);
}

template <typename T>
static size_t
set_const_to_miss(size_t gridsize, Varray<T> &array, T missval, T rconst)
{
  size_t nmiss = 0;
  if (std::isnan(rconst))
    {
      for (size_t i = 0; i < gridsize; ++i)
        if (std::isnan(array[i]))
          {
            array[i] = missval;
            nmiss++;
          }
    }
  else
    {
      for (size_t i = 0; i < gridsize; ++i)
        if (DBL_IS_EQUAL(array[i], rconst) || DBL_IS_EQUAL(array[i], (float) rconst))
          {
            array[i] = missval;
            nmiss++;
          }
    }

  return nmiss;
}

static void
set_const_to_miss(Field &field, double rconst)
{
  if (field.memType == MemType::Float)
    field.nmiss += set_const_to_miss(field.size, field.vec_f, (float) field.missval, (float) rconst);
  else
    field.nmiss += set_const_to_miss(field.size, field.vec_d, field.missval, rconst);
}

template <typename T>
static size_t
set_miss_to_const(size_t gridsize, Varray<T> &array, T missval, T rconst)
{
  for (size_t i = 0; i < gridsize; ++i)
    if (DBL_IS_EQUAL(array[i], missval) || DBL_IS_EQUAL(array[i], (float) missval)) { array[i] = rconst; }

  return 0;
}

static void
set_miss_to_const(Field &field, double rconst)
{
  if (field.memType == MemType::Float)
    field.nmiss = set_miss_to_const(field.size, field.vec_f, (float) field.missval, (float) rconst);
  else
    field.nmiss = set_miss_to_const(field.size, field.vec_d, field.missval, rconst);
}

template <typename T>
static size_t
set_range_to_miss(size_t gridsize, Varray<T> &array, T missval, T rmin, T rmax)
{
  size_t nmiss = 0;
  for (size_t i = 0; i < gridsize; ++i)
    if (array[i] >= rmin && array[i] <= rmax)
      {
        array[i] = missval;
        nmiss++;
      }

  return nmiss;
}

static void
set_range_to_miss(Field &field, double rmin, double rmax)
{
  if (field.memType == MemType::Float)
    field.nmiss += set_range_to_miss(field.size, field.vec_f, (float) field.missval, (float) rmin, (float) rmax);
  else
    field.nmiss += set_range_to_miss(field.size, field.vec_d, field.missval, rmin, rmax);
}

template <typename T>
static size_t
set_valid_range(size_t gridsize, Varray<T> &array, T missval, T rmin, T rmax)
{
  for (size_t i = 0; i < gridsize; ++i)
    if (array[i] < rmin || array[i] > rmax) array[i] = missval;

  const auto nmiss = varray_num_mv(gridsize, array, missval);

  return nmiss;
}

static void
set_valid_range(Field &field, double rmin, double rmax)
{
  if (field.memType == MemType::Float)
    field.nmiss = set_valid_range(field.size, field.vec_f, (float) field.missval, (float) rmin, (float) rmax);
  else
    field.nmiss = set_valid_range(field.size, field.vec_d, field.missval, rmin, rmax);
}

void *
Setmiss(void *process)
{
  cdo_initialize(process);

  // clang-format off
  const auto SETMISSVAL = cdo_operator_add("setmissval", 0, 0, "missing value");
  const auto SETCTOMISS = cdo_operator_add("setctomiss", 0, 0, "constant");
  const auto SETMISSTOC = cdo_operator_add("setmisstoc", 0, 0, "constant");
  const auto SETRTOMISS = cdo_operator_add("setrtomiss", 0, 0, "range (min, max)");
  const auto SETVRANGE  = cdo_operator_add("setvrange",  0, 0, "range (min, max)");
  // clang-format on

  const auto operatorID = cdo_operator_id();

  double new_missval = 0.0;
  double rconst = 0.0, rmin = 0.0, rmax = 0.0;

  if (operatorID == SETMISSVAL)
    {
      operator_check_argc(1);
      new_missval = parameter_to_double(cdo_operator_argv(0));
    }
  else if (operatorID == SETCTOMISS || operatorID == SETMISSTOC)
    {
      operator_check_argc(1);
      rconst = parameter_to_double(cdo_operator_argv(0));
    }
  else
    {
      operator_check_argc(2);
      rmin = parameter_to_double(cdo_operator_argv(0));
      rmax = parameter_to_double(cdo_operator_argv(1));
    }

  const auto streamID1 = cdo_open_read(0);

  const auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  const auto vlistID2 = vlistDuplicate(vlistID1);

  const auto taxisID1 = vlistInqTaxis(vlistID1);
  const auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  VarList varList;
  varListInit(varList, vlistID1);

  Field field;

  if (operatorID == SETMISSVAL)
    {
      const auto nvars = vlistNvars(vlistID2);
      for (int varID = 0; varID < nvars; ++varID) vlistDefVarMissval(vlistID2, varID, new_missval);
    }
  else if (operatorID == SETMISSTOC)
    {
      auto nvars = vlistNvars(vlistID2);
      for (int varID = 0; varID < nvars; ++varID)
        {
          const auto &var = varList[varID];
          if (DBL_IS_EQUAL(rconst, var.missval))
            {
              cdo_warning("Missing value and constant have the same value!");
              break;
            }
        }
    }

  /*
  if (operatorID == SETVRANGE)
    {
      double range[2] = {rmin, rmax};
      nvars = vlistNvars(vlistID2);
      for (varID = 0; varID < nvars; ++varID)
        cdiDefAttFlt(vlistID2, varID, "valid_range", CDI_DATATYPE_FLT64, 2, range);
    }
  */
  const auto streamID2 = cdo_open_write(1);
  cdo_def_vlist(streamID2, vlistID2);

  int tsID = 0;
  while (true)
    {
      const auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
      if (nrecs == 0) break;

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_def_timestep(streamID2, tsID);

      for (int recID = 0; recID < nrecs; ++recID)
        {
          int varID, levelID;
          cdo_inq_record(streamID1, &varID, &levelID);
          const auto &var = varList[varID];
          field.init(var);
          cdo_read_record(streamID1, field);

          // clang-format off
          if      (operatorID == SETMISSVAL) set_missval(field, new_missval);
          else if (operatorID == SETCTOMISS) set_const_to_miss(field, rconst);
          else if (operatorID == SETMISSTOC) set_miss_to_const(field, rconst);
          else if (operatorID == SETRTOMISS) set_range_to_miss(field, rmin, rmax);
          else if (operatorID == SETVRANGE)  set_valid_range(field, rmin, rmax);
          // clang-format on

          cdo_def_record(streamID2, varID, levelID);
          cdo_write_record(streamID2, field);
        }

      tsID++;
    }

  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  cdo_finish();

  return nullptr;
}
