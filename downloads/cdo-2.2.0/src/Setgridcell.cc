/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Setgridcell    setgridcell      Set grid cells to value
*/

#include <limits>
#include <cdi.h>

#include "process_int.h"
#include "param_conversion.h"
#include "pmlist.h"

template <typename T>
static size_t
set_value(size_t gridsize, Varray<T> &array, T missval, T value)
{
  for (size_t i = 0; i < gridsize; ++i) { array[i] = value; }

  return varray_num_mv(gridsize, array, missval);
}

static void
set_value(Field &field, double value)
{
  if (field.memType == MemType::Float)
    field.nmiss = set_value(field.size, field.vec_f, (float) field.missval, (float) value);
  else
    field.nmiss = set_value(field.size, field.vec_d, field.missval, value);
}

template <typename T>
static size_t
set_value(size_t gridsize, Varray<T> &array, T missval, T value, const Varray<size_t> &cells)
{
  const auto numCells = cells.size();
  for (size_t i = 0; i < numCells; ++i)
    {
      const auto index = cells[i];
      array[index - 1] = value;
    }

  return varray_num_mv(gridsize, array, missval);
}

static void
set_value(Field &field, double value, const Varray<size_t> &cells)
{
  if (field.memType == MemType::Float)
    field.nmiss = set_value(field.size, field.vec_f, (float) field.missval, (float) value, cells);
  else
    field.nmiss = set_value(field.size, field.vec_d, field.missval, value, cells);
}

static void
read_index_from_maskfile(const std::string &maskfile, Varray<size_t> &cells)
{
  size_t cdo_read_mask(const char *maskfile, std::vector<bool> &imask);
  std::vector<bool> mask;
  const auto n = cdo_read_mask(maskfile.c_str(), mask);
  size_t nind = 0;
  for (size_t i = 0; i < n; ++i)
    if (mask[i]) nind++;
  if (nind == 0) cdo_abort("Mask is empty!");

  cells.resize(nind);
  nind = 0;
  for (size_t i = 0; i < n; ++i)
    if (mask[i]) cells[nind++] = i + 1;

  if (nind == 0) cdo_abort("Mask file %s generates no input!", cdo_operator_argv(0));
}

static void
setgridcell_get_parameter(double &constant, Varray<size_t> &cells, std::vector<std::string> &varnames, std::string &maskfile)
{
  const auto pargc = cdo_operator_argc();
  if (pargc)
    {
      const auto &pargv = cdo_get_oper_argv();

      KVList kvlist;
      kvlist.name = cdo_module_name();
      if (kvlist.parse_arguments(pargc, pargv) != 0) cdo_abort("Parse error!");
      if (Options::cdoVerbose) kvlist.print();

      for (const auto &kv : kvlist)
        {
          const auto &key = kv.key;
          if (kv.nvalues < 1) cdo_abort("Missing value for parameter key >%s<!", key);
          const auto &values = kv.values;
          const auto &value = kv.values[0];
          int nvalues = kv.nvalues;

          // clang-format off
          if      (key == "value")
            {
              if (nvalues > 1) cdo_abort("Too many values for parameter key >%s<!", key);
              constant = parameter_to_double(value);
            }
          else if (key == "cell")
            {
              cells.resize(nvalues);
              for (int i = 0; i < nvalues; ++i) cells[i] = parameter_to_size_t(values[i]);
            }
          else if (key == "name")
            {
              varnames.resize(nvalues);
              for (int i = 0; i < nvalues; ++i) varnames[i] = values[i];
            }
          else if (key == "mask")
            {
              if (nvalues > 1) cdo_abort("Too many values for parameter key >%s<!", key);
              maskfile = value;
            }
          else cdo_abort("Invalid parameter key >%s<!", key);
          // clang-format on
        }
    }
}

void *
Setgridcell(void *process)
{
  cdo_initialize(process);

  cdo_operator_add("setgridcell", 0, 0, "value=constant[, cell=grid cell indices (1-N)]");

  operator_input_arg(cdo_operator_enter(0));

  const auto nparam = cdo_operator_argc();
  if (nparam == 0) cdo_abort("Parameter missing!");

  std::string maskfile;
  std::vector<std::string> varnames;
  Varray<size_t> cells;
  double value = DBL_MAX;
  setgridcell_get_parameter(value, cells, varnames, maskfile);
  if (DBL_IS_EQUAL(value, DBL_MAX)) cdo_abort("Parameter <values> not set!");

  if (cells.size() && maskfile.size()) cdo_abort("Too many parameter, choose either cell or mask!");

  if (maskfile.size()) read_index_from_maskfile(maskfile, cells);

  const auto numCells = cells.size();
  size_t minIndex = std::numeric_limits<size_t>::max();
  ;
  size_t maxIndex = 0;
  for (size_t i = 0; i < numCells; ++i)
    {
      minIndex = std::min(minIndex, cells[i]);
      maxIndex = std::max(maxIndex, cells[i]);
    }

  if (numCells && minIndex == 0) cdo_abort("Min cell index is 0, muss be >= 1!");

  const auto streamID1 = cdo_open_read(0);

  const auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  const auto vlistID2 = vlistDuplicate(vlistID1);

  const auto taxisID1 = vlistInqTaxis(vlistID1);
  const auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  VarList varList;
  varListInit(varList, vlistID1);

  const int numVars = varList.size();
  Varray<bool> selectVars(numVars, !varnames.size());
  if (varnames.size())
    {
      for (const auto &varname : varnames)
        {
          auto varFound = false;
          for (int varID = 0; varID < numVars; ++varID)
            {
              if (varname == varList[varID].name)
                {
                  selectVars[varID] = true;
                  varFound = true;
                  break;
                }
            }
          if (!varFound) cdo_abort("Variable %s not found!", varname);
        }
    }

  Field field;

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
          field.init(varList[varID]);
          cdo_read_record(streamID1, field);

          if (selectVars[varID])
            {
              if (numCells)
                {
                  if (maxIndex > field.size) cdo_abort("Max cell index (%zu) > gridsize (%zu)!", maxIndex, field.size);
                  set_value(field, value, cells);
                }
              else { set_value(field, value); }
            }

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
