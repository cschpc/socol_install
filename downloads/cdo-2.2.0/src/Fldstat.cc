/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Fldstat    fldrange        Field range (max-min)
      Fldstat    fldmin          Field minimum
      Fldstat    fldmax          Field maximum
      Fldstat    fldsum          Field sum
      Fldstat    fldmean         Field mean
      Fldstat    fldavg          Field average
      Fldstat    fldstd          Field standard deviation
      Fldstat    fldstd1         Field standard deviation [Normalize by (n-1)]
      Fldstat    fldvar          Field variance
      Fldstat    fldvar1         Field variance [Normalize by (n-1)]
      Fldstat    fldpctl         Field percentiles
*/

#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "param_conversion.h"
#include <mpim_grid.h>
#include "pmlist.h"
#include "cdo_zaxis.h"
#include "printinfo.h"
#include "progress.h"
#include "field_functions.h"

void gridcell_areas(int gridID, Varray<double> &array);

template <typename T>
static void
print_location_LL(int operfunc, const VarList &varList, int varID, int levelID, int gridID, double sglval,
                  const Varray<T> &fieldvec, CdiDateTime vDateTime)
{
  static auto printHeader = true;
  auto code = varList[varID].code;

  auto isReg2d = (gridInqType(gridID) == GRID_GAUSSIAN || gridInqType(gridID) == GRID_LONLAT);

  if (isReg2d || gridInqType(gridID) == GRID_CURVILINEAR || gridInqType(gridID) == GRID_UNSTRUCTURED)
    {
      auto zaxisID = varList[varID].zaxisID;
      auto level = cdo_zaxis_inq_level(zaxisID, levelID);
      auto gridsize = gridInqSize(gridID);
      auto nlon = gridInqXsize(gridID);
      const T value = sglval;
      for (size_t ij = 0; ij < gridsize; ++ij)
        {
          if (dbl_is_equal(fieldvec[ij], value))
            {
              auto j = ij / nlon;
              auto i = ij - j * nlon;
              auto xval = gridInqXval(gridID, isReg2d ? i : ij);
              auto yval = gridInqYval(gridID, isReg2d ? j : ij);
              if (printHeader)
                {
                  fprintf(stdout, "  Date     Time     Code  Level   Lon      Lat          %s\n",
                          (operfunc == FieldFunc_Min) ? "Minval" : "Maxval");
                  printHeader = false;
                }

              fprintf(stdout, "%s %s %3d %7g %9.7g %9.7g %12.5g\n", date_to_string(vDateTime.date).c_str(),
                      time_to_string(vDateTime.time).c_str(), code, level, xval, yval, sglval);
              break;
            }
        }
    }
}

template <typename T>
static void
field_mul_weights(size_t len, Varray<T> &v1, const Varray<double> &v2, size_t nmiss, T missval)
{
  if (nmiss)
    {
      for (size_t i = 0; i < len; ++i)
        if (!dbl_is_equal(v1[i], missval)) v1[i] *= v2[i];
    }
  else
    {
      for (size_t i = 0; i < len; ++i) v1[i] *= v2[i];
    }
}

static void
field_mul_weights(Field &field)
{
  if (field.memType == MemType::Float)
    field_mul_weights(field.size, field.vec_f, field.weightv, field.nmiss, (float) field.missval);
  else
    field_mul_weights(field.size, field.vec_d, field.weightv, field.nmiss, field.missval);
}

static void
remove_global_grid_attr(int vlistID)
{
  cdiDelAtt(vlistID, CDI_GLOBAL, "ICON_grid_file_uri");
  cdiDelAtt(vlistID, CDI_GLOBAL, "number_of_grid_used");
  cdiDelAtt(vlistID, CDI_GLOBAL, "uuidOfHGrid");
}

static int
gen_target_gridpoint(int gridID1)
{
  int gridID2 = -1;

  auto gridType = gridInqType(gridID1);
  if (gridType == GRID_UNSTRUCTURED)
    {
      gridID2 = gridCreate(gridType, 1);
      grid_copy_names(gridID1, gridID2);
    }
  else if (gridType == GRID_GENERIC)
    {
      gridID2 = gridCreate(GRID_GENERIC, 1);
      grid_copy_names(gridID1, gridID2);
      gridDefXsize(gridID2, 1);
      gridDefYsize(gridID2, 1);
    }
  else
    {
      gridID2 = gridCreate(GRID_LONLAT, 1);
      gridDefXsize(gridID2, 1);
      gridDefYsize(gridID2, 1);
    }

  auto value = 0.0;
  gridDefXvals(gridID2, &value);
  gridDefYvals(gridID2, &value);

  return gridID2;
}

static void
printWeightsWarning(int ngrids, const std::string &varname)
{
  if (ngrids == 1)
    cdo_warning("Grid cell bounds not available, using constant grid cell area weights!");
  else
    cdo_warning("Grid cell bounds not available, using constant grid cell area weights for variable %s!", varname);
}

static void
fldstatGetParameter(bool &useWeights)
{
  auto pargc = cdo_operator_argc();
  if (pargc)
    {
      auto pargv = cdo_get_oper_argv();

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

          if (key == "weights")
            useWeights = parameter_to_bool(value);
          else
            cdo_abort("Invalid parameter key >%s<!", key);
        }
    }
}

void *
Fldstat(void *process)
{
  cdo_initialize(process);

  // clang-format off
                cdo_operator_add("fldrange",  FieldFunc_Range,  0, nullptr);
                cdo_operator_add("fldmin",    FieldFunc_Min,    0, nullptr);
                cdo_operator_add("fldmax",    FieldFunc_Max,    0, nullptr);
                cdo_operator_add("fldsum",    FieldFunc_Sum,    0, nullptr);
  auto FLDINT = cdo_operator_add("fldint",    FieldFunc_Sum,    0, nullptr);
                cdo_operator_add("fldmean",   FieldFunc_Meanw,  1, nullptr);
                cdo_operator_add("fldavg",    FieldFunc_Avgw,   1, nullptr);
                cdo_operator_add("fldstd",    FieldFunc_Stdw,   1, nullptr);
                cdo_operator_add("fldstd1",   FieldFunc_Std1w,  1, nullptr);
                cdo_operator_add("fldvar",    FieldFunc_Varw,   1, nullptr);
                cdo_operator_add("fldvar1",   FieldFunc_Var1w,  1, nullptr);
                cdo_operator_add("fldskew",   FieldFunc_Skew,   0, nullptr);
                cdo_operator_add("fldkurt",   FieldFunc_Kurt,   0, nullptr);
                cdo_operator_add("fldmedian", FieldFunc_Median, 0, nullptr);
                cdo_operator_add("fldcount",  FieldFunc_Count,  0, nullptr);
                cdo_operator_add("fldpctl",   FieldFunc_Pctl,   0, nullptr);
  // clang-format on

  auto operatorID = cdo_operator_id();
  auto operfunc = cdo_operator_f1(operatorID);
  auto needWeights = (cdo_operator_f2(operatorID) != 0);
  auto needCellarea = (operatorID == FLDINT);
  auto useWeights = true;

  double pn = 0.0;
  if (operfunc == FieldFunc_Pctl)
    {
      operator_input_arg("percentile number");
      pn = parameter_to_double(cdo_operator_argv(0));
    }
  else if (needWeights) { fldstatGetParameter(useWeights); }
  else { operator_check_argc(0); }

  auto streamID1 = cdo_open_read(0);

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto vlistID2 = vlistDuplicate(vlistID1);
  remove_global_grid_attr(vlistID2);

  auto taxisID1 = vlistInqTaxis(vlistID1);
  auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  auto ngrids = vlistNgrids(vlistID1);

  for (int index = 0; index < ngrids; ++index)
    {
      auto gridID1 = vlistGrid(vlistID1, index);
      auto gridID2 = gen_target_gridpoint(gridID1);
      vlistChangeGridIndex(vlistID2, index, gridID2);
    }

  auto streamID2 = cdo_open_write(1);

  cdo_def_vlist(streamID2, vlistID2);

  auto gridsizemax = vlistGridsizeMax(vlistID1);

  Field field;
  if (needWeights || needCellarea)
    {
      field.weightv.resize(gridsizemax);
      if (needWeights && !useWeights)
        {
          cdo_print("Using constant grid cell area weights!");
          for (size_t i = 0; i < gridsizemax; ++i) field.weightv[i] = 1.0;
        }
    }

  VarList varList1;
  varListInit(varList1, vlistID1);

  auto ntsteps = vlistNtsteps(vlistID1);

  progress::init();

  int lastgrid = -1;
  int tsID = 0;
  while (true)
    {
      auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
      if (nrecs == 0) break;

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_def_timestep(streamID2, tsID);

      auto vDateTime = taxisInqVdatetime(taxisID1);

      for (int recID = 0; recID < nrecs; ++recID)
        {
          double fstatus = (ntsteps > 1) ? (tsID + (recID + 1.0) / nrecs) / ntsteps : 1.0;
          if (!Options::cdoVerbose) progress::update(0, 1, fstatus);

          int varID, levelID;
          cdo_inq_record(streamID1, &varID, &levelID);
          field.init(varList1[varID]);
          cdo_read_record(streamID1, field);

          if (needWeights && field.grid != lastgrid)
            {
              lastgrid = field.grid;
              field.weightv[0] = 1;
              if (useWeights && field.size > 1)
                {
                  auto wstatus = (gridcell_weights(field.grid, field.weightv) != 0);
                  if (wstatus && tsID == 0 && levelID == 0) printWeightsWarning(ngrids, varList1[varID].name);
                }
            }

          if (needCellarea && field.grid != lastgrid)
            {
              lastgrid = field.grid;
              gridcell_areas(field.grid, field.weightv);
            }

          if (needCellarea) field_mul_weights(field);

          auto singleValue = (operfunc == FieldFunc_Pctl) ? field_pctl(field, pn) : field_function(field, operfunc);

          if (Options::cdoVerbose && (operfunc == FieldFunc_Min || operfunc == FieldFunc_Max))
            {
              if (field.memType == MemType::Float)
                print_location_LL(operfunc, varList1, varID, levelID, field.grid, singleValue, field.vec_f, vDateTime);
              else
                print_location_LL(operfunc, varList1, varID, levelID, field.grid, singleValue, field.vec_d, vDateTime);
            }

          size_t nmiss = dbl_is_equal(singleValue, field.missval);

          cdo_def_record(streamID2, varID, levelID);
          cdo_write_record(streamID2, &singleValue, nmiss);
        }

      tsID++;
    }

  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  vlistDestroy(vlistID2);

  cdo_finish();

  return nullptr;
}
