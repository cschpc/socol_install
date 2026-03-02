/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Merstat    merrange        Meridional range
      Merstat    mermin          Meridional minimum
      Merstat    mermax          Meridional maximum
      Merstat    mersum          Meridional sum
      Merstat    mermean         Meridional mean
      Merstat    meravg          Meridional average
      Merstat    merstd          Meridional standard deviation
      Merstat    merstd1         Meridional standard deviation [Normalize by (n-1)]
      Merstat    mervar          Meridional variance
      Merstat    mervar1         Meridional variance [Normalize by (n-1)]
      Merstat    merpctl         Meridional percentiles
*/

#include <cdi.h>

#include "process_int.h"
#include "param_conversion.h"
#include <mpim_grid.h>
#include "field_functions.h"

static void
add_operators(void)
{
  // clang-format off
  cdo_operator_add("merrange",  FieldFunc_Range,  0, nullptr);
  cdo_operator_add("mermin",    FieldFunc_Min,    0, nullptr);
  cdo_operator_add("mermax",    FieldFunc_Max,    0, nullptr);
  cdo_operator_add("mersum",    FieldFunc_Sum,    0, nullptr);
  cdo_operator_add("mermean",   FieldFunc_Meanw,  1, nullptr);
  cdo_operator_add("meravg",    FieldFunc_Avgw,   1, nullptr);
  cdo_operator_add("mervar",    FieldFunc_Varw,   1, nullptr);
  cdo_operator_add("mervar1",   FieldFunc_Var1w,  1, nullptr);
  cdo_operator_add("merstd",    FieldFunc_Stdw,   1, nullptr);
  cdo_operator_add("merstd1",   FieldFunc_Std1w,  1, nullptr);
  cdo_operator_add("merskew",   FieldFunc_Skew,   0, nullptr);
  cdo_operator_add("merkurt",   FieldFunc_Kurt,   0, nullptr);
  cdo_operator_add("mermedian", FieldFunc_Median, 0, nullptr);
  cdo_operator_add("merpctl",   FieldFunc_Pctl,   0, nullptr);
  // clang-format on
}

void *
Merstat(void *process)
{
  int gridID1, gridID2 = -1, lastgrid = -1;
  int index;

  cdo_initialize(process);

  add_operators();

  auto operatorID = cdo_operator_id();
  auto operfunc = cdo_operator_f1(operatorID);
  auto needWeights = (cdo_operator_f2(operatorID) != 0);

  double pn = 0.0;
  if (operfunc == FieldFunc_Pctl)
    {
      operator_input_arg("percentile number");
      pn = parameter_to_double(cdo_operator_argv(0));
    }
  else { operator_check_argc(0); }

  auto streamID1 = cdo_open_read(0);

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto vlistID2 = vlistDuplicate(vlistID1);

  auto taxisID1 = vlistInqTaxis(vlistID1);
  auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  auto ngrids = vlistNgrids(vlistID1);
  int ndiffgrids = 0;
  for (index = 1; index < ngrids; ++index)
    if (vlistGrid(vlistID1, 0) != vlistGrid(vlistID1, index)) ndiffgrids++;

  if (ndiffgrids > 0) cdo_abort("Too many different grids!");

  index = 0;
  gridID1 = vlistGrid(vlistID1, index);

  if (gridInqType(gridID1) == GRID_LONLAT || gridInqType(gridID1) == GRID_GAUSSIAN || gridInqType(gridID1) == GRID_GENERIC)
    {
      gridID2 = gridToMeridional(gridID1);
    }
  else { cdo_abort("Unsupported gridtype: %s", gridNamePtr(gridInqType(gridID1))); }

  vlistChangeGridIndex(vlistID2, index, gridID2);

  auto streamID2 = cdo_open_write(1);
  cdo_def_vlist(streamID2, vlistID2);

  gridID1 = vlistInqVarGrid(vlistID1, 0);
  const int nlonmax = gridInqXsize(gridID1);  // max nlon?
  auto gridsizemax = vlistGridsizeMax(vlistID1);

  Field field1, field2;
  if (needWeights) field1.weightv.resize(gridsizemax);

  field2.resize(nlonmax);
  field2.grid = gridID2;
  field2.memType = MemType::Double;

  VarList varList1;
  varListInit(varList1, vlistID1);

  int tsID = 0;
  while (true)
    {
      auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
      if (nrecs == 0) break;

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_def_timestep(streamID2, tsID);

      for (int recID = 0; recID < nrecs; ++recID)
        {
          int varID, levelID;
          cdo_inq_record(streamID1, &varID, &levelID);
          field1.init(varList1[varID]);
          cdo_read_record(streamID1, field1);

          field2.missval = field1.missval;

          bool wstatus = false;
          if (needWeights && field1.grid != lastgrid)
            {
              lastgrid = field1.grid;
              wstatus = gridcell_weights(field1.grid, field1.weightv);
            }

          if (wstatus != 0 && tsID == 0 && levelID == 0)
            cdo_warning("Grid cell bounds not available, using constant grid cell area weights for variable %s!",
                        varList1[varID].name);

          (operfunc == FieldFunc_Pctl) ? meridional_pctl(field1, field2, pn) : meridional_function(field1, field2, operfunc);

          cdo_def_record(streamID2, varID, levelID);
          cdo_write_record(streamID2, field2);
        }

      tsID++;
    }

  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  cdo_finish();

  return nullptr;
}
