/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Setbox     setclonlatbox   Set lon/lat box to constant
      Setbox     setcindexbox    Set index box to constant
*/

#include <cdi.h>

#include "process_int.h"
#include "param_conversion.h"
#include <mpim_grid.h>
#include "selboxinfo.h"

static void
setcbox(double constant, double *array, int gridID, const SelboxInfo &selboxInfo)
{
  const auto &lat1 = selboxInfo.lat1;
  const auto &lat2 = selboxInfo.lat2;
  const auto &lon11 = selboxInfo.lon11;
  const auto &lon12 = selboxInfo.lon12;
  const auto &lon21 = selboxInfo.lon21;
  const auto &lon22 = selboxInfo.lon22;
  const long nlon = gridInqXsize(gridID);
  const long nlat = gridInqYsize(gridID);

  for (long ilat = 0; ilat < nlat; ilat++)
    for (long ilon = 0; ilon < nlon; ilon++)
      if ((lat1 <= ilat && ilat <= lat2 && ((lon11 <= ilon && ilon <= lon12) || (lon21 <= ilon && ilon <= lon22))))
        {
          array[nlon * ilat + ilon] = constant;
        }
}

static int
get_gridID(int vlistID1, bool operIndexBox)
{
  std::vector<int> gridsFound;

  const auto ngrids = vlistNgrids(vlistID1);
  for (int index = 0; index < ngrids; ++index)
    {
      const auto gridID1 = vlistGrid(vlistID1, index);
      if (gridInqSize(gridID1) == 1) continue;

      const auto gridtype = gridInqType(gridID1);
      const auto projtype = gridInqProjType(gridID1);

      const auto isReg2dGeoGrid = (gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN || gridtype == GRID_CURVILINEAR);
      const auto projHasGeoCoords = (gridtype == GRID_PROJECTION && projtype == CDI_PROJ_RLL);

      if (isReg2dGeoGrid || projHasGeoCoords || (operIndexBox && (gridtype == GRID_GENERIC || gridtype == GRID_PROJECTION)))
        {
          gridsFound.push_back(gridID1);
        }
      else
        {
          if (gridInqSize(gridID1) > 2) cdo_warning("Unsupported grid type: %s", gridNamePtr(gridtype));
        }
    }

  if (gridsFound.size() == 0) cdo_abort("No processable grid found!");
  if (gridsFound.size() > 1) cdo_abort("Too many different grids!");

  const auto gridID = gridsFound[0];
  return gridID;
}

static std::vector<bool>
get_processVars(int vlistID1, int gridID)
{
  const auto nvars = vlistNvars(vlistID1);

  std::vector<bool> processVars(nvars, false);

  int varID;
  for (varID = 0; varID < nvars; ++varID)
    if (gridID == vlistInqVarGrid(vlistID1, varID)) processVars[varID] = true;

  for (varID = 0; varID < nvars; ++varID)
    if (processVars[varID]) break;

  if (varID >= nvars) cdo_abort("No processable variable found!");

  return processVars;
}

void *
Setbox(void *process)
{
  cdo_initialize(process);

  // clang-format off
  const auto SETCLONLATBOX = cdo_operator_add("setclonlatbox", 0, 0, "constant, western and eastern longitude and southern and northern latitude");
  const auto SETCINDEXBOX = cdo_operator_add("setcindexbox", 0, 0, "constant, index of first and last longitude and index of first and last latitude");
  // clang-format on

  (void) SETCLONLATBOX;

  const auto operatorID = cdo_operator_id();
  const auto operIndexBox = (operatorID == SETCINDEXBOX);

  operator_input_arg(cdo_operator_enter(operatorID));

  const auto constant = parameter_to_double(cdo_operator_argv(0));

  const auto streamID1 = cdo_open_read(0);
  const auto vlistID1 = cdo_stream_inq_vlist(streamID1);

  const auto gridID = get_gridID(vlistID1, operIndexBox);

  const auto processVars = get_processVars(vlistID1, gridID);

  operator_input_arg(cdo_operator_enter(operatorID));

  auto selboxInfo = operIndexBox ? gen_index_selbox(1, gridID) : gen_lonlat_selbox(1, gridID);

  const auto vlistID2 = vlistDuplicate(vlistID1);

  const auto taxisID1 = vlistInqTaxis(vlistID1);
  const auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  const auto streamID2 = cdo_open_write(1);

  cdo_def_vlist(streamID2, vlistID2);

  const auto gridsize = gridInqSize(gridID);
  Varray<double> array(gridsize);

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

          if (processVars[varID])
            {
              size_t nmiss;
              cdo_read_record(streamID1, array.data(), &nmiss);

              setcbox(constant, array.data(), gridID, selboxInfo);

              const auto missval = vlistInqVarMissval(vlistID1, varID);
              nmiss = varray_num_mv(gridsize, array, missval);
              cdo_def_record(streamID2, varID, levelID);
              cdo_write_record(streamID2, array.data(), nmiss);
            }
        }

      tsID++;
    }

  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  cdo_finish();

  return nullptr;
}
