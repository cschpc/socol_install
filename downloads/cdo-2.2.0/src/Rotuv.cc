/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Rotuv      rotuvb          Backward rotation
*/

#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "param_conversion.h"
#include <mpim_grid.h>

static void
rot_uv_back(const int gridID, Varray<double> &us, Varray<double> &vs)
{
  double xpole = 0, ypole = 0, angle = 0;
  if (gridInqType(gridID) == GRID_PROJECTION && gridInqProjType(gridID) == CDI_PROJ_RLL)
    gridInqParamRLL(gridID, &xpole, &ypole, &angle);

  auto nlon = gridInqXsize(gridID);
  auto nlat = gridInqYsize(gridID);

  Varray<double> xvals(nlon), yvals(nlat);
  gridInqXvals(gridID, xvals.data());
  gridInqYvals(gridID, yvals.data());

  // Convert lat/lon units if required
  cdo_grid_to_degree(gridID, CDI_XAXIS, 1, &angle, "angle");
  cdo_grid_to_degree(gridID, CDI_XAXIS, 1, &xpole, "xpole");
  cdo_grid_to_degree(gridID, CDI_XAXIS, nlon, xvals.data(), "grid center lon");
  cdo_grid_to_degree(gridID, CDI_YAXIS, 1, &ypole, "ypole");
  cdo_grid_to_degree(gridID, CDI_YAXIS, nlat, yvals.data(), "grid center lat");

  if (xpole > 180) xpole -= 360;
  if (angle > 180) angle -= 360;

  for (size_t ilat = 0; ilat < nlat; ilat++)
    for (size_t ilon = 0; ilon < nlon; ilon++)
      {
        auto i = ilat * nlon + ilon;
        auto xval = lamrot_to_lam(yvals[ilat], xvals[ilon], ypole, xpole, angle);
        auto yval = phirot_to_phi(yvals[ilat], xvals[ilon], ypole, angle);
        usvs_to_uv(us[i], vs[i], yval, xval, ypole, xpole, &us[i], &vs[i]);
      }
}

#define MAXARG 16384

void *
Rotuv(void *process)
{
  int chcodes[MAXARG];
  const char *chvars[MAXARG];

  cdo_initialize(process);

  operator_input_arg("pairs of u and v in the rotated system");

  const int nch = cdo_operator_argc();
  if (nch % 2) cdo_abort("Odd number of input arguments!");

  bool lvar = false;  // We have a list of codes
  const int len = (int) cdo_operator_argv(0).size();
  const int ix = (cdo_operator_argv(0)[0] == '-') ? 1 : 0;
  for (int i = ix; i < len; ++i)
    if (!isdigit(cdo_operator_argv(0)[i]))
      {
        lvar = true;  // We have a list of variables
        break;
      }

  if (lvar)
    {
      for (int i = 0; i < nch; ++i) chvars[i] = cdo_operator_argv(i).c_str();
    }
  else
    {
      for (int i = 0; i < nch; ++i) chcodes[i] = parameter_to_int(cdo_operator_argv(i));
    }

  auto streamID1 = cdo_open_read(0);

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto vlistID2 = vlistDuplicate(vlistID1);

  VarList varList1;
  varListInit(varList1, vlistID1);

  auto nvars = vlistNvars(vlistID1);

  auto maxrecs = vlistNrecs(vlistID1);
  std::vector<RecordInfo> recList(maxrecs);

  std::vector<std::vector<size_t>> varnmiss(nvars);
  Varray3D<double> vardata(nvars);

  bool lfound[MAXARG];
  for (int i = 0; i < nch; ++i) lfound[i] = false;

  if (lvar)
    {
      for (int varID = 0; varID < nvars; ++varID)
        {
          for (int i = 0; i < nch; ++i)
            if (varList1[varID].name == chvars[i]) lfound[i] = true;
        }
      for (int i = 0; i < nch; ++i)
        if (!lfound[i]) cdo_abort("Variable %s not found!", chvars[i]);
    }
  else
    {
      for (int varID = 0; varID < nvars; ++varID)
        {
          auto code = varList1[varID].code;
          for (int i = 0; i < nch; ++i)
            if (code == chcodes[i]) lfound[i] = true;
        }
      for (int i = 0; i < nch; ++i)
        if (!lfound[i]) cdo_abort("Code %d not found!", chcodes[i]);
    }

  for (int varID = 0; varID < nvars; ++varID)
    {
      auto gridID = varList1[varID].gridID;
      if (!(gridInqType(gridID) == GRID_PROJECTION && gridInqProjType(gridID) == CDI_PROJ_RLL))
        cdo_abort("Only rotated lon/lat grids supported!");

      auto gridsize = gridInqSize(gridID);
      auto nlevels = varList1[varID].nlevels;
      varnmiss[varID].resize(nlevels);
      vardata[varID].resize(nlevels);
      for (int levelID = 0; levelID < nlevels; ++levelID) vardata[varID][levelID].resize(gridsize);
    }

  auto taxisID1 = vlistInqTaxis(vlistID1);
  auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  auto streamID2 = cdo_open_write(1);
  cdo_def_vlist(streamID2, vlistID2);

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

          recList[recID].set(varID, levelID);

          cdo_read_record(streamID1, vardata[varID][levelID].data(), &varnmiss[varID][levelID]);
          if (varnmiss[varID][levelID]) cdo_abort("Missing values unsupported for this operator!");
        }

      for (int i = 0; i < nch; i += 2)
        {
          int varID;
          for (varID = 0; varID < nvars; ++varID)
            {
              if (lvar)
                {
                  if (varList1[varID].name == chvars[i]) break;
                }
              else
                {
                  if (varList1[varID].code == chcodes[i]) break;
                }
            }

          if (varID == nvars) cdo_abort("u-wind not found!");

          auto usvarID = varID;

          for (varID = 0; varID < nvars; ++varID)
            {
              if (lvar)
                {
                  if (varList1[varID].name == chvars[i + 1]) break;
                }
              else
                {
                  if (varList1[varID].code == chcodes[i + 1]) break;
                }
            }

          if (varID == nvars) cdo_abort("v-wind not found!");

          auto vsvarID = varID;

          if (Options::cdoVerbose)
            {
              if (lvar)
                cdo_print("Using var %s [%s](u) and var %s [%s](v)", varList1[usvarID].name, chvars[i], varList1[vsvarID].name,
                          chvars[i + 1]);
              else
                cdo_print("Using code %d [%d](u) and code %d [%d](v)", varList1[usvarID].code, chcodes[i], varList1[vsvarID].code,
                          chcodes[i + 1]);
            }

          auto gridID = varList1[varID].gridID;
          auto nlevels1 = varList1[usvarID].nlevels;
          auto nlevels2 = varList1[vsvarID].nlevels;
          if (nlevels1 != nlevels2) cdo_abort("u-wind and v-wind have different number of levels!");

          for (int levelID = 0; levelID < nlevels1; ++levelID)
            rot_uv_back(gridID, vardata[usvarID][levelID], vardata[vsvarID][levelID]);
        }

      for (int recID = 0; recID < nrecs; ++recID)
        {
          auto [varID, levelID] = recList[recID].get();
          cdo_def_record(streamID2, varID, levelID);
          cdo_write_record(streamID2, vardata[varID][levelID].data(), varnmiss[varID][levelID]);
        }

      tsID++;
    }

  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  cdo_finish();

  return nullptr;
}
