/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Vertint    gh2hl           Model geometric height level to height level interpolation
*/

#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "cdo_vlist.h"
#include "field_vinterp.h"
#include "stdnametable.h"
#include "util_string.h"
#include "const.h"
#include "cdo_zaxis.h"
#include "param_conversion.h"
#include "vertint_util.h"

static bool
is_height_axis(int zaxisID)
{
  auto isHeight = false;
  if (zaxisInqType(zaxisID) == ZAXIS_REFERENCE)
    {
      auto units = cdo::inq_key_string(zaxisID, CDI_GLOBAL, CDI_KEY_UNITS);
      auto stdname = cdo::inq_key_string(zaxisID, CDI_GLOBAL, CDI_KEY_STDNAME);
      if (stdname == "height" && units.empty()) isHeight = true;
    }
  return isHeight;
}

static int
create_zaxis_height(Varray<double> &heightLevels)
{
  int zaxisID = CDI_UNDEFID;
  const auto &arg1 = cdo_operator_argv(0);
  if (cdo_operator_argc() == 1 && !isdigit(arg1[0]))
    {
      auto zfilename = arg1.c_str();
      auto zfp = std::fopen(zfilename, "r");
      if (zfp)
        {
          zaxisID = zaxis_from_file(zfp, zfilename);
          std::fclose(zfp);
          if (zaxisID == CDI_UNDEFID) cdo_abort("Invalid zaxis description file %s!", zfilename);
          auto nlevels = zaxisInqSize(zaxisID);
          heightLevels.resize(nlevels);
          zaxisInqLevels(zaxisID, heightLevels.data());
        }
      else if (arg1 == "default")
        heightLevels = { 10, 50, 100, 500, 1000, 5000, 10000, 15000, 20000, 25000, 30000 };
      else
        cdo_abort("Open failed on %s", zfilename);
    }
  else { heightLevels = cdo_argv_to_flt(cdo_get_oper_argv()); }

  if (zaxisID == CDI_UNDEFID)
    {
      zaxisID = zaxisCreate(ZAXIS_HEIGHT, heightLevels.size());
      zaxisDefLevels(zaxisID, heightLevels.data());
    }

  return zaxisID;
}

void *
Vertintgh(void *process)
{
  cdo_initialize(process);

  // clang-format off
                cdo_operator_add("gh2hl",   0, 0, "height levels in meter");
  auto GH2HLX = cdo_operator_add("gh2hlx",  0, 0, "height levels in meter");
  // clang-format on

  auto operatorID = cdo_operator_id();

  auto extrapolate = (operatorID == GH2HLX);
  if (extrapolate == false) extrapolate = getenv_extrapolate();

  operator_input_arg(cdo_operator_enter(operatorID));

  Varray<double> heightLevels;
  auto zaxisID2 = create_zaxis_height(heightLevels);

  auto streamID1 = cdo_open_read(0);

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto vlistID2 = vlistDuplicate(vlistID1);

  auto taxisID1 = vlistInqTaxis(vlistID1);
  auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  auto gridsize = vlist_check_gridsize(vlistID1);

  VarList varList1;
  varListInit(varList1, vlistID1);
  varListSetUniqueMemtype(varList1);
  auto memtype = varList1[0].memType;

  auto stdnameHeight_FL = var_stdname(geometric_height_at_full_level_center);
  auto stdnameHeight_HL = var_stdname(geometric_height_at_half_level_center);
  int heightID_FL = -1, heightID_HL = -1;

  auto nvars = vlistNvars(vlistID1);
  for (int varID = 0; varID < nvars; ++varID)
    {
      auto stdname = string_to_lower(cdo::inq_key_string(vlistID1, varID, CDI_KEY_STDNAME));

      // clang-format off
      if (stdname == stdnameHeight_FL) heightID_FL = varID;
      if (stdname == stdnameHeight_HL) heightID_HL = varID;
      // clang-format on
    }

  if (Options::cdoVerbose)
    {
      cdo_print("Found:");
      // clang-format off
      if (-1 != heightID_FL) cdo_print("  %s -> %s", stdnameHeight_FL, varList1[heightID_FL].name);
      if (-1 != heightID_HL) cdo_print("  %s -> %s", stdnameHeight_HL, varList1[heightID_HL].name);
      // clang-format on
    }

  if (-1 == heightID_FL && -1 == heightID_HL) cdo_abort("%s not found!", stdnameHeight_FL);

  auto zaxisID_FL = (-1 == heightID_FL) ? -1 : varList1[heightID_FL].zaxisID;
  auto zaxisID_HL = (-1 == heightID_HL) ? -1 : varList1[heightID_HL].zaxisID;
  auto numFullLevels = (-1 == zaxisID_FL) ? 0 : zaxisInqSize(zaxisID_FL);
  auto numHalfLevels = (-1 == zaxisID_HL) ? 0 : zaxisInqSize(zaxisID_HL);

  auto nzaxis = vlistNzaxis(vlistID1);
  for (int index = 0; index < nzaxis; ++index)
    {
      auto zaxisID = vlistZaxis(vlistID1, index);
      auto nlevels = zaxisInqSize(zaxisID);
      if (zaxisID == zaxisID_FL || zaxisID == zaxisID_HL
          || (is_height_axis(zaxisID) && (nlevels == numHalfLevels || nlevels == numFullLevels)))
        vlistChangeZaxis(vlistID2, zaxisID, zaxisID2);
    }

  VarList varList2;
  varListInit(varList2, vlistID2);
  varListSetMemtype(varList2, memtype);

  Field heightBottom;

  Varray<size_t> pnmiss;
  if (!extrapolate) pnmiss.resize(heightLevels.size());

  std::vector<int> vertIndex_FL, vertIndex_HL;
  if (-1 != heightID_FL) vertIndex_FL.resize(gridsize * heightLevels.size());
  if (-1 != heightID_HL) vertIndex_HL.resize(gridsize * heightLevels.size());

  std::vector<bool> processVars(nvars), interpVars(nvars);
  Varray2D<size_t> varnmiss(nvars);
  Field3DVector vardata1(nvars), vardata2(nvars);

  auto maxlev = std::max(std::max(numFullLevels, numHalfLevels), (int) heightLevels.size());

  for (int varID = 0; varID < nvars; ++varID)
    {
      const auto &var = varList1[varID];
      auto isHeightAxis = is_height_axis(var.zaxisID);

      if (gridInqType(var.gridID) == GRID_SPECTRAL) cdo_abort("Spectral data unsupported!");

      vardata1[varID].init(var);

      interpVars[varID] = (var.zaxisID == zaxisID_FL || var.zaxisID == zaxisID_HL
                          || (isHeightAxis && (var.nlevels == numHalfLevels || var.nlevels == numFullLevels)));

      if (interpVars[varID])
        {
          varnmiss[varID].resize(maxlev, 0);
          vardata2[varID].init(varList2[varID]);
        }
      else
        {
          if (isHeightAxis && var.nlevels > 1)
            {
              if (-1 == heightID_FL && -1 != heightID_HL && var.nlevels == (numHalfLevels - 1))
                cdo_abort("%s not found (needed for %s)!", stdnameHeight_FL, var.name);
              else if (-1 != heightID_FL && -1 == heightID_HL && var.nlevels == (numFullLevels + 1))
                cdo_abort("%s not found (needed for %s)!", stdnameHeight_HL, var.name);
              else
                cdo_warning("Parameter %d has wrong number of levels, skipped! (name=%s nlevel=%d)", varID + 1,
                            var.name, var.nlevels);
            }
          varnmiss[varID].resize(var.nlevels);
        }
    }

  for (int varID = 0; varID < nvars; ++varID)
    {
      const auto &var = varList1[varID];
      if (interpVars[varID] && var.isConstant) vlistDefVarTimetype(vlistID2, varID, TIME_VARYING);
    }

  auto streamID2 = cdo_open_write(1);
  cdo_def_vlist(streamID2, vlistID2);

  int tsID = 0;
  while (true)
    {
      auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
      if (nrecs == 0) break;

      for (int varID = 0; varID < nvars; ++varID)
        {
          processVars[varID] = false;
          const auto &var = varList1[varID];
          for (int levelID = 0; levelID < var.nlevels; ++levelID) varnmiss[varID][levelID] = 0;
        }

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_def_timestep(streamID2, tsID);

      for (int recID = 0; recID < nrecs; ++recID)
        {
          int varID, levelID;
          cdo_inq_record(streamID1, &varID, &levelID);
          cdo_read_record(streamID1, vardata1[varID], levelID, &varnmiss[varID][levelID]);
          processVars[varID] = true;
        }

      for (int varID = 0; varID < nvars; ++varID)
        if (interpVars[varID]) processVars[varID] = true;

      auto lreverse = true;
      if (-1 != heightID_FL && (tsID == 0 || !varList1[heightID_FL].isConstant))
        {
          gen_vert_index(vertIndex_FL, heightLevels, vardata1[heightID_FL], gridsize, lreverse);
          if (!extrapolate)
            {
              heightBottom.init(varList1[heightID_FL]);
              field_copy(vardata1[heightID_FL], numFullLevels - 1, heightBottom);
              gen_vert_index_mv(vertIndex_FL, heightLevels, gridsize, heightBottom, pnmiss, lreverse);
            }
        }

      if (-1 != heightID_HL && (tsID == 0 || !varList1[heightID_HL].isConstant))
        {
          gen_vert_index(vertIndex_HL, heightLevels, vardata1[heightID_HL], gridsize, lreverse);
          if (!extrapolate)
            {
              heightBottom.init(varList1[heightID_HL]);
              field_copy(vardata1[heightID_HL], numHalfLevels - 1, heightBottom);
              gen_vert_index_mv(vertIndex_HL, heightLevels, gridsize, heightBottom, pnmiss, lreverse);
            }
        }

      for (int varID = 0; varID < nvars; ++varID)
        {
          if (processVars[varID])
            {
              const auto &var = varList1[varID];
              if (tsID > 0 && !interpVars[varID] && var.isConstant) continue;

              if (interpVars[varID])
                {
                  if (var.nlevels != numFullLevels && var.nlevels != numHalfLevels)
                    cdo_abort("Number of generalized height level differ from full/half level (param=%s)!", var.name);

                  for (int levelID = 0; levelID < var.nlevels; ++levelID)
                    {
                      if (varnmiss[varID][levelID]) cdo_abort("Missing values unsupported for this operator!");
                    }

                  const auto &height3D = (var.nlevels == numFullLevels) ? vardata1[heightID_FL] : vardata1[heightID_HL];
                  const auto &vertIndex3D = (var.nlevels == numFullLevels) ? vertIndex_FL : vertIndex_HL;
                  vertical_interp_X(height3D, vardata1[varID], vardata2[varID], vertIndex3D, heightLevels, gridsize);

                  if (!extrapolate) varray_copy(heightLevels.size(), pnmiss, varnmiss[varID]);
                }

              for (int levelID = 0; levelID < varList2[varID].nlevels; ++levelID)
                {
                  cdo_def_record(streamID2, varID, levelID);
                  cdo_write_record(streamID2, interpVars[varID] ? vardata2[varID] : vardata1[varID], levelID,
                                   varnmiss[varID][levelID]);
                }
            }
        }

      tsID++;
    }

  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  cdo_finish();

  return nullptr;
}
