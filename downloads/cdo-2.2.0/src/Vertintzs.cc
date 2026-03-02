/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Vertint    zs2zl           Model depth level to depth level interpolation
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

static bool
is_depth_axis(int zaxisID)
{
  return (zaxisInqType(zaxisID) == ZAXIS_DEPTH_BELOW_SEA);
}

static int
create_zaxis_depth(Varray<double> &depthLevels)
{
  int zaxisID = CDI_UNDEFID;
  auto &arg1 = cdo_operator_argv(0);
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
          depthLevels.resize(nlevels);
          zaxisInqLevels(zaxisID, depthLevels.data());
        }
      else if (arg1 == "default")
        depthLevels = { 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000 };
      else
        cdo_abort("Open failed on %s", zfilename);
    }
  else { depthLevels = cdo_argv_to_flt(cdo_get_oper_argv()); }

  if (zaxisID == CDI_UNDEFID)
    {
      zaxisID = zaxisCreate(ZAXIS_DEPTH_BELOW_SEA, depthLevels.size());
      zaxisDefLevels(zaxisID, depthLevels.data());
    }

  return zaxisID;
}

class ModuleVertintzs
{
private:
  int nvars;
  std::vector<bool> processVars;
  std::vector<bool> interpVars;

  CdoStreamID streamID1;
  CdoStreamID streamID2;

  int taxisID1;
  int taxisID2;

  VarList varList1;
  VarList varList2;
  Field3DVector vardata1;
  Field3DVector vardata2;

  Varray2D<size_t> varnmiss;

  int depthID = -1;

  Field3D fullDepth;

  std::vector<int> vertIndexFull;
  Varray<double> depthLevels;

  int gridsize;
  bool extrapolate = true;  // do not use missing values
                            //
  Field depthBottom;
  int numFullLevels;

  Varray<size_t> pnmiss;

public:
  void
  init(void *process)
  {

    cdo_initialize(process);

    cdo_operator_add("zs2zl", 0, 0, "depth levels in meter");

    auto operatorID = cdo_operator_id();

    operator_input_arg(cdo_operator_enter(operatorID));

    auto zaxisID2 = create_zaxis_depth(depthLevels);

    streamID1 = cdo_open_read(0);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = vlistDuplicate(vlistID1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    gridsize = vlist_check_gridsize(vlistID1);

    varListInit(varList1, vlistID1);
    varListSetUniqueMemtype(varList1);
    auto memType = varList1[0].memType;

    nvars = vlistNvars(vlistID1);

    for (int varID = 0; varID < nvars; ++varID)
      {
        if (string_to_lower(varList1[varID].name) == "depth_c") depthID = varID;
      }

    if (Options::cdoVerbose)
      {
        cdo_print("Found:");
        // clang-format off
      if (-1 != depthID)   cdo_print("  %s -> %s", "zstar depth at cell center", varList1[depthID].name);
        // clang-format on
      }

    if (-1 == depthID) cdo_abort("depth_c not found!");

    auto zaxisIDfull = (-1 == depthID) ? -1 : varList1[depthID].zaxisID;
    numFullLevels = (-1 == zaxisIDfull) ? 0 : zaxisInqSize(zaxisIDfull);

    auto nzaxis = vlistNzaxis(vlistID1);
    for (int index = 0; index < nzaxis; ++index)
      {
        auto zaxisID = vlistZaxis(vlistID1, index);
        auto nlevels = zaxisInqSize(zaxisID);
        if (zaxisID == zaxisIDfull || (is_depth_axis(zaxisID) && nlevels == numFullLevels))
          vlistChangeZaxis(vlistID2, zaxisID, zaxisID2);
      }

    varListInit(varList2, vlistID2);
    varListSetMemtype(varList2, memType);

    if (!extrapolate) pnmiss.resize(depthLevels.size());

    vertIndexFull.resize(gridsize * depthLevels.size());

    processVars = std::vector<bool>(nvars);
    interpVars = std::vector<bool>(nvars);
    varnmiss = Varray2D<size_t>(nvars);
    vardata1 = Field3DVector(nvars);
    vardata2 = Field3DVector(nvars);

    auto maxlev = std::max(numFullLevels, (int) depthLevels.size());

    for (int varID = 0; varID < nvars; ++varID)
      {
        const auto &var = varList1[varID];

        if (gridInqType(var.gridID) == GRID_SPECTRAL) cdo_abort("Spectral data unsupported!");

        vardata1[varID].init(var);
        varnmiss[varID].resize(maxlev, 0);

        interpVars[varID] = (var.zaxisID == zaxisIDfull || (is_depth_axis(var.zaxisID) && (var.nlevels == numFullLevels)));

        if (interpVars[varID]) { vardata2[varID].init(varList2[varID]); }
        else if (is_depth_axis(var.zaxisID) && var.nlevels > 1)
          {

            cdo_warning("Parameter %d has wrong number of levels, skipped! (name=%s nlevel=%d)", varID + 1, var.name,
                        var.nlevels);
          }
      }

    fullDepth.init(varList1[depthID]);

    for (int varID = 0; varID < nvars; ++varID)
      {
        if (interpVars[varID] && varList1[varID].isConstant) vlistDefVarTimetype(vlistID2, varID, TIME_VARYING);
      }

    streamID2 = cdo_open_write(1);

    cdo_def_vlist(streamID2, vlistID2);
  }

  void
  run()
  {
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

        if (tsID == 0 || !varList1[depthID].isConstant)
          {
            constexpr auto lreverse = false;
            field_copy(vardata1[depthID], fullDepth);
            gen_vert_index(vertIndexFull, depthLevels, fullDepth, gridsize, lreverse);
            if (!extrapolate)
              {
                depthBottom.init(varList1[depthID]);
                field_copy(fullDepth, numFullLevels - 1, depthBottom);
                gen_vert_index_mv(vertIndexFull, depthLevels, gridsize, depthBottom, pnmiss, lreverse);
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
                    if (var.nlevels != numFullLevels)
                      cdo_abort("Number of depth level differ from full level (param=%s)!", var.name);

                    for (int levelID = 0; levelID < var.nlevels; ++levelID)
                      {
                        if (varnmiss[varID][levelID]) cdo_abort("Missing values unsupported for this operator!");
                      }

                    vertical_interp_X(fullDepth, vardata1[varID], vardata2[varID], vertIndexFull, depthLevels, gridsize);

                    if (!extrapolate) varray_copy(depthLevels.size(), pnmiss, varnmiss[varID]);
                  }

                for (int levelID = 0; levelID < varList2[varID].nlevels; ++levelID)
                  {
                    cdo_def_record(streamID2, varID, levelID);
                    auto varout = (interpVars[varID] ? vardata2[varID] : vardata1[varID]);
                    cdo_write_record(streamID2, varout, levelID, varnmiss[varID][levelID]);
                  }
              }
          }

        tsID++;
      }
  }

  void
  close()
  {
    cdo_stream_close(streamID2);
    cdo_stream_close(streamID1);

    cdo_finish();
  }
};

void *
Vertintzs(void *process)
{
  ModuleVertintzs vertintzs;
  vertintzs.init(process);
  vertintzs.run();
  vertintzs.close();
  return nullptr;
}
