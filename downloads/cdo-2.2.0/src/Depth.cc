/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

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
#include "field_functions.h"

template <typename T>
static void
calc_full_depth(size_t gridsize, size_t nlevels, const Varray<T> &thick_c, const Varray<T> &stretch_c, const Varray<T> &zos,
                Varray<T> &fullDepth)
{
  varray_fill(gridsize, fullDepth, 0.0);

  for (size_t k = 1; k < nlevels; ++k)
    {
      auto depth = &fullDepth[k * gridsize];
      auto depthm1 = &fullDepth[(k - 1) * gridsize];
      auto thickm1 = &thick_c[(k - 1) * gridsize];
      for (size_t i = 0; i < gridsize; ++i) depth[i] = depthm1[i] + stretch_c[i] * thickm1[i];
    }

  for (size_t k = 0; k < nlevels; ++k)
    {
      auto depth = &fullDepth[k * gridsize];
      auto thick = &thick_c[k * gridsize];
      for (size_t i = 0; i < gridsize; ++i) depth[i] += 0.5 * stretch_c[i] * thick[i] - zos[i];
    }
}

static void
calc_full_depth(const Field3D &thick_c, const Field3D &stretch_c, const Field3D &zos, Field3D &fullDepth)
{
  auto gridsize = thick_c.gridsize;
  auto nlevels = thick_c.nlevels;
  auto memType = thick_c.memType;
  if (memType == MemType::Float)
    calc_full_depth(gridsize, nlevels, thick_c.vec_f, stretch_c.vec_f, zos.vec_f, fullDepth.vec_f);
  else
    calc_full_depth(gridsize, nlevels, thick_c.vec_d, stretch_c.vec_d, zos.vec_d, fullDepth.vec_d);
}

void *
Depth(void *process)
{
  cdo_initialize(process);

  cdo_operator_add("zsdepth", 0, 0, nullptr);

  auto streamID1 = cdo_open_read(0);

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);

  auto taxisID1 = vlistInqTaxis(vlistID1);
  auto taxisID2 = taxisDuplicate(taxisID1);

  VarList varList1;
  varListInit(varList1, vlistID1);
  varListSetUniqueMemtype(varList1);

  auto nvars = vlistNvars(vlistID1);

  int thickID = -1, zosID = -1, stretchID = -1, draftaveID = -1;
  for (int varID = 0; varID < nvars; ++varID)
    {
      auto varname = string_to_lower(varList1[varID].name);

      // clang-format off
      if      (varname == "prism_thick_c") thickID = varID;
      else if (varname == "stretch_c")     stretchID = varID;
      else if (varname == "zos")           zosID = varID;
      else if (varname == "draftave")      draftaveID = varID;
      // clang-format on
    }

  if (Options::cdoVerbose)
    {
      cdo_print("Found:");
      // clang-format off

      if (-1 != thickID)    cdo_print("  %s -> %s", "prism thickness at cells", varList1[thickID].name);
      if (-1 != stretchID)  cdo_print("  %s -> %s", "zstar surface stretch at cell center", varList1[stretchID].name);
      if (-1 != zosID)      cdo_print("  %s -> %s", "zstar sfc elevation at cell center", varList1[zosID].name);
      if (-1 != draftaveID) cdo_print("  %s -> %s", "draftave", varList1[draftaveID].name);
      // clang-format on
    }

  if (-1 == thickID) cdo_abort("prism_thick_c not found!");
  if (-1 == stretchID) cdo_abort("stretch_c not found!");
  if (-1 == zosID) cdo_abort("zos not found!");
  if (-1 == draftaveID) cdo_warning("draftave not found, set to zero!");

  auto zaxisID = varList1[thickID].zaxisID;
  auto nlevels = varList1[thickID].nlevels;

  auto ngrids = vlistNgrids(vlistID1);
  if (ngrids > 1) cdo_abort("Too many different grids!");

  int index = 0;
  auto gridID = vlistGrid(vlistID1, index);

  auto vlistID2 = vlistCreate();
  vlistDefNtsteps(vlistID2, vlistNtsteps(vlistID1));
  vlistDefTaxis(vlistID2, taxisID2);

  auto depthID = vlistDefVar(vlistID2, gridID, zaxisID, TIME_VARYING);

  cdiDefKeyString(vlistID2, depthID, CDI_KEY_NAME, "depth_c");
  cdiDefKeyString(vlistID2, depthID, CDI_KEY_STDNAME, "depth");
  cdiDefKeyString(vlistID2, depthID, CDI_KEY_LONGNAME, "depth_below_sea");

  Field3DVector vardata1(nvars);

  for (int varID = 0; varID < nvars; ++varID) vardata1[varID].init(varList1[varID]);

  Field3D fullDepth;
  fullDepth.init(varList1[thickID]);

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
          size_t nmiss;
          cdo_read_record(streamID1, vardata1[varID], levelID, &nmiss);
          if (nmiss) cdo_abort("Missing values unsupported!");
        }

      if (-1 != draftaveID) field2_add(vardata1[zosID], vardata1[draftaveID]);
      calc_full_depth(vardata1[thickID], vardata1[stretchID], vardata1[zosID], fullDepth);

      for (int levelID = 0; levelID < nlevels; ++levelID)
        {
          cdo_def_record(streamID2, depthID, levelID);
          cdo_write_record(streamID2, fullDepth, levelID, 0);
        }

      tsID++;
    }

  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  cdo_finish();

  return nullptr;
}
