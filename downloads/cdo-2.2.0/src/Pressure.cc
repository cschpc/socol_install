/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Pressure    pressure_fl          Pressure on full hybrid levels
      Pressure    pressure_hl          Pressure on half hybrid levels
      Pressure    deltap               Difference of two half hybrid levels
*/

#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "cdo_vlist.h"
#include "vertical_interp.h"
#include "stdnametable.h"
#include "const.h"

void *
Pressure(void *process)
{
  cdo_initialize(process);

  // clang-format off
  auto PRESSURE_FL = cdo_operator_add("pressure_fl", 0, 0, nullptr);
  auto PRESSURE_HL = cdo_operator_add("pressure_hl", 0, 0, nullptr);
  auto DELTAP      = cdo_operator_add("deltap",      0, 0, nullptr);
  // clang-format on

  auto operatorID = cdo_operator_id();

  operator_check_argc(0);

  auto streamID1 = cdo_open_read(0);
  auto vlistID1 = cdo_stream_inq_vlist(streamID1);

  VarList varList1;
  varListInit(varList1, vlistID1);
  // varListSetUniqueMemtype(varList1);
  // auto memType = varList1[0].memType;

  auto gridsize = vlist_check_gridsize(vlistID1);

  int zaxisID_ML = -1;
  int numHybridLevels = 0, numFullLevels = 0, numHalfLevels = 0;
  auto vct = vlist_read_vct(vlistID1, zaxisID_ML, numHybridLevels, numFullLevels, numHalfLevels);

  auto hasVars3D = (zaxisID_ML != -1 && gridsize > 0);
  if (!hasVars3D) cdo_abort("No 3D variable with hybrid sigma pressure coordinate found!");

  Varray<double> psProg(gridsize);
  Varray<double> deltap(gridsize * numFullLevels), fullPress(gridsize * numFullLevels), halfPress(gridsize * numHalfLevels);

  int zaxisID_PL = -1;
  int zaxisType = -1;
  if (operatorID == PRESSURE_FL || operatorID == DELTAP)
    {
      if (numFullLevels == numHybridLevels)
        zaxisID_PL = zaxisID_ML;
      else
        {
          if (Options::cdoVerbose) cdo_print("Creating ZAXIS_HYBRID .. (numFullLevels=%d)", numFullLevels);
          zaxisType = ZAXIS_HYBRID;
          numHybridLevels = numFullLevels;
        }
    }
  else
    {
      if (numHalfLevels == numHybridLevels)
        zaxisID_PL = zaxisID_ML;
      else
        {
          if (Options::cdoVerbose) cdo_print("Creating ZAXIS_HYBRID_HALF .. (numHalfLevels=%d)", numHalfLevels);
          zaxisType = ZAXIS_HYBRID_HALF;
          numHybridLevels = numHalfLevels;
        }
    }

  if (zaxisID_PL == -1)
    {
      zaxisID_PL = zaxisCreate(zaxisType, numHybridLevels);

      Varray<double> level(numHalfLevels);
      for (int l = 0; l < numHalfLevels; ++l) level[l] = l + 1;
      zaxisDefLevels(zaxisID_PL, level.data());
      zaxisDefVct(zaxisID_PL, 2 * numHalfLevels, vct.data());
    }

  VarIDs varIDs = search_varIDs(varList1, vlistID1, numFullLevels);

  if (Options::cdoVerbose)
    {
      cdo_print("Found:");
      // clang-format off
      if (-1 != varIDs.psID)      cdo_print("  %s -> %s", var_stdname(surface_air_pressure), varList1[varIDs.psID].name);
      if (-1 != varIDs.lnpsID)    cdo_print("  LOG(%s) -> %s", var_stdname(surface_air_pressure), varList1[varIDs.lnpsID].name);
      // clang-format on
    }

  auto pvarID = varIDs.lnpsID;
  if (zaxisID_ML != -1 && varIDs.lnpsID != -1)
    {
      auto gridID = varList1[varIDs.lnpsID].gridID;
      if (gridInqType(gridID) == GRID_SPECTRAL)
        {
          varIDs.lnpsID = -1;
          cdo_warning("Spectral LOG(%s) not supported - using %s!", var_stdname(surface_air_pressure),
                      var_stdname(surface_air_pressure));
        }
    }

  if (zaxisID_ML != -1 && varIDs.lnpsID == -1)
    {
      pvarID = varIDs.psID;
      if (varIDs.psID == -1) cdo_abort("%s not found!", var_stdname(surface_air_pressure));
    }

  auto gridID = varList1[pvarID].gridID;
  if (gridInqType(gridID) == GRID_SPECTRAL)
    cdo_abort("%s on spectral representation not supported!", var_stdname(surface_air_pressure));

  Varray<double> array(gridsize);

  auto vlistID2 = vlistCreate();
  vlistDefNtsteps(vlistID2, vlistNtsteps(vlistID1));

  auto ovarID = vlistDefVar(vlistID2, gridID, zaxisID_PL, TIME_VARYING);
  vlistDefVarParam(vlistID2, ovarID, cdiEncodeParam(1, 255, 255));
  cdiDefKeyString(vlistID2, ovarID, CDI_KEY_NAME, "pressure");
  cdiDefKeyString(vlistID2, ovarID, CDI_KEY_STDNAME, "air_pressure");
  cdiDefKeyString(vlistID2, ovarID, CDI_KEY_UNITS, "Pa");

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

          if (varID == pvarID)
            {
              size_t nmiss;
              cdo_read_record(streamID1, array.data(), &nmiss);
              if (nmiss) cdo_abort("Missing valus unsupported!");
            }
        }

      if (zaxisID_ML != -1)
        {
          if (varIDs.lnpsID != -1)
            for (size_t i = 0; i < gridsize; ++i) psProg[i] = std::exp(array[i]);
          else if (varIDs.psID != -1)
            varray_copy(gridsize, array, psProg);

          // check range of psProg
          auto mm = varray_min_max(psProg);
          if (mm.min < MIN_PS || mm.max > MAX_PS) cdo_warning("Surface pressure out of range (min=%g max=%g)!", mm.min, mm.max);

          vct_to_hybrid_pressure(fullPress.data(), halfPress.data(), vct.data(), psProg.data(), numFullLevels, gridsize);
        }

      double *pout = nullptr;
      int nlevels = 0;
      if (operatorID == PRESSURE_FL)
        {
          nlevels = numFullLevels;
          pout = fullPress.data();
        }
      else if (operatorID == DELTAP)
        {
          nlevels = numFullLevels;
          for (int k = 0; k < numFullLevels; ++k)
            for (size_t i = 0; i < gridsize; ++i)
              deltap[k * gridsize + i] = halfPress[(k + 1) * gridsize + i] - halfPress[k * gridsize + i];

          pout = deltap.data();
        }
      else if (operatorID == PRESSURE_HL)
        {
          nlevels = numHalfLevels;
          pout = halfPress.data();
        }

      int varID = 0;
      for (int levelID = 0; levelID < nlevels; ++levelID)
        {
          cdo_def_record(streamID2, varID, levelID);
          cdo_write_record(streamID2, &pout[levelID * gridsize], 0);
        }

      tsID++;
    }

  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  cdo_finish();

  return nullptr;
}
