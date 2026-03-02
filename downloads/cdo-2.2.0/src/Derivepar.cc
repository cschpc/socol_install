/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Derivepar     gheight           Geopotential height
      Derivepar     sealevelpressure  Sea level pressure
*/

#include <cdi.h>

#include "process_int.h"
#include "cdo_vlist.h"
#include "vertical_interp.h"
#include "stdnametable.h"
#include "const.h"
#include "cdo_zaxis.h"
#include "cdo_options.h"

void geopot_height_halflevel(double *gheight, const double *const ta_fl, const double *const hus_fl, const double *const p_hl,
                             const long ngp, const long nlev);
void geopot_height_fulllevel(double *gheight, const double *const ta_fl, const double *const hus_fl, const double *const p_hl,
                             const long ngp, const long nlev);

static void
check_range_var2d(int stepNum, const Varray<double> &var2d, double rMin, double rMax, const char *varname)
{
  auto mm = varray_min_max(var2d);
  if (mm.min < rMin || mm.max > rMax)
    cdo_warning("%s out of range (min=%g max=%g) [timestep:%d]!", varname, mm.min, mm.max, stepNum);
}

static void
check_range_var3d(int stepNum, int nlevels, size_t gridsize, const Varray<double> &var3d, double rMin, double rMax,
                  const char *varname)
{
  static auto printWarning = true;
  if (printWarning)
    {
      double minVal = 1.e33, maxVal = -1.e33;
      for (int levelID = 0; levelID < nlevels; ++levelID)
        {
          auto mm = varray_min_max(gridsize, &var3d[gridsize * levelID]);
          minVal = std::min(minVal, mm.min);
          maxVal = std::max(maxVal, mm.max);
        }

      if (minVal < rMin || maxVal > rMax)
        {
          printWarning = false;
          cdo_warning("%s out of range (min=%g max=%g) [timestep:%d]!", varname, minVal, maxVal, stepNum);
        }
    }
}

void *
Derivepar(void *process)
{
  int surfaceID = -1;
  int presID = -1;

  cdo_initialize(process);

  // clang-format off
  auto GHEIGHT          = cdo_operator_add("gheight",            0, 0, nullptr);
  auto GHEIGHTHALF      = cdo_operator_add("gheighthalf",        0, 0, nullptr);
  auto SEALEVELPRESSURE = cdo_operator_add("sealevelpressure",   0, 0, nullptr);
  // clang-format on

  auto operatorID = cdo_operator_id();

  operator_check_argc(0);

  auto streamID1 = cdo_open_read(0);
  auto vlistID1 = cdo_stream_inq_vlist(streamID1);

  auto gridID0 = vlistGrid(vlistID1, 0);
  if (gridInqType(gridID0) == GRID_SPECTRAL) cdo_abort("Spectral data unsupported!");

  auto gridsize = vlist_check_gridsize(vlistID1);

  int zaxisID_ML = -1;
  int numHybridLevels = 0, numFullLevels = 0, numHalfLevels = 0;
  auto vct = vlist_read_vct(vlistID1, zaxisID_ML, numHybridLevels, numFullLevels, numHalfLevels);
  const int vctSize = vct.size();

  if (Options::cdoVerbose)
    for (int i = 0; i < vctSize / 2; ++i) cdo_print("vct: %5d %25.17f %25.17f", i, vct[i], vct[vctSize / 2 + i]);

  if (zaxisID_ML == -1) cdo_abort("No 3D variable with hybrid sigma pressure coordinate found!");

  VarList varList1;
  varListInit(varList1, vlistID1);
  varListSetUniqueMemtype(varList1);

  auto nvars = vlistNvars(vlistID1);

  VarIDs varIDs = search_varIDs(varList1, vlistID1, numFullLevels);

  if (Options::cdoVerbose)
    {
      cdo_print("Found:");
      // clang-format off
      if (-1 != varIDs.humID)     cdo_print("  %s -> %s", var_stdname(specific_humidity), varList1[varIDs.humID].name);
      if (-1 != varIDs.tempID)    cdo_print("  %s -> %s", var_stdname(air_temperature), varList1[varIDs.tempID].name);
      if (-1 != varIDs.psID)      cdo_print("  %s -> %s", var_stdname(surface_air_pressure), varList1[varIDs.psID].name);
      if (-1 != varIDs.lnpsID)    cdo_print("  LOG(%s) -> %s", var_stdname(surface_air_pressure), varList1[varIDs.lnpsID].name);
      if (-1 != varIDs.sgeopotID) cdo_print("  %s -> %s", var_stdname(surface_geopotential), varList1[varIDs.sgeopotID].name);
      if (-1 != varIDs.geopotID)  cdo_print("  %s -> %s", var_stdname(geopotential), varList1[varIDs.geopotID].name);
      if (-1 != varIDs.gheightID) cdo_print("  %s -> %s", var_stdname(geopotential_height), varList1[varIDs.gheightID].name);
      // clang-format on
    }

  if (varIDs.lnpsID != -1 && varIDs.lnpsID2 != -1)
    cdo_abort("Found LOG(%s) twice: lsp and lnps!", var_stdname(surface_air_pressure));

  if (varIDs.tempID == -1) cdo_abort("%s not found!", var_stdname(air_temperature));

  for (int varID = 0; varID < nvars; ++varID)
    {
      auto gridID = varList1[varID].gridID;
      auto zaxisID = varList1[varID].zaxisID;

      if (operatorID == SEALEVELPRESSURE) varIDs.humID = -1;

      if (gridInqType(gridID) == GRID_SPECTRAL && zaxisInqType(zaxisID) == ZAXIS_HYBRID)
        cdo_abort("Spectral data on model level unsupported!");

      if (gridInqType(gridID) == GRID_SPECTRAL) cdo_abort("Spectral data unsupported!");
    }

  Varray<double> array(gridsize);
  Varray<double> sgeopot(gridsize), ps(gridsize);
  Varray<double> temp(gridsize * numFullLevels);
  Varray<double> halfPress(gridsize * (numFullLevels + 1));
  Varray<double> hum, gheight;
  if (operatorID == GHEIGHT || operatorID == GHEIGHTHALF)
    {
      if (varIDs.humID == -1)
        cdo_warning("%s not found - using algorithm without %s!", var_stdname(specific_humidity), var_stdname(specific_humidity));
      else
        hum.resize(gridsize * numFullLevels);

      gheight.resize(gridsize * (numFullLevels + 1));
    }

  Varray<double> fullPress, sealevelpressure;
  if (operatorID == SEALEVELPRESSURE)
    {
      fullPress.resize(gridsize * numFullLevels);

      surfaceID = zaxis_from_name("surface");
      sealevelpressure.resize(gridsize);
    }

  if (zaxisID_ML != -1 && varIDs.sgeopotID == -1)
    {
      if (varIDs.geopotID == -1)
        cdo_warning("%s not found - set to zero!", var_stdname(surface_geopotential));
      else
        cdo_print("%s not found - using bottom layer of %s!", var_stdname(surface_geopotential), var_stdname(geopotential));

      varray_fill(sgeopot, 0.0);
    }

  presID = varIDs.lnpsID;
  if (zaxisID_ML != -1 && varIDs.lnpsID == -1)
    {
      if (varIDs.psID == -1)
        cdo_abort("%s not found!", var_stdname(surface_air_pressure));
      else
        presID = varIDs.psID;
    }

  if (Options::cdoVerbose)
    {
      if (presID == varIDs.lnpsID)
        cdo_print("using LOG(%s)", var_stdname(surface_air_pressure));
      else
        cdo_print("using %s", var_stdname(surface_air_pressure));
    }

  auto vlistID2 = vlistCreate();
  vlistDefNtsteps(vlistID2, vlistNtsteps(vlistID1));

  {
    int var_id = -1;
    int varID = -1;

    if (operatorID == GHEIGHT)
      {
        var_id = geopotential_height;
        varID = vlistDefVar(vlistID2, gridID0, zaxisID_ML, TIME_VARYING);
      }
    else if (operatorID == GHEIGHTHALF)
      {
        auto zaxisID_ML_Half = zaxisCreate(ZAXIS_HYBRID_HALF, numHalfLevels);
        zaxisDefVct(zaxisID_ML_Half, 2 * numHalfLevels, vct.data());
        Varray<double> levs(numHalfLevels);
        for (int i = 0; i < numHalfLevels; ++i) levs[i] = i + 1;
        zaxisDefLevels(zaxisID_ML_Half, levs.data());
        var_id = geopotential_height;
        varID = vlistDefVar(vlistID2, gridID0, zaxisID_ML_Half, TIME_VARYING);
      }
    else if (operatorID == SEALEVELPRESSURE)
      {
        var_id = air_pressure_at_sea_level;
        varID = vlistDefVar(vlistID2, gridID0, surfaceID, TIME_VARYING);
      }
    else
      cdo_abort("Internal problem, invalid operatorID: %d!", operatorID);

    vlistDefVarParam(vlistID2, varID, cdiEncodeParam(var_echamcode(var_id), 128, 255));
    cdiDefKeyString(vlistID2, varID, CDI_KEY_NAME, var_name(var_id));
    cdiDefKeyString(vlistID2, varID, CDI_KEY_STDNAME, var_stdname(var_id));
    cdiDefKeyString(vlistID2, varID, CDI_KEY_UNITS, var_units(var_id));
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
          size_t nmiss;
          cdo_read_record(streamID1, array.data(), &nmiss);

          auto offset = gridsize * levelID;

          if (zaxisID_ML != -1)
            {
              if (varID == varIDs.sgeopotID) { varray_copy(gridsize, array, sgeopot); }
              else if (varID == varIDs.geopotID && varIDs.sgeopotID == -1 && (levelID + 1) == numFullLevels)
                {
                  varray_copy(gridsize, array, sgeopot);
                }
              else if (varID == presID)
                {
                  if (varIDs.lnpsID != -1)
                    for (size_t i = 0; i < gridsize; ++i) ps[i] = std::exp(array[i]);
                  else if (varIDs.psID != -1)
                    varray_copy(gridsize, array, ps);
                }
              else if (varID == varIDs.tempID)
                array_copy(gridsize, array.data(), &temp[offset]);
              else if (varID == varIDs.humID)
                array_copy(gridsize, array.data(), &hum[offset]);
            }
        }

      if (zaxisID_ML != -1)
        {
          // check range of psProg
          check_range_var2d(tsID + 1, ps, MIN_PS, MAX_PS, "Surface pressure");

          // check range of surface geopot
          check_range_var2d(tsID + 1, sgeopot, MIN_FIS, MAX_FIS, "Orography");
        }

      check_range_var3d(tsID + 1, varList1[varIDs.tempID].nlevels, gridsize, temp, MIN_T, MAX_T, "Temperature");

      if (varIDs.humID != -1) check_range_var3d(tsID + 1, varList1[varIDs.humID].nlevels, gridsize, hum, -0.1, MAX_Q, "Humidity");

      if (operatorID == GHEIGHT || operatorID == GHEIGHTHALF)
        {
          vct_to_hybrid_pressure((double *) nullptr, halfPress.data(), vct.data(), ps.data(), numFullLevels, gridsize);
          array_copy(gridsize, sgeopot.data(), gheight.data() + gridsize * numFullLevels);
          if (operatorID == GHEIGHT)
            geopot_height_fulllevel(gheight.data(), temp.data(), hum.data(), halfPress.data(), gridsize, numFullLevels);
          else
            geopot_height_halflevel(gheight.data(), temp.data(), hum.data(), halfPress.data(), gridsize, numFullLevels);

          int varID = 0;
          auto nlevels = (operatorID == GHEIGHT) ? numFullLevels : numHalfLevels;
          for (int levelID = 0; levelID < nlevels; ++levelID)
            {
              cdo_def_record(streamID2, varID, levelID);
              cdo_write_record(streamID2, gheight.data() + levelID * gridsize, 0);
            }
        }
      else if (operatorID == GHEIGHTHALF)
        {

          int varID = 0;
          auto nlevels = numHalfLevels;
          for (int levelID = 0; levelID < nlevels; ++levelID)
            {
              cdo_def_record(streamID2, varID, levelID);
              cdo_write_record(streamID2, gheight.data() + levelID * gridsize, 0);
            }
        }
      else if (operatorID == SEALEVELPRESSURE)
        {
          vct_to_hybrid_pressure(fullPress.data(), halfPress.data(), vct.data(), ps.data(), numFullLevels, gridsize);

          extrapolate_P(sealevelpressure.data(), &halfPress[gridsize * (numFullLevels)], &fullPress[gridsize * (numFullLevels - 1)],
                        sgeopot.data(), &temp[gridsize * (numFullLevels - 1)], gridsize);

          cdo_def_record(streamID2, 0, 0);
          cdo_write_record(streamID2, sealevelpressure.data(), 0);
        }
      else
        cdo_abort("Internal error");

      tsID++;
    }

  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  vlistDestroy(vlistID2);

  cdo_finish();

  return nullptr;
}
