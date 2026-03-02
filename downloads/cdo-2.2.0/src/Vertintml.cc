/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Vertint    ml2pl           Model to pressure level interpolation
      Vertint    ml2hl           Model to height level interpolation
*/

#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "cdo_vlist.h"
#include "field_vinterp.h"
#include "stdnametable.h"
#include "constants.h"
#include "const.h"
#include "cdo_zaxis.h"
#include "param_conversion.h"
#include "vertint_util.h"

static void
field_copy_div_2d_to_3d(MemType memType, size_t gridsize, int nlevels, const Field &field2d, Field3D &field3d)
{
  if (memType == MemType::Float)
    for (size_t i = 0; i < gridsize; ++i) field3d.vec_f[gridsize * nlevels + i] = field2d.vec_f[i] / PlanetGrav;
  else
    for (size_t i = 0; i < gridsize; ++i) field3d.vec_d[gridsize * nlevels + i] = field2d.vec_d[i] / PlanetGrav;
}

static void
vct_to_hybrid_pressure(MemType memType, Field3D &pressure_FL, Field3D &pressure_HL, const double *vct, const Field &ps,
                       long numHybridLevels, long ngp)
{
  if (memType == MemType::Float)
    vct_to_hybrid_pressure(pressure_FL.vec_f.data(), pressure_HL.vec_f.data(), vct, ps.vec_f.data(), numHybridLevels, ngp);
  else
    vct_to_hybrid_pressure(pressure_FL.vec_d.data(), pressure_HL.vec_d.data(), vct, ps.vec_d.data(), numHybridLevels, ngp);
}

static void
invert_vct(Varray<double> &vct)
{
  Varray<double> vctbuf = vct;
  auto vctSize = vct.size();
  for (size_t i = 0; i < vctSize / 2; ++i)
    {
      vct[vctSize / 2 - 1 - i] = vctbuf[i];
      vct[vctSize - 1 - i] = vctbuf[i + vctSize / 2];
    }
}

static void
check_vct(const Varray<double> &vct, int numHalfLevels)
{
  auto sum = varray_sum(2 * numHalfLevels, vct);
  if (!(sum > 0.0)) cdo_warning("VCT is empty!");
}

static void
check_range_ps(int stepNum, const Field &psProg)
{
  auto mm = field_min_max(psProg);
  if (mm.min < MIN_PS || mm.max > MAX_PS)
    cdo_warning("Surface pressure out of range (min=%g max=%g) [timestep:%d]!", mm.min, mm.max, stepNum);
}

static void
check_range_sgeopot(int stepNum, const Field &sgeopot)
{
  auto mm = field_min_max(sgeopot);
  if (mm.min < MIN_FIS || mm.max > MAX_FIS)
    cdo_warning("Surface geopotential out of range (min=%g max=%g) [timestep:%d]!", mm.min, mm.max, stepNum);
  if (sgeopot.gridsize > 1 && mm.min >= 0.0 && mm.max <= 9000.0 && IS_NOT_EQUAL(mm.min, mm.max))
    cdo_warning("Surface geopotential has an unexpected range (min=%g max=%g) [timestep:%d]!", mm.min, mm.max, stepNum);
}

static bool
zaxis_is_hybrid(int zaxistype)
{
  return (zaxistype == ZAXIS_HYBRID || zaxistype == ZAXIS_HYBRID_HALF);
}

static void
change_hybrid_zaxis(int vlistID1, int vlistID2, int vctSize, double *vct, int zaxisID2, int numFullLevels, int numHalfLevels)
{
  auto nzaxis = vlistNzaxis(vlistID1);
  for (int iz = 0; iz < nzaxis; ++iz)
    {
      auto zaxisID = vlistZaxis(vlistID1, iz);
      auto nlevels = zaxisInqSize(zaxisID);
      auto zaxistype = zaxisInqType(zaxisID);

      if (zaxis_is_hybrid(zaxistype) && (nlevels == numHalfLevels || nlevels == numFullLevels))
        {
          auto vctSize2 = zaxisInqVctSize(zaxisID);
          if (vctSize2 == vctSize && memcmp(vct, zaxisInqVctPtr(zaxisID), vctSize * sizeof(double)) == 0)
            vlistChangeZaxisIndex(vlistID2, iz, zaxisID2);
        }
    }
}

static void
pressure_level_interpolation(Varray<double> &pressureLevels, bool useHeightLevel, bool extrapolate)
{
  int numPL = pressureLevels.size();

  auto streamID1 = cdo_open_read(0);

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto vlistID2 = vlistDuplicate(vlistID1);

  auto taxisID1 = vlistInqTaxis(vlistID1);
  auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  auto gridsize = vlist_check_gridsize(vlistID1);

  auto zaxisID_PL = zaxisCreate(useHeightLevel ? ZAXIS_ALTITUDE : ZAXIS_PRESSURE, numPL);
  zaxisDefLevels(zaxisID_PL, pressureLevels.data());

  int zaxisID_ML = -1;
  int numHybridLevels = 0, numFullLevels = 0, numHalfLevels = 0;
  auto vct = vlist_read_vct(vlistID1, zaxisID_ML, numHybridLevels, numFullLevels, numHalfLevels);
  int vctSize = vct.size();

  // check VCT
  if (zaxisID_ML != -1) check_vct(vct, numHalfLevels);

  change_hybrid_zaxis(vlistID1, vlistID2, vctSize, vct.data(), zaxisID_PL, numFullLevels, numHalfLevels);

  int psvarID = -1;
  auto vctIsInverted = false;
  if (vctSize && vctSize % 2 == 0)
    {
      psvarID = vlist_get_psvarid(vlistID1, zaxisID_ML);

      int i;
      for (i = vctSize / 2 + 1; i < vctSize; ++i)
        if (vct[i] > vct[i - 1]) break;
      if (i == vctSize) vctIsInverted = true;
    }

  if (Options::cdoVerbose) cdo_print("vctIsInverted = %d", static_cast<int>(vctIsInverted));

  if (vctIsInverted) invert_vct(vct);

  VarList varList1;
  varListInit(varList1, vlistID1);
  varListSetUniqueMemtype(varList1);
  auto memType = varList1[0].memType;

  VarList varList2;
  varListInit(varList2, vlistID2);
  varListSetMemtype(varList2, memType);

  auto nvars = vlistNvars(vlistID1);

  std::vector<bool> processVars(nvars), interpVars(nvars);
  Varray2D<size_t> varnmiss(nvars);
  Field3DVector vardata1(nvars), vardata2(nvars);

  auto maxLevels = std::max(numHalfLevels, numPL);

  Varray<size_t> pnmiss;
  if (!extrapolate) pnmiss.resize(numPL);

  // check levels
  if (zaxisID_ML != -1)
    {
      auto nlev = zaxisInqSize(zaxisID_ML);
      if (nlev != numHybridLevels) cdo_abort("Internal error, wrong number of hybrid level!");
    }

  std::vector<int> vertIndex;

  Field3D pressure_FL, pressure_HL;
  if (zaxisID_ML != -1 && gridsize > 0)
    {
      vertIndex.resize(gridsize * numPL);

      CdoVar var3Df, var3Dh;
      var3Df.gridsize = gridsize;
      var3Df.nlevels = numFullLevels;
      var3Df.memType = memType;
      pressure_FL.init(var3Df);

      var3Dh.gridsize = gridsize;
      var3Dh.nlevels = numHalfLevels;
      var3Dh.memType = memType;
      pressure_HL.init(var3Dh);
    }
  else
    cdo_warning("No 3D variable with hybrid sigma pressure coordinate found!");

  if (useHeightLevel)
    {
      std::vector<double> phlev(numPL);
      height_to_pressure(phlev.data(), pressureLevels.data(), numPL);

      if (Options::cdoVerbose)
        for (int i = 0; i < numPL; ++i) cdo_print("level=%d  height=%g  pressure=%g", i + 1, pressureLevels[i], phlev[i]);

      pressureLevels = phlev;
    }

  VarIDs varIDs = search_varIDs(varList1, vlistID1, numFullLevels);

  if (Options::cdoVerbose)
    {
      cdo_print("Found:");
      // clang-format off
      if (-1 != varIDs.tempID)    cdo_print("  %s -> %s", var_stdname(air_temperature), varList1[varIDs.tempID].name);
      if (-1 != varIDs.psID)      cdo_print("  %s -> %s", var_stdname(surface_air_pressure), varList1[varIDs.psID].name);
      if (-1 != varIDs.lnpsID)    cdo_print("  LOG(%s) -> %s", var_stdname(surface_air_pressure), varList1[varIDs.lnpsID].name);
      if (-1 != varIDs.sgeopotID) cdo_print("  %s -> %s", var_stdname(surface_geopotential), varList1[varIDs.sgeopotID].name);
      if (-1 != varIDs.geopotID)  cdo_print("  %s -> %s", var_stdname(geopotential), varList1[varIDs.geopotID].name);
      if (-1 != varIDs.gheightID) cdo_print("  %s -> %s", var_stdname(geopotential_height), varList1[varIDs.gheightID].name);
      // clang-format on
    }

  for (int varID = 0; varID < nvars; ++varID)
    {
      auto &var1 = varList1[varID];
      auto gridID = var1.gridID;
      auto zaxisID = var1.zaxisID;
      auto zaxistype = zaxisInqType(zaxisID);
      auto nlevels = var1.nlevels;

      if (gridInqType(gridID) == GRID_SPECTRAL && zaxis_is_hybrid(zaxistype))
        cdo_abort("Spectral data on model level unsupported!");

      if (gridInqType(gridID) == GRID_SPECTRAL) cdo_abort("Spectral data unsupported!");

      if (varID == varIDs.gheightID) var1.nlevels = nlevels + 1;
      vardata1[varID].init(var1);
      if (varID == varIDs.gheightID) var1.nlevels = nlevels;

      // interpVars[varID] = (zaxis_is_hybrid(zaxistype) && zaxisID_ML != -1 && nlevels == numHybridLevels);
      interpVars[varID]
          = (zaxisID == zaxisID_ML
             || (zaxis_is_hybrid(zaxistype) && zaxisID_ML != -1 && (nlevels == numHalfLevels || nlevels == numFullLevels)));

      if (interpVars[varID])
        {
          varnmiss[varID].resize(maxLevels, 0);
          vardata2[varID].init(varList2[varID]);
        }
      else
        {
          varnmiss[varID].resize(nlevels);
          if (zaxis_is_hybrid(zaxistype) && zaxisID_ML != -1 && nlevels > 1)
            cdo_warning("Parameter %d has wrong number of levels, skipped! (param=%s nlevel=%d)", varID + 1, var1.name, nlevels);
        }
    }

  if (zaxisID_ML != -1 && varIDs.gheightID != -1 && varIDs.tempID == -1)
    cdo_abort("%s not found, needed for vertical interpolation of %s!", var_stdname(air_temperature),
              var_stdname(geopotential_height));

  auto presID = (psvarID != -1) ? psvarID : varIDs.lnpsID;

  if (zaxisID_ML != -1 && presID == -1)
    {
      if (varIDs.psID == -1) cdo_abort("%s not found!", var_stdname(surface_air_pressure));
      presID = varIDs.psID;
    }

  if (Options::cdoVerbose && presID != -1)
    {
      if (presID == varIDs.lnpsID)
        cdo_print("Using LOG(%s) from %s", var_stdname(surface_air_pressure), varList1[presID].name);
      else
        cdo_print("Using %s from %s", var_stdname(surface_air_pressure), varList1[presID].name);
    }

  Field psProg;
  if (zaxisID_ML != -1 && presID != -1) psProg.init(varList1[presID]);

  auto sgeopotNeeded = (varIDs.tempID != -1 || varIDs.gheightID != -1);

  Field sgeopot;
  if (zaxisID_ML != -1 && sgeopotNeeded)
    {
      sgeopot.init(varList1[presID]);
      if (varIDs.sgeopotID == -1)
        {
          if (extrapolate)
            {
              if (varIDs.geopotID == -1)
                cdo_warning("%s not found - set to zero!", var_stdname(surface_geopotential));
              else
                cdo_print("%s not found - using bottom layer of %s!", var_stdname(surface_geopotential), var_stdname(geopotential));
            }
          field_fill(sgeopot, 0.0);
        }
    }

  auto streamID2 = cdo_open_write(1);
  cdo_def_vlist(streamID2, vlistID2);

  int tsID = 0;
  while (true)
    {
      auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
      if (nrecs == 0) break;

      for (int varID = 0; varID < nvars; ++varID) processVars[varID] = false;

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_def_timestep(streamID2, tsID);

      for (int recID = 0; recID < nrecs; ++recID)
        {
          int varID, levelID;
          cdo_inq_record(streamID1, &varID, &levelID);
          const auto &var = varList1[varID];

          if (vctIsInverted && zaxisID_ML != -1 && var.zaxisID == zaxisID_ML) levelID = var.nlevels - 1 - levelID;

          cdo_read_record(streamID1, vardata1[varID], levelID, &varnmiss[varID][levelID]);

          processVars[varID] = true;
        }

      if (zaxisID_ML != -1)
        {
          if (sgeopotNeeded)
            {
              if (varIDs.sgeopotID != -1)
                field_copy(vardata1[varIDs.sgeopotID], sgeopot);
              else if (varIDs.geopotID != -1)
                field_copy(vardata1[varIDs.geopotID], numFullLevels - 1, sgeopot);

              // check range of surface geopot
              if (extrapolate && (varIDs.sgeopotID != -1 || varIDs.geopotID != -1)) check_range_sgeopot(tsID + 1, sgeopot);
            }

          if (presID == varIDs.lnpsID)
            field_transform(vardata1[varIDs.lnpsID], psProg, unary_op_exp);
          else if (presID != -1)
            field_copy(vardata1[presID], psProg);

          // check range of psProg
          check_range_ps(tsID + 1, psProg);

          vct_to_hybrid_pressure(memType, pressure_FL, pressure_HL, vct.data(), psProg, numFullLevels, gridsize);

          gen_vert_index(vertIndex, pressureLevels, pressure_FL, gridsize);

          if (!extrapolate) gen_vert_index_mv(vertIndex, pressureLevels, gridsize, psProg, pnmiss);
        }

      for (int varID = 0; varID < nvars; ++varID)
        {
          if (processVars[varID])
            {
              const auto &var = varList1[varID];

              if (tsID > 0 && var.isConstant) continue;

              if (interpVars[varID])
                {
                  auto nlevels = var.nlevels;
                  if (nlevels != numFullLevels && nlevels != numHalfLevels)
                    cdo_abort("Number of hybrid level differ from full/half level (param=%s)!", var.name);

                  for (int levelID = 0; levelID < nlevels; ++levelID)
                    {
                      if (varnmiss[varID][levelID]) cdo_abort("Missing values unsupported for this operator!");
                    }

                  if (varID == varIDs.tempID)
                    {
                      if (nlevels == numHalfLevels) cdo_abort("Temperature on half level unsupported!");

                      vertical_interp_T(nlevels, pressure_FL, pressure_HL, vardata1[varID], vardata2[varID], sgeopot, vertIndex,
                                        pressureLevels, gridsize);
                    }
                  else if (varID == varIDs.gheightID)
                    {
                      field_copy_div_2d_to_3d(memType, gridsize, nlevels, sgeopot, vardata1[varID]);

                      vertical_interp_Z(nlevels, pressure_FL, pressure_HL, vardata1[varID], vardata2[varID],
                                        vardata1[varIDs.tempID], sgeopot, vertIndex, pressureLevels, gridsize);
                    }
                  else
                    {
                      const auto &levels3D = (nlevels == numFullLevels) ? pressure_FL : pressure_HL;
                      vertical_interp_X(levels3D, vardata1[varID], vardata2[varID], vertIndex, pressureLevels, gridsize);
                    }

                  if (!extrapolate) varray_copy(numPL, pnmiss, varnmiss[varID]);
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
}

static void
height_level_interpolation(Varray<double> &heightLevels, bool extrapolate)
{
  int numHL = heightLevels.size();

  auto streamID1 = cdo_open_read(0);

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto vlistID2 = vlistDuplicate(vlistID1);

  auto taxisID1 = vlistInqTaxis(vlistID1);
  auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  auto gridsize = vlist_check_gridsize(vlistID1);

  auto zaxisID_HL = zaxisCreate(ZAXIS_HEIGHT, numHL);
  zaxisDefLevels(zaxisID_HL, heightLevels.data());

  int zaxisID_ML = -1;
  int numHybridLevels = 0, numFullLevels = 0, numHalfLevels = 0;
  auto vct = vlist_read_vct(vlistID1, zaxisID_ML, numHybridLevels, numFullLevels, numHalfLevels);
  int vctSize = vct.size();

  // check VCT
  if (zaxisID_ML != -1) check_vct(vct, numHalfLevels);

  change_hybrid_zaxis(vlistID1, vlistID2, vctSize, vct.data(), zaxisID_HL, numFullLevels, numHalfLevels);

  auto vctIsInverted = false;
  if (vctSize && vctSize % 2 == 0)
    {
      (void) vlist_get_psvarid(vlistID1, zaxisID_ML);

      int i;
      for (i = vctSize / 2 + 1; i < vctSize; ++i)
        if (vct[i] > vct[i - 1]) break;
      if (i == vctSize) vctIsInverted = true;
    }

  if (Options::cdoVerbose) cdo_print("vctIsInverted = %d", static_cast<int>(vctIsInverted));

  if (vctIsInverted) invert_vct(vct);

  VarList varList1;
  varListInit(varList1, vlistID1);
  varListSetUniqueMemtype(varList1);
  auto memType = varList1[0].memType;

  VarList varList2;
  varListInit(varList2, vlistID2);
  varListSetMemtype(varList2, memType);

  auto nvars = vlistNvars(vlistID1);

  std::vector<bool> processVars(nvars), interpVars(nvars);
  Varray2D<size_t> varnmiss(nvars);
  Field3DVector vardata1(nvars), vardata2(nvars);

  auto maxLevels = std::max(numHalfLevels, numHL);

  Varray<size_t> pnmiss;
  if (!extrapolate) pnmiss.resize(numHL);

  // check levels
  if (zaxisID_ML != -1)
    {
      auto nlev = zaxisInqSize(zaxisID_ML);
      if (nlev != numHybridLevels) cdo_abort("Internal error, wrong number of hybrid level!");
    }

  std::vector<int> vertIndex;

  Field3D pressure_FL, pressure_HL;
  if (zaxisID_ML != -1 && gridsize > 0)
    {
      vertIndex.resize(gridsize * numHL);

      CdoVar var3Df, var3Dh;
      var3Df.gridsize = gridsize;
      var3Df.nlevels = numFullLevels;
      var3Df.memType = memType;
      pressure_FL.init(var3Df);

      var3Dh.gridsize = gridsize;
      var3Dh.nlevels = numHalfLevels;
      var3Dh.memType = memType;
      pressure_HL.init(var3Dh);
    }
  else
    cdo_warning("No 3D variable with hybrid sigma pressure coordinate found!");

  VarIDs varIDs = search_varIDs(varList1, vlistID1, numFullLevels);

  if (Options::cdoVerbose)
    {
      cdo_print("Found:");
      // clang-format off
      if (-1 != varIDs.tempID)    cdo_print("  %s -> %s", var_stdname(air_temperature), varList1[varIDs.tempID].name);
      if (-1 != varIDs.psID)      cdo_print("  %s -> %s", var_stdname(surface_air_pressure), varList1[varIDs.psID].name);
      if (-1 != varIDs.lnpsID)    cdo_print("  LOG(%s) -> %s", var_stdname(surface_air_pressure), varList1[varIDs.lnpsID].name);
      if (-1 != varIDs.sgeopotID) cdo_print("  %s -> %s", var_stdname(surface_geopotential), varList1[varIDs.sgeopotID].name);
      if (-1 != varIDs.geopotID)  cdo_print("  %s -> %s", var_stdname(geopotential), varList1[varIDs.geopotID].name);
      if (-1 != varIDs.gheightID) cdo_print("  %s -> %s", var_stdname(geopotential_height), varList1[varIDs.gheightID].name);
      // clang-format on
    }

  for (int varID = 0; varID < nvars; ++varID)
    {
      auto &var1 = varList1[varID];
      auto gridID = var1.gridID;
      auto zaxisID = var1.zaxisID;
      auto zaxistype = zaxisInqType(zaxisID);
      auto nlevels = var1.nlevels;

      if (gridInqType(gridID) == GRID_SPECTRAL && zaxis_is_hybrid(zaxistype))
        cdo_abort("Spectral data on model level unsupported!");

      if (gridInqType(gridID) == GRID_SPECTRAL) cdo_abort("Spectral data unsupported!");

      vardata1[varID].init(var1);

      // interpVars[varID] = (zaxis_is_hybrid(zaxistype) && zaxisID_ML != -1 && nlevels == numHybridLevels);
      interpVars[varID]
          = (zaxisID == zaxisID_ML
             || (zaxis_is_hybrid(zaxistype) && zaxisID_ML != -1 && (nlevels == numHalfLevels || nlevels == numFullLevels)));

      if (interpVars[varID])
        {
          varnmiss[varID].resize(maxLevels, 0);
          vardata2[varID].init(varList2[varID]);
        }
      else
        {
          varnmiss[varID].resize(nlevels);
          if (zaxis_is_hybrid(zaxistype) && zaxisID_ML != -1 && nlevels > 1)
            cdo_warning("Parameter %d has wrong number of levels, skipped! (param=%s nlevel=%d)", varID + 1, var1.name, nlevels);
        }
    }

  if (zaxisID_ML != -1 && varIDs.gheightID == -1) cdo_abort("%s not found!", var_stdname(geopotential_height));

  auto sgeopotNeeded = (!extrapolate && varIDs.gheightID != -1);

  Field sgeopot;
  if (zaxisID_ML != -1 && sgeopotNeeded)
    {
      if (varIDs.sgeopotID == -1)
        {
          CdoVar var2D;
          var2D.gridsize = gridsize;
          var2D.nlevels = 1;
          var2D.memType = memType;
          sgeopot.init(var2D);

          if (extrapolate) cdo_warning("%s not found - set to zero!", var_stdname(surface_geopotential));
        }
      else { sgeopot.init(varList1[varIDs.sgeopotID]); }

      field_fill(sgeopot, 0.0);
    }

  auto streamID2 = cdo_open_write(1);
  cdo_def_vlist(streamID2, vlistID2);

  int tsID = 0;
  while (true)
    {
      auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
      if (nrecs == 0) break;

      for (int varID = 0; varID < nvars; ++varID) processVars[varID] = false;

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_def_timestep(streamID2, tsID);

      for (int recID = 0; recID < nrecs; ++recID)
        {
          int varID, levelID;
          cdo_inq_record(streamID1, &varID, &levelID);

          cdo_read_record(streamID1, vardata1[varID], levelID, &varnmiss[varID][levelID]);

          processVars[varID] = true;
        }

      if (zaxisID_ML != -1)
        {
          gen_vert_index(vertIndex, heightLevels, vardata1[varIDs.gheightID], gridsize);

          // if (!extrapolate) gen_vert_index_mv(vertIndex, heightLevels, gridsize, psProg, pnmiss);
        }

      for (int varID = 0; varID < nvars; ++varID)
        {
          if (processVars[varID])
            {
              const auto &var = varList1[varID];
              if (tsID > 0 && var.isConstant) continue;

              if (interpVars[varID])
                {
                  if (var.nlevels != numFullLevels && var.nlevels != numHalfLevels)
                    cdo_abort("Number of hybrid level differ from full/half level (param=%s)!", var.name);

                  for (int levelID = 0; levelID < var.nlevels; ++levelID)
                    {
                      if (varnmiss[varID][levelID]) cdo_abort("Missing values unsupported for this operator!");
                    }

                  // const auto &levels3D = (nlevels == numFullLevels) ? vardata1[varIDs.gheightID] : pressure_HL;
                  const auto &levels3D = vardata1[varIDs.gheightID];
                  vertical_interp_X(levels3D, vardata1[varID], vardata2[varID], vertIndex, heightLevels, gridsize);

                  if (!extrapolate) varray_copy(numHL, pnmiss, varnmiss[varID]);
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
}

void *
Vertintml(void *process)
{
  enum
  {
    func_pl,
    func_hl
  };

  cdo_initialize(process);

  // clang-format off
                    cdo_operator_add("ml2pl",      func_pl, 0, "pressure levels in pascal");
  auto ML2PLX     = cdo_operator_add("ml2plx",     func_pl, 0, "pressure levels in pascal");
                    cdo_operator_add("ml2hl",      func_hl, 0, "height levels in meter");
  auto ML2HLX     = cdo_operator_add("ml2hlx",     func_hl, 0, "height levels in meter");
                    cdo_operator_add("ml2height",  func_hl, 1, "height levels in meter");
  auto ML2HEIGHTX = cdo_operator_add("ml2heightx", func_hl, 1, "height levels in meter");
  // clang-format on

  auto operatorID = cdo_operator_id();
  auto useHeightLevel = (cdo_operator_f1(operatorID) == func_hl);
  auto doPressureInterpolation = (cdo_operator_f2(operatorID) == 0);

  auto extrapolate = (operatorID == ML2PLX || operatorID == ML2HLX || operatorID == ML2HEIGHTX);
  if (extrapolate == false) extrapolate = getenv_extrapolate();

  operator_input_arg(cdo_operator_enter(operatorID));

  Varray<double> levels;
  if (cdo_operator_argc() == 1 && cdo_operator_argv(0) == "default")
    {
      if (useHeightLevel)
        levels = { 10, 50, 100, 500, 1000, 5000, 10000, 15000, 20000, 25000, 30000 };
      else
        levels
            = { 100000, 92500, 85000, 70000, 60000, 50000, 40000, 30000, 25000, 20000, 15000, 10000, 7000, 5000, 3000, 2000, 1000 };
    }
  else { levels = cdo_argv_to_flt(cdo_get_oper_argv()); }

  if (doPressureInterpolation)
    pressure_level_interpolation(levels, useHeightLevel, extrapolate);
  else
    height_level_interpolation(levels, extrapolate);

  cdo_finish();

  return nullptr;
}
