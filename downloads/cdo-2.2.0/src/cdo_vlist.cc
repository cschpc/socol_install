/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <array>
#include <algorithm>

#include <cdi.h>
#include "cdo_cdi_wrapper.h"
#include "cdo_options.h"
#include "util_string.h"
#include "cdo_vlist.h"
#include "compare.h"
#include "cdo_output.h"
#include "field_functions.h"

void
vlist_define_timestep_type(int vlistID, int operfunc)
{
  int stepType = -1;
  // clang-format off
  if      (operfunc == FieldFunc_Mean)  stepType = TSTEP_AVG;
  else if (operfunc == FieldFunc_Avg)   stepType = TSTEP_AVG;
  else if (operfunc == FieldFunc_Sum)   stepType = TSTEP_SUM;
  else if (operfunc == FieldFunc_Range) stepType = TSTEP_RANGE;
  else if (operfunc == FieldFunc_Min)   stepType = TSTEP_MIN;
  else if (operfunc == FieldFunc_Max)   stepType = TSTEP_MAX;
  // clang-format on

  if (stepType != -1)
    {
      auto nvars = vlistNvars(vlistID);
      for (int varID = 0; varID < nvars; ++varID) vlistDefVarTsteptype(vlistID, varID, stepType);
    }
}

double
cdo_zaxis_inq_level(int zaxisID, int levelID)
{
  auto zaxistype = zaxisInqType(zaxisID);
  return zaxisInqLevels(zaxisID, nullptr) ? zaxisInqLevel(zaxisID, levelID) : (zaxistype == ZAXIS_SURFACE) ? 0.0 : levelID + 1.0;
}

int
cdo_zaxis_inq_levels(int zaxisID, double *levels)
{
  auto size = zaxisInqLevels(zaxisID, nullptr);

  if (levels)
    {
      if (size)
        zaxisInqLevels(zaxisID, levels);
      else
        {
          size = zaxisInqSize(zaxisID);
          if (size == 1 && zaxisInqType(zaxisID) == ZAXIS_SURFACE)
            levels[0] = 0.0;
          else
            for (int i = 0; i < size; ++i) levels[i] = i + 1.0;
        }
    }

  return size;
}

static void
compare_lat_reg2d(size_t ysize, int gridID1, int gridID2)
{
  if (ysize > 1)
    {
      Varray<double> yvals1(ysize), yvals2(ysize);
      auto ny1 = gridInqYvals(gridID1, &yvals1[0]);
      auto ny2 = gridInqYvals(gridID2, &yvals2[0]);
      if (ny1 == 0 || ny2 == 0) return;

      if (IS_EQUAL(yvals1[0], yvals2[ysize - 1]) && IS_EQUAL(yvals1[ysize - 1], yvals2[0]))
        {
          if (yvals1[0] > yvals2[0])
            cdo_abort("Latitude orientation differ! First grid: N->S; second grid: S->N");
          else
            cdo_abort("Latitude orientation differ! First grid: S->N; second grid: N->S");
        }
      else
        {
          for (size_t i = 0; i < ysize; ++i)
            if (std::fabs(yvals1[i] - yvals2[i]) > 3.e-5)
              {
                cdo_warning("Grid latitudes differ!");
                break;
              }
        }
    }
}

static void
compare_lon_reg2d(size_t xsize, int gridID1, int gridID2)
{
  if (xsize > 1)
    {
      Varray<double> xvals1(xsize), xvals2(xsize);
      auto nx1 = gridInqXvals(gridID1, &xvals1[0]);
      auto nx2 = gridInqXvals(gridID2, &xvals2[0]);
      if (nx1 == 0 || nx2 == 0) return;

      for (size_t i = 0; i < xsize; ++i)
        if (std::fabs(xvals1[i] - xvals2[i]) > 3.e-5)
          {
            cdo_warning("Grid longitudes differ!");
            break;
          }
    }
}

static void
compare_grid_unstructured(int gridID1, int gridID2)
{
  if (gridInqXvals(gridID1, nullptr) && gridInqXvals(gridID1, nullptr) == gridInqXvals(gridID2, nullptr)
      && gridInqYvals(gridID1, nullptr) && gridInqYvals(gridID1, nullptr) == gridInqYvals(gridID2, nullptr))
    {
      auto gridsize = gridInqSize(gridID1);
      Varray<double> xvals1(gridsize), yvals1(gridsize), xvals2(gridsize), yvals2(gridsize);
      gridInqXvals(gridID1, xvals1.data());
      gridInqYvals(gridID1, yvals1.data());
      gridInqXvals(gridID2, xvals2.data());
      gridInqYvals(gridID2, yvals2.data());

      // size_t inc = (gridsize > 10000) ? gridsize / 1000 : 1;
      constexpr size_t inc = 1;
      for (size_t i = 0; i < gridsize; i += inc)
        if (std::fabs(xvals1[i] - xvals2[i]) > 2.e-5 || std::fabs(yvals1[i] - yvals2[i]) > 2.e-5)
          {
            cdo_warning("Geographic location of some grid points differ!");
            if (Options::cdoVerbose)
              printf("cell=%zu x1=%g x2=%g y1=%g y2=%g dx=%g dy=%g\n", i + 1, xvals1[i], xvals2[i], yvals1[i], yvals2[i],
                     xvals1[i] - xvals2[i], yvals1[i] - yvals2[i]);
            break;
          }
    }
}

void
cdo_compare_grids(int gridID1, int gridID2)
{
  if (gridID1 == gridID2) return;

  // compare grids of first variable

  auto gridType1 = gridInqType(gridID1);
  auto gridType2 = gridInqType(gridID2);
  if (gridType1 == gridType2)
    {
      if (gridType1 == GRID_GAUSSIAN || gridType2 == GRID_LONLAT)
        {
          auto xsize = gridInqXsize(gridID1);
          auto ysize = gridInqYsize(gridID1);

          if (ysize == gridInqYsize(gridID2))
            compare_lat_reg2d(ysize, gridID1, gridID2);
          else
            cdo_warning("ysize of input grids differ!");

          if (xsize == gridInqXsize(gridID2))
            compare_lon_reg2d(xsize, gridID1, gridID2);
          else
            cdo_warning("xsize of input grids differ!");
        }
      else if (gridType1 == GRID_CURVILINEAR || gridType2 == GRID_UNSTRUCTURED) { compare_grid_unstructured(gridID1, gridID2); }
    }
  else if (gridInqSize(gridID1) > 1)
    {
      cdo_warning("Grids have different types! First grid: %s; second grid: %s", gridNamePtr(gridType1), gridNamePtr(gridType2));
    }
}

static int
zaxisCheckLevels(int zaxisID1, int zaxisID2)
{
  if (zaxisID1 != zaxisID2)
    {
      auto nlev1 = zaxisInqSize(zaxisID1);
      auto nlev2 = zaxisInqSize(zaxisID2);
      if (nlev1 != nlev2) cdo_abort("Number of levels of the input fields do not match!");

      Varray<double> lev1(nlev1), lev2(nlev1);
      cdo_zaxis_inq_levels(zaxisID1, &lev1[0]);
      cdo_zaxis_inq_levels(zaxisID2, &lev2[0]);

      auto ldiffer = false;
      for (int i = 0; i < nlev1; ++i)
        if (IS_NOT_EQUAL(lev1[i], lev2[i]))
          {
            ldiffer = true;
            break;
          }
      if (ldiffer)
        {
          ldiffer = false;
          for (int i = 0; i < nlev1; ++i)
            if (IS_NOT_EQUAL(lev1[i], lev2[nlev1 - 1 - i]))
              {
                ldiffer = true;
                break;
              }

          if (ldiffer)
            cdo_warning("Input parameters have different levels!");
          else
            cdo_warning("Z-axis orientation differ!");

          return 1;
        }
    }

  return 0;
}

static void
vlistCheckNames(int vlistID1, int vlistID2)
{
  auto nvars = vlistNvars(vlistID1);

  std::vector<std::string> names1(nvars);
  std::vector<std::string> names2(nvars);
  for (int varID = 0; varID < nvars; ++varID) names1[varID] = cdo::inq_var_name(vlistID1, varID);
  for (int varID = 0; varID < nvars; ++varID) names2[varID] = cdo::inq_var_name(vlistID2, varID);

  std::sort(names1.begin(), names1.end());
  std::sort(names2.begin(), names2.end());

  int varID;
  for (varID = 0; varID < nvars; ++varID)
    if (names1[varID] != names2[varID]) break;

  if (varID == nvars) cdo_print("Use CDO option --sortname to sort the parameter by name (NetCDF only)!");
}

static void
vlistPrintMissingVars(int vlistID1, int vlistID2)
{
  auto nvars1 = vlistNvars(vlistID1);
  auto nvars2 = vlistNvars(vlistID2);

  std::vector<std::string> names1(nvars1);
  std::vector<std::string> names2(nvars2);
  for (int varID1 = 0; varID1 < nvars1; ++varID1) names1[varID1] = cdo::inq_var_name(vlistID1, varID1);
  for (int varID2 = 0; varID2 < nvars2; ++varID2) names2[varID2] = cdo::inq_var_name(vlistID2, varID2);

  if (nvars1 > nvars2)
    {
      for (int varID1 = 0; varID1 < nvars1; ++varID1)
        {
          int varID2;
          for (varID2 = 0; varID2 < nvars2; ++varID2)
            {
              if (names1[varID1] == names2[varID2]) break;
            }
          if (varID2 == nvars2) cdo_print("Variable %s not found in second input stream!", names1[varID1]);
        }
    }
  else
    {
      for (int varID2 = 0; varID2 < nvars2; ++varID2)
        {
          int varID1;
          for (varID1 = 0; varID1 < nvars1; ++varID1)
            {
              if (names1[varID1] == names2[varID2]) break;
            }
          if (varID1 == nvars1) cdo_print("Variable %s not found in first input stream!", names2[varID2]);
        }
    }
}

void
vlist_compare(int vlistID1, int vlistID2, int flag)
{
  auto lchecknames = false;

  auto nvars = vlistNvars(vlistID1);

  if (nvars != vlistNvars(vlistID2))
    {
      vlistPrintMissingVars(vlistID1, vlistID2);
      cdo_abort("Input streams have different number of variables per timestep!");
    }

  if (vlistNrecs(vlistID1) != vlistNrecs(vlistID2))
    cdo_abort("Input streams have different number of %s per timestep!", (nvars == 1) ? "layers" : "records");

  for (int varID = 0; varID < nvars; ++varID)
    {
      if (nvars > 1)
        {
          if (flag & CMP_NAME)
            {
              auto name1 = string_to_lower(cdo::inq_var_name(vlistID1, varID));
              auto name2 = string_to_lower(cdo::inq_var_name(vlistID2, varID));
              if (name1 != name2)
                {
                  cdo_warning("Input streams have different parameter names!");
                  lchecknames = true;
                  flag -= CMP_NAME;
                  //    break;
                }
            }
        }

      if (flag & CMP_GRIDSIZE)
        {
          if (gridInqSize(vlistInqVarGrid(vlistID1, varID)) != gridInqSize(vlistInqVarGrid(vlistID2, varID)))
            {
              cdo_abort("Grid size of the input field '%s' do not match!", cdo::inq_var_name(vlistID1, varID));
            }
        }

      if (flag & CMP_NLEVEL)
        {
          auto zaxisID1 = vlistInqVarZaxis(vlistID1, varID);
          auto zaxisID2 = vlistInqVarZaxis(vlistID2, varID);
          if (zaxisCheckLevels(zaxisID1, zaxisID2) != 0) break;
        }
    }

  if (flag & CMP_GRID)
    {
      auto gridID1 = vlistInqVarGrid(vlistID1, 0);
      auto gridID2 = vlistInqVarGrid(vlistID2, 0);
      if (gridID1 != gridID2) cdo_compare_grids(gridID1, gridID2);
    }

  if (lchecknames) vlistCheckNames(vlistID1, vlistID2);
}

int
vlist_compare_x(int vlistID1, int vlistID2, int flag)
{
  auto nvars = vlistNvars(vlistID1);
  auto nvars2 = vlistNvars(vlistID2);
  auto nlevels2 = zaxisInqSize(vlistInqVarZaxis(vlistID2, 0));

  if (nvars2 != 1) cdo_abort("Internal problem, vlist_compare_x() called with unexpected vlistID2 argument!");

  for (int varID = 0; varID < nvars; ++varID)
    {
      if (flag & CMP_GRIDSIZE)
        {
          if (gridInqSize(vlistInqVarGrid(vlistID1, varID)) != gridInqSize(vlistInqVarGrid(vlistID2, 0)))
            cdo_abort("Grid size of the input fields do not match!");
        }

      if (flag & CMP_NLEVEL)
        {
          if ((zaxisInqSize(vlistInqVarZaxis(vlistID1, varID)) != nlevels2) && nlevels2 > 1)
            cdo_abort("Number of levels of the input fields do not match!");
        }
    }

  if (flag & CMP_GRID)
    {
      auto gridID1 = vlistInqVarGrid(vlistID1, 0);
      auto gridID2 = vlistInqVarGrid(vlistID2, 0);
      if (gridID1 != gridID2) cdo_compare_grids(gridID1, gridID2);
    }

  return nlevels2;
}

void
vlist_map(int vlistID1, int vlistID2, int flag, int mapflag, std::map<int, int> &mapOfVarIDs)
{
  auto nvars1 = vlistNvars(vlistID1);
  auto nvars2 = vlistNvars(vlistID2);

  std::vector<std::string> names1(nvars1);
  std::vector<std::string> names2(nvars2);
  for (int varID1 = 0; varID1 < nvars1; ++varID1) names1[varID1] = cdo::inq_var_name(vlistID1, varID1);
  for (int varID2 = 0; varID2 < nvars2; ++varID2) names2[varID2] = cdo::inq_var_name(vlistID2, varID2);

  if (mapflag == 2)
    {
      for (int varID2 = 0; varID2 < nvars2; ++varID2)
        {
          int varID1;
          for (varID1 = 0; varID1 < nvars1; ++varID1)
            {
              if (names1[varID1] == names2[varID2]) break;
            }
          if (varID1 == nvars1) { cdo_abort("Variable %s not found in first input stream!", names2[varID2]); }
          else { mapOfVarIDs[varID1] = varID2; }
        }
    }
  else
    {
      for (int varID1 = 0; varID1 < nvars1; ++varID1)
        {
          int varID2;
          for (varID2 = 0; varID2 < nvars2; ++varID2)
            {
              if (names1[varID1] == names2[varID2]) break;
            }
          if (varID2 == nvars2)
            {
              if (mapflag == 3) continue;
              cdo_abort("Variable %s not found in second input stream!", names1[varID1]);
            }
          else { mapOfVarIDs[varID1] = varID2; }
        }
    }

  if (mapOfVarIDs.empty()) cdo_abort("No variable found that occurs in both streams!");

  if (Options::cdoVerbose)
    for (int varID1 = 0; varID1 < nvars1; ++varID1)
      {
        const auto &it = mapOfVarIDs.find(varID1);
        if (it != mapOfVarIDs.end())
          cdo_print("Variable %d:%s mapped to %d:%s", varID1, names1[varID1], it->second, names2[it->second]);
      }

  if (mapOfVarIDs.size() > 1)
    {
      auto varID2 = mapOfVarIDs.begin()->second;
      for (auto it = ++mapOfVarIDs.begin(); it != mapOfVarIDs.end(); ++it)
        {
          if (it->second < varID2)
            cdo_abort("Variable names must be sorted, use CDO option --sortname to sort the parameter by name (NetCDF only)!");

          varID2 = it->second;
        }
    }

  for (auto it = mapOfVarIDs.begin(); it != mapOfVarIDs.end(); ++it)
    {
      auto varID1 = it->first;
      auto varID2 = it->second;

      if (flag & CMP_GRIDSIZE)
        {
          if (gridInqSize(vlistInqVarGrid(vlistID1, varID1)) != gridInqSize(vlistInqVarGrid(vlistID2, varID2)))
            cdo_abort("Grid size of the input fields do not match!");
        }

      if (flag & CMP_NLEVEL)
        {
          auto zaxisID1 = vlistInqVarZaxis(vlistID1, varID1);
          auto zaxisID2 = vlistInqVarZaxis(vlistID2, varID2);
          if (zaxisCheckLevels(zaxisID1, zaxisID2) != 0) break;
        }

      if (flag & CMP_GRID && varID1 == mapOfVarIDs.begin()->first)
        {
          auto gridID1 = vlistInqVarGrid(vlistID1, varID1);
          auto gridID2 = vlistInqVarGrid(vlistID2, varID2);
          if (gridID1 != gridID2) cdo_compare_grids(gridID1, gridID2);
        }
    }
}

bool
vlist_is_szipped(int vlistID)
{
  auto nvars = vlistNvars(vlistID);
  for (int varID = 0; varID < nvars; ++varID)
    {
      auto comptype = vlistInqVarCompType(vlistID, varID);
      if (comptype == CDI_COMPRESS_SZIP) return true;
    }

  return false;
}

int
vlist_inq_nwpv(int vlistID, int varID)
{
  auto datatype = vlistInqVarDatatype(vlistID, varID);
  // number of words per value; real:1  complex:2
  auto nwpv = (datatype == CDI_DATATYPE_CPX32 || datatype == CDI_DATATYPE_CPX64) ? 2 : 1;

  return nwpv;
}

size_t
vlist_check_gridsize(int vlistID)
{
  auto lerror = false;
  auto ngp = gridInqSize(vlistGrid(vlistID, 0));

  // check gridsize
  auto ngrids = vlistNgrids(vlistID);
  for (int index = 0; index < ngrids; ++index)
    {
      auto gridID = vlistGrid(vlistID, index);
      if (ngp != gridInqSize(gridID))
        {
          lerror = true;
          break;
        }
    }

  if (lerror)
    {
      cdo_print("This operator requires all variables on the same horizontal grid.");
      cdo_print("Horizontal grids found:");
      for (int index = 0; index < ngrids; ++index)
        {
          auto gridID = vlistGrid(vlistID, index);
          cdo_print("  grid=%d  type=%s  points=%zu", index + 1, gridNamePtr(gridInqType(gridID)), gridInqSize(gridID));
        }
      cdo_abort("The input stream contains variables on different horizontal grids!");
    }

  return ngp;
}

Varray<double>
vlist_read_vct(int vlistID, int &zaxisID_ML, int &numHybridLevels, int &numFullLevels, int &numHalfLevels)
{
  Varray<double> vct;
  int vctSize = 0;
  zaxisID_ML = -1;
  numHybridLevels = 0;
  numFullLevels = 0;
  numHalfLevels = 0;

  auto haveVCT = false;
  auto nzaxis = vlistNzaxis(vlistID);
  for (int iz = 0; iz < nzaxis; ++iz)
    {
      // auto monoLevel = false;
      auto monoLevel = true;
      auto zaxisID = vlistZaxis(vlistID, iz);
      auto nlevels = zaxisInqSize(zaxisID);
      auto zaxistype = zaxisInqType(zaxisID);

      if (Options::cdoVerbose)
        cdo_print("ZAXIS_HYBRID=%d ZAXIS_HYBRID_HALF=%d nlevels=%d monoLevel=%d", zaxisInqType(zaxisID) == ZAXIS_HYBRID,
                  zaxisInqType(zaxisID) == ZAXIS_HYBRID_HALF, nlevels, monoLevel);

      if ((zaxistype == ZAXIS_HYBRID || zaxistype == ZAXIS_HYBRID_HALF) && nlevels > 1 && !monoLevel)
        {
          Varray<double> level(nlevels);
          cdo_zaxis_inq_levels(zaxisID, &level[0]);
          int l;
          for (l = 0; l < nlevels; ++l)
            {
              if ((l + 1) != (int) (level[l] + 0.5)) break;
            }
          if (l == nlevels) monoLevel = true;
        }

      if ((zaxistype == ZAXIS_HYBRID || zaxistype == ZAXIS_HYBRID_HALF) && nlevels > 1 && monoLevel)
        {
          vctSize = zaxisInqVctSize(zaxisID);
          if (nlevels == (vctSize / 2 - 1))
            {
              if (!haveVCT)
                {
                  haveVCT = true;
                  zaxisID_ML = zaxisID;
                  numHybridLevels = nlevels;
                  numFullLevels = numHybridLevels;
                  numHalfLevels = numFullLevels + 1;

                  vct.resize(vctSize);
                  zaxisInqVct(zaxisID, vct.data());
                  if (Options::cdoVerbose)
                    cdo_print("Detected half-level model definition : nlevels == (vctSize/2 - 1) (nlevels: %d, vctSize: %d, "
                              "numFullLevels: %d, numHalfLevels: %d) ",
                              nlevels, vctSize, numFullLevels, numHalfLevels);
                }
            }
          else if (nlevels == (vctSize / 2))
            {
              if (!haveVCT)
                {
                  haveVCT = true;
                  zaxisID_ML = zaxisID;
                  numHybridLevels = nlevels;
                  numFullLevels = numHybridLevels - 1;
                  numHalfLevels = numHybridLevels;

                  vct.resize(vctSize);
                  zaxisInqVct(zaxisID, vct.data());
                  if (Options::cdoVerbose)
                    cdo_print("Detected full-level model definition : nlevels == (vctSize/2) (nlevels: %d, vctSize: %d, "
                              "numFullLevels: %d, numHalfLevels: %d) ",
                              nlevels, vctSize, numFullLevels, numHalfLevels);
                }
            }
          else if (nlevels == (vctSize - 4 - 1))
            {
              if (!haveVCT)
                {
                  Varray<double> vctRead(vctSize);
                  zaxisInqVct(zaxisID, vctRead.data());

                  constexpr int voff = 4;
                  if ((int) (vctRead[0] + 0.5) == 100000 && vctRead[voff] < vctRead[voff + 1])
                    {
                      haveVCT = true;
                      zaxisID_ML = zaxisID;
                      numHybridLevels = nlevels;
                      numFullLevels = numHybridLevels;
                      numHalfLevels = numFullLevels + 1;

                      int vctsize = 2 * numHalfLevels;
                      vct.resize(vctsize);

                      // calculate VCT for LM

                      for (int i = 0; i < vctsize / 2; ++i)
                        {
                          if (vctRead[voff + i] >= vctRead[voff] && vctRead[voff + i] <= vctRead[3])
                            {
                              vct[i] = vctRead[0] * vctRead[voff + i];
                              vct[vctsize / 2 + i] = 0;
                            }
                          else
                            {
                              vct[i] = (vctRead[0] * vctRead[3] * (1 - vctRead[voff + i])) / (1 - vctRead[3]);
                              vct[vctsize / 2 + i] = (vctRead[voff + i] - vctRead[3]) / (1 - vctRead[3]);
                            }
                        }

                      if (Options::cdoVerbose)
                        {
                          for (int i = 0; i < vctsize / 2; ++i)
                            fprintf(stdout, "%5d %25.17f %25.17f\n", i, vct[i], vct[vctsize / 2 + i]);
                        }
                    }
                }
            }
        }
    }

  return vct;
}

void
vlist_change_hybrid_zaxis(int vlistID1, int vlistID2, int zaxisID1, int zaxisID2)
{
  int vctSize0 = 0;
  Varray<double> vct;

  auto nzaxis = vlistNzaxis(vlistID1);
  for (int i = 0; i < nzaxis; ++i)
    {
      auto zaxisID = vlistZaxis(vlistID1, i);
      auto nlevels = zaxisInqSize(zaxisID);

      if (zaxisID == zaxisID1 && nlevels > 1)
        {
          auto vctSize = zaxisInqVctSize(zaxisID);
          if (!vct.size())
            {
              vctSize0 = vctSize;
              vct.resize(vctSize);
              zaxisInqVct(zaxisID, vct.data());

              vlistChangeZaxisIndex(vlistID2, i, zaxisID2);
            }
          else
            {
              if (vctSize0 == vctSize && memcmp(vct.data(), zaxisInqVctPtr(zaxisID), vctSize * sizeof(double)) == 0)
                vlistChangeZaxisIndex(vlistID2, i, zaxisID2);
            }
        }
    }
}

int
vlist_get_psvarid(int vlistID, int zaxisID)
{
  auto psname = cdo::inq_key_string(zaxisID, CDI_GLOBAL, CDI_KEY_PSNAME);
  if (psname.size())
    {
      auto nvars = vlistNvars(vlistID);
      for (int varID = 0; varID < nvars; ++varID)
        {
          auto name = cdo::inq_var_name(vlistID, varID);
          if (name == psname) return varID;
        }
      if (Options::cdoVerbose) cdo_warning("Surface pressure variable not found - %s", psname);
    }

  return -1;
}

int
vlist_get_first_spectral_grid(int vlistID)
{
  // find first spectral grid
  auto ngrids = vlistNgrids(vlistID);
  for (int index = 0; index < ngrids; ++index)
    {
      auto gridID = vlistGrid(vlistID, index);
      if (gridInqType(gridID) == GRID_SPECTRAL) return gridID;
    }

  return -1;
}

int
vlist_get_first_gaussian_grid(int vlistID)
{
  // find first gaussian grid
  auto ngrids = vlistNgrids(vlistID);
  for (int index = 0; index < ngrids; ++index)
    {
      auto gridID = vlistGrid(vlistID, index);
      if (gridInqType(gridID) == GRID_GAUSSIAN) return gridID;
    }

  return -1;
}

int
vlist_get_first_fourier_grid(int vlistID)
{
  // find first fourier grid
  auto ngrids = vlistNgrids(vlistID);
  for (int index = 0; index < ngrids; ++index)
    {
      auto gridID = vlistGrid(vlistID, index);
      if (gridInqType(gridID) == GRID_FOURIER) return gridID;
    }

  return -1;
}

void
cdo_check_missval(double missval, const std::string &varname)
{
  if (dbl_is_equal(0.0, missval) || dbl_is_equal(1.0, missval))
    {
      static auto printWarning = true;
      if (printWarning)
        {
          printWarning = false;
          cdo_warning("Variable %s has a missing value of %g, this can lead to incorrect results with this operator!", varname,
                      missval);
        }
    }
}
