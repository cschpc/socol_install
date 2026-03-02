/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <algorithm>
#include <cdi.h>

#include "process_int.h"
#include "cdo_vlist.h"
#include "param_conversion.h"
#include "util_string.h"

/* ================================================= */
/* LayerCloud calculates random overlap cloud cover */
/* ================================================= */

static void
layer_cloud(const Varray<double> &cc, Varray<double> &ll, long maxLevIndex, long minLevIndex, long dimgp)
{
  constexpr double ZEPSEC = 1.0 - 1.0e-12;

  for (long i = 0; i < dimgp; ++i) ll[i] = 1.0 - cc[i + maxLevIndex * dimgp];

  for (long k = maxLevIndex + 1; k <= minLevIndex; ++k)
    {
      for (long i = 0; i < dimgp; ++i)
        {
          double maxval = std::max(cc[i + (k - 1) * dimgp], cc[i + k * dimgp]);
          double minval = std::min(cc[i + (k - 1) * dimgp], ZEPSEC);
          ll[i] *= (1.0 - maxval) / (1.0 - minval);
        }
    }

  for (long i = 0; i < dimgp; ++i) ll[i] = 1.0 - ll[i];
}

static void
vct2plev(const Varray<double> &vct, Varray<double> &plevs, long nlevels)
{
  constexpr double SCALESLP = 101325.0;
  for (long k = 0; k < nlevels; ++k) plevs[k] = vct[k] + vct[k + nlevels] * SCALESLP;
}

static void
hl_index(long &maxLevIndex, long &minLevIndex, double pmax, double pmin, long nlevels, const Varray<double> &levels)
{
  maxLevIndex = -1;
  minLevIndex = -1;

  for (long k = 0; k < nlevels; ++k)
    if (levels[k] > pmax)
      {
        maxLevIndex = k - 1;
        break;
      }

  for (long k = nlevels - 1; k >= 0; --k)
    if (levels[k] < pmin)
      {
        minLevIndex = k;
        break;
      }
}

static void
pl_index(long &maxLevIndex, long &minLevIndex, double pmax, double pmin, long nlevels, const Varray<double> &levels)
{
  maxLevIndex = -1;
  minLevIndex = -1;

  for (long k = 0; k < nlevels; ++k)
    if (levels[k] >= pmax)
      {
        maxLevIndex = k;
        break;
      }

  for (long k = nlevels - 1; k >= 0; --k)
    if (levels[k] < pmin)
      {
        minLevIndex = k;
        break;
      }
}

void *
Cloudlayer(void *process)
{
  constexpr int NumVars = 3;
  int gridID, zaxisID;
  bool zrev = false;
  int aclcacID = -1;
  int nvars2 = 0;
  int aclcac_code_found = 0;
  long kmin[NumVars] = { -1, -1, -1 }, kmax[NumVars] = { -1, -1, -1 };
  double sfclevel = 0;
  double pmin = 0, pmax = 0;

  cdo_initialize(process);

  if (cdo_operator_argc() > 0)
    {
      operator_check_argc(2);
      nvars2 = 1;
      pmin = parameter_to_double(cdo_operator_argv(0));
      pmax = parameter_to_double(cdo_operator_argv(1));
    }
  else { nvars2 = NumVars; }

  auto streamID1 = cdo_open_read(0);
  auto vlistID1 = cdo_stream_inq_vlist(streamID1);

  VarList varList1;
  varListInit(varList1, vlistID1);

  auto gridsize = vlist_check_gridsize(vlistID1);

  auto aclcac_code = 223;

  auto nvars = vlistNvars(vlistID1);
  for (int varID = 0; varID < nvars; ++varID)
    {
      zaxisID = varList1[varID].zaxisID;
      auto code = varList1[varID].code;

      if (code <= 0)
        {
          if (string_to_lower(varList1[varID].name) == "aclcac") code = 223;
        }

      if (code == aclcac_code)
        {
          aclcac_code_found = 1;
          if (zaxisInqType(zaxisID) == ZAXIS_PRESSURE || zaxisInqType(zaxisID) == ZAXIS_HYBRID)
            {
              aclcacID = varID;
              break;
            }
        }
    }

  if (aclcacID == -1)
    {
      if (aclcac_code_found)
        cdo_abort("Cloud cover (parameter 223) not found on pressure or hybrid levels!");
      else
        cdo_abort("Cloud cover (parameter 223) not found!");
    }

  auto missval = varList1[aclcacID].missval;
  gridID = varList1[aclcacID].gridID;
  zaxisID = varList1[aclcacID].zaxisID;

  auto nlevels = varList1[aclcacID].nlevels;
  auto nhlev = nlevels + 1;

  Varray<double> aclcac(gridsize * nlevels);
  Varray<double> cloud[NumVars];
  for (int varID = 0; varID < nvars2; ++varID) cloud[varID].resize(gridsize);

  if (zaxisInqType(zaxisID) == ZAXIS_PRESSURE)
    {
      Varray<double> plevs(nlevels);
      zaxisInqLevels(zaxisID, plevs.data());
      if (plevs[0] > plevs[nlevels - 1])
        {
          zrev = true;
          for (int levelID = 0; levelID < nlevels / 2; ++levelID) std::swap(plevs[levelID], plevs[nlevels - 1 - levelID]);
        }
      /*
      for (int levelID = 0; levelID < nlevels; ++levelID)
        {
          printf("level %d %g\n", levelID, plevs[levelID]);
        }
      */
      if (nvars2 == 1) { pl_index(kmax[0], kmin[0], pmin, pmax, nlevels, plevs); }
      else
        {
          pl_index(kmax[2], kmin[2], 5000., 44000., nlevels, plevs);
          pl_index(kmax[1], kmin[1], 46000., 73000., nlevels, plevs);
          pl_index(kmax[0], kmin[0], 75000., 101300., nlevels, plevs);
        }
    }
  else if (zaxisInqType(zaxisID) == ZAXIS_HYBRID)
    {
      int nvct = zaxisInqVctSize(zaxisID);
      if (nlevels == (nvct / 2 - 1))
        {
          Varray<double> vct(nvct);
          zaxisInqVct(zaxisID, vct.data());

          auto nlevs = nlevels + 1;
          Varray<double> plevs(nlevs);
          vct2plev(vct, plevs, nlevs);

          if (nvars2 == 1) { hl_index(kmax[0], kmin[0], pmin, pmax, nhlev, plevs); }
          else
            {
              hl_index(kmax[2], kmin[2], 5000., 44000., nhlev, plevs);
              hl_index(kmax[1], kmin[1], 46000., 73000., nhlev, plevs);
              hl_index(kmax[0], kmin[0], 75000., 101300., nhlev, plevs);
            }
        }
      else
        cdo_abort("Unsupported vertical coordinate table format!");
    }
  else
    cdo_abort("Unsupported Z-Axis type!");

  auto surfaceID = zaxisCreate(ZAXIS_SURFACE, 1);
  zaxisDefLevels(surfaceID, &sfclevel);

  auto vlistID2 = vlistCreate();
  vlistDefNtsteps(vlistID2, vlistNtsteps(vlistID1));

  if (nvars2 == 1)
    {
      auto varID = vlistDefVar(vlistID2, gridID, surfaceID, TIME_VARYING);
      vlistDefVarParam(vlistID2, varID, cdiEncodeParam(33, 128, 255));
      cdiDefKeyString(vlistID2, varID, CDI_KEY_NAME, "cld_lay");
      cdiDefKeyString(vlistID2, varID, CDI_KEY_LONGNAME, "cloud layer");
      vlistDefVarMissval(vlistID2, varID, missval);
    }
  else
    {
      auto varID = vlistDefVar(vlistID2, gridID, surfaceID, TIME_VARYING);
      vlistDefVarParam(vlistID2, varID, cdiEncodeParam(34, 128, 255));
      cdiDefKeyString(vlistID2, varID, CDI_KEY_NAME, "low_cld");
      cdiDefKeyString(vlistID2, varID, CDI_KEY_LONGNAME, "low cloud");
      vlistDefVarMissval(vlistID2, varID, missval);

      varID = vlistDefVar(vlistID2, gridID, surfaceID, TIME_VARYING);
      vlistDefVarParam(vlistID2, varID, cdiEncodeParam(35, 128, 255));
      cdiDefKeyString(vlistID2, varID, CDI_KEY_NAME, "mid_cld");
      cdiDefKeyString(vlistID2, varID, CDI_KEY_LONGNAME, "mid cloud");
      vlistDefVarMissval(vlistID2, varID, missval);

      varID = vlistDefVar(vlistID2, gridID, surfaceID, TIME_VARYING);
      vlistDefVarParam(vlistID2, varID, cdiEncodeParam(36, 128, 255));
      cdiDefKeyString(vlistID2, varID, CDI_KEY_NAME, "hih_cld");
      cdiDefKeyString(vlistID2, varID, CDI_KEY_LONGNAME, "high cloud");
      vlistDefVarMissval(vlistID2, varID, missval);
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

          size_t offset = zrev ? (nlevels - 1 - levelID) * gridsize : levelID * gridsize;

          if (varID == aclcacID)
            {
              size_t nmiss;
              cdo_read_record(streamID1, aclcac.data() + offset, &nmiss);
              if (nmiss != 0) cdo_abort("Missing values unsupported!");
            }
        }

      for (int varID = 0; varID < nvars2; ++varID)
        {
          for (size_t i = 0; i < gridsize; ++i) cloud[varID][i] = missval;
        }

      for (int varID = 0; varID < nvars2; ++varID)
        {
          if (kmax[varID] != -1 && kmin[varID] != -1) layer_cloud(aclcac, cloud[varID], kmax[varID], kmin[varID], gridsize);
        }

      for (int varID = 0; varID < nvars2; ++varID)
        {
          auto nmiss = varray_num_mv(gridsize, cloud[varID], missval);

          cdo_def_record(streamID2, varID, 0);
          cdo_write_record(streamID2, cloud[varID].data(), nmiss);
        }

      tsID++;
    }

  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  vlistDestroy(vlistID2);

  cdo_finish();

  return nullptr;
}
