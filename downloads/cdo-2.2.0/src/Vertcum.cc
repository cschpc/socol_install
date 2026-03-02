/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Vertcum    vertcum         Vertical cumulative
      Vertcum    vertcumhl       Vertical cumulative on hybrid sigma half level
*/

#include <cdi.h>

#include "process_int.h"
#include "cdo_vlist.h"

#define IS_SURFACE_LEVEL(zaxisID) (zaxisInqType(zaxisID) == ZAXIS_SURFACE && zaxisInqSize(zaxisID) == 1)

static void
add_vars_mv(size_t gridsize, double missval, const Varray<double> &var1, const Varray<double> &var2, Varray<double> &var3)
{
  auto missval1 = missval;
  auto missval2 = missval;
  // for ( size_t i = 0; i < gridsize; ++i ) var3[i] = ADDMN(var2[i], var1[i]);
  for (size_t i = 0; i < gridsize; ++i)
    {
      var3[i] = var2[i];
      if (!DBL_IS_EQUAL(var1[i], missval1))
        {
          if (!DBL_IS_EQUAL(var2[i], missval2))
            var3[i] += var1[i];
          else
            var3[i] = var1[i];
        }
    }
}

void *
Vertcum(void *process)
{
  int nlevshl = 0;

  cdo_initialize(process);

  // clang-format off
                  cdo_operator_add("vertcum",    0,  0, nullptr);
  int VERTCUMHL = cdo_operator_add("vertcumhl",  0,  0, nullptr);
  // clang-format on

  auto operatorID = cdo_operator_id();

  operator_check_argc(0);

  auto streamID1 = cdo_open_read(0);

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto vlistID2 = vlistDuplicate(vlistID1);

  int zaxisIDhl = -1;

  if (operatorID == VERTCUMHL)
    {
      std::vector<double> vct;
      bool lhybrid = false;
      auto nzaxis = vlistNzaxis(vlistID1);
      for (int i = 0; i < nzaxis; ++i)
        {
          auto zaxisID = vlistZaxis(vlistID1, i);
          auto nlevs = zaxisInqSize(zaxisID);

          if (zaxisInqType(zaxisID) == ZAXIS_HYBRID && nlevs > 1)
            {
              int nvct = zaxisInqVctSize(zaxisID);
              if (nlevs == (nvct / 2 - 1))
                {
                  if (lhybrid == false)
                    {
                      lhybrid = true;
                      nlevshl = nlevs + 1;

                      vct.resize(nvct);
                      zaxisInqVct(zaxisID, vct.data());

                      zaxisIDhl = zaxisCreate(ZAXIS_HYBRID_HALF, nlevshl);
                      std::vector<double> levels(nlevshl);
                      for (int levelID = 0; levelID < nlevshl; ++levelID) levels[levelID] = levelID + 1;
                      zaxisDefLevels(zaxisIDhl, levels.data());
                      zaxisDefVct(zaxisIDhl, nvct, vct.data());
                      vlistChangeZaxisIndex(vlistID2, i, zaxisIDhl);
                    }
                  else if (vct.size())
                    {
                      if (memcmp(vct.data(), zaxisInqVctPtr(zaxisID), nvct * sizeof(double)) == 0)
                        vlistChangeZaxisIndex(vlistID2, i, zaxisIDhl);
                    }
                }
            }
        }
    }

  VarList varList1, varList2;
  varListInit(varList1, vlistID1);
  varListInit(varList2, vlistID2);

  auto nvars = vlistNvars(vlistID1);
  std::vector<std::vector<size_t>> varnmiss(nvars);
  Varray3D<double> vardata1(nvars), vardata2(nvars);

  for (int varID = 0; varID < nvars; ++varID)
    {
      auto gridsize = varList1[varID].gridsize;
      auto nlevs = varList1[varID].nlevels;
      auto nlevs2 = varList2[varID].nlevels;

      varnmiss[varID].resize(nlevs);
      vardata1[varID].resize(nlevs);
      vardata2[varID].resize(nlevs2);
      for (int levelID = 0; levelID < nlevs; ++levelID) vardata1[varID][levelID].resize(gridsize);
      for (int levelID = 0; levelID < nlevs2; ++levelID) vardata2[varID][levelID].resize(gridsize);
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
          cdo_read_record(streamID1, &vardata1[varID][levelID][0], &varnmiss[varID][levelID]);
        }

      for (int varID = 0; varID < nvars; ++varID)
        {
          auto missval = varList2[varID].missval;
          auto gridsize = varList2[varID].gridsize;
          auto nlevs2 = varList2[varID].nlevels;

          if (operatorID == VERTCUMHL && nlevs2 == nlevshl)
            {
              for (size_t i = 0; i < gridsize; ++i) vardata2[varID][0][i] = 0;
            }
          else
            {
              for (size_t i = 0; i < gridsize; ++i) vardata2[varID][0][i] = vardata1[varID][0][i];
            }

          for (int levelID = 1; levelID < nlevs2; ++levelID)
            {
              if (operatorID == VERTCUMHL && nlevs2 == nlevshl)
                add_vars_mv(gridsize, missval, vardata1[varID][levelID - 1], vardata2[varID][levelID - 1],
                            vardata2[varID][levelID]);
              else
                add_vars_mv(gridsize, missval, vardata1[varID][levelID], vardata2[varID][levelID - 1], vardata2[varID][levelID]);
            }

          if (operatorID == VERTCUMHL && nlevs2 == nlevshl)
            {
              const auto &var1 = vardata2[varID][nlevs2 - 1];
              for (int levelID = 0; levelID < nlevs2; ++levelID)
                {
                  auto &var2 = vardata2[varID][levelID];
                  for (size_t i = 0; i < gridsize; ++i)
                    {
                      if (IS_NOT_EQUAL(var1[i], 0))
                        var2[i] /= var1[i];
                      else
                        var2[i] = 0;
                    }
                }
            }
        }

      for (int varID = 0; varID < nvars; ++varID)
        {
          auto missval = varList2[varID].missval;
          auto gridsize = varList2[varID].gridsize;
          auto nlevs2 = varList2[varID].nlevels;
          for (int levelID = 0; levelID < nlevs2; ++levelID)
            {
              auto &single = vardata2[varID][levelID];
              auto nmiss = varray_num_mv(gridsize, single, missval);
              cdo_def_record(streamID2, varID, levelID);
              cdo_write_record(streamID2, single.data(), nmiss);
            }
        }

      tsID++;
    }

  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  vlistDestroy(vlistID2);

  cdo_finish();

  return nullptr;
}
