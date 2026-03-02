/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Invertlev     invertlev       Invert level
*/

#include <cdi.h>

#include "process_int.h"
#include "cdo_vlist.h"

static void
invertLevDes(int vlistID)
{
  auto nzaxis = vlistNzaxis(vlistID);
  for (int index = 0; index < nzaxis; ++index)
    {
      auto zaxisID1 = vlistZaxis(vlistID, index);
      auto zaxisID2 = zaxisDuplicate(zaxisID1);
      auto zaxistype = zaxisInqType(zaxisID1);

      auto nlev = zaxisInqSize(zaxisID1);
      if (nlev <= 1) continue;

      if (zaxisInqLevels(zaxisID1, nullptr))
        {
          Varray<double> yv1(nlev), yv2(nlev);
          zaxisInqLevels(zaxisID1, yv1.data());
          for (int ilev = 0; ilev < nlev; ++ilev) yv2[nlev - ilev - 1] = yv1[ilev];
          zaxisDefLevels(zaxisID2, yv2.data());
        }

      if (zaxisInqLbounds(zaxisID1, nullptr) && zaxisInqUbounds(zaxisID1, nullptr))
        {
          Varray<double> yb1(nlev), yb2(nlev);
          zaxisInqLbounds(zaxisID1, yb1.data());
          for (int ilev = 0; ilev < nlev; ++ilev) yb2[nlev - ilev - 1] = yb1[ilev];
          zaxisDefLbounds(zaxisID2, yb2.data());

          zaxisInqUbounds(zaxisID1, yb1.data());
          for (int ilev = 0; ilev < nlev; ++ilev) yb2[nlev - ilev - 1] = yb1[ilev];
          zaxisDefUbounds(zaxisID2, yb2.data());
        }

      if (zaxistype == ZAXIS_HYBRID || zaxistype == ZAXIS_HYBRID_HALF)
        {
          const int vctsize = zaxisInqVctSize(zaxisID1);
          if (vctsize && vctsize % 2 == 0)
            {
              Varray<double> vct1(vctsize), vct2(vctsize);
              zaxisInqVct(zaxisID1, vct1.data());
              for (int i = 0; i < vctsize / 2; ++i)
                {
                  vct2[vctsize / 2 - 1 - i] = vct1[i];
                  vct2[vctsize - 1 - i] = vct1[vctsize / 2 + i];
                }
              zaxisDefVct(zaxisID2, vctsize, vct2.data());
            }
        }

      vlistChangeZaxis(vlistID, zaxisID1, zaxisID2);
    }
}

void *
Invertlev(void *process)
{
  size_t nmiss;

  cdo_initialize(process);

  auto dataIsUnchanged = data_is_unchanged();

  cdo_operator_add("invertlev", 0, 0, nullptr);

  operator_check_argc(0);

  auto streamID1 = cdo_open_read(0);

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto vlistID2 = vlistDuplicate(vlistID1);

  auto taxisID1 = vlistInqTaxis(vlistID1);
  auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  invertLevDes(vlistID2);

  auto streamID2 = cdo_open_write(1);

  cdo_def_vlist(streamID2, vlistID2);

  auto gridsizemax = vlistGridsizeMax(vlistID1);
  Varray<double> array(gridsizemax);

  auto nvars = vlistNvars(vlistID1);

  Varray2D<double> vardata(nvars);
  std::vector<std::vector<size_t>> varnmiss(nvars);

  VarList varList1;
  varListInit(varList1, vlistID1);

  auto has3dVar = false;
  for (int varID = 0; varID < nvars; ++varID)
    {
      auto nlevels = varList1[varID].nlevels;
      if (nlevels > 1)
        {
          has3dVar = true;
          vardata[varID].resize(varList1[varID].gridsize * nlevels);
          varnmiss[varID].resize(nlevels);
        }
    }

  if (!has3dVar) cdo_warning("No variables with invertable levels found!");

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

          if (vardata[varID].size())
            {
              auto offset = varList1[varID].gridsize * levelID;
              cdo_read_record(streamID1, &vardata[varID][offset], &nmiss);
              varnmiss[varID][levelID] = nmiss;
            }
          else
            {
              cdo_def_record(streamID2, varID, levelID);
              if (dataIsUnchanged) { cdo_copy_record(streamID2, streamID1); }
              else
                {
                  cdo_read_record(streamID1, array.data(), &nmiss);
                  cdo_write_record(streamID2, array.data(), nmiss);
                }
            }
        }

      for (int varID = 0; varID < nvars; ++varID)
        {
          if (vardata[varID].size())
            {
              auto gridsize = varList1[varID].gridsize;
              auto nlevels = varList1[varID].nlevels;
              for (int levelID = 0; levelID < nlevels; ++levelID)
                {
                  cdo_def_record(streamID2, varID, levelID);

                  auto offset = gridsize * (nlevels - levelID - 1);
                  nmiss = varnmiss[varID][nlevels - levelID - 1];

                  cdo_write_record(streamID2, &vardata[varID][offset], nmiss);
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
