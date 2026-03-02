/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Replace    replace         Replace variables
*/

#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "cdo_vlist.h"
#include "cdo_zaxis.h"

void *
Replace(void *process)
{
  constexpr int MaxVars = 1024;
  int nchvars = 0;
  int idx;
  size_t nmiss = 0;
  int vars1[MaxVars], vars2[MaxVars];
  std::vector<std::vector<int>> varlevel;
  std::vector<std::vector<size_t>> varnmiss2;
  Varray2D<double> vardata2;

  cdo_initialize(process);

  operator_check_argc(0);

  auto streamID1 = cdo_open_read(0);

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto taxisID1 = vlistInqTaxis(vlistID1);
  auto taxisID3 = taxisDuplicate(taxisID1);

  auto streamID2 = cdo_open_read(1);

  auto vlistID2 = cdo_stream_inq_vlist(streamID2);

  VarList varList1, varList2;
  varListInit(varList1, vlistID1);
  varListInit(varList2, vlistID2);

  // compare all variables in vlistID2

  auto nvars1 = vlistNvars(vlistID1);
  auto nvars2 = vlistNvars(vlistID2);

  for (int varID2 = 0; varID2 < nvars2; ++varID2)
    {
      int varID1;
      for (varID1 = 0; varID1 < nvars1; ++varID1)
        {
          if (varList1[varID1].name == varList2[varID2].name) break;
        }

      if (varID1 < nvars1)
        {
          auto gridsize1 = varList1[varID1].gridsize;
          auto nlevel1 = varList1[varID1].nlevels;

          auto gridsize2 = varList2[varID2].gridsize;
          auto nlevel2 = varList2[varID2].nlevels;

          if (gridsize1 != gridsize2) cdo_abort("Variables have different gridsize!");

          if (nlevel1 < nlevel2) cdo_abort("Variables have different number of levels!");

          if (Options::cdoVerbose) cdo_print("Variable %s replaced.", varList1[varID1].name);

          vars1[nchvars] = varID1;
          vars2[nchvars] = varID2;
          nchvars++;
          if (nchvars > MaxVars) cdo_abort("Internal problem - too many variables!");
        }
      else { cdo_warning("Variable %s not found!", varList2[varID2].name); }
    }

  if (nchvars)
    {
      vardata2.resize(nchvars);
      varnmiss2.resize(nchvars);
      varlevel.resize(nchvars);
      for (idx = 0; idx < nchvars; idx++)
        {
          auto varID1 = vars1[idx];
          auto varID2 = vars2[idx];
          auto nlevel1 = varList1[varID1].nlevels;
          auto nlevel2 = varList2[varID2].nlevels;
          auto gridsize = varList2[varID2].gridsize;
          vardata2[idx].resize(nlevel2 * gridsize);
          varnmiss2[idx].resize(nlevel2);
          varlevel[idx].resize(nlevel1);
          /*
          for ( levelID = 0; levelID < nlevel1; levelID++ )
            varlevel[idx][levelID] = levelID;
          */
          if (nlevel2 <= nlevel1)
            {
              Varray<double> level1(nlevel1), level2(nlevel2);
              cdo_zaxis_inq_levels(vlistInqVarZaxis(vlistID1, varID1), level1.data());
              cdo_zaxis_inq_levels(vlistInqVarZaxis(vlistID2, varID2), level2.data());

              for (int levelID = 0; levelID < nlevel1; ++levelID) varlevel[idx][levelID] = -1;

              for (int l2 = 0; l2 < nlevel2; ++l2)
                {
                  int l1;
                  for (l1 = 0; l1 < nlevel1; ++l1)
                    if (IS_EQUAL(level2[l2], level1[l1]))
                      {
                        varlevel[idx][l1] = l2;
                        break;
                      }

                  if (l1 == nlevel1) cdo_warning("Variable %s on level %g not found!", varList2[varID2].name, level2[l2]);
                }
            }
        }
    }

  auto vlistID3 = vlistDuplicate(vlistID1);

  auto streamID3 = cdo_open_write(2);

  vlistDefTaxis(vlistID3, taxisID3);
  cdo_def_vlist(streamID3, vlistID3);

  auto gridsizemax = vlistGridsizeMax(vlistID1);
  Varray<double> array(gridsizemax);

  auto nts2 = vlistNtsteps(vlistID2);

  int tsID = 0;
  while (true)
    {
      auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
      if (nrecs == 0) break;

      cdo_taxis_copy_timestep(taxisID3, taxisID1);

      if (tsID == 0 || (nts2 != 0 && nts2 != 1))
        {
          auto nrecs2 = cdo_stream_inq_timestep(streamID2, tsID);
          if (nrecs2 == 0) cdo_abort("Input streams have different number of timesteps!");

          for (int recID = 0; recID < nrecs2; ++recID)
            {
              int varID, levelID;
              cdo_inq_record(streamID2, &varID, &levelID);

              for (idx = 0; idx < nchvars; idx++)
                if (vars2[idx] == varID)
                  {
                    auto offset = varList2[varID].gridsize * levelID;
                    cdo_read_record(streamID2, &vardata2[idx][offset], &nmiss);
                    varnmiss2[idx][levelID] = nmiss;
                    break;
                  }
            }
        }

      cdo_def_timestep(streamID3, tsID);

      for (int recID = 0; recID < nrecs; ++recID)
        {
          int varID, levelID;
          cdo_inq_record(streamID1, &varID, &levelID);

          auto parray = array.data();

          for (idx = 0; idx < nchvars; idx++)
            if (vars1[idx] == varID)
              {
                auto levelID2 = varlevel[idx][levelID];
                if (levelID2 != -1)
                  {
                    auto offset = varList1[varID].gridsize * levelID2;
                    parray = &vardata2[idx][offset];
                    nmiss = varnmiss2[idx][levelID2];
                    break;
                  }
              }

          if (idx == nchvars) cdo_read_record(streamID1, parray, &nmiss);

          cdo_def_record(streamID3, varID, levelID);
          cdo_write_record(streamID3, parray, nmiss);
        }

      tsID++;
    }

  cdo_stream_close(streamID3);
  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  cdo_finish();

  return nullptr;
}
