/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

*/

#include <cdi.h>

#include "process_int.h"
#include "printinfo.h"
#include "cdo_zaxis.h"

void *
Pinfo(void *process)
{
  size_t imiss = 0;
  double arrmin, arrmax, arrmean;

  cdo_initialize(process);

  // clang-format off
  auto PINFO  = cdo_operator_add("pinfo",  0, 0, nullptr);
  auto PINFOV = cdo_operator_add("pinfov", 0, 0, nullptr);
  // clang-format on

  (void) (PINFO);  // unused

  auto operatorID = cdo_operator_id();

  operator_check_argc(0);

  auto streamID1 = cdo_open_read(0);

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto vlistID2 = vlistDuplicate(vlistID1);

  auto taxisID1 = vlistInqTaxis(vlistID1);
  auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  auto streamID2 = cdo_open_write(1);
  cdo_def_vlist(streamID2, vlistID2);

  VarList varList1;
  varListInit(varList1, vlistID1);

  auto gridsizemax = vlistGridsizeMax(vlistID1);
  Varray<double> array1(gridsizemax), array2(gridsizemax);

  int indg = 0;
  int tsID = 0;
  while (true)
    {
      auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
      if (nrecs == 0) break;

      auto vDateTime = taxisInqVdatetime(taxisID1);
      auto vdateString = date_to_string(vDateTime.date);
      auto vtimeString = time_to_string(vDateTime.time);

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_def_timestep(streamID2, tsID);

      for (int recID = 0; recID < nrecs; ++recID)
        {
          if (tsID == 0 && recID == 0)
            {
              if (operatorID == PINFOV)
                fprintf(stdout,
                        "   Rec :       Date  Time    Varname     Level    Size    Miss :     Minimum        Mean     Maximum\n");
              else
                fprintf(stdout, "   Rec :       Date  Time    Code  Level    Size    Miss :     Minimum        Mean     Maximum\n");
            }

          int varID, levelID;
          cdo_inq_record(streamID1, &varID, &levelID);
          size_t nmiss;
          cdo_read_record(streamID1, array1.data(), &nmiss);

          const auto &var = varList1[varID];
          indg += 1;
          auto gridsize = var.gridsize;

          fprintf(stdout, "%6d :%s %s ", indg, vdateString.c_str(), vtimeString.c_str());
          if (operatorID == PINFOV)
            fprintf(stdout, "%-8s ", var.name.c_str());
          else
            fprintf(stdout, "%3d", var.code);

          auto level = cdo_zaxis_inq_level(var.zaxisID, levelID);
          fprintf(stdout, " %7g ", level);

          fprintf(stdout, "%7zu %7zu :", gridsize, nmiss);

          if (gridInqType(var.gridID) == GRID_SPECTRAL || (gridsize == 1 && nmiss == 0))
            {
              fprintf(stdout, "            %#12.5g\n", array1[0]);
            }
          else
            {
              if (nmiss)
                {
                  auto mmm = varray_min_max_mean_mv(gridsize, array1, var.missval);
                  arrmin = mmm.min;
                  arrmax = mmm.max;
                  arrmean = mmm.mean;
                  auto ivals = mmm.n;
                  imiss = gridsize - ivals;
                  gridsize = ivals;
                }
              else
                {
                  auto mmm = varray_min_max_mean(gridsize, array1);
                  arrmin = mmm.min;
                  arrmax = mmm.max;
                  arrmean = mmm.mean;
                }

              if (gridsize) { fprintf(stdout, "%#12.5g%#12.5g%#12.5g\n", arrmin, arrmean, arrmax); }
              else { fprintf(stdout, "                     nan\n"); }

              if (imiss != nmiss && nmiss) fprintf(stdout, "Found %zu of %zu missing values!\n", imiss, nmiss);
            }

          varray_copy(gridsize, array1, array2);

          cdo_def_record(streamID2, varID, levelID);
          cdo_write_record(streamID2, array2.data(), nmiss);
        }

      tsID++;
    }

  cdo_stream_close(streamID1);
  cdo_stream_close(streamID2);

  cdo_finish();

  return nullptr;
}
