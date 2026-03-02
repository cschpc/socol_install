/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Michal Koutek, KMNI
          Uwe Schulzweida

*/

/*
   This module "SampleGrid" contains the following operators:

    samplegrid      Resample current grid with given factor, typically 2 (which will half the resolution);
                    tested on curvilinear and LCC grids;
    subgrid         Similar to selindexbox but this operator works for LCC grids (tested on HARMONIE NWP model).
*/

#include <cdi.h>

#include "process_int.h"
#include "param_conversion.h"
#include <mpim_grid.h>
#include "grid_define.h"

static void
sampleData(const double *array1, int gridID1, double *array2, int gridID2, int resampleFactor)
{
  const auto nlon1 = gridInqXsize(gridID1);
  const auto nlat1 = gridInqYsize(gridID1);

  const auto nlon2 = gridInqXsize(gridID2);
  const auto nlat2 = gridInqYsize(gridID2);

  if (cdoDebugExt >= 100)
    cdo_print("%s(): (nlon1: %zu; nlat1: %zu) => (nlon2: %zu; nlat2: %zu); "
              "gridID1: %d; gridID2: %d; resampleFactor: %d)",
              __func__, nlon1, nlat1, nlon2, nlat2, gridID1, gridID2, resampleFactor);

  for (size_t ilat1 = 0; ilat1 < nlat1; ilat1 += resampleFactor)
    for (size_t ilon1 = 0; ilon1 < nlon1; ilon1 += resampleFactor) *array2++ = array1[ilat1 * nlon1 + ilon1];
}

static void
cropData(double *array1, int gridID1, double *array2, int gridID2, int subI0, int subI1, int subJ0, int subJ1)
{
  const long nlon1 = gridInqXsize(gridID1);
  const long nlon2 = gridInqXsize(gridID2);
  const long rowLen = subI1 - subI0 + 1;  // must be same as nlon1

  if (rowLen != nlon2) cdo_abort("cropData() rowLen!= nlon2 [%d != %d]", rowLen, nlon2);

  if (cdoDebugExt >= 10) cdo_print("cropData(%d,%d,%d,%d) ...", subI0, subI1, subJ0, subJ1);

  long array2Idx = 0;
  for (long ilat1 = subJ0; ilat1 <= subJ1; ilat1++)  // copy the last row as well..
    {
      array_copy(rowLen, &array1[ilat1 * nlon1 + subI0], &array2[array2Idx]);
      array2Idx += rowLen;
    }
}

void *
Samplegrid(void *process)
{
  int resampleFactor = 0;
  int subI0 = 0, subI1 = 0, subJ0 = 0, subJ1 = 0;
  struct sbox_t
  {
    int gridSrcID, gridIDsampled;
  };

  cdo_initialize(process);

  const auto SAMPLEGRID = cdo_operator_add("samplegrid", 0, 0, "resample factor, typically 2 (which will half the resolution)");
  const auto SUBGRID = cdo_operator_add("subgrid", 0, 0, " sub-grid indices: i0,i1,j0,j1");

  const auto operatorID = cdo_operator_id();

  const auto nch = cdo_operator_argc();

  if (operatorID == SAMPLEGRID)
    {
      Debug(cdoDebugExt, "samplegrid operator requested..");
      if (nch < 1) cdo_abort("Number of input arguments < 1; At least 1 argument needed: resample-factor (2,3,4, .. etc)");
      resampleFactor = parameter_to_int(cdo_operator_argv(0));

      Debug(cdoDebugExt, "resampleFactor = %d", resampleFactor);
    }
  else if (operatorID == SUBGRID)
    {
      Debug(cdoDebugExt, "subgrid operator requested..");
      if (nch < 4)
        cdo_abort("Number of input arguments < 4; Must specify sub-grid indices: i0,i1,j0,j1; This works only with LCC grid."
                  " For other grids use: selindexbox");
      subI0 = parameter_to_int(cdo_operator_argv(0));
      subI1 = parameter_to_int(cdo_operator_argv(1));
      subJ0 = parameter_to_int(cdo_operator_argv(2));
      subJ1 = parameter_to_int(cdo_operator_argv(3));
    }
  else
    cdo_abort("Unknown operator ...");

  const auto streamID1 = cdo_open_read(0);

  const auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  const auto vlistID2 = vlistDuplicate(vlistID1);

  const auto taxisID1 = vlistInqTaxis(vlistID1);
  const auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  const auto nvars = vlistNvars(vlistID1);
  std::vector<bool> vars(nvars, false);

  auto ngrids = vlistNgrids(vlistID1);

  Debug(cdoDebugExt, "ngrids=%d", ngrids);

  std::vector<sbox_t> sbox(ngrids);

  for (int index = 0; index < ngrids; ++index)
    {
      const auto gridSrcID = vlistGrid(vlistID1, index);
      int gridIDsampled = -1;

      if (gridInqSize(gridSrcID) <= 1) continue;

      int gridtype = gridInqType(gridSrcID);
      if (!(gridtype == GRID_GAUSSIAN || gridtype == GRID_LONLAT || gridtype == GRID_PROJECTION || gridtype == GRID_CURVILINEAR
            || gridtype == GRID_GENERIC))
        cdo_abort("Unsupported gridtype: %s", gridNamePtr(gridtype));

      if (operatorID == SAMPLEGRID) { gridIDsampled = cdo_define_sample_grid(gridSrcID, resampleFactor); }
      else if (operatorID == SUBGRID) { gridIDsampled = cdo_define_subgrid_grid(gridSrcID, subI0, subI1, subJ0, subJ1); }

      sbox[index].gridSrcID = gridSrcID;
      sbox[index].gridIDsampled = gridIDsampled;

      // if ( cdoDebugExt>=10 ) cdo_print_griddes(gridSrcID, 1);
      // if ( cdoDebugExt>=10 ) cdo_print_griddes(gridIDsampled, 1);

      vlistChangeGridIndex(vlistID2, index, gridIDsampled);

      for (int varID = 0; varID < nvars; ++varID)
        if (gridSrcID == vlistInqVarGrid(vlistID1, varID)) vars[varID] = true;
    }

  Debug(cdoDebugExt, [&]() {
    if (operatorID == SAMPLEGRID) Debug("Resampled grid has been created.");
    if (operatorID == SUBGRID) Debug("Sub-grid has been created.");
  });

  const auto streamID2 = cdo_open_write(1);
  cdo_def_vlist(streamID2, vlistID2);

  auto gridsize = vlistGridsizeMax(vlistID1);
  if (vlistNumber(vlistID1) != CDI_REAL) gridsize *= 2;
  Varray<double> array1(gridsize);

  size_t gridsize2 = vlistGridsizeMax(vlistID2);
  if (vlistNumber(vlistID2) != CDI_REAL) gridsize2 *= 2;
  Varray<double> array2(gridsize2);

  Debug(cdoDebugExt, "gridsize = %ld, gridsize2 = %ld", gridsize, gridsize2);

  int tsID = 0;
  while (true)
    {
      const auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
      if (nrecs == 0) break;

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_def_timestep(streamID2, tsID);

      for (int recID = 0; recID < nrecs; ++recID)
        {
          int varID, levelID;
          cdo_inq_record(streamID1, &varID, &levelID);
          size_t nmiss;
          cdo_read_record(streamID1, array1.data(), &nmiss);

          cdo_def_record(streamID2, varID, levelID);

          Debug(cdoDebugExt >= 20, "Processing record (%d) of %d.", recID, nrecs);

          if (vars[varID])
            {
              const auto gridSrcID = vlistInqVarGrid(vlistID1, varID);

              int index;
              for (index = 0; index < ngrids; ++index)
                if (gridSrcID == sbox[index].gridSrcID) break;

              if (index == ngrids) cdo_abort("Internal problem, grid not found!");

              const int gridIDsampled = sbox[index].gridIDsampled;
              gridsize2 = gridInqSize(gridIDsampled);

              if (operatorID == SAMPLEGRID) { sampleData(array1.data(), gridSrcID, array2.data(), gridIDsampled, resampleFactor); }
              else if (operatorID == SUBGRID)
                {
                  cropData(array1.data(), gridSrcID, array2.data(), gridIDsampled, subI0, subI1, subJ0, subJ1);
                }

              if (nmiss)
                {
                  const auto missval = vlistInqVarMissval(vlistID2, varID);
                  nmiss = varray_num_mv(gridsize2, array2, missval);
                }

              cdo_write_record(streamID2, array2.data(), nmiss);
            }
          else { cdo_write_record(streamID2, array1.data(), nmiss); }
        }

      tsID++;
    }

  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  vlistDestroy(vlistID2);

  cdo_finish();

  return nullptr;
}
