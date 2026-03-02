/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Enlarge    enlarge         Enlarge fields
*/

#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "griddes.h"

void *
Enlarge(void *process)
{
  bool linfo = true;

  cdo_initialize(process);

  operator_check_argc(1);

  const auto gridID2 = cdo_define_grid(cdo_operator_argv(0));
  const auto xsize2 = gridInqXsize(gridID2);
  const auto ysize2 = gridInqYsize(gridID2);

  if (Options::cdoVerbose) fprintf(stderr, "gridID2 %d, xsize2 %zu, ysize2 %zu\n", gridID2, xsize2, ysize2);

  const auto streamID1 = cdo_open_read(0);

  const auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  const auto vlistID2 = vlistDuplicate(vlistID1);

  const auto taxisID1 = vlistInqTaxis(vlistID1);
  const auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  const auto gridsize2 = gridInqSize(gridID2);
  if (gridsize2 < vlistGridsizeMax(vlistID1)) cdo_abort("Gridsize of input stream is greater than new gridsize!");

  const auto ngrids = vlistNgrids(vlistID1);
  for (int index = 0; index < ngrids; ++index) vlistChangeGridIndex(vlistID2, index, gridID2);

  const auto streamID2 = cdo_open_write(1);

  cdo_def_vlist(streamID2, vlistID2);

  Varray<double> array1(gridsize2), array2(gridsize2);

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

          const auto missval = vlistInqVarMissval(vlistID1, varID);
          const auto gridID1 = vlistInqVarGrid(vlistID1, varID);
          const auto gridsize1 = gridInqSize(gridID1);

          auto xsize1 = gridInqXsize(gridID1);
          auto ysize1 = gridInqYsize(gridID1);
          if (xsize1 == 0) xsize1 = 1;
          if (ysize1 == 0) ysize1 = 1;

          if (xsize1 == 1 && ysize1 == ysize2 && xsize1 * ysize1 == gridsize1)
            {
              if (linfo)
                {
                  cdo_print("Enlarge zonal");
                  linfo = false;
                }

              for (size_t iy = 0; iy < ysize2; iy++)
                for (size_t ix = 0; ix < xsize2; ix++) array2[ix + iy * xsize2] = array1[iy];

              if (nmiss) nmiss *= xsize2;
            }
          else if (ysize1 == 1 && xsize1 == xsize2 && xsize1 * ysize1 == gridsize1)
            {
              if (linfo)
                {
                  cdo_print("Enlarge meridional");
                  linfo = false;
                }

              for (size_t iy = 0; iy < ysize2; iy++)
                for (size_t ix = 0; ix < xsize2; ix++) array2[ix + iy * xsize2] = array1[ix];

              if (nmiss) nmiss *= ysize2;
            }
          else
            {
              varray_copy(gridsize1, array1, array2);
              for (size_t i = gridsize1; i < gridsize2; ++i) array2[i] = array1[gridsize1 - 1];

              if (nmiss && DBL_IS_EQUAL(array1[gridsize1 - 1], missval)) nmiss += (gridsize2 - gridsize1);
            }

          cdo_def_record(streamID2, varID, levelID);
          cdo_write_record(streamID2, array2.data(), nmiss);
        }

      tsID++;
    }

  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  cdo_finish();

  return nullptr;
}
