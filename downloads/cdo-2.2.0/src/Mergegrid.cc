/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

*/

#include <cdi.h>

#include "process_int.h"
#include "cdo_vlist.h"
#include <mpim_grid.h>

void genGridIndex(int gridID1, int gridID2, std::vector<long> &index);

void *
Mergegrid(void *process)
{
  cdo_initialize(process);

  operator_check_argc(0);

  auto streamID1 = cdo_open_read(0);

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto taxisID1 = vlistInqTaxis(vlistID1);
  auto taxisID3 = taxisDuplicate(taxisID1);

  auto streamID2 = cdo_open_read(1);
  auto vlistID2 = cdo_stream_inq_vlist(streamID2);

  vlist_compare(vlistID1, vlistID2, CMP_NAME | CMP_NLEVEL);

  int ndiffgrids = 0;
  for (int index = 1; index < vlistNgrids(vlistID1); ++index)
    if (vlistGrid(vlistID1, 0) != vlistGrid(vlistID1, index)) ndiffgrids++;

  if (ndiffgrids > 0) cdo_abort("Too many different grids in %s!", cdo_get_stream_name(0));

  ndiffgrids = 0;
  for (int index = 1; index < vlistNgrids(vlistID2); ++index)
    if (vlistGrid(vlistID2, 0) != vlistGrid(vlistID2, index)) ndiffgrids++;

  if (ndiffgrids > 0) cdo_abort("Too many different grids in %s!", cdo_get_stream_name(1));

  auto gridID1 = vlistGrid(vlistID1, 0);
  auto gridID2 = vlistGrid(vlistID2, 0);

  auto gridsize1 = gridInqSize(gridID1);
  auto gridsize2 = gridInqSize(gridID2);

  Varray<double> array1(gridsize1), array2(gridsize2);
  std::vector<long> gindex(gridsize2);

  genGridIndex(gridID1, gridID2, gindex);

  auto vlistID3 = vlistDuplicate(vlistID1);
  auto streamID3 = cdo_open_write(2);

  vlistDefTaxis(vlistID3, taxisID3);
  cdo_def_vlist(streamID3, vlistID3);

  int tsID = 0;
  while (true)
    {
      auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
      if (nrecs == 0) break;

      cdo_taxis_copy_timestep(taxisID3, taxisID1);

      auto nrecs2 = cdo_stream_inq_timestep(streamID2, tsID);
      if (nrecs2 == 0) cdo_abort("Input streams have different number of timesteps!");

      if (nrecs != nrecs2) cdo_abort("Input streams have different number of records!");

      cdo_def_timestep(streamID3, tsID);

      for (int recID = 0; recID < nrecs; ++recID)
        {
          int varID, levelID;
          cdo_inq_record(streamID2, &varID, &levelID);
          size_t nmiss2;
          cdo_read_record(streamID2, array2.data(), &nmiss2);

          auto missval2 = vlistInqVarMissval(vlistID2, varID);

          cdo_inq_record(streamID1, &varID, &levelID);
          size_t nmiss1;
          cdo_read_record(streamID1, array1.data(), &nmiss1);

          auto missval1 = vlistInqVarMissval(vlistID1, varID);

          for (size_t i = 0; i < gridsize2; ++i)
            {
              if (gindex[i] >= 0 && !DBL_IS_EQUAL(array2[i], missval2)) { array1[gindex[i]] = array2[i]; }
            }

          if (nmiss1)
            {
              nmiss1 = 0;
              for (size_t i = 0; i < gridsize1; ++i)
                if (DBL_IS_EQUAL(array1[i], missval1)) nmiss1++;
            }

          cdo_def_record(streamID3, varID, levelID);
          cdo_write_record(streamID3, array1.data(), nmiss1);
        }

      tsID++;
    }

  cdo_stream_close(streamID3);
  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  cdo_finish();

  return nullptr;
}
