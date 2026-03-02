/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Cond2      ifthenelse      If then else
*/

#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "cdo_vlist.h"
#include "cdo_fill.h"

void *
Cond2(void *process)
{
  enum
  {
    FILL_NONE,
    FILL_TS,
    FILL_REC
  };
  int filltype = FILL_NONE;
  double missval1 = -9.E33;
  size_t nmiss1 = 0;
  Varray2D<size_t> varnmiss1;
  Varray2D<double> vardata1;

  cdo_initialize(process);

  cdo_operator_add("ifthenelse", 0, 0, nullptr);

  operator_check_argc(0);

  auto streamID1 = cdo_open_read(0);
  auto streamID2 = cdo_open_read(1);
  auto streamID3 = cdo_open_read(2);

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto vlistID2 = cdo_stream_inq_vlist(streamID2);
  auto vlistID3 = cdo_stream_inq_vlist(streamID3);
  auto vlistID4 = vlistDuplicate(vlistID2);

  auto taxisID2 = vlistInqTaxis(vlistID2);
  auto taxisID4 = taxisDuplicate(taxisID2);
  vlistDefTaxis(vlistID4, taxisID4);

  auto ntsteps1 = vlistNtsteps(vlistID1);
  auto ntsteps2 = vlistNtsteps(vlistID2);
  if (ntsteps1 == 0) ntsteps1 = 1;
  if (ntsteps2 == 0) ntsteps2 = 1;

  if (vlistNrecs(vlistID1) == 1 && vlistNrecs(vlistID2) != 1)
    {
      filltype = FILL_REC;
      cdo_print("Filling up stream1 >%s< by copying the first record.", cdo_get_stream_name(0));
    }

  if (filltype == FILL_NONE) vlist_compare(vlistID1, vlistID2, CMP_DIM);

  vlist_compare(vlistID2, vlistID3, CMP_DIM);

  auto streamID4 = cdo_open_write(3);
  cdo_def_vlist(streamID4, vlistID4);

  VarList varList1, varList2;
  varListInit(varList1, vlistID1);
  varListInit(varList2, vlistID2);

  auto gridsizemax = vlistGridsizeMax(vlistID1);

  if (filltype == FILL_REC && gridsizemax != gridInqSize(vlistGrid(vlistID1, 0)))
    cdo_abort("Stream1 >%s< has wrong gridsize!", cdo_get_stream_name(0));

  Varray<double> array1(gridsizemax), array2(gridsizemax), array3(gridsizemax), array4(gridsizemax);

  if (Options::cdoVerbose)
    cdo_print("Number of timesteps: file1 %d, file2 %d, file3 %d", ntsteps1, ntsteps2, vlistNtsteps(vlistID3));

  if (filltype == FILL_NONE)
    {
      if (ntsteps1 == 1 && ntsteps2 != 1)
        {
          filltype = FILL_TS;
          cdo_print("Filling up stream1 >%s< by copying the first timestep.", cdo_get_stream_name(0));

          cdo_fill_ts(vlistID1, vardata1, varnmiss1);
        }
    }

  int tsID = 0;
  while (true)
    {
      auto nrecs = cdo_stream_inq_timestep(streamID2, tsID);
      if (nrecs == 0) break;

      auto nrecs3 = cdo_stream_inq_timestep(streamID3, tsID);
      if (nrecs3 == 0) cdo_abort("Input streams have different number of timesteps!");

      if (tsID == 0 || filltype == FILL_NONE)
        {
          auto nrecs2 = cdo_stream_inq_timestep(streamID1, tsID);
          if (nrecs2 == 0) cdo_abort("Input streams have different number of timesteps!");
        }

      cdo_taxis_copy_timestep(taxisID4, taxisID2);
      cdo_def_timestep(streamID4, tsID);

      for (int recID = 0; recID < nrecs; ++recID)
        {
          int varID, levelID;
          cdo_inq_record(streamID2, &varID, &levelID);
          size_t nmiss;
          cdo_read_record(streamID2, &array2[0], &nmiss);

          cdo_inq_record(streamID3, &varID, &levelID);
          cdo_read_record(streamID3, &array3[0], &nmiss);

          if (tsID == 0 || filltype == FILL_NONE)
            {
              if (recID == 0 || filltype != FILL_REC)
                {
                  cdo_inq_record(streamID1, &varID, &levelID);
                  cdo_read_record(streamID1, &array1[0], &nmiss1);
                }

              if (filltype == FILL_TS)
                {
                  auto gridsize = varList1[varID].gridsize;
                  auto offset = gridsize * levelID;
                  array_copy(gridsize, &array1[0], &vardata1[varID][offset]);
                  varnmiss1[varID][levelID] = nmiss1;
                }
            }
          else if (filltype == FILL_TS)
            {
              auto gridsize = varList1[varID].gridsize;
              auto offset = gridsize * levelID;
              array_copy(gridsize, &vardata1[varID][offset], &array1[0]);
              nmiss1 = varnmiss1[varID][levelID];
            }

          auto gridsize = varList2[varID].gridsize;
          auto missval2 = varList2[varID].missval;
          if (recID == 0 || filltype != FILL_REC) missval1 = varList1[varID].missval;

          if (nmiss1 > 0) cdo_check_missval(missval1, varList1[varID].name);

          for (size_t i = 0; i < gridsize; ++i)
            array4[i] = DBL_IS_EQUAL(array1[i], missval1) ? missval2 : !DBL_IS_EQUAL(array1[i], 0.) ? array2[i] : array3[i];

          nmiss = varray_num_mv(gridsize, array4, missval2);
          cdo_def_record(streamID4, varID, levelID);
          cdo_write_record(streamID4, array4.data(), nmiss);
        }

      tsID++;
    }

  cdo_stream_close(streamID4);
  cdo_stream_close(streamID3);
  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  cdo_finish();

  return nullptr;
}
