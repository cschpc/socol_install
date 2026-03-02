/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Cond       ifthen          If then
      Cond       ifnotthen       If not then
*/

#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "cdo_vlist.h"
#include "cdo_fill.h"

static void
operator_IFTHEN(size_t n, double mv1, double mv2, const Varray<double> &vIn1, const Varray<double> &vIn2, Varray<double> &vOut)
{
  if (std::isnan(mv1))
    for (size_t i = 0; i < n; ++i) vOut[i] = (!dbl_is_equal(vIn1[i], mv1) && !dbl_is_equal(vIn1[i], 0.0)) ? vIn2[i] : mv2;
  else
    for (size_t i = 0; i < n; ++i) vOut[i] = (!is_equal(vIn1[i], mv1) && !is_equal(vIn1[i], 0.0)) ? vIn2[i] : mv2;
}

static void
operator_IFNOTTHEN(size_t n, double mv1, double mv2, const Varray<double> &vIn1, const Varray<double> &vIn2, Varray<double> &vOut)
{
  if (std::isnan(mv1))
    for (size_t i = 0; i < n; ++i) vOut[i] = (!dbl_is_equal(vIn1[i], mv1) && dbl_is_equal(vIn1[i], 0.0)) ? vIn2[i] : mv2;
  else
    for (size_t i = 0; i < n; ++i) vOut[i] = (!is_equal(vIn1[i], mv1) && is_equal(vIn1[i], 0.0)) ? vIn2[i] : mv2;
}

void *
Cond(void *process)
{
  enum
  {
    FILL_NONE,
    FILL_TS,
    FILL_REC
  };
  int filltype = FILL_NONE;
  double missval1 = -9.E33;
  Varray2D<size_t> varnmiss1;
  Varray2D<double> vardata1;

  cdo_initialize(process);

  // clang-format off
  auto IFTHEN    = cdo_operator_add("ifthen",    0, 0, nullptr);
  auto IFNOTTHEN = cdo_operator_add("ifnotthen", 0, 0, nullptr);
  // clang-format on

  auto operatorID = cdo_operator_id();

  operator_check_argc(0);

  auto streamID1 = cdo_open_read(0);
  auto streamID2 = cdo_open_read(1);

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto vlistID2 = cdo_stream_inq_vlist(streamID2);
  auto vlistID3 = vlistDuplicate(vlistID2);

  auto taxisID2 = vlistInqTaxis(vlistID2);
  auto taxisID3 = taxisDuplicate(taxisID2);
  vlistDefTaxis(vlistID3, taxisID3);

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

  VarList varList1, varList2;
  varListInit(varList1, vlistID1);
  varListInit(varList2, vlistID2);

  auto streamID3 = cdo_open_write(2);
  cdo_def_vlist(streamID3, vlistID3);

  auto gridsizemax = vlistGridsizeMax(vlistID2);

  if (filltype == FILL_REC && gridsizemax != gridInqSize(vlistGrid(vlistID1, 0)))
    cdo_abort("Stream1 >%s< has wrong gridsize!", cdo_get_stream_name(0));

  Varray<double> vaIn1(gridsizemax), vaIn2(gridsizemax), vaOut(gridsizemax);

  if (Options::cdoVerbose) cdo_print("Number of timesteps: file1 %d, file2 %d", ntsteps1, ntsteps2);

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

      if (tsID == 0 || filltype == FILL_NONE)
        {
          auto nrecs2 = cdo_stream_inq_timestep(streamID1, tsID);
          if (nrecs2 == 0) cdo_abort("Input streams have different number of timesteps!");
        }

      cdo_taxis_copy_timestep(taxisID3, taxisID2);
      cdo_def_timestep(streamID3, tsID);

      for (int recID = 0; recID < nrecs; ++recID)
        {
          int varID, levelID;
          cdo_inq_record(streamID2, &varID, &levelID);
          size_t nmiss1 = 0, nmiss2;
          cdo_read_record(streamID2, &vaIn2[0], &nmiss2);

          if (tsID == 0 || filltype == FILL_NONE)
            {
              if (recID == 0 || filltype != FILL_REC)
                {
                  cdo_inq_record(streamID1, &varID, &levelID);
                  cdo_read_record(streamID1, &vaIn1[0], &nmiss1);
                }

              if (filltype == FILL_TS)
                {
                  auto offset = varList1[varID].gridsize * levelID;
                  array_copy(varList1[varID].gridsize, &vaIn1[0], &vardata1[varID][offset]);
                  varnmiss1[varID][levelID] = nmiss1;
                }
            }
          else if (filltype == FILL_TS)
            {
              auto offset = varList1[varID].gridsize * levelID;
              array_copy(varList1[varID].gridsize, &vardata1[varID][offset], &vaIn1[0]);
              nmiss1 = varnmiss1[varID][levelID];
            }

          auto ngp = varList2[varID].gridsize;
          auto missval2 = varList2[varID].missval;
          if (recID == 0 || filltype != FILL_REC) missval1 = varList1[varID].missval;

          if (nmiss1 > 0) cdo_check_missval(missval1, varList1[varID].name);

          // clang-format off
          if      (operatorID == IFTHEN)    operator_IFTHEN(ngp, missval1, missval2, vaIn1, vaIn2, vaOut);
          else if (operatorID == IFNOTTHEN) operator_IFNOTTHEN(ngp, missval1, missval2, vaIn1, vaIn2, vaOut);
          else cdo_abort("Operator not implemented!");
          // clang-format on

          auto nmiss3 = varray_num_mv(ngp, vaOut, missval2);
          cdo_def_record(streamID3, varID, levelID);
          cdo_write_record(streamID3, vaOut.data(), nmiss3);
        }

      tsID++;
    }

  cdo_stream_close(streamID3);
  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  cdo_finish();

  return nullptr;
}
