/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Comp       eq              Equal
      Comp       ne              Not equal
      Comp       le              Less equal
      Comp       lt              Less than
      Comp       ge              Greater equal
      Comp       gt              Greater than
*/

#include <cdi.h>

#include <utility>

#include "cdo_options.h"
#include "process_int.h"
#include "cdo_vlist.h"
#include "cdo_fill.h"
#include "field_functions.h"

static auto func_comp
    = [](auto hasMissvals, auto n, auto mv1, auto mv2, const auto &vIn1, const auto &vIn2, auto &vOut, auto binray_operator) {
        if (hasMissvals)
          {
            if (std::isnan(mv1) || std::isnan(mv2))
              for (size_t i = 0; i < n; ++i)
                vOut[i] = (dbl_is_equal(vIn1[i], mv1) || dbl_is_equal(vIn2[i], mv2)) ? mv1 : binray_operator(vIn1[i], vIn2[i]);
            else
              for (size_t i = 0; i < n; ++i)
                vOut[i] = (is_equal(vIn1[i], mv1) || is_equal(vIn2[i], mv2)) ? mv1 : binray_operator(vIn1[i], vIn2[i]);
          }
        else
          {
            for (size_t i = 0; i < n; ++i) vOut[i] = binray_operator(vIn1[i], vIn2[i]);
          }
      };

void *
Comp(void *process)
{
  enum
  {
    FILL_NONE,
    FILL_TS,
    FILL_REC
  };
  int filltype = FILL_NONE;
  Varray2D<double> vardata;

  cdo_initialize(process);

  cdo_operator_add("eq", FieldFunc_EQ, 0, nullptr);
  cdo_operator_add("ne", FieldFunc_NE, 0, nullptr);
  cdo_operator_add("le", FieldFunc_LE, 0, nullptr);
  cdo_operator_add("lt", FieldFunc_LT, 0, nullptr);
  cdo_operator_add("ge", FieldFunc_GE, 0, nullptr);
  cdo_operator_add("gt", FieldFunc_GT, 0, nullptr);

  auto operatorID = cdo_operator_id();
  auto operFunc = cdo_operator_f1(operatorID);

  operator_check_argc(0);

  auto streamID1 = cdo_open_read(0);
  auto streamID2 = cdo_open_read(1);

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto vlistID2 = cdo_stream_inq_vlist(streamID2);

  auto taxisID1 = vlistInqTaxis(vlistID1);
  auto taxisID2 = vlistInqTaxis(vlistID2);

  auto ntsteps1 = vlistNtsteps(vlistID1);
  auto ntsteps2 = vlistNtsteps(vlistID2);
  if (ntsteps1 == 0) ntsteps1 = 1;
  if (ntsteps2 == 0) ntsteps2 = 1;

  auto fillstream1 = false;

  if (vlistNrecs(vlistID1) != 1 && vlistNrecs(vlistID2) == 1)
    {
      filltype = FILL_REC;
      cdo_print("Filling up stream2 >%s< by copying the first record.", cdo_get_stream_name(1));
      if (ntsteps2 != 1) cdo_abort("stream2 has more than 1 timestep!");
    }
  else if (vlistNrecs(vlistID1) == 1 && vlistNrecs(vlistID2) != 1)
    {
      filltype = FILL_REC;
      cdo_print("Filling up stream1 >%s< by copying the first record.", cdo_get_stream_name(0));
      if (ntsteps1 != 1) cdo_abort("stream1 has more than 1 timestep!");
      fillstream1 = true;
      std::swap(streamID1, streamID2);
      std::swap(vlistID1, vlistID2);
      std::swap(taxisID1, taxisID2);
    }

  if (filltype == FILL_NONE) vlist_compare(vlistID1, vlistID2, CMP_ALL);

  auto gridsizemax = vlistGridsizeMax(vlistID1);
  Varray<double> vaIn1(gridsizemax), vaIn2(gridsizemax), vaOut(gridsizemax);

  double *arrayx1 = vaIn1.data();
  double *arrayx2 = vaIn2.data();

  if (Options::cdoVerbose) cdo_print("Number of timesteps: file1 %d, file2 %d", ntsteps1, ntsteps2);

  if (filltype == FILL_NONE)
    {
      if (ntsteps1 != 1 && ntsteps2 == 1)
        {
          filltype = FILL_TS;
          cdo_print("Filling up stream2 >%s< by copying the first timestep.", cdo_get_stream_name(1));
        }
      else if (ntsteps1 == 1 && ntsteps2 != 1)
        {
          filltype = FILL_TS;
          cdo_print("Filling up stream1 >%s< by copying the first timestep.", cdo_get_stream_name(0));
          fillstream1 = true;
          std::swap(streamID1, streamID2);
          std::swap(vlistID1, vlistID2);
          std::swap(taxisID1, taxisID2);
        }

      if (filltype == FILL_TS) cdo_fill_ts(vlistID2, vardata);
    }

  if (fillstream1)
    {
      arrayx1 = vaIn2.data();
      arrayx2 = vaIn1.data();
    }

  VarList varList1, varList2;
  varListInit(varList1, vlistID1);
  varListInit(varList2, vlistID2);

  auto vlistID3 = vlistDuplicate(vlistID1);

  auto taxisID3 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID3, taxisID3);

  auto streamID3 = cdo_open_write(2);
  cdo_def_vlist(streamID3, vlistID3);

  int tsID = 0;
  while (true)
    {
      auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
      if (nrecs == 0) break;

      if (tsID == 0 || filltype == FILL_NONE)
        {
          auto nrecs2 = cdo_stream_inq_timestep(streamID2, tsID);
          if (nrecs2 == 0) cdo_abort("Input streams have different number of timesteps!");
        }

      cdo_taxis_copy_timestep(taxisID3, taxisID1);
      cdo_def_timestep(streamID3, tsID);

      for (int recID = 0; recID < nrecs; ++recID)
        {
          size_t nmiss1 = 0, nmiss2 = 0;
          int varID, levelID;
          cdo_inq_record(streamID1, &varID, &levelID);
          cdo_read_record(streamID1, arrayx1, &nmiss1);

          if (tsID == 0 || filltype == FILL_NONE)
            {
              if (recID == 0 || filltype != FILL_REC)
                {
                  cdo_inq_record(streamID2, &varID, &levelID);
                  cdo_read_record(streamID2, arrayx2, &nmiss2);
                }

              if (filltype == FILL_TS)
                {
                  auto offset = varList2[varID].gridsize * levelID;
                  array_copy(varList2[varID].gridsize, arrayx2, &vardata[varID][offset]);
                }
            }
          else if (filltype == FILL_TS)
            {
              auto offset = varList2[varID].gridsize * levelID;
              array_copy(varList2[varID].gridsize, &vardata[varID][offset], arrayx2);
            }

          auto datatype1 = varList1[varID].datatype;
          auto gridsize1 = varList1[varID].gridsize;
          auto missval1 = varList1[varID].missval;

          auto xvarID = (filltype == FILL_REC) ? 0 : varID;
          auto datatype2 = varList2[xvarID].datatype;
          auto gridsize2 = varList2[xvarID].gridsize;
          auto missval2 = varList2[xvarID].missval;

          if (gridsize1 != gridsize2)
            cdo_abort("Streams have different gridsize (gridsize1 = %zu; gridsize2 = %zu)!", gridsize1, gridsize2);

          auto ngp = gridsize1;

          if (datatype1 != datatype2)
            {
              if (datatype1 == CDI_DATATYPE_FLT32 && datatype2 == CDI_DATATYPE_FLT64)
                {
                  missval2 = (float) missval2;
                  for (size_t i = 0; i < ngp; ++i) vaIn2[i] = (float) vaIn2[i];
                }
              else if (datatype1 == CDI_DATATYPE_FLT64 && datatype2 == CDI_DATATYPE_FLT32)
                {
                  missval1 = (float) missval1;
                  for (size_t i = 0; i < ngp; ++i) vaIn1[i] = (float) vaIn1[i];
                }
            }

          if (nmiss1 > 0) cdo_check_missval(missval1, varList1[varID].name);
          // if (nmiss2 > 0) cdo_check_missval(missval2, varList2[varID].name);

          auto hasMissvals = (nmiss1 > 0 || nmiss2 > 0);
          // clang-format off
          if      (operFunc == FieldFunc_EQ) func_comp(hasMissvals, ngp, missval1, missval2, vaIn1, vaIn2, vaOut, binary_op_EQ);
          else if (operFunc == FieldFunc_NE) func_comp(hasMissvals, ngp, missval1, missval2, vaIn1, vaIn2, vaOut, binary_op_NE);
          else if (operFunc == FieldFunc_LE) func_comp(hasMissvals, ngp, missval1, missval2, vaIn1, vaIn2, vaOut, binary_op_LE);
          else if (operFunc == FieldFunc_LT) func_comp(hasMissvals, ngp, missval1, missval2, vaIn1, vaIn2, vaOut, binary_op_LT);
          else if (operFunc == FieldFunc_GE) func_comp(hasMissvals, ngp, missval1, missval2, vaIn1, vaIn2, vaOut, binary_op_GE);
          else if (operFunc == FieldFunc_GT) func_comp(hasMissvals, ngp, missval1, missval2, vaIn1, vaIn2, vaOut, binary_op_GT);
          else cdo_abort("Operator not implemented!");
          // clang-format on

          auto nmissOut = varray_num_mv(ngp, vaOut, missval1);
          cdo_def_record(streamID3, varID, levelID);
          cdo_write_record(streamID3, vaOut.data(), nmissOut);
        }

      tsID++;
    }

  cdo_stream_close(streamID3);
  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  cdo_finish();

  return nullptr;
}
