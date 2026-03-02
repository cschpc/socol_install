/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Arith      add             Add two fields
      Arith      sub             Subtract two fields
      Arith      mul             Multiply two fields
      Arith      div             Divide two fields
      Arith      min             Minimum of two fields
      Arith      max             Maximum of two fields
      Arith      atan2           Arc tangent of two fields
*/

#include <cdi.h>

#include <utility>

#include "cdo_options.h"
#include "process_int.h"
#include "cdo_vlist.h"
#include "cdo_options.h"
#include "cdo_fill.h"
#include "field_functions.h"

static void
addOperators(void)
{
  // clang-format off
  cdo_operator_add("add",     FieldFunc_Add,     1, nullptr);
  cdo_operator_add("sub",     FieldFunc_Sub,     1, nullptr);
  cdo_operator_add("mul",     FieldFunc_Mul,     1, nullptr);
  cdo_operator_add("div",     FieldFunc_Div,     1, nullptr);
  cdo_operator_add("min",     FieldFunc_Min,     0, nullptr);
  cdo_operator_add("max",     FieldFunc_Max,     0, nullptr);
  cdo_operator_add("atan2",   FieldFunc_Atan2,   0, nullptr);
  cdo_operator_add("setmiss", FieldFunc_Setmiss, 0, nullptr);
  // clang-format on
}

void *
Arith(void *process)
{
  enum class FillType
  {
    NONE,
    TS,
    VAR,
    VARTS,
    FILE
  };
  auto fillType{ FillType::NONE };
  int nlevels2 = 1;
  int levelID2 = -1;
  Varray2D<size_t> varnmiss;
  Varray2D<double> vardata;
  std::vector<size_t> varnmiss2;
  Varray<double> vardata2;

  cdo_initialize(process);

  addOperators();

  auto operatorID = cdo_operator_id();
  auto operfunc = cdo_operator_f1(operatorID);
  bool opercplx = cdo_operator_f2(operatorID);

  operator_check_argc(0);

  auto streamID1 = cdo_open_read(0);
  auto streamID2 = cdo_open_read(1);

  auto streamID2x = streamID2;

  Field field1, field2;
  auto fieldx1 = &field1;
  auto fieldx2 = &field2;

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto vlistID2 = cdo_stream_inq_vlist(streamID2);
  auto vlistID1x = vlistID1;
  auto vlistID2x = vlistID2;

  if (Options::cdoVerbose) vlistPrint(vlistID1);
  if (Options::cdoVerbose) vlistPrint(vlistID2);

  auto taxisID1 = vlistInqTaxis(vlistID1);
  auto taxisID2 = vlistInqTaxis(vlistID2);

  auto ntsteps1 = vlistNtsteps(vlistID1);
  auto ntsteps2 = vlistNtsteps(vlistID2);
  if (ntsteps1 == 0) ntsteps1 = 1;
  if (ntsteps2 == 0) ntsteps2 = 1;

  auto nvars1 = vlistNvars(vlistID1);
  auto nvars2 = vlistNvars(vlistID2);

  auto fillStream1 = false, fillStream2 = false, fillStream1x = false;
  if (nvars1 == 1 && nvars2 == 1)
    {
      fillStream2 = (vlistNrecs(vlistID1) != 1 && vlistNrecs(vlistID2) == 1);
      fillStream1 = (vlistNrecs(vlistID1) == 1 && vlistNrecs(vlistID2) != 1);
      if (fillStream1 && ntsteps1 != 1 && ntsteps2 == 1)
        {
          fillStream1 = false;
          fillStream2 = true;
          fillStream1x = true;
        }
    }
  else
    {
      fillStream2 = (nvars1 != 1 && nvars2 == 1);
      fillStream1 = (nvars1 == 1 && nvars2 != 1);
    }

  if (fillStream1x)
    {
      nlevels2 = vlist_compare_x(vlistID2, vlistID1, CMP_DIM);

      fillType = FillType::NONE;
      cdo_print("Filling up stream1 >%s< by copying the first variable of each timestep.", cdo_get_stream_name(0));
    }
  else if (fillStream2)
    {
      nlevels2 = vlist_compare_x(vlistID1, vlistID2, CMP_DIM);

      if (ntsteps1 != 1 && ntsteps2 == 1)
        {
          fillType = FillType::VAR;
          cdo_print("Filling up stream2 >%s< by copying the first variable.", cdo_get_stream_name(1));
        }
      else
        {
          fillType = FillType::VARTS;
          cdo_print("Filling up stream2 >%s< by copying the first variable of each timestep.", cdo_get_stream_name(1));
        }
    }
  else if (fillStream1)
    {
      nlevels2 = vlist_compare_x(vlistID2, vlistID1, CMP_DIM);

      if (ntsteps1 == 1 && ntsteps2 != 1)
        {
          fillType = FillType::VAR;
          cdo_print("Filling up stream1 >%s< by copying the first variable.", cdo_get_stream_name(0));
        }
      else
        {
          fillType = FillType::VARTS;
          cdo_print("Filling up stream1 >%s< by copying the first variable of each timestep.", cdo_get_stream_name(0));
        }

      std::swap(streamID1, streamID2);
      vlistID1x = vlistID2;
      vlistID2x = vlistID1;
      std::swap(taxisID1, taxisID2);
      fieldx1 = &field2;
      fieldx2 = &field1;
    }

  if (fillStream1x == false && fillType == FillType::NONE) vlist_compare(vlistID1, vlistID2, CMP_ALL);

  size_t nwpv = (vlistNumber(vlistID1x) == CDI_COMP && vlistNumber(vlistID2x) == CDI_COMP) ? 2 : 1;
  if (nwpv == 2 && !opercplx) cdo_abort("Fields with complex numbers are not supported by this operator!");
  auto gridsizemax = nwpv * vlistGridsizeMax(vlistID1x);

  field1.resize(gridsizemax);
  field2.resize(gridsizemax);

  if (fillStream1x || fillType == FillType::VAR || fillType == FillType::VARTS)
    {
      vardata2.resize(gridsizemax * nlevels2);
      varnmiss2.resize(nlevels2);
    }

  if (Options::cdoVerbose) cdo_print("Number of timesteps: file1 %d, file2 %d", ntsteps1, ntsteps2);

  if (fillType == FillType::NONE)
    {
      if (ntsteps1 != 1 && ntsteps2 == 1)
        {
          fillType = FillType::TS;
          cdo_print("Filling up stream2 >%s< by copying the first timestep.", cdo_get_stream_name(1));
        }
      else if (ntsteps1 == 1 && ntsteps2 != 1)
        {
          fillType = FillType::TS;
          cdo_print("Filling up stream1 >%s< by copying the first timestep.", cdo_get_stream_name(0));

          std::swap(streamID1, streamID2);
          vlistID1x = vlistID2;
          vlistID2x = vlistID1;
          std::swap(taxisID1, taxisID2);
          fieldx1 = &field2;
          fieldx2 = &field1;
        }

      if (fillType == FillType::TS) cdo_fill_ts(vlistID2x, vardata, varnmiss);
    }

  auto streamsSwaped = (fillType == FillType::TS && vlistID1x != vlistID1);

  auto vlistID3 = vlistDuplicate(streamsSwaped ? vlistID1 : vlistID1x);
  if (streamsSwaped)
    {
      auto nvars = vlistNvars(vlistID1);
      for (int varID = 0; varID < nvars; ++varID) vlistDefVarMissval(vlistID3, varID, vlistInqVarMissval(vlistID1, varID));
      for (int varID = 0; varID < nvars; ++varID) vlistDefVarTimetype(vlistID3, varID, vlistInqVarTimetype(vlistID1x, varID));
    }

  if (fillStream1x)
    {
      auto zaxisID2 = vlistZaxis(vlistID2x, 0);
      vlistChangeZaxisIndex(vlistID3, 0, zaxisID2);
    }

  auto taxisID3 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID3, taxisID3);

  auto streamID3 = cdo_open_write(2);
  cdo_def_vlist(streamID3, vlistID3);

  VarList varList1, varList2;
  varListInit(varList1, vlistID1x);
  varListInit(varList2, vlistID2x);

  int nrecs, nrecs2;
  int tsID = 0;
  int tsID2 = 0;
  while (true)
    {
      nrecs = cdo_stream_inq_timestep(streamID1, tsID);

      nrecs2 = 0;
      if (tsID == 0 || fillType == FillType::NONE || fillType == FillType::FILE || fillType == FillType::VARTS)
        {
          nrecs2 = cdo_stream_inq_timestep(streamID2, tsID2);
          if (nrecs2 == 0)
            {
              if (nrecs == 0) break;

              if (fillType == FillType::NONE && streamID2x == streamID2)
                {
                  fillType = FillType::FILE;
                  cdo_print("Filling up stream2 >%s< by copying all timesteps.", cdo_get_stream_name(1));
                }

              if (fillType == FillType::FILE)
                {
                  cdo_stream_close(streamID2);

                  if (cdo_assert_files_only() == false) cdo_abort("infile2 cannot be a pipe in fill mode!");

                  streamID2 = cdo_open_read(1);
                  streamID2x = streamID2;

                  vlistID2 = cdo_stream_inq_vlist(streamID2);

                  vlist_compare(vlistID1, vlistID2, CMP_DIM);

                  tsID2 = 0;
                  nrecs2 = cdo_stream_inq_timestep(streamID2, tsID2);
                  if (nrecs2 == 0) cdo_abort("Empty input stream %s!", cdo_get_stream_name(1));
                }
              else
                cdo_abort("Input streams have different number of timesteps!");
            }
        }

      if (nrecs == 0 || (nrecs2 == 0 && fillType != FillType::TS && fillType != FillType::VAR)) break;

      cdo_taxis_copy_timestep(taxisID3, taxisID1);
      cdo_def_timestep(streamID3, tsID);

      auto numrecs = fillStream1x ? nrecs2 : nrecs;
      for (int recID = 0; recID < numrecs; ++recID)
        {
          auto lread1 = true;
          if (fillStream1x && recID > 0) lread1 = false;
          int varID = -1, levelID;
          if (lread1)
            {
              cdo_inq_record(streamID1, &varID, &levelID);
              cdo_read_record(streamID1, fieldx1->vec_d.data(), &fieldx1->nmiss);

              if (fillStream1x)
                {
                  auto gridsize = nwpv * varList1[varID].gridsize;
                  array_copy(gridsize, fieldx1->vec_d.data(), &vardata2[0]);
                  varnmiss2[0] = fieldx1->nmiss;
                }
            }

          if (fillStream1x) levelID = recID;

          auto varID2 = varID;

          if (tsID == 0 || fillType == FillType::NONE || fillType == FillType::FILE || fillType == FillType::VARTS)
            {
              auto lstatus = (nlevels2 > 1) ? (varID == 0) : (recID == 0);
              if (lstatus || (fillType != FillType::VAR && fillType != FillType::VARTS))
                {
                  cdo_inq_record(streamID2, &varID2, &levelID2);
                  cdo_read_record(streamID2, fieldx2->vec_d.data(), &fieldx2->nmiss);
                  if (varID != varID2) cdo_abort("Internal error, varIDs of input streams differ!");
                  if (fillStream1x == false && levelID != levelID2) cdo_abort("Internal error, levelIDs of input streams differ!");
                }

              if (fillType == FillType::TS)
                {
                  auto gridsize = nwpv * varList2[varID].gridsize;
                  auto offset = gridsize * levelID;
                  array_copy(gridsize, fieldx2->vec_d.data(), &vardata[varID][offset]);
                  varnmiss[varID][levelID] = fieldx2->nmiss;
                }
              else if (lstatus && (fillType == FillType::VAR || fillType == FillType::VARTS))
                {
                  auto gridsize = nwpv * varList2[0].gridsize;
                  auto offset = gridsize * levelID2;
                  array_copy(gridsize, fieldx2->vec_d.data(), &vardata2[offset]);
                  varnmiss2[levelID2] = fieldx2->nmiss;
                }
            }
          else if (fillType == FillType::TS)
            {
              auto gridsize = nwpv * varList2[varID2].gridsize;
              auto offset = gridsize * levelID;
              array_copy(gridsize, &vardata[varID][offset], fieldx2->vec_d.data());
              fieldx2->nmiss = varnmiss[varID][levelID];
            }

          if (fillStream1x)
            {
              auto gridsize = nwpv * varList1[0].gridsize;
              array_copy(gridsize, &vardata2[0], fieldx1->vec_d.data());
              fieldx1->nmiss = varnmiss2[0];
            }

          fieldx1->grid = varList1[varID].gridID;
          fieldx1->missval = varList1[varID].missval;
          fieldx1->nwpv = varList1[varID].nwpv;

          if (fillType == FillType::VAR || fillType == FillType::VARTS)
            {
              levelID2 = (nlevels2 > 1) ? levelID : 0;
              auto gridsize = nwpv * varList2[0].gridsize;
              auto offset = gridsize * levelID2;
              array_copy(gridsize, &vardata2[offset], fieldx2->vec_d.data());
              fieldx2->nmiss = varnmiss2[levelID2];
              fieldx2->grid = varList2[0].gridID;
              fieldx2->missval = varList2[0].missval;
              fieldx2->nwpv = varList2[0].nwpv;
            }
          else
            {
              fieldx2->grid = varList2[varID2].gridID;
              fieldx2->missval = varList2[varID2].missval;
              fieldx2->nwpv = varList2[varID2].nwpv;
            }

          if (nwpv == 2)
            field2_function_complex(field1, field2, operfunc);
          else
            field2_function(field1, field2, operfunc);

          cdo_def_record(streamID3, varID, levelID);
          cdo_write_record(streamID3, field1.vec_d.data(), field1.nmiss);
        }

      tsID++;
      tsID2++;
    }

  if (nrecs == 0 && nrecs2 > 0) cdo_warning("stream2 has more time steps than stream1!");
  // if (nrecs > 0 && nrecs2 == 0) cdo_warning("stream1 has more time steps than stream2!");

  cdo_stream_close(streamID3);
  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  vlistDestroy(vlistID3);

  cdo_finish();

  return nullptr;
}
