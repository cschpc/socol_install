/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Compc      eqc             Equal constant
      Compc      nec             Not equal constant
      Compc      lec             Less equal constant
      Compc      ltc             Less then constant
      Compc      gec             Greater equal constant
      Compc      gtc             Greater then constant
*/

#include <cdi.h>

#include "process_int.h"
#include "cdo_vlist.h"
#include "param_conversion.h"
#include "field_functions.h"

static auto func_compc = [](auto hasMissvals, auto n, auto mv, const auto &vIn, const auto cVal, auto &vOut, auto binray_operator) {
  if (hasMissvals)
    {
      const auto c_is_missval = dbl_is_equal(cVal, mv);
      if (std::isnan(mv))
        for (size_t i = 0; i < n; ++i) vOut[i] = (dbl_is_equal(vIn[i], mv) || c_is_missval) ? mv : binray_operator(vIn[i], cVal);
      else
        for (size_t i = 0; i < n; ++i) vOut[i] = (is_equal(vIn[i], mv) || c_is_missval) ? mv : binray_operator(vIn[i], cVal);
    }
  else
    {
      for (size_t i = 0; i < n; ++i) vOut[i] = binray_operator(vIn[i], cVal);
    }
};

void *
Compc(void *process)
{
  cdo_initialize(process);

  cdo_operator_add("eqc", FieldFunc_EQ, 0, nullptr);
  cdo_operator_add("nec", FieldFunc_NE, 0, nullptr);
  cdo_operator_add("lec", FieldFunc_LE, 0, nullptr);
  cdo_operator_add("ltc", FieldFunc_LT, 0, nullptr);
  cdo_operator_add("gec", FieldFunc_GE, 0, nullptr);
  cdo_operator_add("gtc", FieldFunc_GT, 0, nullptr);

  auto operatorID = cdo_operator_id();
  auto operFunc = cdo_operator_f1(operatorID);

  operator_input_arg("constant value");
  auto rc = parameter_to_double(cdo_operator_argv(0));

  auto streamID1 = cdo_open_read(0);

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto vlistID2 = vlistDuplicate(vlistID1);

  auto taxisID1 = vlistInqTaxis(vlistID1);
  auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  VarList varList;
  varListInit(varList, vlistID1);

  auto gridsizemax = vlistGridsizeMax(vlistID1);
  Varray<double> vaIn(gridsizemax), vaOut(gridsizemax);

  auto streamID2 = cdo_open_write(1);
  cdo_def_vlist(streamID2, vlistID2);

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
          size_t nmiss;
          cdo_read_record(streamID1, vaIn.data(), &nmiss);

          auto missval = varList[varID].missval;
          auto ngp = varList[varID].gridsize;
          auto datatype = varList[varID].datatype;
          double rcv = (datatype == CDI_DATATYPE_FLT32) ? (float) rc : rc;

          if (nmiss > 0) cdo_check_missval(missval, varList[varID].name);

          auto hasMissvals = (nmiss > 0 || dbl_is_equal(rc, missval));
          // clang-format off
          if      (operFunc == FieldFunc_EQ) func_compc(hasMissvals, ngp, missval, vaIn, rcv, vaOut, binary_op_EQ);
          else if (operFunc == FieldFunc_NE) func_compc(hasMissvals, ngp, missval, vaIn, rcv, vaOut, binary_op_NE);
          else if (operFunc == FieldFunc_LE) func_compc(hasMissvals, ngp, missval, vaIn, rcv, vaOut, binary_op_LE);
          else if (operFunc == FieldFunc_LT) func_compc(hasMissvals, ngp, missval, vaIn, rcv, vaOut, binary_op_LT);
          else if (operFunc == FieldFunc_GE) func_compc(hasMissvals, ngp, missval, vaIn, rcv, vaOut, binary_op_GE);
          else if (operFunc == FieldFunc_GT) func_compc(hasMissvals, ngp, missval, vaIn, rcv, vaOut, binary_op_GT);
          else cdo_abort("Operator not implemented!");
          // clang-format on

          nmiss = varray_num_mv(ngp, vaOut, missval);
          cdo_def_record(streamID2, varID, levelID);
          cdo_write_record(streamID2, vaOut.data(), nmiss);
        }

      tsID++;
    }

  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  cdo_finish();

  return nullptr;
}
