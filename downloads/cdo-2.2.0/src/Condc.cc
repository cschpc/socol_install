/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Condc      ifthenc         If then constant
      Condc      ifnotthenc      If not then constant
*/

#include <cdi.h>

#include "process_int.h"
#include "cdo_vlist.h"
#include "param_conversion.h"

static void
operator_IFTHENC(size_t n, double mv, const Varray<double> &vIn, const double cVal, Varray<double> &vOut)
{
  if (std::isnan(mv))
    for (size_t i = 0; i < n; ++i) vOut[i] = (!dbl_is_equal(vIn[i], mv) && !dbl_is_equal(vIn[i], 0.0)) ? cVal : mv;
  else
    for (size_t i = 0; i < n; ++i) vOut[i] = (!is_equal(vIn[i], mv) && !is_equal(vIn[i], 0.0)) ? cVal : mv;
}

static void
operator_IFNOTTHENC(size_t n, double mv, const Varray<double> &vIn, const double cVal, Varray<double> &vOut)
{
  if (std::isnan(mv))
    for (size_t i = 0; i < n; ++i) vOut[i] = (!dbl_is_equal(vIn[i], mv) && dbl_is_equal(vIn[i], 0.0)) ? cVal : mv;
  else
    for (size_t i = 0; i < n; ++i) vOut[i] = (!is_equal(vIn[i], mv) && is_equal(vIn[i], 0.0)) ? cVal : mv;
}

void *
Condc(void *process)
{
  cdo_initialize(process);

  // clang-format off
  const auto IFTHENC    = cdo_operator_add("ifthenc",    0, 0, nullptr);
  const auto IFNOTTHENC = cdo_operator_add("ifnotthenc", 0, 0, nullptr);
  // clang-format on

  const auto operatorID = cdo_operator_id();

  operator_input_arg("constant value");
  const auto rc = parameter_to_double(cdo_operator_argv(0));

  const auto streamID1 = cdo_open_read(0);

  const auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  const auto vlistID2 = vlistDuplicate(vlistID1);

  const auto taxisID1 = vlistInqTaxis(vlistID1);
  const auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  VarList varList;
  varListInit(varList, vlistID1);

  const auto gridsizemax = vlistGridsizeMax(vlistID1);
  Varray<double> vaIn(gridsizemax), vaOut(gridsizemax);

  const auto streamID2 = cdo_open_write(1);
  cdo_def_vlist(streamID2, vlistID2);

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
          cdo_read_record(streamID1, vaIn.data(), &nmiss);

          const auto missval = varList[varID].missval;
          const auto ngp = varList[varID].gridsize;

          if (nmiss > 0) cdo_check_missval(missval, varList[varID].name);

          // clang-format off
          if      (operatorID == IFTHENC)    operator_IFTHENC(ngp, missval, vaIn, rc, vaOut);
          else if (operatorID == IFNOTTHENC) operator_IFNOTTHENC(ngp, missval, vaIn, rc, vaOut);
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
