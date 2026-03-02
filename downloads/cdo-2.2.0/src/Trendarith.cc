/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

   Trendarith    addtrend        Add trend
   Trendarith    subtrend        Subtract trend
*/

#include <cdi.h>

#include "process_int.h"
#include "cdo_vlist.h"
#include "datetime.h"
#include "pmlist.h"
#include "param_conversion.h"
#include "field_functions.h"

template <typename T>
static void
add_trend(double zj, const Varray<T> &v1, const Varray<double> &v2, const Varray<double> &v3, Varray<T> &v4, size_t n, double mv)
{
  auto missval1 = mv;
  auto missval2 = mv;
  for (size_t i = 0; i < n; ++i) v4[i] = ADDMN(v1[i], ADDMN(v2[i], MULMN(v3[i], zj)));
}

static void
add_trend(double zj, const Field &field1, const Field &field2, const Field &field3, Field &field4)
{
  if (field1.memType != field4.memType) cdo_abort("Interal error, memType of field1 and field2 differ!");

  if (field1.memType == MemType::Float)
    add_trend(zj, field1.vec_f, field2.vec_d, field3.vec_d, field4.vec_f, field1.size, field1.missval);
  else
    add_trend(zj, field1.vec_d, field2.vec_d, field3.vec_d, field4.vec_d, field1.size, field1.missval);
}

template <typename T>
static void
sub_trend(double zj, const Varray<T> &v1, const Varray<double> &v2, const Varray<double> &v3, Varray<T> &v4, size_t n, double mv)
{
  auto missval1 = mv;
  auto missval2 = mv;
  for (size_t i = 0; i < n; ++i) v4[i] = SUBMN(v1[i], ADDMN(v2[i], MULMN(v3[i], zj)));
}

static void
sub_trend(double zj, const Field &field1, const Field &field2, const Field &field3, Field &field4)
{
  if (field1.memType != field4.memType) cdo_abort("Interal error, memType of field1 and field2 differ!");

  if (field1.memType == MemType::Float)
    sub_trend(zj, field1.vec_f, field2.vec_d, field3.vec_d, field4.vec_f, field1.size, field1.missval);
  else
    sub_trend(zj, field1.vec_d, field2.vec_d, field3.vec_d, field4.vec_d, field1.size, field1.missval);
}

static void
trendarithGetParameter(bool &tstepIsEqual)
{
  auto pargc = cdo_operator_argc();
  if (pargc)
    {
      const auto &pargv = cdo_get_oper_argv();

      KVList kvlist;
      kvlist.name = cdo_module_name();
      if (kvlist.parse_arguments(pargc, pargv) != 0) cdo_abort("Parse error!");
      if (Options::cdoVerbose) kvlist.print();

      for (const auto &kv : kvlist)
        {
          const auto &key = kv.key;
          if (kv.nvalues > 1) cdo_abort("Too many values for parameter key >%s<!", key);
          if (kv.nvalues < 1) cdo_abort("Missing value for parameter key >%s<!", key);
          const auto &value = kv.values[0];

          // clang-format off
          if (key == "equal") tstepIsEqual = parameter_to_bool(value);
          else cdo_abort("Invalid parameter key >%s<!", key);
          // clang-format on
        }
    }
}

void *
Trendarith(void *process)
{
  cdo_initialize(process);

  cdo_operator_add("addtrend", FieldFunc_Add, 0, nullptr);
  cdo_operator_add("subtrend", FieldFunc_Sub, 0, nullptr);

  auto operatorID = cdo_operator_id();
  auto operfunc = cdo_operator_f1(operatorID);

  auto tstepIsEqual = true;
  trendarithGetParameter(tstepIsEqual);

  auto streamID1 = cdo_open_read(0);
  auto streamID2 = cdo_open_read(1);
  auto streamID3 = cdo_open_read(2);

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto vlistID2 = cdo_stream_inq_vlist(streamID2);
  auto vlistID3 = cdo_stream_inq_vlist(streamID3);
  auto vlistID4 = vlistDuplicate(vlistID1);

  vlist_compare(vlistID1, vlistID2, CMP_DIM);
  vlist_compare(vlistID1, vlistID3, CMP_DIM);

  auto taxisID1 = vlistInqTaxis(vlistID1);
  auto taxisID4 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID4, taxisID4);

  VarList varList;
  varListInit(varList, vlistID1);

  Field field1, field4;

  auto streamID4 = cdo_open_write(3);
  cdo_def_vlist(streamID4, vlistID4);

  FieldVector2D vars2, vars3;
  fields_from_vlist(vlistID1, vars2, FIELD_VEC);
  fields_from_vlist(vlistID1, vars3, FIELD_VEC);

  {
    auto nrecs = cdo_stream_inq_timestep(streamID2, 0);
    for (int recID = 0; recID < nrecs; ++recID)
      {
        int varID, levelID;
        cdo_inq_record(streamID2, &varID, &levelID);
        cdo_read_record(streamID2, vars2[varID][levelID].vec_d.data(), &vars2[varID][levelID].nmiss);
      }
  }

  {
    auto nrecs = cdo_stream_inq_timestep(streamID3, 0);
    for (int recID = 0; recID < nrecs; ++recID)
      {
        int varID, levelID;
        cdo_inq_record(streamID3, &varID, &levelID);
        cdo_read_record(streamID3, vars3[varID][levelID].vec_d.data(), &vars3[varID][levelID].nmiss);
      }
  }

  auto calendar = taxisInqCalendar(taxisID1);

  CheckTimeIncr checkTimeIncr;
  JulianDate julianDate0;
  double deltat1 = 0.0;
  int tsID = 0;
  while (true)
    {
      auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
      if (nrecs == 0) break;

      auto vDateTime = taxisInqVdatetime(taxisID1);

      if (tstepIsEqual) check_time_increment(tsID, calendar, vDateTime, checkTimeIncr);
      auto zj = tstepIsEqual ? (double) tsID : delta_time_step_0(tsID, calendar, vDateTime, julianDate0, deltat1);

      cdo_taxis_copy_timestep(taxisID4, taxisID1);
      cdo_def_timestep(streamID4, tsID);

      for (int recID = 0; recID < nrecs; ++recID)
        {
          int varID, levelID;
          cdo_inq_record(streamID1, &varID, &levelID);
          field1.init(varList[varID]);
          cdo_read_record(streamID1, field1);

          field4.init(varList[varID]);

          if (operfunc == FieldFunc_Add)
            add_trend(zj, field1, vars2[varID][levelID], vars3[varID][levelID], field4);
          else
            sub_trend(zj, field1, vars2[varID][levelID], vars3[varID][levelID], field4);

          field4.nmiss = field1.nmiss;
          cdo_def_record(streamID4, varID, levelID);
          cdo_write_record(streamID4, field4);
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
