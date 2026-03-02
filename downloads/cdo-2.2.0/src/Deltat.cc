/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

     Deltat    deltat         Delta t
*/

#include "cdi.h"
#include "julian_date.h"
#include "process_int.h"
#include "field_functions.h"

template <typename T>
static void
varray_deltat(const size_t len, const Varray<T> &v0, const Varray<T> &v1, Varray<T> &v2, double factor, double mv)
{
  assert(len > 0);
  assert(v0.size() > 0);
  assert(v1.size() > 0);
  assert(v2.size() > 0);
  assert(len <= v0.size());
  assert(len <= v1.size());
  assert(len <= v2.size());

  for (size_t i = 0; i < len; ++i)
    {
      if (DBL_IS_EQUAL(v0[i], mv) || DBL_IS_EQUAL(v1[i], mv))
        v2[i] = mv;
      else
        v2[i] = (v1[i] - v0[i]) * factor;
    }
}

template <typename T>
static void
varray_deltat(const size_t len, const Varray<T> &v0, const Varray<T> &v1, Varray<T> &v2, double factor)
{
  assert(len > 0);

  for (size_t i = 0; i < len; ++i) v2[i] = (v1[i] - v0[i]) * factor;
}

void
field_deltat(const Field &field0, const Field &field1, Field &field2, double factor)
{
  if (field1.memType != field2.memType) cdo_abort("Interal error, memType of field1 and field2 differ!");

  if (field1.nmiss || field0.nmiss)
    {
      if (field1.memType == MemType::Float)
        varray_deltat(field1.size, field0.vec_f, field1.vec_f, field2.vec_f, factor, field1.missval);
      else
        varray_deltat(field1.size, field0.vec_d, field1.vec_d, field2.vec_d, factor, field1.missval);

      field2.nmiss = field_num_mv(field2);
    }
  else
    {
      if (field1.memType == MemType::Float)
        varray_deltat(field1.size, field0.vec_f, field1.vec_f, field2.vec_f, factor);
      else
        varray_deltat(field1.size, field0.vec_d, field1.vec_d, field2.vec_d, factor);
    }
}

class ModuleDeltat
{
  CdoStreamID streamID1;
  CdoStreamID streamID2;

  int taxisID1;
  int taxisID2;
  int calendar;
  VarList varList1;
  FieldVector2D vars;

  Field field1, field2;
  bool ldivdt;

public:
  void
  init(void *process)
  {
    cdo_initialize(process);

    // clang-format off
  cdo_operator_add("deltat",            0,     0, nullptr);
  cdo_operator_add("timederivative",    0,     1, nullptr);
    // clang-format on

    auto operatorID = cdo_operator_id();
    ldivdt = cdo_operator_f2(operatorID);

    operator_check_argc(0);

    streamID1 = cdo_open_read(0);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = vlistDuplicate(vlistID1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    calendar = taxisInqCalendar(taxisID1);

    varListInit(varList1, vlistID1);

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);

    fields_from_vlist(vlistID1, vars, FIELD_VEC | FIELD_NAT);
  }
  void
  run()
  {
    int tsID = 0;
    auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);

    auto julianDate0 = julianDate_encode(calendar, taxisInqVdatetime(taxisID1));

    for (int recID = 0; recID < nrecs; ++recID)
      {
        int varID, levelID;
        cdo_inq_record(streamID1, &varID, &levelID);
        cdo_read_record(streamID1, vars[varID][levelID]);
      }

    tsID++;
    int tsID2 = 0;
    while (true)
      {
        nrecs = cdo_stream_inq_timestep(streamID1, tsID);
        if (nrecs == 0) break;

        auto julianDate1 = julianDate_encode(calendar, taxisInqVdatetime(taxisID1));
        auto idtInSec = ldivdt ? 1.0 / julianDate_to_seconds(julianDate_sub(julianDate1, julianDate0)) : 1.0;
        julianDate0 = julianDate1;

        cdo_taxis_copy_timestep(taxisID2, taxisID1);
        cdo_def_timestep(streamID2, tsID2);

        for (int recID = 0; recID < nrecs; ++recID)
          {
            int varID, levelID;
            cdo_inq_record(streamID1, &varID, &levelID);
            field1.init(varList1[varID]);
            cdo_read_record(streamID1, field1);

            auto &field0 = vars[varID][levelID];

            field2.init(varList1[varID]);
            field_deltat(field0, field1, field2, idtInSec);

            field_copy(field1, field0);

            cdo_def_record(streamID2, varID, levelID);
            cdo_write_record(streamID2, field2);
          }

        tsID++;
        tsID2++;
      }
  }

  void
  close()
  {
    cdo_stream_close(streamID2);
    cdo_stream_close(streamID1);

    cdo_finish();
  }
};
void *
Deltat(void *process)
{
  ModuleDeltat deltat;
  deltat.init(process);
  deltat.run();
  deltat.close();

  return nullptr;
}
