/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Ydayarith  ydayadd         Add multi-year daily time series
      Ydayarith  ydaysub         Subtract multi-year daily time series
      Ydayarith  ydaymul         Multiply multi-year daily time series
      Ydayarith  ydaydiv         Divide multi-year daily time series
*/

#include <cdi.h>

#include "cdo_vlist.h"
#include "datetime.h"
#include "process_int.h"
#include "printinfo.h"
#include "field_functions.h"

static void
add_operators(void)
{
  cdo_operator_add("ydayadd", FieldFunc_Add, 0, nullptr);
  cdo_operator_add("ydaysub", FieldFunc_Sub, 0, nullptr);
  cdo_operator_add("ydaymul", FieldFunc_Mul, 0, nullptr);
  cdo_operator_add("ydaydiv", FieldFunc_Div, 0, nullptr);
}

class ModuleYdayarith
{
private:
  static const int MaxDays = 373;  //~31*12
  Field field;
  FieldVector2D vars2[MaxDays];
  VarList varList1;

  int operfunc;

  CdoStreamID streamID1;
  CdoStreamID streamID2;
  CdoStreamID streamID3;

  int vlistID1;
  int vlistID2;

  int taxisID1;
  int taxisID2;
  int taxisID3;

public:
  void
  init(void *process)
  {
    cdo_initialize(process);

    add_operators();

    auto operatorID = cdo_operator_id();
    operfunc = cdo_operator_f1(operatorID);

    operator_check_argc(0);

    streamID1 = cdo_open_read(0);
    streamID2 = cdo_open_read(1);

    vlistID1 = cdo_stream_inq_vlist(streamID1);
    vlistID2 = cdo_stream_inq_vlist(streamID2);
    auto vlistID3 = vlistDuplicate(vlistID1);

    vlist_compare(vlistID1, vlistID2, CMP_ALL);

    varListInit(varList1, vlistID1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = vlistInqTaxis(vlistID2);
    taxisID3 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID3, taxisID3);

    streamID3 = cdo_open_write(2);
    cdo_def_vlist(streamID3, vlistID3);
  }
  int
  run()
  {
    int tsID = 0;
    while (true)
      {
        auto nrecs = cdo_stream_inq_timestep(streamID2, tsID);
        if (nrecs == 0) break;

        auto vDateTime = taxisInqVdatetime(taxisID2);
        int dayOfYear = decode_day_of_year(vDateTime.date);
        // assert(dayOfYear < 1 || dayOfYear >= MaxDays);
        if (dayOfYear == 0)
          {
            cdo_error("Day of year %d out of range (date=%s)!", dayOfYear, date_to_string(vDateTime.date));
            return -1;
          }
        if (vars2[dayOfYear].size() > 0)
          {
            cdo_error("Day of year index %d already allocated (date=%s)! Each day of year must only exist once", dayOfYear,
                      date_to_string(vDateTime.date));
            return -1;
          }

        fields_from_vlist(vlistID2, vars2[dayOfYear], FIELD_VEC | FIELD_NAT);

        for (int recID = 0; recID < nrecs; ++recID)
          {
            int varID, levelID;
            cdo_inq_record(streamID2, &varID, &levelID);
            cdo_read_record(streamID2, vars2[dayOfYear][varID][levelID]);
          }

        tsID++;
      }

    tsID = 0;
    while (true)
      {
        auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
        if (nrecs == 0) break;

        auto vDateTime = taxisInqVdatetime(taxisID1);
        int dayOfYear = decode_day_of_year(vDateTime.date);
        // assert(dayOfYear < 1 || dayOfYear >= MaxDays);
        if (dayOfYear == 0)
          {
            cdo_error("Day of year %d out of range (date=%s)!", dayOfYear, date_to_string(vDateTime.date));
            return -1;
          }
        if (vars2[dayOfYear].size() == 0)
          {
            cdo_error("Day of year index %d not found (date=%s)!", dayOfYear, date_to_string(vDateTime.date));
            return -1;
          }

        cdo_taxis_copy_timestep(taxisID3, taxisID1);
        cdo_def_timestep(streamID3, tsID);

        for (int recID = 0; recID < nrecs; ++recID)
          {
            int varID, levelID;
            cdo_inq_record(streamID1, &varID, &levelID);
            field.init(varList1[varID]);
            cdo_read_record(streamID1, field);

            field2_function(field, vars2[dayOfYear][varID][levelID], operfunc);

            cdo_def_record(streamID3, varID, levelID);
            cdo_write_record(streamID3, field);
          }
        tsID++;
      }
    return 0;
  }
  void
  close()
  {
    cdo_stream_close(streamID2);
    cdo_stream_close(streamID3);
    cdo_stream_close(streamID1);

    cdo_finish();
  }
};
void *
Ydayarith(void *process)
{
  ModuleYdayarith ydayarith;
  ydayarith.init(process);
  int success = ydayarith.run();
  ydayarith.close();

  return nullptr;
}
