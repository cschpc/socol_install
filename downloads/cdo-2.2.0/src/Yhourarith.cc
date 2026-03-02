/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Yhourarith  yhouradd         Add multi-year hourly time series
      Yhourarith  yhoursub         Subtract multi-year hourly time series
      Yhourarith  yhourmul         Multiply multi-year hourly time series
      Yhourarith  yhourdiv         Divide multi-year hourly time series
*/

#include <cdi.h>

#include "cdo_vlist.h"
#include "datetime.h"
#include "field.h"
#include "process_int.h"
#include "printinfo.h"
#include "field_functions.h"

constexpr int MaxHours = 9301;  // 31*12*25 + 1

static int
getHourOfYearIndex(const CdiDateTime &vDateTime)
{
  auto houroy = decode_hour_of_year(vDateTime);

  if (houroy < 0 || houroy >= MaxHours) cdo_abort("Hour of year %d out of range (%s)!", houroy, datetime_to_string(vDateTime));

  return houroy;
}

static void
add_operators(void)
{
  cdo_operator_add("yhouradd", FieldFunc_Add, 0, nullptr);
  cdo_operator_add("yhoursub", FieldFunc_Sub, 0, nullptr);
  cdo_operator_add("yhourmul", FieldFunc_Mul, 0, nullptr);
  cdo_operator_add("yhourdiv", FieldFunc_Div, 0, nullptr);
}

class ModuleYhourArith
{
  int operfunc;
  CdoStreamID streamID1;
  CdoStreamID streamID2;
  CdoStreamID streamID3;

  int taxisID1;
  int taxisID2;
  int taxisID3;

  int vlistID2;

  Field field;
  FieldVector2D vars2[MaxHours];
  VarList varList1;

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

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
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
  void
  run()
  {
    int tsID = 0;
    while (true)
      {
        auto nrecs = cdo_stream_inq_timestep(streamID2, tsID);
        if (nrecs == 0) break;

        auto houroy = getHourOfYearIndex(taxisInqVdatetime(taxisID2));
        if (vars2[houroy].size() > 0) cdo_abort("Hour of year index %d already allocated!", houroy);

        fields_from_vlist(vlistID2, vars2[houroy], FIELD_VEC | FIELD_NAT);

        for (int recID = 0; recID < nrecs; ++recID)
          {
            int varID, levelID;
            cdo_inq_record(streamID2, &varID, &levelID);
            cdo_read_record(streamID2, vars2[houroy][varID][levelID]);
          }

        tsID++;
      }

    cdo_stream_close(streamID2);

    tsID = 0;
    while (true)
      {
        auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
        if (nrecs == 0) break;

        auto houroy = getHourOfYearIndex(taxisInqVdatetime(taxisID1));
        if (vars2[houroy].size() == 0) cdo_abort("Hour of year index %d not found!", houroy);

        cdo_taxis_copy_timestep(taxisID3, taxisID1);
        cdo_def_timestep(streamID3, tsID);

        for (int recID = 0; recID < nrecs; ++recID)
          {
            int varID, levelID;
            cdo_inq_record(streamID1, &varID, &levelID);
            field.init(varList1[varID]);
            cdo_read_record(streamID1, field);

            field2_function(field, vars2[houroy][varID][levelID], operfunc);

            cdo_def_record(streamID3, varID, levelID);
            cdo_write_record(streamID3, field);
          }

        tsID++;
      }
  }
  void
  close()
  {

    cdo_stream_close(streamID3);
    cdo_stream_close(streamID1);

    cdo_finish();
  }
};

void *
Yhourarith(void *process)
{
  ModuleYhourArith yhourarith;
  yhourarith.init(process);
  yhourarith.run();
  yhourarith.close();

  return nullptr;
}
