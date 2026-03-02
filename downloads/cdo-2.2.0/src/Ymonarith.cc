/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Ymonarith  ymonadd         Add multi-year monthly time series
      Ymonarith  ymonsub         Subtract multi-year monthly time series
      Ymonarith  ymonmul         Multiply multi-year monthly time series
      Ymonarith  ymondiv         Divide multi-year monthly time series
      Ymonarith  yseasadd        Add multi-year seasonal time series
      Ymonarith  yseassub        Subtract multi-year seasonal time series
      Ymonarith  yseasmul        Multiply multi-year seasonal time series
      Ymonarith  yseasdiv        Divide multi-year seasonal time series
*/

#include <cdi.h>

#include "cdo_vlist.h"
#include "cdo_season.h"
#include "process_int.h"
#include "field_functions.h"

constexpr int MaxMonths = 12;

static const char *monthNames[]
    = { "January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December" };

static int
get_month_index(const CdiDate vDate)
{
  int year, mon, day;
  cdiDate_decode(vDate, &year, &mon, &day);
  if (mon < 1 || mon > MaxMonths) cdo_abort("Month %d out of range!", mon);
  mon--;
  return mon;
}

enum
{
  MONTHLY,
  SEASONAL
};

static void
add_operators(void)
{
  // clang-format off
  cdo_operator_add("ymonadd",  FieldFunc_Add, MONTHLY, nullptr);
  cdo_operator_add("ymonsub",  FieldFunc_Sub, MONTHLY, nullptr);
  cdo_operator_add("ymonmul",  FieldFunc_Mul, MONTHLY, nullptr);
  cdo_operator_add("ymondiv",  FieldFunc_Div, MONTHLY, nullptr);
  cdo_operator_add("yseasadd", FieldFunc_Add, SEASONAL, nullptr);
  cdo_operator_add("yseassub", FieldFunc_Sub, SEASONAL, nullptr);
  cdo_operator_add("yseasmul", FieldFunc_Mul, SEASONAL, nullptr);
  cdo_operator_add("yseasdiv", FieldFunc_Div, SEASONAL, nullptr);
  // clang-format on
}

class ModuleYmonArith
{
  int operfunc;
  int opertype;

  CdoStreamID streamID1;
  CdoStreamID streamID2;
  CdoStreamID streamID3;

  int taxisID1;
  int taxisID2;
  int taxisID3;

  int vlistID2;

  VarList varList1;
  Field field;
  FieldVector2D vars2[MaxMonths];

  void
  process_data(int tsID, int nrecs, int mon)
  {
    cdo_taxis_copy_timestep(taxisID3, taxisID1);
    cdo_def_timestep(streamID3, tsID);

    for (int recID = 0; recID < nrecs; ++recID)
      {
        int varID, levelID;
        cdo_inq_record(streamID1, &varID, &levelID);
        field.init(varList1[varID]);
        cdo_read_record(streamID1, field);

        field2_function(field, vars2[mon][varID][levelID], operfunc);

        cdo_def_record(streamID3, varID, levelID);
        cdo_write_record(streamID3, field);
      }
  }

  void
  load_data_from_stream(int mon, int nrecs)
  {
    fields_from_vlist(vlistID2, vars2[mon], FIELD_VEC | FIELD_NAT);

    for (int recID = 0; recID < nrecs; ++recID)
      {
        int varID, levelID;
        cdo_inq_record(streamID2, &varID, &levelID);
        cdo_read_record(streamID2, vars2[mon][varID][levelID]);
      }
  }

  void
  run_monthly()
  {
    for (int nrecs, tsID = 0; (nrecs = cdo_stream_inq_timestep(streamID2, tsID) != 0); tsID++)
      {
        auto mon = get_month_index(taxisInqVdatetime(taxisID2).date);
        if (vars2[mon].size())
          cdo_abort("%s already allocated! The second input file must contain monthly mean values for a maximum of one year.",
                    monthNames[mon]);

        load_data_from_stream(mon, nrecs);
      }

    for (int nrecs, tsID = 0; (nrecs = cdo_stream_inq_timestep(streamID1, tsID) != 0); tsID++)
      {
        auto mon = get_month_index(taxisInqVdatetime(taxisID1).date);
        if (vars2[mon].size() == 0)
          cdo_abort("%s not found! The second input file must contain monthly mean values for a maximum of one year.",
                    monthNames[mon]);

        process_data(tsID, nrecs, mon);
      }
  }

  void
  run_seasonal()
  {
    const char **seasonNames = get_season_name();
    for (int nrecs, tsID = 0; (nrecs = cdo_stream_inq_timestep(streamID2, tsID) != 0); tsID++)
      {
        auto mon = get_month_index(taxisInqVdatetime(taxisID2).date);
        auto season = month_to_season(mon + 1);
        if (vars2[season].size()) cdo_abort("Season %s already allocated!", seasonNames[season]);

        load_data_from_stream(mon, nrecs);
      }

    for (int nrecs, tsID = 0; (nrecs = cdo_stream_inq_timestep(streamID1, tsID) != 0); tsID++)
      {
        auto mon = get_month_index(taxisInqVdatetime(taxisID1).date);
        auto season = month_to_season(mon + 1);
        if (vars2[season].size() == 0) cdo_abort("Season %s not found!", seasonNames[season]);

        process_data(tsID, nrecs, season);
      }
  }

public:
  void
  init(void *process)
  {
    cdo_initialize(process);

    add_operators();

    auto operatorID = cdo_operator_id();
    opertype = cdo_operator_f2(operatorID);
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
    if (opertype == SEASONAL)
      run_seasonal();
    else
      run_monthly();
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
Ymonarith(void *process)
{
  ModuleYmonArith ymonarith;
  ymonarith.init(process);
  ymonarith.run();
  ymonarith.close();

  return nullptr;
}
