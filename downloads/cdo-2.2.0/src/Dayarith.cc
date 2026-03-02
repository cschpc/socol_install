/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Dayarith  dayadd         Add daily time series
      Dayarith  daysub         Subtract daily time series
      Dayarith  daymul         Multiply daily time series
      Dayarith  daydiv         Divide daily time series
*/

#include <cdi.h>

#include "process_int.h"
#include "cdo_vlist.h"
#include "field_functions.h"

static void
add_operators(void)
{
  cdo_operator_add("dayadd", FieldFunc_Add, 0, nullptr);
  cdo_operator_add("daysub", FieldFunc_Sub, 0, nullptr);
  cdo_operator_add("daymul", FieldFunc_Mul, 0, nullptr);
  cdo_operator_add("daydiv", FieldFunc_Div, 0, nullptr);
}

void *
Dayarith(void *process)
{
  cdo_initialize(process);

  add_operators();

  auto operatorID = cdo_operator_id();
  auto operfunc = cdo_operator_f1(operatorID);

  operator_check_argc(0);

  auto streamID1 = cdo_open_read(0);
  auto streamID2 = cdo_open_read(1);

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto vlistID2 = cdo_stream_inq_vlist(streamID2);
  auto vlistID3 = vlistDuplicate(vlistID1);

  vlist_compare(vlistID1, vlistID2, CMP_ALL);

  VarList varList1;
  varListInit(varList1, vlistID1);

  Field field;

  auto taxisID1 = vlistInqTaxis(vlistID1);
  auto taxisID2 = vlistInqTaxis(vlistID2);
  auto taxisID3 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID3, taxisID3);

  auto streamID3 = cdo_open_write(2);
  cdo_def_vlist(streamID3, vlistID3);

  FieldVector2D vars2;
  fields_from_vlist(vlistID2, vars2, FIELD_VEC | FIELD_NAT);

  CdiDate vDate2 = {};
  int tsID = 0;
  int tsID2 = 0;
  while (true)
    {
      auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
      if (nrecs == 0) break;

      auto vDate1 = taxisInqVdatetime(taxisID1).date;

      if (!cdiDate_isEQ(vDate1, vDate2))
        {
          int year1, month1, day1;
          cdiDate_decode(vDate1, &year1, &month1, &day1);
          if (Options::cdoVerbose) cdo_print("Process: Year=%4d  Month=%2d  Day=%2d", year1, month1, day1);

          auto nrecs2 = cdo_stream_inq_timestep(streamID2, tsID2);
          if (nrecs2 == 0) cdo_abort("Missing year=%4d mon=%2d day=%2d in %s!", year1, month1, day1, cdo_get_stream_name(1));

          vDate2 = taxisInqVdatetime(taxisID2).date;
          if (!cdiDate_isEQ(vDate1, vDate2))
            {
              int year2, month2, day2;
              cdiDate_decode(vDate2, &year2, &month2, &day2);
              cdo_abort("Timestep %d in %s has wrong date! Current year=%4d mon=%2d day=%2d, expected year=%4d mon=%2d day=%2d",
                        tsID2 + 1, cdo_get_stream_name(1), year2, month2, day2, year1, month1, day1);
            }

          for (int recID = 0; recID < nrecs2; ++recID)
            {
              int varID, levelID;
              cdo_inq_record(streamID2, &varID, &levelID);
              cdo_read_record(streamID2, vars2[varID][levelID]);
            }

          tsID2++;
        }

      cdo_taxis_copy_timestep(taxisID3, taxisID1);
      cdo_def_timestep(streamID3, tsID);

      for (int recID = 0; recID < nrecs; ++recID)
        {
          int varID, levelID;
          cdo_inq_record(streamID1, &varID, &levelID);
          field.init(varList1[varID]);
          cdo_read_record(streamID1, field);

          field2_function(field, vars2[varID][levelID], operfunc);

          cdo_def_record(streamID3, varID, levelID);
          cdo_write_record(streamID3, field);
        }

      tsID++;
    }

  cdo_stream_close(streamID3);
  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  cdo_finish();

  return nullptr;
}
