/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Yeararith  yearadd         Add yearly time series
      Yeararith  yearsub         Subtract yearly time series
      Yeararith  yearmul         Multiply yearly time series
      Yeararith  yeardiv         Divide yearly time series
*/

#include <cdi.h>

#include <climits>

#include "cdo_vlist.h"
#include "cdo_season.h"
#include "process_int.h"
#include "field_functions.h"

static void
add_operators(void)
{
  // clang-format off
  cdo_operator_add("yearadd",  FieldFunc_Add, 0, nullptr);
  cdo_operator_add("yearsub",  FieldFunc_Sub, 0, nullptr);
  cdo_operator_add("yearmul",  FieldFunc_Mul, 0, nullptr);
  cdo_operator_add("yeardiv",  FieldFunc_Div, 0, nullptr);
  // clang-format on
}
class ModuleYearArith
{
  CdoStreamID streamID1;
  CdoStreamID streamID2;
  CdoStreamID streamID3;

  int taxisID1;
  int taxisID2;
  int taxisID3;

  FieldVector2D vars2;
  Field field;

  int operfunc;
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
    auto vlistID2 = cdo_stream_inq_vlist(streamID2);
    auto vlistID3 = vlistDuplicate(vlistID1);

    vlist_compare(vlistID1, vlistID2, CMP_ALL);

    varListInit(varList1, vlistID1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = vlistInqTaxis(vlistID2);
    taxisID3 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID3, taxisID3);

    streamID3 = cdo_open_write(2);
    cdo_def_vlist(streamID3, vlistID3);

    fields_from_vlist(vlistID2, vars2, FIELD_VEC | FIELD_NAT);
  }
  void
  run()
  {
    int year0 = -INT_MAX + 1;
    int year2last = 0;
    int tsID2 = 0;
    int tsID = 0;
    while (true)
      {
        auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
        if (nrecs == 0) break;

        auto year = taxisInqVdatetime(taxisID1).date.year;
        if (year > year0)
          {
            auto lfound = false;
            while (true)
              {
                auto nrecs2 = cdo_stream_inq_timestep(streamID2, tsID2);
                if (nrecs2 == 0) break;

                tsID2++;
                auto year2 = taxisInqVdatetime(taxisID2).date.year;
                if (year == year2)
                  {
                    lfound = true;
                    year0 = year;
                    for (int recID = 0; recID < nrecs2; ++recID)
                      {
                        int varID, levelID;
                        cdo_inq_record(streamID2, &varID, &levelID);
                        cdo_read_record(streamID2, vars2[varID][levelID]);
                      }
                    break;
                  }

                if (tsID2 > 1 && year2 <= year2last) cdo_abort("stream2 doesn't contain yearly data!");
                year2last = year2;
              }

            if (!lfound) cdo_abort("Data of year %d not found in stream2!", year);
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
  }
  void
  close()
  {
    cdo_stream_close(streamID3);
    cdo_stream_close(streamID2);
    cdo_stream_close(streamID1);

    cdo_finish();
  }
};
void *
Yeararith(void *process)
{
  ModuleYearArith yeararith;
  yeararith.init(process);
  yeararith.run();
  yeararith.close();
  return nullptr;
}
