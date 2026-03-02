/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Copyright (C) 2006 Brockmann Consult

  Author: Ralf Quast
          Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Timcount    timcount          Time counts
      Hourcount   hourcount         Hourly counts
      Daycount    daycount          Daily counts
      Moncount    moncount          Monthly counts
      Yearcount   yearcount         Yearly counts
*/

#include <cdi.h>

#include "process_int.h"
#include "util_date.h"
#include "field_functions.h"

void *
Timcount(void *process)
{
  CdiDateTime vDateTime0{};
  CdiDateTime vDateTimeN{};

  cdo_initialize(process);

  // clang-format off
  cdo_operator_add("timcount",  0, CMP_DATE, nullptr);
  cdo_operator_add("yearcount", 0, CMP_YEAR, nullptr);
  cdo_operator_add("moncount",  0, CMP_MONTH, nullptr);
  cdo_operator_add("daycount",  0, CMP_DAY, nullptr);
  cdo_operator_add("hourcount", 0, CMP_HOUR, nullptr);
  // clang-format on

  auto operatorID = cdo_operator_id();
  auto compareDate = cdo_operator_f2(operatorID);

  operator_check_argc(0);

  auto streamID1 = cdo_open_read(0);

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto vlistID2 = vlistDuplicate(vlistID1);

  auto nvars = vlistNvars(vlistID1);
  for (int varID = 0; varID < nvars; ++varID) cdiDefKeyString(vlistID2, varID, CDI_KEY_UNITS, "No.");

  vlistDefNtsteps(vlistID2, (compareDate == CMP_DATE) ? 1 : -1);

  auto taxisID1 = vlistInqTaxis(vlistID1);
  auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  auto streamID2 = cdo_open_write(1);
  cdo_def_vlist(streamID2, vlistID2);

  auto maxrecs = vlistNrecs(vlistID1);
  std::vector<RecordInfo> recList(maxrecs);

  VarList varList1;
  varListInit(varList1, vlistID1);

  Field field;

  FieldVector2D vars1;
  fields_from_vlist(vlistID1, vars1, FIELD_VEC);

  int tsID = 0;
  int otsID = 0;
  while (true)
    {
      int nrecs = 0;
      int numSets = 0;
      while (true)
        {
          nrecs = cdo_stream_inq_timestep(streamID1, tsID);
          if (nrecs == 0) break;

          auto vDateTime = taxisInqVdatetime(taxisID1);

          if (numSets == 0) vDateTime0 = vDateTime;

          if (date_is_neq(vDateTime, vDateTime0, compareDate))
            {
              cdo_add_steps(-1);
              break;
            }

          for (int recID = 0; recID < nrecs; ++recID)
            {
              int varID, levelID;
              cdo_inq_record(streamID1, &varID, &levelID);

              if (tsID == 0) recList[recID].set(varID, levelID);

              const size_t fieldsize = vars1[varID][levelID].size;

              if (numSets == 0)
                {
                  for (size_t i = 0; i < fieldsize; ++i) vars1[varID][levelID].vec_d[i] = vars1[varID][levelID].missval;
                  vars1[varID][levelID].nmiss = fieldsize;
                }

              field.init(varList1[varID]);
              cdo_read_record(streamID1, field);

              field2_count(vars1[varID][levelID], field);
            }

          vDateTimeN = vDateTime;
          numSets++;
          tsID++;
        }

      if (nrecs == 0 && numSets == 0) break;

      taxisDefVdatetime(taxisID2, vDateTimeN);
      cdo_def_timestep(streamID2, otsID);

      for (int recID = 0; recID < maxrecs; ++recID)
        {
          auto [varID, levelID] = recList[recID].get();
          if (otsID && varList1[varID].isConstant) continue;

          cdo_def_record(streamID2, varID, levelID);
          cdo_write_record(streamID2, vars1[varID][levelID].vec_d.data(), vars1[varID][levelID].nmiss);
        }

      if (nrecs == 0) break;
      otsID++;
    }

  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  cdo_finish();

  return nullptr;
}
