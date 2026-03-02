/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Copyright (C) 2006 Brockmann Consult

  Author: Ralf Quast
          Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Timpctl    timpctl         Time percentiles
      Hourpctl   hourpctl        Hourly percentiles
      Daypctl    daypctl         Daily percentiles
      Monpctl    monpctl         Monthly percentiles
      Yearpctl   yearpctl        Yearly percentiles
*/

#include <cdi.h>

#include "util_date.h"
#include "process_int.h"
#include "cdo_vlist.h"
#include "param_conversion.h"
#include "percentiles_hist.h"
#include "datetime.h"
#include "field_functions.h"

static void
timpctl(int operatorID)
{
  auto timestat_date = TimeStat::MEAN;
  CdiDateTime vDateTime0{};

  operator_input_arg("percentile number");
  auto pn = parameter_to_double(cdo_operator_argv(0));

  auto compareDate = cdo_operator_f2(operatorID);

  auto streamID1 = cdo_open_read(0);
  auto streamID2 = cdo_open_read(1);
  auto streamID3 = cdo_open_read(2);

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto vlistID2 = cdo_stream_inq_vlist(streamID2);
  auto vlistID3 = cdo_stream_inq_vlist(streamID3);
  auto vlistID4 = vlistDuplicate(vlistID1);

  vlist_compare(vlistID1, vlistID2, CMP_ALL);
  vlist_compare(vlistID1, vlistID3, CMP_ALL);

  if (cdo_operator_f2(operatorID) == 16) vlistDefNtsteps(vlistID4, 1);

  auto taxisID1 = vlistInqTaxis(vlistID1);
  auto taxisID2 = vlistInqTaxis(vlistID2);
  auto taxisID3 = vlistInqTaxis(vlistID3);
  // TODO - check that time axes 2 and 3 are equal

  auto taxisID4 = taxisDuplicate(taxisID1);
  taxisWithBounds(taxisID4);
  vlistDefTaxis(vlistID4, taxisID4);

  auto streamID4 = cdo_open_write(3);
  cdo_def_vlist(streamID4, vlistID4);

  auto ntsteps = vlistNtsteps(vlistID1);
  auto nvars = vlistNvars(vlistID1);

  auto maxrecs = vlistNrecs(vlistID1);
  std::vector<RecordInfo> recList(maxrecs);

  FieldVector constFields(maxrecs);

  DateTimeList dtlist;
  dtlist.set_stat(timestat_date);
  dtlist.set_calendar(taxisInqCalendar(taxisID1));

  Field field1, field2;

  VarList varList1;
  varListInit(varList1, vlistID1);

  HistogramSet hset(nvars, ntsteps);

  for (int varID = 0; varID < nvars; ++varID) hset.createVarLevels(varID, varList1[varID].nlevels, varList1[varID].gridsize);

  int tsID = 0;
  int otsID = 0;
  while (true)
    {
      auto nrecs = cdo_stream_inq_timestep(streamID2, otsID);
      if (nrecs != cdo_stream_inq_timestep(streamID3, otsID))
        cdo_abort("Number of records at time step %d of %s and %s differ!", otsID + 1, cdo_get_stream_name(1),
                  cdo_get_stream_name(2));

      auto vDateTime2 = taxisInqVdatetime(taxisID2);
      auto vDateTime3 = taxisInqVdatetime(taxisID3);
      if (cdiDateTime_isNE(vDateTime2, vDateTime3))
        cdo_abort("Verification dates at time step %d of %s and %s differ!", otsID + 1, cdo_get_stream_name(1),
                  cdo_get_stream_name(2));

      for (int recID = 0; recID < nrecs; ++recID)
        {
          int varID, levelID;
          cdo_inq_record(streamID2, &varID, &levelID);
          field1.init(varList1[varID]);
          cdo_read_record(streamID2, field1);

          cdo_inq_record(streamID3, &varID, &levelID);
          field2.init(varList1[varID]);
          cdo_read_record(streamID3, field2);

          hset.defVarLevelBounds(varID, levelID, field1, field2);
        }

      int numSets = 0;
      while (nrecs && (nrecs = cdo_stream_inq_timestep(streamID1, tsID)))
        {
          dtlist.taxis_inq_timestep(taxisID1, numSets);
          auto vDateTime1 = dtlist.get_vDateTime(numSets);
          if (numSets == 0) vDateTime0 = vDateTime1;

          if (date_is_neq(vDateTime1, vDateTime0, compareDate))
            {
              cdo_add_steps(-1);
              break;
            }

          for (int recID = 0; recID < nrecs; ++recID)
            {
              int varID, levelID;
              cdo_inq_record(streamID1, &varID, &levelID);

              if (tsID == 0) recList[recID].set(varID, levelID);

              if (tsID == 0 && varList1[varID].isConstant)
                {
                  constFields[recID].init(varList1[varID]);
                  cdo_read_record(streamID1, constFields[recID]);
                }
              else
                {
                  field1.init(varList1[varID]);
                  cdo_read_record(streamID1, field1);

                  hset.addVarLevelValues(varID, levelID, field1);
                }
            }

          numSets++;
          tsID++;
        }

      if (nrecs == 0 && numSets == 0) break;

      dtlist.stat_taxis_def_timestep(taxisID4, numSets);
      cdo_def_timestep(streamID4, otsID);

      for (int recID = 0; recID < maxrecs; ++recID)
        {
          auto [varID, levelID] = recList[recID].get();
          if (otsID && varList1[varID].isConstant) continue;

          cdo_def_record(streamID4, varID, levelID);

          if (varList1[varID].isConstant) { cdo_write_record(streamID4, constFields[recID]); }
          else
            {
              field1.init(varList1[varID]);
              hset.getVarLevelPercentiles(field1, varID, levelID, pn);
              cdo_write_record(streamID4, field1);
            }
        }

      if (nrecs == 0) break;
      otsID++;
    }

  cdo_stream_close(streamID4);
  cdo_stream_close(streamID3);
  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);
}

void *
Timpctl(void *process)
{
  cdo_initialize(process);

  // clang-format off
  cdo_operator_add("timpctl",  FieldFunc_Pctl, CMP_DATE,  nullptr);
  cdo_operator_add("yearpctl", FieldFunc_Pctl, CMP_YEAR,  nullptr);
  cdo_operator_add("monpctl",  FieldFunc_Pctl, CMP_MONTH, nullptr);
  cdo_operator_add("daypctl",  FieldFunc_Pctl, CMP_DAY,   nullptr);
  cdo_operator_add("hourpctl", FieldFunc_Pctl, CMP_HOUR,  nullptr);
  // clang-format on

  timpctl(cdo_operator_id());

  cdo_finish();

  return nullptr;
}
