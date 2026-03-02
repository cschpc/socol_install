/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Copyright (C) 2006 Brockmann Consult

  Author: Ralf Quast
          Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Ydaypctl   ydaypctl        Multi-year daily percentiles
*/

#include <cdi.h>

#include "cdo_options.h"
#include "cdo_vlist.h"
#include "datetime.h"
#include "process_int.h"
#include "param_conversion.h"
#include "percentiles_hist.h"
#include "field_functions.h"

void *
Ydaypctl(void *process)
{
  constexpr int MaxDays = 373;
  CdiDateTime vDateTimes1[MaxDays]{};
  CdiDateTime vDateTimes2[MaxDays]{};
  long numSets[MaxDays] = { 0 };
  std::vector<bool> vars1(MaxDays, false);
  HistogramSet hsets[MaxDays];

  cdo_initialize(process);
  cdo_operator_add("ydaypctl", FieldFunc_Pctl, 0, nullptr);

  operator_input_arg("percentile number");
  const auto pn = parameter_to_double(cdo_operator_argv(0));

  const auto streamID1 = cdo_open_read(0);
  const auto streamID2 = cdo_open_read(1);
  const auto streamID3 = cdo_open_read(2);

  const auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  const auto vlistID2 = cdo_stream_inq_vlist(streamID2);
  const auto vlistID3 = cdo_stream_inq_vlist(streamID3);
  const auto vlistID4 = vlistDuplicate(vlistID1);

  vlist_compare(vlistID1, vlistID2, CMP_ALL);
  vlist_compare(vlistID1, vlistID3, CMP_ALL);

  const auto taxisID1 = vlistInqTaxis(vlistID1);
  const auto taxisID2 = vlistInqTaxis(vlistID2);
  const auto taxisID3 = vlistInqTaxis(vlistID3);
  // TODO - check that time axes 2 and 3 are equal

  const auto taxisID4 = taxisDuplicate(taxisID1);
  if (taxisHasBounds(taxisID4)) taxisDeleteBounds(taxisID4);
  vlistDefTaxis(vlistID4, taxisID4);

  const auto streamID4 = cdo_open_write(3);
  cdo_def_vlist(streamID4, vlistID4);

  const auto ntsteps = vlistNtsteps(vlistID1);
  const auto nvars = vlistNvars(vlistID1);

  const auto maxrecs = vlistNrecs(vlistID1);
  std::vector<RecordInfo> recList(maxrecs);

  FieldVector constFields(maxrecs);

  Field field1, field2;

  VarList varList1;
  varListInit(varList1, vlistID1);

  int tsID = 0;
  while (true)
    {
      const auto nrecs = cdo_stream_inq_timestep(streamID2, tsID);
      if (nrecs == 0) break;

      if (nrecs != cdo_stream_inq_timestep(streamID3, tsID))
        cdo_abort("Number of records at time step %d of %s and %s differ!", tsID + 1, cdo_get_stream_name(1),
                  cdo_get_stream_name(2));

      const auto vDateTime = taxisInqVdatetime(taxisID2);

      if (cdiDateTime_isNE(vDateTime, taxisInqVdatetime(taxisID3)))
        cdo_abort("Verification dates at time step %d of %s and %s differ!", tsID + 1, cdo_get_stream_name(1),
                  cdo_get_stream_name(2));

      // if (Options::cdoVerbose) cdo_print("process timestep: %d %s", tsID + 1, datetime_to_string(vDateTime));

      const auto dayOfYear = decode_day_of_year(vDateTime.date);
      if (dayOfYear < 0 || dayOfYear >= MaxDays) cdo_abort("Day %d out of range!", dayOfYear);

      vDateTimes2[dayOfYear] = vDateTime;

      if (!vars1[dayOfYear])
        {
          vars1[dayOfYear] = true;
          hsets[dayOfYear].create(nvars, ntsteps);

          for (int varID = 0; varID < nvars; ++varID)
            {
              const auto &var = varList1[varID];
              hsets[dayOfYear].createVarLevels(varID, var.nlevels, var.gridsize);
            }
        }

      for (int recID = 0; recID < nrecs; ++recID)
        {
          int varID, levelID;
          cdo_inq_record(streamID2, &varID, &levelID);
          const auto &var = varList1[varID];
          field1.init(var);
          cdo_read_record(streamID2, field1);

          cdo_inq_record(streamID3, &varID, &levelID);
          field2.init(var);
          cdo_read_record(streamID3, field2);

          hsets[dayOfYear].defVarLevelBounds(varID, levelID, field1, field2);
        }

      tsID++;
    }

  tsID = 0;
  while (true)
    {
      const auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
      if (nrecs == 0) break;

      const auto vDateTime = taxisInqVdatetime(taxisID1);

      // if (Options::cdoVerbose) cdo_print("process timestep: %d %s", tsID + 1, datetime_to_string(vDateTime));

      const auto dayOfYear = decode_day_of_year(vDateTime.date);
      if (dayOfYear < 0 || dayOfYear >= MaxDays) cdo_abort("Day %d out of range!", dayOfYear);

      vDateTimes1[dayOfYear] = vDateTime;

      if (!vars1[dayOfYear])
        cdo_abort("No data for day %d in %s and %s", dayOfYear, cdo_get_stream_name(1), cdo_get_stream_name(2));

      for (int recID = 0; recID < nrecs; ++recID)
        {
          int varID, levelID;
          cdo_inq_record(streamID1, &varID, &levelID);
          const auto &var = varList1[varID];

          if (tsID == 0) recList[recID].set(varID, levelID);

          if (tsID == 0 && var.isConstant)
            {
              constFields[recID].init(var);
              cdo_read_record(streamID1, constFields[recID]);
            }
          else
            {
              field1.init(var);
              cdo_read_record(streamID1, field1);

              hsets[dayOfYear].addVarLevelValues(varID, levelID, field1);
            }
        }

      numSets[dayOfYear]++;
      tsID++;
    }

  int otsID = 0;
  for (int dayOfYear = 0; dayOfYear < MaxDays; dayOfYear++)
    if (numSets[dayOfYear])
      {
        if (decode_month_and_day(vDateTimes1[dayOfYear].date) != decode_month_and_day(vDateTimes2[dayOfYear].date))
          cdo_abort("Verification dates for the day %d of %s and %s are different!", dayOfYear, cdo_get_stream_name(0),
                    cdo_get_stream_name(1));

        taxisDefVdatetime(taxisID4, vDateTimes1[dayOfYear]);
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
                hsets[dayOfYear].getVarLevelPercentiles(field1, varID, levelID, pn);
                cdo_write_record(streamID4, field1);
              }
          }

        otsID++;
      }

  cdo_stream_close(streamID4);
  cdo_stream_close(streamID3);
  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  cdo_finish();

  return nullptr;
}
