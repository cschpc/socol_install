/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Copyright (C) 2006 Brockmann Consult

  Author: Ralf Quast
          Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Ymonpctl   ymonpctl        Multi-year monthly percentiles
*/

#include <cdi.h>

#include "cdo_options.h"
#include "cdo_vlist.h"
#include "field.h"
#include "datetime.h"
#include "process_int.h"
#include "param_conversion.h"
#include "percentiles_hist.h"
#include "field_functions.h"

class ModuleYmonpctl
{
  constexpr static int MaxMonths = 17;

  CdoStreamID streamID1;
  CdoStreamID streamID2;
  CdoStreamID streamID3;
  CdoStreamID streamID4;

  int taxisID1;
  int taxisID2;
  int taxisID3;
  int taxisID4;

  CdiDateTime vDateTimes1[MaxMonths]{};
  CdiDateTime vDateTimes2[MaxMonths]{};

  std::vector<bool> vars1 = std::vector<bool>(MaxMonths, false);
  HistogramSet hsets[MaxMonths];
  long numSets[MaxMonths] = { 0 };
  VarList varList1;
  Field field1, field2;
  FieldVector constFields;
  std::vector<RecordInfo> recList;

  int maxrecs;
  int nvars;
  int ntsteps;

  double pn;

public:
  void
  init(void *process)
  {
    cdo_initialize(process);
    cdo_operator_add("ymonpctl", FieldFunc_Pctl, 0, nullptr);

    operator_input_arg("percentile number");
    pn = parameter_to_double(cdo_operator_argv(0));

    streamID1 = cdo_open_read(0);
    streamID2 = cdo_open_read(1);
    streamID3 = cdo_open_read(2);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = cdo_stream_inq_vlist(streamID2);
    auto vlistID3 = cdo_stream_inq_vlist(streamID3);
    auto vlistID4 = vlistDuplicate(vlistID1);

    vlist_compare(vlistID1, vlistID2, CMP_ALL);
    vlist_compare(vlistID1, vlistID3, CMP_ALL);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = vlistInqTaxis(vlistID2);
    taxisID3 = vlistInqTaxis(vlistID3);
    // TODO - check that time axes 2 and 3 are equal

    taxisID4 = taxisDuplicate(taxisID1);
    if (taxisHasBounds(taxisID4)) taxisDeleteBounds(taxisID4);
    vlistDefTaxis(vlistID4, taxisID4);

    streamID4 = cdo_open_write(3);
    cdo_def_vlist(streamID4, vlistID4);

    ntsteps = vlistNtsteps(vlistID1);
    nvars = vlistNvars(vlistID1);

    maxrecs = vlistNrecs(vlistID1);
    recList = std::vector<RecordInfo>(maxrecs);

    constFields = FieldVector(maxrecs);

    varListInit(varList1, vlistID1);
  }
  void
  run()
  {
    int tsID = 0;
    while (true)
      {
        auto nrecs = cdo_stream_inq_timestep(streamID2, tsID);
        if (nrecs == 0) break;

        if (nrecs != cdo_stream_inq_timestep(streamID3, tsID))
          cdo_abort("Number of records at time step %d of %s and %s differ!", tsID + 1, cdo_get_stream_name(1),
                    cdo_get_stream_name(2));

        auto vDateTime2 = taxisInqVdatetime(taxisID2);
        auto vDateTime3 = taxisInqVdatetime(taxisID3);

        if (cdiDate_get(vDateTime2.date) != cdiDate_get(vDateTime3.date))
          cdo_abort("Verification dates at time step %d of %s and %s differ!", tsID + 1, cdo_get_stream_name(1),
                    cdo_get_stream_name(2));

        // if (Options::cdoVerbose) cdo_print("process timestep: %d %s", tsID + 1, datetime_to_string(vDateTime2));

        auto month = decode_month(vDateTime2.date);
        if (month < 0 || month >= MaxMonths) cdo_abort("Month %d out of range!", month);

        vDateTimes2[month] = vDateTime2;

        if (!vars1[month])
          {
            vars1[month] = true;
            hsets[month].create(nvars, ntsteps);
            for (int varID = 0; varID < nvars; ++varID)
              {
                const auto &var = varList1[varID];
                hsets[month].createVarLevels(varID, var.nlevels, var.gridsize);
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

            hsets[month].defVarLevelBounds(varID, levelID, field1, field2);
          }

        tsID++;
      }

    tsID = 0;
    while (true)
      {
        auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
        if (nrecs == 0) break;

        auto vDateTime = taxisInqVdatetime(taxisID1);

        // if (Options::cdoVerbose) cdo_print("process timestep: %d %s", tsID + 1, datetime_to_string(vDateTime));

        auto month = decode_month(vDateTime.date);
        if (month < 0 || month >= MaxMonths) cdo_abort("Month %d out of range!", month);

        vDateTimes1[month] = vDateTime;

        if (!vars1[month]) cdo_abort("No data for month %d in %s and %s", month, cdo_get_stream_name(1), cdo_get_stream_name(2));

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

                hsets[month].addVarLevelValues(varID, levelID, field1);
              }
          }

        numSets[month]++;
        tsID++;
      }

    int otsID = 0;
    for (int month = 0; month < MaxMonths; ++month)
      if (numSets[month])
        {
          if (decode_month(vDateTimes1[month].date) != decode_month(vDateTimes2[month].date))
            cdo_abort("Verification dates for the month %d of %s and %s are different!", month, cdo_get_stream_name(0),
                      cdo_get_stream_name(1));

          taxisDefVdatetime(taxisID4, vDateTimes1[month]);
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
                  hsets[month].getVarLevelPercentiles(field1, varID, levelID, pn);
                  cdo_write_record(streamID4, field1);
                }
            }

          otsID++;
        }
  }
  void
  close()
  {
    cdo_stream_close(streamID4);
    cdo_stream_close(streamID3);
    cdo_stream_close(streamID2);
    cdo_stream_close(streamID1);

    cdo_finish();
  }
};
void *
Ymonpctl(void *process)
{
  ModuleYmonpctl ymonpctl;
  ymonpctl.init(process);
  ymonpctl.run();
  ymonpctl.close();

  return nullptr;
}
