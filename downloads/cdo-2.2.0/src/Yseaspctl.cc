/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Copyright (C) 2006 Brockmann Consult

  Author: Ralf Quast
          Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Yseaspctl  yseaspctl       Multi-year seasonal percentiles
*/

#include <cdi.h>

#include "cdo_options.h"
#include "cdo_vlist.h"
#include "cdo_season.h"
#include "datetime.h"
#include "process_int.h"
#include "param_conversion.h"
#include "percentiles_hist.h"
#include "printinfo.h"
#include "field_functions.h"

class ModuleYseaspctl
{
  CdoStreamID streamID1;
  CdoStreamID streamID2;
  CdoStreamID streamID3;
  CdoStreamID streamID4;
  FieldVector constFields;

  int taxisID1;
  int taxisID2;
  int taxisID3;
  int taxisID4;

  VarList varList1;
  Field field1, field2;

  static const int MaxSeasons = 4;
  long numSets[MaxSeasons] = { 0 };
  CdiDateTime vDateTimes1[MaxSeasons]{};
  CdiDateTime vDateTimes2[MaxSeasons]{};

  std::vector<bool> vars1;
  HistogramSet hsets[MaxSeasons];
  std::vector<RecordInfo> recList;

  int maxrecs;
  int ntsteps;
  int nvars;
  int pn;

public:
  void
  init(void *process)
  {
    vars1 = std::vector<bool>(MaxSeasons, false);

    cdo_initialize(process);

    cdo_operator_add("yseaspctl", FieldFunc_Pctl, 0, nullptr);

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

        auto vDateTime = taxisInqVdatetime(taxisID2);

        if (cdiDateTime_isNE(vDateTime, taxisInqVdatetime(taxisID3)))
          cdo_abort("Verification dates at time step %d of %s and %s differ!", tsID + 1, cdo_get_stream_name(1),
                    cdo_get_stream_name(2));

        if (Options::cdoVerbose) cdo_print("process timestep: %d %s", tsID + 1, datetime_to_string(vDateTime));

        auto seas = month_to_season(decode_month(vDateTime.date));

        set_date_time(vDateTimes2[seas], vDateTime);

        if (!vars1[seas])
          {
            vars1[seas] = true;
            hsets[seas].create(nvars, ntsteps);
            for (int varID = 0; varID < nvars; ++varID)
              {
                const auto &var = varList1[varID];
                hsets[seas].createVarLevels(varID, var.nlevels, var.gridsize);
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

            hsets[seas].defVarLevelBounds(varID, levelID, field1, field2);
          }

        tsID++;
      }

    tsID = 0;
    while (true)
      {
        auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
        if (nrecs == 0) break;

        auto vDateTime = taxisInqVdatetime(taxisID1);

        auto month = decode_month(vDateTime.date);
        if (month < 0 || month > 16) cdo_abort("Month %d out of range!", month);

        auto seas = month_to_season(month);

        set_date_time(vDateTimes1[seas], vDateTime);

        if (!vars1[seas]) cdo_abort("No data for season %d in %s and %s", seas, cdo_get_stream_name(1), cdo_get_stream_name(2));

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

                hsets[seas].addVarLevelValues(varID, levelID, field1);
              }
          }

        numSets[seas]++;
        tsID++;
      }

    int otsID = 0;
    for (int seas = 0; seas < MaxSeasons; ++seas)
      if (numSets[seas])
        {
          if (decode_month_and_day(vDateTimes1[seas].date) != decode_month_and_day(vDateTimes2[seas].date))
            cdo_abort("Verification dates for the season %d of %s and %s are different!", seas + 1, cdo_get_stream_name(0),
                      cdo_get_stream_name(1));

          taxisDefVdatetime(taxisID4, vDateTimes1[seas]);
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
                  hsets[seas].getVarLevelPercentiles(field1, varID, levelID, pn);
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
Yseaspctl(void *process)
{
  ModuleYseaspctl yseaspctl;
  yseaspctl.init(process);
  yseaspctl.run();
  yseaspctl.close();
  return nullptr;
}
