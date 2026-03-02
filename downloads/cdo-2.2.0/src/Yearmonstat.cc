/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Yearmonstat   yearmonmean        Yearly mean from monthly data
      Yearmonstat   yearmonavg         Yearly average from monthly data
*/

#include "cdi.h"
#include "calendar.h"

#include "cdo_options.h"
#include "cdo_vlist.h"
#include "process_int.h"
#include "datetime.h"
#include "printinfo.h"
#include "field_functions.h"

class YearMonStat
{

private:
  CdoStreamID streamID1;
  CdoStreamID streamID2;

  int vlistID1;
  int vlistID2;

  int taxisID1;
  int taxisID2;

  TimeStat timestat_date = TimeStat::MEAN;

  int year0 = 0, month0 = 0;
  int year, month, day;

  DateTimeList dtlist;

  VarList varList;
  Field field;

  FieldVector2D samp1, vars1;
  CdiDateTime vDateTime0{};

  int calendar;

  std::vector<RecordInfo> recList;
  int operfunc;

  int nvars;
  int maxrecs;

public:
  void
  init(void *process)
  {

    cdo_initialize(process);

    // clang-format off
    cdo_operator_add("yearmonmean",  FieldFunc_Mean, 0, nullptr);
    cdo_operator_add("yearmonavg",   FieldFunc_Avg,  0, nullptr);
    // clang-format on

    auto operatorID = cdo_operator_id();

    operfunc = cdo_operator_f1(operatorID);

    operator_check_argc(0);

    streamID1 = cdo_open_read(0);

    vlistID1 = cdo_stream_inq_vlist(streamID1);
    vlistID2 = vlistDuplicate(vlistID1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    taxisWithBounds(taxisID2);
    vlistDefTaxis(vlistID2, taxisID2);

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);

    nvars = vlistNvars(vlistID1);

    maxrecs = vlistNrecs(vlistID1);
    recList = std::vector<RecordInfo>(maxrecs);

    calendar = taxisInqCalendar(taxisID1);
    dtlist.set_stat(timestat_date);
    dtlist.set_calendar(calendar);

    varListInit(varList, vlistID1);

    fields_from_vlist(vlistID1, samp1);
    fields_from_vlist(vlistID1, vars1, FIELD_VEC);
  }

  void
  run()
  {
    int tsID = 0;
    int otsID = 0;
    while (true)
      {
        int nrecs = 0;
        long numSets = 0;
        double dsets = 0.0;
        while (true)
          {
            nrecs = cdo_stream_inq_timestep(streamID1, tsID);
            if (nrecs == 0) break;

            dtlist.taxis_inq_timestep(taxisID1, numSets);
            auto vDateTime = dtlist.get_vDateTime(numSets);
            cdiDate_decode(vDateTime.date, &year, &month, &day);

            if (numSets == 0) year0 = year;

            if (year != year0)
              {
                cdo_add_steps(-1);
                break;
              }

            if (numSets > 0 && month == month0)
              {
                cdo_warning("   last timestep: %s", datetime_to_string(vDateTime0));
                cdo_warning("current timestep: %s", datetime_to_string(vDateTime));
                cdo_abort("Month does not change!");
              }

            auto dpm = days_per_month(calendar, year, month);

            for (int recID = 0; recID < nrecs; ++recID)
              {
                int varID, levelID;
                cdo_inq_record(streamID1, &varID, &levelID);
                const auto &var = varList[varID];

                if (tsID == 0) recList[recID].set(varID, levelID);

                auto &rsamp1 = samp1[varID][levelID];
                auto &rvars1 = vars1[varID][levelID];

                if (numSets == 0)
                  {
                    cdo_read_record(streamID1, rvars1);

                    fieldc_mul(rvars1, dpm);

                    if (rvars1.nmiss || !rsamp1.empty())
                      {
                        if (rsamp1.empty()) rsamp1.resize(rvars1.size);
                        field2_vinit(rsamp1, rvars1, dpm);
                      }
                  }
                else
                  {
                    field.init(var);
                    cdo_read_record(streamID1, field);

                    fieldc_mul(field, dpm);

                    if (field.nmiss || !rsamp1.empty())
                      {
                        if (rsamp1.empty()) rsamp1.resize(rvars1.size, dsets);
                        field2_vincr(rsamp1, field, dpm);
                      }

                    field2_function(rvars1, field, operfunc);
                  }
              }

            month0 = month;
            vDateTime0 = vDateTime;
            numSets++;
            dsets += dpm;
            tsID++;
          }

        if (nrecs == 0 && numSets == 0) break;

        for (int varID = 0; varID < nvars; ++varID)
          {
            const auto &var = varList[varID];
            if (var.isConstant) continue;
            for (int levelID = 0; levelID < var.nlevels; ++levelID)
              {
                const auto &rsamp1 = samp1[varID][levelID];
                auto &rvars1 = vars1[varID][levelID];
                if (!rsamp1.empty())
                  field2_div(rvars1, rsamp1);
                else
                  fieldc_div(rvars1, dsets);
              }
          }

        if (Options::cdoVerbose) cdo_print("%s  numSets = %ld", datetime_to_string(vDateTime0), numSets);

        dtlist.stat_taxis_def_timestep(taxisID2, numSets);
        cdo_def_timestep(streamID2, otsID);

        for (int recID = 0; recID < maxrecs; ++recID)
          {
            auto [varID, levelID] = recList[recID].get();
            if (otsID && varList[varID].isConstant) continue;

            auto &rvars1 = vars1[varID][levelID];
            cdo_def_record(streamID2, varID, levelID);
            cdo_write_record(streamID2, rvars1);
          }

        if (nrecs == 0) break;
        otsID++;
      }
  }

  void
  close()
  {
    cdo_stream_close(streamID2);
    cdo_stream_close(streamID1);

    cdo_finish();
  }
};

void *
Yearmonstat(void *process)
{
  YearMonStat yearmonstat;
  yearmonstat.init(process);

  yearmonstat.run();
  yearmonstat.close();

  return nullptr;
}
