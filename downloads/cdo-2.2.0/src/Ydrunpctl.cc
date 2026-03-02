/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Copyright (C) 2006 Brockmann Consult

  Author: Ralf Quast
          Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Ydrunpctl    ydrunpctl         Multi-year daily running percentiles
*/

#include <cdi.h>
#include "calendar.h"

#include "cdo_options.h"
#include "cdo_vlist.h"
#include "datetime.h"
#include "process_int.h"
#include "param_conversion.h"
#include "percentiles_hist.h"
#include "percentiles.h"
#include "pmlist.h"
#include "field_functions.h"

constexpr int MaxDays = 373;
class ModuleYdrunpctl
{
private:
  CdoStreamID streamID1;
  CdoStreamID streamID2;
  CdoStreamID streamID3;
  CdoStreamID streamID4;

  int taxisID1;
  int taxisID2;
  int taxisID3;
  int taxisID4;

  int nvars;
  int ntsteps;
  int ndates;
  int dpy;
  double pn;
  int maxrecs;

  std::string percMethod, readMethod;

  CdiDateTime vDateTimes1[MaxDays]{};
  CdiDateTime vDateTimes2[MaxDays]{};
  HistogramSet hsets[MaxDays];
  int numSets[MaxDays] = { 0 };

  std::vector<bool> vars2;

  Field field1, field2;

  FieldVector constFields;
  VarList varList1;
  FieldVector3D vars1;
  std::vector<CdiDateTime> cdiDateTimes;
  std::vector<RecordInfo> recList;

public:
  void
  init(void *process)
  {
    vars2 = std::vector<bool>(MaxDays, false);

    cdo_initialize(process);
    cdo_operator_add("ydrunpctl", FieldFunc_Pctl, 0, nullptr);

    pn = parameter_to_double(cdo_operator_argv(0));
    ndates = parameter_to_int(cdo_operator_argv(1));

    if (cdo_operator_argc() > 2)
      {
        auto params = cdo_get_oper_argv();
        params = std::vector<std::string>(params.begin() + 2, params.end());
        KVList kvlist;
        if (kvlist.parse_arguments(cdo_operator_argc() - 2, params) != 0) cdo_abort("Argument parse error!");
        auto kv = kvlist.search("pm");
        if (kv && kv->nvalues > 0) percMethod = parameter_to_word(kv->values[0]);
        kv = kvlist.search("rm");
        if (kv && kv->nvalues > 0) readMethod = parameter_to_word(kv->values[0]);
      }

    if (percMethod == "r8") percentile_set_method("rtype8");

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

    dpy = calendar_dpy(taxisInqCalendar(taxisID1));

    streamID4 = cdo_open_write(3);
    cdo_def_vlist(streamID4, vlistID4);

    ntsteps = vlistNtsteps(vlistID1);
    nvars = vlistNvars(vlistID1);

    maxrecs = vlistNrecs(vlistID1);
    recList = std::vector<RecordInfo>(maxrecs);

    constFields = FieldVector(maxrecs);

    varListInit(varList1, vlistID1);

    cdiDateTimes = std::vector<CdiDateTime>(ndates + 1);

    vars1 = FieldVector3D(ndates + 1);
    for (int its = 0; its < ndates; its++) fields_from_vlist(vlistID1, vars1[its], FIELD_VEC | FIELD_NAT);
  }
  void
  run()
  {
    int startYear = 0;
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

        // if (Options::cdoVerbose) cdo_print("process timestep: %d %d %d", tsID + 1, vdate, vtime);

        auto dayOfYear = decode_day_of_year(vDateTime.date);
        if (dayOfYear < 0 || dayOfYear >= MaxDays) cdo_abort("Day %d out of range!", dayOfYear);

        vDateTimes2[dayOfYear] = vDateTime;

        if (!vars2[dayOfYear])
          {
            vars2[dayOfYear] = true;
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

    for (tsID = 0; tsID < ndates; ++tsID)
      {
        auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
        if (nrecs == 0) cdo_abort("File has less then %d timesteps!", ndates);

        cdiDateTimes[tsID] = taxisInqVdatetime(taxisID1);

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
            else { cdo_read_record(streamID1, vars1[tsID][varID][levelID]); }
          }
      }

    while (true)
      {
        cdiDateTimes[ndates] = datetime_avg(dpy, ndates, cdiDateTimes);

        auto vDateTime = cdiDateTimes[ndates];

        auto dayOfYear = decode_day_of_year(vDateTime.date);
        if (dayOfYear < 0 || dayOfYear >= MaxDays) cdo_abort("Day %d out of range!", dayOfYear);

        vDateTimes1[dayOfYear] = vDateTime;

        if (!vars2[dayOfYear])
          cdo_abort("No data for day %d in %s and %s", dayOfYear, cdo_get_stream_name(1), cdo_get_stream_name(2));

        for (int varID = 0; varID < nvars; ++varID)
          {
            const auto &var = varList1[varID];
            if (var.isConstant) continue;

            for (int levelID = 0; levelID < var.nlevels; ++levelID)
              for (int inp = 0; inp < ndates; ++inp) hsets[dayOfYear].addVarLevelValues(varID, levelID, vars1[inp][varID][levelID]);
          }

        cdiDateTimes[ndates] = cdiDateTimes[0];
        vars1[ndates] = vars1[0];

        for (int inp = 0; inp < ndates; ++inp)
          {
            cdiDateTimes[inp] = cdiDateTimes[inp + 1];
            vars1[inp] = vars1[inp + 1];
          }

        auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
        if (nrecs == 0) break;

        cdiDateTimes[ndates - 1] = taxisInqVdatetime(taxisID1);

        for (int recID = 0; recID < nrecs; ++recID)
          {
            int varID, levelID;
            cdo_inq_record(streamID1, &varID, &levelID);
            cdo_read_record(streamID1, vars1[ndates - 1][varID][levelID]);
          }

        numSets[dayOfYear] += ndates;
        tsID++;
      }

    if (readMethod == "c" && cdo_assert_files_only())
      {
        auto endYear = cdiDateTimes[ndates - 1].date.year;
        auto cdiStream = streamOpenRead(cdo_get_stream_name(0));
        auto cdiVlistID = streamInqVlist(cdiStream);
        auto cdiTaxisID = vlistInqTaxis(cdiVlistID);
        int missTimes = 0;
        for (missTimes = 0; missTimes < ndates - 1; missTimes++)
          {
            auto nrecs = streamInqTimestep(cdiStream, missTimes);
            if (nrecs == 0) break;

            cdiDateTimes[ndates - 1] = taxisInqVdatetime(cdiTaxisID);
            cdiDateTimes[ndates - 1].date.year = endYear + 1;

            for (int recID = 0; recID < nrecs; ++recID)
              {
                int varID, levelID;
                streamInqRecord(cdiStream, &varID, &levelID);
                auto &pvars1 = vars1[ndates - 1][varID][levelID];
                if (pvars1.memType == MemType::Float)
                  streamReadRecordF(cdiStream, pvars1.vec_f.data(), &pvars1.nmiss);
                else
                  streamReadRecord(cdiStream, pvars1.vec_d.data(), &pvars1.nmiss);
              }

            cdiDateTimes[ndates] = datetime_avg(dpy, ndates, cdiDateTimes);
            auto vDateTime = cdiDateTimes[ndates];
            if (vDateTime.date.year > endYear) vDateTime.date.year = startYear;

            auto dayOfYear = decode_day_of_year(vDateTime.date);
            if (dayOfYear < 0 || dayOfYear >= MaxDays) cdo_abort("Day %d out of range!", dayOfYear);

            numSets[dayOfYear] += ndates;

            for (int varID = 0; varID < nvars; ++varID)
              {
                const auto &var = varList1[varID];
                if (var.isConstant) continue;

                for (int levelID = 0; levelID < var.nlevels; ++levelID)
                  for (int inp = 0; inp < ndates; ++inp)
                    hsets[dayOfYear].addVarLevelValues(varID, levelID, vars1[inp][varID][levelID]);
              }

            cdiDateTimes[ndates] = cdiDateTimes[0];
            vars1[ndates] = vars1[0];

            for (int inp = 0; inp < ndates; ++inp)
              {
                cdiDateTimes[inp] = cdiDateTimes[inp + 1];
                vars1[inp] = vars1[inp + 1];
              }
          }

        if (missTimes != ndates - 1) cdo_abort("Addding the missing values when using the 'readMethod' method was not possible");

        streamClose(cdiStream);
      }
    else if (readMethod == "c")
      cdo_warning("Operators cannot be piped in circular mode");

    /*
    int outyear = 1e9;
    for (dayOfYear = 0; dayOfYear < MaxDays; dayOfYear++)
      if (numSets[dayOfYear])
        {
          int year, month, day;
          cdiDate_decode(vDateTimes1[dayOfYear].date, &year, &month, &day);
          if (year < outyear) outyear = year;
        }

    for (dayOfYear = 0; dayOfYear < MaxDays; dayOfYear++)
      if (numSets[dayOfYear])
        {
          int year, month, day;
          cdiDate_decode(vDateTimes1[dayOfYear].date, &year, &month, &day);
          vDateTimes1[dayOfYear].date = cdiDate_encode(outyear, month, day);
        }
    */
    int otsID = 0;
    for (int dayOfYear = 0; dayOfYear < MaxDays; dayOfYear++)
      if (numSets[dayOfYear])
        {
          if (decode_month_and_day(vDateTimes1[dayOfYear].date) != decode_month_and_day(vDateTimes2[dayOfYear].date))
            cdo_abort("Verification dates for day %d of %s, %s and %s are different!", dayOfYear, cdo_get_stream_name(0),
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
Ydrunpctl(void *process)
{
  ModuleYdrunpctl ydrunstat;
  ydrunstat.init(process);
  ydrunstat.run();
  ydrunstat.close();

  return nullptr;
}
