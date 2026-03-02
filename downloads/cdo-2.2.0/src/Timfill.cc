/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <cdi.h>

#include "process_int.h"
#include "cdo_vlist.h"
#include "param_conversion.h"
#include "statistic.h"
#include "cdo_options.h"
#include "datetime.h"
#include "cimdOmp.h"
#include "field_functions.h"

class ModuleTimFill
{
  DateTimeList dtlist;
  CdoStreamID streamID1;
  CdoStreamID streamID2;
  int vlistID1;
  int vlistID2;

  int taxisID1;
  int taxisID2;

  VarList varList;
  FieldVector3D vars;

  int calendar;
  int nvars;

  Varray2D<double> data2D;

  std::vector<JulianDate> julianDates;
  int nts;

public:
  void
  init(void *process)
  {

    cdo_initialize(process);

    // nn =nearest neighbour
    cdo_operator_add("setmisstonntime", 0, 0, nullptr);

    streamID1 = cdo_open_read(0);
    streamID2 = cdo_open_write(1);

    vlistID1 = cdo_stream_inq_vlist(streamID1);
    vlistID2 = vlistDuplicate(vlistID1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    calendar = taxisInqCalendar(taxisID1);

    varListInit(varList, vlistID1);

    nvars = vlistNvars(vlistID1);
  }

  void
  nearest_neighbour_with_lowest_delta(Varray<double> &data, int current_misval, int next_misval)
  {
    for (int k = current_misval + 1; k < next_misval; ++k)
      {
        // nearest_neighbour
        const auto jdelta1 = julianDate_to_seconds(julianDate_sub(julianDates[k], julianDates[current_misval]));
        const auto jdelta2 = julianDate_to_seconds(julianDate_sub(julianDates[next_misval], julianDates[k]));
        data[k] = data[(jdelta1 <= jdelta2) ? current_misval : next_misval];
      }
  }

public:
  void
  step(int varID)
  {
    auto fieldMemType = varList[varID].memType;
    const auto gridsize = varList[varID].gridsize;
    const auto missval = varList[varID].missval;
    for (int levelID = 0; levelID < varList[varID].nlevels; ++levelID)
      {
#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
        for (size_t i = 0; i < gridsize; ++i)
          {
            const auto ompthID = cdo_omp_get_thread_num();
            auto &data = data2D[ompthID];

            if (fieldMemType == MemType::Float)
              for (int t = 0; t < nts; ++t) data[t] = vars[t][varID][levelID].vec_f[i];
            else
              for (int t = 0; t < nts; ++t) data[t] = vars[t][varID][levelID].vec_d[i];

            constexpr int not_found = -1;
            int current_misval = not_found;  // first_found
            int next_misval = not_found;     // nearest_neighbour
            for (int t = 0; t < nts; ++t)
              {
                // if is_missing_val
                if (dbl_is_equal(data[t], missval))
                  {
                    // search for next non missing value
                    for (int k = t + 1; k < nts; ++k)
                      {
                        if (!dbl_is_equal(data[k], missval))
                          {
                            next_misval = k;
                            break;
                          }
                      }

                    // if first not found and then no other missingvalue
                    if (current_misval == not_found && next_misval != not_found)
                      {
                        for (int k = t; k < next_misval; ++k) data[k] = data[next_misval];
                        t = next_misval;               // advance iterator
                        current_misval = next_misval;  // current_misval is now the first ever found
                        next_misval = not_found;       // reset current is the frist ever found
                      }
                    else if (current_misval != not_found && next_misval != not_found)
                      {
                        nearest_neighbour_with_lowest_delta(data, current_misval, next_misval);
                        t = next_misval;
                        current_misval = next_misval;
                        next_misval = not_found;
                      }
                    // only one was found aka rest, take the first value
                    else if (current_misval != not_found && next_misval == not_found)
                      {
                        for (int k = t; k < nts; ++k) data[k] = data[current_misval];
                        break;
                      }
                    else if (current_misval == not_found && next_misval == not_found) { break; }
                  }
                else { current_misval = t; }
              }

            if (fieldMemType == MemType::Float)
              for (int t = 0; t < nts; ++t) vars[t][varID][levelID].vec_f[i] = data[t];
            else
              for (int t = 0; t < nts; ++t) vars[t][varID][levelID].vec_d[i] = data[t];
          }
      }
  }
  void
  run()
  {
    int tsID = 0;
    while (true)
      {
        const auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
        if (nrecs == 0) break;

        constexpr size_t NALLOC_INC = 1024;
        if ((size_t) tsID >= vars.size()) vars.resize(vars.size() + NALLOC_INC);

        dtlist.taxis_inq_timestep(taxisID1, tsID);

        fields_from_vlist(vlistID1, vars[tsID]);

        for (int recID = 0; recID < nrecs; ++recID)
          {
            int varID, levelID;
            cdo_inq_record(streamID1, &varID, &levelID);
            auto &field = vars[tsID][varID][levelID];
            field.init(varList[varID]);
            cdo_read_record(streamID1, field);
          }

        tsID++;
      }

    nts = tsID;
    if (nts <= 1) cdo_abort("Number of time steps <= 1!");

    julianDates = std::vector<JulianDate>(nts);
    for (tsID = 0; tsID < nts; ++tsID) { julianDates[tsID] = julianDate_encode(calendar, dtlist.get_vDateTime(tsID)); }

    data2D = Varray2D<double>(Threading::ompNumThreads);
    for (auto &data : data2D) data.resize(nts);

    for (int varID = 0; varID < nvars; ++varID) { step(varID); }

    cdo_def_vlist(streamID2, vlistID2);

    for (tsID = 0; tsID < nts; ++tsID)
      {
        dtlist.taxis_def_timestep(taxisID2, tsID);
        cdo_def_timestep(streamID2, tsID);

        for (int varID = 0; varID < nvars; ++varID)
          {
            for (int levelID = 0; levelID < varList[varID].nlevels; ++levelID)
              {
                auto &field = vars[tsID][varID][levelID];
                if (field.hasData())
                  {
                    cdo_def_record(streamID2, varID, levelID);
                    cdo_write_record(streamID2, field);
                  }
              }
          }
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
Timfill(void *process)
{
  ModuleTimFill timfill;
  timfill.init(process);
  timfill.run();
  timfill.close();

  return nullptr;
}
