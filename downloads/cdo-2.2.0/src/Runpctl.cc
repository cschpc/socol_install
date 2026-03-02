/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida
          Ralf Quast

*/

/*
   This module contains the following operators:

      Runpctl    runpctl         Running percentiles
*/

#include <cdi.h>

#include "process_int.h"
#include "param_conversion.h"
#include "percentiles.h"
#include "datetime.h"
#include "field_functions.h"
#include "cimdOmp.h"

template <typename T>
static size_t
runpctl(double pn, int ndates, size_t gridsize, Varray<T> &v2, T missval, const FieldVector3D &vars1, int varID, int levelID,
        MemType memType)
{
  size_t nmiss = 0;
  Varray2D<T> array_2D(Threading::ompNumThreads, Varray<T>(ndates));

#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(dynamic)
#endif
  for (size_t i = 0; i < gridsize; ++i)
    {
      const auto ompthID = cdo_omp_get_thread_num();
      auto &array = array_2D[ompthID];

      int j = 0;
      if (memType == MemType::Float)
        {
          for (int inp = 0; inp < ndates; ++inp)
            {
              const auto val = vars1[inp][varID][levelID].vec_f[i];
              if (!dbl_is_equal(val, missval)) array[j++] = val;
            }
        }
      else
        {
          for (int inp = 0; inp < ndates; ++inp)
            {
              const auto val = vars1[inp][varID][levelID].vec_d[i];
              if (!dbl_is_equal(val, missval)) array[j++] = val;
            }
        }

      if (j > 0) { v2[i] = percentile(array.data(), j, pn); }
      else
        {
          v2[i] = missval;
          nmiss++;
        }
    }

  return nmiss;
}

static void
runpctl(double pn, int ndates, Field &field1, const FieldVector3D &vars1, int varID, int levelID)
{
  if (field1.memType == MemType::Float)
    field1.nmiss
        = runpctl(pn, ndates, field1.gridsize, field1.vec_f, (float) field1.missval, vars1, varID, levelID, field1.memType);
  else
    field1.nmiss = runpctl(pn, ndates, field1.gridsize, field1.vec_d, field1.missval, vars1, varID, levelID, field1.memType);
}

class RunPctl
{
  CdoStreamID streamID1;
  CdoStreamID streamID2;
  int vlistID1;
  int vlistID2;
  VarList varList1;
  int taxisID1;
  int taxisID2;
  FieldVector3D vars1;

  DateTimeList dtlist;
  double pn;
  int ndates;
  int nvars;
  int maxrecs;
  std::vector<RecordInfo> recList;
  int tsID;

public:
  void
  init(void *process)
  {

    const auto timestat_date = TimeStat::MEAN;

    cdo_initialize(process);

    operator_input_arg("percentile number, number of timesteps");
    operator_check_argc(2);
    pn = parameter_to_double(cdo_operator_argv(0));
    ndates = parameter_to_int(cdo_operator_argv(1));

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

    dtlist.set_stat(timestat_date);
    dtlist.set_calendar(taxisInqCalendar(taxisID1));

    Varray<float> array_f(ndates);
    Varray<double> array_d(ndates);
    vars1 = FieldVector3D(ndates + 1);

    for (int its = 0; its < ndates; its++) fields_from_vlist(vlistID1, vars1[its]);

    varListInit(varList1, vlistID1);
  }
  void
  write_records(const int otsID)
  {

    dtlist.stat_taxis_def_timestep(taxisID2, ndates);
    cdo_def_timestep(streamID2, otsID);
    for (int recID = 0; recID < maxrecs; ++recID)
      {
        auto [varID, levelID] = recList[recID].get();
        if (otsID && varList1[varID].isConstant) continue;

        cdo_def_record(streamID2, varID, levelID);
        auto &field1 = vars1[0][varID][levelID];
        cdo_write_record(streamID2, field1);
      }
  }
  void
  run()
  {

    for (tsID = 0; tsID < ndates; ++tsID)
      {
        auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
        if (nrecs == 0) cdo_abort("File has less than %d timesteps!", ndates);

        dtlist.taxis_inq_timestep(taxisID1, tsID);

        for (int recID = 0; recID < nrecs; ++recID)
          {
            int varID, levelID;
            cdo_inq_record(streamID1, &varID, &levelID);

            if (tsID == 0) recList[recID].set(varID, levelID);

            auto &field = vars1[tsID][varID][levelID];
            field.init(varList1[varID]);
            cdo_read_record(streamID1, field);
          }
      }
    int otsID = 0;
    while (true)
      {
        for (int varID = 0; varID < nvars; ++varID)
          {
            if (varList1[varID].isConstant) continue;

            const auto nlevels = varList1[varID].nlevels;
            for (int levelID = 0; levelID < nlevels; ++levelID)
              {
                auto &field1 = vars1[0][varID][levelID];
                runpctl(pn, ndates, field1, vars1, varID, levelID);
              }
          }

        write_records(otsID);
        otsID++;

        dtlist.shift();

        vars1[ndates] = vars1[0];
        for (int inp = 0; inp < ndates; ++inp) vars1[inp] = vars1[inp + 1];

        const auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
        if (nrecs == 0) break;

        dtlist.taxis_inq_timestep(taxisID1, ndates - 1);

        for (int recID = 0; recID < nrecs; ++recID)
          {
            int varID, levelID;
            cdo_inq_record(streamID1, &varID, &levelID);
            auto &fieldN = vars1[ndates - 1][varID][levelID];
            cdo_read_record(streamID1, fieldN);
          }

        tsID++;
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
Runpctl(void *process)
{
  RunPctl pctl;
  pctl.init(process);
  pctl.run();
  pctl.close();

  return nullptr;
}
