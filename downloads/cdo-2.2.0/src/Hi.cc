/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Copyright (C) 2006 Brockmann Consult

  Author: Ralf Quast

*/

/*
   This module contains the following operators:

      Hi      hi           Compute the humidity index
*/

#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "cdo_vlist.h"

static const char HI_NAME[] = "hum_index";
static const char HI_LONGNAME[]
    = "Humindex describes empirically in units of temperature how the temperature and humidity influence the wellness of a human "
      "being. HI = T + 5/9 * (A - 10) with A = e * (6.112 * 10 ** ((7.5 * T)/(237.7 + T)) * R), T  = air temperature in degree "
      "Celsius, R = relative humidity, e = vapour pressure. Humindex is only defined for temperatures of at least 26 degree "
      "Celsius and relative humidity of at least 40 percent.";
static const char HI_UNITS[] = "Celsius";

static const int FIRST_VAR = 0;

static double
humidityIndex(double t, double e, double r, double missval)
{
  static const double tmin = 26.0;
  static const double rmin = 40.0;

  if (t < tmin || r < rmin) return missval;

  return t + (5.0 / 9.0) * ((0.01 * r * e * 6.112 * std::pow(10.0, (7.5 * t) / (237.7 + t))) - 10.0);
}

static void
farexpr(Field &field1, const Field &field2, const Field &field3, double (*expression)(double, double, double, double))
{
  const auto grid1 = field1.grid;
  const auto grid2 = field2.grid;
  const auto grid3 = field3.grid;
  const auto missval1 = field1.missval;
  const auto missval2 = field2.missval;
  const auto missval3 = field3.missval;

  const auto len = gridInqSize(grid1);
  if (len != gridInqSize(grid2) || len != gridInqSize(grid3)) cdo_abort("Fields have different gridsize (%s)", __func__);

  if (field1.nmiss || field2.nmiss || field3.nmiss)
    {
      for (size_t i = 0; i < len; ++i)
        if (DBL_IS_EQUAL(field1.vec_d[i], missval1) || DBL_IS_EQUAL(field2.vec_d[i], missval2)
            || DBL_IS_EQUAL(field3.vec_d[i], missval3))
          field1.vec_d[i] = missval1;
        else
          field1.vec_d[i] = expression(field1.vec_d[i], field2.vec_d[i], field3.vec_d[i], missval1);
    }
  else
    {
      for (size_t i = 0; i < len; ++i) field1.vec_d[i] = expression(field1.vec_d[i], field2.vec_d[i], field3.vec_d[i], missval1);
    }

  field1.nmiss = field_num_miss(field1);
}

class ModuleHi
{
  CdoStreamID streamID1;
  CdoStreamID streamID2;
  CdoStreamID streamID3;
  CdoStreamID streamID4;

  int vlistID1;
  int vlistID2;
  int vlistID3;
  int vlistID4;

  int taxisID1;
  // int taxisID2
  // int taxisID3
  int taxisID4;
  int varID4;

  Field field1, field2, field3;

public:
  void
  init(void *process)
  {
    cdo_initialize(process);

    cdo_operator_add("hi", 0, 0, nullptr);

    operator_check_argc(0);

    streamID1 = cdo_open_read(0);
    streamID2 = cdo_open_read(1);
    streamID3 = cdo_open_read(2);

    vlistID1 = cdo_stream_inq_vlist(streamID1);
    vlistID2 = cdo_stream_inq_vlist(streamID2);
    vlistID3 = cdo_stream_inq_vlist(streamID3);

    taxisID1 = vlistInqTaxis(vlistID1);
    // taxisID2 = vlistInqTaxis(vlistID2);
    // taxisID3 = vlistInqTaxis(vlistID3);

    vlist_compare(vlistID1, vlistID2, CMP_DIM);
    vlist_compare(vlistID1, vlistID3, CMP_DIM);

    const auto gridsizemax = vlistGridsizeMax(vlistID1);

    field1.resize(gridsizemax);
    field2.resize(gridsizemax);
    field3.resize(gridsizemax);

    if (Options::cdoVerbose)
      cdo_print("Number of timesteps: file1 %d, file2 %d, file3 %d", vlistNtsteps(vlistID1), vlistNtsteps(vlistID2),
                vlistNtsteps(vlistID3));

    vlistID4 = vlistCreate();
    const auto gridID = vlistInqVarGrid(vlistID1, FIRST_VAR);
    const auto zaxisID = vlistInqVarZaxis(vlistID1, FIRST_VAR);
    varID4 = vlistDefVar(vlistID4, gridID, zaxisID, TIME_VARYING);

    taxisID4 = cdo_taxis_create(TAXIS_RELATIVE);
    taxisDefTunit(taxisID4, TUNIT_MINUTE);
    taxisDefCalendar(taxisID4, CALENDAR_STANDARD);
    taxisDefRdate(taxisID4, 19550101);
    taxisDefRtime(taxisID4, 0);
    vlistDefTaxis(vlistID4, taxisID4);

    cdiDefKeyString(vlistID4, varID4, CDI_KEY_NAME, HI_NAME);
    cdiDefKeyString(vlistID4, varID4, CDI_KEY_LONGNAME, HI_LONGNAME);
    cdiDefKeyString(vlistID4, varID4, CDI_KEY_UNITS, HI_UNITS);

    streamID4 = cdo_open_write(3);

    cdo_def_vlist(streamID4, vlistID4);
  }
  void
  run()
  {
    int tsID = 0;
    while (true)
      {
        const auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
        if (nrecs == 0) break;

        const auto nrecs2 = cdo_stream_inq_timestep(streamID2, tsID);
        const auto nrecs3 = cdo_stream_inq_timestep(streamID3, tsID);
        if (nrecs2 == 0 || nrecs3 == 0) cdo_abort("Input streams have different number of timesteps!");

        cdo_taxis_copy_timestep(taxisID4, taxisID1);
        cdo_def_timestep(streamID4, tsID);

        for (int recID = 0; recID < nrecs; ++recID)
          {
            int varID1, levelID1;
            cdo_inq_record(streamID1, &varID1, &levelID1);
            cdo_read_record(streamID1, field1.vec_d.data(), &field1.nmiss);

            int varID2, levelID2;
            cdo_inq_record(streamID2, &varID2, &levelID2);
            cdo_read_record(streamID2, field2.vec_d.data(), &field2.nmiss);

            int varID3, levelID3;
            cdo_inq_record(streamID3, &varID3, &levelID3);
            cdo_read_record(streamID3, field3.vec_d.data(), &field3.nmiss);

            if (varID1 != varID2 || varID1 != varID3 || levelID1 != levelID2 || levelID1 != levelID3)
              cdo_abort("Input streams have different structure!");

            if (varID1 != FIRST_VAR) continue;

            field1.grid = vlistInqVarGrid(vlistID1, varID1);
            field1.missval = vlistInqVarMissval(vlistID1, varID1);

            field2.grid = vlistInqVarGrid(vlistID2, varID2);
            field2.missval = vlistInqVarMissval(vlistID2, varID2);

            field3.grid = vlistInqVarGrid(vlistID3, varID3);
            field3.missval = vlistInqVarMissval(vlistID3, varID3);

            farexpr(field1, field2, field3, humidityIndex);

            cdo_def_record(streamID4, varID4, levelID1);
            cdo_write_record(streamID4, field1.vec_d.data(), field1.nmiss);
          }

        tsID++;
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
Hi(void *process)
{
  ModuleHi hi;
  hi.init(process);
  hi.run();
  hi.close();
  return nullptr;
}
