/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Copyright (C) 2006 Brockmann Consult

  Author: Ralf Quast

*/

/*
   This module contains the following operators:

      Wct     wct          Compute the windchill temperature (degree C)
*/

#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "cdo_vlist.h"

static const char WCT_NAME[] = "wind_chill_temperature";
static const char WCT_LONGNAME[] = "Windchill temperature describes the fact that low temperatures are felt "
                                   "to be even lower in case of wind. It is based on the rate of heat loss "
                                   "from exposed skin caused by wind and cold. It is calculated according "
                                   "to the empirical formula: 33 + (T - 33) * (0.478 + 0.237 * ( "
                                   "SQRT(ff*3.6) - 0.0124 * ff * 3.6)) with T  = air temperature in "
                                   "degree Celsius, ff = 10 m wind speed in m/s. Windchill temperature is "
                                   "only defined for temperatures at or below 33 degree Celsius and wind "
                                   "speeds above 1.39 m/s. It is mainly used for freezing temperatures.";
static const char WCT_UNITS[] = "Celsius";

static const int FIRST_VAR = 0;

static double
windchillTemperature(double t, double ff, double missval)
{
  constexpr double tmax = 33.0;
  constexpr double vmin = 1.39; /* minimum wind speed (m/s) */

  return ff < vmin || t > tmax ? missval : tmax + (t - tmax) * (0.478 + 0.237 * (std::sqrt(ff * 3.6) - 0.0124 * ff * 3.6));
}

static void
farexpr(Field &field1, Field &field2, double (*expression)(double, double, double))
{
  const auto missval1 = field1.missval;
  const auto missval2 = field2.missval;

  const auto len = field1.size;
  if (len != field2.size) cdo_abort("Fields have different size (%s)", __func__);

  if (field1.nmiss || field2.nmiss)
    {
      for (size_t i = 0; i < len; ++i)
        if (DBL_IS_EQUAL(field1.vec_d[i], missval1) || DBL_IS_EQUAL(field2.vec_d[i], missval2))
          field1.vec_d[i] = missval1;
        else
          field1.vec_d[i] = expression(field1.vec_d[i], field2.vec_d[i], missval1);
    }
  else
    {
      for (size_t i = 0; i < len; ++i) field1.vec_d[i] = expression(field1.vec_d[i], field2.vec_d[i], missval1);
    }

  field1.nmiss = field_num_miss(field1);
}

class ModuleWct
{
  CdoStreamID streamID1;
  CdoStreamID streamID2;
  CdoStreamID streamID3;

  int taxisID1;
  int taxisID3;

  int vlistID1;
  int vlistID2;

  int varID3;

  Field field1, field2;

public:
  void
  init(void *process)
  {
    cdo_initialize(process);
    cdo_operator_add("wct", 0, 0, nullptr);

    operator_check_argc(0);

    streamID1 = cdo_open_read(0);
    streamID2 = cdo_open_read(1);

    vlistID1 = cdo_stream_inq_vlist(streamID1);
    vlistID2 = cdo_stream_inq_vlist(streamID2);

    taxisID1 = vlistInqTaxis(vlistID1);

    vlist_compare(vlistID1, vlistID2, CMP_DIM);

    const auto gridsizemax = vlistGridsizeMax(vlistID1);

    field1.resize(gridsizemax);
    field2.resize(gridsizemax);

    if (Options::cdoVerbose) cdo_print("Number of timesteps: file1 %d, file2 %d", vlistNtsteps(vlistID1), vlistNtsteps(vlistID2));

    const auto vlistID3 = vlistCreate();
    const auto gridID = vlistInqVarGrid(vlistID1, FIRST_VAR);
    const auto zaxisID = vlistInqVarZaxis(vlistID1, FIRST_VAR);
    varID3 = vlistDefVar(vlistID3, gridID, zaxisID, TIME_VARYING);

    taxisID3 = cdo_taxis_create(TAXIS_RELATIVE);
    taxisDefTunit(taxisID3, TUNIT_MINUTE);
    taxisDefCalendar(taxisID3, CALENDAR_STANDARD);
    taxisDefRdate(taxisID3, 19550101);
    taxisDefRtime(taxisID3, 0);
    vlistDefTaxis(vlistID3, taxisID3);

    cdiDefKeyString(vlistID3, varID3, CDI_KEY_NAME, WCT_NAME);
    cdiDefKeyString(vlistID3, varID3, CDI_KEY_LONGNAME, WCT_LONGNAME);
    cdiDefKeyString(vlistID3, varID3, CDI_KEY_UNITS, WCT_UNITS);

    streamID3 = cdo_open_write(2);
    cdo_def_vlist(streamID3, vlistID3);
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
        if (nrecs2 == 0) cdo_abort("Input streams have different number of timesteps!");

        cdo_taxis_copy_timestep(taxisID3, taxisID1);
        cdo_def_timestep(streamID3, tsID);

        for (int recID = 0; recID < nrecs; ++recID)
          {
            int varID1, levelID1;
            cdo_inq_record(streamID1, &varID1, &levelID1);
            cdo_read_record(streamID1, field1.vec_d.data(), &field1.nmiss);

            int varID2, levelID2;
            cdo_inq_record(streamID2, &varID2, &levelID2);
            cdo_read_record(streamID2, field2.vec_d.data(), &field2.nmiss);

            if (varID1 != varID2 || levelID1 != levelID2) cdo_abort("Input streams have different structure!");

            if (varID1 != FIRST_VAR) continue;

            field1.missval = vlistInqVarMissval(vlistID1, varID1);
            field2.missval = vlistInqVarMissval(vlistID2, varID2);

            farexpr(field1, field2, windchillTemperature);

            cdo_def_record(streamID3, varID3, levelID1);
            cdo_write_record(streamID3, field1.vec_d.data(), field1.nmiss);
          }

        tsID++;
      }
  }

  void
  close()
  {
    cdo_stream_close(streamID3);
    cdo_stream_close(streamID2);
    cdo_stream_close(streamID1);

    cdo_finish();
  }
};
void *
Wct(void *process)
{
  ModuleWct wct;
  wct.init(process);
  wct.run();
  wct.close();
  return nullptr;
}
