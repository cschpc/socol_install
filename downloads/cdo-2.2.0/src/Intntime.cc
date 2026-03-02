/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Intntime   intntime        Time interpolation
*/

#include "cdi.h"
#include "julian_date.h"

#include "cdo_options.h"
#include "process_int.h"
#include "cdo_vlist.h"
#include "param_conversion.h"
#include "datetime.h"
#include "printinfo.h"

void interp_time(double fac1, double fac2, const double *array1, const double *array2, Field &field3, bool withMissval);

void *
Intntime(void *process)
{
  cdo_initialize(process);

  operator_input_arg("number of timesteps between 2 timesteps");
  if (cdo_operator_argc() < 1) cdo_abort("Too few arguments!");

  auto numts = parameter_to_int(cdo_operator_argv(0));
  if (numts < 2) cdo_abort("parameter must be greater than 1!");

  auto streamID1 = cdo_open_read(0);

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto vlistID2 = vlistDuplicate(vlistID1);

  VarList varList;
  varListInit(varList, vlistID1);

  auto nvars = vlistNvars(vlistID1);

  auto maxrecs = vlistNrecs(vlistID1);
  std::vector<RecordInfo> recList(maxrecs);

  Field field3;

  Varray3D<size_t> nmiss(2);
  nmiss[0].resize(nvars);
  nmiss[1].resize(nvars);
  Varray3D<double> vardata(2);
  vardata[0].resize(nvars);
  vardata[1].resize(nvars);

  for (int varID = 0; varID < nvars; ++varID)
    {
      auto gridsize = varList[varID].gridsize;
      auto nlevel = varList[varID].nlevels;
      nmiss[0][varID].resize(nlevel);
      nmiss[1][varID].resize(nlevel);
      vardata[0][varID].resize(gridsize * nlevel);
      vardata[1][varID].resize(gridsize * nlevel);
    }

  auto taxisID1 = vlistInqTaxis(vlistID1);
  auto taxisID2 = taxisDuplicate(taxisID1);
  if (taxisHasBounds(taxisID2)) taxisDeleteBounds(taxisID2);

  vlistDefNtsteps(vlistID2, -1);
  vlistDefTaxis(vlistID2, taxisID2);

  auto streamID2 = cdo_open_write(1);

  cdo_def_vlist(streamID2, vlistID2);

  auto calendar = taxisInqCalendar(taxisID1);

  int curFirst = 0, curSecond = 1;

  int tsID = 0;
  int tsIDo = 0;
  auto nrecs = cdo_stream_inq_timestep(streamID1, tsID++);
  auto julianDate1 = julianDate_encode(calendar, taxisInqVdatetime(taxisID1));

  cdo_taxis_copy_timestep(taxisID2, taxisID1);
  cdo_def_timestep(streamID2, tsIDo++);
  for (int recID = 0; recID < nrecs; ++recID)
    {
      int varID, levelID;
      cdo_inq_record(streamID1, &varID, &levelID);
      auto offset = varList[varID].gridsize * levelID;
      auto single1 = &vardata[curFirst][varID][offset];
      cdo_read_record(streamID1, single1, &nmiss[curFirst][varID][levelID]);

      cdo_def_record(streamID2, varID, levelID);
      cdo_write_record(streamID2, single1, nmiss[curFirst][varID][levelID]);
    }

  while (true)
    {
      nrecs = cdo_stream_inq_timestep(streamID1, tsID++);
      if (nrecs == 0) break;

      auto vDateTime2 = taxisInqVdatetime(taxisID1);
      auto julianDate2 = julianDate_encode(calendar, vDateTime2);

      for (int recID = 0; recID < nrecs; ++recID)
        {
          int varID, levelID;
          cdo_inq_record(streamID1, &varID, &levelID);

          recList[recID].set(varID, levelID);

          auto offset = varList[varID].gridsize * levelID;
          auto single2 = &vardata[curSecond][varID][offset];
          cdo_read_record(streamID1, single2, &nmiss[curSecond][varID][levelID]);
        }

      for (int it = 1; it < numts; it++)
        {
          auto seconds = it * julianDate_to_seconds(julianDate_sub(julianDate2, julianDate1)) / numts;
          auto julianDate = julianDate_add_seconds(julianDate1, lround(seconds));
          auto dt = julianDate_decode(calendar, julianDate);

          if (Options::cdoVerbose) cdo_print("%s %s", date_to_string(dt.date), time_to_string(dt.time));

          taxisDefVdatetime(taxisID2, dt);
          cdo_def_timestep(streamID2, tsIDo++);

          auto diff = julianDate_to_seconds(julianDate_sub(julianDate2, julianDate1));
          auto fac1 = julianDate_to_seconds(julianDate_sub(julianDate2, julianDate)) / diff;
          auto fac2 = julianDate_to_seconds(julianDate_sub(julianDate, julianDate1)) / diff;

          for (int recID = 0; recID < nrecs; ++recID)
            {
              auto [varID, levelID] = recList[recID].get();

              auto offset = varList[varID].gridsize * levelID;
              auto single1 = &vardata[curFirst][varID][offset];
              auto single2 = &vardata[curSecond][varID][offset];

              field3.init(varList[varID]);

              auto withMissval = (nmiss[curFirst][varID][levelID] || nmiss[curSecond][varID][levelID]);
              interp_time(fac1, fac2, single1, single2, field3, withMissval);

              cdo_def_record(streamID2, varID, levelID);
              cdo_write_record(streamID2, field3);
            }
        }

      taxisDefVdatetime(taxisID2, vDateTime2);
      cdo_def_timestep(streamID2, tsIDo++);
      for (int recID = 0; recID < nrecs; ++recID)
        {
          auto [varID, levelID] = recList[recID].get();

          auto offset = varList[varID].gridsize * levelID;
          auto single2 = &vardata[curSecond][varID][offset];

          cdo_def_record(streamID2, varID, levelID);
          cdo_write_record(streamID2, single2, nmiss[curSecond][varID][levelID]);
        }

      julianDate1 = julianDate2;
      std::swap(curFirst, curSecond);
    }

  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  cdo_finish();

  return nullptr;
}
