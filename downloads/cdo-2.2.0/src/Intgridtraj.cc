/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

*/

#include "cdi.h"
#include "julian_date.h"

#include "cdo_vlist.h"
#include "process_int.h"
#include "interpol.h"
#include "datetime.h"
#include "varray.h"
#include "printinfo.h"

static int
read_nextpos(FILE *fp, int calendar, JulianDate &julianDate, double &xpos, double &ypos)
{
  xpos = 0.0;
  ypos = 0.0;

  int year = 0, month = 0, day = 0, hour = 0, minute = 0, second = 0, ms = 0;
  int stat = std::fscanf(fp, "%d-%d-%d %d:%d:%d %lf %lf", &year, &month, &day, &hour, &minute, &second, &xpos, &ypos);

  CdiDateTime dateTime{};
  dateTime.date = cdiDate_set(10101);
  if (stat != EOF)
    {
      dateTime.date = cdiDate_encode(year, month, day);
      dateTime.time = cdiTime_encode(hour, minute, second, ms);
    }

  julianDate = julianDate_encode(calendar, dateTime);

  return stat;
}

void *
Intgridtraj(void *process)
{
  size_t nmiss = 0;
  double xpos, ypos;
  int calendar = CALENDAR_STANDARD;

  cdo_initialize(process);

  operator_input_arg("filename with grid trajectories");
  operator_check_argc(1);

  auto posfile = cdo_operator_argv(0).c_str();
  auto fp = std::fopen(posfile, "r");
  if (fp == nullptr) cdo_abort("Open failed on %s!", posfile);

  JulianDate julianDate;
  read_nextpos(fp, calendar, julianDate, xpos, ypos);

  const auto streamID1 = cdo_open_read(0);
  const auto vlistID1 = cdo_stream_inq_vlist(streamID1);

  VarList varList1;
  varListInit(varList1, vlistID1);

  const auto nvars = vlistNvars(vlistID1);

  const auto maxrecs = vlistNrecs(vlistID1);
  std::vector<RecordInfo> recList(maxrecs);

  const auto gridsizemax = vlistGridsizeMax(vlistID1);
  Field field1, field2;
  field1.resize(gridsizemax);
  field2.resize(1);

  Varray2D<double> vardata1(nvars), vardata2(nvars);

  for (int varID = 0; varID < nvars; ++varID)
    {
      const auto gridsize = varList1[varID].gridsize;
      const auto nlevels = varList1[varID].nlevels;
      vardata1[varID].resize(gridsize * nlevels);
      vardata2[varID].resize(gridsize * nlevels);
    }

  const auto gridID2 = gridCreate(GRID_TRAJECTORY, 1);
  gridDefXsize(gridID2, 1);
  gridDefYsize(gridID2, 1);
  gridDefXvals(gridID2, &xpos);
  gridDefYvals(gridID2, &ypos);

  const auto vlistID2 = vlistDuplicate(vlistID1);

  const auto ngrids = vlistNgrids(vlistID1);
  for (int index = 0; index < ngrids; ++index)
    {
      const auto gridID1 = vlistGrid(vlistID1, index);

      if (gridInqType(gridID1) != GRID_LONLAT && gridInqType(gridID1) != GRID_GAUSSIAN)
        cdo_abort("Unsupported grid type: %s", gridNamePtr(gridInqType(gridID1)));

      vlistChangeGridIndex(vlistID2, index, gridID2);
    }

  const auto taxisID1 = vlistInqTaxis(vlistID1);
  const auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  CdoStreamID streamID2 = CDO_STREAM_UNDEF;

  int tsID = 0;
  auto nrecs = cdo_stream_inq_timestep(streamID1, tsID++);
  auto julianDate1 = julianDate_encode(calendar, taxisInqVdatetime(taxisID1));
  for (int recID = 0; recID < nrecs; ++recID)
    {
      int varID, levelID;
      cdo_inq_record(streamID1, &varID, &levelID);
      const auto gridsize = varList1[varID].gridsize;
      const auto offset = gridsize * levelID;
      auto single1 = &vardata1[varID][offset];
      cdo_read_record(streamID1, single1, &nmiss);
      if (nmiss) cdo_abort("Missing values unsupported for this operator!");
    }

  int tsIDo = 0;
  while (julianDate_to_seconds(julianDate1) <= julianDate_to_seconds(julianDate))
    {
      nrecs = cdo_stream_inq_timestep(streamID1, tsID++);
      if (nrecs == 0) break;
      auto julianDate2 = julianDate_encode(calendar, taxisInqVdatetime(taxisID1));

      for (int recID = 0; recID < nrecs; ++recID)
        {
          int varID, levelID;
          cdo_inq_record(streamID1, &varID, &levelID);

          recList[recID].set(varID, levelID);

          const auto gridsize = varList1[varID].gridsize;
          const auto offset = gridsize * levelID;
          auto single2 = &vardata2[varID][offset];
          cdo_read_record(streamID1, single2, &nmiss);
          if (nmiss) cdo_abort("Missing values unsupported for this operator!");
        }

      while (julianDate_to_seconds(julianDate) < julianDate_to_seconds(julianDate2))
        {
          if (julianDate_to_seconds(julianDate) >= julianDate_to_seconds(julianDate1)
              && julianDate_to_seconds(julianDate) < julianDate_to_seconds(julianDate2))
            {
              if (streamID2 == CDO_STREAM_UNDEF)
                {
                  streamID2 = cdo_open_write(1);
                  cdo_def_vlist(streamID2, vlistID2);
                }

              taxisDefVdatetime(taxisID2, julianDate_decode(calendar, julianDate));
              cdo_def_timestep(streamID2, tsIDo++);

              const auto deltat = julianDate_to_seconds(julianDate_sub(julianDate2, julianDate1));
              const auto fac1 = julianDate_to_seconds(julianDate_sub(julianDate2, julianDate)) / deltat;
              const auto fac2 = julianDate_to_seconds(julianDate_sub(julianDate, julianDate1)) / deltat;
              // printf("      %f %f %f %f %f\n", julianDate_to_seconds(julianDate), julianDate_to_seconds(julianDate1),
              // julianDate_to_seconds(julianDate2), fac1, fac2);
              for (int recID = 0; recID < nrecs; ++recID)
                {
                  auto [varID, levelID] = recList[recID].get();

                  const auto gridsize = varList1[varID].gridsize;
                  const auto offset = gridsize * levelID;
                  auto single1 = &vardata1[varID][offset];
                  auto single2 = &vardata2[varID][offset];

                  for (size_t i = 0; i < gridsize; ++i) field1.vec_d[i] = single1[i] * fac1 + single2[i] * fac2;

                  field1.grid = varList1[varID].gridID;
                  field1.missval = varList1[varID].missval;
                  field1.nmiss = nmiss;
                  field2.grid = gridID2;
                  field2.nmiss = 0;

                  intgridbil(field1, field2);

                  cdo_def_record(streamID2, varID, levelID);
                  cdo_write_record(streamID2, field2.vec_d.data(), field2.nmiss);
                }
            }

          if (read_nextpos(fp, calendar, julianDate, xpos, ypos) == EOF) break;
          gridDefXvals(gridID2, &xpos);
          gridDefYvals(gridID2, &ypos);
        }

      julianDate1 = julianDate2;
      for (int varID = 0; varID < nvars; ++varID)
        {
          auto vardatap = vardata1[varID];
          vardata1[varID] = vardata2[varID];
          vardata2[varID] = vardatap;
        }
    }

  if (tsIDo == 0)
    {
      const auto dt = julianDate_decode(calendar, julianDate);
      cdo_warning("Date/time %s %s not found!", date_to_string(dt.date), time_to_string(dt.time));
    }

  std::fclose(fp);
  if (streamID2 != CDO_STREAM_UNDEF) cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  cdo_finish();

  return nullptr;
}
