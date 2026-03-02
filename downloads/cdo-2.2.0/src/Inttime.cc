/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Inttime    inttime         Time interpolation
*/

#include "cdi.h"
#include "julian_date.h"

#include <utility>

#include "cdo_options.h"
#include "process_int.h"
#include "cdo_vlist.h"
#include "datetime.h"
#include "printinfo.h"

template <typename T>
size_t
interp_time(double fac1, double fac2, size_t gridsize, const double *v1, const double *v2, Varray<T> &v3, bool withMissval,
            double missval)
{
  size_t nmiss3 = 0;

  if (withMissval)
    {
      for (size_t i = 0; i < gridsize; ++i)
        {
          if (!dbl_is_equal(v1[i], missval) && !dbl_is_equal(v2[i], missval))
            v3[i] = v1[i] * fac1 + v2[i] * fac2;
          else if (dbl_is_equal(v1[i], missval) && !dbl_is_equal(v2[i], missval) && fac2 >= 0.5)
            v3[i] = v2[i];
          else if (dbl_is_equal(v2[i], missval) && !dbl_is_equal(v1[i], missval) && fac1 >= 0.5)
            v3[i] = v1[i];
          else
            {
              v3[i] = missval;
              nmiss3++;
            }
        }
    }
  else
    {
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
      for (size_t i = 0; i < gridsize; ++i) v3[i] = v1[i] * fac1 + v2[i] * fac2;
    }

  return nmiss3;
}

void
interp_time(double fac1, double fac2, const double *array1, const double *array2, Field &field3, bool withMissval)
{
  if (field3.memType == MemType::Float)
    field3.nmiss = interp_time(fac1, fac2, field3.gridsize, array1, array2, field3.vec_f, withMissval, field3.missval);
  else
    field3.nmiss = interp_time(fac1, fac2, field3.gridsize, array1, array2, field3.vec_d, withMissval, field3.missval);
}

static void
julianDate_add_increment(JulianDate &julianDate, int64_t ijulinc, int calendar, int timeUnits)
{
  if (timeUnits == TUNIT_MONTH || timeUnits == TUNIT_YEAR)
    {
      auto vDateTime = julianDate_decode(calendar, julianDate);

      int year, month, day;
      cdiDate_decode(vDateTime.date, &year, &month, &day);

      month += (int) ijulinc;
      adjust_month_and_year(month, year);

      vDateTime.date = cdiDate_encode(year, month, day);
      julianDate = julianDate_encode(calendar, vDateTime);
    }
  else
    {
      julianDate = julianDate_add_seconds(julianDate, ijulinc);
    }
}

static void
adjust_time_units(int taxisID, int timeUnitsOut)
{
  auto timeUnitsIn = taxisInqTunit(taxisID);
  if ((timeUnitsIn == TUNIT_MONTH || timeUnitsIn == TUNIT_YEAR) && timeUnitsOut != TUNIT_MONTH && timeUnitsIn != TUNIT_YEAR)
    {
      taxisDefTunit(taxisID, timeUnitsOut);
    }
}

void *
Inttime(void *process)
{
  CdoStreamID streamID2 = CDO_STREAM_UNDEF;
  CdiDateTime sDateTime{};
  int incrPeriod = 0, incrUnits = 3600, timeUnits = TUNIT_HOUR;

  cdo_initialize(process);

  operator_input_arg("date,time<,increment> (format YYYY-MM-DD,hh:mm:ss)");
  if (cdo_operator_argc() < 2) cdo_abort("Too few arguments!");

  sDateTime.date = decode_datestring(cdo_operator_argv(0));
  sDateTime.time = decode_timestring(cdo_operator_argv(1));
  if (cdo_operator_argc() == 3) decode_timeunits(cdo_operator_argv(2), incrPeriod, incrUnits, timeUnits);

  // increment in seconds
  auto ijulinc = (int64_t) incrPeriod * incrUnits;

  auto streamID1 = cdo_open_read(0);

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto vlistID2 = vlistDuplicate(vlistID1);

  VarList varList;
  varListInit(varList, vlistID1);

  if (ijulinc == 0) vlistDefNtsteps(vlistID2, 1);

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
  adjust_time_units(taxisID2, timeUnits);
  if (taxisHasBounds(taxisID2)) taxisDeleteBounds(taxisID2);
  vlistDefTaxis(vlistID2, taxisID2);

  auto calendar = taxisInqCalendar(taxisID1);

  auto julianDate = julianDate_encode(calendar, sDateTime);

  if (Options::cdoVerbose)
    {
      cdo_print("Start date/time %s", datetime_to_string(sDateTime));
      cdo_print("julianDate = %f", julianDate_to_seconds(julianDate));
      cdo_print("ijulinc = %lld", ijulinc);
    }

  int curFirst = 0, curSecond = 1;

  int tsID = 0;
  auto nrecs = cdo_stream_inq_timestep(streamID1, tsID++);
  auto vDateTime1 = taxisInqVdatetime(taxisID1);
  auto julianDate1 = julianDate_encode(calendar, vDateTime1);
  for (int recID = 0; recID < nrecs; ++recID)
    {
      int varID, levelID;
      cdo_inq_record(streamID1, &varID, &levelID);
      auto offset = varList[varID].gridsize * levelID;
      auto single1 = &vardata[curFirst][varID][offset];
      cdo_read_record(streamID1, single1, &nmiss[curFirst][varID][levelID]);
    }

  if (Options::cdoVerbose)
    {
      cdo_print("Dataset begins on %s", datetime_to_string(vDateTime1));
      cdo_print("julianDate1 = %f", julianDate_to_seconds(julianDate1));
    }

  if (julianDate_to_seconds(julianDate1) > julianDate_to_seconds(julianDate))
    {
      cdo_print("Dataset begins on %s", datetime_to_string(vDateTime1));
      cdo_warning("The start time %s is before the beginning of the dataset!", datetime_to_string(sDateTime));
    }

  int tsIDo = 0;
  while (julianDate_to_seconds(julianDate1) <= julianDate_to_seconds(julianDate))
    {
      nrecs = cdo_stream_inq_timestep(streamID1, tsID++);
      if (nrecs == 0) break;

      auto vDateTime = taxisInqVdatetime(taxisID1);
      auto julianDate2 = julianDate_encode(calendar, vDateTime);
      if (Options::cdoVerbose)
        {
          cdo_print("date/time: %s", datetime_to_string(vDateTime));
          cdo_print("julianDate2 = %f", julianDate_to_seconds(julianDate2));
        }

      for (int recID = 0; recID < nrecs; ++recID)
        {
          int varID, levelID;
          cdo_inq_record(streamID1, &varID, &levelID);

          recList[recID].set(varID, levelID);

          auto offset = varList[varID].gridsize * levelID;
          auto single2 = &vardata[curSecond][varID][offset];
          cdo_read_record(streamID1, single2, &nmiss[curSecond][varID][levelID]);
        }

      while (julianDate_to_seconds(julianDate) <= julianDate_to_seconds(julianDate2))
        {
          if (julianDate_to_seconds(julianDate) >= julianDate_to_seconds(julianDate1)
              && julianDate_to_seconds(julianDate) <= julianDate_to_seconds(julianDate2))
            {
              auto dt = julianDate_decode(calendar, julianDate);

              if (Options::cdoVerbose)
                cdo_print("%s %s  %f  %d", date_to_string(dt.date), time_to_string(dt.time), julianDate_to_seconds(julianDate),
                          calendar);

              if (streamID2 == CDO_STREAM_UNDEF)
                {
                  streamID2 = cdo_open_write(1);
                  cdo_def_vlist(streamID2, vlistID2);
                }

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

          if (ijulinc == 0) break;

          julianDate_add_increment(julianDate, ijulinc, calendar, timeUnits);
        }

      julianDate1 = julianDate2;
      std::swap(curFirst, curSecond);
    }

  if (streamID2 != CDO_STREAM_UNDEF) cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  if (tsIDo == 0) cdo_warning("Start date/time %s out of range, no time steps interpolated!", datetime_to_string(sDateTime));

  cdo_finish();

  return nullptr;
}
