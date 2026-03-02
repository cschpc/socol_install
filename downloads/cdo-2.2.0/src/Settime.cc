/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Settime    setdate         Set date
      Settime    settime         Set time
      Settime    setday          Set day
      Settime    setmon          Set month
      Settime    setyear         Set year
      Settime    settunits       Set time units
      Settime    settaxis        Set time axis
      Settime    setreftime      Set reference time
      Settime    setcalendar     Set calendar
      Settime    shifttime       Shift timesteps
*/

#include "cdi.h"
#include "julian_date.h"

#include "cdo_options.h"
#include "process_int.h"
#include "param_conversion.h"
#include "util_string.h"
#include "datetime.h"
#include "printinfo.h"

static void
shift_time(int calendar, int tunit, int64_t ijulinc, CdiDateTime &vDateTime)
{
  if (tunit == TUNIT_MONTH || tunit == TUNIT_YEAR)
    {
      int year, month, day;
      cdiDate_decode(vDateTime.date, &year, &month, &day);

      month += (int) ijulinc;
      adjust_month_and_year(month, year);

      vDateTime.date = cdiDate_encode(year, month, day);
    }
  else
    {
      auto julianDate = julianDate_encode(calendar, vDateTime);
      julianDate = julianDate_add_seconds(julianDate, ijulinc);
      vDateTime = julianDate_decode(calendar, julianDate);

      if (Options::cdoVerbose)
        cdo_print("julianDate, ijulinc, vdate, vtime: %g %lld %s", julianDate_to_seconds(julianDate), ijulinc,
                  datetime_to_string(vDateTime));
    }
}

static void
time_gen_bounds(int calendar, int tunit, int incrPeriod, const CdiDateTime &vDateTime, CdiDateTime *vDateTimeBounds)
{
  cdiDateTime_init(&vDateTimeBounds[0]);
  cdiDateTime_init(&vDateTimeBounds[1]);
  vDateTimeBounds[0].date = vDateTime.date;
  vDateTimeBounds[1].date = vDateTime.date;

  int year, month, day;
  cdiDate_decode(vDateTime.date, &year, &month, &day);

  if (tunit == TUNIT_MONTH)
    {
      vDateTimeBounds[0].date = cdiDate_encode(year, month, 1);
      month++;
      if (month > 12)
        {
          month = 1;
          year++;
        }
      vDateTimeBounds[1].date = cdiDate_encode(year, month, 1);
    }
  else if (tunit == TUNIT_YEAR)
    {
      vDateTimeBounds[0].date = cdiDate_encode(year, 1, 1);
      vDateTimeBounds[1].date = cdiDate_encode(year + 1, 1, 1);
    }
  else if (tunit == TUNIT_DAY)
    {
      auto julianDate = julianDate_encode(calendar, vDateTimeBounds[0]);
      julianDate = julianDate_add_seconds(julianDate, 86400);
      vDateTimeBounds[1] = julianDate_decode(calendar, julianDate);
    }
  else if (tunit == TUNIT_HOUR || tunit == TUNIT_3HOURS || tunit == TUNIT_6HOURS || tunit == TUNIT_12HOURS)
    {
      if (incrPeriod == 0) incrPeriod = 1;
      if (incrPeriod > 24) cdo_abort("Time period must be less equal 24!");

      // clang-format off
      if      (tunit == TUNIT_3HOURS)  incrPeriod = 3;
      else if (tunit == TUNIT_6HOURS)  incrPeriod = 6;
      else if (tunit == TUNIT_12HOURS) incrPeriod = 12;
      // clang-format on

      int hour, minute, second, ms;
      cdiTime_decode(vDateTime.time, &hour, &minute, &second, &ms);
      int h0 = (hour / incrPeriod) * incrPeriod;
      vDateTimeBounds[0].time = cdiTime_encode(h0, 0, 0, 0);
      int h1 = h0 + incrPeriod;
      if (h1 >= 24)
        {
          auto julianDate = julianDate_encode(calendar, vDateTimeBounds[0]);
          julianDate = julianDate_add_seconds(julianDate, incrPeriod * 3600);
          vDateTimeBounds[1] = julianDate_decode(calendar, julianDate);
        }
      else
        vDateTimeBounds[1].time = cdiTime_encode(h1, 0, 0, 0);
    }
}

int
evaluate_calendar_string(int operatorID, const std::string &calendarName)
{
  int calendar = CALENDAR_STANDARD;
  auto calendarString = string_to_lower(calendarName);
  // clang-format off
  if      (calendarString == "standard")  calendar = CALENDAR_STANDARD;
  else if (calendarString == "gregorian") calendar = CALENDAR_GREGORIAN;
  else if (calendarString == "proleptic") calendar = CALENDAR_PROLEPTIC;
  else if (calendarString == "proleptic_gregorian") calendar = CALENDAR_PROLEPTIC;
  else if (calendarString == "360days") calendar = CALENDAR_360DAYS;
  else if (calendarString == "360_day") calendar = CALENDAR_360DAYS;
  else if (calendarString == "365days") calendar = CALENDAR_365DAYS;
  else if (calendarString == "365_day") calendar = CALENDAR_365DAYS;
  else if (calendarString == "366days") calendar = CALENDAR_366DAYS;
  else if (calendarString == "366_day") calendar = CALENDAR_366DAYS;
  else cdo_abort("Calendar >%s< unsupported! Available %s", calendarName, cdo_operator_enter(operatorID));
  // clang-format on

  return calendar;
}

static CdiDateTime
argument2datetimeinc(int &incrPeriod, int &incrUnits, int &timeUnits)
{
  if (cdo_operator_argc() < 1) cdo_abort("Too few arguments!");
  if (cdo_operator_argc() > 3) cdo_abort("Too many arguments!");

  CdiDateTime sDateTime{};
  sDateTime.date = decode_datestring(cdo_operator_argv(0));

  if (cdo_operator_argc() > 1)
    {
      sDateTime.time = decode_timestring(cdo_operator_argv(1));
      if (cdo_operator_argc() == 3) decode_timeunits(cdo_operator_argv(2), incrPeriod, incrUnits, timeUnits);
    }

  return sDateTime;
}

static bool
timeunits_is_valid(int timeUnits)
{
  return (timeUnits == TUNIT_HOUR || timeUnits == TUNIT_3HOURS || timeUnits == TUNIT_6HOURS || timeUnits == TUNIT_12HOURS
          || timeUnits == TUNIT_DAY || timeUnits == TUNIT_MONTH || timeUnits == TUNIT_YEAR);
}

void *
Settime(void *process)
{
  int64_t newval = 0;
  int taxisID2 = CDI_UNDEFID;
  int timeUnits = TUNIT_DAY;
  int64_t ijulinc = 0;
  int incrPeriod = 1, incrUnits = 86400;
  int year = 1, month = 1, day = 1;
  int day0 = 0;
  bool copyTimestep = false;
  int newcalendar = CALENDAR_STANDARD;
  // int nargs;
  CdiDateTime sDateTime{};
  CdiDateTime vDateTimeBounds[2]{};
  JulianDate julianDate;

  cdo_initialize(process);

  // clang-format off
  auto SETYEAR     = cdo_operator_add("setyear",      0,  1, "year");
  auto SETMON      = cdo_operator_add("setmon",       0,  1, "month");
  auto SETDAY      = cdo_operator_add("setday",       0,  1, "day");
  auto SETDATE     = cdo_operator_add("setdate",      0,  1, "date (format: YYYY-MM-DD)");
  auto SETTIME     = cdo_operator_add("settime",      0,  1, "time (format: hh:mm:ss)");
  auto SETTUNITS   = cdo_operator_add("settunits",    0,  1, "time units (seconds, minutes, hours, days, months, years)");
  auto SETTAXIS    = cdo_operator_add("settaxis",     0, -2, "date<,time<,increment>> (format YYYY-MM-DD,hh:mm:ss)");
  auto SETTBOUNDS  = cdo_operator_add("settbounds",   0,  1, "frequency (hour, day, month, year)");
  auto SETREFTIME  = cdo_operator_add("setreftime",   0, -2, "date<,time<,units>> (format YYYY-MM-DD,hh:mm:ss)");
  auto SETCALENDAR = cdo_operator_add("setcalendar",  0,  1, "calendar (standard, proleptic_gregorian, 360_day, 365_day, 366_day)");
  auto SHIFTTIME   = cdo_operator_add("shifttime",    0,  1, "shift value");
  // clang-format on

  auto operatorID = cdo_operator_id();
  // nargs = cdo_operator_f2(operatorID);

  operator_input_arg(cdo_operator_enter(operatorID));

  if (operatorID == SETTAXIS || operatorID == SETREFTIME)
    {
      sDateTime = argument2datetimeinc(incrPeriod, incrUnits, timeUnits);
      // increment in seconds
      ijulinc = (int64_t) incrPeriod * incrUnits;
    }
  else if (operatorID == SETDATE)
    {
      operator_check_argc(1);
      sDateTime.date = decode_datestring(cdo_operator_argv(0));
    }
  else if (operatorID == SETTIME)
    {
      operator_check_argc(1);
      sDateTime.time = decode_timestring(cdo_operator_argv(0));
    }
  else if (operatorID == SHIFTTIME)
    {
      operator_check_argc(1);
      decode_timeunits(cdo_operator_argv(0), incrPeriod, incrUnits, timeUnits);
      // increment in seconds
      ijulinc = (int64_t) incrPeriod * incrUnits;
    }
  else if (operatorID == SETTUNITS || operatorID == SETTBOUNDS)
    {
      operator_check_argc(1);
      decode_timeunits(cdo_operator_argv(0), incrPeriod, incrUnits, timeUnits);

      if (operatorID == SETTBOUNDS && !timeunits_is_valid(timeUnits))
        cdo_abort("Unsupported frequency %s! Use hour, 3hours, 6hours, day, month or year.", cdo_operator_argv(0));
    }
  else if (operatorID == SETCALENDAR)
    {
      operator_check_argc(1);
      auto cname = cdo_operator_argv(0);
      newcalendar = evaluate_calendar_string(operatorID, cname);
    }
  else
    {
      operator_check_argc(1);
      newval = parameter_to_int(cdo_operator_argv(0));
    }

  auto streamID1 = cdo_open_read(0);

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto vlistID2 = vlistDuplicate(vlistID1);

  auto taxisID1 = vlistInqTaxis(vlistID1);
  auto taxis_has_bounds = (taxisHasBounds(taxisID1) > 0);
  auto ntsteps = vlistNtsteps(vlistID1);
  auto nvars = vlistNvars(vlistID1);

  VarList varList1;
  varListInit(varList1, vlistID1);

  if (ntsteps == 1 && varList_numVaryingVars(varList1) == 0) ntsteps = 0;

  if (ntsteps == 0)
    {
      for (int varID = 0; varID < nvars; ++varID) vlistDefVarTimetype(vlistID2, varID, TIME_VARYING);
    }

  auto calendar = taxisInqCalendar(taxisID1);

  if (Options::cdoVerbose) cdo_print("calendar = %d", calendar);

  if (operatorID == SETREFTIME)
    {
      copyTimestep = true;

      if (taxisInqType(taxisID1) == TAXIS_ABSOLUTE)
        {
          cdo_print("Changing absolute to relative time axis!");
          taxisID2 = cdo_taxis_create(TAXIS_RELATIVE);
        }
      else
        taxisID2 = taxisDuplicate(taxisID1);

      if (cdo_operator_argc() != 3) timeUnits = taxisInqTunit(taxisID1);
      taxisDefTunit(taxisID2, timeUnits);
    }
  else if (operatorID == SETTUNITS)
    {
      copyTimestep = true;

      if (taxisInqType(taxisID1) == TAXIS_ABSOLUTE)
        {
          cdo_print("Changing absolute to relative time axis!");
          taxisID2 = cdo_taxis_create(TAXIS_RELATIVE);
          taxisDefTunit(taxisID2, timeUnits);
        }
      else
        taxisID2 = taxisDuplicate(taxisID1);
    }
  else if (operatorID == SETCALENDAR)
    {
      copyTimestep = true;
      // if ( ((char *)argument)[0] == '-' ) cdo_abort("This operator does not work with pipes!");
      if (taxisInqType(taxisID1) == TAXIS_ABSOLUTE)
        {
          // if ( CdoDefault::FileType != CDI_FILETYPE_NC ) cdo_abort("This operator does not work on an absolute time axis!");
          cdo_print("Changing absolute to relative time axis!");
          taxisID2 = cdo_taxis_create(TAXIS_RELATIVE);
        }
      else
        taxisID2 = taxisDuplicate(taxisID1);
    }
  else
    taxisID2 = taxisDuplicate(taxisID1);

  if (operatorID == SETTAXIS)
    {
      taxisDefTunit(taxisID2, timeUnits);
      taxisDefRdatetime(taxisID2, sDateTime);
      julianDate = julianDate_encode(calendar, sDateTime);
    }
  else if (operatorID == SETTUNITS) { taxisDefTunit(taxisID2, timeUnits); }
  else if (operatorID == SETCALENDAR) { taxisDefCalendar(taxisID2, newcalendar); }
  else if (operatorID == SETTBOUNDS) { taxisWithBounds(taxisID2); }

  if (operatorID != SHIFTTIME)
    if (taxis_has_bounds && !copyTimestep)
      {
        cdo_warning("Time bounds unsupported by this operator, removed!");
        taxisDeleteBounds(taxisID2);
        taxis_has_bounds = false;
      }

  vlistDefTaxis(vlistID2, taxisID2);

  CdoStreamID streamID2 = CDO_STREAM_UNDEF;

  Field field;

  int tsID = 0;
  while (true)
    {
      auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
      if (nrecs == 0) break;

      auto vDateTime = taxisInqVdatetime(taxisID1);

      if (operatorID == SETTAXIS)
        {
          if (timeUnits == TUNIT_MONTH || timeUnits == TUNIT_YEAR)
            {
              vDateTime.time = sDateTime.time;
              if (tsID == 0)
                {
                  vDateTime.date = sDateTime.date;
                  cdiDate_decode(vDateTime.date, &year, &month, &day0);
                }
              else
                {
                  month += (int) ijulinc;
                  adjust_month_and_year(month, year);
                  day = (day0 == 31) ? days_per_month(calendar, year, month) : day0;
                  vDateTime.date = cdiDate_encode(year, month, day);
                }
            }
          else
            {
              vDateTime = julianDate_decode(calendar, julianDate);
              julianDate = julianDate_add_seconds(julianDate, ijulinc);
            }
        }
      else if (operatorID == SETTBOUNDS)
        {
          time_gen_bounds(calendar, timeUnits, incrPeriod, vDateTime, vDateTimeBounds);

          if (Options::CMOR_Mode)
            {
              auto julianDate1 = julianDate_encode(calendar, vDateTimeBounds[0]);
              auto julianDate2 = julianDate_encode(calendar, vDateTimeBounds[1]);
              auto seconds = julianDate_to_seconds(julianDate_sub(julianDate2, julianDate1)) / 2;
              auto julianDatem = julianDate_add_seconds(julianDate1, lround(seconds));
              vDateTime = julianDate_decode(calendar, julianDatem);
            }
        }
      else if (operatorID == SHIFTTIME)
        {
          shift_time(calendar, timeUnits, ijulinc, vDateTime);
          if (taxis_has_bounds)
            {
              taxisInqVdatetimeBounds(taxisID1, &vDateTimeBounds[0], &vDateTimeBounds[1]);
              shift_time(calendar, timeUnits, ijulinc, vDateTimeBounds[0]);
              shift_time(calendar, timeUnits, ijulinc, vDateTimeBounds[1]);
            }
        }
      else if (operatorID == SETREFTIME)
        {
          if (ntsteps == 0) vDateTime = sDateTime;
        }
      else if (operatorID == SETCALENDAR || operatorID == SETTUNITS) {}
      else
        {
          cdiDate_decode(vDateTime.date, &year, &month, &day);

          if (operatorID == SETYEAR) year = newval;
          if (operatorID == SETMON) month = newval;
          if (operatorID == SETMON && (month < 0 || month > 16)) cdo_abort("parameter month=%d out of range!", month);
          if (operatorID == SETDAY) day = newval;
          if (operatorID == SETDAY && (day < 0 || day > 31)) cdo_abort("parameter day=%d out of range!", day);

          vDateTime.date = cdiDate_encode(year, month, day);

          if (operatorID == SETDATE) vDateTime.date = sDateTime.date;
          if (operatorID == SETTIME) vDateTime.time = sDateTime.time;
        }

      if (copyTimestep)
        {
          cdo_taxis_copy_timestep(taxisID2, taxisID1);
          if (operatorID == SETREFTIME) taxisDefRdatetime(taxisID2, sDateTime);
          if (operatorID == SETREFTIME && ntsteps == 0) taxisDefVdatetime(taxisID2, vDateTime);
        }
      else
        {
          auto numavg = taxisInqNumavg(taxisID1);
          taxisDefNumavg(taxisID2, numavg);

          taxisDefVdatetime(taxisID2, vDateTime);

          if (taxis_has_bounds || operatorID == SETTBOUNDS)
            taxisDefVdatetimeBounds(taxisID2, vDateTimeBounds[0], vDateTimeBounds[1]);
        }

      if (streamID2 == CDO_STREAM_UNDEF)
        {
          streamID2 = cdo_open_write(1);
          cdo_def_vlist(streamID2, vlistID2);
        }

      cdo_def_timestep(streamID2, tsID);

      for (int recID = 0; recID < nrecs; ++recID)
        {
          int varID, levelID;
          cdo_inq_record(streamID1, &varID, &levelID);
          cdo_def_record(streamID2, varID, levelID);

          field.init(varList1[varID]);
          cdo_read_record(streamID1, field);
          cdo_write_record(streamID2, field);
        }

      tsID++;
    }

  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  cdo_finish();

  return nullptr;
}
