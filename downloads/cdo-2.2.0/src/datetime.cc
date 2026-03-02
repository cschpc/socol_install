/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include <utility>

#include "cdo_options.h"
#include "process_int.h"
#include "datetime.h"
#include "param_conversion.h"
#include "util_string.h"

TimeStat CDO_Timestat_Date = TimeStat::UNDEF;
bool CDO_Ignore_Time_Bounds = false;
bool CDO_Use_Time_Bounds = false;
static bool dateTimeInit = false;

const char *
time_units_cstr(TimeUnits timeUnit)
{
  // clang-format off
  if      (timeUnit == TimeUnits::SECONDS) return "second";
  else if (timeUnit == TimeUnits::MINUTES) return "minute";
  else if (timeUnit == TimeUnits::HOURS)   return "hour";
  else if (timeUnit == TimeUnits::DAYS)    return "day";
  else if (timeUnit == TimeUnits::MONTHS)  return "month";
  else if (timeUnit == TimeUnits::YEARS)   return "year";
  // clang-format on

  return NULL;
}

void
set_timestat_date(const std::string &optarg)
{
  TimeStat timestatdate = TimeStat::UNDEF;

  // clang-format off
  if      (optarg == "first")   timestatdate = TimeStat::FIRST;
  else if (optarg == "last")    timestatdate = TimeStat::LAST;
  else if (optarg == "middle")  timestatdate = TimeStat::MEAN;
  else if (optarg == "midhigh") timestatdate = TimeStat::MIDHIGH;
  // clang-format on

  if (timestatdate == TimeStat::UNDEF) cdo_abort("option --%s: unsupported argument: %s", "timestat_date", optarg);

  CDO_Timestat_Date = timestatdate;
}

static void
get_timestat_date(TimeStat *tstat_date)
{
  char *envstr = getenv("CDO_TIMESTAT_DATE");
  if (envstr == nullptr) envstr = getenv("RUNSTAT_DATE");
  if (envstr)
    {
      TimeStat env_date = TimeStat::UNDEF;
      auto envstrl = string_to_lower(envstr);

      // clang-format off
      if      (envstrl == "first")   env_date = TimeStat::FIRST;
      else if (envstrl == "last")    env_date = TimeStat::LAST;
      else if (envstrl == "middle")  env_date = TimeStat::MEAN;
      else if (envstrl == "midhigh") env_date = TimeStat::MIDHIGH;
      // clang-format on

      if (env_date != TimeStat::UNDEF)
        {
          *tstat_date = env_date;
          if (Options::cdoVerbose) cdo_print("Set CDO_TIMESTAT_DATE to %s", envstr);
        }
    }
}

void
DateTimeList::init()
{
  if (!dateTimeInit) get_timestat_date(&CDO_Timestat_Date);
  dateTimeInit = true;
}

CdiDateTime
DateTimeList::get_vDateTime(int tsID)
{
  if (tsID < 0 || (size_t) tsID >= this->size) cdo_abort("Internal error; tsID out of bounds!");

  return this->dtInfo[tsID].c;
}

void
DateTimeList::shift()
{
  for (size_t inp = 0; inp < this->size - 1; ++inp) this->dtInfo[inp] = this->dtInfo[inp + 1];
}

void
DateTimeList::taxis_inq_timestep(int taxisID, int tsID)
{
  constexpr size_t NALLOC = 128;

  if ((size_t) tsID >= this->nalloc)
    {
      this->nalloc += NALLOC;
      this->dtInfo.resize(this->nalloc);
    }

  if ((size_t) tsID >= this->size) this->size = (size_t) tsID + 1;

  this->dtInfo[tsID].v = taxisInqVdatetime(taxisID);
  this->dtInfo[tsID].c = this->dtInfo[tsID].v;

  if (tsID == 0)
    {
      if (this->hasBounds == -1) this->hasBounds = CDO_Ignore_Time_Bounds ? 0 : taxisHasBounds(taxisID);
      if (this->calendar == -1) this->calendar = taxisInqCalendar(taxisID);
    }

  if (this->hasBounds)
    {
      taxisInqVdatetimeBounds(taxisID, &(this->dtInfo[tsID].b[0]), &(this->dtInfo[tsID].b[1]));

      auto time = cdiTime_get(this->dtInfo[tsID].b[1].time);
      if (CDO_Use_Time_Bounds && time == 0 && cdiDateTime_isEQ(this->dtInfo[tsID].v, this->dtInfo[tsID].b[1]))
        {
          auto julianDate1 = julianDate_encode(this->calendar, this->dtInfo[tsID].b[0]);
          auto julianDate2 = julianDate_encode(this->calendar, this->dtInfo[tsID].b[1]);

          if (julianDate_to_seconds(julianDate1) < julianDate_to_seconds(julianDate2))
            {
              auto julianDate = julianDate_add_seconds(julianDate2, -1);
              this->dtInfo[tsID].c = julianDate_decode(this->calendar, julianDate);
            }
        }
    }
  else
    {
      cdiDateTime_init(&this->dtInfo[tsID].b[0]);
      cdiDateTime_init(&this->dtInfo[tsID].b[1]);
    }
}

void
DateTimeList::taxis_set_next_timestep(int taxisID)
{
  int tsID = this->size;
  this->taxis_inq_timestep(taxisID, tsID);
}

void
DateTimeList::taxis_def_timestep(int taxisID, int tsID)
{
  if (tsID < 0 || (size_t) tsID >= this->size) cdo_abort("Internal error; tsID out of bounds!");

  taxisDefVdatetime(taxisID, this->dtInfo[tsID].v);
  if (this->hasBounds) taxisDefVdatetimeBounds(taxisID, this->dtInfo[tsID].b[0], this->dtInfo[tsID].b[1]);
}

void
DateTimeList::mean(int nsteps)
{
  if (nsteps % 2 == 0)
    {
#ifdef TEST_DTLIST_MEAN
      auto julianDate0 = julianDate_encode(this->calendar, this->dtInfo[0].v);

      double seconds = 0.0;
      for (int i = 1; i < nsteps; ++i)
        {
          auto julianDate = julianDate_encode(this->calendar, this->dtInfo[i].v);
          seconds += julianDate_to_seconds(julianDate_sub(julianDate, julianDate0));
        }

      auto julianDate = julianDate_add_seconds(julianDate0, lround(seconds / nsteps));
      this->timestat.v = julianDate_decode(this->calendar, julianDate);
#else
      auto julianDate1 = julianDate_encode(this->calendar, this->dtInfo[nsteps / 2 - 1].v);
      auto julianDate2 = julianDate_encode(this->calendar, this->dtInfo[nsteps / 2].v);

      auto seconds = julianDate_to_seconds(julianDate_sub(julianDate2, julianDate1)) / 2;
      auto julianDatem = julianDate_add_seconds(julianDate1, lround(seconds));
      this->timestat.v = julianDate_decode(this->calendar, julianDatem);
#endif
    }
  else
    {
      this->timestat.v = this->dtInfo[nsteps / 2].v;
    }
}

void
DateTimeList::midhigh(int nsteps)
{
  this->timestat.v = this->dtInfo[nsteps / 2].v;
}

void
DateTimeList::stat_taxis_def_timestep(int taxisID, int nsteps)
{
  if ((size_t) nsteps > this->size) cdo_abort("Internal error; unexpected nsteps=%d (limit=%ld)!", nsteps, this->size);

  if (CDO_Timestat_Date != TimeStat::UNDEF) this->stat = CDO_Timestat_Date;

  // clang-format off
  if      (this->stat == TimeStat::MEAN)    this->mean(nsteps);
  else if (this->stat == TimeStat::MIDHIGH) this->midhigh(nsteps);
  else if (this->stat == TimeStat::FIRST)   this->timestat.v = this->dtInfo[0].v;
  else if (this->stat == TimeStat::LAST)    this->timestat.v = this->dtInfo[nsteps - 1].v;
  else cdo_abort("Internal error; implementation missing for timestat=%d", (int)this->stat);
  // clang-format on

  if (this->hasBounds)
    {
      this->timestat.b[0] = this->dtInfo[0].b[0];
      this->timestat.b[1] = this->dtInfo[nsteps - 1].b[1];
    }
  else
    {
      this->timestat.b[0] = this->dtInfo[0].v;
      this->timestat.b[1] = this->dtInfo[nsteps - 1].v;
    }

  taxisDefVdatetime(taxisID, this->timestat.v);
  // if (this->hasBounds)
  {
    taxisDefVdatetimeBounds(taxisID, this->timestat.b[0], this->timestat.b[1]);
  }
}

void
DateTimeList::stat_taxis_def_timestep(int taxisID)
{
  int nsteps = this->size;
  this->stat_taxis_def_timestep(taxisID, nsteps);
}

CdiDateTime
datetime_avg(int calendar, int ndates, const std::vector<CdiDateTime> &cdiDateTimes)
{
  if (ndates % 2 == 0)
    {
      auto julianDate1 = julianDate_encode(calendar, cdiDateTimes[ndates / 2 - 1]);
      auto julianDate2 = julianDate_encode(calendar, cdiDateTimes[ndates / 2]);

      auto seconds = julianDate_to_seconds(julianDate_sub(julianDate2, julianDate1)) / 2;
      auto julianDatem = julianDate_add_seconds(julianDate1, lround(seconds));
      return julianDate_decode(calendar, julianDatem);
    }
  else
    {
      return cdiDateTimes[ndates / 2];
    }
}

void
adjust_month_and_year(int &month, int &year)
{
  // clang-format off
  while (month > 12) { month -= 12; year++; }
  while (month <  1) { month += 12; year--; }
  // clang-format on
}

double
delta_time_step_0(int tsID, int calendar, const CdiDateTime &vDateTime, JulianDate &julianDate0, double &deltat1)
{
  double zj = 0.0;

  auto julianDate = julianDate_encode(calendar, vDateTime);

  if (tsID == 0) { julianDate0 = julianDate; }
  else
    {
      auto deltat = julianDate_to_seconds(julianDate_sub(julianDate, julianDate0));
      if (tsID == 1) deltat1 = deltat;
      zj = deltat / deltat1;
    }

  return zj;
}

TimeIncrement
get_time_increment(double jdelta, CdiDate vDate0, CdiDate vDate1)
{
  int64_t lperiod = (jdelta < 0) ? (int64_t) (jdelta - 0.5) : (int64_t) (jdelta + 0.5);

  int sign = 1;
  if (lperiod < 0)
    {
      std::swap(vDate0, vDate1);
      lperiod = -lperiod;
      sign = -1;
    }

  int year0, month0, day0;
  cdiDate_decode(vDate0, &year0, &month0, &day0);
  int year1, month1, day1;
  cdiDate_decode(vDate1, &year1, &month1, &day1);

  auto deltay = year1 - year0;
  auto deltam = deltay * 12 + (month1 - month0);
  if (deltay == 0) deltay = 1;
  if (deltam == 0) deltam = 1;

  TimeIncrement timeIncr;
  if (lperiod / 60 >= 1 && lperiod / 60 < 60) { timeIncr = { lperiod / 60, TimeUnits::MINUTES }; }
  else if (lperiod / 3600 >= 1 && lperiod / 3600 < 24)
    {
      timeIncr = { lperiod / 3600, TimeUnits::HOURS };
    }
  else if (lperiod / (3600 * 24) >= 1 && lperiod / (3600 * 24) < 32)
    {
      timeIncr = { lperiod / (3600 * 24), TimeUnits::DAYS };
      if (timeIncr.period >= 27 && deltam == 1) timeIncr = { 1, TimeUnits::MONTHS };
    }
  else if (lperiod / (3600 * 24 * 30) >= 1 && lperiod / (3600 * 24 * 30) < 12)
    {
      timeIncr = { deltam, TimeUnits::MONTHS };
    }
  else if (lperiod / (3600 * 24 * 30 * 12) >= 1)
    {
      timeIncr = { deltay, TimeUnits::YEARS };
    }
  else
    {
      timeIncr = { lperiod, TimeUnits::SECONDS };
    }

  timeIncr.period *= sign;

  return timeIncr;
}

void
check_time_increment(int tsID, int calendar, const CdiDateTime &vDateTime, CheckTimeIncr &checkTimeIncr)
{
  auto julianDate = julianDate_encode(calendar, vDateTime);

  if (tsID)
    {
      auto jdeltat = julianDate_to_seconds(julianDate_sub(julianDate, checkTimeIncr.julianDate0));
      auto timeIncr = get_time_increment(jdeltat, checkTimeIncr.vDate0, vDateTime.date);

      if (tsID == 1) checkTimeIncr.timeIncr = timeIncr;

      if (checkTimeIncr.printWarning
          && (timeIncr.period != checkTimeIncr.timeIncr.period || timeIncr.units != checkTimeIncr.timeIncr.units))
        {
          checkTimeIncr.printWarning = false;
          cdo_warning("Time increment in step %d (%lld%s) differs from step 1 (%lld%s)!"
                      " Set parameter equal=false for unequal time increments!",
                      tsID + 1, timeIncr.period, time_units_cstr(timeIncr.units), checkTimeIncr.timeIncr.period,
                      time_units_cstr(checkTimeIncr.timeIncr.units));
        }

      /*
      if (Options::cdoVerbose)
        fprintf(stdout, "Timestep: %d  increment: %3ld %s%s\n",
                tsID+1, (long) incrPeriod, tunits[(int)incrUnits], (std::abs(incrPeriod) != 1) ? "s" : "");
      */
    }

  checkTimeIncr.vDate0 = vDateTime.date;
  checkTimeIncr.julianDate0 = julianDate;
}

int
decode_month(const CdiDate &date)
{
  return date.month;
}

int
decode_month_and_day(const CdiDate &date)
{
  int month = date.month, day = date.day;
  return (month * 100 + day);
}

int
decode_day_of_year(const CdiDate &date)
{
  int year, mon, day;
  cdiDate_decode(date, &year, &mon, &day);

  if (day < 1 || day > 31) { return 0; }
  if (mon < 1 || mon > 12) { return 0; }

  return (mon - 1) * 31 + day;
}

int
decode_hour_of_year(const CdiDateTime &cdiDateTime)
{
  int year, month, day;
  int hour, minute, second, ms;
  cdiDate_decode(cdiDateTime.date, &year, &month, &day);
  cdiTime_decode(cdiDateTime.time, &hour, &minute, &second, &ms);

  int houroy = 0;
  if (month >= 1 && month <= 12 && day >= 1 && day <= 31 && hour >= 0 && hour < 24)
    houroy = ((month - 1) * 31 + day - 1) * 25 + hour + 1;

  return houroy;
}

int
decode_hour_of_day(const CdiDateTime &cdiDateTime)
{
  int year, month, day;
  int hour, minute, second, ms;
  cdiDate_decode(cdiDateTime.date, &year, &month, &day);
  cdiTime_decode(cdiDateTime.time, &hour, &minute, &second, &ms);

  int hourod = 0;
  if (month >= 1 && month <= 12 && day >= 1 && day <= 31 && hour >= 0 && hour < 24) hourod = hour + 1;

  return hourod;
}

void
set_date_time(CdiDateTime &datetime1, CdiDateTime datetime2)
{
  if (datetime2.date.month == 12) datetime2.date.year -= 1;

  if (cdiDate_get(datetime2.date) > cdiDate_get(datetime1.date)) datetime1 = datetime2;
}

static int
get_timeunits(const char *unitsStr, int &incrPeriod, int &incrUnits, int &timeUnits)
{
  auto len = strlen(unitsStr);

  // clang-format off
  if      (memcmp(unitsStr, "seconds", len) == 0) { incrUnits = 1;     timeUnits = TUNIT_SECOND; }
  else if (memcmp(unitsStr, "minutes", len) == 0) { incrUnits = 60;    timeUnits = TUNIT_MINUTE; }
  else if (memcmp(unitsStr, "hours", len) == 0)   { incrUnits = 3600;  timeUnits = TUNIT_HOUR; }
  else if (memcmp(unitsStr, "3hours", len) == 0)  { incrUnits = 10800; timeUnits = TUNIT_3HOURS; }
  else if (memcmp(unitsStr, "6hours", len) == 0)  { incrUnits = 21600; timeUnits = TUNIT_6HOURS; }
  else if (memcmp(unitsStr, "12hours", len) == 0) { incrUnits = 43200; timeUnits = TUNIT_12HOURS; }
  else if (memcmp(unitsStr, "days", len) == 0)    { incrUnits = 86400; timeUnits = TUNIT_DAY; }
  else if (memcmp(unitsStr, "months", len) == 0)  { incrUnits = 1;     timeUnits = TUNIT_MONTH; }
  else if (memcmp(unitsStr, "years", len) == 0)   { incrUnits = 12;    timeUnits = TUNIT_YEAR; }
  else cdo_abort("Time units >%s< unsupported!", unitsStr);

  if (timeUnits == TUNIT_HOUR)
    {
      if      (incrPeriod ==  3) { incrPeriod = 1; incrUnits = 10800; timeUnits = TUNIT_3HOURS;  }
      else if (incrPeriod ==  6) { incrPeriod = 1; incrUnits = 21600; timeUnits = TUNIT_6HOURS;  }
      else if (incrPeriod == 12) { incrPeriod = 1; incrUnits = 43200; timeUnits = TUNIT_12HOURS; }
    }
  // clang-format on

  return 0;
}

CdiDate
decode_datestring(const std::string &dateString)
{
  if (strchr(dateString.c_str() + 1, '-'))
    {
      int year = 1, month = 1, day = 1;
      std::sscanf(dateString.c_str(), "%d-%d-%d", &year, &month, &day);
      return cdiDate_encode(year, month, day);
    }
  else
    {
      return cdiDate_set(parameter_to_long(dateString));
    }
}

CdiTime
decode_timestring(const std::string &timeString)
{
  if (strchr(timeString.c_str(), ':'))
    {
      int hour = 0, minute = 0, second = 0, ms = 0;
      double fseconds = 0.0;
      std::sscanf(timeString.c_str(), "%d:%d:%lf", &hour, &minute, &fseconds);
      second = (int) fseconds;
      ms = (fseconds - second) * 1000;
      return cdiTime_encode(hour, minute, second, ms);
    }
  else
    {
      return cdiTime_set(parameter_to_int(timeString));
    }
}

void
decode_timeunits(const std::string &timeUnitsString, int &incrPeriod, int &incrUnits, int &timeUnits)
{
  incrPeriod = 0;
  incrUnits = 0;
  timeUnits = 0;

  char *pUnits = nullptr;
  auto fperiod = strtod(timeUnitsString.c_str(), &pUnits);
  if (pUnits != timeUnitsString.c_str()) incrPeriod = lround(fperiod);

  if (pUnits) get_timeunits(pUnits, incrPeriod, incrUnits, timeUnits);
}
