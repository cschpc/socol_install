/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/
#ifndef DATETIME_H
#define DATETIME_H

#include <cstdio>
#include <cstdint>
#include <vector>
#include <string>

#include "julian_date.h"

enum class TimeUnits
{
  SECONDS = 0,
  MINUTES,
  HOURS,
  DAYS,
  MONTHS,
  YEARS
};

enum class TimeStat
{
  UNDEF,
  FIRST,
  LAST,
  MEAN,
  MIDHIGH,
};

struct DateTimeInfo
{
  CdiDateTime c{};     // corrected verification time
  CdiDateTime v{};     // verification time
  CdiDateTime b[2]{};  // time bounds
};

class TimeIncrement
{
public:
  int64_t period = 0;
  TimeUnits units = TimeUnits::SECONDS;

  TimeIncrement() {}
  TimeIncrement(int64_t _period, TimeUnits _units) : period(_period), units(_units) {}

  bool
  operator==(const TimeIncrement &timeIncr) const
  {
    return (period == timeIncr.period && units == timeIncr.units);
  }

  bool
  operator!=(const TimeIncrement &timeIncr) const
  {
    return (period != timeIncr.period || units != timeIncr.units);
  }
};

struct CheckTimeIncr
{
  JulianDate julianDate0;
  TimeIncrement timeIncr;
  CdiDate vDate0{};
  bool printWarning = true;
};

// clang-format off
class  // DateTimeList
#ifdef WARN_UNUSED
[[gnu::warn_unused]]
#endif
DateTimeList
// clang-format on
{
public:
  DateTimeList() { init(); }
  // clang-format off
  void set_stat(TimeStat _stat) { this->stat = _stat; }
  void set_calendar(int _calendar) { this->calendar = _calendar; }
  // clang-format on
  CdiDateTime get_vDateTime(int tsID);
  void shift();

  void taxis_set_next_timestep(int taxisID);
  void taxis_inq_timestep(int taxisID, int tsID);
  void taxis_def_timestep(int taxisID, int tsID);
  void stat_taxis_def_timestep(int taxisID, int nsteps);
  void stat_taxis_def_timestep(int taxisID);

private:
  size_t nalloc = 0;
  size_t size = 0;
  int hasBounds = -1;
  int calendar = -1;
  TimeStat stat = TimeStat::LAST;
  DateTimeInfo timestat;
  std::vector<DateTimeInfo> dtInfo;

  void init();
  void mean(int nsteps);
  void midhigh(int nsteps);
};

CdiDateTime datetime_avg(int calendar, int ndates, const std::vector<CdiDateTime> &cdiDateTimes);
void set_timestat_date(const std::string &optarg);

void adjust_month_and_year(int &month, int &year);

double delta_time_step_0(int tsID, int calendar, const CdiDateTime &vDateTime, JulianDate &juldate0, double &deltat1);

TimeIncrement get_time_increment(double jdelta, CdiDate vDate0, CdiDate vDate1);

void check_time_increment(int tsID, int calendar, const CdiDateTime &vDateTime, CheckTimeIncr &checkTimeIncr);

int decode_month(const CdiDate &date);
int decode_month_and_day(const CdiDate &date);

int decode_day_of_year(const CdiDate &date);
int decode_hour_of_year(const CdiDateTime &cdiDateTime);
int decode_hour_of_day(const CdiDateTime &cdiDateTime);

void set_date_time(CdiDateTime &datetime1, CdiDateTime datetime2);

const char *time_units_cstr(TimeUnits timeUnit);

CdiDate decode_datestring(const std::string &dateString);
CdiTime decode_timestring(const std::string &timeString);
void decode_timeunits(const std::string &timeUnitsString, int &incrPeriod, int &incrUnits, int &timeUnits);

#endif /* DATETIME_H */
