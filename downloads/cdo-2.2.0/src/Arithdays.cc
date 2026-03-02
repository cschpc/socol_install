/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Arithdays  muldpm          Multiply with days per month
      Arithdays  divdpm          Divide by days per month
      Arithdays  muldpy          Multiply with days per year
      Arithdays  divdpy          Divide by days per year
      Arithdays  muldoy          Multiply with day of year
*/

#include <cdi.h>
#include "calendar.h"

#include "cdo_options.h"
#include "process_int.h"
#include "printinfo.h"
#include "field_functions.h"

static double
dayofyear(int calendar, const CdiDateTime &vDateTime)
{
  const int month_360[12] = { 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30 };
  const int month_365[12] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
  const int month_366[12] = { 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };

  int year, month, day;
  int hour, minute, second, ms;
  cdiDate_decode(vDateTime.date, &year, &month, &day);
  cdiTime_decode(vDateTime.time, &hour, &minute, &second, &ms);

  const auto dpy = days_per_year(calendar, year);

  double doy = 0.0;
  for (int im = 1; im < month; ++im)
    {
      auto dpm = (dpy == 360) ? month_360 : ((dpy == 365) ? month_365 : month_366);
      if (im >= 1 && im <= 12) doy += dpm[im - 1];
    }

  doy += (day - 1);
  doy += (second + minute * 60 + hour * 3600) / 86400.0;

  if (Options::cdoVerbose) cdo_print("vDateTime, dpy, doy: %s %d %g", datetime_to_string(vDateTime), dpy, doy);

  return doy;
}

void *
Arithdays(void *process)
{
  enum
  {
    Func_Month = 1,
    Func_Year,
  };

  cdo_initialize(process);

  // clang-format off
               cdo_operator_add("muldpm", FieldFunc_Mul, Func_Month, nullptr);
               cdo_operator_add("divdpm", FieldFunc_Div, Func_Month, nullptr);
               cdo_operator_add("muldpy", FieldFunc_Mul, Func_Year, nullptr);
               cdo_operator_add("divdpy", FieldFunc_Div, Func_Year, nullptr);
  int MULDOY = cdo_operator_add("muldoy", FieldFunc_Mul, 0, nullptr);
  // clang-format on

  const auto operatorID = cdo_operator_id();
  const auto operfunc = cdo_operator_f1(operatorID);
  const auto operfunc2 = cdo_operator_f2(operatorID);

  operator_check_argc(0);

  const auto streamID1 = cdo_open_read(0);

  const auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  const auto vlistID2 = vlistDuplicate(vlistID1);

  const auto taxisID1 = vlistInqTaxis(vlistID1);
  const auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  const auto calendar = taxisInqCalendar(taxisID1);

  const auto streamID2 = cdo_open_write(1);
  cdo_def_vlist(streamID2, vlistID2);

  VarList varList1;
  varListInit(varList1, vlistID1);

  Field field;

  int tsID = 0;
  while (true)
    {
      const auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
      if (nrecs == 0) break;

      const auto vDateTime = taxisInqVdatetime(taxisID1);

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_def_timestep(streamID2, tsID);

      int year, month, day;
      cdiDate_decode(vDateTime.date, &year, &month, &day);

      double rconst;
      if (operatorID == MULDOY)
        rconst = dayofyear(calendar, vDateTime);
      else
        rconst = (operfunc2 == Func_Month) ? days_per_month(calendar, year, month) : days_per_year(calendar, year);

      if (Options::cdoVerbose) cdo_print("calendar %d  year %d  month %d  result %g", calendar, year, month, rconst);

      for (int recID = 0; recID < nrecs; ++recID)
        {
          int varID, levelID;
          cdo_inq_record(streamID1, &varID, &levelID);
          field.init(varList1[varID]);
          cdo_read_record(streamID1, field);

          fieldc_function(field, rconst, operfunc);

          cdo_def_record(streamID2, varID, levelID);
          cdo_write_record(streamID2, field);
        }

      tsID++;
    }

  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  cdo_finish();

  return nullptr;
}
