/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Showinfo   showparam       Show parameters
      Showinfo   showcode        Show code numbers
      Showinfo   showname        Show variable names
      Showinfo   showstdname     Show variable standard names
      Showinfo   showlevel       Show levels
      Showinfo   showyear        Show years
      Showinfo   showmon         Show months
      Showinfo   showdate        Show dates
      Showinfo   showtime        Show timesteps
      Showinfo   showltype       Show level types
      Showinfo   showformat      Show file format
*/

#include <cdi.h>

#include "process_int.h"
#include "Showattribute.h"
#include "printinfo.h"
#include "cdo_zaxis.h"

static void
print_newline(int nout, int maxOut)
{
  if (!Options::silentMode && !(nout % maxOut)) fprintf(stdout, "\n");
}

static void
print_newline_if_missing(int nout, int maxOut)
{
  if (Options::silentMode || (nout % maxOut)) fprintf(stdout, "\n");
}

static void
show_year(CdoStreamID streamID)
{
  auto vlistID = cdo_stream_inq_vlist(streamID);
  auto taxisID = vlistInqTaxis(vlistID);
  auto ntsteps = vlistNtsteps(vlistID);
  if (ntsteps == 0) return;

  constexpr int maxOut = 20;
  int nout = 0;
  int year0 = 0;

  int tsID = 0;
  while (true)
    {
      auto nrecs = cdo_stream_inq_timestep(streamID, tsID);
      if (nrecs == 0) break;

      auto vDateTime = taxisInqVdatetime(taxisID);
      const int year = vDateTime.date.year;

      if (tsID == 0 || year0 != year)
        {
          nout++;
          year0 = year;
          fprintf(stdout, " %4d", year0);
          print_newline(nout, maxOut);
        }

      tsID++;
    }
  print_newline_if_missing(nout, maxOut);
}

static void
show_mon(CdoStreamID streamID)
{
  auto vlistID = cdo_stream_inq_vlist(streamID);
  auto taxisID = vlistInqTaxis(vlistID);
  auto ntsteps = vlistNtsteps(vlistID);
  if (ntsteps == 0) return;

  constexpr int maxOut = 36;
  int nout = 0;
  int month0 = 0;

  int tsID = 0;
  while (true)
    {
      auto nrecs = cdo_stream_inq_timestep(streamID, tsID);
      if (nrecs == 0) break;

      auto vDateTime = taxisInqVdatetime(taxisID);
      const int month = vDateTime.date.month;

      if (tsID == 0 || month0 != month)
        {
          nout++;
          month0 = month;
          fprintf(stdout, " %2d", month0);
          print_newline(nout, maxOut);
        }

      tsID++;
    }
  print_newline_if_missing(nout, maxOut);
}

static void
show_date(CdoStreamID streamID)
{
  auto vlistID = cdo_stream_inq_vlist(streamID);
  auto taxisID = vlistInqTaxis(vlistID);
  auto ntsteps = vlistNtsteps(vlistID);
  if (ntsteps == 0) return;

  constexpr int maxOut = 12;
  int nout = 0;
  int64_t date0 = 0;

  int tsID = 0;
  while (true)
    {
      auto nrecs = cdo_stream_inq_timestep(streamID, tsID);
      if (nrecs == 0) break;

      auto vDateTime = taxisInqVdatetime(taxisID);
      auto vdate = cdiDate_get(vDateTime.date);

      if (tsID == 0 || date0 != vdate)
        {
          nout++;
          date0 = vdate;
          fprintf(stdout, " %s", date_to_string(vDateTime.date).c_str());
          print_newline(nout, maxOut);
        }

      tsID++;
    }
  print_newline_if_missing(nout, maxOut);
}

static void
show_time(CdoStreamID streamID)
{
  auto vlistID = cdo_stream_inq_vlist(streamID);
  auto taxisID = vlistInqTaxis(vlistID);
  auto ntsteps = vlistNtsteps(vlistID);
  if (ntsteps == 0) return;

  constexpr int maxOut = 12;
  int nout = 0;

  int tsID = 0;
  while (true)
    {
      auto nrecs = cdo_stream_inq_timestep(streamID, tsID);
      if (nrecs == 0) break;

      auto vDateTime = taxisInqVdatetime(taxisID);
      nout++;
      fprintf(stdout, " %s", time_to_string(vDateTime.time).c_str());
      print_newline(nout, maxOut);

      tsID++;
    }
  print_newline_if_missing(nout, maxOut);
}

static void
show_timestamp(CdoStreamID streamID)
{
  auto vlistID = cdo_stream_inq_vlist(streamID);
  auto taxisID = vlistInqTaxis(vlistID);
  auto ntsteps = vlistNtsteps(vlistID);
  if (ntsteps == 0) return;

  constexpr int maxOut = 4;
  int nout = 0;

  int tsID = 0;
  while (true)
    {
      auto nrecs = cdo_stream_inq_timestep(streamID, tsID);
      if (nrecs == 0) break;

      nout++;
      fprintf(stdout, " %s", datetime_to_string(taxisInqVdatetime(taxisID)).c_str());
      print_newline(nout, maxOut);

      tsID++;
    }
  print_newline_if_missing(nout, maxOut);
}

static void
show_code(const VarList &varList)
{
  constexpr int maxOut = 25;
  int nout = 0;

  int nvars = varList.size();
  for (int varID = 0; varID < nvars; ++varID)
    {
      nout++;
      const auto &var = varList[varID];
      fprintf(stdout, " %d", var.code);
      print_newline(nout, maxOut);
    }
  if (Options::silentMode || (nout % maxOut)) fprintf(stdout, "\n");
}

static void
show_grid(int vlistID)
{
  auto nvars = vlistNvars(vlistID);

  fprintf(stdout, "# param nr | grid nr | z-axis nr:   /* Use in combination with operatores: griddes and zaxisdes */\n");
  for (int varID = 0; varID < nvars; ++varID)
    {
      auto gridID = vlistInqVarGrid(vlistID, varID);
      auto zaxisID = vlistInqVarZaxis(vlistID, varID);

      fprintf(stdout, "      %3d     %3d      %3d\n", vlistInqVarCode(vlistID, varID), vlistGridIndex(vlistID, gridID) + 1,
              vlistZaxisIndex(vlistID, zaxisID) + 1);
    }
}

static void
show_unit(const VarList &varList)
{
  constexpr int maxOut = 10;
  int nout = 0;

  int nvars = varList.size();
  for (int varID = 0; varID < nvars; ++varID)
    {
      nout++;
      const auto &var = varList[varID];
      if (var.units.size()) fprintf(stdout, " %s", var.units.c_str());
      print_newline(nout, maxOut);
    }
  print_newline_if_missing(nout, maxOut);
}

static void
show_param(const VarList &varList)
{
  constexpr int maxOut = 10;
  int nout = 0;

  char paramstr[32];
  int nvars = varList.size();
  for (int varID = 0; varID < nvars; ++varID)
    {
      nout++;
      const auto &var = varList[varID];
      cdiParamToString(var.param, paramstr, sizeof(paramstr));

      fprintf(stdout, " %s", paramstr);
      print_newline(nout, maxOut);
    }
  print_newline_if_missing(nout, maxOut);
}

static void
show_name(const VarList &varList)
{
  constexpr int maxOut = 10;
  int nout = 0;

  int nvars = varList.size();
  for (int varID = 0; varID < nvars; ++varID)
    {
      nout++;
      const auto &var = varList[varID];
      fprintf(stdout, " %s", var.name.c_str());
      print_newline(nout, maxOut);
    }
  print_newline_if_missing(nout, maxOut);
}

static void
show_stdname(int vlistID)
{
  auto nvars = vlistNvars(vlistID);

  constexpr int maxOut = 1;
  int nout = 0;

  for (int varID = 0; varID < nvars; ++varID)
    {
      nout++;
      auto stdname = cdo::inq_key_string(vlistID, varID, CDI_KEY_STDNAME);
      fprintf(stdout, " %s", stdname.size() ? stdname.c_str() : "unknown");
      print_newline(nout, maxOut);
    }
  print_newline_if_missing(nout, maxOut);
}

static void
show_level(const VarList &varList)
{
  int nvars = varList.size();
  for (int varID = 0; varID < nvars; ++varID)
    {
      const auto &var = varList[varID];
      for (int levelID = 0; levelID < var.nlevels; ++levelID) fprintf(stdout, " %.9g", cdo_zaxis_inq_level(var.zaxisID, levelID));
      fprintf(stdout, "\n");
    }
}

static void
show_ltype(int vlistID)
{
  auto nzaxis = vlistNzaxis(vlistID);
  for (int index = 0; index < nzaxis; ++index)
    {
      auto zaxisID = vlistZaxis(vlistID, index);
      auto ltype = zaxis_to_ltype(zaxisID);

      if (ltype != -1) fprintf(stdout, " %d", ltype);
    }
  fprintf(stdout, "\n");
}

static void
show_atts(int vlistID, const VarList &varList)
{
  auto nvars = vlistNvars(vlistID);

  for (int varID = 0; varID < nvars; ++varID)
    {
      const auto &var = varList[varID];
      fprintf(stdout, "%s:\n", var.name.c_str());

      int nattsvar;
      cdiInqNatts(vlistID, varID, &nattsvar);
      print_attributes(varList, vlistID, varID, nattsvar, nullptr);
    }
}

void *
Showinfo(void *process)
{
  cdo_initialize(process);

  // clang-format off
  auto SHOWYEAR      = cdo_operator_add("showyear",      0, 0, nullptr);
  auto SHOWMON       = cdo_operator_add("showmon",       0, 0, nullptr);
  auto SHOWDATE      = cdo_operator_add("showdate",      0, 0, nullptr);
  auto SHOWTIME      = cdo_operator_add("showtime",      0, 0, nullptr);
  auto SHOWTIMESTAMP = cdo_operator_add("showtimestamp", 0, 0, nullptr);
  auto SHOWCODE      = cdo_operator_add("showcode",      0, 0, nullptr);
  auto SHOWUNIT      = cdo_operator_add("showunit",      0, 0, nullptr);
  auto SHOWPARAM     = cdo_operator_add("showparam",     0, 0, nullptr);
  auto SHOWNAME      = cdo_operator_add("showname",      0, 0, nullptr);
  auto SHOWSTDNAME   = cdo_operator_add("showstdname",   0, 0, nullptr);
  auto SHOWLEVEL     = cdo_operator_add("showlevel",     0, 0, nullptr);
  auto SHOWLTYPE     = cdo_operator_add("showltype",     0, 0, nullptr);
  auto SHOWFORMAT    = cdo_operator_add("showformat",    0, 0, nullptr);
  auto SHOWGRID      = cdo_operator_add("showgrid",      0, 0, nullptr); 
  auto SHOWATTS      = cdo_operator_add("showatts",      0, 0, nullptr);
  auto SHOWATTSGLOB  = cdo_operator_add("showattsglob",  0, 0, nullptr);

  auto operatorID = cdo_operator_id();

  operator_check_argc(0);

  auto streamID = cdo_open_read(0);
  auto vlistID = cdo_stream_inq_vlist(streamID);
  VarList varList;
  varListInit(varList, vlistID);

  if      (operatorID == SHOWYEAR)      show_year(streamID);
  else if (operatorID == SHOWMON)       show_mon(streamID);
  else if (operatorID == SHOWDATE)      show_date(streamID);
  else if (operatorID == SHOWTIME)      show_time(streamID);
  else if (operatorID == SHOWTIMESTAMP) show_timestamp(streamID);
  else if (operatorID == SHOWCODE)      show_code(varList);
  else if (operatorID == SHOWGRID)      show_grid(vlistID);
  else if (operatorID == SHOWUNIT)      show_unit(varList);
  else if (operatorID == SHOWPARAM)     show_param(varList);
  else if (operatorID == SHOWNAME)      show_name(varList);
  else if (operatorID == SHOWSTDNAME)   show_stdname(vlistID);
  else if (operatorID == SHOWLEVEL)     show_level(varList);
  else if (operatorID == SHOWLTYPE)     show_ltype(vlistID);
  else if (operatorID == SHOWFORMAT)
    {
      print_filetype(streamID, vlistID);
    }
  else if (operatorID == SHOWATTS || operatorID == SHOWATTSGLOB)
    {
      if (operatorID == SHOWATTS) show_atts(vlistID, varList);

      fprintf(stdout, "Global:\n");
      int natts;
      cdiInqNatts(vlistID, CDI_GLOBAL, &natts);
      print_attributes(varList, vlistID, CDI_GLOBAL, natts, nullptr);
    }
  // clang-format on

  cdo_stream_close(streamID);

  cdo_finish();

  return nullptr;
}
