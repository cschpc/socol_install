/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Splittime  splithour       Split hours
      Splittime  splitday        Split days
      Splittime  splitmon        Split months
      Splittime  splitseas       Split seasons
*/

#include <time.h>
#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "cdo_season.h"
#include "util_files.h"

constexpr int MaxStreams = 32;

struct tm
datetime_to_tm(CdiDateTime vDateTime)
{
  int year, month, day, hour, minute, second, ms;
  cdiDate_decode(vDateTime.date, &year, &month, &day);
  cdiTime_decode(vDateTime.time, &hour, &minute, &second, &ms);

  struct tm stime;
  memset(&stime, 0, sizeof(struct tm));

  stime.tm_sec = second;
  stime.tm_min = minute;
  stime.tm_hour = hour;
  stime.tm_mday = day;
  stime.tm_mon = month - 1;
  stime.tm_year = year - 1900;

  return stime;
}

void *
Splittime(void *process)
{
  CdoStreamID streamID2;
  CdoStreamID streamIDs[MaxStreams];
  int tsIDs[MaxStreams];
  int index = 0;

  cdo_initialize(process);

  if (process_self().m_ID != 0) cdo_abort("This operator can't be combined with other operators!");

  auto dataIsUnchanged = data_is_unchanged();

  enum
  {
    func_time,
    func_date
  };

  // clang-format off
                   cdo_operator_add("splithour", func_time, 10000, nullptr);
                   cdo_operator_add("splitday",  func_date,     1, nullptr);
  auto SPLITMON  = cdo_operator_add("splitmon",  func_date,   100, nullptr);
  auto SPLITSEAS = cdo_operator_add("splitseas", func_date,   100, nullptr);
  // clang-format on

  auto operatorID = cdo_operator_id();
  auto operfunc = cdo_operator_f1(operatorID);
  auto operintval = cdo_operator_f2(operatorID);

  const char *format = nullptr;
  if (operatorID == SPLITMON && cdo_operator_argc() == 1)
    format = cdo_operator_argv(0).c_str();
  else
    operator_check_argc(0);

  auto seasonNames = get_season_name();

  for (int i = 0; i < MaxStreams; ++i) streamIDs[i] = CDO_STREAM_UNDEF;
  for (int i = 0; i < MaxStreams; ++i) tsIDs[i] = 0;

  auto streamID1 = cdo_open_read(0);

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto vlistID2 = vlistDuplicate(vlistID1);

  auto taxisID1 = vlistInqTaxis(vlistID1);
  auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  char filename[8192];
  strcpy(filename, cdo_get_obase().c_str());
  const int nchars = strlen(filename);

  auto refname = cdo_get_stream_name(0);
  char filesuffix[32] = { 0 };
  FileUtils::gen_suffix(filesuffix, sizeof(filesuffix), cdo_inq_filetype(streamID1), vlistID1, refname);

  VarList varList1;
  varListInit(varList1, vlistID1);

  Varray<double> array;
  //  if (! dataIsUnchanged)
  {
    auto gridsizemax = vlistGridsizeMax(vlistID1);
    if (vlistNumber(vlistID1) != CDI_REAL) gridsizemax *= 2;
    array.resize(gridsizemax);
  }

  auto haveConstVars = (varList_numConstVars(varList1) > 0);

  FieldVector2D vars;
  if (haveConstVars)
    {
      int numVars = varList1.size();
      vars.resize(numVars);

      for (int varID = 0; varID < numVars; ++varID)
        {
          const auto &var = varList1[varID];
          if (var.isConstant)
            {
              vars[varID].resize(var.nlevels);

              for (int levelID = 0; levelID < var.nlevels; ++levelID)
                {
                  vars[varID][levelID].grid = var.gridID;
                  vars[varID][levelID].resize(var.gridsize);
                }
            }
        }
    }

  int tsID = 0;
  while (true)
    {
      auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
      if (nrecs == 0) break;

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      auto vDateTime = taxisInqVdatetime(taxisID1);

      if (operfunc == func_date)
        {
          index = (cdiDate_get(vDateTime.date) / operintval) % 100;
          if (index < 0) index = -index;

          if (operatorID == SPLITSEAS) index = month_to_season(index);
        }
      else if (operfunc == func_time) { index = (cdiTime_get(vDateTime.time) / operintval) % 100; }

      if (index < 0 || index >= MaxStreams) cdo_abort("Index out of range!");

      streamID2 = streamIDs[index];
      if (streamID2 == CDO_STREAM_UNDEF)
        {
          if (operatorID == SPLITSEAS)
            {
              sprintf(filename + nchars, "%3s", seasonNames[index]);
              if (filesuffix[0]) sprintf(filename + nchars + 3, "%s", filesuffix);
            }
          else
            {
              char oformat[32];
              strcpy(oformat, "%02d");

              if (operatorID == SPLITMON && format)
                {
                  char sbuf[32];
                  auto stime = datetime_to_tm(vDateTime);
                  auto slen = strftime(sbuf, sizeof(sbuf), format, &stime);
                  if (slen) strcpy(oformat, sbuf);
                }

              auto slen = sprintf(filename + nchars, oformat, index);
              if (filesuffix[0]) sprintf(filename + nchars + slen, "%s", filesuffix);
            }

          if (Options::cdoVerbose) cdo_print("create file %s", filename);

          streamID2 = cdo_open_write(filename);
          cdo_def_vlist(streamID2, vlistID2);
          streamIDs[index] = streamID2;
        }

      cdo_def_timestep(streamID2, tsIDs[index]);

      if (tsID > 0 && tsIDs[index] == 0 && haveConstVars)
        {
          int numVars = varList1.size();
          for (int varID = 0; varID < numVars; ++varID)
            {
              const auto &var = varList1[varID];
              if (var.isConstant)
                {
                  for (int levelID = 0; levelID < var.nlevels; ++levelID)
                    {
                      cdo_def_record(streamID2, varID, levelID);
                      auto nmiss = vars[varID][levelID].nmiss;
                      cdo_write_record(streamID2, vars[varID][levelID].vec_d.data(), nmiss);
                    }
                }
            }
        }

      for (int recID = 0; recID < nrecs; ++recID)
        {
          int varID, levelID;
          cdo_inq_record(streamID1, &varID, &levelID);
          cdo_def_record(streamID2, varID, levelID);

          if (dataIsUnchanged && !(tsID == 0 && haveConstVars)) { cdo_copy_record(streamID2, streamID1); }
          else
            {
              size_t nmiss;
              cdo_read_record(streamID1, array.data(), &nmiss);
              cdo_write_record(streamID2, array.data(), nmiss);

              if (tsID == 0 && haveConstVars)
                {
                  const auto &var = varList1[varID];
                  if (var.isConstant)
                    {
                      varray_copy(var.gridsize, array, vars[varID][levelID].vec_d);
                      vars[varID][levelID].nmiss = nmiss;
                    }
                }
            }
        }

      tsIDs[index]++;
      tsID++;
    }

  cdo_stream_close(streamID1);

  for (index = 0; index < MaxStreams; ++index)
    {
      if (streamIDs[index] != CDO_STREAM_UNDEF) cdo_stream_close(streamIDs[index]);
    }

  vlistDestroy(vlistID2);

  cdo_finish();

  return nullptr;
}
