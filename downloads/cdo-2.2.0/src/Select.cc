/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Select      select         Select fields
*/

#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "cdo_vlist.h"
#include "cdo_zaxis.h"
#include "datetime.h"
#include "sellist.h"
#include "printinfo.h"
#include "param_conversion.h"
#include "progress.h"
#include "cdi_lockedIO.h"

int cdo_read_timestepmask(const char *maskfile, std::vector<bool> &imask);

std::string
dom_to_string(CdiDate date)
{
  int year, month, day;
  cdiDate_decode(date, &year, &month, &day);

  const char *cmons[] = { "", "JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC" };
  if (month < 0 || month > 12) month = 0;
  char cstr[32];
  std::snprintf(cstr, sizeof(cstr), "%d%s", day, cmons[month]);

  return std::string(cstr);
}

static void
write_const_vars(CdoStreamID streamID2, const VarList &varList2, const int nvars, Varray2D<double> &vardata2)
{
  for (int varID2 = 0; varID2 < nvars; ++varID2)
    {
      if (vardata2[varID2].size())
        {
          const auto &var = varList2[varID2];
          for (int levelID2 = 0; levelID2 < var.nlevels; ++levelID2)
            {
              auto pdata = &vardata2[varID2][var.gridsize * levelID2];
              auto nmiss = array_num_mv(var.gridsize, pdata, var.missval);
              cdo_def_record(streamID2, varID2, levelID2);
              cdo_write_record(streamID2, pdata, nmiss);
            }
          vardata2[varID2].clear();
          vardata2[varID2].shrink_to_fit();
        }
    }
}

static void
eval_timestepmask(const char *maskfile, KVList &kvlist)
{
  std::vector<bool> imask;
  auto n = cdo_read_timestepmask(maskfile, imask);
  if (n < 0) cdo_abort("Read of timestep mask failed!");

  int nvalues = 0;
  for (int i = 0; i < n; ++i)
    if (imask[i]) nvalues++;
  if (nvalues == 0)
    cdo_print("timestepmask has no values!");
  else
    {
      KeyValues kv;
      kv.key = "timestep";
      kv.nvalues = nvalues;
      kv.values.resize(nvalues);

      std::vector<char> value(32);
      int j = 0;
      for (int i = 0; i < n; ++i)
        if (imask[i])
          {
            std::snprintf(value.data(), value.size(), "%d", i + 1);
            kv.values[j++] = value.data();
          }

      kvlist.push_back(kv);
    }
}

const char *
get_steptype_name(const int tsteptype)
{
  // clang-format off
  if      (tsteptype == TSTEP_INSTANT)  return "instant";
  else if (tsteptype == TSTEP_INSTANT2) return "instant";
  else if (tsteptype == TSTEP_INSTANT3) return "instant";
  else if (tsteptype == TSTEP_MIN)      return "min";
  else if (tsteptype == TSTEP_MAX)      return "max";
  else if (tsteptype == TSTEP_AVG)      return "avg";
  else if (tsteptype == TSTEP_ACCUM)    return "accum";
  else if (tsteptype == TSTEP_RANGE)    return "range";
  else if (tsteptype == TSTEP_DIFF)     return "diff";
  else if (tsteptype == TSTEP_SUM)      return "sum";
  // clang-format on
  return "unknown";
}

static bool
has_selected_params(int nvars, const VarList &varList, int vlistID)
{
  for (int varID = 0; varID < nvars; ++varID)
    {
      const auto &var = varList[varID];
      for (int levID = 0; levID < var.nlevels; ++levID)
        if (vlistInqFlag(vlistID, varID, levID) == true) return true;
    }

  return false;
}

void *
Select(void *process)
{
  auto hasConstVars = true;
  CdoStreamID streamID2 = CDO_STREAM_UNDEF;
  int nvars, nvars2 = 0;
  int last_year = -999999999;
  char paramstr[32];
  char gname[CDI_MAX_NAME];
  char zname[CDI_MAX_NAME];
  int vlistID0 = -1, vlistID2 = -1;
  int taxisID2 = CDI_UNDEFID;
  int numStepsOut = 0;
  auto doTimeSel = false;
  std::vector<bool> processVars;
  Varray2D<double> vardata2;
  Field field;
  double fstartdate = -99999999999.0;
  double fenddate = -99999999999.0;

  cdo_initialize(process);

  auto SELECT = cdo_operator_add("select", 0, 0, "parameter list");
  auto DELETE = cdo_operator_add("delete", 0, 0, "parameter list");

  auto dataIsUnchanged = data_is_unchanged();

  auto operatorID = cdo_operator_id();

  operator_input_arg(cdo_operator_enter(operatorID));

  auto argc = cdo_operator_argc();
  const auto &argnames = cdo_get_oper_argv();

  if (argc == 0) cdo_abort("Parameter missing!");

  KVList kvlist;
  kvlist.name = cdo_module_name();
  if (kvlist.parse_arguments(argc, argnames) != 0) cdo_abort("Parse error!");
  if (Options::cdoVerbose) kvlist.print();

  auto kv = kvlist.search("timestepmask");
  if (kv && kv->nvalues > 0)
    {
      if (kvlist.search("timestep")) cdo_abort("Parameter timestep and timestepmask can't be combined!");
      eval_timestepmask(kv->values[0].c_str(), kvlist);
    }

  SelectInfo selInfo(kvlist);

  // clang-format off
  SELINFO_ADD_INT(timestep_of_year, "Timestep of year");
  SELINFO_ADD_INT(timestep,         "Timestep");
  SELINFO_ADD_INT(year,             "Year");
  SELINFO_ADD_INT(month,            "Month");
  SELINFO_ADD_INT(day,              "Day");
  SELINFO_ADD_INT(hour,             "Hour");
  SELINFO_ADD_INT(minute,           "Minute");
  SELINFO_ADD_INT(code,             "Code number");
  SELINFO_ADD_INT(levidx,           "Level index");
  SELINFO_ADD_INT(ltype,            "Level type");
  SELINFO_ADD_INT(zaxisnum,         "Zaxis number");
  SELINFO_ADD_INT(gridnum,          "Grid number");
  SELINFO_ADD_FLT(level,            "Level");
  SELINFO_ADD_FLT(levrange,         "Level range");
  SELINFO_ADD_WORD(name,            "Variable name");
  SELINFO_ADD_WORD(param,           "Parameter");
  SELINFO_ADD_WORD(zaxisname,       "Zaxis name");
  SELINFO_ADD_WORD(gridname,        "Grid name");
  SELINFO_ADD_WORD(steptype,        "Time step type");
  SELINFO_ADD_WORD(startdate,       "Start date");
  SELINFO_ADD_WORD(enddate,         "End date");
  SELINFO_ADD_WORD(season,          "Season");
  SELINFO_ADD_WORD(date,            "Date");
  SELINFO_ADD_WORD(timestepmask,    "Timestep mask");
  SELINFO_ADD_WORD(dom,             "Day of month");
  // clang-format on

  if (Options::cdoVerbose) selInfo.print();

  selInfo.verify();

  if (SELINFO_NVAL(levrange) > 0 && SELINFO_NVAL(levrange) != 2) cdo_abort("Key levrange needs two values!");
  if (SELINFO_NVAL(timestepmask) > 1) cdo_abort("Key timestepmask has too many values!");
  (void) (levrange);      // unused
  (void) (timestepmask);  // unused

  auto nfiles = cdo_stream_cnt() - 1;

  DateTimeList dtlist;

  if (!Options::cdoVerbose && nfiles > 1) progress::init();

  int tsmax = -1;

  timestep = 0;
  int tsID2 = 0;
  for (int indf = 0; indf < nfiles; ++indf)
    {
      if (!Options::cdoVerbose && nfiles > 1) progress::update(0, 1, (indf + 1.) / nfiles);
      if (Options::cdoVerbose) cdo_print("Process file: %s", cdo_get_stream_name(indf));

      auto streamID1 = cdo_open_read(indf);

      auto vlistID1 = cdo_stream_inq_vlist(streamID1);
      auto taxisID1 = vlistInqTaxis(vlistID1);

      VarList varList1, varList2;
      varListInit(varList1, vlistID1);

      auto copyConstVars = false;

      if (indf == 0)
        {
          auto xresult = true;

          // vlistID0 = vlistDuplicate(vlistID1);

          vlistClearFlag(vlistID1);
          nvars = vlistNvars(vlistID1);
          processVars.resize(nvars);

          if (operatorID == DELETE)
            {
              xresult = false;
              for (int varID = 0; varID < nvars; ++varID)
                {
                  const auto &var = varList1[varID];
                  for (int levID = 0; levID < var.nlevels; ++levID) vlistDefFlag(vlistID1, varID, levID, true);
                }
            }

          auto findVariable = SELINFO_NVAL(code) || SELINFO_NVAL(name) || SELINFO_NVAL(param);

          auto doVarSel = findVariable || SELINFO_NVAL(ltype) || SELINFO_NVAL(zaxisnum) || SELINFO_NVAL(gridnum)
                          || SELINFO_NVAL(zaxisname) || SELINFO_NVAL(gridname) || SELINFO_NVAL(steptype);

          auto doLevSel = SELINFO_NVAL(level) || SELINFO_NVAL(levrange) || SELINFO_NVAL(levidx);

          doTimeSel = SELINFO_NVAL(date) || SELINFO_NVAL(startdate) || SELINFO_NVAL(enddate) || SELINFO_NVAL(season)
                      || SELINFO_NVAL(timestep_of_year) || SELINFO_NVAL(timestep) || SELINFO_NVAL(year) || SELINFO_NVAL(month)
                      || SELINFO_NVAL(day) || SELINFO_NVAL(hour) || SELINFO_NVAL(minute) || SELINFO_NVAL(dom);

          for (int varID = 0; varID < nvars; ++varID)
            {
              const auto &var = varList1[varID];

              code = var.code;

              cdiParamToString(var.param, paramstr, sizeof(paramstr));

              name = var.name.c_str();
              // stdname = var.stdname;
              param = paramstr;

              auto zaxisID = var.zaxisID;
              ltype = zaxis_to_ltype(zaxisID);

              zaxisnum = vlistZaxisIndex(vlistID1, zaxisID) + 1;
              zaxisName(zaxisInqType(zaxisID), zname);
              zaxisname = zname;

              gridnum = vlistGridIndex(vlistID1, var.gridID) + 1;
              gridName(gridInqType(var.gridID), gname);
              gridname = gname;

              steptype = var.isConstant ? "constant" : "varying";

              auto found_code = SELINFO_CHECK(code);
              auto found_name = SELINFO_CHECK(name);
              auto found_param = SELINFO_CHECK(param);
              auto found_grid = SELINFO_CHECK(gridnum);
              auto found_gname = SELINFO_CHECK(gridname);
              auto found_ltype = SELINFO_CHECK(ltype);
              auto found_zaxis = SELINFO_CHECK(zaxisnum);
              auto found_zname = SELINFO_CHECK(zaxisname);
              auto found_stype = SELINFO_CHECK(steptype);

              if (SELINFO_NVAL(steptype) && !found_stype)
                {
                  steptype = get_steptype_name(var.tsteptype);
                  found_stype = SELINFO_CHECK(steptype);
                }

              bool lstep = SELINFO_NVAL(steptype) ? found_stype : true;
              bool lvar = (found_code || found_name || found_param);
              bool lgrid = (SELINFO_NVAL(gridnum) || SELINFO_NVAL(gridname)) ? (found_grid || found_gname) : true;
              bool lvert = (SELINFO_NVAL(ltype) || SELINFO_NVAL(zaxisnum) || SELINFO_NVAL(zaxisname))
                               ? (found_ltype || found_zaxis || found_zname)
                               : true;

              processVars[varID] = (lvar && lgrid && lvert && lstep);

              if (!processVars[varID] && !lvar && !findVariable)
                {
                  if (found_grid || found_gname)
                    processVars[varID] = true;
                  else if (found_stype)
                    processVars[varID] = true;
                  else if (found_ltype || found_zaxis || found_zname)
                    processVars[varID] = true;
                  else if (!doVarSel && (SELINFO_NVAL(levidx) || SELINFO_NVAL(level) || SELINFO_NVAL(levrange)))
                    {
                      auto nlevels = var.nlevels;
                      for (int levID = 0; levID < nlevels; ++levID)
                        {
                          levidx = levID + 1;
                          level = cdo_zaxis_inq_level(zaxisID, levID);
                          if (!processVars[varID] && SELINFO_CHECK(levidx)) processVars[varID] = true;
                          if (!processVars[varID] && SELINFO_CHECK(level)) processVars[varID] = true;
                          if (!processVars[varID] && SELINFO_CHECK_RANGE(levrange, level)) processVars[varID] = true;
                        }
                    }
                }
            }

          for (int varID = 0; varID < nvars; ++varID)
            {
              if (processVars[varID])
                {
                  const auto &var = varList1[varID];
                  if (zaxisInqType(var.zaxisID) == ZAXIS_HYBRID)
                    {
                      auto psvarid = vlist_get_psvarid(vlistID1, var.zaxisID);
                      if (psvarid != -1 && !processVars[psvarid]) processVars[psvarid] = true;
                    }
                }
            }

          for (int varID = 0; varID < nvars; ++varID)
            {
              if (processVars[varID])
                {
                  const auto &var = varList1[varID];

                  for (int levID = 0; levID < var.nlevels; ++levID)
                    {
                      levidx = levID + 1;
                      level = cdo_zaxis_inq_level(var.zaxisID, levID);

                      if (var.nlevels == 1 && IS_EQUAL(level, 0))
                        {
                          SELINFO_CHECK(level);
                          vlistDefFlag(vlistID1, varID, levID, xresult);
                        }
                      else
                        {
                          if (SELINFO_NVAL(levidx))
                            {
                              if (SELINFO_CHECK(levidx)) vlistDefFlag(vlistID1, varID, levID, xresult);
                            }
                          else if (SELINFO_NVAL(level))
                            {
                              if (SELINFO_CHECK(level)) vlistDefFlag(vlistID1, varID, levID, xresult);
                            }
                          else if (SELINFO_NVAL(levrange))
                            {
                              if (SELINFO_CHECK_RANGE(levrange, level)) vlistDefFlag(vlistID1, varID, levID, xresult);
                            }
                          else { vlistDefFlag(vlistID1, varID, levID, xresult); }
                        }
                    }
                }
            }

          SELINFO_CHECK_FLAG(code);
          SELINFO_CHECK_FLAG(levidx);
          SELINFO_CHECK_FLAG(ltype);
          SELINFO_CHECK_FLAG(zaxisnum);
          SELINFO_CHECK_FLAG(gridnum);
          SELINFO_CHECK_FLAG(level);
          SELINFO_CHECK_FLAG(name);
          SELINFO_CHECK_FLAG(param);
          SELINFO_CHECK_FLAG(zaxisname);
          SELINFO_CHECK_FLAG(gridname);
          SELINFO_CHECK_FLAG(steptype);
          SELINFO_CHECK_RANGE_FLAG(levrange);

          if (has_selected_params(nvars, varList1, vlistID1))
            {
              if (doVarSel && doTimeSel)
                {
                  for (int varID = 0; varID < nvars; ++varID)
                    {
                      const auto &var = varList1[varID];
                      if (processVars[varID] && var.isConstant)
                        {
                          copyConstVars = true;
                          break;
                        }
                    }
                }
            }
          else
            {
              if ((!doVarSel) && (!doLevSel) && doTimeSel)
                {
                  copyConstVars = true;

                  for (int varID = 0; varID < nvars; ++varID)
                    {
                      processVars[varID] = true;
                      const auto &var = varList1[varID];
                      for (int levID = 0; levID < var.nlevels; ++levID) vlistDefFlag(vlistID1, varID, levID, true);
                    }
                }
              else { cdo_abort("No variable selected!"); }
            }

          // if (Options::cdoVerbose) vlistPrint(vlistID1);

          vlistID0 = vlistDuplicate(vlistID1);
          for (int varID = 0; varID < nvars; ++varID)
            {
              const auto &var = varList1[varID];
              for (int levID = 0; levID < var.nlevels; ++levID)
                vlistDefFlag(vlistID0, varID, levID, vlistInqFlag(vlistID1, varID, levID));
            }

          // if (Options::cdoVerbose) vlistPrint(vlistID0);

          vlistID2 = vlistCreate();
          cdo_vlist_copy_flag(vlistID2, vlistID0);

          varListInit(varList2, vlistID2);

          // if (Options::cdoVerbose) vlistPrint(vlistID2);

          taxisID2 = taxisDuplicate(taxisID1);

          auto ntsteps = vlistNtsteps(vlistID1);

          nvars2 = vlistNvars(vlistID2);

          if (ntsteps == 1 && nfiles == 1 && varList_numVaryingVars(varList2) == 0) ntsteps = 0;

          numStepsOut = (nfiles == 1 && doTimeSel == false) ? ntsteps : -1;
          if (operatorID == SELECT && SELINFO_NVAL(timestep) > 0) numStepsOut = SELINFO_NVAL(timestep);

          if (numStepsOut == 0 && nfiles > 1)
            {
              hasConstVars = false;
              for (int varID = 0; varID < nvars2; ++varID) vlistDefVarTimetype(vlistID2, varID, TIME_VARYING);
            }

          auto nsel = SELINFO_NVAL(timestep);

          // support for negative timestep values
          if (nsel > 0)
            {
              for (int i = 0; i < nsel; ++i)
                {
                  int ptimestep;
                  SELINFO_GET_VAL(timestep, i, &ptimestep);
                  if (ptimestep < 0)
                    {
                      if (nfiles != 1) cdo_abort("Negative timesteps only supported with one input stream!");
                      if (ntsteps < 0 && cdo_assert_files_only())
                        {
                          int tsID = 0;
                          while (cdo_stream_inq_timestep(streamID1, tsID)) tsID++;
                          ntsteps = tsID;
                          if (Options::cdoVerbose) cdo_print("Found %d timesteps", ntsteps);
                        }
                      if (ntsteps > 0)
                        {
                          if (Options::cdoVerbose) cdo_print("timestep %d changed to %d", ptimestep, ptimestep + ntsteps + 1);
                          ptimestep += ntsteps + 1;
                          SELINFO_DEF_VAL(timestep, i, &ptimestep);
                        }
                    }
                }

              for (int i = 0; i < nsel; ++i)
                {
                  int ptimestep;
                  SELINFO_GET_VAL(timestep, i, &ptimestep);
                  tsmax = std::max(tsmax, ptimestep);
                }
            }

          SELINFO_GET_VAL(startdate, 0, &startdate);
          SELINFO_GET_VAL(enddate, 0, &enddate);
          if (SELINFO_NVAL(startdate)) fstartdate = date_str_to_double(startdate, 0);
          if (SELINFO_NVAL(enddate)) fenddate = date_str_to_double(enddate, 1);
        }
      else { vlist_compare(vlistID0, vlistID1, CMP_ALL); }

      if (nvars2 == 0)
        {
          cdo_warning("No variable selected!");
          goto END_LABEL;
        }

      if (copyConstVars) vardata2.resize(nvars2);

      auto stopReading = false;
      int tsID1 = 0;
      while (true)
        {
          auto nrecs = cdo_stream_inq_timestep(streamID1, tsID1);
          if (nrecs == 0) break;

          timestep++;
          auto copyTimestep = true;

          if (doTimeSel)
            {
              copyTimestep = false;

              if (operatorID == SELECT && SELINFO_NVAL(timestep) > 0)
                {
                  if (timestep > tsmax)
                    {
                      stopReading = true;
                      break;
                    }
                }

              dtlist.taxis_inq_timestep(taxisID1, 0);
              auto vDateTime = dtlist.get_vDateTime(0);
              int second, ms;
              cdiDate_decode(vDateTime.date, &year, &month, &day);
              cdiTime_decode(vDateTime.time, &hour, &minute, &second, &ms);
              (void) (season);  // unused

              if (year != last_year)
                {
                  timestep_of_year = 0;
                  last_year = year;
                }

              timestep_of_year++;

              if (SELINFO_CHECK(timestep)) copyTimestep = true;
              if (SELINFO_CHECK(timestep_of_year)) copyTimestep = true;

              if (!copyTimestep && SELINFO_NVAL(date) == 0 && SELINFO_NVAL(timestep) == 0 && SELINFO_NVAL(timestep_of_year) == 0
                  && SELINFO_NVAL(dom) == 0)
                {
                  auto lseason = (SELINFO_NVAL(season) == 0 || SELINFO_CHECK_SEASON(season, month));
                  auto lyear = (SELINFO_NVAL(year) == 0 || SELINFO_CHECK(year));
                  auto lmonth = (SELINFO_NVAL(month) == 0 || SELINFO_CHECK(month));
                  auto lday = (SELINFO_NVAL(day) == 0 || SELINFO_CHECK(day));
                  auto lhour = (SELINFO_NVAL(hour) == 0 || SELINFO_CHECK(hour));
                  auto lminute = (SELINFO_NVAL(minute) == 0 || SELINFO_CHECK(minute));

                  if (lseason && lyear && lmonth && lday && lhour && lminute) copyTimestep = true;
                }

              auto vdate = cdiDate_get(vDateTime.date);
              auto vtime = cdiTime_get(vDateTime.time);
              const double fdate = ((double) vdate) + ((double) vtime) / 1000000.0;

              if (SELINFO_NVAL(enddate))
                {
                  copyTimestep = (fdate <= fenddate);
                  if (fdate > fenddate)
                    {
                      SELINFO_DEF_FLAG(enddate, 0, true);
                      if (operatorID == SELECT)
                        {
                          stopReading = true;
                          break;
                        }
                    }
                }

              if (SELINFO_NVAL(startdate))
                {
                  copyTimestep = (fdate >= fstartdate);
                  if (fdate >= fstartdate) SELINFO_DEF_FLAG(startdate, 0, true);
                }

              if (SELINFO_NVAL(date))
                {
                  auto datetimeString = datetime_to_string(vDateTime);
                  date = datetimeString.c_str();
                  if (SELINFO_CHECK_DATE(date)) copyTimestep = true;
                }

              if (SELINFO_NVAL(dom))
                {
                  auto domString = dom_to_string(vDateTime.date);
                  dom = domString.c_str();
                  if (SELINFO_CHECK_DATE(dom)) copyTimestep = true;
                }

              if (operatorID == DELETE) copyTimestep = !copyTimestep;

              if (copyTimestep && indf == 0 && tsID1 == 0) copyConstVars = false;
            }

          if (copyTimestep)
            {
              cdo_taxis_copy_timestep(taxisID2, taxisID1);
              if (streamID2 == CDO_STREAM_UNDEF)
                {
                  auto isLastTimestep = ((nfiles == 1) && (numStepsOut > 1) && (numStepsOut == (tsID1 + 1)));
                  if (isLastTimestep && tsID2 == 0) numStepsOut = 1;
                  vlistDefNtsteps(vlistID2, numStepsOut);
                  streamID2 = cdo_open_write(nfiles);
                  vlistDefTaxis(vlistID2, taxisID2);
                  cdo_def_vlist(streamID2, vlistID2);
                }
              cdo_def_timestep(streamID2, tsID2);

              for (int recID = 0; recID < nrecs; ++recID)
                {
                  int varID, levelID;
                  cdo_inq_record(streamID1, &varID, &levelID);
                  if (vlistInqFlag(vlistID0, varID, levelID) == true)
                    {
                      const auto &var = varList1[varID];

                      if (hasConstVars && tsID2 > 0 && tsID1 == 0)
                        if (var.isConstant) continue;

                      auto varID2 = vlistFindVar(vlistID2, varID);
                      auto levelID2 = vlistFindLevel(vlistID2, varID, levelID);
                      if (copyConstVars && tsID2 == 0) write_const_vars(streamID2, varList2, varID2, vardata2);

                      cdo_def_record(streamID2, varID2, levelID2);
                      if (dataIsUnchanged) { cdo_copy_record(streamID2, streamID1); }
                      else
                        {
                          field.init(var);
                          cdo_read_record(streamID1, field);
                          cdo_write_record(streamID2, field);
                        }
                    }
                }

              if (copyConstVars && tsID2 == 0) write_const_vars(streamID2, varList2, nvars2, vardata2);

              tsID2++;
            }
          else if (copyConstVars && indf == 0 && tsID1 == 0)
            {
              for (int recID = 0; recID < nrecs; ++recID)
                {
                  int varID, levelID;
                  cdo_inq_record(streamID1, &varID, &levelID);
                  if (vlistInqFlag(vlistID0, varID, levelID) == true)
                    {
                      auto varID2 = vlistFindVar(vlistID2, varID);
                      const auto &var = varList2[varID2];
                      if (var.isConstant)
                        {
                          auto levelID2 = vlistFindLevel(vlistID2, varID, levelID);
                          if (levelID == 0) vardata2[varID2].resize(var.gridsize * var.nlevels);
                          size_t nmiss;
                          cdo_read_record(streamID1, &vardata2[varID2][var.gridsize * levelID2], &nmiss);
                        }
                    }
                }
            }

          tsID1++;
        }

      cdo_stream_close(streamID1);

      if (stopReading) break;
    }

END_LABEL:

  if (!Options::cdoVerbose && nfiles > 1) progress::update(0, 1, 1);

  SELINFO_CHECK_FLAG(timestep_of_year);
  SELINFO_CHECK_FLAG(timestep);
  SELINFO_CHECK_FLAG(year);
  SELINFO_CHECK_FLAG(month);
  SELINFO_CHECK_FLAG(day);
  SELINFO_CHECK_FLAG(hour);
  SELINFO_CHECK_FLAG(minute);
  SELINFO_CHECK_FLAG(startdate);
  //  SELINFO_CHECK_FLAG(enddate);
  SELINFO_CHECK_FLAG(season);
  SELINFO_CHECK_FLAG(date);
  SELINFO_CHECK_FLAG(dom);

  if (streamID2 != CDO_STREAM_UNDEF) cdo_stream_close(streamID2);

  vlistDestroy(vlistID0);
  vlistDestroy(vlistID2);

  if (tsID2 == 0) cdo_abort("No timesteps selected!");

  cdo_finish();

  return nullptr;
}
