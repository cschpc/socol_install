/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Diff       diff            Compare two datasets
*/

#include <map>
#include <algorithm>

#include <cdi.h>

#include "process_int.h"
#include "cdo_vlist.h"
#include "param_conversion.h"
#include "mpmo_color.h"
#include "cdo_options.h"
#include "printinfo.h"
#include "pmlist.h"
#include "cdo_zaxis.h"

struct DiffParam
{
  size_t ndiff = 0;
  double absm = 0.0;
  double relm = 0.0;
  bool dsgn = false;
  bool zero = false;
};

static inline void
diff_kernel(double v1, double v2, DiffParam &dp)
{
  auto absdiff = std::fabs(v1 - v2);
  if (absdiff > 0.0) dp.ndiff++;

  dp.absm = std::max(dp.absm, absdiff);

  auto vv = v1 * v2;
  if (vv < 0.0)
    dp.dsgn = true;
  else if (is_equal(vv, 0.0))
    dp.zero = true;
  else
    dp.relm = std::max(dp.relm, absdiff / std::max(std::fabs(v1), std::fabs(v2)));
}

static void
diff_kernel_mv(double v1, double v2, double missval1, double missval2, DiffParam &dp)
{
  auto v1IsNan = std::isnan(v1);
  auto v2IsNan = std::isnan(v2);
  auto v1IsMissval = dbl_is_equal(v1, missval1);
  auto v2IsMissval = dbl_is_equal(v2, missval2);
  if (v1IsNan != v2IsNan)
    {
      dp.ndiff++;
      dp.relm = 1.0;
    }
  else if (!v1IsMissval && !v2IsMissval) { diff_kernel(v1, v2, dp); }
  else if (v1IsMissval != v2IsMissval)
    {
      dp.ndiff++;
      dp.relm = 1.0;
    }
}

static DiffParam
diff(size_t gridsize, const Field &field1, const Field &field2)
{
  DiffParam diffParam;
  auto hasMissvals = (field1.nmiss || field2.nmiss);
  if (hasMissvals)
    {
      if (memtype_is_float_float(field1.memType, field2.memType))
        {
          for (size_t i = 0; i < gridsize; ++i)
            diff_kernel_mv(field1.vec_f[i], field2.vec_f[i], field1.missval, field2.missval, diffParam);
        }
      else if (memtype_is_float_double(field1.memType, field2.memType))
        {
          for (size_t i = 0; i < gridsize; ++i)
            diff_kernel_mv(field1.vec_f[i], field2.vec_d[i], field1.missval, field2.missval, diffParam);
        }
      else if (memtype_is_double_float(field1.memType, field2.memType))
        {
          for (size_t i = 0; i < gridsize; ++i)
            diff_kernel_mv(field1.vec_d[i], field2.vec_f[i], field1.missval, field2.missval, diffParam);
        }
      else
        {
          for (size_t i = 0; i < gridsize; ++i)
            diff_kernel_mv(field1.vec_d[i], field2.vec_d[i], field1.missval, field2.missval, diffParam);
        }
    }
  else
    {
      if (memtype_is_float_float(field1.memType, field2.memType))
        {
          for (size_t i = 0; i < gridsize; ++i) diff_kernel(field1.vec_f[i], field2.vec_f[i], diffParam);
        }
      else if (memtype_is_float_double(field1.memType, field2.memType))
        {
          for (size_t i = 0; i < gridsize; ++i) diff_kernel(field1.vec_f[i], field2.vec_d[i], diffParam);
        }
      else if (memtype_is_double_float(field1.memType, field2.memType))
        {
          for (size_t i = 0; i < gridsize; ++i) diff_kernel(field1.vec_d[i], field2.vec_f[i], diffParam);
        }
      else
        {
          for (size_t i = 0; i < gridsize; ++i) diff_kernel(field1.vec_d[i], field2.vec_d[i], diffParam);
        }
    }

  return diffParam;
}

static void
use_real_part(size_t gridsize, Field &field)
{
  if (field.memType == MemType::Float)
    for (size_t i = 0; i < gridsize; ++i) field.vec_f[i] = field.vec_f[i * 2];
  else
    for (size_t i = 0; i < gridsize; ++i) field.vec_d[i] = field.vec_d[i * 2];
}

static void
diff_get_parameter(double &abslim, double &abslim2, double &rellim, int &mapflag, int &maxcount)
{
  auto pargc = cdo_operator_argc();
  if (pargc)
    {
      const auto &pargv = cdo_get_oper_argv();

      KVList kvlist;
      kvlist.name = cdo_module_name();
      if (kvlist.parse_arguments(pargc, pargv) != 0) cdo_abort("Parse error!");
      if (Options::cdoVerbose) kvlist.print();

      for (const auto &kv : kvlist)
        {
          const auto &key = kv.key;
          if (kv.nvalues > 1) cdo_abort("Too many values for parameter key >%s<!", key);
          if (kv.nvalues < 1) cdo_abort("Missing value for parameter key >%s<!", key);
          const auto &value = kv.values[0];

          // clang-format off
          if      (key == "abslim")   abslim = parameter_to_double(value);
          else if (key == "abslim2")  abslim2 = parameter_to_double(value);
          else if (key == "rellim")   rellim = parameter_to_double(value);
          else if (key == "maxcount") maxcount = parameter_to_int(value);
          else if (key == "names")
            {
              if      (value == "left")      mapflag = 1;
              else if (value == "right")     mapflag = 2;
              else if (value == "intersect") mapflag = 3;
              else cdo_abort("Invalid value for key >%s< (names=<left/right/intersect>)", key, value);
            }
          else cdo_abort("Invalid parameter key >%s<!", key);
          // clang-format on
        }
    }
}

void *
Diff(void *process)
{
  auto printHeader = true;
  int varID1, varID2 = -1;
  int levelID;
  int ndrec = 0, nd2rec = 0, ngrec = 0;
  char paramstr[32];

  cdo_initialize(process);

  // clang-format off
  auto DIFF  = cdo_operator_add("diff",  0, 0, nullptr);
  auto DIFFP = cdo_operator_add("diffp", 0, 0, nullptr);
  auto DIFFN = cdo_operator_add("diffn", 0, 0, nullptr);
  auto DIFFC = cdo_operator_add("diffc", 0, 0, nullptr);
  // clang-format on

  auto operatorID = cdo_operator_id();

  int mapflag = 0, maxcount = 0;
  double abslim = 0.0, abslim2 = 1.e-3, rellim = 1.0;
  diff_get_parameter(abslim, abslim2, rellim, mapflag, maxcount);

  constexpr double rangeMin = -1.e33;
  constexpr double rangeMax = 1.e33;
  if (rellim < rangeMin || rellim > rangeMax) cdo_abort("Rel. limit out of range!");
  if (abslim < rangeMin || abslim > rangeMax) cdo_abort("Abs. limit out of range!");
  if (abslim2 < rangeMin || abslim2 > rangeMax) cdo_abort("Abs2. limit out of range!");

  auto streamID1 = cdo_open_read(0);
  auto streamID2 = cdo_open_read(1);

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto vlistID2 = cdo_stream_inq_vlist(streamID2);

  auto nvars = vlistNvars(vlistID1);
  std::map<int, int> mapOfVarIDs;

  if (mapflag == 0)
    {
      vlist_compare(vlistID1, vlistID2, CMP_ALL);
      for (int varID = 0; varID < nvars; ++varID) mapOfVarIDs[varID] = varID;
    }
  else { vlist_map(vlistID1, vlistID2, CMP_ALL, mapflag, mapOfVarIDs); }

  VarList varList1, varList2;
  varListInit(varList1, vlistID1);
  varListInit(varList2, vlistID2);

  Field field1, field2;

  auto taxisID = vlistInqTaxis(vlistID1);

  int nrecs, nrecs2;
  int indg = 0;
  int tsID = 0;
  while (true)
    {
      auto stopLoop = false;

      nrecs = cdo_stream_inq_timestep(streamID1, tsID);
      auto vDateTime = taxisInqVdatetime(taxisID);
      auto vdateString = date_to_string(vDateTime.date);
      auto vtimeString = time_to_string(vDateTime.time);

      nrecs2 = cdo_stream_inq_timestep(streamID2, tsID);

      if (nrecs == 0 || nrecs2 == 0) break;

      int recID2next = 0;

      for (int recID = 0; recID < nrecs; ++recID)
        {
          cdo_inq_record(streamID1, &varID1, &levelID);

          auto it = mapOfVarIDs.find(varID1);
          if (it == mapOfVarIDs.end())
            {
              if (mapflag == 2 || mapflag == 3) continue;
              cdo_abort("Internal problem (tsID=%d recID=%d): varID1=%d not found!", tsID + 1, recID + 1, varID1);
            }

          for (; recID2next < nrecs2; ++recID2next)
            {
              cdo_inq_record(streamID2, &varID2, &levelID);
              if (it->second == varID2)
                {
                  ++recID2next;
                  break;
                }
            }

          if (it->second != varID2 && recID2next == nrecs2)
            cdo_abort("Internal problem (tsID=%d recID=%d): varID2=%d not found in second stream!", tsID + 1, recID + 1,
                      it->second);

          indg += 1;

          auto gridsize = varList1[varID1].gridsize;

          // checkrel = gridInqType(gridID) != GRID_SPECTRAL;
          auto checkrel = true;

          cdiParamToString(varList1[varID1].param, paramstr, sizeof(paramstr));

          field1.init(varList1[varID1]);
          cdo_read_record(streamID1, field1);
          if (varList1[varID1].nwpv == CDI_COMP) use_real_part(gridsize, field1);

          field2.init(varList2[varID2]);
          cdo_read_record(streamID2, field2);
          if (varList2[varID2].nwpv == CDI_COMP) use_real_part(gridsize, field2);

          auto dp = diff(gridsize, field1, field2);

          if (!Options::silentMode || Options::cdoVerbose)
            {
              if (dp.absm > abslim || (checkrel && dp.relm >= rellim) || Options::cdoVerbose)
                {
                  if (printHeader)
                    {
                      printHeader = false;

                      fprintf(stdout, "               Date     Time   Level Gridsize    Miss ");
                      fprintf(stdout, "   Diff ");
                      fprintf(stdout, ": S Z  Max_Absdiff Max_Reldiff : ");

                      if (operatorID == DIFFN)
                        fprintf(stdout, "Parameter name");
                      else if (operatorID == DIFF || operatorID == DIFFP)
                        fprintf(stdout, "Parameter ID");
                      else if (operatorID == DIFFC)
                        fprintf(stdout, "Code number");

                      fprintf(stdout, "\n");
                    }

                  fprintf(stdout, "%6d ", indg);
                  fprintf(stdout, ":");

                  set_text_color(stdout, MAGENTA);
                  fprintf(stdout, "%s %s ", vdateString.c_str(), vtimeString.c_str());
                  reset_text_color(stdout);
                  set_text_color(stdout, GREEN);
                  fprintf(stdout, "%7g ", cdo_zaxis_inq_level(varList1[varID1].zaxisID, levelID));
                  fprintf(stdout, "%8zu %7zu ", gridsize, std::max(field1.nmiss, field2.nmiss));
                  fprintf(stdout, "%7zu ", dp.ndiff);
                  reset_text_color(stdout);

                  fprintf(stdout, ":");
                  fprintf(stdout, " %c %c ", dp.dsgn ? 'T' : 'F', dp.zero ? 'T' : 'F');
                  set_text_color(stdout, BLUE);
                  fprintf(stdout, "%#12.5g%#12.5g", dp.absm, dp.relm);
                  reset_text_color(stdout);
                  fprintf(stdout, " : ");

                  set_text_color(stdout, BRIGHT, GREEN);
                  if (operatorID == DIFFN)
                    fprintf(stdout, "%-11s", varList1[varID1].name.c_str());
                  else if (operatorID == DIFF || operatorID == DIFFP)
                    fprintf(stdout, "%-11s", paramstr);
                  else if (operatorID == DIFFC)
                    fprintf(stdout, "%4d", varList1[varID1].code);
                  reset_text_color(stdout);

                  fprintf(stdout, "\n");
                }
            }

          ngrec++;
          if (dp.absm > abslim || (checkrel && dp.relm >= rellim)) ndrec++;
          if (dp.absm > abslim2 || (checkrel && dp.relm >= rellim)) nd2rec++;

          if (maxcount > 0 && ndrec >= maxcount)
            {
              stopLoop = true;
              break;
            }
        }

      if (stopLoop) break;

      tsID++;
    }

  if (ndrec > 0)
    {
      Options::cdoExitStatus = 1;

      set_text_color(stdout, BRIGHT, RED);
      fprintf(stdout, "  %d of %d records differ", ndrec, ngrec);
      reset_text_color(stdout);
      fprintf(stdout, "\n");

      if (ndrec != nd2rec && abslim < abslim2) fprintf(stdout, "  %d of %d records differ more than %g\n", nd2rec, ngrec, abslim2);
      //  fprintf(stdout, "  %d of %d records differ more then one thousandth\n", nprec, ngrec);
    }

  if (nrecs == 0 && nrecs2 > 0) cdo_warning("stream2 has more time steps than stream1!");
  if (nrecs > 0 && nrecs2 == 0) cdo_warning("stream1 has more time steps than stream2!");

  cdo_stream_close(streamID1);
  cdo_stream_close(streamID2);

  cdo_finish();

  return nullptr;
}
