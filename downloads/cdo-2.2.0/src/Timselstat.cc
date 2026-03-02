/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Timselstat    timselrange        Time selection range
      Timselstat    timselmin          Time selection minimum
      Timselstat    timselmax          Time selection maximum
      Timselstat    timselsum          Time selection sum
      Timselstat    timselmean         Time selection mean
      Timselstat    timselavg          Time selection average
      Timselstat    timselvar          Time selection variance
      Timselstat    timselvar1         Time selection variance [Normalize by (n-1)]
      Timselstat    timselstd          Time selection standard deviation
      Timselstat    timselstd1         Time selection standard deviation [Normalize by (n-1)]
*/

#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "param_conversion.h"
#include "datetime.h"
#include "field_functions.h"

static void
addOperators(void)
{
  // clang-format off
  cdo_operator_add("timselrange", FieldFunc_Range, 0, nullptr);
  cdo_operator_add("timselmin",   FieldFunc_Min,   0, nullptr);
  cdo_operator_add("timselmax",   FieldFunc_Max,   0, nullptr);
  cdo_operator_add("timselsum",   FieldFunc_Sum,   0, nullptr);
  cdo_operator_add("timselmean",  FieldFunc_Mean,  0, nullptr);
  cdo_operator_add("timselavg",   FieldFunc_Avg,   0, nullptr);
  cdo_operator_add("timselvar",   FieldFunc_Var,   0, nullptr);
  cdo_operator_add("timselvar1",  FieldFunc_Var1,  0, nullptr);
  cdo_operator_add("timselstd",   FieldFunc_Std,   0, nullptr);
  cdo_operator_add("timselstd1",  FieldFunc_Std1,  0, nullptr);
  // clang-format on
}

void *
Timselstat(void *process)
{
  TimeStat timestat_date = TimeStat::MEAN;

  cdo_initialize(process);

  addOperators();

  const auto operatorID = cdo_operator_id();
  const auto operfunc = cdo_operator_f1(operatorID);

  const auto lrange = (operfunc == FieldFunc_Range);
  const auto lmean = (operfunc == FieldFunc_Mean || operfunc == FieldFunc_Avg);
  const auto lstd = (operfunc == FieldFunc_Std || operfunc == FieldFunc_Std1);
  const auto lvarstd = (lstd || operfunc == FieldFunc_Var || operfunc == FieldFunc_Var1);
  const auto lvars2 = (lvarstd || lrange);
  const int divisor = (operfunc == FieldFunc_Std1 || operfunc == FieldFunc_Var1);

  auto field2_stdvar_func = lstd ? field2_std : field2_var;
  auto fieldc_stdvar_func = lstd ? fieldc_std : fieldc_var;

  operator_input_arg("numSets <noffset <nskip>>");

  const auto nargc = cdo_operator_argc();
  const auto ndates = parameter_to_int(cdo_operator_argv(0));
  const auto noffset = (nargc > 1) ? parameter_to_int(cdo_operator_argv(1)) : 0;
  const auto nskip = (nargc > 2) ? parameter_to_int(cdo_operator_argv(2)) : 0;

  if (Options::cdoVerbose) cdo_print("numSets = %d, noffset = %d, nskip = %d", ndates, noffset, nskip);

  const auto streamID1 = cdo_open_read(0);

  const auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  const auto vlistID2 = vlistDuplicate(vlistID1);

  const auto taxisID1 = vlistInqTaxis(vlistID1);
  const auto taxisID2 = taxisDuplicate(taxisID1);
  taxisWithBounds(taxisID2);
  vlistDefTaxis(vlistID2, taxisID2);

  const auto streamID2 = cdo_open_write(1);
  cdo_def_vlist(streamID2, vlistID2);

  const auto maxrecs = vlistNrecs(vlistID1);
  std::vector<RecordInfo> recList(maxrecs);

  DateTimeList dtlist;
  dtlist.set_stat(timestat_date);
  dtlist.set_calendar(taxisInqCalendar(taxisID1));

  VarList varList;
  varListInit(varList, vlistID1);

  int VARS_MEMTYPE = 0;
  if ((operfunc == FieldFunc_Min) || (operfunc == FieldFunc_Max)) VARS_MEMTYPE = FIELD_NAT;

  Field field;

  FieldVector2D samp1, vars1, vars2;
  fields_from_vlist(vlistID1, samp1);
  fields_from_vlist(vlistID1, vars1, FIELD_VEC | VARS_MEMTYPE);
  if (lvars2) fields_from_vlist(vlistID1, vars2, FIELD_VEC);

  int tsID;
  for (tsID = 0; tsID < noffset; ++tsID)
    {
      const auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
      if (nrecs == 0) break;

      for (int recID = 0; recID < nrecs; ++recID)
        {
          int varID, levelID;
          cdo_inq_record(streamID1, &varID, &levelID);

          if (tsID == 0) recList[recID].set(varID, levelID);
        }
    }

  int otsID = 0;
  if (tsID < noffset)
    {
      cdo_warning("noffset is larger than number of timesteps!");
      goto LABEL_END;
    }

  while (true)
    {
      int nrecs = 0;
      int numSets;
      for (numSets = 0; numSets < ndates; numSets++)
        {
          nrecs = cdo_stream_inq_timestep(streamID1, tsID);
          if (nrecs == 0) break;

          dtlist.taxis_inq_timestep(taxisID1, numSets);

          for (int recID = 0; recID < nrecs; ++recID)
            {
              int varID, levelID;
              cdo_inq_record(streamID1, &varID, &levelID);

              if (tsID == 0) recList[recID].set(varID, levelID);

              auto &rsamp1 = samp1[varID][levelID];
              auto &rvars1 = vars1[varID][levelID];

              if (numSets == 0)
                {
                  cdo_read_record(streamID1, rvars1);
                  if (lrange)
                    {
                      vars2[varID][levelID].nmiss = rvars1.nmiss;
                      vars2[varID][levelID].vec_d = rvars1.vec_d;
                    }

                  if (rvars1.nmiss || !rsamp1.empty())
                    {
                      if (rsamp1.empty()) rsamp1.resize(rvars1.size);
                      field2_vinit(rsamp1, rvars1);
                    }
                }
              else
                {
                  field.init(varList[varID]);
                  cdo_read_record(streamID1, field);

                  if (field.nmiss || !rsamp1.empty())
                    {
                      if (rsamp1.empty()) rsamp1.resize(rvars1.size, numSets);
                      field2_vincr(rsamp1, field);
                    }

                  // clang-format off
                  if      (lvarstd) field2_sumsumq(rvars1, vars2[varID][levelID], field);
                  else if (lrange)  field2_maxmin(rvars1, vars2[varID][levelID], field);
                  else              field2_function(rvars1, field, operfunc);
                  // clang-format on
                }
            }

          if (numSets == 0 && lvarstd)
            for (int recID = 0; recID < maxrecs; ++recID)
              {
                auto [varID, levelID] = recList[recID].get();
                if (varList[varID].isConstant) continue;

                field2_moq(vars2[varID][levelID], vars1[varID][levelID]);
              }

          tsID++;
        }

      if (nrecs == 0 && numSets == 0) break;

      for (int recID = 0; recID < maxrecs; ++recID)
        {
          auto [varID, levelID] = recList[recID].get();
          if (varList[varID].isConstant) continue;

          const auto &rsamp1 = samp1[varID][levelID];
          auto &rvars1 = vars1[varID][levelID];

          if (lmean)
            {
              if (!rsamp1.empty())
                field2_div(rvars1, rsamp1);
              else
                fieldc_div(rvars1, (double) numSets);
            }
          else if (lvarstd)
            {
              if (!rsamp1.empty())
                field2_stdvar_func(rvars1, vars2[varID][levelID], rsamp1, divisor);
              else
                fieldc_stdvar_func(rvars1, vars2[varID][levelID], numSets, divisor);
            }
          else if (lrange) { field2_sub(rvars1, vars2[varID][levelID]); }
        }

      dtlist.stat_taxis_def_timestep(taxisID2, numSets);
      cdo_def_timestep(streamID2, otsID);

      for (int recID = 0; recID < maxrecs; ++recID)
        {
          auto [varID, levelID] = recList[recID].get();
          if (otsID && varList[varID].isConstant) continue;

          auto &rvars1 = vars1[varID][levelID];

          cdo_def_record(streamID2, varID, levelID);
          cdo_write_record(streamID2, rvars1);
        }

      if (nrecs == 0) break;
      otsID++;

      for (int i = 0; i < nskip; ++i)
        {
          nrecs = cdo_stream_inq_timestep(streamID1, tsID);
          if (nrecs == 0) break;
          tsID++;
        }

      if (nrecs == 0) break;
    }

LABEL_END:

  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  cdo_finish();

  return nullptr;
}
