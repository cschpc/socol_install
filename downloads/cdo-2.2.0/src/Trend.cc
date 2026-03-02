/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Trend      trend           Trend
*/

#include <cdi.h>

#include "process_int.h"
#include "cdo_vlist.h"
#include "cdo_options.h"
#include "datetime.h"
#include "pmlist.h"
#include "param_conversion.h"
#include "field_functions.h"

static void
trendGetParameter(bool &tstepIsEqual)
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
          if (key == "equal") tstepIsEqual = parameter_to_bool(value);
          else cdo_abort("Invalid parameter key >%s<!", key);
          // clang-format on
        }
    }
}

void *
Trend(void *process)
{
  cdo_initialize(process);

  auto tstepIsEqual = true;
  trendGetParameter(tstepIsEqual);

  auto streamID1 = cdo_open_read(0);

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto vlistID2 = vlistDuplicate(vlistID1);

  vlistDefNtsteps(vlistID2, 1);

  auto taxisID1 = vlistInqTaxis(vlistID1);
  auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  VarList varList1;
  varListInit(varList1, vlistID1);

  auto nvars = vlistNvars(vlistID1);

  auto maxrecs = vlistNrecs(vlistID1);
  std::vector<RecordInfo> recList(maxrecs);

  for (int varID = 0; varID < nvars; ++varID) vlistDefVarDatatype(vlistID2, varID, CDI_DATATYPE_FLT64);

  auto streamID2 = cdo_open_write(1);
  auto streamID3 = cdo_open_write(2);

  cdo_def_vlist(streamID2, vlistID2);
  cdo_def_vlist(streamID3, vlistID2);

  auto gridsizemax = vlistGridsizeMax(vlistID1);

  Field field1, field2;
  field1.resize(gridsizemax);
  field2.resize(gridsizemax);

  constexpr int nwork = 5;
  FieldVector2D work[nwork];
  for (auto &w : work) fields_from_vlist(vlistID1, w, FIELD_VEC, 0);

  auto calendar = taxisInqCalendar(taxisID1);

  CheckTimeIncr checkTimeIncr;
  JulianDate julianDate0;
  double deltat1 = 0.0;
  CdiDateTime vDateTime{};
  int tsID = 0;
  while (true)
    {
      auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
      if (nrecs == 0) break;

      vDateTime = taxisInqVdatetime(taxisID1);

      if (tstepIsEqual) check_time_increment(tsID, calendar, vDateTime, checkTimeIncr);
      auto zj = tstepIsEqual ? (double) tsID : delta_time_step_0(tsID, calendar, vDateTime, julianDate0, deltat1);

      for (int recID = 0; recID < nrecs; ++recID)
        {
          int varID, levelID;
          cdo_inq_record(streamID1, &varID, &levelID);

          recList[recID].set(varID, levelID);

          cdo_read_record(streamID1, field1.vec_d.data(), &field1.nmiss);

          auto gridsize = varList1[varID].gridsize;
          auto missval = varList1[varID].missval;

          auto &sumj = work[0][varID][levelID].vec_d;
          auto &sumjj = work[1][varID][levelID].vec_d;
          auto &sumjx = work[2][varID][levelID].vec_d;
          auto &sumx = work[3][varID][levelID].vec_d;
          auto &zn = work[4][varID][levelID].vec_d;

          for (size_t i = 0; i < gridsize; ++i)
            if (!dbl_is_equal(field1.vec_d[i], missval))
              {
                sumj[i] += zj;
                sumjj[i] += zj * zj;
                sumjx[i] += zj * field1.vec_d[i];
                sumx[i] += field1.vec_d[i];
                zn[i]++;
              }
        }

      tsID++;
    }

  taxisDefVdatetime(taxisID2, vDateTime);
  cdo_def_timestep(streamID2, 0);
  cdo_def_timestep(streamID3, 0);

  for (int recID = 0; recID < maxrecs; ++recID)
    {
      auto [varID, levelID] = recList[recID].get();

      auto gridsize = varList1[varID].gridsize;
      auto missval = varList1[varID].missval;
      auto missval1 = missval;
      auto missval2 = missval;
      field1.size = gridsize;
      field1.missval = missval;
      field2.size = gridsize;
      field2.missval = missval;

      const auto &sumj = work[0][varID][levelID].vec_d;
      const auto &sumjj = work[1][varID][levelID].vec_d;
      const auto &sumjx = work[2][varID][levelID].vec_d;
      const auto &sumx = work[3][varID][levelID].vec_d;
      const auto &zn = work[4][varID][levelID].vec_d;

      for (size_t i = 0; i < gridsize; ++i)
        {
          auto temp1 = SUBMN(sumjx[i], DIVMN(MULMN(sumj[i], sumx[i]), zn[i]));
          auto temp2 = SUBMN(sumjj[i], DIVMN(MULMN(sumj[i], sumj[i]), zn[i]));

          field2.vec_d[i] = DIVMN(temp1, temp2);
          field1.vec_d[i] = SUBMN(DIVMN(sumx[i], zn[i]), MULMN(DIVMN(sumj[i], zn[i]), field2.vec_d[i]));
        }

      cdo_def_record(streamID2, varID, levelID);
      cdo_write_record(streamID2, field1.vec_d.data(), field_num_miss(field1));

      cdo_def_record(streamID3, varID, levelID);
      cdo_write_record(streamID3, field2.vec_d.data(), field_num_miss(field2));
    }

  cdo_stream_close(streamID3);
  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  cdo_finish();

  return nullptr;
}
