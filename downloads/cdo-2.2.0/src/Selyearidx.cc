/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Selyearidx    selyearidx         Select index of year
*/

#include <cdi.h>

#include "process_int.h"
#include "cdo_vlist.h"
#include "field_functions.h"

void *
Selyearidx(void *process)
{
  cdo_initialize(process);

  cdo_operator_add("selyearidx", 0, 0, nullptr);
  cdo_operator_add("seltimeidx", 1, 0, nullptr);

  auto ltime = cdo_operator_f1(cdo_operator_id());

  auto streamID1 = cdo_open_read(0);
  auto vlistID1 = cdo_stream_inq_vlist(streamID1);

  auto taxisID1 = vlistInqTaxis(vlistID1);

  auto streamID2 = cdo_open_read(1);
  auto vlistID2 = cdo_stream_inq_vlist(streamID2);
  auto taxisID2 = vlistInqTaxis(vlistID2);

  vlist_compare(vlistID1, vlistID2, CMP_ALL);

  auto vlistID3 = vlistDuplicate(vlistID2);
  auto taxisID3 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID3, taxisID3);

  auto streamID3 = cdo_open_write(2);
  cdo_def_vlist(streamID3, vlistID3);

  VarList varList1, varList2;
  varListInit(varList1, vlistID1);
  varListInit(varList2, vlistID2);

  auto maxrecs = vlistNrecs(vlistID1);
  std::vector<RecordInfo> recList(maxrecs);

  auto gridsizemax = vlistGridsizeMax(vlistID1);

  Varray<double> array(gridsizemax);

  FieldVector2D vars1, vars2;
  fields_from_vlist(vlistID1, vars1, FIELD_VEC);
  fields_from_vlist(vlistID2, vars2, FIELD_VEC);

  const int nvars = vlistNvars(vlistID1);
  for (int varID = 0; varID < nvars; ++varID)
    {
      const auto &var = varList1[varID];
      auto missval = varList2[varID].missval;
      for (int levelID = 0; levelID < var.nlevels; ++levelID)
        {
          for (size_t i = 0; i < var.gridsize; ++i) vars2[varID][levelID].vec_d[i] = missval;
        }
    }

  int tsID = 0;
  int tsID2 = 0;
  int tsID3 = 0;
  while (true)
    {
      auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
      if (nrecs == 0) break;

      auto year1 = taxisInqVdatetime(taxisID1).date.year;

      auto lexactdate = (gridsizemax == 1 && nrecs == 1);

      for (int recID = 0; recID < nrecs; ++recID)
        {
          int varID, levelID;
          cdo_inq_record(streamID1, &varID, &levelID);
          size_t nmiss;
          cdo_read_record(streamID1, vars1[varID][levelID].vec_d.data(), &nmiss);
          vars1[varID][levelID].nmiss = nmiss;

          if (tsID == 0) recList[recID].set(varID, levelID);
        }

      int nrecs2;
      int numSets = 0;
      while ((nrecs2 = cdo_stream_inq_timestep(streamID2, tsID2)))
        {
          auto year = taxisInqVdatetime(taxisID2).date.year;

          if (ltime == false && year1 != year) break;

          for (int recID = 0; recID < nrecs2; ++recID)
            {
              int varID, levelID;
              cdo_inq_record(streamID2, &varID, &levelID);
              size_t nmiss;
              cdo_read_record(streamID2, array.data(), &nmiss);

              const auto &var = varList2[varID];
              for (size_t i = 0; i < var.gridsize; ++i)
                if (numSets == (int) std::lround(vars1[varID][levelID].vec_d[i]))
                  {
                    if (lexactdate) cdo_taxis_copy_timestep(taxisID3, taxisID2);
                    vars2[varID][levelID].vec_d[i] = array[i];
                  }
            }

          numSets++;
          tsID2++;
        }

      if (numSets)
        {
          if (!lexactdate) cdo_taxis_copy_timestep(taxisID3, taxisID1);
          cdo_def_timestep(streamID3, tsID3);

          for (int recID = 0; recID < maxrecs; ++recID)
            {
              auto [varID, levelID] = recList[recID].get();
              if (tsID && varList1[varID].isConstant) continue;

              auto &var2 = vars2[varID][levelID];
              cdo_def_record(streamID3, varID, levelID);
              cdo_write_record(streamID3, var2.vec_d.data(), field_num_miss(var2));
            }

          tsID3++;
        }

      if (numSets == 0)
        {
          cdo_warning("First input stream has more timesteps than the second input stream!");
          break;
        }

      if (nrecs2 == 0) break;

      tsID++;
    }

  cdo_stream_close(streamID3);
  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  cdo_finish();

  return nullptr;
}
