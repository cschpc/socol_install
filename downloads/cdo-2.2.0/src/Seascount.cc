/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Copyright (C) 2007 Brockmann Consult

  Author: Ralf Quast
          Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Seascount   seascount         Seasonal counts
*/

#include <cdi.h>

#include "cdo_season.h"
#include "datetime.h"
#include "process_int.h"
#include "field_functions.h"

void *
Seascount(void *process)
{
  CdiDateTime vDateTime0{};
  int seas0 = 0;
  int oldmon = 0;

  cdo_initialize(process);

  cdo_operator_add("seascount", 0, 0, nullptr);

  operator_check_argc(0);

  const auto seasonStart = get_season_start();

  const auto streamID1 = cdo_open_read(0);

  const auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  const auto vlistID2 = vlistDuplicate(vlistID1);

  const auto taxisID1 = vlistInqTaxis(vlistID1);
  const auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  const auto streamID2 = cdo_open_write(1);
  cdo_def_vlist(streamID2, vlistID2);

  const auto maxrecs = vlistNrecs(vlistID1);
  std::vector<RecordInfo> recList(maxrecs);

  VarList varList1;
  varListInit(varList1, vlistID1);

  Field field;

  FieldVector2D vars1;
  fields_from_vlist(vlistID1, vars1, FIELD_VEC);

  int tsID = 0;
  int otsID = 0;
  while (true)
    {
      int nrecs = 0;
      int numSets = 0;
      auto newseas = false;
      while (true)
        {
          nrecs = cdo_stream_inq_timestep(streamID1, tsID);
          if (nrecs == 0) break;

          const auto vDateTime = taxisInqVdatetime(taxisID1);

          const auto month = decode_month(vDateTime.date);
          auto newmon = month;
          if (seasonStart == SeasonStart::DEC && newmon == 12) newmon = 0;

          const auto seas = month_to_season(month);

          if (numSets == 0)
            {
              seas0 = seas;
              oldmon = newmon;
            }

          if (newmon < oldmon) newseas = true;

          if ((seas != seas0) || newseas)
            {
              cdo_add_steps(-1);
              break;
            }

          oldmon = newmon;

          for (int recID = 0; recID < nrecs; ++recID)
            {
              int varID, levelID;
              cdo_inq_record(streamID1, &varID, &levelID);
              const auto &var = varList1[varID];

              if (tsID == 0) recList[recID].set(varID, levelID);

              const auto fieldsize = vars1[varID][levelID].size;

              if (numSets == 0)
                {
                  for (size_t i = 0; i < fieldsize; ++i) vars1[varID][levelID].vec_d[i] = vars1[varID][levelID].missval;
                  vars1[varID][levelID].nmiss = fieldsize;
                }

              field.init(var);
              cdo_read_record(streamID1, field);

              field2_count(vars1[varID][levelID], field);
            }

          vDateTime0 = vDateTime;
          numSets++;
          tsID++;
        }

      if (nrecs == 0 && numSets == 0) break;

      taxisDefVdatetime(taxisID2, vDateTime0);
      cdo_def_timestep(streamID2, otsID);

      for (int recID = 0; recID < maxrecs; ++recID)
        {
          auto [varID, levelID] = recList[recID].get();
          if (otsID && varList1[varID].isConstant) continue;

          cdo_def_record(streamID2, varID, levelID);
          cdo_write_record(streamID2, vars1[varID][levelID].vec_d.data(), vars1[varID][levelID].nmiss);
        }

      if (nrecs == 0) break;
      otsID++;
    }

  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  cdo_finish();

  return nullptr;
}
