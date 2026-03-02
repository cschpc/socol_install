/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Timcumsum    timcumsum         Cumulative sum over time
*/

#include <cdi.h>

#include "process_int.h"
#include "field_functions.h"

void *
Timcumsum(void *process)
{
  int varID, levelID;

  cdo_initialize(process);

  operator_check_argc(0);

  auto streamID1 = cdo_open_read(0);

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto vlistID2 = vlistDuplicate(vlistID1);

  auto taxisID1 = vlistInqTaxis(vlistID1);
  auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  auto streamID2 = cdo_open_write(1);
  cdo_def_vlist(streamID2, vlistID2);

  auto gridsizemax = vlistGridsizeMax(vlistID1);
  if (vlistNumber(vlistID1) != CDI_REAL) gridsizemax *= 2;

  Field field;
  field.resize(gridsizemax);

  FieldVector2D vars1;
  fields_from_vlist(vlistID1, vars1, FIELD_VEC);

  int tsID = 0;
  while (true)
    {
      auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
      if (nrecs == 0) break;

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_def_timestep(streamID2, tsID);

      for (int recID = 0; recID < nrecs; ++recID)
        {
          cdo_inq_record(streamID1, &varID, &levelID);

          auto &rvars1 = vars1[varID][levelID];

          auto fieldsize = rvars1.size;

          if (tsID == 0)
            {
              cdo_read_record(streamID1, rvars1.vec_d.data(), &rvars1.nmiss);
              if (rvars1.nmiss)
                for (size_t i = 0; i < fieldsize; ++i)
                  if (DBL_IS_EQUAL(rvars1.vec_d[i], rvars1.missval)) rvars1.vec_d[i] = 0;
            }
          else
            {
              cdo_read_record(streamID1, field.vec_d.data(), &field.nmiss);
              field.size = fieldsize;
              field.grid = rvars1.grid;
              field.missval = rvars1.missval;

              if (field.nmiss)
                for (size_t i = 0; i < fieldsize; ++i)
                  if (DBL_IS_EQUAL(field.vec_d[i], rvars1.missval)) field.vec_d[i] = 0;

              field2_sum(rvars1, field);
            }

          cdo_def_record(streamID2, varID, levelID);
          cdo_write_record(streamID2, rvars1.vec_d.data(), rvars1.nmiss);
        }

      tsID++;
    }

  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  cdo_finish();

  return nullptr;
}
