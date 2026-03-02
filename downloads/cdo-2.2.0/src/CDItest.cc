/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:
*/

#include <cdi.h>

#include "param_conversion.h"
#include "cdo_wtime.h"
#include "process_int.h"

void *
CDItest(void *process)
{
  cdo_initialize(process);

  const auto dataIsUnchanged = false;
  // const auto dataIsUnchanged = data_is_unchanged();

  const auto NCOPY = cdo_operator_add("ncopy", 0, 0, nullptr);
  (void) (NCOPY);  // unused

  const auto operatorID = cdo_operator_id();
  (void) (operatorID);  // unused

  //  operator_input_arg("Number of copies");
  const auto max_copy = (cdo_operator_argc() == 1) ? parameter_to_int(cdo_operator_argv(0)) : 3;

  int n = 0;
  while (true)
    {
      const auto startTime = cdo_get_wtime();

      const auto streamID1 = cdo_open_read(0);

      const auto vlistID1 = cdo_stream_inq_vlist(streamID1);
      const auto taxisID1 = vlistInqTaxis(vlistID1);

      const auto vlistID2 = vlistDuplicate(vlistID1);
      const auto taxisID2 = taxisDuplicate(taxisID1);
      vlistDefTaxis(vlistID2, taxisID2);

      const auto streamID2 = cdo_open_write(1);
      cdo_def_vlist(streamID2, vlistID2);

      Field field;

      VarList varList1;
      varListInit(varList1, vlistID1);

      int tsID = 0;
      while (true)
        {
          const auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
          if (nrecs == 0) break;

          cdo_taxis_copy_timestep(taxisID2, taxisID1);

          cdo_def_timestep(streamID2, tsID);

          for (int recID = 0; recID < nrecs; ++recID)
            {
              int varID, levelID;
              cdo_inq_record(streamID1, &varID, &levelID);
              cdo_def_record(streamID2, varID, levelID);

              if (dataIsUnchanged) { cdo_copy_record(streamID2, streamID1); }
              else
                {
                  field.init(varList1[varID]);
                  cdo_read_record(streamID1, field);
                  cdo_write_record(streamID2, field);
                }
            }

          tsID++;
        }

      cdo_stream_close(streamID1);
      cdo_stream_close(streamID2);

      vlistDestroy(vlistID2);
      taxisDestroy(taxisID2);

      n++;

      cdo_print("Copy number %d: %.2fs", n, cdo_get_wtime() - startTime);

      if (n == max_copy) break;
    }

  cdo_finish();

  return nullptr;
}
