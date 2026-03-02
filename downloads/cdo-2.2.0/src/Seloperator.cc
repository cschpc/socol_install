/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include "process_int.h"
#include "cdo_zaxis.h"
#include "param_conversion.h"
#include "cdi_lockedIO.h"

void *
Seloperator(void *process)
{
  int varID, levelID;
  bool selfound = false;

  cdo_initialize(process);

  auto dataIsUnchanged = data_is_unchanged();

  operator_input_arg("code, ltype, level");

  auto scode = parameter_to_int(cdo_operator_argv(0));
  auto sltype = parameter_to_int(cdo_operator_argv(1));

  auto slevel = (cdo_operator_argc() == 3) ? parameter_to_double(cdo_operator_argv(2)) : 0.0;

  auto streamID1 = cdo_open_read(0);
  auto vlistID1 = cdo_stream_inq_vlist(streamID1);

  auto nvars = vlistNvars(vlistID1);
  for (varID = 0; varID < nvars; ++varID)
    {
      auto code = vlistInqVarCode(vlistID1, varID);
      auto zaxisID = vlistInqVarZaxis(vlistID1, varID);
      auto nlevels = zaxisInqSize(zaxisID);
      auto ltype = zaxis_to_ltype(zaxisID);

      for (int levID = 0; levID < nlevels; levID++)
        {
          auto level = cdo_zaxis_inq_level(zaxisID, levID);
          auto sellevel = (cdo_operator_argc() == 3) ? IS_EQUAL(level, slevel) : true;
          auto selcode = (scode == -1 || scode == code);
          auto selltype = (sltype == -1 || sltype == ltype);

          if (selcode && selltype && sellevel)
            {
              vlistDefFlag(vlistID1, varID, levID, true);
              selfound = true;
            }
        }
    }

  if (!selfound) cdo_warning("Code %d, ltype %d, level %g not found!", scode, sltype, slevel);

  auto vlistID2 = vlistCreate();
  cdo_vlist_copy_flag(vlistID2, vlistID1);
  vlistDefNtsteps(vlistID2, vlistNtsteps(vlistID1));

  auto taxisID1 = vlistInqTaxis(vlistID1);
  auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  auto streamID2 = cdo_open_write(1);
  cdo_def_vlist(streamID2, vlistID2);

  Field field;

  VarList varList1;
  varListInit(varList1, vlistID1);

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
          if (vlistInqFlag(vlistID1, varID, levelID) == true)
            {
              auto varID2 = vlistFindVar(vlistID2, varID);
              auto levelID2 = vlistFindLevel(vlistID2, varID, levelID);
              cdo_def_record(streamID2, varID2, levelID2);

              if (dataIsUnchanged) { cdo_copy_record(streamID2, streamID1); }
              else
                {
                  field.init(varList1[varID]);
                  cdo_read_record(streamID1, field);
                  cdo_write_record(streamID2, field);
                }
            }
        }

      tsID++;
    }

  cdo_stream_close(streamID1);
  cdo_stream_close(streamID2);

  vlistDestroy(vlistID2);

  cdo_finish();

  return nullptr;
}
