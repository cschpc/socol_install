/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "cdo_vlist.h"
#include "param_conversion.h"
#include "field_functions.h"

void *
Duplicate(void *process)
{
  cdo_initialize(process);

  if (cdo_operator_argc() > 1) cdo_abort("Too many arguments!");

  auto ndup = (cdo_operator_argc() == 1) ? parameter_to_int(cdo_operator_argv(0)) : 2;
  if (Options::cdoVerbose) cdo_print("ndup = %d", ndup);

  auto streamID1 = cdo_open_read(0);

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto vlistID2 = vlistDuplicate(vlistID1);

  auto taxisID1 = vlistInqTaxis(vlistID1);
  auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  VarList varList;
  varListInit(varList, vlistID1);

  auto nvars = vlistNvars(vlistID1);
  auto ntsteps = vlistNtsteps(vlistID1);

  if (ntsteps == 1 && varList_numVaryingVars(varList) == 0) ntsteps = 0;

  if (ntsteps == 0)
    {
      for (int varID = 0; varID < nvars; ++varID) vlistDefVarTimetype(vlistID2, varID, TIME_VARYING);
      ntsteps = 1;
    }

  auto numStepsOut = (ntsteps > 0) ? ndup * ntsteps : -1;
  vlistDefNtsteps(vlistID2, numStepsOut);

  auto streamID2 = cdo_open_write(1);
  cdo_def_vlist(streamID2, vlistID2);

  FieldVector3D vars;
  std::vector<CdiDateTime> vDateTimes;

  int tsID = 0;
  while (true)
    {
      auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
      if (nrecs == 0) break;

      constexpr size_t NALLOC_INC = 1024;
      if ((size_t) tsID >= vDateTimes.size()) vDateTimes.resize(vDateTimes.size() + NALLOC_INC);
      if ((size_t) tsID >= vars.size()) vars.resize(vars.size() + NALLOC_INC);

      vDateTimes[tsID] = taxisInqVdatetime(taxisID1);

      fields_from_vlist(vlistID1, vars[tsID]);

      for (int recID = 0; recID < nrecs; ++recID)
        {
          int varID, levelID;
          cdo_inq_record(streamID1, &varID, &levelID);
          auto &field = vars[tsID][varID][levelID];
          field.init(varList[varID]);
          cdo_read_record(streamID1, field);
        }

      tsID++;
    }

  auto nts = tsID;
  Field fieldOut;

  for (int idup = 0; idup < ndup; idup++)
    {
      for (tsID = 0; tsID < nts; ++tsID)
        {
          taxisDefVdatetime(taxisID2, vDateTimes[tsID]);
          cdo_def_timestep(streamID2, idup * nts + tsID);

          for (int varID = 0; varID < nvars; ++varID)
            {
              fieldOut.init(varList[varID]);
              for (int levelID = 0; levelID < varList[varID].nlevels; ++levelID)
                {
                  const auto &field = vars[tsID][varID][levelID];
                  if (field.hasData())
                    {
                      field_copy(field, fieldOut);
                      cdo_def_record(streamID2, varID, levelID);
                      cdo_write_record(streamID2, fieldOut);
                    }
                }
            }
        }
    }

  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  cdo_finish();

  return nullptr;
}
