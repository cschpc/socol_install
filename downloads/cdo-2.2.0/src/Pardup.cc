/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Pardup     pardup          Duplicate parameters
      Pardup     parmul          Multiply parameters
*/

#include <cdi.h>

#include "cdo_vlist.h"
#include "process_int.h"
#include "param_conversion.h"

void *
Pardup(void *process)
{
  cdo_initialize(process);

  const auto PARDUP = cdo_operator_add("pardup", 0, 0, nullptr);
  const auto PARMUL = cdo_operator_add("parmul", 0, 0, nullptr);

  const auto operatorID = cdo_operator_id();

  int nmul = 0;
  if (operatorID == PARDUP) { nmul = 2; }
  else if (operatorID == PARMUL)
    {
      operator_input_arg("number of multiply");
      nmul = parameter_to_int(cdo_operator_argv(0));
    }
  else
    cdo_abort("operator not implemented!");

  const auto streamID1 = cdo_open_read(0);

  const auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  const auto vlistID2 = vlistDuplicate(vlistID1);

  const auto taxisID1 = vlistInqTaxis(vlistID1);
  const auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  VarList varList1;
  varListInit(varList1, vlistID1);

  const auto nvars = vlistNvars(vlistID1);

  const auto maxrecs = vlistNrecs(vlistID1);
  std::vector<RecordInfo> recList(maxrecs);

  const auto gridsizemax = vlistGridsizeMax(vlistID1);
  Varray<double> array(gridsizemax);
  Varray2D<double> vardata(nvars);
  std::vector<std::vector<size_t>> varnmiss(nvars);

  for (int varID = 0; varID < nvars; ++varID)
    {
      const auto gridsize = varList1[varID].gridsize;
      const auto nlevels = varList1[varID].nlevels;
      vardata[varID].resize(gridsize * nlevels);
      varnmiss[varID].resize(nlevels);
    }

  for (int i = 1; i < nmul; ++i)
    {
      vlistCat(vlistID2, vlistID1);
      for (int varID = 0; varID < nvars; ++varID)
        vlistDefVarParam(vlistID2, varID + nvars * i, cdiEncodeParam(-(varID + nvars * i + 1), 255, 255));
    }

  const auto streamID2 = cdo_open_write(1);
  cdo_def_vlist(streamID2, vlistID2);

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

          recList[recID].set(varID, levelID);

          const auto gridsize = varList1[varID].gridsize;
          const auto offset = gridsize * levelID;
          auto single = &vardata[varID][offset];

          size_t nmiss;
          cdo_read_record(streamID1, single, &nmiss);
          varnmiss[varID][levelID] = nmiss;
        }

      for (int i = 0; i < nmul; ++i)
        for (int recID = 0; recID < nrecs; ++recID)
          {
            auto [varID, levelID] = recList[recID].get();

            const auto varID2 = varID + i * nvars;

            const auto gridsize = varList1[varID].gridsize;
            const auto offset = gridsize * levelID;
            auto single = &vardata[varID][offset];
            const auto nmiss = varnmiss[varID][levelID];

            array_copy(gridsize, single, array.data());
            cdo_def_record(streamID2, varID2, levelID);
            cdo_write_record(streamID2, array.data(), nmiss);
          }

      tsID++;
    }

  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  cdo_finish();

  return nullptr;
}
