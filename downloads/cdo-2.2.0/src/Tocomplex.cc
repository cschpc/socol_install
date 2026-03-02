/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include "process_int.h"
#include "cdo_default_values.h"  // Namespace CdoDefault

void *
Tocomplex(void *process)
{
  cdo_initialize(process);

  auto RETOCOMPLEX = cdo_operator_add("retocomplex", 0, 0, nullptr);
  auto IMTOCOMPLEX = cdo_operator_add("imtocomplex", 0, 0, nullptr);

  auto operatorID = cdo_operator_id();

  operator_check_argc(0);

  auto streamID1 = cdo_open_read(0);

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto vlistID2 = vlistDuplicate(vlistID1);

  auto nvars = vlistNvars(vlistID2);
  for (int varID = 0; varID < nvars; ++varID)
    {
      auto datatype = vlistInqVarDatatype(vlistID2, varID);
      datatype = (datatype == CDI_DATATYPE_FLT64) ? CDI_DATATYPE_CPX64 : CDI_DATATYPE_CPX32;
      vlistDefVarDatatype(vlistID2, varID, datatype);
    }

  auto taxisID1 = vlistInqTaxis(vlistID1);
  auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  // if (CdoDefault::FileType != CDI_FILETYPE_EXT) cdo_abort("Complex numbers need EXTRA format; used CDO option -f ext!");
  auto streamID2 = cdo_open_write(1);
  cdo_def_vlist(streamID2, vlistID2);

  auto gridsizemax = vlistGridsizeMax(vlistID1);
  Varray<double> array1(gridsizemax), array2(2 * gridsizemax);

  VarList varList1;
  varListInit(varList1, vlistID1);

  int tsID = 0;
  int tsID2 = 0;
  while (true)
    {
      auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
      if (nrecs == 0) break;

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_def_timestep(streamID2, tsID2++);

      for (int recID = 0; recID < nrecs; ++recID)
        {
          size_t nmiss;
          int varID, levelID;
          cdo_inq_record(streamID1, &varID, &levelID);
          cdo_def_record(streamID2, varID, levelID);

          cdo_read_record(streamID1, array1.data(), &nmiss);

          auto gridsize = varList1[varID].gridsize;
          if (operatorID == RETOCOMPLEX)
            {
              for (size_t i = 0; i < gridsize; ++i)
                {
                  array2[2 * i] = array1[i];
                  array2[2 * i + 1] = 0;
                }
            }
          else if (operatorID == IMTOCOMPLEX)
            {
              for (size_t i = 0; i < gridsize; ++i)
                {
                  array2[2 * i] = 0;
                  array2[2 * i + 1] = array1[i];
                }
            }

          cdo_write_record(streamID2, array2.data(), nmiss);
        }

      tsID++;
    }

  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  vlistDestroy(vlistID2);

  cdo_finish();

  return nullptr;
}
