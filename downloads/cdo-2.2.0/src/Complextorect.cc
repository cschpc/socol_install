/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include "process_int.h"

void *
Complextorect(void *process)
{
  cdo_initialize(process);

  auto COMPLEXTORECT = cdo_operator_add("complextorect", 0, 0, nullptr);
  auto COMPLEXTOPOL = cdo_operator_add("complextopol", 0, 0, nullptr);

  auto operatorID = cdo_operator_id();

  operator_check_argc(0);

  auto streamID1 = cdo_open_read(0);

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto vlistID2 = vlistDuplicate(vlistID1);
  auto vlistID3 = vlistDuplicate(vlistID1);

  auto nvars = vlistNvars(vlistID2);
  for (int varID = 0; varID < nvars; ++varID)
    {
      auto datatype = vlistInqVarDatatype(vlistID2, varID);
      datatype = (datatype == CDI_DATATYPE_CPX64) ? CDI_DATATYPE_FLT64 : CDI_DATATYPE_FLT32;
      vlistDefVarDatatype(vlistID2, varID, datatype);
      vlistDefVarDatatype(vlistID3, varID, datatype);
    }

  auto taxisID1 = vlistInqTaxis(vlistID1);
  auto taxisID2 = taxisDuplicate(taxisID1);
  auto taxisID3 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);
  vlistDefTaxis(vlistID3, taxisID3);

  auto streamID2 = cdo_open_write(1);
  auto streamID3 = cdo_open_write(2);

  cdo_def_vlist(streamID2, vlistID2);
  cdo_def_vlist(streamID3, vlistID3);

  VarList varList1;
  varListInit(varList1, vlistID1);

  auto gridsizemax = vlistGridsizeMax(vlistID1);
  Varray<double> array1(2 * gridsizemax);
  Varray<double> array2(gridsizemax), array3(gridsizemax);

  int tsID = 0;
  while (true)
    {
      auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
      if (nrecs == 0) break;

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_taxis_copy_timestep(taxisID3, taxisID1);

      cdo_def_timestep(streamID2, tsID);
      cdo_def_timestep(streamID3, tsID);

      for (int recID = 0; recID < nrecs; ++recID)
        {
          int varID, levelID;
          cdo_inq_record(streamID1, &varID, &levelID);
          cdo_def_record(streamID2, varID, levelID);
          cdo_def_record(streamID3, varID, levelID);

          auto gridsize = varList1[varID].gridsize;

          size_t nmiss;
          cdo_read_record(streamID1, array1.data(), &nmiss);

          auto missval1 = varList1[varID].missval;
          auto missval2 = missval1;

          if (operatorID == COMPLEXTORECT)
            {
              for (size_t i = 0; i < gridsize; ++i)
                {
                  array2[i] = array1[2 * i];
                  array3[i] = array1[2 * i + 1];
                }
            }
          else if (operatorID == COMPLEXTOPOL)
            {
              for (size_t i = 0; i < gridsize; ++i)
                {
                  array2[i] = SQRTMN(ADDMN(MULMN(array1[2 * i], array1[2 * i]), MULMN(array1[2 * i + 1], array1[2 * i + 1])));
                  array3[i] = (DBL_IS_EQUAL(array1[2 * i], missval1) || DBL_IS_EQUAL(array1[2 * i + 1], missval1))
                                  ? missval1
                                  : std::atan2(array1[2 * i + 1], array1[2 * i]);
                }
            }

          cdo_write_record(streamID2, array2.data(), nmiss);
          cdo_write_record(streamID3, array3.data(), nmiss);
        }

      tsID++;
    }

  cdo_stream_close(streamID3);
  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  vlistDestroy(vlistID2);
  vlistDestroy(vlistID3);

  cdo_finish();

  return nullptr;
}
