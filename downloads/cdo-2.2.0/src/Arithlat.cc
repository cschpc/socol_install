/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Arithlat   mulcoslat       Multiply with cos(lat)
      Arithlat   divcoslat       Divide by cos(lat)
*/

#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include <mpim_grid.h>
#include "field_functions.h"

void *
Arithlat(void *process)
{
  cdo_initialize(process);

  cdo_operator_add("mulcoslat", FieldFunc_Mul, 0, nullptr);
  cdo_operator_add("divcoslat", FieldFunc_Div, 0, nullptr);

  const auto operatorID = cdo_operator_id();
  const auto operfunc = cdo_operator_f1(operatorID);

  operator_check_argc(0);

  const auto streamID1 = cdo_open_read(0);

  const auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  const auto vlistID2 = vlistDuplicate(vlistID1);

  const auto taxisID1 = vlistInqTaxis(vlistID1);
  const auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  const auto streamID2 = cdo_open_write(1);
  cdo_def_vlist(streamID2, vlistID2);

  VarList varList1;
  varListInit(varList1, vlistID1);

  const auto gridsizemax = vlistGridsizeMax(vlistID1);
  Varray<double> array(gridsizemax);
  Varray<double> scale;

  int gridID0 = -1;
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
          size_t nmiss;
          cdo_read_record(streamID1, array.data(), &nmiss);

          auto gridID = varList1[varID].gridID;
          const auto gridsize = varList1[varID].gridsize;
          const auto missval = varList1[varID].missval;

          if (gridID != gridID0)
            {
              gridID0 = gridID;
              gridID = generate_full_point_grid(gridID);
              if (!gridHasCoordinates(gridID)) cdo_abort("Cell center coordinates missing!");

              scale.resize(gridsize);
              gridInqYvals(gridID, scale.data());

              // Convert lat/lon units if required
              cdo_grid_to_radian(gridID, CDI_XAXIS, gridsize, &scale[0], "grid latitudes");

              if (operfunc == FieldFunc_Mul)
                for (size_t i = 0; i < gridsize; ++i) scale[i] = std::cos(scale[i]);
              else
                for (size_t i = 0; i < gridsize; ++i) scale[i] = 1.0 / std::cos(scale[i]);

              if (Options::cdoVerbose)
                for (int i = 0; i < 10; ++i) cdo_print("coslat  %3d  %g", i + 1, scale[i]);
            }

          if (nmiss)
            {
              for (size_t i = 0; i < gridsize; ++i)
                if (!DBL_IS_EQUAL(array[i], missval)) array[i] *= scale[i];
            }
          else
            {
              for (size_t i = 0; i < gridsize; ++i) array[i] *= scale[i];
            }

          cdo_def_record(streamID2, varID, levelID);
          cdo_write_record(streamID2, array.data(), nmiss);
        }

      tsID++;
    }

  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  cdo_finish();

  return nullptr;
}
