/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Change_e5slm      change_e5slm          Change ECHAM5 sea land mask
*/

#include <cdi.h>

#include "process_int.h"
#include "cdi_lockedIO.h"

void *
Change_e5slm(void *process)
{
  int nrecs;
  int varID, levelID;
  size_t nmiss;

  cdo_initialize(process);

  const auto streamID1 = cdo_open_read(0);

  const auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  const auto taxisID1 = vlistInqTaxis(vlistID1);

  const auto vlistID2 = vlistDuplicate(vlistID1);
  const auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  const auto streamID2 = cdo_open_write(1);
  cdo_def_vlist(streamID2, vlistID2);

  /* get filename of SLM */
  operator_input_arg("filename of the sea land mask");
  operator_check_argc(1);
  const char *fn_slm = cdo_operator_argv(0).c_str();

  /* read SLM */
  const auto streamIDslm = stream_open_read_locked(fn_slm);
  const auto vlistIDslm = streamInqVlist(streamIDslm);

  const auto gridsize = gridInqSize(vlistInqVarGrid(vlistIDslm, 0));

  Varray<double> array(gridsize);
  Varray<double> cland(gridsize);
  Varray<bool> lsea(gridsize);

  streamInqTimestep(streamIDslm, 0);

  streamInqRecord(streamIDslm, &varID, &levelID);
  streamReadRecord(streamIDslm, cland.data(), &nmiss);

  if (nmiss) cdo_abort("SLM with missing values are unsupported!");

  const auto mm = varray_min_max(cland);
  if (mm.min < 0 || mm.max > 1) cdo_warning("Values of SLM out of bounds! (minval=%g, maxval=%g)", mm.min, mm.max);

  streamClose(streamIDslm);

  for (size_t i = 0; i < gridsize; ++i) lsea[i] = cland[i] <= 0;

  const auto nvars = vlistNvars(vlistID1);
  std::vector<short> codes(nvars);

  for (varID = 0; varID < nvars; ++varID)
    {
      if (gridsize != gridInqSize(vlistInqVarGrid(vlistID1, varID))) cdo_abort("gridsize differ!");

      auto code = vlistInqVarCode(vlistID1, varID);
      auto varname = cdo::inq_var_name(vlistID1, varID);

      if (code < 0)
        {
          // clang-format off
          if      (varname == "SLM")       code = 172;
          else if (varname == "ALAKE")     code = 99;
          else if (varname == "WS")        code = 140;
          else if (varname == "AZ0")       code = 173;
          else if (varname == "ALB")       code = 174;
          else if (varname == "VGRAT")     code = 198;
          else if (varname == "FOREST")    code = 212;
          else if (varname == "FAO")       code = 226;
          else if (varname == "WSMX")      code = 229;
          else if (varname == "GLAC")      code = 232;
          else if (varname == "VLTCLIM")   code = 71;
          else if (varname == "VGRATCLIM") code = 70;
          // clang-format on
        }

      codes[varID] = code;
    }

  int tsID = 0;
  while ((nrecs = cdo_stream_inq_timestep(streamID1, tsID)))
    {
      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_def_timestep(streamID2, tsID);

      for (int recID = 0; recID < nrecs; ++recID)
        {
          cdo_inq_record(streamID1, &varID, &levelID);
          cdo_read_record(streamID1, array.data(), &nmiss);

          auto code = codes[varID];
          if (code == 172)
            {
              cdo_print("SLM changed!");
              for (size_t i = 0; i < gridsize; ++i) array[i] = cland[i];
            }
          else if (code == 99)
            {
              cdo_print("ALAKE set all values to zero!");
              for (size_t i = 0; i < gridsize; ++i) array[i] = 0;
            }
          else if (code == 232)
            {
              cdo_print("GLAC set sea points to %g!", array[0]);
              for (size_t i = 0; i < gridsize; ++i)
                if (cland[i] < 0.5) array[i] = array[0];
            }
          else if (code == 70 || code == 71 || code == 140 || code == 173 || code == 174 || code == 198 || code == 200
                   || code == 212 || code == 226 || code == 229)
            {
              cdo_print("Code %d set sea points to %g!", code, array[0]);
              for (size_t i = 0; i < gridsize; ++i)
                if (lsea[i]) array[i] = array[0];
            }

          cdo_def_record(streamID2, varID, levelID);
          cdo_write_record(streamID2, array.data(), nmiss);
        }

      tsID++;
    }

  cdo_stream_close(streamID1);
  cdo_stream_close(streamID2);

  cdo_finish();

  return nullptr;
}
