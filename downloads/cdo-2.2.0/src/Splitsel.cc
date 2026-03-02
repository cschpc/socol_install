/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Splitsel   splitsel        Split time selection
*/

#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "param_conversion.h"
#include "util_files.h"

void *
Splitsel(void *process)
{
  cdo_initialize(process);

  if (process_self().m_ID != 0) cdo_abort("This operator can't be combined with other operators!");

  auto dataIsUnchanged = data_is_unchanged();

  cdo_operator_add("splitsel", 0, 0, nullptr);

  // operator_input_arg("numSets <noffset <nskip>>");

  auto nargc = cdo_operator_argc();
  if (nargc < 1) cdo_abort("Too few arguments! Need %d found %d.", 1, nargc);

  auto numSets = parameter_to_int(cdo_operator_argv(0));
  auto noffset = (nargc > 1) ? parameter_to_int(cdo_operator_argv(1)) : 0;
  auto nskip = (nargc > 2) ? parameter_to_int(cdo_operator_argv(2)) : 0;

  if (Options::cdoVerbose) cdo_print("numSets = %d, noffset = %d, nskip = %d", numSets, noffset, nskip);

  if (numSets < 1) cdo_abort("numSets must be greater than 0!");
  if (noffset < 0) cdo_abort("noffset must be greater or equal 0!");
  if (nskip < 0)
    {
      if (-nskip >= numSets) cdo_abort("Absolute value of negative nskip must be less than numSets!");
      if (cdo_assert_files_only() == false) cdo_abort("Negative nskip not allowed in combination with other operators!");
    }

  auto streamID1 = cdo_open_read(0);

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto vlistID2 = vlistDuplicate(vlistID1);

  auto taxisID1 = vlistInqTaxis(vlistID1);
  // int taxisID2 = cdo_taxis_create(TAXIS_ABSOLUTE);
  auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  VarList varList1;
  varListInit(varList1, vlistID1);

  char filename[8192];
  strcpy(filename, cdo_get_obase().c_str());
  auto nchars = strlen(filename);

  auto refname = cdo_get_stream_name(0);
  char filesuffix[32] = { 0 };
  FileUtils::gen_suffix(filesuffix, sizeof(filesuffix), cdo_inq_filetype(streamID1), vlistID1, refname);

  Varray<double> array;
  //  if (! dataIsUnchanged)
  {
    auto gridsizemax = vlistGridsizeMax(vlistID1);
    if (vlistNumber(vlistID1) != CDI_REAL) gridsizemax *= 2;
    array.resize(gridsizemax);
  }

  auto haveConstVars = (varList_numConstVars(varList1) > 0);

  FieldVector2D vars;
  if (haveConstVars)
    {
      int numVars = varList1.size();
      vars.resize(numVars);

      for (int varID = 0; varID < numVars; ++varID)
        {
          const auto &var = varList1[varID];
          if (var.isConstant)
            {
              vars[varID].resize(var.nlevels);

              for (int levelID = 0; levelID < var.nlevels; ++levelID)
                {
                  vars[varID][levelID].grid = var.gridID;
                  vars[varID][levelID].resize(var.gridsize);
                }
            }
        }
    }

  int index = 1;
  int tsID;
  int nrecs = 0;

  // offset
  for (tsID = 0; tsID < noffset; ++tsID)
    {
      nrecs = cdo_stream_inq_timestep(streamID1, tsID);
      if (nrecs == 0)
        {
          cdo_warning("noffset is larger than number of timesteps!");
          goto LABEL_END;
        }

      if (tsID == 0 && haveConstVars)
        for (int recID = 0; recID < nrecs; ++recID)
          {
            int varID, levelID;
            cdo_inq_record(streamID1, &varID, &levelID);

            const auto &var = varList1[varID];
            if (var.isConstant)
              {
                size_t nmiss;
                cdo_read_record(streamID1, vars[varID][levelID].vec_d.data(), &nmiss);
                vars[varID][levelID].nmiss = nmiss;
              }
          }
    }

  while (true)
    {
      sprintf(filename + nchars, "%06d", index);
      sprintf(filename + nchars + 6, "%s", filesuffix);

      if (Options::cdoVerbose) cdo_print("create file %s", filename);

      CdoStreamID streamID2 = CDO_STREAM_UNDEF;
      int tsID2 = 0;

      for (int i = 0; i < numSets; ++i)
        {
          nrecs = cdo_stream_inq_timestep(streamID1, tsID);
          if (nrecs == 0) break;

          cdo_taxis_copy_timestep(taxisID2, taxisID1);

          if (streamID2 == CDO_STREAM_UNDEF)
            {
              streamID2 = cdo_open_write(filename);
              cdo_def_vlist(streamID2, vlistID2);
            }

          cdo_def_timestep(streamID2, tsID2);

          if (tsID > 0 && tsID2 == 0 && haveConstVars)
            {
              int numVars = varList1.size();
              for (int varID = 0; varID < numVars; ++varID)
                {
                  const auto &var = varList1[varID];
                  if (var.isConstant)
                    {
                      for (int levelID = 0; levelID < var.nlevels; ++levelID)
                        {
                          cdo_def_record(streamID2, varID, levelID);
                          auto nmiss = vars[varID][levelID].nmiss;
                          cdo_write_record(streamID2, vars[varID][levelID].vec_d.data(), nmiss);
                        }
                    }
                }
            }

          for (int recID = 0; recID < nrecs; ++recID)
            {
              int varID, levelID;
              cdo_inq_record(streamID1, &varID, &levelID);
              cdo_def_record(streamID2, varID, levelID);

              if (dataIsUnchanged && !(tsID == 0 && haveConstVars)) { cdo_copy_record(streamID2, streamID1); }
              else
                {
                  size_t nmiss;
                  cdo_read_record(streamID1, array.data(), &nmiss);
                  cdo_write_record(streamID2, array.data(), nmiss);

                  if (tsID == 0 && haveConstVars)
                    {
                      const auto &var = varList1[varID];
                      if (var.isConstant)
                        {
                          varray_copy(var.gridsize, array, vars[varID][levelID].vec_d);
                          vars[varID][levelID].nmiss = nmiss;
                        }
                    }
                }
            }

          tsID++;
          tsID2++;
        }

      cdo_stream_close(streamID2);
      if (nrecs == 0) break;

      if (cdo_stream_inq_timestep(streamID1, tsID) == 0) break;

      // skip
      for (int i = 0; i < nskip; ++i)
        if (cdo_stream_inq_timestep(streamID1, tsID + i) == 0) break;

      tsID += nskip;

      if (cdo_stream_inq_timestep(streamID1, tsID) == 0) break;

      index++;
    }

LABEL_END:

  cdo_stream_close(streamID1);

  vlistDestroy(vlistID2);

  cdo_finish();

  return nullptr;
}
