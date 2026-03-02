/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Splitdate   splitdate        Split into dates
*/

#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "param_conversion.h"
#include "util_files.h"

void *
Splitdate(void *process)
{
  cdo_initialize(process);

  if (process_self().m_ID != 0) cdo_abort("This operator can't be combined with other operators!");

  auto dataIsUnchanged = data_is_unchanged();

  cdo_operator_add("splitdate", 0, 0, nullptr);

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

  int64_t vDate0 = -1;

  CdoStreamID streamID2 = CDO_STREAM_UNDEF;
  int tsID2 = 0;
  int tsID = 0;

  while (true)
    {
      auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
      if (nrecs == 0) break;
      
      cdo_taxis_copy_timestep(taxisID2, taxisID1);

      auto vDateTime = taxisInqVdatetime(taxisID1);
      auto vDate = cdiDate_get(vDateTime.date);

      if (vDate != vDate0)
        {
          if (streamID2 != CDO_STREAM_UNDEF) cdo_stream_close(streamID2);

          vDate0 = vDate;
          tsID2 = 0;

          int year, month, day;
          cdiDate_decode(vDateTime.date, &year, &month, &day);
          sprintf(filename + nchars, "%04d-%02d-%02d", year, month, day);
          sprintf(filename + nchars + 10, "%s", filesuffix);

          if (Options::cdoVerbose) cdo_print("create file %s", filename);

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

  if (streamID2 != CDO_STREAM_UNDEF) cdo_stream_close(streamID2);


  cdo_stream_close(streamID1);

  vlistDestroy(vlistID2);

  cdo_finish();

  return nullptr;
}
