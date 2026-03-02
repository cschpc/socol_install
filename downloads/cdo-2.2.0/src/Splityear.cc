/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

     Splityear  splityear       Split in years
     Splityear  splityearmon    Split in years and month
*/

#include <climits>
#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "util_files.h"

#define MAX_YEARS 99999

void *
Splityear(void *process)
{
  CdoStreamID streamID2 = CDO_STREAM_UNDEF;

  cdo_initialize(process);

  if (process_self().m_ID != 0) cdo_abort("This operator can't be combined with other operators!");

  auto dataIsUnchanged = data_is_unchanged();

  // clang-format off
  auto SPLITYEAR    = cdo_operator_add("splityear",     0, 0, nullptr);
  auto SPLITYEARMON = cdo_operator_add("splityearmon",  0, 0, nullptr);
  // clang-format on

  auto operatorID = cdo_operator_id();

  operator_check_argc(0);

  Varray<int> cyear(MAX_YEARS, 0);

  auto streamID1 = cdo_open_read(0);

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto vlistID2 = vlistDuplicate(vlistID1);

  auto taxisID1 = vlistInqTaxis(vlistID1);
  auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  char filename[8192];
  strcpy(filename, cdo_get_obase().c_str());
  const int nchars = strlen(filename);

  auto refname = cdo_get_stream_name(0);
  char filesuffix[32] = { 0 };
  FileUtils::gen_suffix(filesuffix, sizeof(filesuffix), cdo_inq_filetype(streamID1), vlistID1, refname);

  Varray<double> array;
  // if (! dataIsUnchanged)
  {
    auto gridsizemax = vlistGridsizeMax(vlistID1);
    if (vlistNumber(vlistID1) != CDI_REAL) gridsizemax *= 2;
    array.resize(gridsizemax);
  }

  VarList varList1;
  varListInit(varList1, vlistID1);

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

  int ic = 0;
  int index1 = -INT_MAX;
  int index2;
  int year1 = -1, year2;
  int mon1 = -1, mon2;
  int tsID = 0;
  int tsID2 = 0;
  while (true)
    {
      auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
      if (nrecs == 0) break;

      cdo_taxis_copy_timestep(taxisID2, taxisID1);

      auto vDate = taxisInqVdatetime(taxisID1).date;
      int day;
      cdiDate_decode(vDate, &year2, &mon2, &day);

      if (operatorID == SPLITYEAR)
        {
          if (tsID == 0 || year1 != year2 || mon1 > mon2)
            {
              tsID2 = 0;

              ic = (year1 != year2) ? 0 : ic + 1;
              if (year2 >= 0 && year2 < MAX_YEARS)
                {
                  ic = cyear[year2];
                  cyear[year2]++;
                }

              year1 = year2;

              if (streamID2 != CDO_STREAM_UNDEF) cdo_stream_close(streamID2);

              sprintf(filename + nchars, "%04d", year1);
              if (ic > 0) sprintf(filename + strlen(filename), "_%d", ic + 1);
              if (filesuffix[0]) sprintf(filename + strlen(filename), "%s", filesuffix);

              if (Options::cdoVerbose) cdo_print("create file %s", filename);

              streamID2 = cdo_open_write(filename);
              cdo_def_vlist(streamID2, vlistID2);
            }
          mon1 = mon2;
        }
      else if (operatorID == SPLITYEARMON)
        {
          index2 = year2 * 100 + mon2;

          if (tsID == 0 || index1 != index2)
            {
              tsID2 = 0;

              index1 = index2;

              if (streamID2 != CDO_STREAM_UNDEF) cdo_stream_close(streamID2);

              sprintf(filename + nchars, "%04d", index1);
              // if ( ic > 0 ) sprintf(filename+strlen(filename), "_%d", ic+1);
              if (filesuffix[0]) sprintf(filename + strlen(filename), "%s", filesuffix);

              if (Options::cdoVerbose) cdo_print("create file %s", filename);

              streamID2 = cdo_open_write(filename);
              cdo_def_vlist(streamID2, vlistID2);
            }
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
                      auto nmiss = vars[varID][levelID].nmiss;
                      cdo_def_record(streamID2, varID, levelID);
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

      tsID2++;
      tsID++;
    }

  cdo_stream_close(streamID1);
  cdo_stream_close(streamID2);

  cdo_finish();

  return nullptr;
}
