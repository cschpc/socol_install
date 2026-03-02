/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Cedrick Ansorge

*/

/*
  This module contains the following operators:

  Eofcoeff3d             eofcoeff3d             process eof coefficients
*/

#include "cdi.h"

#include "cdo_options.h"
#include "process_int.h"
#include "cdo_vlist.h"
#include "util_files.h"

// NO MISSING VALUE SUPPORT ADDED SO FAR

void *
Eofcoeff3d(void *process)
{
  char eof_name[16], oname[1024];
  double missval1 = -999, missval2 = -999;
  int varID, levelID;
  int nrecs;
  size_t nmiss;

  cdo_initialize(process);

  if (process_self().m_ID != 0) cdo_abort("This operator can't be combined with other operators!");

  cdo_operator_add("eofcoeff3d", 0, 0, nullptr);

  const auto streamID1 = cdo_open_read(0);
  const auto streamID2 = cdo_open_read(1);

  const auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  const auto vlistID2 = cdo_stream_inq_vlist(streamID2);

  // taxisID1 = vlistInqTaxis(vlistID1);
  const auto taxisID2 = vlistInqTaxis(vlistID2);
  const auto taxisID3 = taxisDuplicate(taxisID2);

  const auto gridID1 = vlistInqVarGrid(vlistID1, 0);
  const auto gridID2 = vlistInqVarGrid(vlistID2, 0);

  const auto gridsize = vlistGridsizeMax(vlistID1);
  if (gridsize != vlistGridsizeMax(vlistID2)) cdo_abort("Gridsize of input files does not match!");

  if (vlistNgrids(vlistID2) > 1 || vlistNgrids(vlistID1) > 1) cdo_abort("Too many different grids in input");

  const auto nvars = (vlistNvars(vlistID1) == vlistNvars(vlistID2)) ? vlistNvars(vlistID1) : -1;
  const auto nlevs = zaxisInqSize(vlistInqVarZaxis(vlistID1, 0));

  if (gridID1 != gridID2) cdo_compare_grids(gridID1, gridID2);

  strcpy(oname, cdo_get_obase().c_str());
  int nchars = strlen(oname);

  auto refname = cdo_get_stream_name(0);
  char filesuffix[32] = { 0 };
  FileUtils::gen_suffix(filesuffix, sizeof(filesuffix), cdo_inq_filetype(streamID1), vlistID1, refname);

  FieldVector3D eof(nvars);
  for (varID = 0; varID < nvars; ++varID) eof[varID].resize(nlevs);

  int eofID = 0;
  while (1)
    {
      nrecs = cdo_stream_inq_timestep(streamID1, eofID);
      if (nrecs == 0) break;

      for (int recID = 0; recID < nrecs; ++recID)
        {
          cdo_inq_record(streamID1, &varID, &levelID);
          missval1 = vlistInqVarMissval(vlistID1, varID);
          eof[varID][levelID].resize((eofID == 0) ? 1 : eofID + 1);
          eof[varID][levelID][eofID].grid = gridID1;
          eof[varID][levelID][eofID].missval = missval1;
          eof[varID][levelID][eofID].resize(gridsize);
          varray_fill(gridsize, eof[varID][levelID][eofID].vec_d, missval1);

          if (varID >= nvars) cdo_abort("Internal error - varID >= nvars");
          if (levelID >= nlevs) cdo_abort("Internal error - levelID >= nlevs");

          cdo_read_record(streamID1, eof[varID][levelID][eofID].vec_d.data(), &nmiss);
          eof[varID][levelID][eofID].nmiss = nmiss;
        }
      eofID++;
    }

  const int neof = eofID;

  if (Options::cdoVerbose) cdo_print("%s contains %i eof's", cdo_get_stream_name(0), neof);
  // Create 1x1 Grid for output
  const auto gridID3 = gridCreate(GRID_LONLAT, 1);
  gridDefXsize(gridID3, 1);
  gridDefYsize(gridID3, 1);
  double xvals = 0, yvals = 0;
  gridDefXvals(gridID3, &xvals);
  gridDefYvals(gridID3, &yvals);

  double zvals = 0.;
  const auto zaxisID3 = zaxisCreate(ZAXIS_GENERIC, 1);
  zaxisDefLevels(zaxisID3, &zvals);
  cdiDefKeyString(zaxisID3, CDI_GLOBAL, CDI_KEY_NAME, "zaxis_Reduced");
  cdiDefKeyString(zaxisID3, CDI_GLOBAL, CDI_KEY_LONGNAME,
                  "Reduced zaxis from EOF3D - only one coefficient per 3D eigenvector and time step");

  // Create var-list and time-axis for output

  const auto vlistID3 = vlistDuplicate(vlistID2);

  const auto ngrids = vlistNgrids(vlistID3);
  for (int i = 0; i < ngrids; ++i) vlistChangeGridIndex(vlistID3, i, gridID3);

  const auto nzaxis = vlistNzaxis(vlistID3);
  for (int i = 0; i < nzaxis; ++i) vlistChangeZaxisIndex(vlistID3, i, zaxisID3);

  vlistDefTaxis(vlistID3, taxisID3);
  for (varID = 0; varID < nvars; ++varID) vlistDefVarTimetype(vlistID3, varID, TIME_VARYING);

  // open streams for eofcoeff output
  std::vector<CdoStreamID> streamIDs(neof);
  for (eofID = 0; eofID < neof; eofID++)
    {
      oname[nchars] = '\0';

      sprintf(eof_name, "%5.5i", eofID);
      strcat(oname, eof_name);
      if (filesuffix[0]) strcat(oname, filesuffix);

      streamIDs[eofID] = cdo_open_write(oname);

      if (Options::cdoVerbose) cdo_print("opened %s ('w')  as stream%i for %i. eof", oname, streamIDs[eofID]->get_id(), eofID + 1);

      cdo_def_vlist(streamIDs[eofID], vlistID3);
    }

  // ALLOCATE temporary fields for data read and write
  Field in;
  in.resize(gridsize);
  in.grid = gridID1;
  FieldVector2D out(nvars);
  for (varID = 0; varID < nvars; ++varID)
    {
      out[varID].resize(neof);
      for (eofID = 0; eofID < neof; eofID++)
        {
          out[varID][eofID].missval = missval1;
          out[varID][eofID].nmiss = 0;
          out[varID][eofID].resize(1);
        }
    }

  int tsID = 0;
  while (1)
    {
      nrecs = cdo_stream_inq_timestep(streamID2, tsID);
      if (nrecs == 0) break;

      for (varID = 0; varID < nvars; ++varID)
        for (eofID = 0; eofID < neof; eofID++)
          {
            out[varID][eofID].vec_d[0] = 0;
            out[varID][eofID].grid = gridID3;
            out[varID][eofID].missval = missval2;
          }

      cdo_taxis_copy_timestep(taxisID3, taxisID2);

      for (int recID = 0; recID < nrecs; ++recID)
        {
          cdo_inq_record(streamID2, &varID, &levelID);
          cdo_read_record(streamID2, in.vec_d.data(), &in.nmiss);
          missval2 = vlistInqVarMissval(vlistID2, varID);

          for (eofID = 0; eofID < neof; eofID++)
            {
              if (recID == 0) cdo_def_timestep(streamIDs[eofID], tsID);

              nmiss = 0;
              for (size_t i = 0; i < gridsize; ++i)
                {
                  if (!DBL_IS_EQUAL(in.vec_d[i], missval2) && !DBL_IS_EQUAL(eof[varID][levelID][eofID].vec_d[i], missval1))
                    {
                      // out[varID][eofID].vec_d[0] += w[i]*in.vec_d[i]*eof[varID][levelID][eofID].vec_d[i];
                      out[varID][eofID].vec_d[0] += in.vec_d[i] * eof[varID][levelID][eofID].vec_d[i];
                    }
                  else
                    nmiss += 1;
                }
              /*
              if ( nmiss ) {
                out[varID][eofID].nmiss=1;
                out[varID][eofID].vec_d[0] = missval2;
              }
              */
            }

          if (varID >= nvars) cdo_abort("Internal error - varID >= nvars");
          if (levelID >= nlevs) cdo_abort("Internal error - levelID >= nlevs");
        }

      for (eofID = 0; eofID < neof; eofID++)
        {
          for (varID = 0; varID < nvars; ++varID)
            {
              cdo_def_record(streamIDs[eofID], varID, 0);
              cdo_write_record(streamIDs[eofID], out[varID][eofID].vec_d.data(), out[varID][eofID].nmiss);
            }
        }

      tsID++;
    }

  for (eofID = 0; eofID < neof; eofID++) cdo_stream_close(streamIDs[eofID]);

  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  cdo_finish();

  return nullptr;
}
