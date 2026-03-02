/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Ralf Müller
          Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Intlevel   intlevel3d      Linear level interpolation on a 3d vertical coordinates variable
      Intlevel   intlevelx3d     Linear level interpolation on a 3d vertical coordinates variable with extrapolation
*/

#include <cdi.h>

#include "cimdOmp.h"
#include "cdo_options.h"
#include "process_int.h"
#include "cdo_vlist.h"
#include "cdi_lockedIO.h"

void vert_interp_lev3d(size_t gridsize, int nlev1, double missval, const Field3D &field1, Field3D &field2, int nlev2,
                       const Varray<int> &lev_idx, const Varray<float> &lev_wgt);
void vert_gen_weights(int expol, int nlev1, const Varray<double> &lev1, int nlev2, const Varray<double> &lev2, Varray<int> &lev_idx,
                      Varray<float> &lev_wgt);
bool levelDirUp(const int nlev, const double *const lev);
bool levelDirDown(const int nlev, const double *const lev);

/*
 * Create weights for the 3d vertical coordinate
 *
 * The resulting index sets lev_idx1 and lev_idx2 contain absolute numbers,i.e.
 * wrt. the given gridsize. They can directly be used to read values from 3d data fields.
 *
 * 3d version of vert_gen_weights() (src/Intlevel.cc)
 */
static void
vert_gen_weights3d(bool expol, size_t gridsize, int nlev1, const Varray<float> &xlev1, int nlev2, const Varray<float> &xlev2,
                   Varray<int> &xlev_idx, Varray<float> &xlev_wgt)
{
  auto nthreads = Threading::ompNumThreads;
  Varray2D<double> lev1p2(nthreads, Varray<double>(nlev1 + 2));
  Varray2D<double> lev2(nthreads, Varray<double>(nlev2));
  Varray2D<float> lev_wgt(nthreads, Varray<float>(nlev2));
  Varray2D<int> lev_idx(nthreads, Varray<int>(nlev2));

  // Check monotony of vertical levels
  for (int k = 0; k < nlev1; ++k) lev1p2[0][k] = xlev1[k * gridsize];
  auto lup = levelDirUp(nlev1, lev1p2[0].data());
  auto ldown = levelDirDown(nlev1, lev1p2[0].data());
  if (!lup && !ldown) cdo_abort("Non monotonic zaxis!");
  double level_0 = lup ? -1.e33 : 1.e33;
  double level_N = lup ? 1.e33 : -1.e33;

#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
  for (size_t i = 0; i < gridsize; ++i)
    {
      auto ompthID = cdo_omp_get_thread_num();

      lev1p2[ompthID][0] = level_0;
      lev1p2[ompthID][nlev1 + 1] = level_N;
      for (int k = 0; k < nlev1; ++k) lev1p2[ompthID][k + 1] = xlev1[k * gridsize + i];
      for (int k = 0; k < nlev2; ++k) lev2[ompthID][k] = xlev2[k * gridsize + i];

      vert_gen_weights(expol, nlev1 + 2, lev1p2[ompthID], nlev2, lev2[ompthID], lev_idx[ompthID], lev_wgt[ompthID]);

      for (int k = 0; k < nlev2; ++k) xlev_idx[k * gridsize + i] = lev_idx[ompthID][k];
      for (int k = 0; k < nlev2; ++k) xlev_wgt[k * gridsize + i] = lev_wgt[ompthID][k];
    }
}

static void
vlist_copy_var_attributes(const VarList &varList0, int varID0, int vlistID3, int oz3dvarID)
{
  const auto &var0 = varList0[varID0];
  cdiDefKeyString(vlistID3, oz3dvarID, CDI_KEY_NAME, var0.name.c_str());
  if (var0.longname.size()) cdiDefKeyString(vlistID3, oz3dvarID, CDI_KEY_LONGNAME, var0.longname.c_str());
  if (var0.units.size()) cdiDefKeyString(vlistID3, oz3dvarID, CDI_KEY_UNITS, var0.units.c_str());
}

static void
read_source_coordinate(int streamNumber, VarList &varList2, Varray<float> &zlevelsIn)
{
  auto streamID2 = cdo_open_read(streamNumber);  // 3d vertical source coordinate
  auto vlistID2 = cdo_stream_inq_vlist(streamID2);

  varListInit(varList2, vlistID2);

  if (vlistNvars(vlistID2) != 1) cdo_abort("infile2: Only one single variable is allowed!");

  auto gridsize = varList2[0].gridsize;
  auto nlevels = varList2[0].nlevels;

  zlevelsIn.resize(gridsize * nlevels);

  auto nrecs = cdo_stream_inq_timestep(streamID2, 0);
  if (Options::cdoVerbose) cdo_print("%d records input 3d vertical height", nrecs);

  for (int recID = 0; recID < nrecs; ++recID)
    {
      int varID, levelID;
      cdo_inq_record(streamID2, &varID, &levelID);
      auto offset = gridsize * levelID;
      size_t nmiss;
      cdo_read_record_f(streamID2, &zlevelsIn[offset], &nmiss);
      if (0 != nmiss) cdo_abort("Input vertical coordinate variables are not allowed to contain missing values.");
    }

  cdo_stream_close(streamID2);
}

static void
read_target_coordinate(const std::string &fileName, VarList &varList0, Varray<float> &zlevelsOut)
{
  auto streamID0 = stream_open_read_locked(fileName.c_str());  // 3d vertical target coordinate
  auto vlistID0 = streamInqVlist(streamID0);

  varListInit(varList0, vlistID0);

  if (vlistNvars(vlistID0) != 1) cdo_abort("tgtcoordinate: Only one single variable is allowed!");

  auto gridsize = varList0[0].gridsize;
  auto nlevels = varList0[0].nlevels;

  zlevelsOut.resize(gridsize * nlevels);

  auto nrecs = streamInqTimestep(streamID0, 0);
  if (Options::cdoVerbose) cdo_print("%d records target 3d vertical height and gridsize %zu", nrecs, gridsize);

  for (int recID = 0; recID < nrecs; ++recID)
    {
      int varID, levelID;
      streamInqRecord(streamID0, &varID, &levelID);
      auto offset = gridsize * levelID;
      size_t nmiss;
      streamReadRecordF(streamID0, &zlevelsOut[offset], &nmiss);
      if (0 != nmiss) cdo_abort("Output vertical coordinate variables are not allowd to contain missing values.");
    }

  streamClose(streamID0);
}

void *
Intlevel3d(void *process)
{
  cdo_initialize(process);

  // clang-format off
  auto INTLEVEL3D  = cdo_operator_add("intlevel3d",  0, 0, nullptr);
  auto INTLEVELX3D = cdo_operator_add("intlevelx3d", 0, 0, nullptr);
  // clang-format on

  (void) INTLEVEL3D;  // unused

  auto operatorID = cdo_operator_id();
  auto expol = (operatorID == INTLEVELX3D);

  operator_input_arg("tgtcoordinate");

  auto streamID1 = cdo_open_read(0);  // input data

  // Read filename from Parameter
  operator_input_arg("filename for vertical target coordinates variable");
  operator_check_argc(1);

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto taxisID1 = vlistInqTaxis(vlistID1);

  auto vlistID3 = vlistDuplicate(vlistID1);
  auto taxisID3 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID3, taxisID1);

  VarList varList1;
  varListInit(varList1, vlistID1);
  varListSetUniqueMemtype(varList1);
  auto memType = varList1[0].memType;

  // Read 3d source coordinate (streamID2)
  VarList varList2;
  Varray<float> zlevelsIn;

  read_source_coordinate(1, varList2, zlevelsIn);

  auto gridsizei = varList2[0].gridsize;  // horizontal gridsize of input z coordinate
  auto nlevi = varList2[0].nlevels;       // number of input levels for later use

  // Read 3d target coordinate (streamID0)
  VarList varList0;
  Varray<float> zlevelsOut;

  read_target_coordinate(cdo_operator_argv(0), varList0, zlevelsOut);

  auto gridsizeo = varList0[0].gridsize;  // horizontal gridsize of output z coordinate
  auto nlevo = varList0[0].nlevels;       // number of output levels for later use
  auto gridID3 = varList0[0].gridID;
  auto zaxisID3 = varList0[0].zaxisID;

  // gridsize of input and output vertical coordinate must be equal (later use of gridsizeo ONLY)
  if (gridsizei != gridsizeo) cdo_abort("Input and output vertical coordinate must have the same gridsize!");

  auto gridSize = gridsizeo;

  /*
   * Check for the correct vertical axis in the input: Variables with the same
   * number of levels as the input vertical levels from operators parameter (streamID0).
   * Variables with a different z-axis should be copied into output.
   */
  int zaxisID1 = -1;
  auto nzaxis = vlistNzaxis(vlistID1);
  int i;
  for (i = 0; i < nzaxis; ++i)
    {
      auto zaxisID = vlistZaxis(vlistID1, i);
      auto nlevel = zaxisInqSize(zaxisID);
      if (nlevel == nlevi)
        {
          zaxisID1 = zaxisID;
          break;
        }
    }
  if (i == nzaxis) cdo_abort("No processable variable found (vertical coordinate differ)!");

  auto ngrids = vlistNgrids(vlistID1);
  for (i = 0; i < ngrids; ++i)
    {
      auto gridsize = gridInqSize(vlistGrid(vlistID1, i));
      if (gridsize == gridSize) break;
    }
  if (i == nzaxis) cdo_abort("No processable variable found (grid coordinate differ)!");

  // Create weights for later interpolation - assumption: input vertical correct is constant in time
  Varray<int> lev_idx(nlevo * gridSize);
  Varray<float> lev_wgt(nlevo * gridSize);

  vert_gen_weights3d(expol, gridSize, nlevi, zlevelsIn, nlevo, zlevelsOut, lev_idx, lev_wgt);
  varray_free(zlevelsIn);

  for (int index = 0; index < nzaxis; ++index)
    if (zaxisID1 == vlistZaxis(vlistID1, index)) vlistChangeZaxisIndex(vlistID3, index, zaxisID3);

  // add the vertical output field to the output stream
  auto oz3dvarID = vlistDefVar(vlistID3, gridID3, zaxisID3, TIME_VARYING);
  vlist_copy_var_attributes(varList0, 0, vlistID3, oz3dvarID);

  auto streamID3 = cdo_open_write(2);  // output stream
  cdo_def_vlist(streamID3, vlistID3);

  VarList varList3;
  varListInit(varList3, vlistID3);
  varListSetMemtype(varList3, memType);

  auto maxlev = std::max(nlevi, nlevo);
  auto nvars = vlistNvars(vlistID1);

  std::vector<bool> processVars(nvars);
  std::vector<bool> interpVars(nvars, false);         // marker for variables to be interpolated
  std::vector<std::vector<size_t>> varnmiss(nvars);  // can for missing values of arbitrary variables
  Field3DVector vardata1(nvars), vardata2(nvars);

  for (int varID = 0; varID < nvars; ++varID)
    {
      auto zaxisID = varList1[varID].zaxisID;
      auto gridsize = varList1[varID].gridsize;
      auto nlevel = varList1[varID].nlevels;

      vardata1[varID].init(varList1[varID]);

      /*  variabls for interpolation:
       *  * have the required vertical axis, i.e. the correct number of levels (nlevi)
       *  * have the same number of horizontal grid points (i.e. same gridSize) like the two vertical coordinates
       *  * are NOT the output vertical coordinates itself
       */
      interpVars[varID] = (zaxisID == zaxisID1 && varID != oz3dvarID && gridsize == gridSize && gridsizeo == gridsize);
      if (interpVars[varID])
        {
          varnmiss[varID].resize(maxlev, 0);
          vardata2[varID].init(varList3[varID]);
        }
      else
        {
          varnmiss[varID].resize(nlevel);
          if (Options::cdoVerbose)
            cdo_print("Ignore variable %s (levels=%d gridsize=%zu)!", varList1[varID].name, nlevel, gridsize);
        }
    }

  {
    int varID;
    for (varID = 0; varID < nvars; ++varID)
      if (interpVars[varID]) break;
    if (varID == nvars) cdo_abort("No processable variable found!");
  }

  int tsID = 0;
  while (true)
    {
      auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
      if (nrecs == 0) break;

      for (int varID = 0; varID < nvars; ++varID) processVars[varID] = false;

      cdo_taxis_copy_timestep(taxisID3, taxisID1);
      cdo_def_timestep(streamID3, tsID);

      // Read the whole 3d data field
      for (int recID = 0; recID < nrecs; ++recID)
        {
          int varID, levelID;
          cdo_inq_record(streamID1, &varID, &levelID);
          cdo_read_record(streamID1, vardata1[varID], levelID, &varnmiss[varID][levelID]);
          processVars[varID] = true;
        }

      // Perform the interpolation on all valid data variables
      for (int varID = 0; varID < nvars; ++varID)
        {
          if (processVars[varID] && interpVars[varID])
            {
              auto gridsize = varList1[varID].gridsize;
              auto missval = varList1[varID].missval;

              vert_interp_lev3d(gridsize, nlevi, missval, vardata1[varID], vardata2[varID], nlevo, lev_idx, lev_wgt);

              for (int levelID = 0; levelID < nlevo; ++levelID)
                {
                  auto offset = gridsize * levelID;
                  if (memType == MemType::Float)
                    varnmiss[varID][levelID] = array_num_mv(gridsize, &vardata2[varID].vec_f[offset], (float) missval);
                  else
                    varnmiss[varID][levelID] = array_num_mv(gridsize, &vardata2[varID].vec_d[offset], missval);
                }
            }
          else
            {
              if (Options::cdoVerbose && tsID <= 1) cdo_print("Perform no interpolation on variable %s", varList1[varID].name);
            }
        }

      // write the output
      for (int varID = 0; varID < nvars; ++varID)
        {
          if (processVars[varID])
            {
              for (int levelID = 0; levelID < varList3[varID].nlevels; ++levelID)
                {
                  cdo_def_record(streamID3, varID, levelID);
                  cdo_write_record(streamID3, interpVars[varID] ? vardata2[varID] : vardata1[varID], levelID,
                                   varnmiss[varID][levelID]);
                }
            }
        }

      // copy output z coordinate to output stream
      for (int levelID = 0; levelID < varList3[oz3dvarID].nlevels; ++levelID)
        {
          auto offset = varList3[oz3dvarID].gridsize * levelID;
          cdo_def_record(streamID3, oz3dvarID, levelID);
          cdo_write_record_f(streamID3, &zlevelsOut[offset], 0);
        }

      tsID++;
    }

  cdo_stream_close(streamID1);
  cdo_stream_close(streamID3);

  vlistDestroy(vlistID3);

  cdo_finish();

  return nullptr;
}
