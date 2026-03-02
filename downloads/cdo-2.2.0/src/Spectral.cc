/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Spectral   sp2gp           Spectral to gridpoint
      Spectral   sp2gpl          Spectral to gridpoint linear (sp2gp,linear)
      Spectral   gp2sp           Gridpoint to spectral
      Spectral   gp2spl          Gridpoint to spectral linear (gp2sp,linear)
      Spectral   sp2sp           Spectral to spectral
      Spectral   spcut           Cut spectral wave number
*/

#include <cdi.h>

#include "cdo_vlist.h"
#include "process_int.h"
#include "param_conversion.h"
#include <mpim_grid.h>
#include "griddes.h"
#include "specspace.h"

static int
gp2sp_init(SP_Transformation &spTrans, int gridID1, int gridIDsp, int defaultTrunc, int (*nlat2ntr)(int))
{
  long nlon = gridInqXsize(gridID1);
  long nlat = gridInqYsize(gridID1);

  long ntr = nlat2ntr(nlat);
  if (defaultTrunc > 0)
    {
      if (defaultTrunc > ntr) cdo_abort("Output trunctation=%d muss be lower than input trunctation=%d", defaultTrunc, ntr);
      ntr = defaultTrunc;
    }
  if (Options::cdoVerbose) cdo_print("trunc=%ld\n", ntr);

  if (gridIDsp != -1)
    if (ntr != gridInqTrunc(gridIDsp)) gridIDsp = -1;

  if (gridIDsp == -1)
    {
      gridIDsp = gridCreate(GRID_SPECTRAL, (ntr + 1) * (ntr + 2));
      gridDefTrunc(gridIDsp, ntr);
      gridDefComplexPacking(gridIDsp, 1);
    }

  if (gridIDsp == -1) cdo_abort("Computation of spherical harmonics failed!");

  int gridID2 = gridIDsp;

  ntr = gridInqTrunc(gridID2);
  spTrans.init(nlon, nlat, ntr, PolFlag::FC2SP);

  return gridID2;
}

static int
sp2gp_init(SP_Transformation &spTrans, int gridID1, int gridIDsp, int gridIDgp, int defaultTrunc, int (*nlat2ntr)(int), const char *ctype)
{
  if (defaultTrunc > 0) gridIDgp = -1;

  if (gridIDgp != -1)
    {
      long nlat = gridInqYsize(gridIDgp);
      long ntr = nlat2ntr(nlat);
      if (gridInqTrunc(gridIDsp) != ntr) gridIDgp = -1;
    }

  if (gridIDgp == -1)
    {
      int ntr = gridInqTrunc(gridIDsp);
      if (defaultTrunc > 0)
        {
          if (defaultTrunc < ntr) cdo_abort("Output trunctation=%d muss be greater than input trunctation=%d", defaultTrunc, ntr);
          ntr = defaultTrunc;
        }
      char gridname[20];
      std::snprintf(gridname, sizeof(gridname), "t%s%dgrid", ctype, ntr);
      gridIDgp = grid_from_name(gridname);
    }

  int gridID2 = gridIDgp;

  long ntr = gridInqTrunc(gridID1);
  long nlon = gridInqXsize(gridID2);
  long nlat = gridInqYsize(gridID2);
  spTrans.init(nlon, nlat, ntr, PolFlag::SP2FC);

  return gridID2;
}

void *
Spectral(void *process)
{
  int gridID1 = -1, gridID2 = -1;
  int defaultTrunc = 0;
  Varray<int> waves;
  SP_Transformation spTrans;

  cdo_initialize(process);

  auto dataIsUnchanged = data_is_unchanged();

  // clang-format off
  auto GP2SP  = cdo_operator_add("gp2sp",  0, 0, nullptr);
  auto GP2SPL = cdo_operator_add("gp2spl", 0, 0, nullptr);
  auto SP2GP  = cdo_operator_add("sp2gp",  0, 0, nullptr);
  auto SP2GPL = cdo_operator_add("sp2gpl", 0, 0, nullptr);
  auto SP2SP  = cdo_operator_add("sp2sp",  0, 0, nullptr);
  auto SPCUT  = cdo_operator_add("spcut",  0, 0, nullptr);

  auto operatorID = cdo_operator_id();

  auto lgp2sp = (operatorID == GP2SP || operatorID == GP2SPL);
  auto lsp2gp = (operatorID == SP2GP || operatorID == SP2GPL);
  auto linear = (operatorID == GP2SPL || operatorID == SP2GPL);

  int (*nlat2ntr)(int) = linear ? nlat_to_ntr_linear : nlat_to_ntr;
  const char *ctype = linear ? "l" : "";

  auto paramArgc = cdo_operator_argc();
  if ((lgp2sp || lsp2gp) && paramArgc == 1)
    {
      const auto &parg = cdo_operator_argv(0);
      auto pos = parg.find("=");
      if (pos > 0 && parg.substr(0, pos) == "trunc")
        {
          defaultTrunc = parameter_to_int(parg.substr(pos + 1));
        }
      else
        {
          auto type = parameter_to_word((pos > 0 && parg.substr(0, pos) == "type") ? parg.substr(pos + 1) : parg);
          if      (type == "linear")    { nlat2ntr = nlat_to_ntr_linear; ctype = "l"; }
          else if (type == "cubic")     { nlat2ntr = nlat_to_ntr_cubic; ctype = "c"; }
          else if (type == "quadratic") { nlat2ntr = nlat_to_ntr; }
          else cdo_abort("Unsupported type: %s\n", type);
        }
    }
  else if (paramArgc > 0 && operatorID != SP2SP && operatorID != SPCUT)
    {
      cdo_abort("Too many parameters");
    }
  // clang-format on

  auto streamID1 = cdo_open_read(0);

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto vlistID2 = vlistDuplicate(vlistID1);

  auto taxisID1 = vlistInqTaxis(vlistID1);
  auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  auto gridIDsp = vlist_get_first_spectral_grid(vlistID1);
  auto gridIDgp = vlist_get_first_gaussian_grid(vlistID1);

  // define output grid
  if (lgp2sp)
    {
      if (gridIDgp == -1) cdo_warning("No data on regular Gaussian grid found!");

      gridID1 = gridIDgp;
      if (gridID1 != -1) gridID2 = gp2sp_init(spTrans, gridID1, gridIDsp, defaultTrunc, nlat2ntr);
    }
  else if (lsp2gp)
    {
      if (gridIDsp == -1) cdo_warning("No spectral data found!");

      gridID1 = gridIDsp;
      if (gridID1 != -1) gridID2 = sp2gp_init(spTrans, gridID1, gridIDsp, gridIDgp, defaultTrunc, nlat2ntr, ctype);
    }
  else if (operatorID == SP2SP)
    {
      gridID1 = gridIDsp;

      operator_input_arg("truncation");
      if (gridID1 != -1)
        {
          if (!isdigit(cdo_operator_argv(0)[0])) cdo_abort("parameter truncation must comprise only digits [0-9]!");
          long ntr = parameter_to_int(cdo_operator_argv(0));
          long nsp = (ntr + 1) * (ntr + 2);
          gridIDsp = gridCreate(GRID_SPECTRAL, nsp);
          gridDefTrunc(gridIDsp, ntr);
          gridDefComplexPacking(gridIDsp, 1);
        }
      else
        cdo_abort("No spectral data found!");

      gridID2 = gridIDsp;
    }
  else if (operatorID == SPCUT)
    {
      gridID1 = gridIDsp;

      operator_input_arg("wave numbers");
      if (gridID1 != -1)
        {
          long maxntr = 1 + gridInqTrunc(gridID1);
          auto wnums = cdo_argv_to_int(cdo_get_oper_argv());
          long ncut = wnums.size();
          waves.resize(maxntr);
          for (long i = 0; i < maxntr; ++i) waves[i] = 1;
          for (long i = 0; i < ncut; ++i)
            {
              long j = wnums[i] - 1;
              if (j < 0 || j >= maxntr) cdo_abort("wave number %ld out of range (min=1, max=%l qd)!", wnums[i], maxntr);
              waves[j] = 0;
            }
        }
      else
        cdo_abort("No spectral data found!");

      gridID2 = gridIDsp;
    }

  auto nvars = vlistNvars(vlistID2);
  std::vector<bool> processVars(nvars);
  for (int varID = 0; varID < nvars; ++varID) processVars[varID] = (gridID1 == vlistInqVarGrid(vlistID1, varID));

  if (gridID1 != -1) vlistChangeGrid(vlistID2, gridID1, gridID2);

  auto streamID2 = cdo_open_write(1);
  cdo_def_vlist(streamID2, vlistID2);

  auto gridsizemax = vlistGridsizeMax(vlistID1);
  Varray<double> array1(gridsizemax), array2;
  if (gridID2 != -1) array2.resize(gridInqSize(gridID2));

  int tsID = 0;
  while (true)
    {
      auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
      if (nrecs == 0) break;

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_def_timestep(streamID2, tsID);

      for (int recID = 0; recID < nrecs; ++recID)
        {
          int varID, levelID;
          cdo_inq_record(streamID1, &varID, &levelID);

          if (processVars[varID])
            {
              size_t nmiss;
              cdo_read_record(streamID1, array1.data(), &nmiss);
              if (nmiss) cdo_abort("Missing values unsupported for spectral data!");

              gridID1 = vlistInqVarGrid(vlistID1, varID);
              // clang-format off
              if      (lgp2sp) grid2spec(spTrans, gridID1, array1, gridID2, array2);
              else if (lsp2gp) spec2grid(spTrans, gridID1, array1, gridID2, array2);
              else if (operatorID == SP2SP) spec2spec(gridID1, array1, gridID2, array2);
              else if (operatorID == SPCUT) speccut(gridID1, array1, array2, waves);
              // clang-format on

              cdo_def_record(streamID2, varID, levelID);
              cdo_write_record(streamID2, array2.data(), nmiss);
            }
          else
            {
              cdo_def_record(streamID2, varID, levelID);
              if (dataIsUnchanged) { cdo_copy_record(streamID2, streamID1); }
              else
                {
                  size_t nmiss;
                  cdo_read_record(streamID1, array1.data(), &nmiss);
                  cdo_write_record(streamID2, array1.data(), nmiss);
                }
            }
        }

      tsID++;
    }

  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  cdo_finish();

  return nullptr;
}
