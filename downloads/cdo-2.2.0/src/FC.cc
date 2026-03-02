/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      FC         fc2sp           Fourier to spectral
      FC         sp2fc           Spectral to fourier
      FC         fc2gp           Fourier to gridpoint
      FC         gp2fc           Gridpoint to fourier
*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_LIBFFTW3
#include <fftw3.h>
#endif

#include <mutex>

#include <cdi.h>

#include "cdo_vlist.h"
#include "process_int.h"
#include <mpim_grid.h>
#include "specspace.h"
#include "griddes.h"
#include "cimdOmp.h"

#ifdef HAVE_LIBFFTW3
static std::mutex fftwMutex;
#endif

static int
vlistGetFirstReg2DGrid(int vlistID)
{
  // find first gaussian grid
  auto ngrids = vlistNgrids(vlistID);
  for (int index = 0; index < ngrids; ++index)
    {
      auto gridID = vlistGrid(vlistID, index);
      if (gridInqType(gridID) == GRID_GAUSSIAN || gridInqType(gridID) == GRID_LONLAT || gridInqType(gridID) == GRID_CURVILINEAR)
        return gridID;
    }

  return -1;
}

static void
fourier2grid(int gridID1, const Varray<double> &array1, int gridID2, Varray<double> &array2)
{
  (void) gridID2;
  auto nlon = gridInqXsize(gridID1);
  auto nlat = gridInqYsize(gridID1);

#ifdef HAVE_LIBFFTW3
  struct FourierMemory
  {
    fftw_complex *in_fft;
    double *out_fft;
    fftw_plan plan;
  };

  std::vector<FourierMemory> ompmem(Threading::ompNumThreads);

  for (int i = 0; i < Threading::ompNumThreads; ++i)
    {
      ompmem[i].in_fft = fftw_alloc_complex(nlon);
      ompmem[i].out_fft = (double *) fftw_malloc(nlon * sizeof(double));
      std::scoped_lock lock(fftwMutex);
      ompmem[i].plan = fftw_plan_dft_c2r_1d(nlon, ompmem[i].in_fft, ompmem[i].out_fft, FFTW_ESTIMATE);
    }

  if (Options::cdoVerbose) fftw_print_plan(ompmem[0].plan);

#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
  for (size_t ilat = 0; ilat < nlat; ++ilat)
    {
      auto ompthID = cdo_omp_get_thread_num();
      auto in_fft = ompmem[ompthID].in_fft;
      auto out_fft = ompmem[ompthID].out_fft;

      for (size_t ifc = 0; ifc < nlon; ++ifc)
        {
          in_fft[ifc][0] = array1[2 * (ilat * nlon + ifc)];
          in_fft[ifc][1] = array1[2 * (ilat * nlon + ifc) + 1];
        }

      fftw_execute(ompmem[ompthID].plan);

      for (size_t ilon = 0; ilon < nlon; ++ilon) array2[ilat * nlon + ilon] = out_fft[ilon];
    }

  for (int i = 0; i < Threading::ompNumThreads; ++i)
    {
      fftw_free(ompmem[i].in_fft);
      fftw_free(ompmem[i].out_fft);
      fftw_destroy_plan(ompmem[i].plan);
    }
#else
  cdo_abort("FFTW support not compiled in!");
#endif
}

static void
grid2fourier(int gridID1, const Varray<double> &array1, int gridID2, Varray<double> &array2)
{
  (void) gridID2;
  auto nlon = gridInqXsize(gridID1);
  auto nlat = gridInqYsize(gridID1);

#ifdef HAVE_LIBFFTW3
  double norm = 1.0 / nlon;
  struct FourierMemory
  {
    double *in_fft;
    fftw_complex *out_fft;
    fftw_plan plan;
  };

  std::vector<FourierMemory> ompmem(Threading::ompNumThreads);

  for (int i = 0; i < Threading::ompNumThreads; ++i)
    {
      ompmem[i].in_fft = (double *) fftw_malloc(nlon * sizeof(double));
      ompmem[i].out_fft = fftw_alloc_complex(nlon);
      std::scoped_lock lock(fftwMutex);
      ompmem[i].plan = fftw_plan_dft_r2c_1d(nlon, ompmem[i].in_fft, ompmem[i].out_fft, FFTW_ESTIMATE);
    }

  if (Options::cdoVerbose) fftw_print_plan(ompmem[0].plan);

#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
  for (size_t ilat = 0; ilat < nlat; ++ilat)
    {
      auto ompthID = cdo_omp_get_thread_num();
      auto in_fft = ompmem[ompthID].in_fft;
      auto out_fft = ompmem[ompthID].out_fft;

      for (size_t ilon = 0; ilon < nlon; ++ilon) in_fft[ilon] = array1[ilat * nlon + ilon];

      fftw_execute(ompmem[ompthID].plan);

      for (size_t ifc = 0; ifc < nlon; ++ifc)
        {
          array2[2 * (ilat * nlon + ifc)] = norm * out_fft[ifc][0];
          array2[2 * (ilat * nlon + ifc) + 1] = norm * out_fft[ifc][1];
        }
    }

  for (int i = 0; i < Threading::ompNumThreads; ++i)
    {
      fftw_free(ompmem[i].in_fft);
      fftw_free(ompmem[i].out_fft);
      fftw_destroy_plan(ompmem[i].plan);
    }
#else
  cdo_abort("FFTW support not compiled in!");
#endif
}

void *
FC(void *process)
{
  int gridID1 = -1, gridID2 = -1;
  size_t nlon = 0, nlat = 0;
  int ntr = 0;
  int nsp = 0, nfc = 0;

  cdo_initialize(process);

  operator_check_argc(0);

  auto dataIsUnchanged = data_is_unchanged();

  auto FC2SP = cdo_operator_add("fc2sp", 0, 0, nullptr);
  auto SP2FC = cdo_operator_add("sp2fc", 0, 0, nullptr);
  auto FC2GP = cdo_operator_add("fc2gp", 0, 0, nullptr);
  auto GP2FC = cdo_operator_add("gp2fc", 0, 0, nullptr);
  auto GRID2FOURIER = cdo_operator_add("grid2fourier", 1, 0, nullptr);
  auto FOURIER2GRID = cdo_operator_add("fourier2grid", 1, 0, nullptr);

  auto operatorID = cdo_operator_id();
  auto operfunc = cdo_operator_f1(operatorID);

  auto streamID1 = cdo_open_read(0);

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto vlistID2 = vlistDuplicate(vlistID1);

  auto taxisID1 = vlistInqTaxis(vlistID1);
  auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int gridIDsp = -1, gridIDgp = -1, gridIDfc = -1;
  if (operfunc == 0)
    {
      gridIDsp = vlist_get_first_spectral_grid(vlistID1);
      gridIDgp = vlist_get_first_gaussian_grid(vlistID1);
      gridIDfc = vlist_get_first_fourier_grid(vlistID1);
    }

  // define output grid
  if (operatorID == FC2SP)
    {
      if (gridIDfc == -1) cdo_warning("No fourier data found!");

      gridID1 = gridIDfc;

      if (gridID1 != -1)
        {
          nfc = gridInqSize(gridID1);
          ntr = gridInqTrunc(gridID1);
          nlat = nfc_to_nlat(nfc, ntr);

          if (gridIDsp != -1)
            if (ntr != gridInqTrunc(gridIDsp)) gridIDsp = -1;

          if (gridIDsp == -1)
            {
              nsp = (ntr + 1) * (ntr + 2);
              gridIDsp = gridCreate(GRID_SPECTRAL, nsp);
              gridDefTrunc(gridIDsp, ntr);
              gridDefComplexPacking(gridIDsp, 1);
            }

          gridID2 = gridIDsp;
          nlon = 2 * nlat;
          ntr = gridInqTrunc(gridID2);
        }
    }
  else if (operatorID == SP2FC)
    {
      if (gridIDsp == -1) cdo_warning("No spectral data found!");

      gridID1 = gridIDsp;

      if (gridID1 != -1)
        {
          ntr = gridInqTrunc(gridID1);
          nlat = ntr_to_nlat(ntr);

          if (gridIDfc != -1)
            {
              if (ntr != gridInqTrunc(gridIDfc)) gridIDfc = -1;
            }

          if (gridIDfc == -1)
            {
              nfc = 2 * nlat * (ntr + 1);
              gridIDfc = gridCreate(GRID_FOURIER, nfc);
              gridDefTrunc(gridIDfc, ntr);
            }

          gridID2 = gridIDfc;
          nlon = 2 * nlat;
        }
    }
  else if (operatorID == GP2FC)
    {
      if (gridIDgp == -1) cdo_warning("No Gaussian grid data found!");

      gridID1 = gridIDgp;

      if (gridID1 != -1)
        {
          nlon = gridInqXsize(gridID1);
          nlat = gridInqYsize(gridID1);
          ntr = nlat_to_ntr(nlat);

          if (gridIDfc != -1)
            if (ntr != gridInqTrunc(gridIDfc)) gridIDfc = -1;

          if (gridIDfc == -1)
            {
              nfc = 2 * nlat * (ntr + 1);
              gridIDfc = gridCreate(GRID_FOURIER, nfc);
              gridDefTrunc(gridIDfc, ntr);
            }

          gridID2 = gridIDfc;
        }
    }
  else if (operatorID == FC2GP)
    {
      if (gridIDfc == -1) cdo_warning("No fourier data found!");

      gridID1 = gridIDfc;

      if (gridID1 != -1)
        {
          nfc = gridInqSize(gridID1);
          ntr = gridInqTrunc(gridID1);
          nlat = nfc_to_nlat(nfc, ntr);

          if (gridIDgp != -1)
            {
              if (nlat != gridInqYsize(gridIDgp)) gridIDgp = -1;
            }

          if (gridIDgp == -1)
            {
              char gridname[20];
              std::snprintf(gridname, sizeof(gridname), "t%dgrid", ntr);

              gridIDgp = grid_from_name(gridname);
            }

          gridID2 = gridIDgp;
          nlon = gridInqXsize(gridID2);
          nlat = gridInqYsize(gridID2);
        }
    }
  else if (operatorID == FOURIER2GRID)
    {
      gridID1 = vlistGetFirstReg2DGrid(vlistID1);
      if (gridID1 == -1) cdo_warning("No regular 2D data found!");
      if (gridID1 != -1)
        {
          nlon = gridInqXsize(gridID1);
          nlat = gridInqYsize(gridID1);
          gridID2 = gridID1;
        }
    }
  else if (operatorID == GRID2FOURIER)
    {
      gridID1 = vlistGetFirstReg2DGrid(vlistID1);
      if (gridID1 == -1) cdo_warning("No regular 2D data found!");

      if (gridID1 != -1)
        {
          nlon = gridInqXsize(gridID1);
          nlat = gridInqYsize(gridID1);
          gridID2 = gridID1;
        }
    }

  FC_Transformation fcTrans;
  SP_Transformation spTrans;
  if (operfunc == 0)
    {
      if (nlon > 0)
        {
          if (operatorID == GP2FC || operatorID == FC2GP)
            fcTrans.init(nlon, nlat, ntr);
          else if (operatorID == SP2FC)
            spTrans.init(nlon, nlat, ntr, PolFlag::SP2FC);
          else if (operatorID == FC2SP)
            spTrans.init(nlon, nlat, ntr, PolFlag::FC2SP);
        }
    }

  // printf("nfc %d, ntr %d, nlat %zu, nlon %zu\n", nfc, ntr, nlat, nlon);

  auto nvars = vlistNvars(vlistID2);
  std::vector<bool> vars(nvars);
  for (int varID = 0; varID < nvars; ++varID) vars[varID] = gridID1 == vlistInqVarGrid(vlistID1, varID);

  if (gridID1 != -1) vlistChangeGrid(vlistID2, gridID1, gridID2);
  if (operatorID == GRID2FOURIER)
    {
      for (int varID = 0; varID < nvars; ++varID)
        if (vars[varID]) vlistDefVarDatatype(vlistID2, varID, CDI_DATATYPE_CPX64);
    }
  else if (operatorID == FOURIER2GRID)
    {
      for (int varID = 0; varID < nvars; ++varID)
        if (vars[varID]) vlistDefVarDatatype(vlistID2, varID, CDI_DATATYPE_FLT32);
    }

  auto streamID2 = cdo_open_write(1);

  cdo_def_vlist(streamID2, vlistID2);

  auto gridsizemax = vlistGridsizeMax(vlistID1);
  if (operatorID == FOURIER2GRID) gridsizemax *= 2;
  Varray<double> array1(gridsizemax), array2;

  if (gridID2 != -1)
    {
      auto gridsize = gridInqSize(gridID2);
      if (operatorID == GRID2FOURIER) gridsize *= 2;
      array2.resize(gridsize);
    }

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

          if (vars[varID])
            {
              size_t nmiss;
              cdo_read_record(streamID1, array1.data(), &nmiss);
              if (nmiss) cdo_abort("Missing values unsupported for spectral/fourier data!");

              gridID1 = vlistInqVarGrid(vlistID1, varID);
              if (operatorID == FC2SP)
                four2spec(spTrans, gridID1, array1, gridID2, array2);
              else if (operatorID == SP2FC)
                spec2four(spTrans, gridID1, array1, gridID2, array2);
              else if (operatorID == FC2GP)
                four2grid(fcTrans, gridID1, array1, gridID2, array2);
              else if (operatorID == GP2FC)
                grid2four(fcTrans, gridID1, array1, gridID2, array2);
              else if (operatorID == FOURIER2GRID)
                fourier2grid(gridID1, array1, gridID2, array2);
              else if (operatorID == GRID2FOURIER)
                grid2fourier(gridID1, array1, gridID2, array2);

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
