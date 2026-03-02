/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_LIBFFTW3
#include <fftw3.h>
#endif

#include <cdi.h>

#include "process_int.h"
#include "cdo_vlist.h"
#include "param_conversion.h"
#include "cdo_fft.h"
#include "cdo_options.h"
#include "cimdOmp.h"
#include "field_functions.h"

struct FourierMemory
{
  Varray<double> real;
  Varray<double> imag;
  Varray<double> work_r;
  Varray<double> work_i;
#ifdef HAVE_LIBFFTW3
  fftw_complex *in_fft;
  fftw_complex *out_fft;
  fftw_plan plan;
#endif
};

static void
fourier_fftw(int sign, int varID, int levelID, int nts, size_t gridsize, double missval, FieldVector3D &vars,
             std::vector<FourierMemory> &fourierMemory)
{
  (void) sign;
#ifdef HAVE_LIBFFTW3
  auto norm = 1.0 / std::sqrt(nts);

#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
  for (size_t i = 0; i < gridsize; ++i)
    {
      auto ompthID = cdo_omp_get_thread_num();
      auto &fm = fourierMemory[ompthID];

      auto hasMissvals = false;
      for (int tsID = 0; tsID < nts; ++tsID)
        {
          auto real = vars[tsID][varID][levelID].vec_d[2 * i];
          auto imag = vars[tsID][varID][levelID].vec_d[2 * i + 1];
          fm.in_fft[tsID][0] = real;
          fm.in_fft[tsID][1] = imag;
          if (DBL_IS_EQUAL(real, missval) || DBL_IS_EQUAL(imag, missval)) hasMissvals = true;
        }

      if (hasMissvals)
        {
          for (int tsID = 0; tsID < nts; ++tsID)
            {
              vars[tsID][varID][levelID].vec_d[2 * i] = missval;
              vars[tsID][varID][levelID].vec_d[2 * i + 1] = missval;
            }
        }
      else
        {
          fftw_execute(fm.plan);

          for (int tsID = 0; tsID < nts; ++tsID)
            {
              vars[tsID][varID][levelID].vec_d[2 * i] = fm.out_fft[tsID][0] * norm;
              vars[tsID][varID][levelID].vec_d[2 * i + 1] = fm.out_fft[tsID][1] * norm;
            }
        }
    }
#endif
}

static void
fourier_intrinsic(int sign, int varID, int levelID, int nts, size_t gridsize, double missval, FieldVector3D &vars,
                  std::vector<FourierMemory> &fourierMemory)
{
  auto isPower2 = ((nts & (nts - 1)) == 0);

#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
  for (size_t i = 0; i < gridsize; ++i)
    {
      auto ompthID = cdo_omp_get_thread_num();
      auto &fm = fourierMemory[ompthID];

      auto hasMissvals = false;
      for (int tsID = 0; tsID < nts; ++tsID)
        {
          auto real = vars[tsID][varID][levelID].vec_d[2 * i];
          auto imag = vars[tsID][varID][levelID].vec_d[2 * i + 1];
          fm.real[tsID] = real;
          fm.imag[tsID] = imag;
          if (DBL_IS_EQUAL(real, missval) || DBL_IS_EQUAL(imag, missval)) hasMissvals = true;
        }

      if (hasMissvals)
        {
          for (int tsID = 0; tsID < nts; ++tsID)
            {
              vars[tsID][varID][levelID].vec_d[2 * i] = missval;
              vars[tsID][varID][levelID].vec_d[2 * i + 1] = missval;
            }
        }
      else
        {
          if (isPower2)  // nts is a power of 2
            cdo::fft(fm.real.data(), fm.imag.data(), nts, sign);
          else
            cdo::ft_r(fm.real.data(), fm.imag.data(), nts, sign, fm.work_r.data(), fm.work_i.data());

          for (int tsID = 0; tsID < nts; ++tsID)
            {
              vars[tsID][varID][levelID].vec_d[2 * i] = fm.real[tsID];
              vars[tsID][varID][levelID].vec_d[2 * i + 1] = fm.imag[tsID];
            }
        }
    }
}

void *
Fourier(void *process)
{
  size_t nalloc = 0;

  cdo_initialize(process);

  bool use_fftw = false;
  if (Options::Use_FFTW)
    {
#ifdef HAVE_LIBFFTW3
      if (Options::cdoVerbose) cdo_print("Using fftw3 lib");
      use_fftw = true;
#else
      if (Options::cdoVerbose) cdo_print("LIBFFTW3 support not compiled in!");
#endif
    }

  if (Options::cdoVerbose && !use_fftw) cdo_print("Using intrinsic FFT function!");

  operator_input_arg("the sign of the exponent (-1 for normal or 1 for reverse transformation)!");
  auto sign = parameter_to_int(cdo_operator_argv(0));

  auto streamID1 = cdo_open_read(0);

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto vlistID2 = vlistDuplicate(vlistID1);

  auto taxisID1 = vlistInqTaxis(vlistID1);
  auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  auto streamID2 = cdo_open_write(1);
  cdo_def_vlist(streamID2, vlistID2);

  VarList varList;
  varListInit(varList, vlistID1);

  auto nvars = vlistNvars(vlistID1);
  FieldVector3D vars;
  std::vector<CdiDateTime> vDateTimes;

  int tsID = 0;
  while (true)
    {
      auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
      if (nrecs == 0) break;

      if ((size_t) tsID >= nalloc)
        {
          constexpr size_t NALLOC_INC = 1024;
          nalloc += NALLOC_INC;
          vDateTimes.resize(nalloc);
          vars.resize(nalloc);
        }

      vDateTimes[tsID] = taxisInqVdatetime(taxisID1);

      fields_from_vlist(vlistID1, vars[tsID]);

      for (int recID = 0; recID < nrecs; ++recID)
        {
          int varID, levelID;
          cdo_inq_record(streamID1, &varID, &levelID);
          auto gridsize = varList[varID].gridsize;
          vars[tsID][varID][levelID].resize(2 * gridsize);
          size_t nmiss;
          cdo_read_record(streamID1, vars[tsID][varID][levelID].vec_d.data(), &nmiss);
          vars[tsID][varID][levelID].nmiss = nmiss;
        }

      tsID++;
    }

  int nts = tsID;

  std::vector<FourierMemory> fourierMemory(Threading::ompNumThreads);

  if (use_fftw)
    {
#ifdef HAVE_LIBFFTW3
      for (auto &fm : fourierMemory)
        {
          fm.in_fft = fftw_alloc_complex(nts);
          fm.out_fft = fftw_alloc_complex(nts);
          fm.plan = fftw_plan_dft_1d(nts, fm.in_fft, fm.out_fft, sign, FFTW_ESTIMATE);
        }
      if (Options::cdoVerbose) fftw_print_plan(fourierMemory[0].plan);
#endif
    }
  else
    {
      const bool isPower2 = ((nts & (nts - 1)) == 0);
      for (auto &fm : fourierMemory)
        {
          fm.real.resize(nts);
          fm.imag.resize(nts);
          if (!isPower2) fm.work_r.resize(nts);
          if (!isPower2) fm.work_i.resize(nts);
        }
    }

  for (int varID = 0; varID < nvars; ++varID)
    {
      auto missval = varList[varID].missval;
      auto gridsize = varList[varID].gridsize;
      auto nlevels = varList[varID].nlevels;
      for (int levelID = 0; levelID < nlevels; ++levelID)
        {
          if (use_fftw)
            fourier_fftw(sign, varID, levelID, nts, gridsize, missval, vars, fourierMemory);
          else
            fourier_intrinsic(sign, varID, levelID, nts, gridsize, missval, vars, fourierMemory);
        }
    }

#ifdef HAVE_LIBFFTW3
  if (use_fftw)
    {
      for (auto &fm : fourierMemory)
        {
          fftw_free(fm.in_fft);
          fftw_free(fm.out_fft);
          fftw_destroy_plan(fm.plan);
        }
      fftw_cleanup();
    }
#endif

  for (tsID = 0; tsID < nts; ++tsID)
    {
      taxisDefVdatetime(taxisID2, vDateTimes[tsID]);
      cdo_def_timestep(streamID2, tsID);

      for (int varID = 0; varID < nvars; ++varID)
        {
          auto nlevels = varList[varID].nlevels;
          for (int levelID = 0; levelID < nlevels; ++levelID)
            {
              if (!vars[tsID][varID][levelID].empty())
                {
                  auto nmiss = vars[tsID][varID][levelID].nmiss;
                  cdo_def_record(streamID2, varID, levelID);
                  cdo_write_record(streamID2, vars[tsID][varID][levelID].vec_d.data(), nmiss);
                }
            }
        }
    }

  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  cdo_finish();

  return nullptr;
}
