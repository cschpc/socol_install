/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida
          Ralf Quast

*/

#include <climits>
#include <algorithm>
#include <stdlib.h>

#include <cdi.h>

#include "cdo_output.h"
#include "percentiles.h"
#include "percentiles_hist.h"

#define FLT_CAPACITY(n, s) ((int) (((n) * (s)) / sizeof(float)))
#define FLT_PTR(p) ((float *) (p))
#define INT_PTR(p) ((unsigned int *) (p))
#define SHR_PTR(p) ((unsigned short *) (p))

static int
histGetEnvNBins()
{
  constexpr int NBINS_DEFAULT = 101;
  constexpr int NBINS_MINIMUM = 11;

  auto str = getenv("CDO_PCTL_NBINS");

  return (str != nullptr) ? std::max(atoi(str), NBINS_MINIMUM) : NBINS_DEFAULT;
}

template <typename T>
static void
hist_init_bins(int nbins, T *bins)
{
  for (int i = 0; i < nbins; ++i) bins[i] = 0;
}

static void
histDefBounds(Histogram &hist, float a, float b)
{
  assert(hist.nbins > 0);

  hist.nsamp = 0;
  hist.min = std::min(a, b);
  hist.max = std::max(a, b);
  hist.step = (hist.max - hist.min) / hist.nbins;

  if (hist.isUint32)
    hist_init_bins(hist.nbins, INT_PTR(hist.ptr));
  else
    hist_init_bins(hist.nbins, SHR_PTR(hist.ptr));
}

static inline int
calc_bin(int nbins, float histMin, float histStep, float value)
{
  //assert(histStep > 0.0f);
  //return std::min((int) ((value - histMin) / histStep), nbins - 1);
  return (histStep > 0.0f) ? std::min((int) ((value - histMin) / histStep), nbins - 1) : 0;
}

template <typename T>
static void
histBinAddValue(Histogram &hist, T *bins, float value)
{
  auto bin = calc_bin(hist.nbins, hist.min, hist.step, value);
  if (bin >= 0 && bin < hist.nbins) bins[bin]++;
}

template <typename T>
static void
histBinSubValue(Histogram &hist, T *bins, float value)
{
  auto bin = calc_bin(hist.nbins, hist.min, hist.step, value);
  if (bin >= 0 && bin < hist.nbins && bins[bin] > 0) bins[bin]--;
}

static void
histBin(Histogram &hist)
{
  assert(hist.nsamp == hist.capacity);

  std::vector<float> values(hist.nsamp);

  auto fltptr = FLT_PTR(hist.ptr);
  for (int i = 0; i < hist.nsamp; ++i) values[i] = fltptr[i];

  if (hist.isUint32)
    hist_init_bins(hist.nbins, INT_PTR(hist.ptr));
  else
    hist_init_bins(hist.nbins, SHR_PTR(hist.ptr));

  if (hist.isUint32)
    for (int i = 0; i < hist.nsamp; ++i) histBinAddValue(hist, INT_PTR(hist.ptr), values[i]);
  else
    for (int i = 0; i < hist.nsamp; ++i) histBinAddValue(hist, SHR_PTR(hist.ptr), values[i]);
}
/* unused
static int
histReset(Histogram &hist)
{
  assert(hist.nbins > 0);

  if (hist.nsamp < hist.capacity)
    {
      for (int i = 0; i < hist.nsamp; ++i) FLT_PTR(hist.ptr)[i] = 0.;
    }
  else
    {
      if (hist.isUint32)
        for (int i = 0; i < hist.nbins; ++i) INT_PTR(hist.ptr)[i] = 0;
      else
        for (int i = 0; i < hist.nbins; ++i) SHR_PTR(hist.ptr)[i] = 0;
    }

  hist.nsamp = 0;

  return 0;
}
*/

static void
histCheckValue(const Histogram &hist, float &value)
{
  // 2011-08-01 Uwe Schulzweida: added check for rounding errors
  if (value < hist.min && (hist.min - value) < 1.e5f) value = hist.min;
  if (value > hist.max && (value - hist.max) < 1.e5f) value = hist.max;
}

static int
histAddValue(Histogram &hist, float value)
{
  assert(hist.nbins > 0);

  // if (is_equal(hist.min, hist.max)) return 0;

  histCheckValue(hist, value);
  if (value < hist.min || value > hist.max) return 1;

  if (hist.nsamp < hist.capacity) { FLT_PTR(hist.ptr)[hist.nsamp] = value; }
  else
    {
      if (hist.nsamp == hist.capacity) histBin(hist);

      if (hist.isUint32)
        histBinAddValue(hist, INT_PTR(hist.ptr), value);
      else
        histBinAddValue(hist, SHR_PTR(hist.ptr), value);
    }

  hist.nsamp++;

  return 0;
}

static void
histRemoveValue(Histogram &hist, float value)
{
  auto fltptr = FLT_PTR(hist.ptr);

  int i = 0;
  for (i = 0; i < hist.nsamp; ++i)
    {
      if (is_equal(fltptr[i], value))
        {
          if (i != hist.nsamp - 1) fltptr[i] = fltptr[hist.nsamp - 1];
          break;
        }
    }
  if (i == hist.nsamp)
    cdo_warning("'%f' not found in histogram!", value);
  else
    hist.nsamp--;
}

static int
histSubValue(Histogram &hist, float value)
{
  assert(hist.nbins > 0);

  if (is_equal(hist.min, hist.max)) return 0;

  histCheckValue(hist, value);
  if (value < hist.min || value > hist.max) return 1;

  if (hist.nsamp < hist.capacity) { histRemoveValue(hist, value); }
  else if (hist.nsamp > hist.capacity)
    {
      if (hist.isUint32)
        histBinSubValue(hist, INT_PTR(hist.ptr), value);
      else
        histBinSubValue(hist, SHR_PTR(hist.ptr), value);

      hist.nsamp--;
    }
  else
    return 1;

  return 0;
}

template <typename T>
static double
histGetBin(int nbins, double s, const T *ptr)
{
  int i = 0, count = 0;

  do count += ptr[i++];
  while (count < s);

  assert(i > 0);
  assert(i - 1 < nbins);
  assert(ptr[i - 1] > 0);

  double t = (count - s) / ptr[i - 1];

  return (i - t);
}

static double
histGetPercentile(const Histogram &hist, double p)
{
  assert(hist.nsamp > 0);
  assert(hist.nbins > 0);
  assert(p >= 0.0);
  assert(p <= 100.0);

  if (hist.nsamp > hist.capacity)
    {
      static auto lprint = true;
      if (lprint && Options::cdoVerbose)
        {
          lprint = false;
          cdo_print("Using percentile method: histogram with %d bins", hist.nbins);
        }

      double s = hist.nsamp * (p / 100.0);

      auto bin = hist.isUint32 ? histGetBin(hist.nbins, s, INT_PTR(hist.ptr)) : histGetBin(hist.nbins, s, SHR_PTR(hist.ptr));

      //assert(hist.step > 0.0f);

      return hist.min + bin * hist.step;
    }
  else { return percentile(FLT_PTR(hist.ptr), hist.nsamp, p); }
}

void
HistogramSet::createVarLevels(int varID, int nlevels, size_t nhists)
{
  auto nbins = histGetEnvNBins();

  assert(nbins > 0);
  assert(nlevels > 0);
  assert(nhists > 0);

  if (varID < 0 || varID >= nvars) cdo_abort("Illegal argument: varID %d is undefined (%s)", varID, __func__);

  this->var_nlevels[varID] = nlevels;
  this->var_nhists[varID] = nhists;
  this->histograms[varID].resize(nlevels);

  for (int levelID = 0; levelID < nlevels; ++levelID)
    {
      this->histograms[varID][levelID].resize(nhists);
      auto &hists = this->histograms[varID][levelID];

      for (size_t histID = 0; histID < nhists; histID++)
        {
          hists[histID].min = 0.0f;
          hists[histID].max = 0.0f;
          hists[histID].step = 0.0f;
          hists[histID].nbins = nbins;
          hists[histID].nsamp = 0;
          auto isUint32 = (this->nsteps >= USHRT_MAX);
          auto vsize = isUint32 ? sizeof(unsigned int) : sizeof(unsigned short);
          hists[histID].isUint32 = isUint32;
          hists[histID].capacity = FLT_CAPACITY(nbins, vsize);
          hists[histID].ptr = malloc(nbins * vsize);
          if (hists[histID].ptr == nullptr) cdo_abort("Not enough memory (%s)", __func__);
        }
    }
}

template <typename T1, typename T2>
static void
def_bounds(size_t nhists, std::vector<Histogram> &hists, const Varray<T1> &v1, const Varray<T2> &v2, float mv1, float mv2)
{
  assert(!v1.empty());
  assert(!v2.empty());

  for (size_t i = 0; i < nhists; ++i)
    {
      float a = v1[i];
      float b = v2[i];

      if (dbl_is_equal(a, mv1) || dbl_is_equal(b, mv2) /*|| dbl_is_equal(a, b)*/)
        histDefBounds(hists[i], 0.0f, 0.0f);
      else
        histDefBounds(hists[i], a, b);
    }
}

void
HistogramSet::defVarLevelBounds(int varID, int levelID, const Field &field1, const Field &field2)
{
  if (varID < 0 || varID >= nvars) cdo_abort("Illegal argument: varID %d is undefined (%s)", varID, __func__);

  auto nlevels = this->var_nlevels[varID];
  if (levelID < 0 || levelID >= nlevels) cdo_abort("Illegal argument: levelID %d is undefined (%s)", levelID, __func__);

  auto nhists = this->var_nhists[varID];
  if (nhists != gridInqSize(field1.grid) || nhists != gridInqSize(field2.grid)) cdo_abort("Grids are different (%s)", __func__);

  auto &hists = this->histograms[varID][levelID];

  float mv1 = (float) field1.missval;
  float mv2 = (float) field2.missval;

  // clang-format off
  if      (field1.memType == MemType::Float && field2.memType == MemType::Float)
    def_bounds(nhists, hists, field1.vec_f, field2.vec_f, mv1, mv2);
  else if (field1.memType == MemType::Float && field2.memType == MemType::Double)
    def_bounds(nhists, hists, field1.vec_f, field2.vec_d, mv1, mv2);
  else if (field1.memType == MemType::Double && field2.memType == MemType::Float)
    def_bounds(nhists, hists, field1.vec_d, field2.vec_f, mv1, mv2);
  else if (field1.memType == MemType::Double && field2.memType == MemType::Double)
    def_bounds(nhists, hists, field1.vec_d, field2.vec_d, mv1, mv2);
  // clang-format on
}

template <typename T>
static int
histAddVarLevelValues(size_t nhists, std::vector<Histogram> &hists, const Varray<T> &v, size_t nmiss, T mv)
{
  assert(!v.empty());

  int nign = 0;
  if (nmiss)
    {
      for (size_t i = 0; i < nhists; ++i)
        if (!dbl_is_equal(v[i], mv)) nign += histAddValue(hists[i], v[i]);
    }
  else
    {
      for (size_t i = 0; i < nhists; ++i) nign += histAddValue(hists[i], v[i]);
    }

  return nign;
}

template <typename T>
static int
histSubVarLevelValues(size_t nhists, std::vector<Histogram> &hists, const Varray<T> &v, size_t nmiss, T mv)
{
  assert(!v.empty());

  int nign = 0;
  if (nmiss)
    {
      for (size_t i = 0; i < nhists; ++i)
        if (!dbl_is_equal(v[i], mv)) nign += histSubValue(hists[i], v[i]);
    }
  else
    {
      for (size_t i = 0; i < nhists; ++i) nign += histSubValue(hists[i], v[i]);
    }

  return nign;
}

int
HistogramSet::addVarLevelValues(int varID, int levelID, const Field &field)
{
  if (varID < 0 || varID >= nvars) cdo_abort("Illegal argument: varID %d is undefined (%s)", varID, __func__);

  auto nlevels = this->var_nlevels[varID];
  if (levelID < 0 || levelID >= nlevels) cdo_abort("Illegal argument: levelID %d is undefined (%s)", levelID, __func__);

  auto nhists = this->var_nhists[varID];
  if (nhists != gridInqSize(field.grid)) cdo_abort("Grids are different (%s)", __func__);

  auto &hists = this->histograms[varID][levelID];

  int nign = 0;

  if (field.memType == MemType::Float)
    nign = histAddVarLevelValues(nhists, hists, field.vec_f, field.nmiss, (float) field.missval);
  else
    nign = histAddVarLevelValues(nhists, hists, field.vec_d, field.nmiss, field.missval);

  if (nign)
    {
      cdo_warning("%d out of %d grid values are out of bounds and have been ignored (%s)", nign, nhists, __func__);
      return 1;
    }

  return 0;
}

int
HistogramSet::subVarLevelValues(int varID, int levelID, const Field &field)
{
  if (varID < 0 || varID >= nvars) cdo_abort("Illegal argument: varID %d is undefined (%s)", varID, __func__);

  auto nlevels = this->var_nlevels[varID];
  if (levelID < 0 || levelID >= nlevels) cdo_abort("Illegal argument: levelID %d is undefined (%s)", levelID, __func__);

  auto nhists = this->var_nhists[varID];
  if (nhists != gridInqSize(field.grid)) cdo_abort("Grids are different (%s)", __func__);

  auto &hists = this->histograms[varID][levelID];

  int nign = 0;

  if (field.memType == MemType::Float)
    nign = histSubVarLevelValues(nhists, hists, field.vec_f, field.nmiss, (float) field.missval);
  else
    nign = histSubVarLevelValues(nhists, hists, field.vec_d, field.nmiss, field.missval);

  if (nign)
    {
      cdo_warning("%d out of %d grid values are out of bounds and have been ignored (%s)", nign, nhists, __func__);
      return 1;
    }

  return 0;
}
/* unused
void
HistogramSet::reset(int varID, int levelID)
{
  assert(nvars > 0);

  if (varID < 0 || varID >= nvars) cdo_abort("Illegal argument: varID %d is undefined (%s)", varID, __func__);

  auto nlevels = this->var_nlevels[varID];
  assert(nlevels > 0);

  if (levelID < 0 || levelID >= nlevels) cdo_abort("Illegal argument: levelID %d is undefined (%s)", levelID, __func__);

  auto nhists = this->var_nhists[varID];
  assert(nhists > 0);

  auto &hists = this->histograms[varID][levelID];
  for (size_t i = 0; i < nhists; ++i) histReset(hists[i]);
}
*/

template <typename T>
static size_t
calcPercentile(size_t nhists, const std::vector<Histogram> &hists, double p, Varray<T> &v, T mv)
{
  assert(!v.empty());

  size_t nmiss = 0;

  for (size_t i = 0; i < nhists; ++i)
    {
      if (hists[i].nsamp)
        {
          v[i] = histGetPercentile(hists[i], p);
        }
      else
        {
          v[i] = mv;
          nmiss++;
        }
    }

  return nmiss;
}

void
HistogramSet::getVarLevelPercentiles(Field &field, int varID, int levelID, double p)
{
  if (varID < 0 || varID >= nvars) cdo_abort("Illegal argument: varID %d is undefined (%s)", varID, __func__);

  auto nlevels = this->var_nlevels[varID];
  if (levelID < 0 || levelID >= nlevels) cdo_abort("Illegal argument: levelID %d is undefined (%s)", levelID, __func__);

  auto nhists = this->var_nhists[varID];
  if (nhists != gridInqSize(field.grid)) cdo_abort("Grids are different (%s)", __func__);

  const auto &hists = this->histograms[varID][levelID];

  if (field.memType == MemType::Float)
    field.nmiss = calcPercentile(nhists, hists, p, field.vec_f, (float) field.missval);
  else
    field.nmiss = calcPercentile(nhists, hists, p, field.vec_d, field.missval);
}
