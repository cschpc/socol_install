/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cmath>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "cdo_options.h"
#include "cimdOmp.h"
#include "constants.h"
#include "varray.h"
#include "gaussian_latitudes.h"

static void
jspleg1(double *pleg, double plat, long ktrunc, double *work)
{
  /*
     jspleg1 - Routine to calculate legendre functions

     Purpose
     --------

     This routine calculates the legendre functions for one latitude.
     (but not their derivatives)

     Interface
     ----------

     jspleg1( pleg, plat, ktrunc)

     Input parameters
     ----------------

     plat      - Latitude in radians
     ktrunc    - Spectral truncation

     Output parameters
     -----------------

     pleg      - Array of legendre functions for one latitude.
                 The array must be at least (KTRUNC+1)*(KTRUNC+4)/2 words long.

     Method
     ------

     Recurrence relation with explicit relations for P(m,m) and P(m,m+1)

     AUTHOR
     ------

     J.D.Chambers         ECMWF        9 November 1993

     Modifications
     -------------

     None
  */

  // Initialization

  const auto itout1 = ktrunc + 1;
  // const auto zsin   = std::sin(plat);
  const auto zsin = plat;
  const auto zcos = std::sqrt(1.0 - zsin * zsin);

  double *zhlp1 = work;
  double *zhlp2 = work + itout1;
  double *zhlp3 = work + itout1 + itout1;

  //  Step 1.        M = 0, N = 0 and N = 1

  long ilm = 1;
  pleg[0] = 1.0;
  auto zf1m = std::sqrt(3.0);
  pleg[1] = zf1m * zsin;

  //  Step 2.       Sum for M = 0 to T (T = truncation)

  for (long jm = 1; jm < itout1; jm++)
    {
      zhlp1[jm] = std::sqrt(2.0 * jm + 3.0);
      zhlp2[jm] = 1.0 / std::sqrt(2.0 * jm);
    }

  zhlp1[0] = std::sqrt(3.0);

  for (long jm = 0; jm < itout1; jm++)
    {
      const auto i1m = jm - 1;
      const auto zre1 = zhlp1[jm];
      auto ze1 = 1.0 / zre1;

      //   Step 3.       M > 0 only

      if (i1m != -1)
        {
          const auto zf2m = zf1m * zcos * zhlp2[jm];
          zf1m = zf2m * zre1;

          //  Step 4.       N = M and N = M+1

          pleg[++ilm] = zf2m;
          pleg[++ilm] = zf1m * zsin;

          // When output truncation is reached, return to calling program

          if (jm == (itout1 - 1)) break;
        }

      //  Step 5.       Sum for N = M+2 to T+1

      const double zjmsqr = jm * jm;
      const auto im2 = i1m + 2;

      for (long jcn = im2; jcn < itout1; jcn++)
        {
          const double znsqr = (jcn + 1) * (jcn + 1);
          zhlp3[jcn] = std::sqrt((4.0 * znsqr - 1.0) / (znsqr - zjmsqr));
        }

      for (long jcn = im2; jcn < itout1; jcn++)
        {
          const auto ze2 = zhlp3[jcn];
          ilm++;
          pleg[ilm] = ze2 * (zsin * pleg[ilm - 1] - ze1 * pleg[ilm - 2]);
          ze1 = 1.0 / ze2;
        }
    }
}

// phcs - Compute values of Legendre polynomials and their meridional derivatives
static void
phcs(bool needHnm, double *pnm, double *hnm, long waves, double pmu, double *ztemp1, double *ztemp2)
{
  long jnmjk;
  double zcospar, zsinpar, zcosfak, zsinfak;
  double zq, zwm2, zw, zwq, zq2m1, zwm2q2, z2q2, zcnm, zdnm, zenm;

  const auto twowaves = waves << 1;

  const auto zcos2 = std::sqrt(1.0 - pmu * pmu);
  const auto lat = std::acos(pmu);
  auto zan = 1.0;

  ztemp1[0] = 0.5;

  for (long jn = 1; jn < twowaves; jn++)
    {
      const auto zsqp = 1.0 / std::sqrt((double) (jn + jn * jn));
      zan *= std::sqrt(1.0 - 1.0 / (4 * jn * jn));

      zcospar = std::cos(lat * jn);
      zsinpar = std::sin(lat * jn) * jn * zsqp;
      zcosfak = 1.0;

      for (long jk = 2; jk < jn; jk += 2)
        {
          jnmjk = jn - jk;
          zcosfak *= (jk - 1.0) * (jn + jnmjk + 2.0) / (jk * (jn + jnmjk + 1.0));
          zsinfak = zcosfak * zsqp * (double) jnmjk;
          zcospar += zcosfak * std::cos(lat * jnmjk);
          zsinpar += zsinfak * std::sin(lat * jnmjk);
        }

      // code for jk == jn

      if ((jn & 1) == 0)
        {
          zcosfak *= (double) ((jn - 1) * (jn + 2)) / (double) (jn * (jn + 1));
          zcospar += zcosfak * 0.5;
        }
      ztemp1[jn] = zan * zcospar;
      ztemp2[jn - 1] = zan * zsinpar;
    }

  memcpy(pnm, ztemp1, waves * sizeof(double));
  pnm += waves;
  memcpy(pnm, ztemp2, waves * sizeof(double));
  pnm += waves;

  if (needHnm)
    {
      hnm[0] = 0.0;
      for (long jn = 1; jn < waves; jn++)
        hnm[jn] = jn * (pmu * ztemp1[jn] - std::sqrt((jn + jn + 1.0) / (jn + jn - 1.0)) * ztemp1[jn - 1]);

      hnm += waves;

      hnm[0] = pmu * ztemp2[0];

      for (long jn = 1; jn < waves; jn++)
        hnm[jn] = (jn + 1) * pmu * ztemp2[jn]
                  - std::sqrt(((jn + jn + 3.0) * ((jn + 1) * (jn + 1) - 1.0)) / (jn + jn + 1.0)) * ztemp2[jn - 1];

      hnm += waves;
    }

  for (long jm = 2; jm < waves; jm++)
    {
      pnm[0] = std::sqrt(1.0 + 1.0 / (jm + jm)) * zcos2 * ztemp2[0];
      if (needHnm) hnm[0] = jm * pmu * pnm[0];
#if defined(CRAY)
#pragma _CRI novector
#endif
#if defined(__uxp__)
#pragma loop scalar
#endif
      for (long jn = 1; jn < (twowaves - jm); jn++)
        {
          zq = jm + jm + jn - 1;
          zwm2 = zq + jn;
          zw = zwm2 + 2.0;
          zwq = zw * zq;
          zq2m1 = zq * zq - 1.0;
          zwm2q2 = zwm2 * zq2m1;
          z2q2 = zq2m1 * 2;
          zcnm = std::sqrt((zwq * (zq - 2.0)) / (zwm2q2 - z2q2));
          zdnm = std::sqrt((zwq * (jn + 1.0)) / zwm2q2);
          zenm = std::sqrt(zw * jn / ((zq + 1.0) * zwm2));
          pnm[jn] = zcnm * ztemp1[jn] - pmu * (zdnm * ztemp1[jn + 1] - zenm * pnm[jn - 1]);
          if (needHnm) hnm[jn] = (jm + jn) * pmu * pnm[jn] - std::sqrt(zw * jn * (zq + 1) / zwm2) * pnm[jn - 1];
        }
      memcpy(ztemp1, ztemp2, twowaves * sizeof(double));
      memcpy(ztemp2, pnm, twowaves * sizeof(double));
      pnm += waves;
      if (needHnm) hnm += waves;
    }
}

void
after_legini_full(long ntr, long nlat, double *restrict poli, double *restrict pold, double *restrict pdev, double *restrict pol2,
                  double *restrict pol3, double *restrict coslat)
{
  auto dimsp = (ntr + 1) * (ntr + 2);
  auto waves = ntr + 1;
  const auto twowaves = waves << 1;

  Varray<double> gmu(nlat), gwt(nlat);
  gaussian_latitudes(nlat, gmu.data(), gwt.data());

  auto needHnm = (pdev != nullptr) || (pol2 != nullptr);

#ifdef _OPENMP
  Varray2D<double> pnm_2(Threading::ompNumThreads);
  Varray2D<double> hnm_2(Threading::ompNumThreads);
  Varray2D<double> ztemp1_2(Threading::ompNumThreads);
  Varray2D<double> ztemp2_2(Threading::ompNumThreads);
  for (long i = 0; i < Threading::ompNumThreads; ++i) pnm_2[i].resize(dimsp);
  if (needHnm)
    for (long i = 0; i < Threading::ompNumThreads; ++i) hnm_2[i].resize(dimsp);
  for (long i = 0; i < Threading::ompNumThreads; ++i) ztemp1_2[i].resize(twowaves);
  for (long i = 0; i < Threading::ompNumThreads; ++i) ztemp2_2[i].resize(twowaves);
#else
  Varray<double> pnm(dimsp), hnm, ztemp1(twowaves), ztemp2(twowaves);
  if (needHnm) hnm.resize(dimsp);
#endif

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
  for (long jgl = 0; jgl < nlat; ++jgl)
    {
#ifdef _OPENMP
      const auto ompthID = cdo_omp_get_thread_num();
      auto &pnm = pnm_2[ompthID];
      auto &hnm = hnm_2[ompthID];
      auto &ztemp1 = ztemp1_2[ompthID];
      auto &ztemp2 = ztemp2_2[ompthID];
#endif
      const auto gmusq = 1.0 - gmu[jgl] * gmu[jgl];
      coslat[jgl] = std::sqrt(gmusq);

      phcs(needHnm, pnm.data(), hnm.data(), waves, gmu[jgl], ztemp1.data(), ztemp2.data());

      const auto zgwt = gwt[jgl];
      const auto zrafgmusqr = 1.0 / (PlanetRadius * gmusq);
      const auto zradsqrtgmusqr = 1.0 / (-PlanetRadius * std::sqrt(gmusq));

      auto jsp = jgl;
      for (long jm = 0; jm < waves; ++jm)
        for (long jn = 0; jn < (waves - jm); ++jn)
          {
            if (poli) poli[jsp] = pnm[jm * waves + jn] * 2.0;
            if (pold) pold[jsp] = pnm[jm * waves + jn] * zgwt;
            if (pdev) pdev[jsp] = hnm[jm * waves + jn] * 2.0 * zradsqrtgmusqr;
            if (pol2) pol2[jsp] = hnm[jm * waves + jn] * zgwt * zrafgmusqr;
            if (pol3) pol3[jsp] = pnm[jm * waves + jn] * zgwt * jm * zrafgmusqr;
            jsp += nlat;
          }
    }
}
/*
void
after_legini(long ntr, long nlat, double *restrict poli, double *restrict pold, double *restrict coslat)
{
  long waves = ntr + 1;                     // used in omp loop
  long dimpnm = (ntr + 1) * (ntr + 4) / 2;  // used in omp loop

#ifdef _OPENMP
  Varray2D<double> pnm_2(Threading::ompNumThreads);
  Varray2D<double> work_2(Threading::ompNumThreads);
  for (long i = 0; i < Threading::ompNumThreads; ++i) pnm_2[i].resize(dimpnm);
  for (long i = 0; i < Threading::ompNumThreads; ++i) work_2[i].resize(3 * waves);
#else
  Varray<double> pnm(dimpnm), work(3 * waves);
#endif

  Varray<double> gmu(nlat), gwt(nlat);
  gaussian_latitudes(nlat, gmu.data(), gwt.data());
  for (long jgl = 0; jgl < nlat; ++jgl) gwt[jgl] *= 0.5;

  for (long jgl = 0; jgl < nlat; ++jgl) coslat[jgl] = std::sqrt(1.0 - gmu[jgl] * gmu[jgl]);

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
  for (long jgl = 0; jgl < nlat / 2; jgl++)
    {
#ifdef _OPENMP
      const auto ompthID = cdo_omp_get_thread_num();
      auto pnm = pnm_2[ompthID].data();
      auto work = work_2[ompthID].data();
#endif
      const auto zgwt = gwt[jgl];

      jspleg1(&pnm[0], gmu[jgl], ntr, &work[0]);

      long latn = jgl;
      long lats;
      long isp = 0;
      double is;
      for (long jm = 0; jm < waves; ++jm)
        {
          for (long jn = 0; jn < waves - jm; ++jn)
            {
              is = (jn + 1) % 2 * 2 - 1;
              lats = latn - jgl + nlat - jgl - 1;
              poli[latn] = pnm[isp];
              pold[latn] = pnm[isp] * zgwt;
              poli[lats] = pnm[isp] * is;
              pold[lats] = pnm[isp] * zgwt * is;
              latn += nlat;
              isp++;
            }
          isp++;
        }
    }
}
*/
static inline void
sp2fc_kernel(long nlat, const double *pol, const double *sal, double *restrict far, double *restrict fai)
{
  const auto sar = sal[0];
  const auto sai = sal[1];
  for (long lat = 0; lat < nlat; lat++)
    {
      far[lat] += pol[lat] * sar;
      fai[lat] += pol[lat] * sai;
    }
}

void
sp2fc(const double *sa, double *fa, const double *poli, long nlev, long nlat, long nfc, long nt)
{
  long ntp1 = nt + 1;
  long nsp2 = (nt + 1) * (nt + 2);

  std::vector<long> cumindex(ntp1);
  cumindex[0] = 0;
  for (long jmm = 1; jmm < ntp1; jmm++) cumindex[jmm] = cumindex[jmm - 1] + (ntp1 - jmm + 1);

  if (nlev >= Threading::ompNumThreads)
    {
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
      for (long lev = 0; lev < nlev; lev++)
        {
          auto sal = sa + lev * nsp2;
          auto fal = fa + lev * nfc * nlat;
          memset(fal, 0, nfc * nlat * sizeof(double));

          for (long jmm = 0; jmm < ntp1; jmm++)
            {
              auto polt = poli + cumindex[jmm] * nlat;
              auto salt = sal + cumindex[jmm] * 2;
              auto far = fal + jmm * 2 * nlat;
              auto fai = far + nlat;
              for (long jfc = 0; jfc < (ntp1 - jmm); jfc++) { sp2fc_kernel(nlat, polt + jfc * nlat, salt + jfc * 2, far, fai); }
            }
        }
    }
  else
    {
      for (long lev = 0; lev < nlev; lev++)
        {
          auto sal = sa + lev * nsp2;
          auto fal = fa + lev * nfc * nlat;
          memset(fal, 0, nfc * nlat * sizeof(double));

#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(dynamic)
#endif
          for (long jmm = 0; jmm < ntp1; jmm++)
            {
              auto polt = poli + cumindex[jmm] * nlat;
              auto salt = sal + cumindex[jmm] * 2;
              auto far = fal + jmm * 2 * nlat;
              auto fai = far + nlat;
              for (long jfc = 0; jfc < (ntp1 - jmm); jfc++) { sp2fc_kernel(nlat, polt + jfc * nlat, salt + jfc * 2, far, fai); }
            }
        }
    }
}

static inline void
fc2sp_kernel(long nlat, const double *pol, const double *far, const double *fai, double *sal)
{
  double sar = 0.0;
  double sai = 0.0;
#ifdef HAVE_OPENMP4
#pragma omp simd reduction(+ : sar) reduction(+ : sai)
#endif
  for (long lat = 0; lat < nlat; lat++)
    {
      sar += pol[lat] * far[lat];
      sai += pol[lat] * fai[lat];
    }
  sal[0] = sar;
  sal[1] = sai;
}

void
fc2sp(const double *fa, double *sa, const double *poli, long nlev, long nlat, long nfc, long nt)
{
  long ntp1 = nt + 1;
  long nsp2 = (nt + 1) * (nt + 2);

  std::vector<long> cumindex(ntp1);
  cumindex[0] = 0;
  for (long jmm = 1; jmm < ntp1; jmm++) cumindex[jmm] = cumindex[jmm - 1] + (ntp1 - jmm + 1);

  if (nlev >= Threading::ompNumThreads)
    {
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
      for (long lev = 0; lev < nlev; lev++)
        {
          auto fal = fa + lev * nfc * nlat;
          auto sal = sa + lev * nsp2;

          for (long jmm = 0; jmm < ntp1; jmm++)
            {
              auto polt = poli + cumindex[jmm] * nlat;
              auto salt = sal + cumindex[jmm] * 2;
              auto far = fal + jmm * 2 * nlat;
              auto fai = far + nlat;
              for (long jfc = 0; jfc < (ntp1 - jmm); jfc++) { fc2sp_kernel(nlat, polt + jfc * nlat, far, fai, salt + jfc * 2); }
            }
        }
    }
  else
    {
      for (long lev = 0; lev < nlev; lev++)
        {
          auto fal = fa + lev * nfc * nlat;
          auto sal = sa + lev * nsp2;

#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(dynamic)
#endif
          for (long jmm = 0; jmm < ntp1; jmm++)
            {
              auto polt = poli + cumindex[jmm] * nlat;
              auto salt = sal + cumindex[jmm] * 2;
              auto far = fal + jmm * 2 * nlat;
              auto fai = far + nlat;
              for (long jfc = 0; jfc < (ntp1 - jmm); jfc++) { fc2sp_kernel(nlat, polt + jfc * nlat, far, fai, salt + jfc * 2); }
            }
        }
    }
}

/* ======================================== */
/* Convert Spectral Array to new truncation */
/* ======================================== */

void
sp2sp(const double *arrayIn, long truncIn, double *arrayOut, long truncOut)
{
  if (truncOut <= truncIn)
    {
      for (long n = 0; n <= truncOut; ++n)
        {
          for (long m = n; m <= truncOut; ++m)
            {
              *arrayOut++ = *arrayIn++;
              *arrayOut++ = *arrayIn++;
            }
          arrayIn += 2 * (truncIn - truncOut);
        }
    }
  else
    {
      for (long n = 0; n <= truncIn; ++n)
        {
          for (long m = n; m <= truncIn; ++m)
            {
              *arrayOut++ = *arrayIn++;
              *arrayOut++ = *arrayIn++;
            }
          for (long m = truncIn + 1; m <= truncOut; ++m)
            {
              *arrayOut++ = 0.0;
              *arrayOut++ = 0.0;
            }
        }
      for (long n = truncIn + 1; n <= truncOut; ++n)
        for (long m = n; m <= truncOut; ++m)
          {
            *arrayOut++ = 0.0;
            *arrayOut++ = 0.0;
          }
    }
}

/* ======================================== */
/* Cut spectral wave numbers                */
/* ======================================== */

void
spcut(const double *arrayIn, double *arrayOut, long trunc, const int *waves)
{
  for (long n = 0; n <= trunc; ++n)
    {
      for (long m = n; m <= trunc; ++m)
        {
          if (waves[m])
            {
              *arrayOut++ = *arrayIn++;
              *arrayOut++ = *arrayIn++;
            }
          else
            {
              *arrayOut++ = 0.0;
              *arrayOut++ = 0.0;
              arrayIn++;
              arrayIn++;
            }
        }
    }
}
