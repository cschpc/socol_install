/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Specinfo specinfo  Spectral information
*/

#include "process_int.h"
#include <mpim_grid.h>

#define NTR2NSP(ntr) ((ntr + 1) * (ntr + 2))
#define NSP2NTR(nsp) ((long) ((((std::sqrt((double) (4 * nsp + 1))) - 3) / 2)))
#define NGP2NLEVEL(ngp) ((long) (log10(((double) ngp) / 80.) / std::log10(4.)))
#define NGP_ICON(nrooti, nlevel) ((long) (20 * nrooti * nrooti * ipow(4, nlevel)))
/*#define NGP_GME(ni)           ((ni+1)*(ni+1)*10)*/
#define NGP_GME(ni) (2 + ni * ni * 10)
#define NGP2NI(ngp) ((long) std::sqrt((double) ngp / 10.) - 1)

static void
fac(long nlonin, long *nlonout, int *ierr)
{
  long n2 = 0;
  long n3 = 0;
  long n5 = 0;

  long m = nlonin;

  while (m % 2 == 0)
    {
      m = m / 2;
      n2++;
    }
  while (m % 3 == 0)
    {
      m = m / 3;
      n3++;
    }
  while (m % 5 == 0)
    {
      m = m / 5;
      n5++;
    }

  if (m == 1)
    {
      *nlonout = nlonin;
      *ierr = 0;
    }
  else
    {
      *nlonout = nlonin + 1;
      *ierr = 1;
    }

  return;
}

static long
nlat2nlon(long nlat)
{
  if (nlat == 0) cdo_abort("nlat = 0!");

  long nlon = 2 * nlat;

  long m;
  int ierr;
  fac(nlon, &m, &ierr);
  /* adjust till fft is possible */
  while (ierr != 0)
    {
      nlon = m;
      /* correct here nlon so that nlat keeps always even */
      while (nlon % 4 != 0) nlon++;
      fac(nlon, &m, &ierr);
    }

  return nlon;
}

long
ngp2ntr(long ngp)
{
  long ntr = (long) std::lround(std::sqrt(0.25 + ngp) - 1.5);
  long nlonl = nlat_to_nlon(ntr_to_nlat_linear(ntr));
  long nlatl = nlonl / 2;

  ntr = (2 * nlatl - 1) / 2;

  return ntr;
}

static long
ipow(long i1, long i2)
{
  long i3 = 1;

  for (long i = 0; i < i2; ++i) i3 *= i1;

  return i3;
}

static constexpr long NiMax = 12;

static void
lookup_ni(long nsp, long *nroot, long *ni)
{
  long tbl2[NiMax], tbl3[NiMax], tbl5[NiMax];
  long d2 = 0, n2 = 0, d3 = 0, n3 = 0, d5 = 0, n5 = 0;

  for (long i = 0; i < NiMax; ++i)
    {
      tbl2[i] = 10 * 2 * 2 * ipow(4, (i + 1)) + 2;
      tbl3[i] = 10 * 3 * 3 * ipow(4, (i + 1)) + 2;
      tbl5[i] = 10 * 5 * 5 * ipow(4, (i + 1)) + 2;
    }

  for (long i = 0; i < NiMax; ++i)
    if (tbl2[i] >= nsp)
      {
        n2 = i;
        d2 = tbl2[n2] - nsp;
        break;
      }

  for (long i = 0; i < NiMax; ++i)
    if (tbl3[i] >= nsp)
      {
        n3 = i;
        d3 = tbl3[n3] - nsp;
        break;
      }

  for (long i = 0; i < NiMax; ++i)
    if (tbl5[i] >= nsp)
      {
        n5 = i;
        d5 = tbl5[n5] - nsp;
        break;
      }

  long d = d2;
  if (d3 < d) d = d3;
  if (d5 < d) d = d5;

  if (d == d2)
    {
      *nroot = 2;
      *ni = 2 * ipow(2, n2 + 1);
    }
  else if (d == d3)
    {
      *nroot = 3;
      *ni = 3 * ipow(2, n3 + 1);
    }
  else if (d == d5)
    {
      *nroot = 5;
      *ni = 5 * ipow(2, n5 + 1);
    }
}

static void
lookup_rl(long nsp, long *nroot, long *nlevel)
{
  long tbl2[NiMax], tbl3[NiMax], tbl5[NiMax];
  long d2 = 0, n2 = 0, d3 = 0, n3 = 0, d5 = 0, n5 = 0;

  for (long i = 0; i < NiMax; ++i)
    {
      tbl2[i] = 20 * 2 * 2 * ipow(4, (i + 1));
      tbl3[i] = 20 * 3 * 3 * ipow(4, (i + 1));
      tbl5[i] = 20 * 5 * 5 * ipow(4, (i + 1));
    }

  for (long i = 0; i < NiMax; ++i)
    if (tbl2[i] >= nsp)
      {
        n2 = i;
        d2 = tbl2[n2] - nsp;
        break;
      }

  for (long i = 0; i < NiMax; ++i)
    if (tbl3[i] >= nsp)
      {
        n3 = i;
        d3 = tbl3[n3] - nsp;
        break;
      }

  for (long i = 0; i < NiMax; ++i)
    if (tbl5[i] >= nsp)
      {
        n5 = i;
        d5 = tbl5[n5] - nsp;
        break;
      }

  long d = d2;
  if (d3 < d) d = d3;
  if (d5 < d) d = d5;

  if (d == d2)
    {
      *nroot = 2;
      *nlevel = n2 + 1;
    }
  else if (d == d3)
    {
      *nroot = 3;
      *nlevel = n3 + 1;
    }
  else if (d == d5)
    {
      *nroot = 5;
      *nlevel = n5 + 1;
    }
}

void *
Specinfo(void *process)
{
  char arg[128], *parg;
  bool nout1 = false, nout2 = false, nout3 = false;
  long i;
  long ntr1 = 0, nsp1 = 0, nlat1 = 0, nlon1 = 0, ngp1 = 0, ni1 = 0, ngp_gme1 = 0;
  long ntr2 = 0, nsp2 = 0, nlat2 = 0, nlon2 = 0, ngp2 = 0, ni2 = 0, ngp_gme2 = 0;
  long ntr3 = 0, nsp3 = 0, nlat3 = 0, nlon3 = 0, ngp3 = 0, ni3 = 0, ngp_gme3 = 0;
  long nlevel1 = 0, nlevel2 = 0, nlevel3 = 0;
  long ngp_icon1 = 0, ngp_icon2 = 0, ngp_icon3 = 0;
  long nrootg1 = 0, nrootg2 = 0, nrootg3 = 0;
  long nrooti1 = 0, nrooti2 = 0, nrooti3 = 0;

  cdo_initialize(process);

  operator_input_arg("Txx, TLxx, NLON=xx, NLAT=xx, NIxx or ICONRyyLxx");

  const long len = cdo_operator_argv(0).size();

  if ((len + 1) >= 128) cdo_abort("Parameter string too large!");

  for (i = 0; i < len; ++i) arg[i] = toupper(cdo_operator_argv(0)[i]);
  arg[len] = 0;

  if (arg[0] == 'T' && arg[1] == 'L')
    {
      parg = &arg[2];
      if (*parg == '=') parg++;
      if (!isdigit((int) *parg)) cdo_abort("Wrong parameter: %s", arg);
      ntr2 = atoi(parg);
      nsp2 = NTR2NSP(ntr2);
      nlat2 = ntr_to_nlat_linear(ntr2);
      nlon2 = nlat_to_nlon(nlat2);
      ngp2 = nlon2 * nlat2;

      lookup_ni(nsp2, &nrootg2, &ni2);
      lookup_rl(nsp2, &nrooti2, &nlevel2);

      nout2 = true;
    }
  else if (arg[0] == 'T' && arg[1] == 'C')
    {
      parg = &arg[2];
      if (*parg == '=') parg++;
      if (!isdigit((int) *parg)) cdo_abort("Wrong parameter: %s", arg);
      ntr3 = atoi(parg);
      nsp3 = NTR2NSP(ntr3);
      nlat3 = ntr_to_nlat_cubic(ntr3);
      nlon3 = nlat_to_nlon(nlat3);
      ngp3 = nlon3 * nlat3;

      lookup_ni(nsp3, &nrootg3, &ni3);
      lookup_rl(nsp3, &nrooti3, &nlevel3);

      nout3 = true;
    }
  else if (arg[0] == 'T')
    {
      parg = &arg[1];
      if (*parg == '=') parg++;
      if (!isdigit((int) *parg)) cdo_abort("Wrong parameter: %s", arg);
      ntr1 = atoi(parg);
      nsp1 = NTR2NSP(ntr1);
      nlat1 = ntr_to_nlat(ntr1);
      nlon1 = nlat_to_nlon(nlat1);
      ngp1 = nlon1 * nlat1;

      lookup_ni(nsp1, &nrootg1, &ni1);
      lookup_rl(nsp1, &nrooti1, &nlevel1);

      nout1 = true;
    }
  else if (arg[0] == 'N' && arg[1] == 'I')
    {
      parg = &arg[2];
      if (*parg == '=') parg++;
      if (!isdigit((int) *parg)) cdo_abort("Wrong parameter: %s", arg);
      ni1 = atoi(parg);
      ni2 = ni1;
      ngp_gme1 = NGP_GME(ni1);
      ngp_gme2 = NGP_GME(ni2);

      ntr1 = ngp2ntr(ngp_gme1);

      nsp1 = NTR2NSP(ntr1);
      ntr1 = NSP2NTR(nsp1);
      ntr2 = ntr1;

      nlat1 = ntr_to_nlat(ntr1);
      nlon1 = nlat_to_nlon(nlat1);
      nlat1 = nlon1 / 2;

      nlat2 = ntr_to_nlat_linear(ntr2);
      nlon2 = nlat_to_nlon(nlat2);
      nlat2 = nlon2 / 2;

      /* lookup_ni(nsp1, &nrootg1, &ni1); */
      lookup_rl(nsp1, &nrooti1, &nlevel1);

      nrootg2 = nrootg1;
      ni2 = ni1;
      nrooti2 = nrooti1;
      nlevel2 = nlevel1;

      nout1 = true;
      nout2 = true;
    }
  else if (arg[0] == 'N' && arg[1] == 'L' && arg[2] == 'O' && arg[3] == 'N')
    {
      parg = &arg[4];
      if (*parg == '=') parg++;
      if (!isdigit((int) *parg)) cdo_abort("Wrong parameter: %s", arg);
      nlon1 = atoi(parg);
      nlon2 = nlon1;
      nlon3 = nlon1;
      nlat1 = nlon1 / 2;
      nlat2 = nlon2 / 2;
      nlat3 = nlon3 / 2;
      nlon1 = nlat2nlon(nlat1);
      nlon2 = nlat2nlon(nlat2);
      nlon3 = nlat2nlon(nlat3);
      nlat1 = nlon1 / 2;
      nlat2 = nlon2 / 2;
      nlat3 = nlon3 / 2;
      ntr1 = (nlat1 * 2 - 1) / 3;
      ntr2 = (nlat2 * 2 - 1) / 2;
      ntr3 = (nlat3 * 2 - 1) / 4;
      ngp1 = nlon1 * nlat1;
      ngp2 = nlon2 * nlat2;
      ngp3 = nlon3 * nlat3;

      nsp1 = NTR2NSP(ntr1);
      nsp2 = NTR2NSP(ntr2);
      nsp3 = NTR2NSP(ntr3);

      lookup_ni(nsp1, &nrootg1, &ni1);
      lookup_rl(nsp1, &nrooti1, &nlevel1);

      lookup_ni(nsp2, &nrootg2, &ni2);
      lookup_rl(nsp2, &nrooti2, &nlevel2);

      lookup_ni(nsp3, &nrootg3, &ni3);
      lookup_rl(nsp3, &nrooti3, &nlevel3);

      nout1 = true;
      nout2 = true;
      nout3 = true;
    }
  else if (arg[0] == 'N' && arg[1] == 'L' && arg[2] == 'A' && arg[3] == 'T')
    {
      parg = &arg[4];
      if (*parg == '=') parg++;
      if (!isdigit((int) *parg)) cdo_abort("Wrong parameter: %s", arg);
      nlat1 = atoi(parg);
      nlat2 = nlat1;
      nlat3 = nlat1;
      nlon1 = nlat2nlon(nlat1);
      nlon2 = nlat2nlon(nlat2);
      nlon3 = nlat2nlon(nlat3);
      nlat1 = nlon1 / 2;
      nlat2 = nlon2 / 2;
      nlat3 = nlon3 / 2;
      ntr1 = (nlat1 * 2 - 1) / 3;
      ntr2 = (nlat2 * 2 - 1) / 2;
      ntr3 = (nlat3 * 2 - 1) / 4;
      ngp1 = nlon1 * nlat1;
      ngp2 = nlon2 * nlat2;
      ngp3 = nlon3 * nlat3;

      nsp1 = NTR2NSP(ntr1);
      nsp2 = NTR2NSP(ntr2);
      nsp3 = NTR2NSP(ntr3);

      lookup_ni(nsp1, &nrootg1, &ni1);
      lookup_rl(nsp1, &nrooti1, &nlevel1);

      lookup_ni(nsp2, &nrootg2, &ni2);
      lookup_rl(nsp2, &nrooti2, &nlevel2);

      lookup_ni(nsp3, &nrootg3, &ni3);
      lookup_rl(nsp3, &nrooti3, &nlevel3);

      nout1 = true;
      nout2 = true;
      nout3 = true;
    }
  else if (arg[0] == 'N')
    {
      parg = &arg[1];
      if (*parg == '=') parg++;
      if (!isdigit((int) *parg)) cdo_abort("Wrong parameter: %s", arg);
      nlat1 = 2 * atoi(parg);
      nlat2 = nlat1;
      nlat3 = nlat1;
      nlon1 = nlat2nlon(nlat1);
      nlon2 = nlat2nlon(nlat2);
      nlon3 = nlat2nlon(nlat3);
      nlat1 = nlon1 / 2;
      nlat2 = nlon2 / 2;
      nlat3 = nlon3 / 2;
      ntr1 = (nlat1 * 2 - 1) / 3;
      ntr2 = (nlat2 * 2 - 1) / 2;
      ntr3 = (nlat3 * 2 - 1) / 4;
      ngp1 = nlon1 * nlat1;
      ngp2 = nlon2 * nlat2;
      ngp3 = nlon3 * nlat3;

      nsp1 = NTR2NSP(ntr1);
      nsp2 = NTR2NSP(ntr2);
      nsp3 = NTR2NSP(ntr3);

      lookup_ni(nsp1, &nrootg1, &ni1);
      lookup_rl(nsp1, &nrooti1, &nlevel1);

      lookup_ni(nsp2, &nrootg2, &ni2);
      lookup_rl(nsp2, &nrooti2, &nlevel2);

      lookup_ni(nsp3, &nrootg3, &ni3);
      lookup_rl(nsp3, &nrooti3, &nlevel3);

      nout1 = true;
      nout2 = true;
      nout3 = true;
    }
  else if (arg[0] == 'O')
    {
      parg = &arg[1];
      if (*parg == '=') parg++;
      if (!isdigit((int) *parg)) cdo_abort("Wrong parameter: %s", arg);
      nlat1 = 2 * atoi(parg);
      nlat2 = nlat1;
      nlat3 = nlat1;
      nlon1 = nlat2nlon(nlat1) + 16;
      nlon2 = nlat2nlon(nlat2) + 16;
      nlon3 = nlat2nlon(nlat3) + 16;
      nlat1 = nlon1 / 2;
      nlat2 = nlon2 / 2;
      nlat3 = nlon3 / 2;
      ntr1 = (nlat1 * 2 - 1) / 3;
      ntr2 = (nlat2 * 2 - 1) / 2;
      ntr3 = (nlat3 * 2 - 1) / 4;
      ngp1 = nlon1 * nlat1;
      ngp2 = nlon2 * nlat2;
      ngp3 = nlon3 * nlat3;

      nsp1 = NTR2NSP(ntr1);
      nsp2 = NTR2NSP(ntr2);
      nsp3 = NTR2NSP(ntr3);

      lookup_ni(nsp1, &nrootg1, &ni1);
      lookup_rl(nsp1, &nrooti1, &nlevel1);

      lookup_ni(nsp2, &nrootg2, &ni2);
      lookup_rl(nsp2, &nrooti2, &nlevel2);

      lookup_ni(nsp3, &nrootg3, &ni3);
      lookup_rl(nsp3, &nrooti3, &nlevel3);

      nout1 = true;
      nout2 = true;
      nout3 = true;
    }
  else if (arg[0] == 'I' && arg[1] == 'C' && arg[2] == 'O' && arg[3] == 'N')
    {
      parg = &arg[4];
      if (*parg != 'R') cdo_abort("Wrong parameter: %s", arg);
      parg++;
      if (!isdigit((int) *parg)) cdo_abort("Wrong parameter: %s", arg);
      nrooti1 = atoi(parg);
      nrooti2 = nrooti1;
      while (isdigit((int) *parg)) parg++;
      if (*parg != 'B') cdo_abort("Wrong parameter: %s", arg);
      parg++;
      if (!isdigit((int) *parg)) cdo_abort("Wrong parameter: %s", arg);
      nlevel1 = atoi(parg);
      nlevel2 = nlevel1;
      ngp_icon1 = NGP_ICON(nrooti1, nlevel1);
      ngp_icon2 = NGP_ICON(nrooti1, nlevel2);

      ntr1 = ngp2ntr(ngp_icon1);
      nsp1 = NTR2NSP(ntr1);
      ntr1 = NSP2NTR(nsp1);
      ntr2 = ntr1;

      nlat1 = ntr_to_nlat(ntr1);
      nlon1 = nlat_to_nlon(nlat1);
      nlat1 = nlon1 / 2;

      nlat2 = ntr_to_nlat_linear(ntr2);
      nlon2 = nlat_to_nlon(nlat2);
      nlat2 = nlon2 / 2;

      lookup_ni(nsp1, &nrootg1, &ni1);
      /* lookup_rl(nsp1, &nrooti1, &nlevel1);*/

      nrootg2 = nrootg1;
      ni2 = ni1;
      nrooti2 = nrooti1;
      nlevel2 = nlevel1;

      nout1 = true;
      nout2 = true;
    }
  else
    cdo_abort("Unsupported parameter: %s", arg);

  nsp1 = NTR2NSP(ntr1);
  nsp2 = NTR2NSP(ntr2);
  nsp3 = NTR2NSP(ntr3);
  ngp1 = nlon1 * nlat1;
  ngp2 = nlon2 * nlat2;
  ngp3 = nlon3 * nlat3;
  ngp_gme1 = NGP_GME(ni1);
  ngp_gme2 = NGP_GME(ni2);
  ngp_gme3 = NGP_GME(ni3);
  ngp_icon1 = NGP_ICON(nrooti1, nlevel1);
  ngp_icon2 = NGP_ICON(nrooti2, nlevel2);
  ngp_icon3 = NGP_ICON(nrooti3, nlevel3);

  fprintf(stdout, "truncation     nsp  nlon  nlat      ngp  gme    ngp_gme  icon   ngp_icon\n");

  if (nout2)
    fprintf(stdout, "   TL%-4ld %8ld %5ld %5ld %8ld  ni%ld %8ld  R%ldB%02ld  %8ld\n", ntr2, nsp2, nlon2, nlat2, ngp2, ni2, ngp_gme2,
            nrooti2, nlevel2, ngp_icon2);

  if (nout1)
    fprintf(stdout, "   TQ%-4ld %8ld %5ld %5ld %8ld  ni%ld %8ld  R%ldB%02ld  %8ld\n", ntr1, nsp1, nlon1, nlat1, ngp1, ni1, ngp_gme1,
            nrooti1, nlevel1, ngp_icon1);

  if (nout3)
    fprintf(stdout, "   TC%-4ld %8ld %5ld %5ld %8ld  ni%ld %8ld  R%ldB%02ld  %8ld\n", ntr3, nsp3, nlon3, nlat3, ngp3, ni3, ngp_gme3,
            nrooti3, nlevel3, ngp_icon3);

  cdo_finish();

  return nullptr;
}
