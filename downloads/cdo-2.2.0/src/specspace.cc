/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include "cdo_output.h"
#include "cdo_fctrans.h"
#include "specspace.h"
#include <mpim_grid.h>

void
grid2spec(const SP_Transformation &spTrans, int gridIDin, const Varray<double> &arrayIn, int gridIDout, Varray<double> &arrayOut)
{
  const auto &fcTrans = spTrans.fcTrans;
  const long nlev = 1;
  const long ntr = gridInqTrunc(gridIDout);
  const long nlon = gridInqXsize(gridIDin);
  const long nlat = gridInqYsize(gridIDin);
  const long waves = ntr + 1;
  const long nfc = waves * 2;

  Varray<double> fpwork(nlat * nfc * nlev);

  if (fcTrans.use_fftw)
    gp2fc(arrayIn.data(), fpwork.data(), nlat, nlon, nlev, nfc);
  else
    gp2fc(fcTrans.vtrig.data(), fcTrans.ifax, arrayIn.data(), fpwork.data(), nlat, nlon, nlev, nfc);

  fc2sp(fpwork.data(), arrayOut.data(), spTrans.pold.data(), nlev, nlat, nfc, ntr);
}

void
spec2grid(const SP_Transformation &spTrans, int gridIDin, const Varray<double> &arrayIn, int gridIDout, Varray<double> &arrayOut)
{
  const auto &fcTrans = spTrans.fcTrans;
  const long nlev = 1;
  const long ntr = gridInqTrunc(gridIDin);
  const long nlon = gridInqXsize(gridIDout);
  const long nlat = gridInqYsize(gridIDout);
  const long waves = ntr + 1;
  const long nfc = waves * 2;

  Varray<double> fpwork(nlat * nfc * nlev);

  sp2fc(arrayIn.data(), fpwork.data(), spTrans.poli.data(), nlev, nlat, nfc, ntr);

  if (fcTrans.use_fftw)
    fc2gp(fpwork.data(), arrayOut.data(), nlat, nlon, nlev, nfc);
  else
    fc2gp(fcTrans.vtrig.data(), fcTrans.ifax, fpwork.data(), arrayOut.data(), nlat, nlon, nlev, nfc);
}

void
four2spec(const SP_Transformation &spTrans, int gridIDin, const Varray<double> &arrayIn, int gridIDout, Varray<double> &arrayOut)
{
  (void) gridIDin;
  const long nlev = 1;
  const long ntr = gridInqTrunc(gridIDout);
  const long nlat = spTrans.nlat;
  const long waves = ntr + 1;
  const long nfc = waves * 2;

  fc2sp(arrayIn.data(), arrayOut.data(), spTrans.pold.data(), nlev, nlat, nfc, ntr);
}

void
spec2four(const SP_Transformation &spTrans, int gridIDin, const Varray<double> &arrayIn, int gridIDout, Varray<double> &arrayOut)
{
  const long nlev = 1;
  const long ntr = gridInqTrunc(gridIDin);
  long nfc = gridInqSize(gridIDout);
  const long nlat = nfc_to_nlat(nfc, ntr);
  const long waves = ntr + 1;
  nfc = waves * 2;

  sp2fc(arrayIn.data(), arrayOut.data(), spTrans.poli.data(), nlev, nlat, nfc, ntr);
}

void
four2grid(const FC_Transformation &fcTrans, int gridIDin, const Varray<double> &arrayIn, int gridIDout, Varray<double> &arrayOut)
{
  const long nlev = 1;
  const long ntr = gridInqTrunc(gridIDin);
  const long nlon = gridInqXsize(gridIDout);
  const long nlat = gridInqYsize(gridIDout);
  const long waves = ntr + 1;
  const long nfc = waves * 2;

  if (fcTrans.use_fftw)
    fc2gp(arrayIn.data(), arrayOut.data(), nlat, nlon, nlev, nfc);
  else
    fc2gp(fcTrans.vtrig.data(), fcTrans.ifax, arrayIn.data(), arrayOut.data(), nlat, nlon, nlev, nfc);
}

void
grid2four(const FC_Transformation &fcTrans, int gridIDin, const Varray<double> &arrayIn, int gridIDout, Varray<double> &arrayOut)
{
  const long nlev = 1;
  const long ntr = gridInqTrunc(gridIDout);
  const long nlon = gridInqXsize(gridIDin);
  const long nlat = gridInqYsize(gridIDin);
  const long waves = ntr + 1;
  const long nfc = waves * 2;

  if (fcTrans.use_fftw)
    gp2fc(arrayIn.data(), arrayOut.data(), nlat, nlon, nlev, nfc);
  else
    gp2fc(fcTrans.vtrig.data(), fcTrans.ifax, arrayIn.data(), arrayOut.data(), nlat, nlon, nlev, nfc);
}

void
spec2spec(int gridIDin, const Varray<double> &arrayIn, int gridIDout, Varray<double> &arrayOut)
{
  const long ntrIn = gridInqTrunc(gridIDin);
  const long ntrOut = gridInqTrunc(gridIDout);

  sp2sp(arrayIn.data(), ntrIn, arrayOut.data(), ntrOut);
}

void
speccut(int gridIDin, const Varray<double> &arrayIn, Varray<double> &arrayOut, Varray<int> &waves)
{
  const long ntr = gridInqTrunc(gridIDin);

  spcut(arrayIn.data(), arrayOut.data(), ntr, waves.data());
}

void
trans_uv2dv(const SP_Transformation &spTrans, long nlev, int gridID1, double *gu, double *gv, int gridID2, double *sd, double *svo)
{
  if (gridInqType(gridID1) != GRID_GAUSSIAN)
    cdo_abort("unexpected grid1 type: %s instead of Gaussian", gridNamePtr(gridInqType(gridID1)));

  if (gridInqType(gridID2) != GRID_SPECTRAL)
    cdo_abort("unexpected grid2 type: %s instead of spectral", gridNamePtr(gridInqType(gridID2)));

  const auto &fcTrans = spTrans.fcTrans;
  const long ntr = gridInqTrunc(gridID2);
  const long nlon = gridInqXsize(gridID1);
  const long nlat = gridInqYsize(gridID1);
  const long waves = ntr + 1;
  const long nfc = waves * 2;

  Varray<double> fpwork1(nlat * nfc * nlev);
  Varray<double> fpwork2(nlat * nfc * nlev);

  if (fcTrans.use_fftw)
    {
      gp2fc(gu, fpwork1.data(), nlat, nlon, nlev, nfc);
      gp2fc(gv, fpwork2.data(), nlat, nlon, nlev, nfc);
    }
  else
    {
      gp2fc(fcTrans.vtrig.data(), fcTrans.ifax, gu, fpwork1.data(), nlat, nlon, nlev, nfc);
      gp2fc(fcTrans.vtrig.data(), fcTrans.ifax, gv, fpwork2.data(), nlat, nlon, nlev, nfc);
    }

  scaluv(fpwork1.data(), spTrans.coslat.data(), nlat, nfc * nlev);
  scaluv(fpwork2.data(), spTrans.coslat.data(), nlat, nfc * nlev);

  uv2dv(fpwork1.data(), fpwork2.data(), sd, svo, spTrans.pol2.data(), spTrans.pol3.data(), nlev, nlat, ntr);
}

void
trans_dv2uv(const SP_Transformation &spTrans, const DV_Transformation &dvTrans, long nlev, int gridID1, double *sd, double *svo,
            int gridID2, double *gu, double *gv)
{
  if (gridInqType(gridID1) != GRID_SPECTRAL)
    cdo_warning("unexpected grid1 type: %s instead of spectral", gridNamePtr(gridInqType(gridID1)));
  if (gridInqType(gridID2) != GRID_GAUSSIAN)
    cdo_warning("unexpected grid2 type: %s instead of Gaussian", gridNamePtr(gridInqType(gridID2)));

  const auto &fcTrans = spTrans.fcTrans;
  const long ntr = gridInqTrunc(gridID1);
  const long nlon = gridInqXsize(gridID2);
  const long nlat = gridInqYsize(gridID2);
  const long waves = ntr + 1;
  const long nfc = waves * 2;
  const long dimsp = (ntr + 1) * (ntr + 2);

  double *su = gu;
  double *sv = gv;

  dv2uv(sd, svo, su, sv, dvTrans.f1.data(), dvTrans.f2.data(), ntr, dimsp, nlev);

  Varray<double> fpwork(nlat * nfc * nlev);

  sp2fc(su, fpwork.data(), spTrans.poli.data(), nlev, nlat, nfc, ntr);
  scaluv(fpwork.data(), spTrans.rcoslat.data(), nlat, nfc * nlev);

  if (fcTrans.use_fftw)
    fc2gp(fpwork.data(), gu, nlat, nlon, nlev, nfc);
  else
    fc2gp(fcTrans.vtrig.data(), fcTrans.ifax, fpwork.data(), gu, nlat, nlon, nlev, nfc);

  sp2fc(sv, fpwork.data(), spTrans.poli.data(), nlev, nlat, nfc, ntr);
  scaluv(fpwork.data(), spTrans.rcoslat.data(), nlat, nfc * nlev);

  if (fcTrans.use_fftw)
    fc2gp(fpwork.data(), gv, nlat, nlon, nlev, nfc);
  else
    fc2gp(fcTrans.vtrig.data(), fcTrans.ifax, fpwork.data(), gv, nlat, nlon, nlev, nfc);
}
