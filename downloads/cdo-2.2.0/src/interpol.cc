/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <atomic>
#include <algorithm>

#include <cdi.h>

#include "process_int.h"
#include <mpim_grid.h>
#include "cdo_options.h"
#include "interpol.h"
#include "progress.h"
#include "cimdOmp.h"
#include "matrix_view.h"

/**
* Find the interval i-1 .. i in which an element x fits and return i, the
* bigger one of the interval borders or x itself if it is an interval border.
*
* If no interval can be found return the length of the array.

* @param *array ascending or descending sorted list
* @param nelem  length of the sorted list
* @param x      the element to find a position for
*/
static long
find_element(double x, long nelem, const Varray<double> &v)
{
  long ii;
  long mid = 0;
  long first = 1;
  long last = nelem;

  if (v[0] < v[nelem - 1])  // ascending order
    {
      // return the length of the array if x is out of bounds
      if (x < v[0] || x > v[nelem - 1]) return nelem;

      // search for the interval in which x fits
      // implementation: binary search algorithm
      for (ii = 1; ii < nelem; ++ii)
        {
          // binary search: divide search room in the middle
          mid = (first + last) >> 1;

          // return the bigger interval border of the interval in which x fits
          if (!(x < v[mid - 1] || x > v[mid])) break;

          // binary search: ignore half of the search room
          if (x > v[mid])
            first = mid;
          else
            last = mid;
        }
    }
  else
    {
      // return the length of the array if x is out of bounds
      if (x < v[nelem - 1] || x > v[0]) return nelem;

      // search for the interval in which x fits
      // implementation: binary search algorithm
      for (ii = 1; ii < nelem; ++ii)
        {
          // binary search: divide search room in the middle
          mid = (first + last) >> 1;

          // return the bigger interval border of the interval in which x fits
          if (!(x < v[mid] || x > v[mid - 1])) break;

          // binary search: ignore half of the search room
          if (x < v[mid])
            first = mid;
          else
            last = mid;
        }
    }

  if (mid > 1 && IS_EQUAL(x, v[mid - 1])) mid--;

  return mid;
}
/*
static
long find_element(double x, long nelem, const double *array)
{
  long ii;

  if ( array[0] < array[nelem-1] )
    {
      for ( ii = 1; ii < nelem; ii++ )
        if ( x >= array[ii-1] && x <= array[ii] ) break;
    }
  else
    {
      for ( ii = 1; ii < nelem; ii++ )
        if ( x >= array[ii] && x <= array[ii-1] ) break;
    }

  return ii;
}
*/

int
rect_grid_search(size_t &ii, size_t &jj, double x, double y, size_t nxm, size_t nym, const Varray<double> &xm,
                 const Varray<double> &ym)
{
  constexpr double rtol = 1.e-12;
  int lfound = 0;

  jj = find_element(y, nym, ym);
  if (jj >= nym && std::fabs(ym[0] - y) < rtol) jj = 1;  // fix rounding errors

  if (jj < nym)
    {
      ii = find_element(x, nxm, xm);
      if (ii >= nxm && std::fabs(xm[0] - x) < rtol) ii = 1;  // fix rounding errors

      if (ii < nxm) lfound = 1;
    }

  return lfound;
}

int
rect_grid_search2(long &imin, long &imax, double xmin, double xmax, long nxm, const Varray<double> &xm)
{
  int lfound = 0;
  imin = nxm;
  imax = -1;

  const bool lascend = (xm[0] < xm[nxm - 1]);

  long i1 = find_element(xmin, nxm, xm);
  long i2 = find_element(xmax, nxm, xm);

  if (i1 > 0 && i1 < nxm)
    {
      lfound = 1;

      if (lascend)
        {
          if (i1 > 1 && xmin <= xm[i1 - 1]) i1--;
          imin = i1 - 1;
          imax = i1 - 1;
        }
      else
        {
          if (i1 < nxm - 1 && xmin <= xm[i1]) i1++;
          imin = i1 - 1;
          imax = i1 - 1;
        }
    }

  if (i2 > 0 && i2 < nxm)
    {
      lfound = 1;

      if (lascend)
        {
          if (i2 < nxm - 1 && xmax >= xm[i2]) i2++;
          imax = i2 - 1;
          if (imin == nxm) imin = imax;
        }
      else
        {
          if (i2 > 1 && xmax >= xm[i2 - 1]) i2--;
          imin = i2 - 1;
          if (imax == -1) imax = imin;
        }
    }

  return lfound;
}

static double
intlinarr2p(long nxm, long nym, double **fieldm, const Varray<double> &xm, const Varray<double> &ym, double x, double y)
{
  long ii, jj;
  double value = 0;

  for (jj = 1; jj < nym; ++jj)
    if (y >= std::min(ym[jj - 1], ym[jj]) && y <= std::max(ym[jj - 1], ym[jj])) break;

  for (ii = 1; ii < nxm; ++ii)
    if (x >= xm[ii - 1] && x <= xm[ii]) break;

  if (jj < nym && ii < nxm)
    {
      value = fieldm[jj - 1][ii - 1] * (x - xm[ii]) * (y - ym[jj]) / ((xm[ii - 1] - xm[ii]) * (ym[jj - 1] - ym[jj]))
              + fieldm[jj - 1][ii] * (x - xm[ii - 1]) * (y - ym[jj]) / ((xm[ii] - xm[ii - 1]) * (ym[jj - 1] - ym[jj]))
              + fieldm[jj][ii - 1] * (x - xm[ii]) * (y - ym[jj - 1]) / ((xm[ii - 1] - xm[ii]) * (ym[jj] - ym[jj - 1]))
              + fieldm[jj][ii] * (x - xm[ii - 1]) * (y - ym[jj - 1]) / ((xm[ii] - xm[ii - 1]) * (ym[jj] - ym[jj - 1]));
    }

  return value;
}

template <typename T>
static void
intlinarr2(T missval, int lon_is_circular, size_t nxm, size_t nym, const Varray<T> &varray1, const Varray<double> &xm,
           const Varray<double> &ym, size_t gridsize2, Varray<T> &varray2, const Varray<double> &x, const Varray<double> &y)
{
  auto nlon1 = nxm;
  std::atomic<size_t> atomicCount{ 0 };

  if (lon_is_circular) nlon1--;
  const size_t gridsize1 = nlon1 * nym;

  std::vector<char> grid1_mask(gridsize1);
  for (size_t jj = 0; jj < nym; ++jj)
    for (size_t ii = 0; ii < nlon1; ++ii) grid1_mask[jj * nlon1 + ii] = !DBL_IS_EQUAL(varray1[jj * nlon1 + ii], missval);

  progress::init();

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
  for (size_t i = 0; i < gridsize2; ++i)
    {
      size_t src_add[4];  // address for the four source points

      varray2[i] = missval;

      atomicCount++;
      if (cdo_omp_get_thread_num() == 0) progress::update(0, 1, (double) atomicCount / gridsize2);

      size_t ii, jj;
      auto lfound = rect_grid_search(ii, jj, x[i], y[i], nxm, nym, xm, ym);

      if (lfound)
        {
          size_t iix = ii;
          if (lon_is_circular && iix == (nxm - 1)) iix = 0;
          src_add[0] = (jj - 1) * nlon1 + (ii - 1);
          src_add[1] = (jj - 1) * nlon1 + (iix);
          src_add[2] = (jj) *nlon1 + (ii - 1);
          src_add[3] = (jj) *nlon1 + (iix);

          // Check to see if points are missing values
          for (int n = 0; n < 4; ++n)
            if (!grid1_mask[src_add[n]]) lfound = 0;
        }

      if (lfound)
        {
          double wgts[4];
          wgts[0] = (x[i] - xm[ii]) * (y[i] - ym[jj]) / ((xm[ii - 1] - xm[ii]) * (ym[jj - 1] - ym[jj]));
          wgts[1] = (x[i] - xm[ii - 1]) * (y[i] - ym[jj]) / ((xm[ii] - xm[ii - 1]) * (ym[jj - 1] - ym[jj]));
          wgts[3] = (x[i] - xm[ii - 1]) * (y[i] - ym[jj - 1]) / ((xm[ii] - xm[ii - 1]) * (ym[jj] - ym[jj - 1]));
          wgts[2] = (x[i] - xm[ii]) * (y[i] - ym[jj - 1]) / ((xm[ii - 1] - xm[ii]) * (ym[jj] - ym[jj - 1]));

          // printf("%2ld %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n",
          // tgt_add, plon, plat, wgts[0], wgts[1], wgts[2], wgts[3], iw, jw);

          double fwsum = 0.0;
          for (int n = 0; n < 4; ++n) fwsum += varray1[src_add[n]] * wgts[n];
          varray2[i] = fwsum;
        }
    }

  progress::update(0, 1, 1);
}

void
intlinarr(long nxm, double *ym, double *xm, int nx, double *y, double *x)
{
  /*
    intlinarr - lineare interpolation over 1D array

    Uwe Schulzweida  04/05/1995
  */
  for (long jj = 1; jj < nxm; ++jj)
    for (long j = 0; j < nx; ++j)
      if (x[j] >= xm[jj - 1] && x[j] <= xm[jj]) y[j] = intlin(x[j], ym[jj - 1], xm[jj - 1], ym[jj], xm[jj]);
}

void
intgridbil(Field &field1, Field &field2)
{
  auto gridID1 = field1.grid;
  auto gridID2 = field2.grid;
  if (gridID1 == -1) cdo_abort("Source grid undefined!");
  if (gridID2 == -1) cdo_abort("Target grid undefined!");

  const auto missval = field1.missval;

  auto nlon1 = gridInqXsize(gridID1);
  auto nlat1 = gridInqYsize(gridID1);

  int lon_is_circular = 0;

  bool lgeorefgrid = true;
  if (grid_is_distance_generic(gridID1) && grid_is_distance_generic(gridID2)) lgeorefgrid = false;

  if (lgeorefgrid)
    {
      if (!gridHasCoordinates(gridID1)) cdo_abort("Source grid has no coordinate values!");

      lon_is_circular = gridIsCircular(gridID1);
      if (lon_is_circular) nlon1 += 1;
    }

  Varray<double> lons1(nlon1), lats1(nlat1);
  gridInqXvals(gridID1, lons1.data());
  gridInqYvals(gridID1, lats1.data());

  if (lgeorefgrid)
    {
      if (lon_is_circular) lons1[nlon1 - 1] = 0;

      cdo_grid_to_radian(gridID1, CDI_XAXIS, nlon1, lons1.data(), "grid1 center lon");
      cdo_grid_to_radian(gridID1, CDI_YAXIS, nlat1, lats1.data(), "grid1 center lat");

      if (lon_is_circular) lons1[nlon1 - 1] = lons1[0] + 2 * M_PI;
    }

  const auto xsize2 = gridInqXsize(gridID2);
  const auto ysize2 = gridInqYsize(gridID2);

  if (lgeorefgrid)
    {
      gridID2 = generate_full_point_grid(gridID2);
      if (!gridHasCoordinates(gridID2)) cdo_abort("Target cell center coordinates missing!");
    }

  const auto gridsize2 = gridInqSize(gridID2);

  Varray<double> xvals2(gridsize2), yvals2(gridsize2);

  if (lgeorefgrid)
    {
      gridInqXvals(gridID2, xvals2.data());
      gridInqYvals(gridID2, yvals2.data());

      cdo_grid_to_radian(gridID2, CDI_XAXIS, gridsize2, xvals2.data(), "grid2 center lon");
      cdo_grid_to_radian(gridID2, CDI_YAXIS, gridsize2, yvals2.data(), "grid2 center lat");

      for (size_t i = 0; i < gridsize2; ++i)
        {
          if (xvals2[i] < lons1[0]) xvals2[i] += 2 * M_PI;
          if (xvals2[i] > lons1[nlon1 - 1]) xvals2[i] -= 2 * M_PI;
        }
    }
  else
    {
      Varray<double> xcoord(xsize2), ycoord(ysize2);
      gridInqXvals(gridID2, xcoord.data());
      gridInqYvals(gridID2, ycoord.data());

      for (size_t j = 0; j < ysize2; ++j)
        for (size_t i = 0; i < xsize2; ++i)
          {
            xvals2[j * xsize2 + i] = xcoord[i];
            yvals2[j * xsize2 + i] = ycoord[j];
          }
    }

  if (field2.grid != gridID2) gridDestroy(gridID2);

  if (field1.memType != field2.memType) cdo_abort("Interal error, memType of field1 and field2 differ!");

  if (field1.memType == MemType::Float)
    intlinarr2((float) missval, lon_is_circular, nlon1, nlat1, field1.vec_f, lons1, lats1, gridsize2, field2.vec_f, xvals2, yvals2);
  else
    intlinarr2(missval, lon_is_circular, nlon1, nlat1, field1.vec_d, lons1, lats1, gridsize2, field2.vec_d, xvals2, yvals2);

  field_num_mv(field2);
}

// source code from pingo
void
interpolate(Field &field1, Field &field2)
{
  long i;
  long ilon, ilat, olon, olat;
  long l11, l12, l21, l22, l1, l2;
  double volon1, volon2, volat1, volat2;
  double vilat1, vilat2;
  double vlon1, vlon2, vlat1, vlat2;
  long ilon1, ilon2;
  long k, n;
  int wrap_around, xlat_is_ascending;
  double a11, a12, a21, a22, b11, b12, b21, b22;
  double faclon1, faclon2, faclat1, faclat2;

  const auto gridIDi = field1.grid;
  const auto gridIDo = field2.grid;
  auto &arrayIn = field1.vec_d;
  auto &arrayOut = field2.vec_d;
  const auto missval = field1.missval;

  // gridsize_i = gridInqSize(gridIDi);
  size_t gridsize_o = gridInqSize(gridIDo);

  long nlon = gridInqXsize(gridIDi);
  long nlat = gridInqYsize(gridIDi);
  long out_nlon = gridInqXsize(gridIDo);
  long out_nlat = gridInqYsize(gridIDo);

  Varray<double> lon_array(nlon + 2);
  Varray<double> lat_array(nlat + 2);
  double *lon = lon_array.data() + 1;
  double *lat = lat_array.data() + 1;

  if (!gridHasCoordinates(gridIDi)) cdo_abort("Source grid has no values");
  if (!gridHasCoordinates(gridIDo)) cdo_abort("Target grid has no values");

  gridInqXvals(gridIDi, lon);
  gridInqYvals(gridIDi, lat);

  // Convert lat/lon units if required
  cdo_grid_to_degree(gridIDi, CDI_XAXIS, nlon, lon, "grid1 center lon");
  cdo_grid_to_degree(gridIDi, CDI_YAXIS, nlat, lat, "grid1 center lat");

  if (nlon > 1)
    {
      lon[-1] = (lon[nlon - 1] - 360 > 2 * lon[0] - lon[1]) ? lon[nlon - 1] - 360 : 2 * lon[0] - lon[1];
      lon[nlon] = (lon[0] + 360 < 2 * lon[nlon - 1] - lon[nlon - 2]) ? lon[0] + 360 : 2 * lon[nlon - 1] - lon[nlon - 2];
    }
  else
    {
      lon[-1] = lon[0] - 360;
      lon[1] = lon[0] + 360;
    }

  if (nlat > 1)
    {
      lat[-1] = 2 * lat[0] - lat[1];
      lat[nlat] = 2 * lat[nlat - 1] - lat[nlat - 2];
    }
  else
    {
      lat[-1] = lat[0] - 10;
      lat[1] = lat[nlat - 1] + 10;
    }

  if (lat[-1] < -90) lat[-1] = -99;
  if (lat[-1] > 90) lat[-1] = 99;
  if (lat[nlat] < -90) lat[nlat] = -99;
  if (lat[nlat] > 90) lat[nlat] = 99;

  Varray<double> lono_array((out_nlon < 2) ? 4 : out_nlon + 2);
  Varray<double> lato_array((out_nlat < 2) ? 4 : out_nlat + 2);
  double *lono = lono_array.data() + 1;
  double *lato = lato_array.data() + 1;

  gridInqXvals(gridIDo, lono);
  gridInqYvals(gridIDo, lato);

  // Convert lat/lon units if required
  cdo_grid_to_degree(gridIDo, CDI_XAXIS, out_nlon, lono, "grid2 center lon");
  cdo_grid_to_degree(gridIDo, CDI_YAXIS, out_nlat, lato, "grid2 center lat");

  for (i = 0; i < out_nlon - 1; ++i)
    if (lono[i + 1] <= lono[i]) break;

  for (i++; i < out_nlon; ++i)
    {
      lono[i] += 360;
      if (i < out_nlon - 1 && lono[i + 1] + 360 <= lono[i]) cdo_abort("Longitudes of output grid are not in ascending order!");
    }

  if (lono[out_nlon - 1] - lono[0] >= 360) cdo_abort("The area covered by the longitudes of output grid must not overlap!");

  if (lato[0] > 90.001 || lato[out_nlat - 1] > 90.001 || lato[0] < -90.001 || lato[out_nlat - 1] < -90.001)
    {
      cdo_abort("Latitudes of output grid must be between 90 and -90!");
    }

  for (i = 0; i < out_nlat - 1; ++i)
    if (IS_EQUAL(lato[i + 1], lato[i]) || (i < out_nlat - 2 && ((lato[i + 1] > lato[i]) != (lato[i + 2] > lato[i + 1]))))
      {
        cdo_abort("Latitudes of output grid must be in descending or ascending order!");
      }

  if (out_nlon > 1)
    {
      lono[-1] = (lono[out_nlon - 1] - 360 > 2 * lono[0] - lono[1]) ? lono[out_nlon - 1] - 360 : 2 * lono[0] - lono[1];
      lono[out_nlon] = (lono[0] + 360 < 2 * lono[out_nlon - 1] - lono[out_nlon - 2]) ? lono[0] + 360
                                                                                     : 2 * lono[out_nlon - 1] - lono[out_nlon - 2];
    }
  else
    {
      lono[-1] = lono[0] - 360;
      lono[1] = lono[0] + 360;
    }

  if (out_nlat > 1)
    {
      lato[-1] = 2 * lato[0] - lato[1];
      lato[out_nlat] = 2 * lato[out_nlat - 1] - lato[out_nlat - 2];
    }
  else
    {
      lato[-1] = lato[0] - 10;
      lato[1] = lato[out_nlat - 1] + 10;
    }

  if (lato[-1] < -90) lato[-1] = -99;
  if (lato[-1] > 90) lato[-1] = 99;
  if (lato[out_nlat] < -90) lato[out_nlat] = -99;
  if (lato[out_nlat] > 90) lato[out_nlat] = 99;

  const auto nxlon = 2 * nlon + 1;
  const auto nxlat = 2 * nlat + 1;
  Varray<double> xin_array(nxlon * nxlat);
  MatrixView<double> xin(xin_array.data(), nxlat, nxlon);

  Varray<double> xlon(nxlon);
  for (ilon = 0; ilon < nlon; ilon++)
    {
      xlon[2 * ilon + 1] = lon[ilon];
      xlon[2 * ilon] = (lon[ilon - 1] + lon[ilon]) / 2;
    }
  xlon[2 * nlon] = (lon[nlon - 1] + lon[nlon]) / 2;

  Varray<double> xlat(2 * nlat + 1);
  for (ilat = 0; ilat < nlat; ilat++)
    {
      xlat[2 * ilat + 1] = lat[ilat];
      xlat[2 * ilat] = (lat[ilat - 1] + lat[ilat]) / 2;
    }
  xlat[2 * nlat] = (lat[nlat - 1] + lat[nlat]) / 2;

  MatrixView<double> in0(arrayIn.data(), nlat, nlon);

  Varray<long> ilon11(out_nlon), ilon12(out_nlon), ilon21(out_nlon), ilon22(out_nlon);
  Varray<double> volon11(out_nlon), volon12(out_nlon), volon21(out_nlon), volon22(out_nlon);

  for (olon = 0; olon < out_nlon; olon++)
    {
      volon1 = (lono[olon - 1] + lono[olon]) / 2;
      volon2 = (lono[olon] + lono[olon + 1]) / 2;
      if (IS_EQUAL(volon1, volon2)) volon2 += 360;
      volon2 -= 360 * std::floor((volon1 - xlon[0]) / 360);
      volon1 -= 360 * std::floor((volon1 - xlon[0]) / 360);
      volon21[olon] = volon1;
      volon22[olon] = volon2;
      for (l21 = 0; l21 < nxlon && xlon[l21] < volon1; l21++)
        ;
      for (l22 = l21; l22 < nxlon && xlon[l22] < volon2; l22++)
        ;
      volon1 -= 360;
      volon2 -= 360;
      volon11[olon] = volon1;
      volon12[olon] = volon2;
      for (l11 = 0; xlon[l11] < volon1; l11++)
        ;
      for (l12 = l11; l12 < nxlon && xlon[l12] < volon2; l12++)
        ;
      ilon11[olon] = l11;
      ilon12[olon] = l12;
      ilon21[olon] = l21;
      ilon22[olon] = l22;
    }

  Varray<long> ilat1(out_nlat), ilat2(out_nlat);

  xlat_is_ascending = xlat[0] <= xlat[nxlat - 1];
  for (olat = 0; olat < out_nlat; olat++)
    {
      volat1 = (lato[olat - 1] + lato[olat]) / 2;
      volat2 = (lato[olat] + lato[olat + 1]) / 2;
      if (!xlat_is_ascending)
        {
          if (volat1 > volat2)
            {
              for (l1 = 0; l1 < nxlat && xlat[l1] > volat1; ++l1)
                ;
              for (l2 = l1; l2 < nxlat && xlat[l2] > volat2; ++l2)
                ;
            }
          else
            {
              for (l1 = 0; l1 < nxlat && xlat[l1] > volat2; ++l1)
                ;
              for (l2 = l1; l2 < nxlat && xlat[l2] > volat1; ++l2)
                ;
            }
        }
      else
        {
          if (volat1 < volat2)
            {
              for (l1 = 0; l1 < nxlat && xlat[l1] < volat1; ++l1)
                ;
              for (l2 = l1; l2 < nxlat && xlat[l2] < volat2; ++l2)
                ;
            }
          else
            {
              for (l1 = 0; l1 < nxlat && xlat[l1] < volat2; ++l1)
                ;
              for (l2 = l1; l2 < nxlat && xlat[l2] < volat1; ++l2)
                ;
            }
        }

      ilat1[olat] = l1;
      ilat2[olat] = l2;
    }

  MatrixView<double> xout(arrayOut.data(), out_nlat, out_nlon);

  wrap_around = nlon > 1 && (lon[nlon - 1] >= lon[-1] + 360 - 0.001 || lon[nlon] >= lon[0] + 360 - 0.001);

  for (ilat = 0; ilat < nlat; ilat++)
    for (ilon = 0; ilon < nlon; ilon++) xin[2 * ilat + 1][2 * ilon + 1] = in0[ilat][ilon];

  for (ilat = 0; ilat < nxlat; ilat += 2)
    for (ilon = 1; ilon < nxlon; ilon += 2)
      {
        double sum = 0.0;
        n = 0;
        if (ilat > 0 && !DBL_IS_EQUAL(xin[ilat - 1][ilon], missval))
          {
            sum += xin[ilat - 1][ilon];
            n++;
          }
        if (ilat < nxlat - 1 && !DBL_IS_EQUAL(xin[ilat + 1][ilon], missval))
          {
            sum += xin[ilat + 1][ilon];
            n++;
          }
        xin[ilat][ilon] = n ? sum / n : missval;
      }

  for (ilat = 1; ilat < nxlat; ilat += 2)
    for (ilon = 0; ilon < nxlon; ilon += 2)
      {
        double sum = 0.0;
        n = 0;
        if (ilon > 0 && !DBL_IS_EQUAL(xin[ilat][ilon - 1], missval))
          {
            sum += xin[ilat][ilon - 1];
            n++;
          }
        if (ilon == 0 && wrap_around && !DBL_IS_EQUAL(xin[ilat][2 * nlon - 1], missval))
          {
            sum += xin[ilat][2 * nlon - 1];
            n++;
          }
        if (ilon < nxlon - 1 && !DBL_IS_EQUAL(xin[ilat][ilon + 1], missval))
          {
            sum += xin[ilat][ilon + 1];
            n++;
          }
        if (ilon == nxlon - 1 && wrap_around && !DBL_IS_EQUAL(xin[ilat][1], missval))
          {
            sum += xin[ilat][1];
            n++;
          }
        xin[ilat][ilon] = n ? sum / n : missval;
      }

  for (ilat = 0; ilat < nxlat; ilat += 2)
    for (ilon = 0; ilon < nxlon; ilon += 2)
      {
        double sum = 0.0;
        n = 0;
        if (ilon > 0 && !DBL_IS_EQUAL(xin[ilat][ilon - 1], missval))
          {
            sum += xin[ilat][ilon - 1];
            n++;
          }
        if (ilon == 0 && wrap_around && !DBL_IS_EQUAL(xin[ilat][2 * nlon - 1], missval))
          {
            sum += xin[ilat][2 * nlon - 1];
            n++;
          }
        if (ilon < nxlon - 1 && !DBL_IS_EQUAL(xin[ilat][ilon + 1], missval))
          {
            sum += xin[ilat][ilon + 1];
            n++;
          }
        if (ilon == nxlon - 1 && wrap_around && !DBL_IS_EQUAL(xin[ilat][1], missval))
          {
            sum += xin[ilat][1];
            n++;
          }
        if (ilat > 0 && !DBL_IS_EQUAL(xin[ilat - 1][ilon], missval))
          {
            sum += xin[ilat - 1][ilon];
            n++;
          }
        if (ilat < nxlat - 1 && !DBL_IS_EQUAL(xin[ilat + 1][ilon], missval))
          {
            sum += xin[ilat + 1][ilon];
            n++;
          }
        xin[ilat][ilon] = n ? sum / n : missval;
      }

  for (olat = 0; olat < out_nlat; olat++)
    {
      if (lato[-1] < lato[out_nlat])
        {
          volat1 = (lato[olat - 1] + lato[olat]) / 2;
          volat2 = (lato[olat] + lato[olat + 1]) / 2;
        }
      else
        {
          volat2 = (lato[olat - 1] + lato[olat]) / 2;
          volat1 = (lato[olat] + lato[olat + 1]) / 2;
        }

      for (olon = 0; olon < out_nlon; olon++)
        {
          double sum = 0.0;
          double wsum = 0.0;
          for (k = 0; k < 2; ++k)
            {
              if (k == 0)
                {
                  ilon1 = ilon11[olon];
                  ilon2 = ilon12[olon];
                  volon1 = volon11[olon];
                  volon2 = volon12[olon];
                }
              else
                {
                  ilon1 = ilon21[olon];
                  ilon2 = ilon22[olon];
                  volon1 = volon21[olon];
                  volon2 = volon22[olon];
                }

              for (ilon = ilon1; ilon <= ilon2; ilon++)
                {
                  if (ilon == 0 || ilon == nxlon) continue;
                  const auto vilon1 = xlon[ilon - 1];
                  const auto vilon2 = xlon[ilon];
                  for (ilat = ilat1[olat]; ilat <= ilat2[olat]; ilat++)
                    {
                      if (ilat == 0 || ilat == nxlat) continue;
                      if (xlat_is_ascending)
                        {
                          vilat1 = xlat[ilat - 1];
                          vilat2 = xlat[ilat];
                          a11 = xin[ilat - 1][ilon - 1];
                          a12 = xin[ilat - 1][ilon];
                          a21 = xin[ilat][ilon - 1];
                          a22 = xin[ilat][ilon];
                        }
                      else
                        {
                          vilat1 = xlat[ilat];
                          vilat2 = xlat[ilat - 1];
                          a11 = xin[ilat][ilon - 1];
                          a12 = xin[ilat][ilon];
                          a21 = xin[ilat - 1][ilon - 1];
                          a22 = xin[ilat - 1][ilon];
                        }
                      if (DBL_IS_EQUAL(a11, missval) || DBL_IS_EQUAL(a12, missval) || DBL_IS_EQUAL(a21, missval)
                          || DBL_IS_EQUAL(a22, missval))
                        {
                          continue;
                        }
                      if (volon1 <= vilon1 && vilon2 <= volon2 && volat1 <= vilat1 && vilat2 <= volat2)
                        {
                          vlon1 = vilon1 * M_PI / 180;
                          vlon2 = vilon2 * M_PI / 180;
                          vlat1 = vilat1 * M_PI / 180;
                          vlat2 = vilat2 * M_PI / 180;
                          b11 = a11;
                          b12 = a12;
                          b21 = a21;
                          b22 = a22;
                        }
                      else
                        {
                          vlon1 = ((volon1 <= vilon1) ? vilon1 : volon1);
                          vlon2 = ((vilon2 <= volon2) ? vilon2 : volon2);
                          vlat1 = ((volat1 <= vilat1) ? vilat1 : volat1);
                          vlat2 = ((vilat2 <= volat2) ? vilat2 : volat2);
                          if (vlon1 >= vlon2 - (volon2 - volon1) * 1e-5 || vlat1 >= vlat2 - (volat2 - volat1) * 1e-5) { continue; }
                          faclon1 = (vlon1 - vilon1) / (vilon2 - vilon1);
                          faclon2 = (vlon2 - vilon1) / (vilon2 - vilon1);
                          faclat1 = (vlat1 - vilat1) / (vilat2 - vilat1);
                          faclat2 = (vlat2 - vilat1) / (vilat2 - vilat1);
                          vlon1 *= M_PI / 180;
                          vlon2 *= M_PI / 180;
                          vlat1 *= M_PI / 180;
                          vlat2 *= M_PI / 180;
                          b11 = a11 + (a12 - a11) * faclon1 + (a21 - a11) * faclat1 + (a22 - a12 - a21 + a11) * faclon1 * faclat1;
                          b12 = a11 + (a12 - a11) * faclon2 + (a21 - a11) * faclat1 + (a22 - a12 - a21 + a11) * faclon2 * faclat1;
                          b21 = a11 + (a12 - a11) * faclon1 + (a21 - a11) * faclat2 + (a22 - a12 - a21 + a11) * faclon1 * faclat2;
                          b22 = a11 + (a12 - a11) * faclon2 + (a21 - a11) * faclat2 + (a22 - a12 - a21 + a11) * faclon2 * faclat2;
                        }
                      wsum += (vlon2 - vlon1) * (std::sin(vlat2) - std::sin(vlat1));
                      const auto t = 2.0 * std::sin((vlat2 + vlat1) / 2.0) * std::sin((vlat2 - vlat1) / 2.0) / (vlat2 - vlat1);
                      sum += (vlon2 - vlon1) / 2 * ((b11 + b12) * (t - std::sin(vlat1)) + (b21 + b22) * (std::sin(vlat2) - t));
                    }
                }
            }
          xout[olat][olon] = IS_NOT_EQUAL(wsum, 0) ? sum / wsum : missval;
        }
    }

  field2.nmiss = varray_num_mv(gridsize_o, arrayOut, missval);
}
