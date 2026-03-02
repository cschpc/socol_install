/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include "process_int.h"
#include "cdo_wtime.h"
#include "remap.h"
#include "remap_store_link_cnsrv.h"
#include "cdo_options.h"
#include "progress.h"
#include "cimdOmp.h"

#define ZERO 0.0
#define ONE 1.0
#define TWO 2.0
#define THREE 3.0
#define HALF 0.5
#define QUART 0.25

#define BABY_STEP 0.001 /* original value */

/* static double north_thresh =  1.45;  */ /* threshold for coord transformation */
/* static double south_thresh = -2.00;  */ /* threshold for coord transformation */
static double north_thresh = 2.00;         /* threshold for coord transformation */
static double south_thresh = -2.00;        /* threshold for coord transformation */

/* threshold for coord transformation */
void
remap_set_threshhold(double threshhold)
{
  north_thresh = threshhold;
  south_thresh = -threshhold;
}

/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
/*                                                                         */
/*      CONSERVATIVE INTERPOLATION                                         */
/*                                                                         */
/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */

/*
    This routine is identical to the intersection routine except
    that a coordinate transformation (using a Lambert azimuthal
    equivalent projection) is performed to treat polar cells more
    accurately.
*/
static void
pole_intersection(long *location, double *intrsct_lat, double *intrsct_lon, bool *lcoinc, bool *lthresh, double beglat,
                  double beglon, double endlat, double endlon, double *begseg, bool lrevers, long numSearchCells, long srch_corners,
                  const size_t *srch_add, const double *srch_corner_lat, const double *srch_corner_lon, bool *luse_last,
                  double *intrsct_x, double *intrsct_y, int *avoid_pole_count, double *avoid_pole_offset)
{
  /*
    Intent(in):
    double beglat, beglon,  ! beginning lat/lon endpoints for segment
    double endlat, endlon   ! ending    lat/lon endpoints for segment
    bool   lrevers          ! flag true if segment integrated in reverse

    Intent(inout) :
    double begseg[2] ! begin lat/lon of full segment
    int *location    ! address in destination array containing this
                     ! segment -- also may contain last location on entry
    bool *lthresh    ! flag segment crossing threshold boundary

    intent(out):
    *int lcoinc      ! flag segment coincident with grid line
    double *intrsct_lat, *intrsct_lon ! lat/lon coords of next intersect.
  */
  /* Local variables */
  long n, next_n, cell;
  long ioffset;
  int loutside; /* flags points outside grid */

  double pi4, rns;                           /*  north/south conversion */
  double x1, x2;                             /*  local x variables for segment */
  double y1, y2;                             /*  local y variables for segment */
  double begx, begy;                         /*  beginning x,y variables for segment */
  double endx, endy;                         /*  beginning x,y variables for segment */
  double begsegx, begsegy;                   /*  beginning x,y variables for segment */
  double grdx1, grdx2;                       /*  local x variables for grid cell */
  double grdy1, grdy2;                       /*  local y variables for grid cell */
  double vec1_y, vec1_x;                     /*  vectors and cross products used */
  double vec2_y, vec2_x;                     /*  during grid search */
  double cross_product, eps;                 /*  eps=small offset away from intersect */
  double s1, s2, determ;                     /*  variables used for linear solve to   */
  double mat1, mat2, mat3, mat4, rhs1, rhs2; /* find intersection */

  /*printf("pole_intersection: %g %g %g %g\n", beglat, beglon, endlat,
   * endlon);*/

  /* Initialize defaults, flags, etc. */

  if (!*lthresh) *location = -1;
  *lcoinc = false;
  *intrsct_lat = endlat;
  *intrsct_lon = endlon;

  loutside = false;
  s1 = ZERO;

  /* Convert coordinates */

  std::vector<double> srch_corner_x(srch_corners * numSearchCells);
  std::vector<double> srch_corner_y(srch_corners * numSearchCells);

  if (beglat > ZERO)
    {
      pi4 = QUART * PI;
      rns = ONE;
    }
  else
    {
      pi4 = -QUART * PI;
      rns = -ONE;
    }

  if (*luse_last)
    {
      x1 = *intrsct_x;
      y1 = *intrsct_y;
    }
  else
    {
      x1 = rns * TWO * std::sin(pi4 - HALF * beglat) * std::cos(beglon);
      y1 = TWO * std::sin(pi4 - HALF * beglat) * std::sin(beglon);
      *luse_last = true;
    }

  x2 = rns * TWO * std::sin(pi4 - HALF * endlat) * std::cos(endlon);
  y2 = TWO * std::sin(pi4 - HALF * endlat) * std::sin(endlon);

  for (n = 0; n < srch_corners * numSearchCells; ++n)
    {
      srch_corner_x[n] = rns * TWO * std::sin(pi4 - HALF * srch_corner_lat[n]) * std::cos(srch_corner_lon[n]);
      srch_corner_y[n] = TWO * std::sin(pi4 - HALF * srch_corner_lat[n]) * std::sin(srch_corner_lon[n]);
    }

  begx = x1;
  begy = y1;
  endx = x2;
  endy = y2;
  begsegx = rns * TWO * std::sin(pi4 - HALF * begseg[0]) * std::cos(begseg[1]);
  begsegy = TWO * std::sin(pi4 - HALF * begseg[0]) * std::sin(begseg[1]);
  *intrsct_x = endx;
  *intrsct_y = endy;

  /*
     Search for location of this segment in ocean grid using cross
     product method to determine whether a point is enclosed by a cell
  */
  while (true) /* srch_loop */
    {
      /* If last segment crossed threshold, use that location */

      if (*lthresh)
        {
          for (cell = 0; cell < numSearchCells; ++cell)
            if (srch_add[cell] == (size_t) *location)
              {
                eps = TINY;
                goto after_srch_loop;
              }
        }

      /* Otherwise normal search algorithm */

      for (cell = 0; cell < numSearchCells; ++cell) /* cell_loop  */
        {
          ioffset = cell * srch_corners;
          for (n = 0; n < srch_corners; ++n) /* corner_loop */
            {
              next_n = (n + 1) % srch_corners;
              /*
                Here we take the cross product of the vector making
                up each cell side with the vector formed by the vertex
                and search point.  If all the cross products are
                positive, the point is contained in the cell.
              */
              vec1_x = srch_corner_x[ioffset + next_n] - srch_corner_x[ioffset + n];
              vec1_y = srch_corner_y[ioffset + next_n] - srch_corner_y[ioffset + n];
              vec2_x = x1 - srch_corner_x[ioffset + n];
              vec2_y = y1 - srch_corner_y[ioffset + n];

              /* If endpoint coincident with vertex, offset the endpoint */

              if (IS_EQUAL(vec2_x, 0) && IS_EQUAL(vec2_y, 0))
                {
                  x1 += 1.e-10 * (x2 - x1);
                  y1 += 1.e-10 * (y2 - y1);
                  vec2_x = x1 - srch_corner_x[ioffset + n];
                  vec2_y = y1 - srch_corner_y[ioffset + n];
                }

              cross_product = vec1_x * vec2_y - vec2_x * vec1_y;

              /*
                If the cross product for a side is ZERO, the point
                  lies exactly on the side or the length of a side
                  is ZERO.  If the length is ZERO set det > 0.
                  otherwise, perform another cross
                  product between the side and the segment itself.
                If this cross product is also ZERO, the line is
                  coincident with the cell boundary - perform the
                  dot product and only choose the cell if the dot
                  product is positive (parallel vs anti-parallel).
              */
              if (IS_EQUAL(cross_product, 0))
                {
                  if (IS_NOT_EQUAL(vec1_x, 0) || IS_NOT_EQUAL(vec1_y, 0))
                    {
                      vec2_x = x2 - x1;
                      vec2_y = y2 - y1;
                      cross_product = vec1_x * vec2_y - vec2_x * vec1_y;
                    }
                  else
                    cross_product = ONE;

                  if (IS_EQUAL(cross_product, 0))
                    {
                      *lcoinc = true;
                      cross_product = vec1_x * vec2_x + vec1_y * vec2_y;
                      if (lrevers) cross_product = -cross_product;
                    }
                }

              /* If cross product is less than ZERO, this cell doesn't work */

              if (cross_product < ZERO) break; /* corner_loop */

            } /* corner_loop */

          /* If cross products all positive, we found the location */

          if (n >= srch_corners)
            {
              *location = srch_add[cell];
              /*
                If the beginning of this segment was outside the
                grid, invert the segment so the intersection found
                will be the first intersection with the grid
              */
              if (loutside)
                {
                  x2 = begx;
                  y2 = begy;
                  *location = -1;
                  eps = -TINY;
                }
              else
                eps = TINY;

              goto after_srch_loop;
            }

          /* Otherwise move on to next cell */

        } /* cell_loop */

      /*
        If no cell found, the point lies outside the grid.
        take some baby steps along the segment to see if any
        part of the segment lies inside the grid.
      */
      loutside = true;
      s1 = s1 + BABY_STEP;
      x1 = begx + s1 * (x2 - begx);
      y1 = begy + s1 * (y2 - begy);

      /* Reached the end of the segment and still outside the grid return no
       * intersection */

      if (s1 >= ONE)
        {
          *luse_last = false;
          return;
        }
    } /* srch_loop */

after_srch_loop:

  /*
    Now that a cell is found, search for the next intersection.
    Loop over sides of the cell to find intersection with side
    must check all sides for coincidences or intersections
  */

  ioffset = cell * srch_corners;

  for (n = 0; n < srch_corners; ++n) /* intrsct_loop */
    {
      next_n = (n + 1) % srch_corners;

      grdy1 = srch_corner_y[ioffset + n];
      grdy2 = srch_corner_y[ioffset + next_n];
      grdx1 = srch_corner_x[ioffset + n];
      grdx2 = srch_corner_x[ioffset + next_n];

      /* Set up linear system to solve for intersection */

      mat1 = x2 - x1;
      mat2 = grdx1 - grdx2;
      mat3 = y2 - y1;
      mat4 = grdy1 - grdy2;
      rhs1 = grdx1 - x1;
      rhs2 = grdy1 - y1;

      determ = mat1 * mat4 - mat2 * mat3;

      /*
         If the determinant is ZERO, the segments are either
           parallel or coincident.  Coincidences were detected
           above so do nothing.
         If the determinant is non-ZERO, solve for the linear
           parameters s for the intersection point on each line
           segment.
         If 0<s1,s2<1 then the segment intersects with this side.
           Return the point of intersection (adding a small
           number so the intersection is off the grid line).
      */
      if (std::fabs(determ) > 1.e-30)
        {
          s1 = (rhs1 * mat4 - mat2 * rhs2) / determ;
          s2 = (mat1 * rhs2 - rhs1 * mat3) / determ;

          /* Uwe Schulzweida: s1 >= ZERO! (bug fix) */
          if (s2 >= ZERO && s2 <= ONE && s1 >= ZERO && s1 <= ONE)
            {
              /*
                Recompute intersection based on full segment
                so intersections are consistent for both sweeps
              */
              if (!loutside)
                {
                  mat1 = x2 - begsegx;
                  mat3 = y2 - begsegy;
                  rhs1 = grdx1 - begsegx;
                  rhs2 = grdy1 - begsegy;
                }
              else
                {
                  mat1 = x2 - endx;
                  mat3 = y2 - endy;
                  rhs1 = grdx1 - endx;
                  rhs2 = grdy1 - endy;
                }

              determ = mat1 * mat4 - mat2 * mat3;

              /*
                Sometimes due to roundoff, the previous
                determinant is non-ZERO, but the lines
                are actually coincident.  If this is the
                case, skip the rest.
              */
              if (IS_NOT_EQUAL(determ, 0))
                {
                  s1 = (rhs1 * mat4 - mat2 * rhs2) / determ;
                  s2 = (mat1 * rhs2 - rhs1 * mat3) / determ;

                  if (!loutside)
                    {
                      *intrsct_x = begsegx + s1 * mat1;
                      *intrsct_y = begsegy + s1 * mat3;
                    }
                  else
                    {
                      *intrsct_x = endx + s1 * mat1;
                      *intrsct_y = endy + s1 * mat3;
                    }

                  /* Convert back to lat/lon coordinates */

                  *intrsct_lon = rns * std::atan2(*intrsct_y, *intrsct_x);
                  if (*intrsct_lon < ZERO) *intrsct_lon = *intrsct_lon + PI2;

                  if (std::fabs(*intrsct_x) > 1.e-10)
                    *intrsct_lat = (pi4 - std::asin(rns * HALF * (*intrsct_x) / std::cos(*intrsct_lon))) * TWO;
                  else if (std::fabs(*intrsct_y) > 1.e-10)
                    *intrsct_lat = (pi4 - std::asin(HALF * (*intrsct_y) / std::sin(*intrsct_lon))) * TWO;
                  else
                    *intrsct_lat = TWO * pi4;

                  /* Add offset in transformed space for next pass. */

                  if (s1 - eps / determ < ONE)
                    {
                      *intrsct_x = *intrsct_x - mat1 * (eps / determ);
                      *intrsct_y = *intrsct_y - mat3 * (eps / determ);
                    }
                  else
                    {
                      if (!loutside)
                        {
                          *intrsct_x = endx;
                          *intrsct_y = endy;
                          *intrsct_lat = endlat;
                          *intrsct_lon = endlon;
                        }
                      else
                        {
                          *intrsct_x = begsegx;
                          *intrsct_y = begsegy;
                          *intrsct_lat = begseg[0];
                          *intrsct_lon = begseg[1];
                        }
                    }

                  break; /* intrsct_loop */
                }
            }
        }

      /* No intersection this side, move on to next side */

    } /* intrsct_loop */

  /*
     If segment manages to cross over pole, shift the beginning
     endpoint in order to avoid hitting pole directly
     (it is ok for endpoint to be pole point)
  */

  if (std::fabs(*intrsct_x) < 1.e-10 && std::fabs(*intrsct_y) < 1.e-10 && (IS_NOT_EQUAL(endx, 0) && IS_NOT_EQUAL(endy, 0)))
    {
      if (*avoid_pole_count > 2)
        {
          *avoid_pole_count = 0;
          *avoid_pole_offset = 10. * (*avoid_pole_offset);
        }

      cross_product = begsegx * (endy - begsegy) - begsegy * (endx - begsegx);
      *intrsct_lat = begseg[0];
      if (cross_product * (*intrsct_lat) > ZERO)
        {
          *intrsct_lon = beglon + *avoid_pole_offset;
          begseg[1] = begseg[1] + *avoid_pole_offset;
        }
      else
        {
          *intrsct_lon = beglon - *avoid_pole_offset;
          begseg[1] = begseg[1] - *avoid_pole_offset;
        }

      *avoid_pole_count = *avoid_pole_count + 1;
      *luse_last = false;
    }
  else
    {
      *avoid_pole_count = 0;
      *avoid_pole_offset = TINY;
    }

  /*
     If the segment crosses a pole threshold, reset the intersection
     to be the threshold latitude and do not reuse x,y intersect
     on next entry.  Only check if did not cross threshold last
     time - sometimes the coordinate transformation can place a
     segment on the other side of the threshold again
  */
  if (*lthresh)
    {
      if (*intrsct_lat > north_thresh || *intrsct_lat < south_thresh) *lthresh = false;
    }
  else if (beglat > ZERO && *intrsct_lat < north_thresh)
    {
      mat4 = endlat - begseg[0];
      mat3 = endlon - begseg[1];
      if (mat3 > PI) mat3 = mat3 - PI2;
      if (mat3 < -PI) mat3 = mat3 + PI2;
      *intrsct_lat = north_thresh - TINY;
      s1 = (north_thresh - begseg[0]) / mat4;
      *intrsct_lon = begseg[1] + s1 * mat3;
      *luse_last = false;
      *lthresh = true;
    }
  else if (beglat<ZERO && * intrsct_lat> south_thresh)
    {
      mat4 = endlat - begseg[0];
      mat3 = endlon - begseg[1];
      if (mat3 > PI) mat3 = mat3 - PI2;
      if (mat3 < -PI) mat3 = mat3 + PI2;
      *intrsct_lat = south_thresh + TINY;
      s1 = (south_thresh - begseg[0]) / mat4;
      *intrsct_lon = begseg[1] + s1 * mat3;
      *luse_last = false;
      *lthresh = true;
    }

  /* If reached end of segment, do not use x,y intersect on next entry */

  if (IS_EQUAL(*intrsct_lat, endlat) && IS_EQUAL(*intrsct_lon, endlon)) *luse_last = false;

} /* pole_intersection */

/*
   This routine finds the next intersection of a destination grid line with
   the line segment given by beglon, endlon, etc.
   A coincidence flag is returned if the segment is entirely coincident with
   an ocean grid line.  The cells in which to search for an intersection must
   have already been restricted in the calling routine.
*/
static void
intersection(long *location, double *intrsct_lat, double *intrsct_lon, bool *lcoinc, double beglat, double beglon, double endlat,
             double endlon, double *begseg, bool lbegin, bool lrevers, long numSearchCells, long srch_corners,
             const size_t *srch_add, const double *srch_corner_lat, const double *srch_corner_lon, long *last_loc, bool *lthresh,
             double *intrsct_lat_off, double *intrsct_lon_off, bool *luse_last, double *intrsct_x, double *intrsct_y,
             int *avoid_pole_count, double *avoid_pole_offset)
{
  /*
    Intent(in):
    bool lbegin,            ! flag for first integration along this segment
    bool lrevers            ! flag whether segment integrated in reverse
    double beglat, beglon,  ! beginning lat/lon endpoints for segment
    double endlat, endlon   ! ending    lat/lon endpoints for segment

    Intent(inout) ::
    double *begseg          ! begin lat/lon of full segment

    intent(out):
    int *location           ! address in destination array containing this
    segment bool *lcoinc            ! flag segments which are entirely
    coincident with a grid line double *intrsct_lat, *intrsct_lon ! lat/lon
    coords of next intersect.
  */
  /* Local variables */
  long n, next_n, cell;
  long ioffset;

  int loutside; /* flags points outside grid */

  double lon1, lon2;         /* local longitude variables for segment */
  double lat1, lat2;         /* local latitude  variables for segment */
  double grdlon1, grdlon2;   /* local longitude variables for grid cell */
  double grdlat1, grdlat2;   /* local latitude  variables for grid cell */
  double vec1_lat, vec1_lon; /* vectors and cross products used */
  double vec2_lat, vec2_lon; /* during grid search */
  double cross_product;
  double eps, offset;                                /* small offset away from intersect */
  double s1, s2, determ;                             /* variables used for linear solve to */
  double mat1 = 0, mat2, mat3 = 0, mat4, rhs1, rhs2; /* find intersection */

  /* Initialize defaults, flags, etc. */

  *location = -1;
  *lcoinc = false;
  *intrsct_lat = endlat;
  *intrsct_lon = endlon;

  if (numSearchCells == 0) return;

  if (beglat > north_thresh || beglat < south_thresh)
    {
      if (*lthresh) *location = *last_loc;
      pole_intersection(location, intrsct_lat, intrsct_lon, lcoinc, lthresh, beglat, beglon, endlat, endlon, begseg, lrevers,
                        numSearchCells, srch_corners, srch_add, srch_corner_lat, srch_corner_lon, luse_last, intrsct_x, intrsct_y,
                        avoid_pole_count, avoid_pole_offset);

      if (*lthresh)
        {
          *last_loc = *location;
          *intrsct_lat_off = *intrsct_lat;
          *intrsct_lon_off = *intrsct_lon;
        }
      return;
    }

  loutside = false;
  if (lbegin)
    {
      lat1 = beglat;
      lon1 = beglon;
    }
  else
    {
      lat1 = *intrsct_lat_off;
      lon1 = *intrsct_lon_off;
    }

  lat2 = endlat;
  lon2 = endlon;
  if ((lon2 - lon1) > THREE * PIH)
    lon2 -= PI2;
  else if ((lon2 - lon1) < -THREE * PIH)
    lon2 += PI2;

  s1 = ZERO;

  /*
     Search for location of this segment in ocean grid using cross
     product method to determine whether a point is enclosed by a cell
  */
  while (true) /* srch_loop */
    {
      /* If last segment crossed threshold, use that location */

      if (*lthresh)
        {
          for (cell = 0; cell < numSearchCells; ++cell)
            if (srch_add[cell] == (size_t) *last_loc)
              {
                *location = *last_loc;
                eps = TINY;
                goto after_srch_loop;
              }
        }

      /* Otherwise normal search algorithm */

      for (cell = 0; cell < numSearchCells; ++cell) /* cell_loop  */
        {
          ioffset = cell * srch_corners;
          for (n = 0; n < srch_corners; ++n) /* corner_loop */
            {
              next_n = (n + 1) % srch_corners;
              /*
                Here we take the cross product of the vector making
                up each cell side with the vector formed by the vertex
                and search point.  If all the cross products are
                positive, the point is contained in the cell.
              */
              vec1_lat = srch_corner_lat[ioffset + next_n] - srch_corner_lat[ioffset + n];
              vec1_lon = srch_corner_lon[ioffset + next_n] - srch_corner_lon[ioffset + n];
              vec2_lat = lat1 - srch_corner_lat[ioffset + n];
              vec2_lon = lon1 - srch_corner_lon[ioffset + n];

              /* If endpoint coincident with vertex, offset the endpoint */

              if (IS_EQUAL(vec2_lat, 0) && IS_EQUAL(vec2_lon, 0))
                {
                  lat1 += 1.e-10 * (lat2 - lat1);
                  lon1 += 1.e-10 * (lon2 - lon1);
                  vec2_lat = lat1 - srch_corner_lat[ioffset + n];
                  vec2_lon = lon1 - srch_corner_lon[ioffset + n];
                }

              /* Check for 0,2pi crossings */

              if (vec1_lon > PI)
                vec1_lon -= PI2;
              else if (vec1_lon < -PI)
                vec1_lon += PI2;

              if (vec2_lon > PI)
                vec2_lon -= PI2;
              else if (vec2_lon < -PI)
                vec2_lon += PI2;

              cross_product = vec1_lon * vec2_lat - vec2_lon * vec1_lat;

              /*
               If the cross product for a side is ZERO, the point
                 lies exactly on the side or the side is degenerate
                 (ZERO length).  If degenerate, set the cross
                 product to a positive number.  Otherwise perform
                 another cross product between the side and the
                 segment itself.
               If this cross product is also ZERO, the line is
                 coincident with the cell boundary - perform the
                 dot product and only choose the cell if the dot
                 product is positive (parallel vs anti-parallel).
              */
              if (IS_EQUAL(cross_product, 0))
                {
                  if (IS_NOT_EQUAL(vec1_lat, 0) || IS_NOT_EQUAL(vec1_lon, 0))
                    {
                      vec2_lat = lat2 - lat1;
                      vec2_lon = lon2 - lon1;

                      if (vec2_lon > PI)
                        vec2_lon -= PI2;
                      else if (vec2_lon < -PI)
                        vec2_lon += PI2;

                      cross_product = vec1_lon * vec2_lat - vec2_lon * vec1_lat;
                    }
                  else
                    cross_product = ONE;

                  if (IS_EQUAL(cross_product, 0))
                    {
                      *lcoinc = true;
                      cross_product = vec1_lon * vec2_lon + vec1_lat * vec2_lat;
                      if (lrevers) cross_product = -cross_product;
                    }
                }

              /* If cross product is less than ZERO, this cell doesn't work */

              if (cross_product < ZERO) break; /* corner_loop */

            } /* corner_loop */

          /* If cross products all positive, we found the location */

          if (n >= srch_corners)
            {
              *location = srch_add[cell];
              /*
                If the beginning of this segment was outside the
                grid, invert the segment so the intersection found
                will be the first intersection with the grid
              */
              if (loutside)
                {
                  lat2 = beglat;
                  lon2 = beglon;
                  *location = -1;
                  eps = -TINY;
                }
              else
                eps = TINY;

              goto after_srch_loop;
            }

          /* Otherwise move on to next cell */

        } /* cell_loop */

      /*
        If still no cell found, the point lies outside the grid.
        Take some baby steps along the segment to see if any
        part of the segment lies inside the grid.
      */
      loutside = true;
      s1 = s1 + BABY_STEP;
      lat1 = beglat + s1 * (endlat - beglat);
      lon1 = beglon + s1 * (lon2 - beglon);

      /* Reached the end of the segment and still outside the grid return no
       * intersection */

      if (s1 >= ONE) return;

    } /* srch_loop */

after_srch_loop:

  /*
    Now that a cell is found, search for the next intersection.
    Loop over sides of the cell to find intersection with side
    must check all sides for coincidences or intersections
  */

  ioffset = cell * srch_corners;

  for (n = 0; n < srch_corners; ++n) /* intrsct_loop */
    {
      next_n = (n + 1) % srch_corners;

      grdlon1 = srch_corner_lon[ioffset + n];
      grdlon2 = srch_corner_lon[ioffset + next_n];
      grdlat1 = srch_corner_lat[ioffset + n];
      grdlat2 = srch_corner_lat[ioffset + next_n];

      /* Set up linear system to solve for intersection */

      mat1 = lat2 - lat1;
      mat2 = grdlat1 - grdlat2;
      mat3 = lon2 - lon1;
      mat4 = grdlon1 - grdlon2;
      rhs1 = grdlat1 - lat1;
      rhs2 = grdlon1 - lon1;

      if (mat3 > PI)
        mat3 -= PI2;
      else if (mat3 < -PI)
        mat3 += PI2;

      if (mat4 > PI)
        mat4 -= PI2;
      else if (mat4 < -PI)
        mat4 += PI2;

      if (rhs2 > PI)
        rhs2 -= PI2;
      else if (rhs2 < -PI)
        rhs2 += PI2;

      determ = mat1 * mat4 - mat2 * mat3;

      /*
         If the determinant is ZERO, the segments are either
           parallel or coincident.  Coincidences were detected
           above so do nothing.
         If the determinant is non-ZERO, solve for the linear
           parameters s for the intersection point on each line
           segment.
         If 0<s1,s2<1 then the segment intersects with this side.
           Return the point of intersection (adding a small
           number so the intersection is off the grid line).
      */
      if (std::fabs(determ) > 1.e-30)
        {
          s1 = (rhs1 * mat4 - mat2 * rhs2) / determ;
          s2 = (mat1 * rhs2 - rhs1 * mat3) / determ;

          if (s2 >= ZERO && s2 <= ONE && s1 >= ZERO && s1 <= ONE)
            {
              /*
                Recompute intersection based on full segment
                so intersections are consistent for both sweeps
              */
              if (!loutside)
                {
                  mat1 = lat2 - begseg[0];
                  mat3 = lon2 - begseg[1];
                  rhs1 = grdlat1 - begseg[0];
                  rhs2 = grdlon1 - begseg[1];
                }
              else
                {
                  mat1 = begseg[0] - endlat;
                  mat3 = begseg[1] - endlon;
                  rhs1 = grdlat1 - endlat;
                  rhs2 = grdlon1 - endlon;
                }

              if (mat3 > PI)
                mat3 -= PI2;
              else if (mat3 < -PI)
                mat3 += PI2;

              if (rhs2 > PI)
                rhs2 -= PI2;
              else if (rhs2 < -PI)
                rhs2 += PI2;

              determ = mat1 * mat4 - mat2 * mat3;

              /*
                Sometimes due to roundoff, the previous
                determinant is non-ZERO, but the lines
                are actually coincident.  If this is the
                case, skip the rest.
              */
              if (IS_NOT_EQUAL(determ, 0))
                {
                  s1 = (rhs1 * mat4 - mat2 * rhs2) / determ;
                  s2 = (mat1 * rhs2 - rhs1 * mat3) / determ;

                  offset = s1 + eps / determ;
                  if (offset > ONE) offset = ONE;

                  if (!loutside)
                    {
                      *intrsct_lat = begseg[0] + mat1 * s1;
                      *intrsct_lon = begseg[1] + mat3 * s1;
                      *intrsct_lat_off = begseg[0] + mat1 * offset;
                      *intrsct_lon_off = begseg[1] + mat3 * offset;
                    }
                  else
                    {
                      *intrsct_lat = endlat + mat1 * s1;
                      *intrsct_lon = endlon + mat3 * s1;
                      *intrsct_lat_off = endlat + mat1 * offset;
                      *intrsct_lon_off = endlon + mat3 * offset;
                    }
                  break; /* intrsct_loop */
                }
            }
        }

      /* No intersection this side, move on to next side */

    } /* intrsct_loop */

  /*
     If the segment crosses a pole threshold, reset the intersection
     to be the threshold latitude.  Only check if this was not a
     threshold segment since sometimes coordinate transform can end
     up on other side of threshold again.
  */
  if (*lthresh)
    {
      if (*intrsct_lat < north_thresh || *intrsct_lat > south_thresh) *lthresh = false;
    }
  else if (lat1 > ZERO && *intrsct_lat > north_thresh)
    {
      *intrsct_lat = north_thresh + TINY;
      *intrsct_lat_off = north_thresh + eps * mat1;
      s1 = (*intrsct_lat - begseg[0]) / mat1;
      *intrsct_lon = begseg[1] + s1 * mat3;
      *intrsct_lon_off = begseg[1] + (s1 + eps) * mat3;
      *last_loc = *location;
      *lthresh = true;
    }
  else if (lat1 < ZERO && *intrsct_lat < south_thresh)
    {
      *intrsct_lat = south_thresh - TINY;
      *intrsct_lat_off = south_thresh + eps * mat1;
      s1 = (*intrsct_lat - begseg[0]) / mat1;
      *intrsct_lon = begseg[1] + s1 * mat3;
      *intrsct_lon_off = begseg[1] + (s1 + eps) * mat3;
      *last_loc = *location;
      *lthresh = true;
    }

} /* intersection */

/*
   This routine computes the line integral of the flux function
   that results in the interpolation weights.  The line is defined
   by the input lat/lon of the endpoints.
*/
static double
phi_gradient(double in_phi1, double in_phi2, double dphi, double f1, double f2, double grid_lon)
{
  double phi1 = in_phi1 - grid_lon;
  if (phi1 > PI)
    phi1 -= PI2;
  else if (phi1 < -PI)
    phi1 += PI2;

  double phi2 = in_phi2 - grid_lon;
  if (phi2 > PI)
    phi2 -= PI2;
  else if (phi2 < -PI)
    phi2 += PI2;

  double weight;
  if ((phi2 - phi1) < PI && (phi2 - phi1) > -PI)
    weight = dphi * (phi1 * f1 + phi2 * f2);
  else
    {
      const double fac = (phi1 > ZERO) ? PI : -PI;
      const double fint = f1 + (f2 - f1) * (fac - phi1) / std::fabs(dphi);
      weight = HALF * phi1 * (phi1 - fac) * f1 - HALF * phi2 * (phi2 + fac) * f2 + HALF * fac * (phi1 + phi2) * fint;
    }

  return weight;
}

static void
line_integral(double *weights, double in_phi1, double in_phi2, double theta1, double theta2, double grid1_lon, double grid2_lon)
{
  /*
    Intent(in):
    double in_phi1, in_phi2,     ! Longitude endpoints for the segment
    double theta1, theta2,       ! Latitude  endpoints for the segment
    double grid1_lon,            ! Reference coordinates for each
    double grid2_lon             ! Grid (to ensure correct 0,2pi interv.)

    Intent(out):
    double weights[6]            ! Line integral contribution to weights
  */

  /*  Weights for the general case based on a trapezoidal approx to the integrals. */

  const double sinth1 = std::sin(theta1);
  const double sinth2 = std::sin(theta2);
  const double costh1 = std::cos(theta1);
  const double costh2 = std::cos(theta2);

  double dphi = in_phi1 - in_phi2;
  if (dphi > PI)
    dphi -= PI2;
  else if (dphi < -PI)
    dphi += PI2;

  dphi = HALF * dphi;

  /*
     The first weight is the area overlap integral. The second and
     fourth are second-order latitude gradient weights.
  */
  weights[0] = dphi * (sinth1 + sinth2);
  weights[1] = dphi * (costh1 + costh2 + (theta1 * sinth1 + theta2 * sinth2));
  weights[3] = weights[0];
  weights[4] = weights[1];

  /*
     The third and fifth weights are for the second-order phi gradient
     component.  Must be careful of longitude range.
  */
  const double f1 = HALF * (costh1 * sinth1 + theta1);
  const double f2 = HALF * (costh2 * sinth2 + theta2);

  weights[2] = phi_gradient(in_phi1, in_phi2, dphi, f1, f2, grid1_lon);
  weights[5] = phi_gradient(in_phi1, in_phi2, dphi, f1, f2, grid2_lon);

} /* line_integral */

static void
correct_pole(RemapGrid *srcGrid, RemapGrid *tgtGrid, RemapVars &rv, std::vector<double> &src_centroid_lat,
             std::vector<double> &src_centroid_lon, std::vector<double> &tgt_centroid_lat, std::vector<double> &tgt_centroid_lon,
             GridStore &grid_store)
{
  /*
     Correct for situations where N/S pole not explicitly included in
     grid (i.e. as a grid corner point). If pole is missing from only
     one grid, need to correct only the area and centroid of that
     grid.  If missing from both, do complete weight calculation.
  */
  long n;
  long num_wts;
  long srcGridSize;
  long tgtGridSize;
  long srcCellIndex; /* current linear address for source grid cell   */
  long tgtCellIndex; /* current linear address for target grid cell   */
  double weights[6]; /* local wgt array */

  num_wts = rv.num_wts;

  srcGridSize = srcGrid->size;
  tgtGridSize = tgtGrid->size;

  /* North Pole */
  weights[0] = PI2;
  weights[1] = PI * PI;
  weights[2] = ZERO;
  weights[3] = PI2;
  weights[4] = PI * PI;
  weights[5] = ZERO;

  srcCellIndex = srcGridSize;
  /* pole_loop1 */
  for (n = 0; n < srcGridSize; ++n)
    if (srcGrid->cell_area[n] < -THREE * PIH && srcGrid->cell_center_lat[n] > ZERO)
      {
        srcCellIndex = n;
        break;
      }

  tgtCellIndex = tgtGridSize;
  /* pole_loop2 */
  for (n = 0; n < tgtGridSize; ++n)
    if (tgtGrid->cell_area[n] < -THREE * PIH && tgtGrid->cell_center_lat[n] > ZERO)
      {
        tgtCellIndex = n;
        break;
      }

  if (srcCellIndex != srcGridSize)
    {
      srcGrid->cell_area[srcCellIndex] += weights[0];
      src_centroid_lat[srcCellIndex] += weights[1];
      src_centroid_lon[srcCellIndex] += weights[2];
    }

  if (tgtCellIndex != tgtGridSize)
    {
      tgtGrid->cell_area[tgtCellIndex] += weights[3];
      tgt_centroid_lat[tgtCellIndex] += weights[4];
      tgt_centroid_lon[tgtCellIndex] += weights[5];
    }

  if (srcCellIndex != srcGridSize && tgtCellIndex != tgtGridSize)
    {
      store_link_cnsrv(rv, srcCellIndex, tgtCellIndex, num_wts, weights, grid_store);

      srcGrid->cell_frac[srcCellIndex] += weights[0];
      tgtGrid->cell_frac[tgtCellIndex] += weights[3];
    }

  /* South Pole */
  weights[0] = PI2;
  weights[1] = -PI * PI;
  weights[2] = ZERO;
  weights[3] = PI2;
  weights[4] = -PI * PI;
  weights[5] = ZERO;

  srcCellIndex = srcGridSize;
  /* pole_loop3 */
  for (n = 0; n < srcGridSize; ++n)
    if (srcGrid->cell_area[n] < -THREE * PIH && srcGrid->cell_center_lat[n] < ZERO)
      {
        srcCellIndex = n;
        break;
      }

  tgtCellIndex = tgtGridSize;
  /* pole_loop4 */
  for (n = 0; n < tgtGridSize; ++n)
    if (tgtGrid->cell_area[n] < -THREE * PIH && tgtGrid->cell_center_lat[n] < ZERO)
      {
        tgtCellIndex = n;
        break;
      }

  if (srcCellIndex != srcGridSize)
    {
      srcGrid->cell_area[srcCellIndex] += weights[0];
      src_centroid_lat[srcCellIndex] += weights[1];
      src_centroid_lon[srcCellIndex] += weights[2];
    }

  if (tgtCellIndex != tgtGridSize)
    {
      tgtGrid->cell_area[tgtCellIndex] += weights[3];
      tgt_centroid_lat[tgtCellIndex] += weights[4];
      tgt_centroid_lon[tgtCellIndex] += weights[5];
    }

  if (srcCellIndex != srcGridSize && tgtCellIndex != tgtGridSize)
    {
      store_link_cnsrv(rv, srcCellIndex, tgtCellIndex, num_wts, weights, grid_store);

      srcGrid->cell_frac[srcCellIndex] += weights[0];
      tgtGrid->cell_frac[tgtCellIndex] += weights[3];
    }
}

static inline void
norm_weight(double norm_factor, double *weights, double src_centroid_lat, double src_centroid_lon)
{
  const auto weight0 = weights[0];
  weights[0] = weight0 * norm_factor;
  weights[1] = (weights[1] - weight0 * src_centroid_lat) * norm_factor;
  weights[2] = (weights[2] - weight0 * src_centroid_lon) * norm_factor;
  constexpr double weps = 1.e-15;
  if (fabs(weights[0]) < weps) weights[0] = 0.0;
  if (fabs(weights[1]) < weps) weights[1] = 0.0;
  if (fabs(weights[2]) < weps) weights[2] = 0.0;
}

static void
normalize_weights(RemapGrid *tgtGrid, RemapVars &rv, std::vector<double> &src_centroid_lat, std::vector<double> &src_centroid_lon)
{
  /* Include centroids in weights and normalize using destination area if requested */
  long numLinks = rv.numLinks;
  // srcCellIndex: current linear address for source grid cell
  // tgtCellIndex: current linear address for target grid cell
  double *weights = &rv.wts[0];

  if (rv.normOpt == NormOpt::DESTAREA)
    {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(numLinks, rv, weights, tgtGrid, src_centroid_lat, src_centroid_lon)
#endif
      for (long n = 0; n < numLinks; ++n)
        {
          const auto src_add = rv.srcCellIndices[n];
          const auto tgt_add = rv.tgtCellIndices[n];
          const auto tgt_cell_area = tgtGrid->cell_area[tgt_add];
          const auto norm_factor = IS_NOT_EQUAL(tgt_cell_area, 0) ? ONE / tgt_cell_area : ZERO;
          norm_weight(norm_factor, &weights[n * 3], src_centroid_lat[src_add], src_centroid_lon[src_add]);
        }
    }
  else if (rv.normOpt == NormOpt::FRACAREA)
    {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(numLinks, rv, weights, tgtGrid, src_centroid_lat, src_centroid_lon)
#endif
      for (long n = 0; n < numLinks; ++n)
        {
          const auto src_add = rv.srcCellIndices[n];
          const auto tgt_add = rv.tgtCellIndices[n];
          const auto tgt_cell_frac = tgtGrid->cell_frac[tgt_add];
          const auto norm_factor = IS_NOT_EQUAL(tgt_cell_frac, 0) ? ONE / tgt_cell_frac : ZERO;
          // printf("i: %ld %g  %g %g  %g %g %g\n", n+1, norm_factor, src_centroid_lon[src_add], src_centroid_lat[src_add],
          // weights[n * 3], weights[n * 3+1], weights[n * 3+2]);
          norm_weight(norm_factor, &weights[n * 3], src_centroid_lat[src_add], src_centroid_lon[src_add]);
          // printf("o: %ld %g  %g %g  %g %g %g\n", n+1, norm_factor, src_centroid_lon[src_add], src_centroid_lat[src_add],
          // weights[n * 3], weights[n * 3+1], weights[n * 3+2]);
        }
    }
  else if (rv.normOpt == NormOpt::NONE)
    {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(numLinks, rv, weights, tgtGrid, src_centroid_lat, src_centroid_lon)
#endif
      for (long n = 0; n < numLinks; ++n)
        {
          const auto src_add = rv.srcCellIndices[n];
          const auto norm_factor = ONE;
          norm_weight(norm_factor, &weights[n * 3], src_centroid_lat[src_add], src_centroid_lon[src_add]);
        }
    }
}

static void
srvout(const char *fname, int nlon, int nlat, double *data)
{
  int ihead[8];

  ihead[0] = 99;
  ihead[1] = 0;
  ihead[2] = 0;
  ihead[3] = 0;
  ihead[4] = nlon;
  ihead[5] = nlat;
  ihead[6] = 0;
  ihead[7] = 0;

  FILE *fp = std::fopen(fname, "w");
  fprintf(fp, "%d %d %d %d %d %d %d %d\n", ihead[0], ihead[1], ihead[2], ihead[3], ihead[4], ihead[5], ihead[6], ihead[7]);

  for (int i = 0; i < nlon * nlat; ++i) fprintf(fp, "%g\n", data[i]);

  std::fclose(fp);
}

/*
  -----------------------------------------------------------------------

   This routine traces the perimeters of every grid cell on each
   grid checking for intersections with the other grid and computing
   line integrals for each subsegment.

  -----------------------------------------------------------------------
*/
void
remap_conserv_weights_scrip(RemapSearch &rsearch, RemapVars &rv)
{
  RemapGrid *srcGrid = rsearch.srcGrid;
  RemapGrid *tgtGrid = rsearch.tgtGrid;

  bool lcheck = true;

  long max_subseg = 100000; /* max number of subsegments per segment to prevent infinite loop */
  /* 1000 is too small!!! */

  long srch_corners; /* num of corners of srch cells           */

  double findex = 0;

  if (Options::cdoVerbose) cdo_print("Called %s()", __func__);

  progress::init();

  long num_wts = rv.num_wts;

  GridStore grid_store;
  grid_store_init(grid_store, tgtGrid->size);

  if (Options::cdoVerbose)
    {
      cdo_print("North threshhold: %g", north_thresh);
      cdo_print("South threshhold: %g", south_thresh);
    }

  double start = Options::cdoVerbose ? cdo_get_wtime() : 0;

  long srcGridSize = srcGrid->size;
  long tgtGridSize = tgtGrid->size;

  long srcNumCorners = srcGrid->num_cell_corners;
  long tgtNumCorners = tgtGrid->num_cell_corners;

  /* Initialize centroid arrays */

  std::vector<double> src_centroid_lat(srcGridSize);
  std::vector<double> src_centroid_lon(srcGridSize);
  std::vector<double> tgt_centroid_lat(tgtGridSize);
  std::vector<double> tgt_centroid_lon(tgtGridSize);

  for (long n = 0; n < srcGridSize; ++n)
    {
      src_centroid_lat[n] = 0;
      src_centroid_lon[n] = 0;
    }

  for (long n = 0; n < tgtGridSize; ++n)
    {
      tgt_centroid_lat[n] = 0;
      tgt_centroid_lon[n] = 0;
    }

  std::vector<double *> srch_corner_lat(Threading::ompNumThreads);  // lat of each corner of srch cells
  std::vector<double *> srch_corner_lon(Threading::ompNumThreads);  // lon of each corner of srch cells
  std::vector<long> maxSearchCells(Threading::ompNumThreads);       // num cells in restricted search arrays

  /* Integrate around each cell on source grid */

  for (long i = 0; i < Threading::ompNumThreads; ++i)
    {
      srch_corner_lat[i] = nullptr;
      srch_corner_lon[i] = nullptr;
      maxSearchCells[i] = 0;
    }

  Varray<Varray<size_t>> srch_add(Threading::ompNumThreads);  // global address of cells in srch arrays
  for (int i = 0; i < Threading::ompNumThreads; ++i) srch_add[i].resize(tgtGridSize);

  srch_corners = tgtNumCorners;

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(findex, rsearch, num_wts, src_centroid_lon, src_centroid_lat, grid_store, rv,      \
                                              Options::cdoVerbose, max_subseg, srch_corner_lat, srch_corner_lon, maxSearchCells, \
                                              srcNumCorners, srch_corners, srcGrid, tgtGrid, tgtGridSize, srcGridSize, srch_add)
#endif
  for (long srcCellIndex = 0; srcCellIndex < srcGridSize; ++srcCellIndex)
    {
      const auto ompthID = cdo_omp_get_thread_num();

#ifdef _OPENMP
#pragma omp atomic
#endif
      findex++;
      if (ompthID == 0) progress::update(0, 0.5, findex / srcGridSize);

      bool lcoinc;  // flag for coincident segments
      bool lthresh = false;
      bool luse_last = false;
      int avoid_pole_count = 0;                         // count attempts to avoid pole
      double avoid_pole_offset = TINY;                  // endpoint offset to avoid pole
      double intrsct_lat, intrsct_lon;                  // lat/lon of next intersect
      double intrsct_lat_off = 0, intrsct_lon_off = 0;  // lat/lon coords offset for next search
      double intrsct_x, intrsct_y;                      // x,y for intersection
      long last_loc = -1;                               // save location when crossing threshold
      long tgtCellIndex;

      // Get search cells
      long numSearchCells = get_srch_cells(srcCellIndex, rsearch.srcBins, rsearch.tgtBins,
                                           &rsearch.srcBins.cell_bound_box[srcCellIndex * 4], srch_add[ompthID]);

      if (numSearchCells == 0) continue;

      // Create search arrays
      if (numSearchCells > maxSearchCells[ompthID])
        {
          srch_corner_lat[ompthID] = (double *) realloc(srch_corner_lat[ompthID], srch_corners * numSearchCells * sizeof(double));
          srch_corner_lon[ompthID] = (double *) realloc(srch_corner_lon[ompthID], srch_corners * numSearchCells * sizeof(double));
          maxSearchCells[ompthID] = numSearchCells;
        }

      // gather1
      long ioffset;
      for (long n = 0; n < numSearchCells; ++n)
        {
          tgtCellIndex = srch_add[ompthID][n];
          ioffset = tgtCellIndex * srch_corners;

          long nsrch_corners = n * srch_corners;
          for (long k = 0; k < srch_corners; ++k)
            {
              srch_corner_lat[ompthID][nsrch_corners + k] = tgtGrid->cell_corner_lat[ioffset + k];
              srch_corner_lon[ompthID][nsrch_corners + k] = tgtGrid->cell_corner_lon[ioffset + k];
            }
        }

      // Integrate around this cell
      ioffset = srcCellIndex * srcNumCorners;
      for (long corner = 0; corner < srcNumCorners; ++corner)
        {
          long next_corn = (corner + 1) % srcNumCorners;

          /* Define endpoints of the current segment */

          double beglat = srcGrid->cell_corner_lat[ioffset + corner];
          double beglon = srcGrid->cell_corner_lon[ioffset + corner];
          double endlat = srcGrid->cell_corner_lat[ioffset + next_corn];
          double endlon = srcGrid->cell_corner_lon[ioffset + next_corn];
          bool lrevers = false;

          /* To ensure exact path taken during both sweeps, always integrate
           * segments in the same direction (SW to NE). */
          if ((endlat < beglat) || (IS_EQUAL(endlat, beglat) && endlon < beglon))
            {
              beglat = srcGrid->cell_corner_lat[ioffset + next_corn];
              beglon = srcGrid->cell_corner_lon[ioffset + next_corn];
              endlat = srcGrid->cell_corner_lat[ioffset + corner];
              endlon = srcGrid->cell_corner_lon[ioffset + corner];
              lrevers = true;
            }

          /*
            If this is a constant-longitude segment, skip the rest
            since the line integral contribution will be ZERO.
          */
          if (IS_EQUAL(endlon, beglon)) continue;

          double weights[6];  // local wgt array
          double begseg[2];   // begin lat/lon for full segment
          begseg[0] = beglat;
          begseg[1] = beglon;
          bool lbegin = true;

          long num_subseg = 0;
          /*
            Integrate along this segment, detecting intersections
            and computing the line integral for each sub-segment
          */
          while (IS_NOT_EQUAL(beglat, endlat) || IS_NOT_EQUAL(beglon, endlon))
            {
              /* Prevent infinite loops if integration gets stuck near cell or threshold boundary */
              num_subseg++;
              if (num_subseg >= max_subseg)
                cdo_abort("Integration stalled: num_subseg exceeded limit (grid1[%zu]: lon1=%g lon2=%g lat1=%g lat2=%g)!",
                          srcCellIndex, beglon, endlon, beglat, endlat);

              /* Uwe Schulzweida: skip very small regions */
              if (num_subseg % 1000 == 0)
                {
                  if (std::fabs(beglat - endlat) < 1.e-10 || std::fabs(beglon - endlon) < 1.e-10)
                    {
                      if (Options::cdoVerbose)
                        cdo_print("Skip very small region (grid1[%ld]): lon=%g dlon=%g lat=%g dlat=%g", srcCellIndex, beglon,
                                  endlon - beglon, beglat, endlat - beglat);
                      break;
                    }
                }

              /* Find next intersection of this segment with a gridline on grid 2. */

              intersection(&tgtCellIndex, &intrsct_lat, &intrsct_lon, &lcoinc, beglat, beglon, endlat, endlon, begseg, lbegin,
                           lrevers, numSearchCells, srch_corners, srch_add[ompthID].data(), srch_corner_lat[ompthID],
                           srch_corner_lon[ompthID], &last_loc, &lthresh, &intrsct_lat_off, &intrsct_lon_off, &luse_last,
                           &intrsct_x, &intrsct_y, &avoid_pole_count, &avoid_pole_offset);

              lbegin = false;

              /* Compute line integral for this subsegment. */

              if (tgtCellIndex != -1)
                line_integral(weights, beglon, intrsct_lon, beglat, intrsct_lat, srcGrid->cell_center_lon[srcCellIndex],
                              tgtGrid->cell_center_lon[tgtCellIndex]);
              else
                line_integral(weights, beglon, intrsct_lon, beglat, intrsct_lat, srcGrid->cell_center_lon[srcCellIndex],
                              srcGrid->cell_center_lon[srcCellIndex]);

              /* If integrating in reverse order, change sign of weights */

              if (lrevers)
                for (long k = 0; k < 6; ++k) weights[k] = -weights[k];

              /*
                Store the appropriate addresses and weights.
                Also add contributions to cell areas and centroids.
              */
              if (tgtCellIndex != -1)
                if (srcGrid->mask[srcCellIndex])
                  {
#ifdef _OPENMP
#pragma omp critical
#endif
                    {
                      store_link_cnsrv(rv, srcCellIndex, tgtCellIndex, num_wts, weights, grid_store);

                      tgtGrid->cell_frac[tgtCellIndex] += weights[3];
                    }
                    srcGrid->cell_frac[srcCellIndex] += weights[0];
                  }

              srcGrid->cell_area[srcCellIndex] += weights[0];
              src_centroid_lat[srcCellIndex] += weights[1];
              src_centroid_lon[srcCellIndex] += weights[2];

              /* Reset beglat and beglon for next subsegment. */
              beglat = intrsct_lat;
              beglon = intrsct_lon;
            }
          /* End of segment */
        }
    }

  /* Finished with all cells: deallocate search arrays */

  for (int i = 0; i < Threading::ompNumThreads; ++i) free(srch_corner_lon[i]);
  for (int i = 0; i < Threading::ompNumThreads; ++i) free(srch_corner_lat[i]);

  /* Integrate around each cell on target grid */

  for (int i = 0; i < Threading::ompNumThreads; ++i) srch_corner_lat[i] = nullptr;
  for (int i = 0; i < Threading::ompNumThreads; ++i) srch_corner_lon[i] = nullptr;

  for (int i = 0; i < Threading::ompNumThreads; ++i) maxSearchCells[i] = 0;
  for (int i = 0; i < Threading::ompNumThreads; ++i) srch_add[i].resize(srcGridSize);

  srch_corners = srcNumCorners;

  findex = 0;

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(findex, rsearch, num_wts, tgt_centroid_lon, tgt_centroid_lat, grid_store, rv,      \
                                              Options::cdoVerbose, max_subseg, srch_corner_lat, srch_corner_lon, maxSearchCells, \
                                              tgtNumCorners, srch_corners, srcGrid, tgtGrid, tgtGridSize, srcGridSize, srch_add)
#endif
  for (long tgtCellIndex = 0; tgtCellIndex < tgtGridSize; ++tgtCellIndex)
    {
      const auto ompthID = cdo_omp_get_thread_num();

#ifdef _OPENMP
#pragma omp atomic
#endif
      findex++;
      if (ompthID == 0) progress::update(0.5, 0.5, findex / tgtGridSize);

      bool lcoinc;  // flag for coincident segments
      bool lthresh = false;
      bool luse_last = false;
      int avoid_pole_count = 0;                         // count attempts to avoid pole
      double avoid_pole_offset = TINY;                  // endpoint offset to avoid pole
      double intrsct_lat, intrsct_lon;                  // lat/lon of next intersect
      double intrsct_lat_off = 0, intrsct_lon_off = 0;  // lat/lon coords offset for next search
      double intrsct_x, intrsct_y;                      // x,y for intersection
      long last_loc = -1;                               // save location when crossing threshold
      long srcCellIndex;

      // Get search cells
      long numSearchCells = get_srch_cells(tgtCellIndex, rsearch.tgtBins, rsearch.srcBins,
                                           &rsearch.tgtBins.cell_bound_box[tgtCellIndex * 4], srch_add[ompthID]);

      if (numSearchCells == 0) continue;

      // Create search arrays
      if (numSearchCells > maxSearchCells[ompthID])
        {
          srch_corner_lat[ompthID] = (double *) realloc(srch_corner_lat[ompthID], srch_corners * numSearchCells * sizeof(double));
          srch_corner_lon[ompthID] = (double *) realloc(srch_corner_lon[ompthID], srch_corners * numSearchCells * sizeof(double));
          maxSearchCells[ompthID] = numSearchCells;
        }

      // gather2
      long ioffset;
      for (long n = 0; n < numSearchCells; ++n)
        {
          srcCellIndex = srch_add[ompthID][n];
          ioffset = srcCellIndex * srch_corners;

          long nsrch_corners = n * srch_corners;
          for (long k = 0; k < srch_corners; ++k)
            {
              srch_corner_lat[ompthID][nsrch_corners + k] = srcGrid->cell_corner_lat[ioffset + k];
              srch_corner_lon[ompthID][nsrch_corners + k] = srcGrid->cell_corner_lon[ioffset + k];
            }
        }

      // Integrate around this cell
      ioffset = tgtCellIndex * tgtNumCorners;
      for (long corner = 0; corner < tgtNumCorners; ++corner)
        {
          long next_corn = (corner + 1) % tgtNumCorners;

          /* Define endpoints of the current segment */
          double beglat = tgtGrid->cell_corner_lat[ioffset + corner];
          double beglon = tgtGrid->cell_corner_lon[ioffset + corner];
          double endlat = tgtGrid->cell_corner_lat[ioffset + next_corn];
          double endlon = tgtGrid->cell_corner_lon[ioffset + next_corn];
          bool lrevers = false;

          /* To ensure exact path taken during both sweeps, always integrate in the same direction */
          if ((endlat < beglat) || (IS_EQUAL(endlat, beglat) && endlon < beglon))
            {
              beglat = tgtGrid->cell_corner_lat[ioffset + next_corn];
              beglon = tgtGrid->cell_corner_lon[ioffset + next_corn];
              endlat = tgtGrid->cell_corner_lat[ioffset + corner];
              endlon = tgtGrid->cell_corner_lon[ioffset + corner];
              lrevers = true;
            }

          /*
            If this is a constant-longitude segment, skip the rest
            since the line integral contribution will be ZERO.
          */
          if (IS_EQUAL(endlon, beglon)) continue;

          double weights[6];  // local wgt array
          double begseg[2];   // begin lat/lon for full segment
          begseg[0] = beglat;
          begseg[1] = beglon;
          bool lbegin = true;

          long num_subseg = 0;
          /*
            Integrate along this segment, detecting intersections
            and computing the line integral for each sub-segment
          */
          while (IS_NOT_EQUAL(beglat, endlat) || IS_NOT_EQUAL(beglon, endlon))
            {
              /* Prevent infinite loops if integration gets stuck near cell or threshold boundary */
              num_subseg++;
              if (num_subseg >= max_subseg)
                cdo_abort("Integration stalled: num_subseg exceeded limit (grid2[%zu]: lon1=%g lon2=%g lat1=%g lat2=%g)!",
                          tgtCellIndex, beglon, endlon, beglat, endlat);

              /* Uwe Schulzweida: skip very small regions */
              if (num_subseg % 1000 == 0)
                {
                  if (std::fabs(beglat - endlat) < 1.e-10 || std::fabs(beglon - endlon) < 1.e-10)
                    {
                      if (Options::cdoVerbose)
                        cdo_print("Skip very small region (grid2[%ld]): lon=%g dlon=%g lat=%g dlat=%g", tgtCellIndex, beglon,
                                  endlon - beglon, beglat, endlat - beglat);
                      break;
                    }
                }

              /* Find next intersection of this segment with a gridline on grid 2. */

              intersection(&srcCellIndex, &intrsct_lat, &intrsct_lon, &lcoinc, beglat, beglon, endlat, endlon, begseg, lbegin,
                           lrevers, numSearchCells, srch_corners, srch_add[ompthID].data(), srch_corner_lat[ompthID],
                           srch_corner_lon[ompthID], &last_loc, &lthresh, &intrsct_lat_off, &intrsct_lon_off, &luse_last,
                           &intrsct_x, &intrsct_y, &avoid_pole_count, &avoid_pole_offset);

              lbegin = false;

              /* Compute line integral for this subsegment. */

              if (srcCellIndex != -1)
                line_integral(weights, beglon, intrsct_lon, beglat, intrsct_lat, srcGrid->cell_center_lon[srcCellIndex],
                              tgtGrid->cell_center_lon[tgtCellIndex]);
              else
                line_integral(weights, beglon, intrsct_lon, beglat, intrsct_lat, tgtGrid->cell_center_lon[tgtCellIndex],
                              tgtGrid->cell_center_lon[tgtCellIndex]);

              /* If integrating in reverse order, change sign of weights */

              if (lrevers)
                for (long k = 0; k < 6; ++k) weights[k] = -weights[k];

              /*
                Store the appropriate addresses and weights.
                Also add contributions to cell areas and centroids.
                If there is a coincidence, do not store weights
                because they have been captured in the previous loop.
                The source grid mask is the master mask
              */
              if (!lcoinc && srcCellIndex != -1)
                if (srcGrid->mask[srcCellIndex])
                  {
#ifdef _OPENMP
#pragma omp critical
#endif
                    {
                      store_link_cnsrv(rv, srcCellIndex, tgtCellIndex, num_wts, weights, grid_store);

                      srcGrid->cell_frac[srcCellIndex] += weights[0];
                    }
                    tgtGrid->cell_frac[tgtCellIndex] += weights[3];
                  }

              tgtGrid->cell_area[tgtCellIndex] += weights[3];
              tgt_centroid_lat[tgtCellIndex] += weights[4];
              tgt_centroid_lon[tgtCellIndex] += weights[5];

              /* Reset beglat and beglon for next subsegment. */
              beglat = intrsct_lat;
              beglon = intrsct_lon;
            }
          /* End of segment */
        }
    }

  progress::update(0, 1, 1);

  /* Finished with all cells: deallocate search arrays */

  for (int i = 0; i < Threading::ompNumThreads; ++i) free(srch_corner_lon[i]);
  for (int i = 0; i < Threading::ompNumThreads; ++i) free(srch_corner_lat[i]);

  /*
     Correct for situations where N/S pole not explicitly included in
     grid (i.e. as a grid corner point). If pole is missing from only
     one grid, need to correct only the area and centroid of that
     grid.  If missing from both, do complete weight calculation.
  */
  correct_pole(srcGrid, tgtGrid, rv, src_centroid_lat, src_centroid_lon, tgt_centroid_lat, tgt_centroid_lon, grid_store);

  grid_store_delete(grid_store);

  if (rv.numLinks != rv.maxLinks) remap_vars_resize(rv, rv.numLinks);

  /* Finish centroid computation */

  // int nlon = srcGrid->dims[0];
  // int nlat = srcGrid->dims[1];
  // srvout("clat1.srv", nlon, nlat, src_centroid_lat);
  // srvout("clon1.srv", nlon, nlat, src_centroid_lon);
  for (long n = 0; n < srcGridSize; ++n)
    if (IS_NOT_EQUAL(srcGrid->cell_area[n], 0))
      {
        src_centroid_lat[n] /= srcGrid->cell_area[n];
        src_centroid_lon[n] /= srcGrid->cell_area[n];
      }
  // srvout("clat2.srv", nlon, nlat, src_centroid_lat);
  // srvout("clon2.srv", nlon, nlat, src_centroid_lon);

  for (long n = 0; n < tgtGridSize; ++n)
    if (IS_NOT_EQUAL(tgtGrid->cell_area[n], 0))
      {
        tgt_centroid_lat[n] /= tgtGrid->cell_area[n];
        tgt_centroid_lon[n] /= tgtGrid->cell_area[n];
      }

  /* 2010-10-08 Uwe Schulzweida: remove all links with weights < 0 */

  /*
  if ( 1 )
    {
      numLinks = rv.numLinks;

      if ( Options::cdoVerbose )
        for ( long n = 0; n < numLinks; n++ )
          printf("wts1: %d %g\n", n, rv.wts[3*n]);

      for ( long n = 0; n < numLinks; n++ )
        {
          if ( rv.wts[3*n] < 0 )
            {
              int i, n2, nd;

              for ( n2 = n+1; n2 < numLinks; n2++ )
                if ( rv.wts[3*n2] >= 0 ) break;

              nd = n2-n;
              numLinks -= nd;
              for ( i = n; i < numLinks; i++ )
                {
                  rv.wts[3*i]   = rv.wts[3*(i+nd)];
                  rv.wts[3*i+1] = rv.wts[3*(i+nd)+1];
                  rv.wts[3*i+2] = rv.wts[3*(i+nd)+2];

                  rv.srcCellIndices[i] = rv.srcCellIndices[i+nd];
                  rv.tgtCellIndices[i] = rv.tgtCellIndices[i+nd];
                }
            }
        }

     if ( Options::cdoVerbose ) cdo_print("Removed number of links = %zu", rv.numLinks - numLinks);

      rv.numLinks = numLinks;
    }
  */

  /* Include centroids in weights and normalize using destination area if requested */
  normalize_weights(tgtGrid, rv, src_centroid_lat, src_centroid_lon);

  long numLinks = rv.numLinks;

  if (Options::cdoVerbose) cdo_print("Total number of links = %zu", rv.numLinks);

  for (long n = 0; n < srcGridSize; ++n)
    if (IS_NOT_EQUAL(srcGrid->cell_area[n], 0)) srcGrid->cell_frac[n] /= srcGrid->cell_area[n];

  for (long n = 0; n < tgtGridSize; ++n)
    if (IS_NOT_EQUAL(tgtGrid->cell_area[n], 0)) tgtGrid->cell_frac[n] /= tgtGrid->cell_area[n];

  /* Perform some error checking on final weights  */

  if (lcheck)
    {
      remap_check_area(srcGridSize, srcGrid->cell_area, "Source");
      remap_check_area(tgtGridSize, tgtGrid->cell_area, "Target");

      for (long n = 0; n < srcGridSize; ++n)
        {
          if (src_centroid_lat[n] < (-PIH - .01) || src_centroid_lat[n] > (PIH + .01))
            cdo_print("Source grid centroid lat error: %ld %g", n, src_centroid_lat[n]);

          src_centroid_lat[n] = 0;
          src_centroid_lon[n] = 0;
        }

      for (long n = 0; n < tgtGridSize; ++n)
        {
          if (tgt_centroid_lat[n] < (-PIH - .01) || tgt_centroid_lat[n] > (PIH + .01))
            cdo_print("Target grid centroid lat error: %ld %g", n, tgt_centroid_lat[n]);

          tgt_centroid_lat[n] = 0;
          tgt_centroid_lon[n] = 0;
        }

      remap_vars_check_weights(rv);

      for (long n = 0; n < numLinks; ++n)
        {
          long tgtCellIndex = rv.tgtCellIndices[n];
          tgt_centroid_lat[tgtCellIndex] += rv.wts[3 * n];
        }

      /* 2012-01-24 Uwe Schulzweida: changed [tgtCellIndex] to [n] (bug fix) */
      double norm_factor = 0;  // factor for normalizing wts
      for (long n = 0; n < tgtGridSize; ++n)
        {
          if (rv.normOpt == NormOpt::DESTAREA)
            norm_factor = tgtGrid->cell_frac[n];
          else if (rv.normOpt == NormOpt::FRACAREA)
            norm_factor = ONE;
          else if (rv.normOpt == NormOpt::NONE)
            norm_factor = tgtGrid->cell_area[n];

          if (tgt_centroid_lat[n] > 0 && std::fabs(tgt_centroid_lat[n] - norm_factor) > .01)
            cdo_print("Error: sum of wts for map1 %ld %g %g", n, tgt_centroid_lat[n], norm_factor);
        }
    }  // lcheck

  if (Options::cdoVerbose) cdo_print("Cells search: %.2f seconds", cdo_get_wtime() - start);
}  // remap_conserv_weights_scrip
