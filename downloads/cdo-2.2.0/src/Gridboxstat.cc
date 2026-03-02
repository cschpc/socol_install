/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Gridboxstat    gridboxrange        Gridbox range
      Gridboxstat    gridboxmin          Gridbox minimum
      Gridboxstat    gridboxmax          Gridbox maximum
      Gridboxstat    gridboxsum          Gridbox sum
      Gridboxstat    gridboxmean         Gridbox mean
      Gridboxstat    gridboxavg          Gridbox average
      Gridboxstat    gridboxstd          Gridbox standard deviation
      Gridboxstat    gridboxstd1         Gridbox standard deviation [Normalize by (n-1)]
      Gridboxstat    gridboxvar          Gridbox variance
      Gridboxstat    gridboxvar1         Gridbox variance [Normalize by (n-1)]
*/

#include <cdi.h>

#include "process_int.h"
#include "param_conversion.h"
#include <mpim_grid.h>
#include "cdo_options.h"
#include "progress.h"
#include "cimdOmp.h"
#include "field_functions.h"

static void
genBoxGridReg2D(int gridID1, size_t xinc, size_t yinc, int gridID2)
{
  auto nlon1 = gridInqXsize(gridID1);
  auto nlat1 = gridInqYsize(gridID1);
  auto nlon2 = gridInqXsize(gridID2);
  auto nlat2 = gridInqYsize(gridID2);

  {
    Varray<double> xvals1(nlon1), yvals1(nlat1);
    Varray<double> xvals2(nlon2), yvals2(nlat2);
    gridInqXvals(gridID1, xvals1.data());
    gridInqYvals(gridID1, yvals1.data());

    size_t j = 0;
    for (size_t i = 0; i < nlon1; i += xinc)
      {
        auto i1 = i + (xinc - 1);
        if (i1 >= nlon1 - 1) i1 = nlon1 - 1;
        xvals2[j] = xvals1[i] + (xvals1[i1] - xvals1[i]) / 2.;
        j++;
      }

    j = 0;
    for (size_t i = 0; i < nlat1; i += yinc)
      {
        auto i1 = i + (yinc - 1);
        if (i1 >= nlat1 - 1) i1 = nlat1 - 1;
        yvals2[j] = yvals1[i] + (yvals1[i1] - yvals1[i]) / 2;
        j++;
      }

    gridDefXvals(gridID2, xvals2.data());
    gridDefYvals(gridID2, yvals2.data());
  }

  if (gridHasBounds(gridID1))
    {
      Varray<double> grid1_corner_lon(2 * nlon1), grid1_corner_lat(2 * nlat1);
      Varray<double> grid2_corner_lon(2 * nlon2), grid2_corner_lat(2 * nlat2);
      gridInqXbounds(gridID1, grid1_corner_lon.data());
      gridInqYbounds(gridID1, grid1_corner_lat.data());

      size_t j = 0;
      for (size_t i = 0; i < nlon1; i += xinc)
        {
          auto i1 = i + (xinc - 1);
          if (i1 >= nlon1 - 1) i1 = nlon1 - 1;
          grid2_corner_lon[2 * j] = grid1_corner_lon[2 * i];
          grid2_corner_lon[2 * j + 1] = grid1_corner_lon[2 * i1 + 1];
          j++;
        }

      j = 0;
      for (size_t i = 0; i < nlat1; i += yinc)
        {
          auto i1 = i + (yinc - 1);
          if (i1 >= nlat1 - 1) i1 = nlat1 - 1;
          grid2_corner_lat[2 * j] = grid1_corner_lat[2 * i];
          grid2_corner_lat[2 * j + 1] = grid1_corner_lat[2 * i1 + 1];
          j++;
        }

      gridDefNvertex(gridID2, 2);
      gridDefXbounds(gridID2, grid2_corner_lon.data());
      gridDefYbounds(gridID2, grid2_corner_lat.data());
    }
}

static int
find_corner(size_t g1_add, size_t off1, size_t off2, const Varray<double> &grid1_corner_lon, const Varray<double> &grid1_corner_lat)
{
  int c_flag[4] = { 0, 0, 0, 0 };
  for (int corner = 0; corner < 4; corner++)
    {
      auto lon = grid1_corner_lon[4 * g1_add + corner];
      auto lat = grid1_corner_lat[4 * g1_add + corner];
      auto g1_add2 = g1_add + off1;
      auto g1_add3 = g1_add + off2;
      for (int corner2 = 0; corner2 < 4; corner2++)
        {
          auto lon2 = grid1_corner_lon[4 * g1_add2 + corner2];
          auto lat2 = grid1_corner_lat[4 * g1_add2 + corner2];
          auto lon3 = grid1_corner_lon[4 * g1_add3 + corner2];
          auto lat3 = grid1_corner_lat[4 * g1_add3 + corner2];
          if ((IS_EQUAL(lon2, lon) && IS_EQUAL(lat2, lat)) || (IS_EQUAL(lon3, lon) && IS_EQUAL(lat3, lat))) c_flag[corner] = 1;
        }
    }

  if (c_flag[0] + c_flag[1] + c_flag[2] + c_flag[3] < 3) cdo_warning("Found two matching corners!");

  int corner = 0;
  for (corner = 0; corner < 4; corner++)
    if (!c_flag[corner]) break;

  return corner;
}

static void
genBoxGridCurv2D(int gridID1, size_t xinc, size_t yinc, int gridID2)
{
  auto nlon1 = gridInqXsize(gridID1);
  auto nlat1 = gridInqYsize(gridID1);
  auto nlon2 = gridInqXsize(gridID2);
  auto nlat2 = gridInqYsize(gridID2);
  auto gridsize1 = gridInqSize(gridID1);
  auto gridsize2 = gridInqSize(gridID2);

  auto circular = gridIsCircular(gridID1);
  double xvals2_0 = 0.0;

  bool lGridHasBounds = gridHasBounds(gridID1);

  Varray<double> xvals1(gridsize1), yvals1(gridsize1);
  Varray<double> xvals2(gridsize2), yvals2(gridsize2);
  gridInqXvals(gridID1, xvals1.data());
  gridInqYvals(gridID1, yvals1.data());

  // Convert lat/lon units if required
  cdo_grid_to_degree(gridID1, CDI_XAXIS, gridsize1, xvals1.data(), "grid center lon");
  cdo_grid_to_degree(gridID1, CDI_YAXIS, gridsize1, yvals1.data(), "grid center lat");

  Varray<double> grid1_corner_lon, grid1_corner_lat;
  Varray<double> grid2_corner_lon, grid2_corner_lat;
  if (lGridHasBounds)
    {
      grid1_corner_lon.resize(4 * gridsize1);
      grid1_corner_lat.resize(4 * gridsize1);
      grid2_corner_lon.resize(4 * gridsize2);
      grid2_corner_lat.resize(4 * gridsize2);
      gridInqXbounds(gridID1, grid1_corner_lon.data());
      gridInqYbounds(gridID1, grid1_corner_lat.data());

      // Convert lat/lon units if required
      cdo_grid_to_degree(gridID1, CDI_XAXIS, 4 * gridsize1, grid1_corner_lon.data(), "grid corner lon");
      cdo_grid_to_degree(gridID1, CDI_YAXIS, 4 * gridsize1, grid1_corner_lat.data(), "grid corner lat");
    }

  // Process grid2 bounds
  double area_norm = xinc * yinc;
  for (size_t y2 = 0; y2 < nlat2; y2++)
    {
      for (size_t x2 = 0; x2 < nlon2; x2++)
        {
          auto g2_add = (y2 * nlon2 + x2);
          double on_up, on_lo, ox_up, ox_lo, an_le, an_ri, ax_le, ax_ri;
          on_up = on_lo = 360.0;
          ox_up = ox_lo = -360.0;
          an_ri = an_le = 90.0;
          ax_ri = ax_le = -90.0;

          for (size_t y1 = y2 * yinc; y1 < yinc * (y2 + 1); y1++)
            {
              auto use_y1 = (y1 >= nlat1) ? nlat1 - 1 : y1;
              for (size_t x1 = x2 * xinc; x1 < xinc * (x2 + 1); x1++)
                {
                  auto use_x1 = x1;
                  if (x1 >= nlon1)
                    {
                      if (circular && use_y1 == y1)
                        use_y1 -= 1;
                      else
                        use_x1 = nlon1 - 1;
                    }

                  auto g1_add = (use_y1 * nlon1) + use_x1;
                  auto xval1 = xvals1[g1_add];
                  auto yval1 = yvals1[g1_add];

                  if (y1 == y2 * yinc && x1 == x2 * xinc)
                    {
                      xvals2_0 = xval1;
                      xvals2[g2_add] = xval1 / area_norm;
                      yvals2[g2_add] = yval1 / area_norm;
                    }
                  else if (std::fabs(xval1 - xvals2_0) > 270.)
                    {
                      if ((xval1 - xvals2_0) > 270.)
                        xvals2[g2_add] += (xval1 - 360.) / area_norm;
                      else if ((xval1 - xvals2_0) < -270.)
                        xvals2[g2_add] += (xval1 + 360.) / area_norm;
                      yvals2[g2_add] += yval1 / area_norm;
                    }
                  else
                    {
                      xvals2[g2_add] += xval1 / area_norm;
                      yvals2[g2_add] += yval1 / area_norm;
                    }

                  if (lGridHasBounds)
                    {
                      // upper left cell
                      if (y1 == y2 * yinc && x1 == x2 * xinc)
                        {
                          int corner = 0;
                          if (g1_add + nlon1 > gridsize1)
                            cdo_warning("Can't find cell below upper left!");
                          else
                            corner = find_corner(g1_add, 1, nlon1, grid1_corner_lon, grid1_corner_lat);
                          on_up = grid1_corner_lon[4 * g1_add + corner];
                          ax_le = grid1_corner_lat[4 * g1_add + corner];
                        }

                      // upper right cell
                      if ((y1 == y2 * yinc) && (x1 == (x2 + 1) * xinc - 1))
                        {
                          int corner = 0;
                          if (g1_add + nlon1 > gridsize1)
                            cdo_warning("Can't find cell below upper right!");
                          else
                            corner = find_corner(g1_add, -1, nlon1, grid1_corner_lon, grid1_corner_lat);
                          ox_up = grid1_corner_lon[4 * g1_add + corner];
                          ax_ri = grid1_corner_lat[4 * g1_add + corner];
                        }

                      // lower right cell
                      if ((y1 == (y2 + 1) * yinc - 1) && (x1 == (x2 + 1) * xinc - 1))
                        {
                          int corner = 0;
                          if (g1_add < nlon1)
                            cdo_warning("Can't find cell above lower right!");
                          else
                            corner = find_corner(g1_add, -1, -nlon1, grid1_corner_lon, grid1_corner_lat);
                          ox_lo = grid1_corner_lon[4 * g1_add + corner];
                          an_ri = grid1_corner_lat[4 * g1_add + corner];
                        }

                      // lower left cell
                      if ((y1 == (y2 + 1) * yinc - 1) && (x1 == x2 * xinc))
                        {
                          int corner = 0;
                          if (g1_add < nlon1)
                            cdo_warning("Can't find cell above lower left!");
                          else
                            corner = find_corner(g1_add, 1, -nlon1, grid1_corner_lon, grid1_corner_lat);
                          on_lo = grid1_corner_lon[4 * g1_add + corner];
                          an_le = grid1_corner_lat[4 * g1_add + corner];
                        }
                    }
                }  // x1
            }      // y1

          if (lGridHasBounds)
            {
              // upper left corner
              grid2_corner_lon[4 * g2_add + 3] = on_up;
              grid2_corner_lat[4 * g2_add + 3] = ax_le;
              // upper right corner
              grid2_corner_lon[4 * g2_add + 2] = ox_up;
              grid2_corner_lat[4 * g2_add + 2] = ax_ri;
              // lower right corner
              grid2_corner_lon[4 * g2_add + 1] = ox_lo;
              grid2_corner_lat[4 * g2_add + 1] = an_ri;
              // lower left corner
              grid2_corner_lon[4 * g2_add + 0] = on_lo;
              grid2_corner_lat[4 * g2_add + 0] = an_le;
            }

          //  while ( xvals2[g2_add] >  180. ) xvals2[g2_add] -= 360.;
          //  while ( xvals2[g2_add] < -180. ) xvals2[g2_add] += 360.;
        }  // x2
    }      // y2

  gridDefXvals(gridID2, xvals2.data());
  gridDefYvals(gridID2, yvals2.data());

  if (lGridHasBounds)
    {
      gridDefNvertex(gridID2, 4);
      gridDefXbounds(gridID2, grid2_corner_lon.data());
      gridDefYbounds(gridID2, grid2_corner_lat.data());
    }
}

static int
genBoxGrid(int gridID1, size_t xinc, size_t yinc)
{
  if (xinc < 1 || yinc < 1) cdo_abort("xinc and yinc must not be smaller than 1!");

  int gridID2 = -1;
  auto gridtype = gridInqType(gridID1);
  if (gridtype == GRID_GAUSSIAN || gridtype == GRID_LONLAT || gridtype == GRID_CURVILINEAR || gridtype == GRID_GENERIC)
    {
      auto nlon1 = gridInqXsize(gridID1);
      auto nlat1 = gridInqYsize(gridID1);
      if (xinc > nlon1 || yinc > nlat1) cdo_abort("xinc and/or yinc exceeds gridsize!");

      auto nlon2 = nlon1 / xinc;
      auto nlat2 = nlat1 / yinc;
      if (nlon1 % xinc) nlon2++;
      if (nlat1 % yinc) nlat2++;
      auto gridsize2 = nlon2 * nlat2;

      gridID2 = gridCreate(gridtype, gridsize2);
      gridDefXsize(gridID2, nlon2);
      gridDefYsize(gridID2, nlat2);
    }
  else { cdo_abort("Unsupported grid: %s", gridNamePtr(gridtype)); }

  if (gridtype == GRID_GAUSSIAN || gridtype == GRID_LONLAT) { genBoxGridReg2D(gridID1, xinc, yinc, gridID2); }
  else if (gridtype == GRID_GENERIC) {}
  else if (gridtype == GRID_CURVILINEAR) { genBoxGridCurv2D(gridID1, xinc, yinc, gridID2); }
  else { cdo_abort("Unsupported grid: %s", gridNamePtr(gridtype)); }

  return gridID2;
}

static void
gridboxstat(const Field &field1, Field &field2, size_t xinc, size_t yinc, int statfunc)
{
  auto useWeight = !field1.weightv.empty();

  auto boxsize = xinc * yinc;
  FieldVector fields(Threading::ompNumThreads);
  for (int i = 0; i < Threading::ompNumThreads; ++i)
    {
      fields[i].resize(boxsize);
      if (useWeight) fields[i].weightv.resize(boxsize);
      fields[i].missval = field1.missval;
    }

  auto gridID1 = field1.grid;
  auto gridID2 = field2.grid;

  auto nlon1 = gridInqXsize(gridID1);
  auto nlat1 = gridInqYsize(gridID1);

  auto nlon2 = gridInqXsize(gridID2);
  auto nlat2 = gridInqYsize(gridID2);

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
  for (size_t ig = 0; ig < nlat2 * nlon2; ++ig)
    {
      auto ompthID = cdo_omp_get_thread_num();

      auto &field = fields[ompthID];
      field.nmiss = 0;

      auto ilat = ig / nlon2;
      auto ilon = ig - ilat * nlon2;

      size_t isize = 0;
      for (size_t j = 0; j < yinc; ++j)
        {
          auto jj = ilat * yinc + j;
          if (jj >= nlat1) break;
          for (size_t i = 0; i < xinc; ++i)
            {
              auto ii = ilon * xinc + i;
              auto index = jj * nlon1 + ii;
              if (ii >= nlon1) break;
              field.vec_d[isize] = (field1.memType == MemType::Float) ? (double) field1.vec_f[index] : field1.vec_d[index];
              if (dbl_is_equal(field.vec_d[isize], field.missval)) field.nmiss++;
              if (useWeight) field.weightv[isize] = field1.weightv[index];
              isize++;
            }
        }

      field.size = isize;
      auto result = field_function(field, statfunc);
      if (field2.memType == MemType::Float)
        field2.vec_f[ig] = result;
      else
        field2.vec_d[ig] = result;
    }

  field_num_mv(field2);
}

static void
add_operators(void)
{
  // clang-format off
  cdo_operator_add("gridboxrange",  FieldFunc_Range,  0, nullptr);
  cdo_operator_add("gridboxmin",    FieldFunc_Min,    0, nullptr);
  cdo_operator_add("gridboxmax",    FieldFunc_Max,    0, nullptr);
  cdo_operator_add("gridboxsum",    FieldFunc_Sum,    0, nullptr);
  cdo_operator_add("gridboxmean",   FieldFunc_Meanw,  1, nullptr);
  cdo_operator_add("gridboxavg",    FieldFunc_Avgw,   1, nullptr);
  cdo_operator_add("gridboxvar",    FieldFunc_Varw,   1, nullptr);
  cdo_operator_add("gridboxvar1",   FieldFunc_Var1w,  1, nullptr);
  cdo_operator_add("gridboxstd",    FieldFunc_Stdw,   1, nullptr);
  cdo_operator_add("gridboxstd1",   FieldFunc_Std1w,  1, nullptr);
  cdo_operator_add("gridboxskew",   FieldFunc_Skew,   0, nullptr);
  cdo_operator_add("gridboxkurt",   FieldFunc_Kurt,   0, nullptr);
  cdo_operator_add("gridboxmedian", FieldFunc_Median, 0, nullptr);
  // clang-format on
}

void *
Gridboxstat(void *process)
{
  int lastgrid = -1;
  auto wstatus = false;

  cdo_initialize(process);

  operator_input_arg("xinc, yinc");
  operator_check_argc(2);
  auto xinc = parameter_to_int(cdo_operator_argv(0));
  auto yinc = parameter_to_int(cdo_operator_argv(1));

  add_operators();

  auto operatorID = cdo_operator_id();
  auto operfunc = cdo_operator_f1(operatorID);
  auto needWeights = (cdo_operator_f2(operatorID) != 0);

  auto streamID1 = cdo_open_read(0);

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto vlistID2 = vlistDuplicate(vlistID1);

  auto taxisID1 = vlistInqTaxis(vlistID1);
  auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  auto ngrids = vlistNgrids(vlistID1);
  if (ngrids > 1) cdo_abort("Too many different grids!");

  auto gridID1 = vlistGrid(vlistID1, 0);

  auto gridID2 = genBoxGrid(gridID1, xinc, yinc);
  for (int index = 0; index < ngrids; ++index) vlistChangeGridIndex(vlistID2, index, gridID2);

  auto streamID2 = cdo_open_write(1);
  cdo_def_vlist(streamID2, vlistID2);

  VarList varList1, varList2;
  varListInit(varList1, vlistID1);
  varListInit(varList2, vlistID2);

  auto gridsize1 = gridInqSize(gridID1);

  Field field1, field2;
  if (needWeights) field1.weightv.resize(gridsize1);

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
          field1.init(varList1[varID]);
          cdo_read_record(streamID1, field1);

          if (needWeights && field1.grid != lastgrid)
            {
              lastgrid = field1.grid;
              wstatus = gridcell_weights(field1.grid, field1.weightv);
            }
          if (wstatus != 0 && tsID == 0 && levelID == 0)
            cdo_warning("Grid cell bounds not available, using constant grid cell area weights for variable %s!",
                        varList1[varID].name);

          field2.init(varList2[varID]);
          gridboxstat(field1, field2, xinc, yinc, operfunc);

          cdo_def_record(streamID2, varID, levelID);
          cdo_write_record(streamID2, field2);
        }
      tsID++;
    }

  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  cdo_finish();

  return nullptr;
}
