/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Intgrid    interpolate     PINGO grid interpolation
      Intgrid    intgridbil      Bilinear grid interpolation
*/

#include <cdi.h>

#include "process_int.h"
#include "param_conversion.h"
#include "interpol.h"
#include <mpim_grid.h>
#include "griddes.h"
#include "matrix_view.h"

template <typename T>
static void
thinout(const Varray<T> &varray1, Varray<T> &varray2, int gridID1, int gridID2, size_t xinc, size_t yinc)
{
  const auto nlon1 = gridInqXsize(gridID1);
  const auto nlat1 = gridInqYsize(gridID1);

  const auto nlon2 = gridInqXsize(gridID2);
  const auto nlat2 = gridInqYsize(gridID2);

  MatrixView<const T> xfield1(varray1.data(), nlat1, nlon1);
  MatrixView<T> xfield2(varray2.data(), nlat2, nlon2);

  size_t olat = 0;
  for (size_t ilat = 0; ilat < nlat1; ilat += yinc)
    {
      size_t olon = 0;
      for (size_t ilon = 0; ilon < nlon1; ilon += xinc)
        {
          xfield2[olat][olon] = xfield1[ilat][ilon];
          olon++;
        }
      olat++;
    }
}

static void
thinout(const Field &field1, Field &field2, size_t xinc, size_t yinc)
{
  if (field1.memType != field2.memType) cdo_abort("Interal error, memType of field1 and field2 differ!");

  if (field1.memType == MemType::Float)
    thinout(field1.vec_f, field2.vec_f, field1.grid, field2.grid, xinc, yinc);
  else
    thinout(field1.vec_d, field2.vec_d, field1.grid, field2.grid, xinc, yinc);

  field_num_mv(field2);
}

static int
gen_thinout_grid(int gridID1, size_t xinc, size_t yinc)
{
  const auto nlon1 = gridInqXsize(gridID1);
  const auto nlat1 = gridInqYsize(gridID1);

  auto nlon2 = nlon1 / xinc;
  auto nlat2 = nlat1 / yinc;
  if (nlon1 % xinc) nlon2++;
  if (nlat1 % yinc) nlat2++;
  const auto gridsize2 = nlon2 * nlat2;

  const auto gridtype = gridInqType(gridID1);
  const auto gridID2 = gridCreate(gridtype, gridsize2);
  gridDefXsize(gridID2, nlon2);
  gridDefYsize(gridID2, nlat2);

  grid_copy_names(gridID1, gridID2);

  if (gridtype == GRID_GAUSSIAN || gridtype == GRID_LONLAT)
    {
      Varray<double> xvals1(nlon1), yvals1(nlat1);
      Varray<double> xvals2(nlon2), yvals2(nlat2);
      gridInqXvals(gridID1, xvals1.data());
      gridInqYvals(gridID1, yvals1.data());

      size_t olat = 0;
      for (size_t ilat = 0; ilat < nlat1; ilat += yinc) yvals2[olat++] = yvals1[ilat];

      size_t olon = 0;
      for (size_t ilon = 0; ilon < nlon1; ilon += xinc) xvals2[olon++] = xvals1[ilon];

      gridDefXvals(gridID2, xvals2.data());
      gridDefYvals(gridID2, yvals2.data());
    }
  else if (gridtype == GRID_CURVILINEAR)
    {
      Varray<double> xvals1(nlon1 * nlat1), yvals1(nlon1 * nlat1);
      Varray<double> xvals2(nlon2 * nlat2), yvals2(nlon2 * nlat2);
      gridInqXvals(gridID1, xvals1.data());
      gridInqYvals(gridID1, yvals1.data());

      thinout(xvals1, xvals2, gridID1, gridID2, xinc, yinc);
      thinout(yvals1, yvals2, gridID1, gridID2, xinc, yinc);

      gridDefXvals(gridID2, xvals2.data());
      gridDefYvals(gridID2, yvals2.data());
    }
  else { cdo_abort("Unsupported grid: %s", gridNamePtr(gridtype)); }

  return gridID2;
}

static int
gen_boxavg_grid(int gridID1, size_t xinc, size_t yinc)
{
  const auto nlon1 = gridInqXsize(gridID1);
  const auto nlat1 = gridInqYsize(gridID1);

  auto nlon2 = nlon1 / xinc;
  auto nlat2 = nlat1 / yinc;
  if (nlon1 % xinc) nlon2++;
  if (nlat1 % yinc) nlat2++;
  const auto gridsize2 = nlon2 * nlat2;

  const auto gridID2 = gridCreate(GRID_LONLAT, gridsize2);
  gridDefXsize(gridID2, nlon2);
  gridDefYsize(gridID2, nlat2);

  const auto gridtype = gridInqType(gridID1);
  if (gridtype == GRID_GAUSSIAN || gridtype == GRID_LONLAT)
    {
      Varray<double> xvals1(nlon1), yvals1(nlat1);
      Varray<double> xvals2(nlon2), yvals2(nlat2);
      gridInqXvals(gridID1, &xvals1[0]);
      gridInqYvals(gridID1, &yvals1[0]);

      Varray<double> grid1_corner_lon, grid1_corner_lat;
      Varray<double> grid2_corner_lon, grid2_corner_lat;
      if (gridHasBounds(gridID1))
        {
          grid1_corner_lon.resize(2 * nlon1);
          grid1_corner_lat.resize(2 * nlat1);
          grid2_corner_lon.resize(2 * nlon2);
          grid2_corner_lat.resize(2 * nlat2);
          gridInqXbounds(gridID1, &grid1_corner_lon[0]);
          gridInqYbounds(gridID1, &grid1_corner_lat[0]);
        }

      size_t j = 0;
      for (size_t i = 0; i < nlon1; i += xinc)
        {
          auto i1 = i + (xinc - 1);
          if (i1 >= nlon1 - 1) i1 = nlon1 - 1;
          xvals2[j] = xvals1[i] + (xvals1[i1] - xvals1[i]) / 2;
          if (!grid2_corner_lon.empty())
            {
              grid2_corner_lon[2 * j] = grid1_corner_lon[2 * i];
              grid2_corner_lon[2 * j + 1] = grid1_corner_lon[2 * i1 + 1];
            }
          j++;
        }
      j = 0;
      for (size_t i = 0; i < nlat1; i += yinc)
        {
          auto i1 = i + (yinc - 1);
          if (i1 >= nlat1 - 1) i1 = nlat1 - 1;
          yvals2[j] = yvals1[i] + (yvals1[i1] - yvals1[i]) / 2;
          if (!grid2_corner_lat.empty())
            {
              grid2_corner_lat[2 * j] = grid1_corner_lat[2 * i];
              grid2_corner_lat[2 * j + 1] = grid1_corner_lat[2 * i1 + 1];
            }
          j++;
        }

      gridDefXvals(gridID2, &xvals2[0]);
      gridDefYvals(gridID2, &yvals2[0]);

      if (!grid2_corner_lon.empty() && !grid2_corner_lat.empty())
        {
          gridDefNvertex(gridID2, 2);
          gridDefXbounds(gridID2, &grid2_corner_lon[0]);
          gridDefYbounds(gridID2, &grid2_corner_lat[0]);
        }
    }
  else { cdo_abort("Unsupported grid: %s", gridNamePtr(gridtype)); }

  return gridID2;
}

template <typename T>
static void
boxavg(const Varray<T> &varray1, Varray<T> &varray2, int gridID1, int gridID2, size_t xinc, size_t yinc)
{
  const auto nlon1 = gridInqXsize(gridID1);
  const auto nlat1 = gridInqYsize(gridID1);

  const auto nlon2 = gridInqXsize(gridID2);
  const auto nlat2 = gridInqYsize(gridID2);

  MatrixView<const T> xfield1(varray1.data(), nlat1, nlon1);
  MatrixView<T> xfield2(varray2.data(), nlat2, nlon2);

  for (size_t ilat = 0; ilat < nlat2; ilat++)
    for (size_t ilon = 0; ilon < nlon2; ilon++)
      {
        double xsum = 0.0;

        size_t in = 0;
        for (size_t j = 0; j < yinc; ++j)
          {
            const auto jj = ilat * yinc + j;
            if (jj >= nlat1) break;
            for (size_t i = 0; i < xinc; ++i)
              {
                const auto ii = ilon * xinc + i;
                if (ii >= nlon1) break;
                in++;
                xsum += xfield1[jj][ii];
              }
          }

        xfield2[ilat][ilon] = xsum / in;
      }
}

static void
boxavg(const Field &field1, Field &field2, size_t xinc, size_t yinc)
{
  if (field1.memType != field2.memType) cdo_abort("Interal error, memType of field1 and field2 differ!");

  if (field1.memType == MemType::Float)
    boxavg(field1.vec_f, field2.vec_f, field1.grid, field2.grid, xinc, yinc);
  else
    boxavg(field1.vec_d, field2.vec_d, field1.grid, field2.grid, xinc, yinc);

  field_num_mv(field2);
}

void *
Intgrid(void *process)
{
  int gridID2 = -1;
  int xinc = 1, yinc = 1;

  cdo_initialize(process);

  // clang-format off
  const auto INTGRIDBIL  = cdo_operator_add("intgridbil",  0, 0, nullptr);
  const auto INTGRIDDIS  = cdo_operator_add("intgriddis",  0, 0, nullptr);
  const auto INTGRIDNN   = cdo_operator_add("intgridnn",   0, 0, nullptr);
  const auto INTERPOLATE = cdo_operator_add("interpolate", 0, 0, nullptr);
  const auto BOXAVG      = cdo_operator_add("boxavg",      0, 0, nullptr);
  const auto THINOUT     = cdo_operator_add("thinout",     0, 0, nullptr);
  // clang-format on

  const auto operatorID = cdo_operator_id();

  const auto lfieldmem = (operatorID == INTGRIDBIL || operatorID == THINOUT || operatorID == BOXAVG);

  if (operatorID == INTGRIDBIL || operatorID == INTERPOLATE || operatorID == INTGRIDDIS || operatorID == INTGRIDNN)
    {
      operator_input_arg("grid description file or name");
      gridID2 = cdo_define_grid(cdo_operator_argv(0));
    }
  else if (operatorID == THINOUT || operatorID == BOXAVG)
    {
      operator_input_arg("xinc, yinc");
      operator_check_argc(2);
      xinc = parameter_to_int(cdo_operator_argv(0));
      yinc = parameter_to_int(cdo_operator_argv(1));
    }

  const auto streamID1 = cdo_open_read(0);

  const auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  const auto vlistID2 = vlistDuplicate(vlistID1);

  const auto taxisID1 = vlistInqTaxis(vlistID1);
  const auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  const auto ngrids = vlistNgrids(vlistID1);
  for (int index = 0; index < ngrids; ++index)
    {
      const auto gridID1 = vlistGrid(vlistID1, index);
      const auto gridtype = gridInqType(gridID1);

      if (operatorID == BOXAVG)
        {
          if (index == 0)
            {
              if (gridtype != GRID_LONLAT && gridtype != GRID_GAUSSIAN)
                cdo_abort("Interpolation of %s data unsupported!", gridNamePtr(gridtype));

              gridID2 = gen_boxavg_grid(gridID1, xinc, yinc);
            }
          else
            cdo_abort("Too many different grids!");
        }
      if (operatorID == THINOUT)
        {
          if (index == 0)
            {
              if (gridtype != GRID_LONLAT && gridtype != GRID_GAUSSIAN && gridtype != GRID_CURVILINEAR)
                cdo_abort("Interpolation of %s data unsupported!", gridNamePtr(gridtype));

              gridID2 = gen_thinout_grid(gridID1, xinc, yinc);
            }
          else
            cdo_abort("Too many different grids!");
        }
      else if (operatorID == INTGRIDBIL || operatorID == INTERPOLATE)
        {
          const auto ldistgen = (grid_is_distance_generic(gridID1) && grid_is_distance_generic(gridID2));
          if (!ldistgen && gridtype != GRID_LONLAT && gridtype != GRID_GAUSSIAN)
            cdo_abort("Interpolation of %s data unsupported!", gridNamePtr(gridtype));
        }
      else if (operatorID == INTGRIDNN || operatorID == INTGRIDDIS)
        {
          const auto hasProjParams = ((gridtype == GRID_PROJECTION) && grid_has_proj_params(gridID1));
          if (!gridProjIsSupported(gridID1) && !hasProjParams && gridtype != GRID_LONLAT && gridtype != GRID_GAUSSIAN
              && gridtype != GRID_GME && gridtype != GRID_CURVILINEAR && gridtype != GRID_UNSTRUCTURED)
            cdo_abort("Interpolation of %s data unsupported!", gridNamePtr(gridtype));
        }

      vlistChangeGridIndex(vlistID2, index, gridID2);
    }

  VarList varList1, varList2;
  varListInit(varList1, vlistID1);
  varListInit(varList2, vlistID2);

  Field field1, field2;
  const auto gridsizemax = vlistGridsizeMax(vlistID1);
  if (!lfieldmem) field1.resize(gridsizemax);

  const auto gridsize = gridInqSize(gridID2);
  if (!lfieldmem) field2.resize(gridsize);

  const auto streamID2 = cdo_open_write(1);
  cdo_def_vlist(streamID2, vlistID2);

  int tsID = 0;
  while (true)
    {
      const auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
      if (nrecs == 0) break;

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_def_timestep(streamID2, tsID);

      for (int recID = 0; recID < nrecs; ++recID)
        {
          int varID, levelID;
          cdo_inq_record(streamID1, &varID, &levelID);
          if (lfieldmem)
            {
              field1.init(varList1[varID]);
              cdo_read_record(streamID1, field1);

              field2.init(varList2[varID]);
            }
          else
            {
              cdo_read_record(streamID1, field1.vec_d.data(), &field1.nmiss);

              field1.grid = varList1[varID].gridID;
              field1.missval = varList1[varID].missval;
              field2.grid = gridID2;
              field2.missval = field1.missval;
              field2.nmiss = 0;
            }

          // clang-format off
	  if      (operatorID == INTGRIDBIL)  intgridbil(field1, field2);
	  else if (operatorID == INTGRIDNN)   intgridnn(field1, field2);
	  else if (operatorID == INTGRIDDIS)  intgriddis(field1, field2, 4);
	  else if (operatorID == INTERPOLATE) interpolate(field1, field2);
	  else if (operatorID == BOXAVG)      boxavg(field1, field2, xinc, yinc);
	  else if (operatorID == THINOUT)     thinout(field1, field2, xinc, yinc);
          // clang-format on

          cdo_def_record(streamID2, varID, levelID);
          if (lfieldmem)
            cdo_write_record(streamID2, field2);
          else
            cdo_write_record(streamID2, field2.vec_d.data(), field2.nmiss);
        }
      tsID++;
    }

  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  cdo_finish();

  return nullptr;
}
