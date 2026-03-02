/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Invert     invertlat       Invert latitude
      Invert     invertlon       Invert longitude
      Invert     invertlatdes    Invert latitude description
      Invert     invertlondes    Invert longitude description
      Invert     invertlatdata   Invert latitude data
      Invert     invertlondata   Invert longitude data
*/

#include <cdi.h>

#include <utility>

#include "process_int.h"
#include "matrix_view.h"

template <typename T>
static void
invert_lon_data(Varray<T> &v, size_t nlon, size_t nlat)
{
  if (nlat > 0 && nlon > 0)
    {
      Varray<T> vtmp(nlon);
      MatrixView<T> mv(v.data(), nlat, nlon);

      for (size_t ilat = 0; ilat < nlat; ilat++)
        {
          for (size_t ilon = 0; ilon < nlon; ilon++) vtmp[ilon] = mv[ilat][ilon];
          for (size_t ilon = 0; ilon < nlon / 2; ilon++) std::swap(vtmp[ilon], vtmp[nlon - ilon - 1]);
          for (size_t ilon = 0; ilon < nlon; ilon++) mv[ilat][ilon] = vtmp[ilon];
        }
    }
}

static void
invert_lon_data(Field &field)
{
  const auto nlon = gridInqXsize(field.grid);
  const auto nlat = gridInqYsize(field.grid);

  if (field.memType == MemType::Float)
    invert_lon_data(field.vec_f, nlon, nlat);
  else
    invert_lon_data(field.vec_d, nlon, nlat);
}

template <typename T>
static void
invert_lat_data(Varray<T> &v, size_t nlon, size_t nlat)
{
  if (nlat > 0 && nlon > 0)
    {
      Varray<T> vtmp(nlon);
      MatrixView<T> mv(v.data(), nlat, nlon);

      for (size_t ilat = 0; ilat < nlat / 2; ilat++)
        {
          for (size_t ilon = 0; ilon < nlon; ilon++) vtmp[ilon] = mv[ilat][ilon];
          for (size_t ilon = 0; ilon < nlon; ilon++) mv[ilat][ilon] = mv[nlat - ilat - 1][ilon];
          for (size_t ilon = 0; ilon < nlon; ilon++) mv[nlat - ilat - 1][ilon] = vtmp[ilon];
        }
    }
}

static void
invert_lat_data(Field &field)
{
  const auto nlon = gridInqXsize(field.grid);
  const auto nlat = gridInqYsize(field.grid);

  if (field.memType == MemType::Float)
    invert_lat_data(field.vec_f, nlon, nlat);
  else
    invert_lat_data(field.vec_d, nlon, nlat);
}

static void
invert_lon_des(int vlistID)
{
  const auto ngrids = vlistNgrids(vlistID);
  for (int index = 0; index < ngrids; ++index)
    {
      const auto gridID1 = vlistGrid(vlistID, index);
      const auto gridID2 = gridDuplicate(gridID1);

      const auto gridtype = gridInqType(gridID1);

      if (!(gridtype == GRID_GENERIC || gridtype == GRID_GAUSSIAN || gridtype == GRID_PROJECTION || gridtype == GRID_LONLAT
            || gridtype == GRID_CURVILINEAR))
        cdo_abort("Unsupported gridtype: %s!", gridNamePtr(gridtype));

      if (gridInqXvals(gridID1, nullptr))
        {
          const auto nlon = gridInqXsize(gridID1);
          const auto nlat = gridInqYsize(gridID1);
          const auto size = (gridtype == GRID_CURVILINEAR) ? nlon * nlat : nlon;

          Varray<double> coords(size);

          if (gridtype == GRID_CURVILINEAR)
            {
              gridInqXvals(gridID1, coords.data());
              invert_lon_data(coords, nlon, nlat);
              gridDefXvals(gridID2, coords.data());

              if (gridInqYvals(gridID1, nullptr))
                {
                  gridInqYvals(gridID1, coords.data());
                  invert_lon_data(coords, nlon, nlat);
                  gridDefYvals(gridID2, coords.data());
                }
            }
          else
            {
              gridInqXvals(gridID1, coords.data());
              for (size_t ilon = 0; ilon < nlon / 2; ilon++) std::swap(coords[ilon], coords[nlon - ilon - 1]);
              gridDefXvals(gridID2, coords.data());
            }
        }

      if (gridInqXbounds(gridID1, nullptr))
        {
          const auto nlon = gridInqXsize(gridID1);
          const auto nlat = gridInqYsize(gridID1);
          const auto nv = gridInqNvertex(gridID1);
          const auto size = (gridtype == GRID_CURVILINEAR) ? nv * nlon * nlat : nv * nlon;

          Varray<double> bounds(size);

          if (gridtype == GRID_CURVILINEAR)
            {
              gridInqXbounds(gridID1, bounds.data());
              invert_lon_data(bounds, nlon * nv, nlat);
              gridDefXbounds(gridID2, bounds.data());

              if (gridInqYbounds(gridID1, nullptr))
                {
                  gridInqYbounds(gridID1, bounds.data());
                  invert_lon_data(bounds, nlon * nv, nlat);
                  gridDefYbounds(gridID2, bounds.data());
                }
            }
          else
            {
              gridInqXbounds(gridID1, bounds.data());
              for (size_t ilon = 0; ilon < nlon / 2; ilon++)
                {
                  std::swap(bounds[nlon * 2 - ilon * 2 - 1], bounds[ilon * 2]);
                  std::swap(bounds[nlon * 2 - ilon * 2 - 2], bounds[ilon * 2 + 1]);
                }
              gridDefXbounds(gridID2, bounds.data());
            }
        }

      vlistChangeGrid(vlistID, gridID1, gridID2);
    }
}

static void
invert_lat_coord(int gridID)
{
  const auto gridtype = gridInqType(gridID);

  if (gridInqYvals(gridID, nullptr))
    {
      const auto nlon = gridInqXsize(gridID);
      const auto nlat = gridInqYsize(gridID);
      const auto size = (gridtype == GRID_CURVILINEAR) ? nlon * nlat : nlat;

      Varray<double> coords(size);

      if (gridtype == GRID_CURVILINEAR)
        {
          if (gridInqXvals(gridID, nullptr))
            {
              gridInqXvals(gridID, coords.data());
              invert_lat_data(coords, nlon, nlat);
              gridDefXvals(gridID, coords.data());
            }

          gridInqYvals(gridID, coords.data());
          invert_lat_data(coords, nlon, nlat);
          gridDefYvals(gridID, coords.data());
        }
      else
        {
          gridInqYvals(gridID, coords.data());
          for (size_t ilat = 0; ilat < nlat / 2; ilat++) std::swap(coords[ilat], coords[nlat - ilat - 1]);
          gridDefYvals(gridID, coords.data());
        }
    }

  if (gridInqYbounds(gridID, nullptr))
    {
      const auto nlon = gridInqXsize(gridID);
      const auto nlat = gridInqYsize(gridID);
      const auto nv = gridInqNvertex(gridID);
      const auto size = (gridtype == GRID_CURVILINEAR) ? nv * nlon * nlat : nv * nlat;

      Varray<double> bounds(size);

      if (gridtype == GRID_CURVILINEAR)
        {
          if (gridInqXbounds(gridID, nullptr))
            {
              gridInqXbounds(gridID, bounds.data());
              invert_lat_data(bounds, nlon * nv, nlat);
              gridDefXbounds(gridID, bounds.data());
            }

          gridInqYbounds(gridID, bounds.data());
          invert_lat_data(bounds, nlon * nv, nlat);
          gridDefYbounds(gridID, bounds.data());
        }
      else
        {
          gridInqYbounds(gridID, bounds.data());
          for (size_t ilat = 0; ilat < nlat / 2; ilat++)
            {
              std::swap(bounds[nlat * 2 - ilat * 2 - 1], bounds[ilat * 2]);
              std::swap(bounds[nlat * 2 - ilat * 2 - 2], bounds[ilat * 2 + 1]);
            }
          gridDefYbounds(gridID, bounds.data());
        }
    }
}

static void
invert_lat_des(int vlistID)
{
  const auto ngrids = vlistNgrids(vlistID);
  for (int index = 0; index < ngrids; ++index)
    {
      const auto gridID1 = vlistGrid(vlistID, index);
      const auto gridID2 = gridDuplicate(gridID1);

      const auto gridtype = gridInqType(gridID1);

      if (!(gridtype == GRID_GENERIC || gridtype == GRID_GAUSSIAN || gridtype == GRID_PROJECTION || gridtype == GRID_LONLAT
            || gridtype == GRID_CURVILINEAR))
        cdo_abort("Unsupported gridtype: %s!", gridNamePtr(gridtype));

      invert_lat_coord(gridID2);

      const auto projID = gridInqProj(gridID2);
      if (projID != CDI_UNDEFID) invert_lat_coord(projID);

      vlistChangeGrid(vlistID, gridID1, gridID2);
    }
}

void *
Invert(void *process)
{
  enum
  {
    func_fld,
    func_all,
    func_hrd,
    func_lon,
    func_lat
  };

  cdo_initialize(process);

  // clang-format off
  cdo_operator_add("invertlat",     func_all, func_lat, nullptr);
  cdo_operator_add("invertlon",     func_all, func_lon, nullptr);
  cdo_operator_add("invertlatdes",  func_hrd, func_lat, nullptr);
  cdo_operator_add("invertlondes",  func_hrd, func_lon, nullptr);
  cdo_operator_add("invertlatdata", func_fld, func_lat, nullptr);
  cdo_operator_add("invertlondata", func_fld, func_lon, nullptr);
  // clang-format on

  const auto operatorID = cdo_operator_id();
  const auto operfunc1 = cdo_operator_f1(operatorID);
  const auto operfunc2 = cdo_operator_f2(operatorID);

  operator_check_argc(0);

  const auto streamID1 = cdo_open_read(0);

  const auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  const auto vlistID2 = vlistDuplicate(vlistID1);

  const auto taxisID1 = vlistInqTaxis(vlistID1);
  const auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  if (operfunc1 == func_all || operfunc1 == func_hrd)
    {
      if (operfunc2 == func_lat)
        invert_lat_des(vlistID2);
      else
        invert_lon_des(vlistID2);
    }

  const auto streamID2 = cdo_open_write(1);

  cdo_def_vlist(streamID2, vlistID2);

  VarList varList1;
  varListInit(varList1, vlistID1);

  Field field;

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
          field.init(varList1[varID]);
          cdo_read_record(streamID1, field);

          cdo_def_record(streamID2, varID, levelID);

          if (operfunc1 == func_all || operfunc1 == func_fld)
            {
              if (operfunc2 == func_lat)
                invert_lat_data(field);
              else
                invert_lon_data(field);

              cdo_write_record(streamID2, field);
            }
          else { cdo_write_record(streamID2, field); }
        }

      tsID++;
    }

  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  cdo_finish();

  return nullptr;
}
