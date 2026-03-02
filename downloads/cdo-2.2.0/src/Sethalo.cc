/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

*/

#include <cdi.h>

#include <utility>

#include "process_int.h"
#include "param_conversion.h"
#include "util_string.h"
#include "pmlist.h"
#include <mpim_grid.h>

constexpr double UndefValue = DBL_MAX;

struct Halo
{
  long east = 0;
  long west = 0;
  long south = 0;
  long north = 0;
  double value = UndefValue;
};

static int
gen_tripolar_grid(int gridID1)
{
  auto nlon1 = gridInqXsize(gridID1);
  auto nlat1 = gridInqYsize(gridID1);

  auto nlon2 = nlon1;
  auto nlat2 = nlat1 + 2;

  auto gridtype = gridInqType(gridID1);

  auto gridID2 = gridCreate(gridtype, nlon2 * nlat2);
  gridDefXsize(gridID2, nlon2);
  gridDefYsize(gridID2, nlat2);

  cdiCopyKey(gridID1, CDI_GLOBAL, CDI_KEY_DATATYPE, gridID2);

  grid_copy_names(gridID1, gridID2);

  // auto xunits = cdo::inq_key_string(gridID1, CDI_XAXIS, CDI_KEY_UNITS);
  // auto yunits = cdo::inq_key_string(gridID1, CDI_YAXIS, CDI_KEY_UNITS);

  if (gridHasCoordinates(gridID1))
    {
      if (gridtype == GRID_CURVILINEAR)
        {
          Varray<double> xvals1(nlon1 * nlat1), yvals1(nlon1 * nlat1);
          Varray<double> xvals2(nlon2 * nlat2), yvals2(nlon2 * nlat2);

          gridInqXvals(gridID1, xvals1.data());
          gridInqYvals(gridID1, yvals1.data());

          for (size_t ilat = 0; ilat < nlat1; ilat++)
            {
              for (size_t ilon = 0; ilon < nlon1; ilon++)
                {
                  xvals2[(ilat + 2) * nlon1 + ilon] = xvals1[ilat * nlon1 + ilon];
                  yvals2[(ilat + 2) * nlon1 + ilon] = yvals1[ilat * nlon1 + ilon];
                }
            }

          for (size_t ilon = 0; ilon < nlon1; ilon++)
            {
              auto ilonr = nlon1 - ilon - 1;
              xvals2[1 * nlon1 + ilon] = xvals2[2 * nlon1 + ilonr];  // syncronise line 2 with line 3
              xvals2[0 * nlon1 + ilon] = xvals2[3 * nlon1 + ilonr];  // syncronise line 1 with line 4
              yvals2[1 * nlon1 + ilon] = yvals2[2 * nlon1 + ilonr];  // syncronise line 2 with line 3
              yvals2[0 * nlon1 + ilon] = yvals2[3 * nlon1 + ilonr];  // syncronise line 1 with line 4
            }

          gridDefXvals(gridID2, xvals2.data());
          gridDefYvals(gridID2, yvals2.data());
        }
    }

  if (gridHasBounds(gridID1))
    {
      if (gridtype == GRID_CURVILINEAR)
        {
          Varray<double> xbounds1(4 * nlon1 * nlat1), ybounds1(4 * nlon1 * nlat1);
          Varray<double> xbounds2(4 * nlon2 * nlat2), ybounds2(4 * nlon2 * nlat2);

          gridInqXbounds(gridID1, xbounds1.data());
          gridInqYbounds(gridID1, ybounds1.data());

          if (gridtype == GRID_CURVILINEAR)
            {
              gridDefNvertex(gridID2, 4);

              auto nlon4 = 4 * nlon1;
              for (size_t ilat = 0; ilat < nlat1; ilat++)
                {
                  for (size_t ilon = 0; ilon < nlon4; ilon++)
                    {
                      xbounds2[(ilat + 2) * nlon4 + ilon] = xbounds1[ilat * nlon4 + ilon];
                      ybounds2[(ilat + 2) * nlon4 + ilon] = ybounds1[ilat * nlon4 + ilon];
                    }
                }

              for (size_t ilon = 0; ilon < nlon1; ilon++)
                {
                  auto ilonr = nlon1 - ilon - 1;
                  for (size_t k = 0; k < 4; ++k)
                    {
                      auto kr = 3 - k;
                      xbounds2[1 * nlon4 + 4 * ilon + k] = xbounds2[2 * nlon4 + 4 * ilonr + kr];
                      xbounds2[0 * nlon4 + 4 * ilon + k] = xbounds2[3 * nlon4 + 4 * ilonr + kr];
                      ybounds2[1 * nlon4 + 4 * ilon + k] = ybounds2[2 * nlon4 + 4 * ilonr + kr];
                      ybounds2[0 * nlon4 + 4 * ilon + k] = ybounds2[3 * nlon4 + 4 * ilonr + kr];
                    }
                }
              /*
              for (size_t ilon = 0; ilon < 4*nlon1; ilon++)
                {
                  auto ilonr = nlon4 - ilon - 1;
                  xbounds2[1*nlon4 + ilon] = xbounds2[2*nlon4 + ilonr]; xbounds2[0*nlon4 + ilon] = xbounds2[3*nlon4 + ilonr];
                  ybounds2[1*nlon4 + ilon] = ybounds2[2*nlon4 + ilonr]; ybounds2[0*nlon4 + ilon] = ybounds2[3*nlon4 + ilonr];
                }
              */
            }

          gridDefXbounds(gridID2, xbounds2.data());
          gridDefYbounds(gridID2, ybounds2.data());
        }
    }

  return gridID2;
}

static int
gen_regular_grid(int gridID1, Halo &halo)
{
  double cpi2 = M_PI * 2;

  auto isCircularGrid = gridIsCircular(gridID1);

  long nlon1 = gridInqXsize(gridID1);
  long nlat1 = gridInqYsize(gridID1);

  long nlon2 = nlon1 + halo.east + halo.west;
  long nlat2 = nlat1 + halo.south + halo.north;

  long lonMinIdx = 0;
  long lonMaxIdx = nlon1;
  if (halo.east < 0) lonMinIdx = -halo.east;
  if (halo.west < 0) lonMaxIdx += halo.west;
  long latMinIdx = 0;
  long latMaxIdx = nlat1;
  if (halo.south < 0) latMinIdx = -halo.south;
  if (halo.north < 0) latMaxIdx += halo.north;

  if (Options::cdoVerbose)
    printf("nlon1=%ld/nlat1=%ld nlon2=%ld/nlat2=%ld east=%ld/west=%ld/south=%ld/north=%ld\n", nlon1, nlat1, nlon2, nlat2, halo.east,
           halo.west, halo.south, halo.north);

  auto gridtype = gridInqType(gridID1);
  auto isCurvilinearGrid = (gridtype == GRID_CURVILINEAR);

  auto xinc = isCurvilinearGrid ? 0.0 : gridInqXinc(gridID1);
  auto yinc = isCurvilinearGrid ? 0.0 : gridInqYinc(gridID1);

  auto gridID2 = gridCreate(gridtype, nlon2 * nlat2);
  gridDefXsize(gridID2, nlon2);
  gridDefYsize(gridID2, nlat2);

  cdiCopyKey(gridID1, CDI_GLOBAL, CDI_KEY_DATATYPE, gridID2);

  grid_copy_names(gridID1, gridID2);

  auto xunits = cdo::inq_key_string(gridID1, CDI_XAXIS, CDI_KEY_UNITS);
  // auto yunits = cdo::inq_key_string(gridID1, CDI_YAXIS, CDI_KEY_UNITS);

  if (xunits.rfind("degree", 0) == 0) cpi2 *= RAD2DEG;

  if (gridHasCoordinates(gridID1))
    {
      Varray<double> xvals1(isCurvilinearGrid ? nlon1 * nlat1 : nlon1);
      Varray<double> yvals1(isCurvilinearGrid ? nlon1 * nlat1 : nlat1);
      Varray<double> xvals2(isCurvilinearGrid ? nlon2 * nlat2 : nlon2);
      Varray<double> yvals2(isCurvilinearGrid ? nlon2 * nlat2 : nlat2);

      auto pxvals2 = xvals2.data();
      auto pyvals2 = yvals2.data();

      gridInqXvals(gridID1, xvals1.data());
      gridInqYvals(gridID1, yvals1.data());

      if (isCurvilinearGrid)
        {
          if (yvals1[0] > yvals1[nlon1 * nlat1 - 1]) std::swap(halo.south, halo.north);

          if (halo.south > 0)
            {
              pxvals2 += nlon2 * halo.south;
              pyvals2 += nlon2 * halo.south;
            }

          for (long ilat = latMinIdx; ilat < latMaxIdx; ilat++)
            {
              const auto pxvals1 = &xvals1[ilat * nlon1];
              const auto pyvals1 = &yvals1[ilat * nlon1];
              if (isCircularGrid)
                {
                  for (long ilon = nlon1 - halo.east; ilon < nlon1; ilon++) *pxvals2++ = pxvals1[ilon];
                  for (long ilon = nlon1 - halo.east; ilon < nlon1; ilon++) *pyvals2++ = pyvals1[ilon];
                }
              else
                {
                  for (long ilon = nlon1 - halo.east; ilon < nlon1; ilon++) *pxvals2++ = pxvals1[0];
                  for (long ilon = nlon1 - halo.east; ilon < nlon1; ilon++) *pyvals2++ = pyvals1[0];
                }

              for (long ilon = lonMinIdx; ilon < lonMaxIdx; ilon++) *pxvals2++ = pxvals1[ilon];
              for (long ilon = lonMinIdx; ilon < lonMaxIdx; ilon++) *pyvals2++ = pyvals1[ilon];

              if (isCircularGrid)
                {
                  for (long ilon = 0; ilon < halo.west; ilon++) *pxvals2++ = pxvals1[ilon];
                  for (long ilon = 0; ilon < halo.west; ilon++) *pyvals2++ = pyvals1[ilon];
                }
              else
                {
                  for (long ilon = 0; ilon < halo.west; ilon++) *pxvals2++ = pxvals1[nlon1 - 1];
                  for (long ilon = 0; ilon < halo.west; ilon++) *pyvals2++ = pyvals1[nlon1 - 1];
                }
            }
          for (long ilat = 0; ilat < halo.north; ilat++)
            {
              long offset = nlon2 * (nlat1 + halo.south - 1);
              for (long ilon = 0; ilon < nlon2; ilon++) *pxvals2++ = xvals2[offset + ilon];
              for (long ilon = 0; ilon < nlon2; ilon++) *pyvals2++ = yvals2[offset + ilon];
            }
          for (long ilat = 0; ilat < halo.south; ilat++)
            {
              long offset = nlon2 * halo.south;
              for (long ilon = 0; ilon < nlon2; ilon++) xvals2[nlon2 * ilat + ilon] = xvals2[offset + ilon];
              for (long ilon = 0; ilon < nlon2; ilon++) yvals2[nlon2 * ilat + ilon] = yvals2[offset + ilon];
            }
        }
      else
        {
          if (yvals1[0] > yvals1[nlat1 - 1]) std::swap(halo.south, halo.north);

          if (isCircularGrid)
            {
              // clang-format off
              for (long ilon = nlon1 - halo.east; ilon < nlon1; ilon++) *pxvals2++ = xvals1[ilon] - cpi2;
              for (long ilon = lonMinIdx; ilon < lonMaxIdx; ilon++)     *pxvals2++ = xvals1[ilon];
              for (long ilon = 0; ilon < halo.west; ilon++)             *pxvals2++ = xvals1[ilon] + cpi2;
              // clang-format on
            }
          else
            {
              // clang-format off
              for (long ilon = halo.east; ilon > 0; ilon--)             *pxvals2++ = xvals1[0] - ilon * xinc;
              for (long ilon = lonMinIdx; ilon < lonMaxIdx; ilon++)     *pxvals2++ = xvals1[ilon];
              for (long ilon = 1; ilon <= halo.west; ilon++)            *pxvals2++ = xvals1[nlon1 - 1] + ilon * xinc;
              // clang-format on
            }

          // clang-format off
          for (long ilat = halo.south; ilat > 0; ilat--)        *pyvals2++ = yvals1[0] - ilat * yinc;
          for (long ilat = latMinIdx; ilat < latMaxIdx; ilat++) *pyvals2++ = yvals1[ilat];
          for (long ilat = 1; ilat <= halo.north; ilat++)       *pyvals2++ = yvals1[nlat1 - 1] + ilat * yinc;
          // clang-format on
        }

      gridDefXvals(gridID2, xvals2.data());
      gridDefYvals(gridID2, yvals2.data());
    }

  if (isCircularGrid && gridHasBounds(gridID1))
    {
      auto nv = isCurvilinearGrid ? 4 : 2;
      gridDefNvertex(gridID2, nv);

      Varray<double> xbounds1(isCurvilinearGrid ? nv * nlon1 * nlat1 : nv * nlon1);
      Varray<double> ybounds1(isCurvilinearGrid ? nv * nlon1 * nlat1 : nv * nlat1);
      Varray<double> xbounds2(isCurvilinearGrid ? nv * nlon2 * nlat2 : nv * nlon2);
      Varray<double> ybounds2(isCurvilinearGrid ? nv * nlon2 * nlat2 : nv * nlat2);

      double *pxbounds2 = xbounds2.data();
      double *pybounds2 = ybounds2.data();

      gridInqXbounds(gridID1, xbounds1.data());
      gridInqYbounds(gridID1, ybounds1.data());

      if (isCurvilinearGrid)
        {
          if (isCircularGrid)
            {
              for (long ilat = 0; ilat < nlat1; ilat++)
                {
                  const auto pxbounds1 = &xbounds1[nv * ilat * nlon1];
                  const auto pybounds1 = &ybounds1[nv * ilat * nlon1];
                  for (long ilon = nv * (nlon1 - halo.east); ilon < nv * nlon1; ilon++) *pxbounds2++ = pxbounds1[ilon];
                  for (long ilon = nv * (nlon1 - halo.east); ilon < nv * nlon1; ilon++) *pybounds2++ = pybounds1[ilon];

                  for (long ilon = nv * lonMinIdx; ilon < nv * lonMaxIdx; ilon++) *pxbounds2++ = pxbounds1[ilon];
                  for (long ilon = nv * lonMinIdx; ilon < nv * lonMaxIdx; ilon++) *pybounds2++ = pybounds1[ilon];

                  for (long ilon = 0; ilon < nv * halo.west; ilon++) *pxbounds2++ = pxbounds1[ilon];
                  for (long ilon = 0; ilon < nv * halo.west; ilon++) *pybounds2++ = pybounds1[ilon];
                }
            }
        }
      else
        {
          if (isCircularGrid)
            {
              // clang-format off
              for (long ilon = nv * (nlon1 - halo.east); ilon < nv * nlon1; ilon++) *pxbounds2++ = xbounds1[ilon] - cpi2;
              for (long ilon = nv * lonMinIdx; ilon < nv * lonMaxIdx; ilon++)       *pxbounds2++ = xbounds1[ilon];
              for (long ilon = 0; ilon < nv * halo.west; ilon++)                    *pxbounds2++ = xbounds1[ilon] + cpi2;
              // clang-format on
            }

          for (long ilat = 0; ilat < nv * nlat2; ++ilat) ybounds2[ilat] = ybounds1[ilat];
        }

      gridDefXbounds(gridID2, xbounds2.data());
      gridDefYbounds(gridID2, ybounds2.data());
    }

  return gridID2;
}

static Halo
get_parameter(void)
{
  Halo halo;

  auto pargc = cdo_operator_argc();
  if (pargc)
    {
      auto &pargv = cdo_get_oper_argv();

      KVList kvlist;
      kvlist.name = cdo_module_name();
      if (kvlist.parse_arguments(pargc, pargv) != 0) cdo_abort("Parse error!");
      if (Options::cdoVerbose) kvlist.print();

      for (const auto &kv : kvlist)
        {
          const auto &key = kv.key;
          if (kv.nvalues > 1) cdo_abort("Too many values for parameter key >%s<!", key);
          if (kv.nvalues < 1) cdo_abort("Missing value for parameter key >%s<!", key);
          const auto &value = kv.values[0];

          // clang-format off
          if      (key == "east")   halo.east = parameter_to_long(value);
          else if (key == "west")   halo.west = parameter_to_long(value);
          else if (key == "south")  halo.south = parameter_to_long(value);
          else if (key == "north")  halo.north = parameter_to_long(value);
          else if (key == "value")  halo.value = parameter_to_double(value);
          else cdo_abort("Invalid parameter key >%s<!", key);
          // clang-format on
        }
    }

  return halo;
}

static int
gen_index_grid(int gridID1, Halo &halo)
{
  if (cdo_operator_argc() == 2 && string_is_int(cdo_operator_argv(0)) && string_is_int(cdo_operator_argv(1)))
    {
      halo.east = parameter_to_long(cdo_operator_argv(0));
      halo.west = parameter_to_long(cdo_operator_argv(1));
    }
  else
    {
      halo = get_parameter();
    }

  auto isCircularGrid = gridIsCircular(gridID1);
  long nlon = gridInqXsize(gridID1);
  long nlat = gridInqYsize(gridID1);

  if (isCircularGrid && halo.east > nlon) cdo_abort("east halo out of range (max=%ld).", nlon);
  if (isCircularGrid && halo.west > nlon) cdo_abort("west halo out of range (max=%ld).", nlon);
  if (halo.east < 0 && -halo.east > nlon - 1) cdo_abort("negative east halo out of range (max=%ld).", -nlon + 1);
  if (halo.west < 0 && -halo.west > nlon - 1) cdo_abort("negative west halo out of range (max=%ld).", -nlon + 1);
  if (halo.south < 0 && -halo.south > nlat - 1) cdo_abort("negative south halo out of range (max=%ld).", -nlat + 1);
  if (halo.north < 0 && -halo.north > nlat - 1) cdo_abort("negative north halo out of range (max=%ld).", -nlat + 1);
  if (halo.east + halo.west < -nlon + 1) cdo_abort("sum of negative east and west halo out of range (max=%ld).", -nlon + 1);
  if (halo.south + halo.north < -nlat + 1) cdo_abort("sum of negative south and north halo out of range (max=%ld).", -nlat + 1);

  auto gridID2 = gen_regular_grid(gridID1, halo);

  return gridID2;
}

static bool
regular_halo(const Varray<double> &array1, int gridID1, Varray<double> &array2, Halo &halo, double missval)
{
  auto recalcNumMiss = false;

  long nlon1 = gridInqXsize(gridID1);
  long nlat1 = gridInqYsize(gridID1);
  long nlon2 = nlon1 + halo.east + halo.west;

  auto isCircularGrid = gridIsCircular(gridID1);
  auto fillValueDefined = is_not_equal(halo.value, UndefValue);
  auto fillValue = fillValueDefined ? halo.value : missval;
  auto useFillValue = (!isCircularGrid || fillValueDefined);
  if (!fillValueDefined && (halo.east > 0 || halo.west > 0 || halo.south > 0 || halo.north > 0)) recalcNumMiss = true;

  long lonMinIdx = 0;
  long lonMaxIdx = nlon1;
  if (halo.east < 0) lonMinIdx = -halo.east;
  if (halo.west < 0) lonMaxIdx += halo.west;
  long latMinIdx = 0;
  long latMaxIdx = nlat1;
  if (halo.south < 0) latMinIdx = -halo.south;
  if (halo.north < 0) latMaxIdx += halo.north;

  auto parray2 = array2.data();

  for (long i = 0; i < nlon2 * halo.south; i++) *parray2++ = fillValue;

  for (long ilat = latMinIdx; ilat < latMaxIdx; ilat++)
    {
      const auto parray1 = &array1[ilat * nlon1];
      if (useFillValue)
        {
          // clang-format off
          for (long ilon = 0; ilon < halo.east; ilon++)             *parray2++ = fillValue;
          for (long ilon = lonMinIdx; ilon < lonMaxIdx; ilon++)     *parray2++ = parray1[ilon];
          for (long ilon = 0; ilon < halo.west; ilon++)             *parray2++ = fillValue;
          // clang-format on
        }
      else
        {
          // clang-format off
          for (long ilon = nlon1 - halo.east; ilon < nlon1; ilon++) *parray2++ = parray1[ilon];
          for (long ilon = lonMinIdx; ilon < lonMaxIdx; ilon++)     *parray2++ = parray1[ilon];
          for (long ilon = 0; ilon < halo.west; ilon++)             *parray2++ = parray1[ilon];
          // clang-format on
        }
    }

  for (long i = 0; i < nlon2 * halo.north; i++) *parray2++ = fillValue;

  return recalcNumMiss;
}

static void
tripolar_halo(const Varray<double> &array1, int gridID1, Varray<double> &array2)
{
  auto nlon = gridInqXsize(gridID1);
  auto nlat = gridInqYsize(gridID1);

  for (size_t ilat = 0; ilat < nlat; ilat++)
    for (size_t ilon = 0; ilon < nlon; ilon++) array2[(ilat + 2) * nlon + ilon] = array1[ilat * nlon + ilon];

  for (size_t ilon = 0; ilon < nlon; ilon++)
    {
      size_t ilonr = nlon - ilon - 1;
      array2[1 * nlon + ilon] = array2[2 * nlon + ilonr];  // syncronise line 2 with line 3
      array2[0 * nlon + ilon] = array2[3 * nlon + ilonr];  // syncronise line 1 with line 4
    }
}

static int
get_gridID(int vlistID1)
{
  int gridID1 = -1;

  auto ngrids = vlistNgrids(vlistID1);
  int ndiffgrids = 0;
  int index;
  for (index = 1; index < ngrids; ++index)
    if (vlistGrid(vlistID1, 0) != vlistGrid(vlistID1, index)) ndiffgrids++;

  for (index = 0; index < ngrids; ++index)
    {
      gridID1 = vlistGrid(vlistID1, index);
      auto gridtype = gridInqType(gridID1);
      if (gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN) break;
      if (gridtype == GRID_CURVILINEAR) break;
      if (gridtype == GRID_GENERIC && gridInqXsize(gridID1) > 0 && gridInqYsize(gridID1) > 0) break;
    }

  if (index == ngrids) cdo_abort("No regular 2d lon/lat grid found!");
  if (ndiffgrids > 0) cdo_abort("Too many different grids!");

  return gridID1;
}

class SetHalo
{
  int SETHALO;
  int operatorID;

  Halo halo;

  CdoStreamID streamID1;
  CdoStreamID streamID2;

  int taxisID1;
  int taxisID2;

  int gridID1 = -1;

  VarList varList;

  std::vector<bool> vars;
  Varray<double> array1;
  Varray<double> array2;

public:
  void
  init(void *process)
  {
    cdo_initialize(process);

    // clang-format off
    SETHALO = cdo_operator_add("sethalo", 0, 0, nullptr);
              cdo_operator_add("tpnhalo", 0, 0, nullptr);
    // clang-format on

    operatorID = cdo_operator_id();

    streamID1 = cdo_open_read(0);
    auto vlistID1 = cdo_stream_inq_vlist(streamID1);

    gridID1 = get_gridID(vlistID1);

    auto gridID2 = (operatorID == SETHALO) ? gen_index_grid(gridID1, halo) : gen_tripolar_grid(gridID1);

    auto vlistID2 = vlistDuplicate(vlistID1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    auto ngrids = vlistNgrids(vlistID1);
    for (int index = 0; index < ngrids; ++index)
      {
        if (gridID1 == vlistGrid(vlistID1, index))
          {
            vlistChangeGridIndex(vlistID2, index, gridID2);
            break;
          }
      }

    varListInit(varList, vlistID1);

    auto nvars = vlistNvars(vlistID1);
    vars = std::vector<bool>(nvars);
    for (int varID = 0; varID < nvars; ++varID) vars[varID] = (gridID1 == varList[varID].gridID);

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);

    array1 = Varray<double>(gridInqSize(gridID1));
    array2 = Varray<double>(gridInqSize(gridID2));
  }

  void
  run()
  {
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

            if (vars[varID])
              {
                size_t nmiss;
                cdo_read_record(streamID1, array1.data(), &nmiss);

                auto recalcNumMiss = false;
                auto missval = varList[varID].missval;
                if (operatorID == SETHALO)
                  recalcNumMiss = regular_halo(array1, gridID1, array2, halo, missval);
                else
                  tripolar_halo(array1, gridID1, array2);

                if (nmiss || recalcNumMiss) nmiss = varray_num_mv(array2.size(), array2, missval);

                cdo_def_record(streamID2, varID, levelID);
                cdo_write_record(streamID2, array2.data(), nmiss);
              }
          }

        tsID++;
      }
  }

  void
  close()
  {
    cdo_stream_close(streamID2);
    cdo_stream_close(streamID1);

    cdo_finish();
  }
};

void *
Sethalo(void *process)
{
  SetHalo sethalo;
  sethalo.init(process);
  sethalo.run();
  sethalo.close();

  return nullptr;
}
