/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include "process_int.h"
#include "param_conversion.h"
#include "pmlist.h"
#include "grid_point_search.h"
#include "mpim_grid.h"
#include "region.h"

int gengridcell(const int gridID1, const size_t gridsize2, const std::vector<long> &cellidx);
void window_cell(const Field &field1, Field &field2, const std::vector<long> &cellidx);
double radiusDegToKm(const double radiusInDeg);

struct CirclePoint
{
  double radius = 1.0;
  double lon = 0.0;
  double lat = 0.0;
  size_t maxpoints = SIZE_MAX;
};

struct RegionInfo
{
  std::vector<long> cellidx;
  long nvals = 0;
  int gridtype = -1;
  int gridID1 = -1, gridID2 = -1;
};

static inline bool
is_point_inside(double xval, double yval, double xi, double xj, double yi, double yj)
{
  return (((yval >= yi && yval < yj) || (yval > yj && yval <= yi)) && (xval < ((xj - xi) * (yval - yi) / (yj - yi) + xi)));
}

static bool
point_is_inside(double xval, double yval, size_t n, const double *xcoords, const double *ycoords)
{
  auto c = false;

  for (size_t i = 0, j = n - 1; i < n; j = i++)
    {
      if (is_point_inside(xval, yval, xcoords[i], xcoords[j], ycoords[i], ycoords[j])) c = !c;
    }

  return c;
}

static bool
point_is_inside(double xval, double yval, double xmin, double xmax, const double *xcoords, const double *ycoords, size_t nofcoords)
{
  auto c = false;

  // clang-format off
  if      (xval >= xmin && xval <= xmax)
    c = point_is_inside(xval,         yval, nofcoords, xcoords, ycoords);
  else if (xval > 180.0 && xval - 360.0 >= xmin && xval - 360.0 <= xmax)
    c = point_is_inside(xval - 360.0, yval, nofcoords, xcoords, ycoords);
  else if (xval <   0.0 && xval + 360.0 >= xmin && xval + 360.0 <= xmax)
    c = point_is_inside(xval + 360.0, yval, nofcoords, xcoords, ycoords);
  // clang-format on

  return c;
}

static void
sel_region_cell(std::vector<char> &mask, size_t gridsize, const Varray<double> &xvals, const Varray<double> &yvals,
                const double *xcoords, const double *ycoords, size_t segmentSize, std::vector<long> &cellidx)
{
  auto xmm = varray_min_max(segmentSize, xcoords);
  auto ymm = varray_min_max(segmentSize, ycoords);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) default(shared)
#endif
  for (size_t i = 0; i < gridsize; ++i)
    {
      if (mask[i]) continue;

      auto yval = yvals[i];
      if (yval > ymm.min && yval < ymm.max)
        {
          if (point_is_inside(xvals[i], yval, xmm.min, xmm.max, xcoords, ycoords, segmentSize)) mask[i] = true;
        }
    }

  for (size_t i = 0; i < gridsize; ++i)
    {
      if (mask[i]) cellidx.push_back(i);
    }
}

static int
generate_region_grid(int gridID1, long &gridsize2, std::vector<long> &cellidx, int numFiles)
{
  auto gridID0 = gridID1;

  gridID1 = generate_full_grid(gridID1);
  if (!gridHasCoordinates(gridID1)) cdo_abort("Cell center coordinates missing!");

  auto gridsize = gridInqSize(gridID1);
  Varray<double> xvals(gridsize), yvals(gridsize);

  gridInqXvals(gridID1, xvals.data());
  gridInqYvals(gridID1, yvals.data());

  // Convert lat/lon units if required
  cdo_grid_to_degree(gridID1, CDI_XAXIS, xvals.size(), xvals.data(), "grid center lon");
  cdo_grid_to_degree(gridID1, CDI_YAXIS, yvals.size(), yvals.data(), "grid center lat");

  std::vector<char> mask(gridsize, false);
  for (int i = 0; i < numFiles; ++i)
    {
      Regions regions;
      auto param = cdo_operator_argv(i).c_str();
      if (strncmp(param, "dcw:", 4) == 0)
        read_regions_from_dcw(param + 4, regions);
      else
        read_regions_from_file(param, regions);

      for (size_t k = 0; k < regions.numSegments; ++k)
        {
          auto segmentSize = regions.segmentSize[k];
          if (segmentSize < 3) continue;
          auto offset = regions.segmentOffset[k];
          auto xcoords = &regions.x[offset];
          auto ycoords = &regions.y[offset];
          sel_region_cell(mask, gridsize, xvals, yvals, xcoords, ycoords, segmentSize, cellidx);
        }
    }

  gridsize2 = cellidx.size();
  if (gridsize2 == 0) cdo_abort("No grid points found!");

  auto gridID2 = gridsize2 ? gengridcell(gridID1, gridsize2, cellidx) : CDI_UNDEFID;

  if (gridID0 != gridID1) gridDestroy(gridID1);

  return gridID2;
}

static int
generate_circle_grid(int gridID1, long &gridsize2, std::vector<long> &cellidx, const CirclePoint &cpoint)
{
  auto gridID0 = gridID1;

  gridID1 = generate_full_grid(gridID1);
  if (!gridHasCoordinates(gridID1)) cdo_abort("Cell center coordinates missing!");

  {
    auto gridsize1 = gridInqSize(gridID1);

    Varray<double> xvals(gridsize1), yvals(gridsize1);
    gridInqXvals(gridID1, xvals.data());
    gridInqYvals(gridID1, yvals.data());

    // Convert lat/lon units if required
    cdo_grid_to_radian(gridID1, CDI_XAXIS, gridsize1, xvals.data(), "grid center lon");
    cdo_grid_to_radian(gridID1, CDI_YAXIS, gridsize1, yvals.data(), "grid center lat");

    GridPointSearch gps;
    // grid_point_search_create(gps, xvals, yvals, PointSearchMethod::spherepart);
    grid_point_search_create(gps, xvals, yvals);

    grid_point_search_set_arc_radius(gps, cpoint.radius);

    auto numNeighbors = cpoint.maxpoints;
    if (numNeighbors > gridsize1) numNeighbors = gridsize1;

    knnWeightsType knnWeights(numNeighbors);
    grid_search_point_smooth(gps, cpoint.lon, cpoint.lat, knnWeights);

    auto nvals = knnWeights.m_numNeighbors;
    cellidx.resize(nvals);

    for (size_t i = 0; i < nvals; ++i) cellidx[i] = knnWeights.m_addr[i];

    grid_point_search_delete(gps);

    if (nvals == 0) cdo_abort("No grid points found!");

    gridsize2 = nvals;
  }

  auto gridID2 = gengridcell(gridID1, gridsize2, cellidx);

  if (gridID0 != gridID1) gridDestroy(gridID1);

  return gridID2;
}

static void
selcircle_get_parameter(CirclePoint &cpoint)
{
  auto pargc = cdo_operator_argc();
  if (pargc)
    {
      const auto &pargv = cdo_get_oper_argv();

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
          if      (key == "maxpoints") cpoint.maxpoints = parameter_to_size_t(value);
          else if (key == "lon")       cpoint.lon = parameter_to_double(value);
          else if (key == "lat")       cpoint.lat = parameter_to_double(value);
          else if (key == "radius")    cpoint.radius = radius_str_to_deg(value);
          else cdo_abort("Invalid parameter key >%s<!", key);
          // clang-format on
        }
    }
}

void *
Selregion(void *process)
{
  cdo_initialize(process);

  // clang-format off
  auto SELREGION  = cdo_operator_add("selregion",   0,   0, nullptr);
  auto SELCIRCLE  = cdo_operator_add("selcircle",   0,   0, nullptr);
  // clang-format on

  auto operatorID = cdo_operator_id();

  auto streamID1 = cdo_open_read(0);

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto vlistID2 = vlistDuplicate(vlistID1);

  auto taxisID1 = vlistInqTaxis(vlistID1);
  auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  auto nvars = vlistNvars(vlistID1);
  std::vector<bool> varIDs(nvars, false);

  auto ngrids = vlistNgrids(vlistID1);
  std::vector<RegionInfo> regions(ngrids);

  int numFiles = 0;
  CirclePoint cpoint;

  if (operatorID == SELREGION)
    {
      numFiles = cdo_operator_argc();
      if (numFiles == 0) cdo_abort("Region parameter missing!");
    }
  else if (operatorID == SELCIRCLE)
    {
      selcircle_get_parameter(cpoint);
      if (cpoint.radius < 0.0 || cpoint.radius > 180.0) cdo_abort("%s=%g out of bounds (0-180 deg)!", "radius", cpoint.radius);

      auto gridsizemax = vlistGridsizeMax(vlistID1);
      if (gridsizemax < cpoint.maxpoints) cpoint.maxpoints = gridsizemax;
      if (Options::cdoVerbose && operatorID == SELCIRCLE)
        cdo_print("lon = %g, lat = %g, radius = %gdeg(%gkm)", cpoint.lon, cpoint.lat, cpoint.radius, radiusDegToKm(cpoint.radius));

      cpoint.radius *= DEG2RAD;
      cpoint.lon *= DEG2RAD;
      cpoint.lat *= DEG2RAD;
    }

  for (int index = 0; index < ngrids; ++index)
    {
      auto &region = regions[index];
      auto gridID1 = vlistGrid(vlistID1, index);
      auto gridtype = gridInqType(gridID1);
      if (is_point_grid(gridID1))
        {
          auto gridsize = gridInqSize(gridID1);
          if (gridsize == 1) continue;

          region.cellidx.reserve(gridsize);

          int gridID2 = CDI_UNDEFID;
          if (operatorID == SELREGION)
            gridID2 = generate_region_grid(gridID1, region.nvals, region.cellidx, numFiles);
          else if (operatorID == SELCIRCLE)
            gridID2 = generate_circle_grid(gridID1, region.nvals, region.cellidx, cpoint);

          region.cellidx.shrink_to_fit();

          if (gridID2 != CDI_UNDEFID)
            {
              region.gridtype = gridtype;
              region.gridID1 = gridID1;
              region.gridID2 = gridID2;

              vlistChangeGridIndex(vlistID2, index, gridID2);

              for (int varID = 0; varID < nvars; ++varID)
                if (gridID1 == vlistInqVarGrid(vlistID1, varID)) varIDs[varID] = true;
            }
        }
      else { cdo_abort("Unsupported grid type: %s", gridNamePtr(gridtype)); }
    }

  VarList varList1, varList2;
  varListInit(varList1, vlistID1);
  varListInit(varList2, vlistID2);

  Field field1, field2;

  auto streamID2 = cdo_open_write(1);
  cdo_def_vlist(streamID2, vlistID2);

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
          const auto &var = varList1[varID];
          field1.init(var);
          cdo_read_record(streamID1, field1);

          cdo_def_record(streamID2, varID, levelID);

          if (varIDs[varID])
            {
              auto gridID1 = var.gridID;
              int index;
              for (index = 0; index < ngrids; ++index)
                if (gridID1 == regions[index].gridID1) break;
              if (index == ngrids) cdo_abort("Internal problem, grid not found!");

              field2.init(varList2[varID]);
              window_cell(field1, field2, regions[index].cellidx);

              if (field1.nmiss) field_num_mv(field2);

              cdo_write_record(streamID2, field2);
            }
          else { cdo_write_record(streamID2, field1); }
        }

      tsID++;
    }

  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  cdo_finish();

  return nullptr;
}
