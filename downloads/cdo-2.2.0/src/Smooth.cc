/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Smooth        smooth          Smooth grid points
      Smooth        smooth9         9 point smoothing
*/

#include <atomic>
#include <sstream>

#include <cdi.h>

#include "process_int.h"
#include "param_conversion.h"
#include "cdo_wtime.h"
#include <mpim_grid.h>
#include "constants.h"  // planet radius
#include "pmlist.h"
#include "cdo_options.h"
#include "progress.h"
#include "cimdOmp.h"
#include "grid_point_search.h"

enum class CurveForm
{
  AVG = 0,
  LINEAR,
};

static const char *Form[] = { "avg", "linear" };

struct SmoothPoint
{
  double arc_radius = 0.0;
  double radius = 1.0;
  double weight0 = 0.25;
  double weightR = 0.25;
  size_t maxpoints = SIZE_MAX;
  CurveForm form = CurveForm::LINEAR;
};

template <typename T>
static size_t
smooth(int gridID, T missval, const Varray<T> &array1, Varray<T> &array2, const SmoothPoint &spoint)
{
  auto gridID0 = gridID;
  auto gridsize = gridInqSize(gridID);
  auto numNeighbors = spoint.maxpoints;
  if (numNeighbors > gridsize) numNeighbors = gridsize;

  Varray<uint8_t> mask(gridsize);
  for (size_t i = 0; i < gridsize; ++i) mask[i] = !dbl_is_equal(array1[i], missval);

  gridID = generate_full_point_grid(gridID);
  if (!gridHasCoordinates(gridID)) cdo_abort("Cell center coordinates missing!");

  Varray<double> xvals(gridsize), yvals(gridsize);
  gridInqXvals(gridID, xvals.data());
  gridInqYvals(gridID, yvals.data());

  // Convert lat/lon units if required
  cdo_grid_to_radian(gridID, CDI_XAXIS, gridsize, xvals.data(), "grid center lon");
  cdo_grid_to_radian(gridID, CDI_YAXIS, gridsize, yvals.data(), "grid center lat");

  std::vector<knnWeightsType> knnWeights;
  for (int i = 0; i < Threading::ompNumThreads; ++i) knnWeights.push_back(knnWeightsType(numNeighbors));

  auto start = Options::cdoVerbose ? cdo_get_wtime() : 0.0;

  bool xIsCyclic = false;
  size_t dims[2] = { gridsize, 0 };
  GridPointSearch gps;
  grid_point_search_create(gps, xIsCyclic, dims, gridsize, xvals, yvals);

  if (spoint.arc_radius > 0.0)
    grid_point_search_set_arc_radius(gps, spoint.arc_radius);
  else
    grid_point_search_set_chord_radius(gps, spoint.radius);

  if (Options::cdoVerbose) cdo_print("Point search created: %.2f seconds (%zu points)", cdo_get_wtime() - start, gridsize);

  if (Options::cdoVerbose) progress::init();

  start = Options::cdoVerbose ? cdo_get_wtime() : 0.0;

  size_t naddsMin = gridsize, naddsMax = 0;
  std::atomic<size_t> atomicCount{ 0 }, atomicSum{ 0 }, atomicNumMiss{ 0 };

#ifdef HAVE_OPENMP4
#pragma omp parallel for default(shared) schedule(dynamic) reduction(min : naddsMin) reduction(max : naddsMax)
#endif
  for (size_t i = 0; i < gridsize; ++i)
    {
      const auto ompthID = cdo_omp_get_thread_num();

      atomicCount++;
      if (Options::cdoVerbose && cdo_omp_get_thread_num() == 0) progress::update(0, 1, (double) atomicCount / gridsize);

      grid_search_point_smooth(gps, xvals[i], yvals[i], knnWeights[ompthID]);

      // Compute weights based on inverse distance if mask is false, eliminate those points
      const auto nadds = knnWeights[ompthID].computeWeights(mask, spoint.radius, spoint.weight0, spoint.weightR);
      naddsMin = std::min(naddsMin, nadds);
      naddsMax = std::max(naddsMax, nadds);

      array2[i] = nadds ? knnWeights[ompthID].arrayWeightsSum(array1) : missval;
      atomicSum += nadds;
      if (nadds == 0) atomicNumMiss++;
    }

  progress::update(0, 1, 1);

  size_t nmissx = atomicNumMiss;
  size_t numPoints = atomicSum;

  if (Options::cdoVerbose) cdo_print("Point search nearest: %.2f seconds (%zu points)", cdo_get_wtime() - start, numPoints);
  if (Options::cdoVerbose) cdo_print("Min/Max points found: %zu/%zu", naddsMin, naddsMax);

  grid_point_search_delete(gps);

  if (gridID0 != gridID) gridDestroy(gridID);

  return nmissx;
}

static void
smooth(const Field &field1, Field &field2, const SmoothPoint &spoint)
{
  if (field1.memType != field2.memType) cdo_abort("Interal error, memType of field1 and field2 differ!");

  if (field1.memType == MemType::Float)
    field2.nmiss = smooth(field1.grid, (float) field1.missval, field1.vec_f, field2.vec_f, spoint);
  else
    field2.nmiss = smooth(field1.grid, field1.missval, field1.vec_d, field2.vec_d, spoint);
}

template <typename T>
static inline void
smooth9_sum(size_t ij, const std::vector<uint8_t> &mask, double sfac, const Varray<T> &array, double &avg, double &divavg)
{
  if (mask[ij])
    {
      avg += sfac * array[ij];
      divavg += sfac;
    }
}

template <typename T>
static size_t
smooth9(int gridID, T missval, const Varray<T> &array1, Varray<T> &array2)
{
  const auto gridsize = gridInqSize(gridID);
  const auto nlon = gridInqXsize(gridID);
  const auto nlat = gridInqYsize(gridID);
  const auto gridIsCyclic = gridIsCircular(gridID);

  std::vector<uint8_t> mask(gridsize);

  for (size_t i = 0; i < gridsize; ++i) mask[i] = !dbl_is_equal(missval, array1[i]);

  size_t nmiss = 0;
  for (size_t i = 0; i < nlat; ++i)
    {
      for (size_t j = 0; j < nlon; ++j)
        {
          double avg = 0, divavg = 0;

          if ((i == 0) || (j == 0) || (i == (nlat - 1)) || (j == (nlon - 1)))
            {
              const auto ij = j + nlon * i;
              if (mask[ij])
                {
                  avg += array1[ij];
                  divavg += 1;
                  // upper left corner
                  if ((i != 0) && (j != 0))
                    smooth9_sum(((i - 1) * nlon) + j - 1, mask, 0.3, array1, avg, divavg);
                  else if (i != 0 && gridIsCyclic)
                    smooth9_sum((i - 1) * nlon + j - 1 + nlon, mask, 0.3, array1, avg, divavg);

                  // upper cell
                  if (i != 0) smooth9_sum(((i - 1) * nlon) + j, mask, 0.5, array1, avg, divavg);

                  // upper right corner
                  if ((i != 0) && (j != (nlon - 1)))
                    smooth9_sum(((i - 1) * nlon) + j + 1, mask, 0.3, array1, avg, divavg);
                  else if ((i != 0) && gridIsCyclic)
                    smooth9_sum((i - 1) * nlon + j + 1 - nlon, mask, 0.3, array1, avg, divavg);

                  // left cell
                  if (j != 0)
                    smooth9_sum(i * nlon + j - 1, mask, 0.5, array1, avg, divavg);
                  else if (gridIsCyclic)
                    smooth9_sum(i * nlon - 1 + nlon, mask, 0.5, array1, avg, divavg);

                  // right cell
                  if (j != (nlon - 1))
                    smooth9_sum((i * nlon) + j + 1, mask, 0.5, array1, avg, divavg);
                  else if (gridIsCyclic)
                    smooth9_sum(i * nlon + j + 1 - nlon, mask, 0.5, array1, avg, divavg);

                  // lower left corner
                  if ((i != (nlat - 1)) && (j != 0))
                    smooth9_sum(((i + 1) * nlon + j - 1), mask, 0.3, array1, avg, divavg);
                  else if ((i != (nlat - 1)) && gridIsCyclic)
                    smooth9_sum((i + 1) * nlon - 1 + nlon, mask, 0.3, array1, avg, divavg);

                  // lower cell
                  if (i != (nlat - 1)) smooth9_sum(((i + 1) * nlon) + j, mask, 0.5, array1, avg, divavg);

                  // lower right corner
                  if ((i != (nlat - 1)) && (j != (nlon - 1)))
                    smooth9_sum(((i + 1) * nlon) + j + 1, mask, 0.3, array1, avg, divavg);
                  else if ((i != (nlat - 1)) && gridIsCyclic)
                    smooth9_sum(((i + 1) * nlon) + j + 1 - nlon, mask, 0.3, array1, avg, divavg);
                }
            }
          else if (mask[j + nlon * i])
            {
              avg += array1[j + nlon * i];
              divavg += 1;

              smooth9_sum(((i - 1) * nlon) + j - 1, mask, 0.3, array1, avg, divavg);
              smooth9_sum(((i - 1) * nlon) + j, mask, 0.5, array1, avg, divavg);
              smooth9_sum(((i - 1) * nlon) + j + 1, mask, 0.3, array1, avg, divavg);
              smooth9_sum(((i) *nlon) + j - 1, mask, 0.5, array1, avg, divavg);
              smooth9_sum((i * nlon) + j + 1, mask, 0.5, array1, avg, divavg);
              smooth9_sum(((i + 1) * nlon + j - 1), mask, 0.3, array1, avg, divavg);
              smooth9_sum(((i + 1) * nlon) + j, mask, 0.5, array1, avg, divavg);
              smooth9_sum(((i + 1) * nlon) + j + 1, mask, 0.3, array1, avg, divavg);
            }

          if (std::fabs(divavg) > 0) { array2[i * nlon + j] = avg / divavg; }
          else
            {
              array2[i * nlon + j] = missval;
              nmiss++;
            }
        }
    }

  return nmiss;
}

static void
smooth9(const Field &field1, Field &field2)
{
  if (field1.memType != field2.memType) cdo_abort("Interal error, memType of field1 and field2 differ!");

  if (field1.memType == MemType::Float)
    field2.nmiss = smooth9(field1.grid, (float) field1.missval, field1.vec_f, field2.vec_f);
  else
    field2.nmiss = smooth9(field1.grid, field1.missval, field1.vec_d, field2.vec_d);
}

double
radiusDegToKm(const double radiusInDeg)
{
  return radiusInDeg * (2.0 * PlanetRadius * M_PI) / (360.0 * 1000.0);
}

static CurveForm
convert_curveform(const std::string &formstr)
{
  CurveForm form = CurveForm::LINEAR;

  if (formstr == "linear")
    form = CurveForm::LINEAR;
  else
    cdo_abort("form=%s unsupported!", formstr);

  return form;
}

static void
get_parameter(int &xnsmooth, SmoothPoint &spoint)
{
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
          if      (key == "nsmooth")    xnsmooth = parameter_to_int(value);
          else if (key == "maxpoints")  spoint.maxpoints = parameter_to_size_t(value);
          else if (key == "weight0")    spoint.weight0 = parameter_to_double(value);
          else if (key == "weightR")    spoint.weightR = parameter_to_double(value);
          else if (key == "radius")     spoint.radius = radius_str_to_deg(value);
          else if (key == "arc_radius") spoint.arc_radius = radius_str_to_deg(value);
          else if (key == "form")       spoint.form = convert_curveform(value);
          else cdo_abort("Invalid parameter key >%s<!", key);
          // clang-format on
        }
    }
}

static void
print_parameter(int xnsmooth, const SmoothPoint &sp)
{
  std::stringstream outbuffer;

  outbuffer << "nsmooth=" << xnsmooth;
  outbuffer << ", maxpoints=" << sp.maxpoints;

  if (sp.arc_radius > 0.0)
    outbuffer << ", arc_radius=" << sp.arc_radius << "deg(" << radiusDegToKm(sp.arc_radius) << "km)";
  else
    outbuffer << ", radius=" << sp.radius << "deg(" << radiusDegToKm(sp.radius) << "km)";

  outbuffer << ", form=" << Form[(int) sp.form];
  if (sp.form == CurveForm::LINEAR) outbuffer << ", weight0=" << sp.weight0 << ", weightR=" << sp.weightR;

  cdo_print("%s", outbuffer.str());
}

static void
check_radius_range(double radius, const char *name)
{
  if (radius < 0.0 || radius > 180.0) cdo_abort("%s=%g out of bounds (0-180 deg)!", name, radius);
}

void *
Smooth(void *process)
{
  int xnsmooth = 1;

  cdo_initialize(process);

  // clang-format off
  const auto SMOOTH  = cdo_operator_add("smooth",   0,   0, nullptr);
  const auto SMOOTH9 = cdo_operator_add("smooth9",  0,   0, nullptr);
  // clang-format on

  const auto operatorID = cdo_operator_id();

  SmoothPoint spoint;
  if (operatorID == SMOOTH) get_parameter(xnsmooth, spoint);

  check_radius_range(spoint.radius, "radius");
  check_radius_range(spoint.arc_radius, "arc_radius");

  const auto streamID1 = cdo_open_read(0);

  const auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  const auto vlistID2 = vlistDuplicate(vlistID1);

  const auto taxisID1 = vlistInqTaxis(vlistID1);
  const auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  VarList varList1;
  varListInit(varList1, vlistID1);

  const auto nvars = vlistNvars(vlistID1);
  std::vector<bool> varIDs(nvars, false);

  for (int varID = 0; varID < nvars; ++varID)
    {
      const auto &var = varList1[varID];
      const auto gridID = var.gridID;
      const auto gridtype = gridInqType(gridID);
      if (gridtype == GRID_GAUSSIAN || gridtype == GRID_LONLAT || gridtype == GRID_CURVILINEAR || gridtype == GRID_PROJECTION
          || (operatorID == SMOOTH9 && gridtype == GRID_GENERIC && gridInqXsize(gridID) > 0 && gridInqYsize(gridID) > 0))
        {
          varIDs[varID] = true;
        }
      else if (operatorID == SMOOTH && gridtype == GRID_UNSTRUCTURED) { varIDs[varID] = true; }
      else { cdo_warning("Unsupported grid for variable %s", var.name); }
    }

  const auto gridsizemax = vlistGridsizeMax(vlistID1);
  if (gridsizemax < spoint.maxpoints) spoint.maxpoints = gridsizemax;
  if (Options::cdoVerbose && operatorID == SMOOTH) print_parameter(xnsmooth, spoint);

  spoint.radius *= DEG2RAD;
  spoint.arc_radius *= DEG2RAD;

  Field field1, field2;

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
          const auto &var = varList1[varID];
          field1.init(var);
          field2.init(var);
          cdo_read_record(streamID1, field1);

          if (varIDs[varID])
            {
              for (int i = 0; i < xnsmooth; ++i)
                {
                  if (operatorID == SMOOTH)
                    smooth(field1, field2, spoint);
                  else if (operatorID == SMOOTH9)
                    smooth9(field1, field2);

                  field_copy(field2, field1);
                }
            }

          cdo_def_record(streamID2, varID, levelID);
          cdo_write_record(streamID2, field1);
        }

      tsID++;
    }

  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  cdo_finish();

  return nullptr;
}
