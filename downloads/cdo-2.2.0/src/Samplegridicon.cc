/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <algorithm>  // sort

#include <cdi.h>

#include "cdo_options.h"
#include "dmemory.h"
#include "process_int.h"
#include "cdi_lockedIO.h"
#include <mpim_grid.h>
#include "grid_point_search.h"
#include "verifygrid.h"
#include "field_functions.h"

constexpr int MAX_CHILDS = 9;

struct CellIndex
{
  long ncells;
  long *neighbor;  // neighbor cell index
  long *parent;    // parent cell index
  long *child;     // child cell index
  const char *filename;
};

static void
copy_data_to_index(long ncells, const Varray<double> &data, long *cellindex)
{
  for (long i = 0; i < ncells; ++i) cellindex[i] = std::lround(data[i]);
}

static void
free_cellindex(CellIndex *cellindex)
{
  if (cellindex->neighbor) Free(cellindex->neighbor);
  if (cellindex->parent) Free(cellindex->parent);
  if (cellindex->child) Free(cellindex->child);
  Free(cellindex);
}

static CellIndex *
read_cellindex(const std::string &filename)
{
  const auto streamID = stream_open_read_locked(filename.c_str());
  const auto vlistID = streamInqVlist(streamID);
  const auto ngrids = vlistNgrids(vlistID);
  int gridID = -1;
  for (int index = 0; index < ngrids; ++index)
    {
      gridID = vlistGrid(vlistID, index);
      if (gridInqType(gridID) == GRID_UNSTRUCTURED && gridInqNvertex(gridID) == 3) break;
    }

  if (gridID == -1) cdo_abort("No ICON grid found in %s!", filename);

  // int nid = CDI_UNDEFID;
  int pid = CDI_UNDEFID;
  // int cid = CDI_UNDEFID;
  const auto nvars = vlistNvars(vlistID);
  for (int varID = 0; varID < nvars; ++varID)
    {
      auto varname = cdo::inq_var_name(vlistID, varID);
      /*
      if (varname == "neighbor_cell_index")
        {
          nid = varID;
          break;
        }
      */
      if (varname == "parent_cell_index")
        {
          pid = varID;
          break;
        }
      /*
      if (varname == "child_cell_index")
        {
          cid = varID;
          break;
        }
      */
    }

  // if (nid == CDI_UNDEFID) cdo_abort("neighbor_cell_index not found in %s!", filename);
  // if (pid == CDI_UNDEFID) cdo_abort("parent_cell_index not found in %s!", filename);
  // if (cid == CDI_UNDEFID) cdo_abort("child_cell_index not found in %s!", filename);

  const long ncells = gridInqSize(gridID);

  CellIndex *cellindex = (CellIndex *) Malloc(sizeof(CellIndex));
  cellindex->ncells = ncells;

  cellindex->neighbor = nullptr;
  // cellindex->neighbor = (long*) Malloc(3*ncells*sizeof(long));
  cellindex->parent = (long *) Malloc(ncells * sizeof(long));
  cellindex->child = nullptr;
  // cellindex->child    = (cid != CDI_UNDEFID) ? (int*) Malloc(MAX_CHILDS*ncells*sizeof(int)) : nullptr;
  Varray<double> data(ncells);

  for (long i = 0; i < ncells; ++i) cellindex->parent[i] = 0;

  const auto nrecs = streamInqTimestep(streamID, 0);
  for (int recID = 0; recID < nrecs; ++recID)
    {
      int varID, levelID;
      size_t nmiss;
      streamInqRecord(streamID, &varID, &levelID);
      if (varID == pid /* || varID == nid || varID == cid */)
        {
          streamReadRecord(streamID, data.data(), &nmiss);
          // if (varID == pid)
          {
            if (Options::cdoVerbose) cdo_print("Read parent_cell_index");
            copy_data_to_index(ncells, data, cellindex->parent);
          }
          // else if ( varID == nid ) copy_data_to_index(ncells, data, cellindex->neighbor+levelID*ncells);
          // else if ( varID == cid ) copy_data_to_index(ncells, data, cellindex->child+levelID*ncells);
        }
    }

  // Fortran to C index
  for (long i = 0; i < ncells; ++i) cellindex->parent[i] -= 1;
  // for ( long i = 0; i < 3*ncells; ++i ) cellindex->neighbor[i] -= 1;

  streamClose(streamID);

  return cellindex;
}

static int
read_grid(const char *filename)
{
  const auto streamID = stream_open_read_locked(filename);
  const auto vlistID = streamInqVlist(streamID);
  const auto ngrids = vlistNgrids(vlistID);
  int gridID = -1;
  for (int index = 0; index < ngrids; ++index)
    {
      gridID = vlistGrid(vlistID, index);
      if (gridInqType(gridID) == GRID_UNSTRUCTURED && gridInqNvertex(gridID) == 3) break;
    }

  if (gridID == -1) cdo_abort("No ICON grid found in %s!", filename);

  const auto gridID2 = gridDuplicate(gridID);

  streamClose(streamID);

  return gridID2;
}

/**
* Return the first index of element x fits.
*
* If no interval can be found return -1.

* @param *array ascending sorted list
* @param n      length of the sorted list
* @param search the element to find a position for
*/
static long
find_index(int search, long n, const long *array)
{
  long first = 0;
  long last = n - 1;
  long middle = (first + last) / 2;

  while (first <= last)
    {
      if (array[middle] < search)
        first = middle + 1;
      else if (array[middle] == search)
        {
          for (long i = middle; i >= 0; i--)
            {
              if (array[i] == search)
                middle = i;
              else
                break;
            }
          return middle;
        }
      else
        last = middle - 1;

      middle = (first + last) / 2;
    }

  return -1;
}

struct SortInfo
{
  int p, i;
};

static bool
cmpsinfo(const SortInfo &a, const SortInfo &b)
{
  return a.p < b.p;
}

static void
compute_child_from_parent(CellIndex *cellindex1, CellIndex *cellindex2)
{
  if (Options::cdoVerbose) cdo_print("%s", __func__);

  const auto ncells1 = cellindex1->ncells;
  long *parent1 = cellindex1->parent;

  std::vector<long> idx1(ncells1);
  for (long i = 0; i < ncells1; ++i) idx1[i] = i;
  for (long i = 1; i < ncells1; ++i)
    if (parent1[i] < parent1[i - 1])
      {
        if (Options::cdoVerbose) cdo_print("Sort parent index of %s!", cellindex1->filename);
        std::vector<SortInfo> sinfo(ncells1);
        for (long j = 0; j < ncells1; ++j)
          {
            sinfo[j].p = parent1[j];
            sinfo[j].i = idx1[j];
          }
        std::sort(sinfo.begin(), sinfo.end(), cmpsinfo);
        for (long j = 0; j < ncells1; ++j)
          {
            parent1[j] = sinfo[j].p;
            idx1[j] = sinfo[j].i;
          }
        break;
      }

  const auto ncells2 = cellindex2->ncells;
  long *child2 = (long *) Malloc(MAX_CHILDS * ncells2 * sizeof(long));
  cellindex2->child = child2;
  for (long i = 0; i < ncells2; ++i)
    {
      for (long k = 0; k < MAX_CHILDS; ++k) child2[i * MAX_CHILDS + k] = -1;
      long j = find_index(i, ncells1, parent1);
      if (j < 0) continue;
      for (long k = 0; k < MAX_CHILDS; ++k)
        {
          if (i != parent1[j + k]) break;
          //  child2[i*MAX_CHILDS+k] = j+k;
          child2[i * MAX_CHILDS + k] = idx1[j + k];
        }
      // if ( i%10000 == 0 ) printf("%d %d %d %d %d %d\n", i, j, parent1[j], parent1[j+1], parent1[j+2], parent1[j+3]);
    }
}

static void
read_coordinates(const char *filename, long n, double *lon, double *lat, int nv, double *lon_bnds, double *lat_bnds)
{
  const auto streamID = streamOpenRead(filename);
  const auto vlistID = streamInqVlist(streamID);
  const auto ngrids = vlistNgrids(vlistID);
  int gridID = -1;
  for (int index = 0; index < ngrids; ++index)
    {
      gridID = vlistGrid(vlistID, index);
      if (gridInqType(gridID) == GRID_UNSTRUCTURED && (long) gridInqSize(gridID) == n && gridInqNvertex(gridID) == 3) break;
    }

  if (gridID == -1) cdo_abort("No ICON grid with %ld cells found in %s!", n, filename);

  gridInqXvals(gridID, lon);
  gridInqYvals(gridID, lat);

  // Convert lat/lon units if required
  cdo_grid_to_radian(gridID, CDI_XAXIS, n, lon, "grid center lon");
  cdo_grid_to_radian(gridID, CDI_YAXIS, n, lat, "grid center lat");

  if (nv == 3 && lon_bnds && lat_bnds)
    {
      gridInqXbounds(gridID, lon_bnds);
      gridInqYbounds(gridID, lat_bnds);

      cdo_grid_to_radian(gridID, CDI_XAXIS, n * 3, lon_bnds, "grid corner lon");
      cdo_grid_to_radian(gridID, CDI_YAXIS, n * 3, lat_bnds, "grid corner lat");
    }

  streamClose(streamID);
}

constexpr int MAX_SEARCH = 128;  // the triangles are distorted!

static void
compute_child_from_bounds(CellIndex *cellindex2, long ncells2, double *grid_center_lon2, double *grid_center_lat2,
                          double *grid_corner_lon2, double *grid_corner_lat2, long ncells1, const Varray<double> &grid_center_lon1,
                          const Varray<double> &grid_center_lat1)
{
  if (Options::cdoVerbose) cdo_print("%s", __func__);

  bool xIsCyclic = false;
  size_t dims[2];
  dims[0] = ncells1;
  dims[1] = 0;
  GridPointSearch gps;
  grid_point_search_create(gps, xIsCyclic, dims, ncells1, grid_center_lon1, grid_center_lat1);
  knnWeightsType knnWeights(MAX_SEARCH);
  size_t *nbr_addr = &knnWeights.m_addr[0];

  int ncorner = 3;
  Point3D centerPoint3D;
  Varray<Point3D> cellCorners3D(4);
  Varray<Point> cellCornersPlaneProjection(4);

  long *child2 = (long *) Malloc(MAX_CHILDS * ncells2 * sizeof(long));
  cellindex2->child = child2;
  for (long cellNo2 = 0; cellNo2 < ncells2; ++cellNo2)
    {
      for (int k = 0; k < MAX_CHILDS; ++k) child2[cellNo2 * MAX_CHILDS + k] = -1;

      set_cell_corners_3D(ncorner, &grid_corner_lon2[cellNo2 * ncorner], &grid_corner_lat2[cellNo2 * ncorner], cellCorners3D);
      cellCorners3D[ncorner] = cellCorners3D[0];

      const auto coordinateToIgnore = find_coordinate_to_ignore(cellCorners3D);

      const auto cval
          = (coordinateToIgnore == 1) ? cellCorners3D[0].X : ((coordinateToIgnore == 2) ? cellCorners3D[0].Y : cellCorners3D[0].Z);
      const auto invertResult = (cval < 0.0);

      set_cell_corners_plane_projection(coordinateToIgnore, ncorner, cellCorners3D, cellCornersPlaneProjection);

      auto isClockwise = are_polygon_vertices_arranged_in_clockwise_order(cellCornersPlaneProjection, ncorner + 1);

      if (invertResult) isClockwise = !isClockwise;
      if (isClockwise) continue;

      grid_search_point(gps, grid_center_lon2[cellNo2], grid_center_lat2[cellNo2], knnWeights);

      int k = 0;
      double centerCoordinates[3];
      for (int i = 0; i < MAX_SEARCH; ++i)
        {
          const auto cellNo1 = nbr_addr[i];
          if (cellNo1 < SIZE_MAX)
            {
              gcLLtoXYZ(grid_center_lon1[cellNo1], grid_center_lat1[cellNo1], centerCoordinates);
              centerPoint3D.X = centerCoordinates[0];
              centerPoint3D.Y = centerCoordinates[1];
              centerPoint3D.Z = centerCoordinates[2];

              const auto centerPoint2D = set_center_point_plane_projection(coordinateToIgnore, centerPoint3D);

              auto winding_number = winding_numbers_algorithm(cellCornersPlaneProjection, ncorner + 1, centerPoint2D);

              if (winding_number != 0)
                {
                  if (k >= MAX_CHILDS) cdo_abort("Internal problem, limit of MAX_CHILDS reached (limit=9).");
                  child2[cellNo2 * MAX_CHILDS + k++] = (long) cellNo1;
                }
            }
        }
    }

  grid_point_search_delete(gps);
}

static void
compute_child_from_coordinates(CellIndex *cellindex1, CellIndex *cellindex2)
{
  if (Options::cdoVerbose) cdo_print("%s", __func__);

  const auto ncells1 = cellindex1->ncells;
  const auto ncells2 = cellindex2->ncells;

  Varray<double> lon1(ncells1), lat1(ncells1), lon2(ncells2), lat2(ncells2);
  Varray<double> lon2_bnds(3 * ncells2), lat2_bnds(3 * ncells2);

  read_coordinates(cellindex1->filename, ncells1, lon1.data(), lat1.data(), 0, nullptr, nullptr);
  read_coordinates(cellindex2->filename, ncells2, lon2.data(), lat2.data(), 3, lon2_bnds.data(), lat2_bnds.data());

  compute_child_from_bounds(cellindex2, ncells2, lon2.data(), lat2.data(), lon2_bnds.data(), lat2_bnds.data(), ncells1, lon1, lat1);
}

static void
compute_child(CellIndex *cellindex1, CellIndex *cellindex2)
{
  bool lparent = true;
  const auto ncells1 = cellindex1->ncells;
  long *parent1 = cellindex1->parent;
  {
    long i;
    for (i = 0; i < ncells1; ++i)
      if (parent1[i] >= 0) break;
    if (i == ncells1) lparent = false;
  }

  if (lparent)
    compute_child_from_parent(cellindex1, cellindex2);
  else
    {
      compute_child_from_coordinates(cellindex1, cellindex2);
      // cdo_abort("Missing parent index of %s!", cellindex1->filename);
    }
}

static void
compute_sum(long i, long &n, double &sum, double &sumq, long kci, CellIndex **cellindex, const Varray<double> &array)
{
  // printf("compute: i, kci %d %d\n", i, kci);
  const auto ncells2 = cellindex[kci]->ncells;
  if (i < 0 || i > ncells2) cdo_abort("Child grid cell index %ld out of bounds %ld!", i, ncells2);

  for (int k = 0; k < MAX_CHILDS; ++k)
    {
      long index = cellindex[kci]->child[i * MAX_CHILDS + k];
      if (index == -1) break;
      if (kci == 1)
        {
          sum += array[index];
          sumq += array[index] * array[index];
          n += 1;
        }
      else
        compute_sum(index, n, sum, sumq, kci - 1, cellindex, array);
    }
}

static void
samplegrid(double missval, long nci, CellIndex **cellindex, const Varray<double> &array1, Varray<double> &array2,
           Varray<double> &array3)
{
  static bool lstat = true;
  long kci = nci - 1;
  const auto ncells2 = cellindex[kci]->ncells;
  long nx = 0;
  double x = 0.0;
#ifdef _OPENMP
//#pragma omp parallel for default(shared)
#endif
  for (long i = 0; i < ncells2; ++i)
    {
      long n = 0;
      double sum = 0, sumq = 0;
      compute_sum(i, n, sum, sumq, kci, cellindex, array1);
      array2[i] = n ? sum / n : missval;  // mean
      double var1 = (n * n > n) ? (sumq * n - sum * sum) / (n * n - n) : missval;
      if (var1 < 0 && var1 > -1.e-5) var1 = 0;
      array3[i] = var_to_std(var1, missval);  // std1
      if (lstat && n)
        {
          nx++;
          x += n;
        }
    }
  if (Options::cdoVerbose && lstat)
    {
      lstat = false;
      cdo_print("Mean number of childs %g", nx ? x / nx : 0);
    }
}

void *
Samplegridicon(void *process)
{
  cdo_initialize(process);

  cdo_operator_add("samplegridicon", 0, 0, "sample grids");

  const auto nsamplegrids = cdo_operator_argc();
  if (nsamplegrids < 2) cdo_abort("Parameter missing!");

  std::vector<CellIndex *> cellindex(nsamplegrids);

  for (int i = 0; i < nsamplegrids; ++i)
    {
      cellindex[i] = read_cellindex(cdo_operator_argv(i));
      cellindex[i]->filename = cdo_operator_argv(i).c_str();
      if (Options::cdoVerbose) cdo_print("Found %ld grid cells in %s", cellindex[i]->ncells, cellindex[i]->filename);
    }

  for (int i = 0; i < nsamplegrids - 1; ++i) compute_child(cellindex[i], cellindex[i + 1]);

  const auto gridID2 = read_grid(cdo_operator_argv(nsamplegrids - 1).c_str());

  const auto streamID1 = cdo_open_read(0);
  const auto vlistID1 = cdo_stream_inq_vlist(streamID1);

  long gridsize = vlistGridsizeMax(vlistID1);
  if (Options::cdoVerbose) cdo_print("Source gridsize = %zu", gridsize);
  if (gridsize != cellindex[0]->ncells)
    cdo_abort("Gridsize (%ld) of input stream and first grid (%ld) differ!", gridsize, cellindex[0]->ncells);
  if (vlistNumber(vlistID1) != CDI_REAL) gridsize *= 2;
  Varray<double> array1(gridsize);

  const auto vlistID2 = vlistDuplicate(vlistID1);
  const auto vlistID3 = vlistDuplicate(vlistID1);

  const auto taxisID1 = vlistInqTaxis(vlistID1);
  const auto taxisID2 = taxisDuplicate(taxisID1);
  const auto taxisID3 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);
  vlistDefTaxis(vlistID3, taxisID3);

  const auto ngrids = vlistNgrids(vlistID1);
  for (int index = 0; index < ngrids; ++index)
    {
      const auto gridID = vlistGrid(vlistID1, index);
      const auto gridtype = gridInqType(gridID);
      if (!(gridtype == GRID_UNSTRUCTURED && gridInqNvertex(gridID) == 3))
        cdo_abort("Unsupported gridtype: %s with %d corners", gridNamePtr(gridtype), gridInqNvertex(gridID));

      vlistChangeGridIndex(vlistID2, index, gridID2);
      vlistChangeGridIndex(vlistID3, index, gridID2);
    }

  const auto streamID2 = cdo_open_write(1);
  cdo_def_vlist(streamID2, vlistID2);

  const auto streamID3 = cdo_open_write(2);
  cdo_def_vlist(streamID3, vlistID3);

  auto gridsize2 = gridInqSize(gridID2);
  if (Options::cdoVerbose) cdo_print("Target gridsize = %ld", gridsize2);
  if (vlistNumber(vlistID2) != CDI_REAL) gridsize2 *= 2;
  Varray<double> array2(gridsize2), array3(gridsize2);

  int tsID = 0;
  while (true)
    {
      const auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
      if (nrecs == 0) break;

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_taxis_copy_timestep(taxisID3, taxisID1);
      cdo_def_timestep(streamID2, tsID);
      cdo_def_timestep(streamID3, tsID);

      for (int recID = 0; recID < nrecs; ++recID)
        {
          int varID, levelID;
          cdo_inq_record(streamID1, &varID, &levelID);
          size_t nmiss;
          cdo_read_record(streamID1, array1.data(), &nmiss);

          const auto missval = vlistInqVarMissval(vlistID1, varID);

          samplegrid(missval, nsamplegrids, &cellindex[0], array1, array2, array3);

          nmiss = varray_num_mv(gridsize2, array2, missval);
          cdo_def_record(streamID2, varID, levelID);
          cdo_write_record(streamID2, array2.data(), nmiss);

          nmiss = varray_num_mv(gridsize2, array3, missval);
          cdo_def_record(streamID3, varID, levelID);
          cdo_write_record(streamID3, array3.data(), nmiss);
        }

      tsID++;
    }

  cdo_stream_close(streamID3);
  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  vlistDestroy(vlistID2);
  gridDestroy(gridID2);

  for (int i = 0; i < nsamplegrids; ++i) free_cellindex(cellindex[i]);

  cdo_finish();

  return nullptr;
}
