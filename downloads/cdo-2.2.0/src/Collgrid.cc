/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <algorithm>  // sort

#include <cdi.h>

#include "cdo_rlimit.h"
#include "process_int.h"
#include "cdo_vlist.h"
#include "param_conversion.h"
#include <mpim_grid.h>
#include "util_files.h"
#include "cdo_options.h"
#include "cdi_lockedIO.h"

static int globalGridType = CDI_UNDEFID;

struct GridInfo2
{
  int globalIndicesID;
  bool needed;
};

struct CollgridInfo
{
  CdoStreamID streamID;
  int vlistID;
  VarList varList;
  Field field;
  std::vector<std::vector<long>> cellIndex;
};

struct xyinfoType
{
  double x = 0.0, y = 0.0;
  int id = -1;
};

static bool
cmpx(const xyinfoType &a, const xyinfoType &b)
{
  return (a.x < b.x);
}

static bool
cmpxy_lt(const xyinfoType &a, const xyinfoType &b)
{
  return (a.y < b.y || (std::fabs(a.y - b.y) <= 0 && a.x < b.x));
}

static bool
cmpxy_gt(const xyinfoType &a, const xyinfoType &b)
{
  return (a.y > b.y || (std::fabs(a.y - b.y) <= 0 && a.x < b.x));
}

static int
gen_coll_grid(int ngrids, int nfiles, std::vector<CollgridInfo> &collgridInfo, int gindex, long nxblocks)
{
  auto isSouthNorth = true;
  auto isRegular = false;
  auto isCurvilinear = false;

  long nx = (nxblocks != -1) ? nxblocks : -1;

  auto gridID = vlistGrid(collgridInfo[0].vlistID, gindex);
  auto gridtype0 = (globalGridType != CDI_UNDEFID) ? globalGridType : gridInqType(gridID);
  if (ngrids > 1 && gridtype0 == GRID_GENERIC && gridInqXsize(gridID) == 0 && gridInqYsize(gridID) == 0) return -1;

  auto isUnstructured = (gridtype0 == GRID_UNSTRUCTURED);
  auto nv = isUnstructured ? gridInqNvertex(gridID) : 0;
  auto withCenter = (globalGridType == CDI_UNDEFID && gridHasCoordinates(gridID));
  auto withBounds = (isUnstructured && globalGridType == CDI_UNDEFID && gridHasBounds(gridID));

  std::vector<xyinfoType> xyinfo(nfiles);
  std::vector<long> xsize(nfiles), ysize(nfiles);
  Varray2D<double> xvals(nfiles), yvals(nfiles);
  Varray2D<double> xbounds(nfiles), ybounds(nfiles);

  for (int fileID = 0; fileID < nfiles; ++fileID)
    {
      gridID = vlistGrid(collgridInfo[fileID].vlistID, gindex);
      auto gridtype = (globalGridType != CDI_UNDEFID) ? globalGridType : gridInqType(gridID);
      if (gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN || gridtype == GRID_PROJECTION)
        isRegular = true;
      else if (gridtype == GRID_CURVILINEAR)
        isCurvilinear = true;
      else if (gridtype == GRID_UNSTRUCTURED)
        isUnstructured = true;
      else if (gridtype == GRID_GENERIC /*&& gridInqXsize(gridID) > 0 && gridInqYsize(gridID) > 0*/)
        isRegular = withCenter;
      else
        cdo_abort("Unsupported grid type: %s!", gridNamePtr(gridtype));

      xsize[fileID] = isUnstructured ? gridInqSize(gridID) : gridInqXsize(gridID);
      ysize[fileID] = isUnstructured ? 1 : gridInqYsize(gridID);
      if (xsize[fileID] == 0) xsize[fileID] = 1;
      if (ysize[fileID] == 0) ysize[fileID] = 1;

      if (isRegular)
        {
          xvals[fileID].resize(xsize[fileID]);
          yvals[fileID].resize(ysize[fileID]);
        }
      else if (isCurvilinear || isUnstructured)
        {
          if (withCenter) xvals[fileID].resize(xsize[fileID] * ysize[fileID]);
          if (withCenter) yvals[fileID].resize(xsize[fileID] * ysize[fileID]);
          if (withBounds) xbounds[fileID].resize(nv * xsize[fileID] * ysize[fileID]);
          if (withBounds) ybounds[fileID].resize(nv * xsize[fileID] * ysize[fileID]);
        }

      if (isRegular || isCurvilinear || isUnstructured)
        {
          if (withCenter) gridInqXvals(gridID, xvals[fileID].data());
          if (withCenter) gridInqYvals(gridID, yvals[fileID].data());
          if (withBounds) gridInqXbounds(gridID, xbounds[fileID].data());
          if (withBounds) gridInqYbounds(gridID, ybounds[fileID].data());
        }
      // printf("fileID %d, gridID %d\n", fileID, gridID);

      xyinfo[fileID].id = fileID;
      if (isRegular)
        {
          xyinfo[fileID].x = xvals[fileID][0];
          xyinfo[fileID].y = yvals[fileID][0];
          if (ysize[fileID] > 1 && yvals[fileID][0] > yvals[fileID][ysize[fileID] - 1]) isSouthNorth = false;
        }
    }

  if (Options::cdoVerbose && isRegular)
    for (int fileID = 0; fileID < nfiles; ++fileID) printf("1 %d %g %g \n", xyinfo[fileID].id, xyinfo[fileID].x, xyinfo[fileID].y);

  if (isRegular)
    {
      std::sort(xyinfo.begin(), xyinfo.end(), cmpx);

      if (Options::cdoVerbose)
        for (int fileID = 0; fileID < nfiles; ++fileID)
          printf("2 %d %g %g \n", xyinfo[fileID].id, xyinfo[fileID].x, xyinfo[fileID].y);

      auto cmpxy = isSouthNorth ? cmpxy_lt : cmpxy_gt;
      std::sort(xyinfo.begin(), xyinfo.end(), cmpxy);

      if (Options::cdoVerbose)
        for (int fileID = 0; fileID < nfiles; ++fileID)
          printf("3 %d %g %g \n", xyinfo[fileID].id, xyinfo[fileID].x, xyinfo[fileID].y);

      if (nx <= 0)
        {
          nx = 1;
          for (int fileID = 1; fileID < nfiles; ++fileID)
            {
              if (DBL_IS_EQUAL(xyinfo[0].y, xyinfo[fileID].y))
                nx++;
              else
                break;
            }
        }
    }
  else
    {
      if (nx <= 0) nx = nfiles;
    }

  const long ny = nfiles / nx;
  if (nx * ny != nfiles) cdo_abort("Number of input files (%ld) and number of blocks (%ldx%ld) differ!", nfiles, nx, ny);

  long xsize2 = 0;
  for (long i = 0; i < nx; ++i) xsize2 += xsize[xyinfo[i].id];
  long ysize2 = 0;
  for (long j = 0; j < ny; ++j) ysize2 += ysize[xyinfo[j * nx].id];
  if (Options::cdoVerbose) cdo_print("xsize2 %ld  ysize2 %ld", xsize2, ysize2);

  {  // verify size of data
    const long xs = xsize[xyinfo[0].id];
    for (long j = 1; j < ny; ++j)
      if (xsize[xyinfo[j * nx].id] != xs) cdo_abort("xsize=%ld differ from first file (xsize=%ld)!", xsize[xyinfo[j * nx].id], xs);
    const long ys = ysize[xyinfo[0].id];
    for (long i = 1; i < nx; ++i)
      if (ysize[xyinfo[i].id] != ys) cdo_abort("ysize=%ld differ from first file (ysize=%ld)!", ysize[xyinfo[i].id], ys);
  }

  Varray<double> xvals2, yvals2;
  Varray<double> xbounds2, ybounds2;
  if (isRegular)
    {
      xvals2.resize(xsize2);
      yvals2.resize(ysize2);
    }
  else if (isCurvilinear || isUnstructured)
    {
      if (withCenter) xvals2.resize(xsize2 * ysize2);
      if (withCenter) yvals2.resize(xsize2 * ysize2);
      if (withBounds) xbounds2.resize(nv * xsize2 * ysize2);
      if (withBounds) ybounds2.resize(nv * xsize2 * ysize2);
    }

  std::vector<long> xoff(nx + 1), yoff(ny + 1);

  xoff[0] = 0;
  for (long i = 0; i < nx; ++i)
    {
      const long idx = xyinfo[i].id;
      if (isRegular) array_copy(xsize[idx], xvals[idx].data(), &xvals2[xoff[i]]);
      xoff[i + 1] = xoff[i] + xsize[idx];
    }

  yoff[0] = 0;
  for (long j = 0; j < ny; ++j)
    {
      const long idx = xyinfo[j * nx].id;
      if (isRegular) array_copy(ysize[idx], yvals[idx].data(), &yvals2[yoff[j]]);
      yoff[j + 1] = yoff[j] + ysize[idx];
    }

  for (int fileID = 0; fileID < nfiles; ++fileID)
    {
      const long idx = xyinfo[fileID].id;
      const long iy = fileID / nx;
      const long ix = fileID - iy * nx;

      const long offset = yoff[iy] * xsize2 + xoff[ix];
      // printf("fileID %d %d, iy %d, ix %d, offset %d\n", fileID, xyinfo[fileID].id, iy, ix, offset);

      long ij = 0;
      for (long j = 0; j < ysize[idx]; ++j)
        for (long i = 0; i < xsize[idx]; ++i)
          {
            if (isCurvilinear || isUnstructured)
              {
                if (withCenter) { xvals2[offset + j * xsize2 + i] = xvals[idx][ij]; }
                if (withCenter) { yvals2[offset + j * xsize2 + i] = yvals[idx][ij]; }
                if (withBounds)
                  {
                    for (long k = 0; k < nv; ++k) xbounds2[(offset + j * xsize2 + i) * nv + k] = xbounds[idx][ij * nv + k];
                  }
                if (withBounds)
                  {
                    for (long k = 0; k < nv; ++k) ybounds2[(offset + j * xsize2 + i) * nv + k] = ybounds[idx][ij * nv + k];
                  }
              }
            collgridInfo[idx].cellIndex[gindex][ij++] = offset + j * xsize2 + i;
          }
    }

  auto gridID2 = gridCreate(gridtype0, xsize2 * ysize2);
  if (!isUnstructured)
    {
      gridDefXsize(gridID2, xsize2);
      gridDefYsize(gridID2, ysize2);
    }
  else if (nv > 0) { gridDefNvertex(gridID2, nv); }

  if (isRegular || isCurvilinear || isUnstructured)
    {
      if (withCenter) gridDefXvals(gridID2, xvals2.data());
      if (withCenter) gridDefYvals(gridID2, yvals2.data());
      if (withBounds) gridDefXbounds(gridID2, xbounds2.data());
      if (withBounds) gridDefYbounds(gridID2, ybounds2.data());
    }

  gridID = vlistGrid(collgridInfo[0].vlistID, gindex);

  grid_copy_names(gridID, gridID2);

  if (gridtype0 == GRID_PROJECTION) grid_copy_mapping(gridID, gridID2);

  return gridID2;
}
/*
static void
coll_cells_reg2d(const Field &field1, Field &field2, const CollgridInfo &collgridInfo, size_t nlon)
{
  auto nx = collgridInfo.nx;
  auto ny = collgridInfo.ny;

  for (size_t j = 0; j < ny; ++j)
    {
      auto offset1 = j * nx;
      auto offset2 = collgridInfo.offset + j * nlon;

      if (field1.memType == MemType::Float)
        for (size_t i = 0; i < nx; ++i) field2.vec_f[offset2 + i] = field1.vec_f[offset1 + i];
      else
        for (size_t i = 0; i < nx; ++i) field2.vec_d[offset2 + i] = field1.vec_d[offset1 + i];
    }
}
*/
static void
collect_cells(const Field &field1, Field &field2, const std::vector<long> &cellIndex)
{
  if (field1.memType != field2.memType) cdo_abort("Interal error, memType of field1 and field2 differ!");

  auto patchSize = field1.size;

  if (field1.memType == MemType::Float)
    for (size_t i = 0; i < patchSize; ++i) field2.vec_f[cellIndex[i]] = field1.vec_f[i];
  else
    for (size_t i = 0; i < patchSize; ++i) field2.vec_d[cellIndex[i]] = field1.vec_d[i];
}

static std::vector<int>
get_var_gridindex(int vlistID)
{
  auto nvars = vlistNvars(vlistID);
  auto ngrids = vlistNgrids(vlistID);

  std::vector<int> varGridIndex(nvars, 0);
  for (int varID = 0; varID < nvars; ++varID)
    {
      auto gridID = vlistInqVarGrid(vlistID, varID);
      for (int index = 0; index < ngrids; ++index)
        {
          if (gridID == vlistGrid(vlistID, index))
            {
              varGridIndex[varID] = index;
              break;
            }
        }
    }

  return varGridIndex;
}

static std::vector<GridInfo2>
get_gridinfo(int vlistID, const VarList &varList, const std::vector<int> &varGridIndex, std::vector<bool> &selectedVars)
{
  auto nvars = vlistNvars(vlistID);
  auto ngrids = vlistNgrids(vlistID);

  std::vector<GridInfo2> gridInfo(ngrids);
  for (int index = 0; index < ngrids; ++index)
    {
      gridInfo[index].globalIndicesID = -1;
      gridInfo[index].needed = false;
    }

  int globalCellIndicesID = -1;
  int globalVertIndicesID = -1;
  int globalEdgeIndicesID = -1;
  for (int varID = 0; varID < nvars; ++varID)
    {
      // clang-format off
      if      (varList[varID].name == "global_cell_indices") globalCellIndicesID = varID;
      else if (varList[varID].name == "global_vert_indices") globalVertIndicesID = varID;
      else if (varList[varID].name == "global_edge_indices") globalEdgeIndicesID = varID;
      // clang-format on
    }
  if (globalCellIndicesID != -1) selectedVars[globalCellIndicesID] = false;
  if (globalVertIndicesID != -1) selectedVars[globalVertIndicesID] = false;
  if (globalEdgeIndicesID != -1) selectedVars[globalEdgeIndicesID] = false;

  if (globalCellIndicesID != -1) gridInfo[varGridIndex[globalCellIndicesID]].globalIndicesID = globalCellIndicesID;
  if (globalVertIndicesID != -1) gridInfo[varGridIndex[globalVertIndicesID]].globalIndicesID = globalVertIndicesID;
  if (globalEdgeIndicesID != -1) gridInfo[varGridIndex[globalEdgeIndicesID]].globalIndicesID = globalEdgeIndicesID;

  for (int varID = 0; varID < nvars; ++varID)
    if (selectedVars[varID]) gridInfo[varGridIndex[varID]].needed = true;

  return gridInfo;
}

static std::vector<bool>
get_selected_vars(int nsel, int noff, int vlistID1, const VarList &varList1)
{
  auto nvars = vlistNvars(vlistID1);
  std::vector<bool> selectedVars(nvars, false);

  if (nsel == 0)
    {
      for (int varID = 0; varID < nvars; ++varID) selectedVars[varID] = true;
    }
  else
    {
      if (Options::cdoVerbose)
        for (int i = 0; i < nsel; ++i) cdo_print("name %d = %s", i + 1, cdo_operator_argv(noff + i));

      std::vector<bool> selfound(nsel);
      for (int i = 0; i < nsel; ++i) selfound[i] = false;

      for (int varID = 0; varID < nvars; ++varID)
        {
          for (int isel = 0; isel < nsel; isel++)
            {
              if (cdo_operator_argv(noff + isel) == varList1[varID].name)
                {
                  selfound[isel] = true;
                  selectedVars[varID] = true;
                }
            }
        }

      int err = 0;
      for (int isel = 0; isel < nsel; isel++)
        {
          if (selfound[isel] == false)
            {
              err++;
              cdo_warning("Variable name %s not found!", cdo_operator_argv(noff + isel));
            }
        }
      if (err) cdo_abort("Could not find all requested variables: (%d/%d)", nsel - err, nsel);
    }

  return selectedVars;
}

static void
select_vars(const std::vector<bool> &selectedVars, int vlistID1, const VarList &varList1)
{
  auto nvars = vlistNvars(vlistID1);
  int numVars = 0;

  for (int varID = 0; varID < nvars; ++varID)
    {
      if (selectedVars[varID])
        {
          numVars++;
          auto nlevels = varList1[varID].nlevels;
          for (int levID = 0; levID < nlevels; levID++) vlistDefFlag(vlistID1, varID, levID, true);
        }
    }

  if (numVars == 0) cdo_abort("No variables selected!");
}

void *
Collgrid(void *process)
{
  int nxblocks = -1;

  cdo_initialize(process);

  auto nfiles = cdo_stream_cnt() - 1;
  auto ofilename = cdo_get_stream_name(nfiles);

  if (!Options::cdoOverwriteMode && FileUtils::file_exists(ofilename) && !FileUtils::user_file_overwrite(ofilename))
    cdo_abort("Outputfile %s already exists!", ofilename);

  std::vector<CollgridInfo> collgridInfo(nfiles);

  cdo::set_numfiles(nfiles + 8);

  for (int fileID = 0; fileID < nfiles; ++fileID)
    {
      auto streamID = cdo_open_read(fileID);
      auto vlistID = cdo_stream_inq_vlist(streamID);
      collgridInfo[fileID].streamID = streamID;
      collgridInfo[fileID].vlistID = vlistID;
      varListInit(collgridInfo[fileID].varList, vlistID);
    }

  const auto &varList1 = collgridInfo[0].varList;
  auto vlistID1 = collgridInfo[0].vlistID;
  vlistClearFlag(vlistID1);

  // check that the contents is always the same
  for (int fileID = 1; fileID < nfiles; ++fileID) vlist_compare(vlistID1, collgridInfo[fileID].vlistID, CMP_NAME | CMP_NLEVEL);

  auto nsel = cdo_operator_argc();
  int noff = 0;

  if (nsel > 0)
    {
      auto argument = cdo_operator_argv(0).c_str();
      if (strcmp(argument, "gridtype=unstructured") == 0)
        {
          nsel--;
          noff++;
          globalGridType = GRID_UNSTRUCTURED;
        }
      else
        {
          auto len = (int) strlen(argument);
          while (--len >= 0 && isdigit(argument[len]))
            ;

          if (len == -1)
            {
              nsel--;
              noff++;
              nxblocks = parameter_to_int(argument);
            }
        }
    }

  auto selectedVars = get_selected_vars(nsel, noff, vlistID1, varList1);
  auto varGridIndex = get_var_gridindex(vlistID1);
  auto gridInfo = get_gridinfo(vlistID1, varList1, varGridIndex, selectedVars);
  select_vars(selectedVars, vlistID1, varList1);

  auto vlistID2 = vlistCreate();
  cdo_vlist_copy_flag(vlistID2, vlistID1);
  vlistDefNtsteps(vlistID2, vlistNtsteps(vlistID1));

  auto taxisID1 = vlistInqTaxis(vlistID1);
  auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  // if ( Options::cdoVerbose ) vlistPrint(vlistID1);
  // if ( Options::cdoVerbose ) vlistPrint(vlistID2);

  auto ngrids1 = vlistNgrids(vlistID1);
  auto ngrids2 = vlistNgrids(vlistID2);

  for (int fileID = 0; fileID < nfiles; ++fileID)
    {
      collgridInfo[fileID].cellIndex.resize(ngrids1);
      for (int gindex = 0; gindex < ngrids1; ++gindex)
        {
          if (gridInfo[gindex].needed)
            {
              auto patchSize = gridInqSize(vlistGrid(collgridInfo[fileID].vlistID, gindex));
              collgridInfo[fileID].cellIndex[gindex].resize(patchSize);
            }
        }
    }

  std::vector<int> gridID2s(ngrids2);

  for (int i2 = 0; i2 < ngrids2; ++i2)
    {
      int i1;
      for (i1 = 0; i1 < ngrids1; ++i1)
        if (vlistGrid(vlistID1, i1) == vlistGrid(vlistID2, i2)) break;

      gridID2s[i2] = gen_coll_grid(ngrids2, nfiles, collgridInfo, i1, nxblocks);
    }

  for (int i = 0; i < ngrids2; ++i)
    {
      if (gridID2s[i] != -1) vlistChangeGridIndex(vlistID2, i, gridID2s[i]);
    }

  VarList varList2;
  varListInit(varList2, vlistID2);

  auto nvars2 = vlistNvars(vlistID2);
  std::vector<bool> collectVars2(nvars2, false);
  for (int varID = 0; varID < nvars2; ++varID)
    {
      auto gridID = varList2[varID].gridID;
      for (int i = 0; i < ngrids2; ++i)
        {
          if (gridID2s[i] != -1 && gridID == vlistGrid(vlistID2, i))
            {
              collectVars2[varID] = true;
              break;
            }
        }
    }

  Field field2;

  auto streamID2 = cdo_open_write(nfiles);
  cdo_def_vlist(streamID2, vlistID2);

  int nrecs0 = 0;
  int tsID = 0;
  do {
      nrecs0 = cdo_stream_inq_timestep(collgridInfo[0].streamID, tsID);
      for (int fileID = 1; fileID < nfiles; ++fileID)
        {
          auto nrecs = cdo_stream_inq_timestep(collgridInfo[fileID].streamID, tsID);
          if (nrecs != nrecs0)
            cdo_abort("Number of records at time step %d of %s and %s differ!", tsID + 1, cdo_get_stream_name(0),
                      cdo_get_stream_name(fileID));
        }

      cdo_taxis_copy_timestep(taxisID2, taxisID1);

      if (nrecs0 > 0) cdo_def_timestep(streamID2, tsID);

      for (int recID = 0; recID < nrecs0; ++recID)
        {
          int varID = 0, levelID = 0;
          for (int fileID = nfiles - 1; fileID >= 0; fileID--) cdo_inq_record(collgridInfo[fileID].streamID, &varID, &levelID);

          // if (Options::cdoVerbose && tsID == 0) printf(" tsID, recID, varID, levelID %d %d %d %d\n", tsID, recID, varID,
          // levelID);

          auto gindex = varGridIndex[varID];
          if (gridInfo[gindex].needed && gridInfo[gindex].globalIndicesID == varID)
            {
              Varray<double> cellIndex;
              for (int fileID = 0; fileID < nfiles; ++fileID)
                {
                  auto patchSize = collgridInfo[fileID].varList[varID].gridsize;
                  if (cellIndex.size() < patchSize) cellIndex.resize(patchSize);
                  size_t nmiss;
                  cdo_read_record(collgridInfo[fileID].streamID, cellIndex.data(), &nmiss);
                  for (size_t i = 0; i < patchSize; ++i) collgridInfo[fileID].cellIndex[gindex][i] = std::lround(cellIndex[i]) - 1;
                }
            }

          if (vlistInqFlag(vlistID1, varID, levelID) == true)
            {
              auto varID2 = vlistFindVar(vlistID2, varID);
              auto levelID2 = vlistFindLevel(vlistID2, varID, levelID);
              // if (Options::cdoVerbose && tsID == 0) printf("varID %d %d levelID %d %d\n", varID, varID2, levelID, levelID2);

              field2.init(varList2[varID2]);
              field_fill(field2, field2.missval);

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
              for (int fileID = 0; fileID < nfiles; ++fileID)
                {
                  auto &field1 = collgridInfo[fileID].field;
                  field1.init(collgridInfo[fileID].varList[varID]);
                  cdo_read_record(collgridInfo[fileID].streamID, field1);

                  if (collectVars2[varID2]) collect_cells(field1, field2, collgridInfo[fileID].cellIndex[gindex]);
                }

              cdo_def_record(streamID2, varID2, levelID2);

              if (collectVars2[varID2])
                {
                  field_num_mv(field2);
                  cdo_write_record(streamID2, field2);
                }
              else { cdo_write_record(streamID2, collgridInfo[0].field); }
            }
        }

      tsID++;
    }
  while (nrecs0 > 0);

  for (int fileID = 0; fileID < nfiles; ++fileID) cdo_stream_close(collgridInfo[fileID].streamID);

  cdo_stream_close(streamID2);

  cdo_finish();

  return nullptr;
}
