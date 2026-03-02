/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Selgridcell     selgridcell    Select grid cells by indices
*/

#include <algorithm>
#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include <mpim_grid.h>
#include "gridreference.h"
#include "util_files.h"
#include "param_conversion.h"

int gengridcell(const int gridID1, const size_t gridsize2, const std::vector<long> &cellidx);

static int
genindexgrid(int gridID1, const size_t gridsize2, const std::vector<long> &cellidx)
{
  const auto gridID0 = gridID1;
  auto gridtype1 = gridInqType(gridID1);

  if (gridtype1 == GRID_LONLAT || gridtype1 == GRID_GAUSSIAN || gridtype1 == GRID_PROJECTION)
    {
      gridID1 = gridToCurvilinear(gridID1);
      gridtype1 = GRID_CURVILINEAR;
    }
  else if (gridtype1 == GRID_UNSTRUCTURED && !gridHasCoordinates(gridID1))
    {
      const auto reference = dereferenceGrid(gridID1);
      if (reference.isValid) gridID1 = reference.gridID;
      if (reference.notFound) cdo_warning("Reference to source grid not found!");
    }

  int gridID2 = -1;
  if (gridtype1 == GRID_UNSTRUCTURED || gridtype1 == GRID_CURVILINEAR)
    gridID2 = gengridcell(gridID1, gridsize2, cellidx);
  else if (gridtype1 == GRID_GENERIC && gridInqYsize(gridID1) == 0)
    gridID2 = gengridcell(gridID1, gridsize2, cellidx);

  if (gridID0 != gridID1) gridDestroy(gridID1);

  return gridID2;
}

template <typename T>
static void
select_index(const Varray<T> &v1, Varray<T> &v2, long n, const std::vector<long> &cellidx)
{
  for (long i = 0; i < n; ++i) v2[i] = v1[cellidx[i]];
}

static void
select_index(const Field &field1, Field &field2, long nind, const std::vector<long> &cellidx)
{
  if (field1.memType != field2.memType) cdo_abort("Interal error, memType of field1 and field2 differ!");

  if (field1.memType == MemType::Float)
    select_index(field1.vec_f, field2.vec_f, nind, cellidx);
  else
    select_index(field1.vec_d, field2.vec_d, nind, cellidx);
}

void *
Selgridcell(void *process)
{
  int gridID1 = -1, gridID2;
  int index, gridtype = -1;
  struct sindex_t
  {
    int gridID1, gridID2;
  };
  std::vector<int> indarr;

  cdo_initialize(process);

  cdo_operator_add("selgridcell", 0, 0, "grid cell indices (1-N)");
  const auto DELGRIDCELL = cdo_operator_add("delgridcell", 0, 0, "grid cell indices (1-N)");

  operator_input_arg(cdo_operator_enter(0));

  const auto operatorID = cdo_operator_id();

  if (cdo_operator_argc() < 1) cdo_abort("Too few arguments!");

  int nind = 0;
  if (cdo_operator_argc() == 1)
    {
      bool maskfile = true, indexfile = false;
      const char *filename = cdo_operator_argv(0).c_str();
      if (strncmp(filename, "index=", 6) == 0)
        {
          filename += 6;
          indexfile = true;
        }
      else if (strncmp(filename, "mask=", 5) == 0)
        {
          filename += 5;
          maskfile = true;
        }

      if (FileUtils::file_exists(filename))
        {
          if (indexfile)
            {
              size_t cdo_read_index(const char *indexfile, std::vector<int> &indarr);
              const auto n = cdo_read_index(filename, indarr);
              nind = n;
              if (nind == 0) cdo_abort("Index file %s generates no input!", cdo_operator_argv(0));
            }
          else if (maskfile)
            {
              size_t cdo_read_mask(const char *maskfile, std::vector<bool> &imask);
              std::vector<bool> mask;
              const auto n = cdo_read_mask(filename, mask);
              nind = 0;
              for (size_t i = 0; i < n; ++i)
                if (mask[i]) nind++;
              if (nind == 0) cdo_abort("Mask is empty!");

              indarr.resize(nind);
              nind = 0;
              for (size_t i = 0; i < n; ++i)
                if (mask[i]) indarr[nind++] = i;
              if (nind == 0) cdo_abort("Mask file %s generates no input!", cdo_operator_argv(0));
            }
        }
    }

  if (nind == 0)
    {
      indarr = cdo_argv_to_int(cdo_get_oper_argv());
      nind = indarr.size();

      if (Options::cdoVerbose)
        for (int i = 0; i < nind; ++i) cdo_print("int %d = %d", i + 1, indarr[i]);

      for (int i = 0; i < nind; ++i) indarr[i] -= 1;
    }

  if (nind == 0) cdo_abort("Argument %s generates no input!", cdo_operator_argv(0));

  auto minmax = std::minmax_element(indarr.begin(), indarr.end());
  const auto indmin = *minmax.first;
  const auto indmax = *minmax.second;
  if (indmin < 0) cdo_abort("Index < 1 not allowed!");

  const auto streamID1 = cdo_open_read(0);

  const auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  const auto vlistID2 = vlistDuplicate(vlistID1);

  const auto taxisID1 = vlistInqTaxis(vlistID1);
  const auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  VarList varList1;
  varListInit(varList1, vlistID1);

  const auto nvars = vlistNvars(vlistID1);
  std::vector<bool> vars(nvars, false);

  const auto ngrids = vlistNgrids(vlistID1);
  std::vector<sindex_t> sindex(ngrids);

  long ncells = nind;
  std::vector<long> cellidx;
  if (operatorID == DELGRIDCELL)
    {
      const auto gridsize = vlistGridsizeMax(vlistID1);
      ncells = gridsize - nind;
      cellidx.resize(gridsize, 1);
      for (long i = 0; i < nind; ++i) cellidx[indarr[i]] = 0;
      long j = 0;
      for (size_t i = 0; i < gridsize; ++i)
        if (cellidx[i] == 1) cellidx[j++] = i;
      if (j != ncells) cdo_abort("Internal error; number of cells differ");
    }
  else
    {
      cellidx.resize(nind);
      for (int i = 0; i < nind; ++i) cellidx[i] = indarr[i];
    }

  if (ncells == 0) cdo_abort("Mask is empty!");

  for (index = 0; index < ngrids; ++index)
    {
      gridID1 = vlistGrid(vlistID1, index);
      gridtype = gridInqType(gridID1);

      const auto gridsize = gridInqSize(gridID1);
      if (gridsize == 1) continue;
      if ((size_t) indmax >= gridsize)
        {
          cdo_warning("Max grid index is greater than grid size, skipped grid %d!", index + 1);
          continue;
        }

      gridID2 = genindexgrid(gridID1, ncells, cellidx);
      if (gridID2 == -1)
        {
          cdo_warning("Unsupported grid type >%s<, skipped grid %d!", gridNamePtr(gridtype), index + 1);
          continue;
        }

      sindex[index].gridID1 = gridID1;
      sindex[index].gridID2 = gridID2;

      vlistChangeGridIndex(vlistID2, index, gridID2);

      for (int varID = 0; varID < nvars; ++varID)
        if (gridID1 == varList1[varID].gridID) vars[varID] = true;
    }

  {
    int varID;
    for (varID = 0; varID < nvars; ++varID)
      if (vars[varID]) break;

    if (varID >= nvars) cdo_abort("No variables selected!");
  }

  VarList varList2;
  varListInit(varList2, vlistID2);

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
          field1.init(varList1[varID]);
          cdo_read_record(streamID1, field1);

          cdo_def_record(streamID2, varID, levelID);

          if (vars[varID])
            {
              field2.init(varList2[varID]);
              select_index(field1, field2, ncells, cellidx);

              if (field1.nmiss) field2.nmiss = field_num_mv(field2);

              cdo_write_record(streamID2, field2);
            }
          else { cdo_write_record(streamID2, field1); }
        }

      tsID++;
    }

  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  vlistDestroy(vlistID2);

  cdo_finish();

  return nullptr;
}
