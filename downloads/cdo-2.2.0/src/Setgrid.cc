/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Setgrid    setgrid         Set grid
      Setgrid    setgridtype     Set grid type
      Setgrid    setgridarea     Set grid area
      Setgrid    setgridmask     Set grid mask
*/

#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "cdi_lockedIO.h"
#include "param_conversion.h"
#include <mpim_grid.h>
#include "griddes.h"
#include "gridreference.h"
#include "util_files.h"

static void
change_grid(int vlistID1, int vlistID2, int gridID2)
{
  auto gridsize2 = gridInqSize(gridID2);
  int found = 0;
  auto ngrids = vlistNgrids(vlistID1);
  for (int index = 0; index < ngrids; ++index)
    {
      if (gridsize2 == gridInqSize(vlistGrid(vlistID1, index)))
        {
          vlistChangeGridIndex(vlistID2, index, gridID2);
          found++;
        }
    }
  if (!found) cdo_warning("No horizontal grid with %zu cells found!", gridsize2);
}

static void
grid_set_type(int vlistID1, int vlistID2, int gridtype, const std::string &gridname, bool lregular, bool lregularnn,
              bool ldereference, bool withBounds, std::vector<int> &grid2_vgpm)
{
  auto needCorners = withBounds ? NeedCorners::Yes : NeedCorners::No;
  int gridID2;
  auto lrgrid = false;
  auto ngrids = vlistNgrids(vlistID1);
  for (int index = 0; index < ngrids; ++index)
    {
      auto gridID1 = vlistGrid(vlistID1, index);
      auto gridtype1 = gridInqType(gridID1);
      gridID2 = -1;

      if (gridtype1 == GRID_GENERIC && gridInqSize(gridID1) == 1) continue;

      if (lregular || lregularnn)
        {
          if (gridtype1 == GRID_GAUSSIAN_REDUCED) gridID2 = gridToRegular(gridID1);
        }
      else if (ldereference)
        {
          auto reference = dereferenceGrid(gridID1);
          if (reference.isValid) gridID2 = reference.gridID;
          if (reference.notFound) cdo_abort("Reference to horizontal grid not found!");
        }
      else
        {
          if (gridtype == GRID_CURVILINEAR)
            {
              gridID2 = (gridtype1 == GRID_CURVILINEAR) ? gridID1 : gridToCurvilinear(gridID1, needCorners);
            }
          else if (gridtype == GRID_UNSTRUCTURED)
            {
              gridID2 = gridToUnstructured(gridID1, NeedCorners::Yes);
              if (gridtype1 == GRID_GME)
                {
                  auto grid2_nvgp = gridInqSize(gridID2);
                  grid2_vgpm.resize(grid2_nvgp);
                  gridInqMaskGME(gridID2, grid2_vgpm.data());
                  gridCompress(gridID2);
                }
            }
          else if (gridtype == GRID_LONLAT && gridtype1 == GRID_CURVILINEAR)
            {
              gridID2 = gridCurvilinearToRegular(gridID1);
              if (gridID2 == -1) cdo_warning("Conversion of curvilinear grid to regular grid failed!");
            }
          else if (gridtype == GRID_LONLAT && gridtype1 == GRID_UNSTRUCTURED)
            {
              gridID2 = -1;
              cdo_warning("Conversion of unstructured grid to regular grid failed!");
            }
          else if (gridtype == GRID_LONLAT && gridtype1 == GRID_GENERIC)
            {
              gridID2 = -1;
              cdo_warning("Conversion of generic grid to regular grid failed!");
            }
          else if (gridtype == GRID_LONLAT && gridtype1 == GRID_LONLAT) { gridID2 = gridID1; }
          else if (gridtype == GRID_PROJECTION)
            {
              if (gridInqType(gridID1) == GRID_PROJECTION)
                gridID2 = gridID1;
              else
                {
                  int projID = gridInqProj(gridID1);
                  if (projID != CDI_UNDEFID && gridInqType(projID) == GRID_PROJECTION) gridID2 = projID;
                }
              if (gridID2 == -1) cdo_warning("Projection not found!");
            }
          else
            cdo_abort("Unsupported grid name: %s", gridname);
        }

      if (gridID2 == -1)
        {
          if (!(lregular || lregularnn)) cdo_abort("Unsupported grid type!");
        }

      if (gridID2 != -1)
        {
          if (lregular || lregularnn) lrgrid = true;
          vlistChangeGridIndex(vlistID2, index, gridID2);
        }
    }

  if ((lregular || lregularnn) && !lrgrid) cdo_warning("No reduced Gaussian grid found!");
}

static void
grid_set_cellarea(int vlistID1, int vlistID2, Varray<double> &gridcellArea)
{
  auto areasize = gridcellArea.size();
  int numGridsChanges = 0;
  auto ngrids = vlistNgrids(vlistID1);
  for (int index = 0; index < ngrids; ++index)
    {
      auto gridID1 = vlistGrid(vlistID1, index);
      auto gridsize = gridInqSize(gridID1);
      if (gridsize == areasize)
        {
          auto gridID2 = gridDuplicate(gridID1);
          gridDefArea(gridID2, gridcellArea.data());
          vlistChangeGridIndex(vlistID2, index, gridID2);
          numGridsChanges++;
        }
    }
  if (!numGridsChanges) cdo_warning("No grid with %zu cells found!", areasize);
}

static void
grid_set_mask(int vlistID1, int vlistID2, Varray<double> &gridmask)
{
  auto masksize = gridmask.size();
  auto ngrids = vlistNgrids(vlistID1);
  for (int index = 0; index < ngrids; ++index)
    {
      auto gridID1 = vlistGrid(vlistID1, index);
      auto gridsize = gridInqSize(gridID1);
      if (gridsize == masksize)
        {
          std::vector<int> mask(masksize);
          for (size_t i = 0; i < masksize; ++i)
            {
              mask[i] = (gridmask[i] < 0 || gridmask[i] > 255) ? 0 : (int) std::lround(gridmask[i]);
            }
          auto gridID2 = gridDuplicate(gridID1);
          gridDefMask(gridID2, mask.data());
          vlistChangeGridIndex(vlistID2, index, gridID2);
        }
    }
}

static void
grid_unset_mask(int vlistID1, int vlistID2)
{
  auto ngrids = vlistNgrids(vlistID1);
  for (int index = 0; index < ngrids; ++index)
    {
      auto gridID1 = vlistGrid(vlistID1, index);
      auto gridID2 = gridDuplicate(gridID1);
      gridDefMask(gridID2, nullptr);
      vlistChangeGridIndex(vlistID2, index, gridID2);
    }
}

static void
grid_set_proj_params(int vlistID1, int vlistID2, const std::string &projparams)
{
  auto ngrids = vlistNgrids(vlistID1);
  for (int index = 0; index < ngrids; ++index)
    {
      auto gridID1 = vlistGrid(vlistID1, index);
      if (gridInqType(gridID1) == GRID_PROJECTION)
        {
          auto gridID2 = gridDuplicate(gridID1);
          cdiDefAttTxt(gridID2, CDI_GLOBAL, "proj4_params", (int) projparams.size(), projparams.c_str());
          vlistChangeGridIndex(vlistID2, index, gridID2);
        }
    }
}

static void
grid_read_cellarea(const std::string &areafile, Varray<double> &gridcellArea)
{
  auto searchName = false;
  std::string filename = areafile;
  std::string varname;

  if (!FileUtils::file_exists(areafile.c_str()))
    {
      auto pos = filename.find_last_of(':');
      if (pos > 1 && pos < (filename.size() - 1))
        {
          varname = filename.substr(pos + 1);
          filename = filename.substr(0, pos);
          searchName = true;
        }
    }

  auto streamID = stream_open_read_locked(filename.c_str());
  auto vlistID = streamInqVlist(streamID);

  VarList varList;
  varListInit(varList, vlistID);

  int svarID = 0;
  if (searchName)
    {
      auto nvars = vlistNvars(vlistID);
      int varID;
      for (varID = 0; varID < nvars; ++varID)
        {
          if (varList[varID].name == varname)
            {
              svarID = varID;
              break;
            }
        }
      if (varID == nvars) cdo_abort("Variable %s not found in %s\n", varname, filename);
    }

  auto nrecs = streamInqTimestep(streamID, 0);
  for (int recID = 0; recID < nrecs; ++recID)
    {
      int varID, levelID;
      streamInqRecord(streamID, &varID, &levelID);
      if (varID == svarID)
        {
          auto areasize = gridInqSize(varList[varID].gridID);
          gridcellArea.resize(areasize);

          size_t nmiss;
          streamReadRecord(streamID, gridcellArea.data(), &nmiss);
          break;
        }
    }

  streamClose(streamID);
}

class ModuleSetgrid
{
  CdoStreamID streamID1;
  CdoStreamID streamID2;

  int taxisID1;
  int taxisID2;

  bool lregular = false;
  bool lregularnn = false;

  std::vector<int> grid2_vgpm;
  Varray<double> gridmask, gridcellArea;
  VarList varList1, varList2;
  Field field;

public:
  void
  init(void *process)
  {
    int number = 0, position = 0;
    std::string griduri;
    std::string projparams;
    std::string gridname;
    std::string gridfile;

    cdo_initialize(process);

    // clang-format off
    auto SETGRID       = cdo_operator_add("setgrid",       0, 0, "grid description file or name");
    auto SETGRIDTYPE   = cdo_operator_add("setgridtype",   0, 0, "grid type");
    auto SETGRIDAREA   = cdo_operator_add("setgridarea",   0, 0, "filename with area weights");
    auto SETGRIDMASK   = cdo_operator_add("setgridmask",   0, 0, "filename with grid mask");
    auto UNSETGRIDMASK = cdo_operator_add("unsetgridmask", 0, 0, nullptr);
    auto SETGRIDNUMBER = cdo_operator_add("setgridnumber", 0, 0, "grid number and optionally grid position");
    auto SETGRIDURI    = cdo_operator_add("setgriduri",    0, 0, "reference URI of the horizontal grid");
    auto USEGRIDNUMBER = cdo_operator_add("usegridnumber", 0, 0, "use existing grid identified by grid number");
    auto SETPROJPARAMS = cdo_operator_add("setprojparams", 0, 0, "proj library parameter (e.g.: +init=EPSG:3413)");
    // clang-format on

    auto operatorID = cdo_operator_id();

    if (operatorID != UNSETGRIDMASK) operator_input_arg(cdo_operator_enter(operatorID));

    if (operatorID == SETGRID)
      {
        operator_check_argc(1);
        gridfile = cdo_operator_argv(0);
      }
    else if (operatorID == SETGRIDTYPE)
      {
        operator_check_argc(1);
        gridname = cdo_operator_argv(0);
      }
    else if (operatorID == SETGRIDAREA)
      {
        operator_check_argc(1);

        grid_read_cellarea(cdo_operator_argv(0), gridcellArea);

        if (Options::cdoVerbose)
          {
            auto areasize = gridcellArea.size();
            auto mmm = varray_min_max_mean(areasize, gridcellArea);
            cdo_print("gridcellAreas: %zu %#12.5g%#12.5g%#12.5g", areasize, mmm.min, mmm.mean, mmm.max);
          }
      }
    else if (operatorID == SETGRIDMASK)
      {
        operator_check_argc(1);
        auto maskfile = cdo_operator_argv(0).c_str();
        auto streamID = stream_open_read_locked(maskfile);
        auto vlistID = streamInqVlist(streamID);

        (void) streamInqTimestep(streamID, 0);
        int varID, levelID;
        streamInqRecord(streamID, &varID, &levelID);

        auto missval = vlistInqVarMissval(vlistID, varID);
        auto gridID = vlistInqVarGrid(vlistID, varID);
        auto masksize = gridInqSize(gridID);
        gridmask.resize(masksize);

        size_t nmiss;
        streamReadRecord(streamID, gridmask.data(), &nmiss);
        streamClose(streamID);

        for (size_t i = 0; i < masksize; ++i)
          if (dbl_is_equal(gridmask[i], missval)) gridmask[i] = 0;
      }
    else if (operatorID == USEGRIDNUMBER)
      {
        operator_check_argc(1);
        number = parameter_to_int(cdo_operator_argv(0));
      }
    else if (operatorID == SETGRIDNUMBER)
      {
        if (cdo_operator_argc() >= 1 && cdo_operator_argc() <= 2)
          {
            number = parameter_to_int(cdo_operator_argv(0));
            if (cdo_operator_argc() == 2) position = parameter_to_int(cdo_operator_argv(1));
          }
        else { operator_check_argc(1); }
      }
    else if (operatorID == SETGRIDURI)
      {
        operator_check_argc(1);
        griduri = cdo_operator_argv(0);
      }
    else if (operatorID == SETPROJPARAMS)
      {
        operator_check_argc(1);
        projparams = cdo_operator_argv(0);
      }

    streamID1 = cdo_open_read(0);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = vlistDuplicate(vlistID1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    if (operatorID == SETGRID)
      {
        auto gridID2 = cdo_define_grid(gridfile);
        change_grid(vlistID1, vlistID2, gridID2);
      }
    else if (operatorID == SETGRIDNUMBER || operatorID == SETGRIDURI || operatorID == USEGRIDNUMBER)
      {
        int gridID2x = -1;
        if (operatorID == SETGRIDNUMBER)
          {
            auto gridID1 = vlistGrid(vlistID1, 0);
            gridID2x = gridCreate(GRID_UNSTRUCTURED, gridInqSize(gridID1));
            cdiDefKeyInt(gridID2x, CDI_GLOBAL, CDI_KEY_NUMBEROFGRIDUSED, number);
            cdiDefKeyInt(gridID2x, CDI_GLOBAL, CDI_KEY_NUMBEROFGRIDINREFERENCE, position);
          }
        else if (operatorID == USEGRIDNUMBER)
          {
            if (number < 1 || number > vlistNgrids(vlistID1))
              cdo_abort("Invalid grid number: %d (max = %d)!", number, vlistNgrids(vlistID1));

            gridID2x = vlistGrid(vlistID1, number - 1);
          }
        else
          {
            auto gridID1 = vlistGrid(vlistID1, 0);
            gridID2x = gridDuplicate(gridID1);
            cdiDefKeyString(gridID2x, CDI_GLOBAL, CDI_KEY_REFERENCEURI, griduri.c_str());
          }

        change_grid(vlistID1, vlistID2, gridID2x);
      }
    else if (operatorID == SETGRIDTYPE)
      {
        int gridtype = -1;
        bool ldereference = false;
        bool withBounds = true;

        // clang-format off
        if      (gridname == "curvilinear0")  {gridtype = GRID_CURVILINEAR; withBounds = false;}
        else if (gridname == "curvilinear")   {gridtype = GRID_CURVILINEAR; withBounds = true;}
        else if (gridname == "cell")           gridtype = GRID_UNSTRUCTURED;
        else if (gridname == "unstructured0") {gridtype = GRID_UNSTRUCTURED; withBounds = false;}
        else if (gridname == "unstructured")  {gridtype = GRID_UNSTRUCTURED; withBounds = true;}
        else if (gridname == "generic")        gridtype = GRID_GENERIC;
        else if (gridname == "dereference")    ldereference = true;
        else if (gridname == "lonlat")         gridtype = GRID_LONLAT;
        else if (gridname == "gaussian")       gridtype = GRID_GAUSSIAN;
        else if (gridname == "regularnn")     {gridtype = GRID_GAUSSIAN; lregularnn = true;}
        else if (gridname == "regular")       {gridtype = GRID_GAUSSIAN; lregular = true;}
        else if (gridname == "projection")    {gridtype = GRID_PROJECTION;}
        else cdo_abort("Unsupported grid name: %s", gridname);
        // clang-format on

        grid_set_type(vlistID1, vlistID2, gridtype, gridname, lregular, lregularnn, ldereference, withBounds, grid2_vgpm);
      }
    else if (operatorID == SETGRIDAREA) { grid_set_cellarea(vlistID1, vlistID2, gridcellArea); }
    else if (operatorID == SETGRIDMASK) { grid_set_mask(vlistID1, vlistID2, gridmask); }
    else if (operatorID == UNSETGRIDMASK) { grid_unset_mask(vlistID1, vlistID2); }
    else if (operatorID == SETPROJPARAMS) { grid_set_proj_params(vlistID1, vlistID2, projparams); }

    streamID2 = cdo_open_write(1);

    cdo_def_vlist(streamID2, vlistID2);
    // vlistPrint(vlistID2);

    varListInit(varList1, vlistID1);
    varListInit(varList2, vlistID2);

    auto nvars = vlistNvars(vlistID1);
    if (lregular || lregularnn)
      for (int varID = 0; varID < nvars; ++varID) varList1[varID].memType = MemType::Double;
    if (lregular || lregularnn)
      for (int varID = 0; varID < nvars; ++varID) varList2[varID].memType = MemType::Double;
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
            cdo_def_record(streamID2, varID, levelID);

            field.init((lregular || lregularnn) ? varList2[varID] : varList1[varID]);
            cdo_read_record(streamID1, field);
            auto nmiss = field.nmiss;

            auto gridID1 = varList1[varID].gridID;
            if (lregular || lregularnn)
              {
                if (gridInqType(gridID1) == GRID_GAUSSIAN_REDUCED)
                  {
                    auto missval = varList1[varID].missval;
                    int lnearst = lregularnn ? 1 : 0;
                    field2regular(gridID1, varList2[varID].gridID, missval, field.vec_d.data(), nmiss, lnearst);
                  }
              }
            else if (gridInqType(gridID1) == GRID_GME)
              {
                auto gridsize = varList1[varID].gridsize;
                size_t j = 0;
                if (field.memType == MemType::Float)
                  {
                    for (size_t i = 0; i < gridsize; ++i)
                      if (grid2_vgpm[i]) field.vec_f[j++] = field.vec_f[i];
                  }
                else
                  {
                    for (size_t i = 0; i < gridsize; ++i)
                      if (grid2_vgpm[i]) field.vec_d[j++] = field.vec_d[i];
                  }
              }

            cdo_write_record(streamID2, field);
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
Setgrid(void *process)
{
  ModuleSetgrid setgrid;
  setgrid.init(process);
  setgrid.run();
  setgrid.close();
  return nullptr;
}
