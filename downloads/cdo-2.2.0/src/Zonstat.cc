/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Zonstat    zonmin          Zonal minimum
      Zonstat    zonmax          Zonal maximum
      Zonstat    zonrange        Zonal range
      Zonstat    zonsum          Zonal sum
      Zonstat    zonmean         Zonal mean
      Zonstat    zonavg          Zonal average
      Zonstat    zonstd          Zonal standard deviation
      Zonstat    zonstd1         Zonal standard deviation [Normalize by (n-1)]
      Zonstat    zonvar          Zonal variance
      Zonstat    zonvar1         Zonal variance [Normalize by (n-1)]
      Zonstat    zonpctl         Zonal percentiles
*/

#include <cdi.h>

#include <utility>

#include "cdo_options.h"
#include "cdo_vlist.h"
#include "process_int.h"
#include "param_conversion.h"
#include <mpim_grid.h>
#include "griddes.h"
#include "field_functions.h"

void remap_weights_zonal_mean(const int gridID1, const int gridID2, Varray2D<size_t> &remapIndices, Varray2D<double> &remapWeights);
void remap_zonal_mean(const Varray2D<size_t> &remapIndices, const Varray2D<double> &remapWeights, const Field &field1,
                      Field &field2);

template <typename T>
static void
reorder_field(Varray<T> &v, const std::vector<int> &hpRingIndices, size_t numIndices)
{
  Varray<T> vtmp(numIndices);
  for (size_t i = 0; i < numIndices; ++i) vtmp[i] = v[hpRingIndices[i]];
  for (size_t i = 0; i < numIndices; ++i) v[i] = vtmp[i];
}

static void
reorder_field(Field &field, const std::vector<int> &hpRingIndices)
{
  auto numIndices = hpRingIndices.size();
  if (numIndices)
    {
      if (field.memType == MemType::Float)
        reorder_field(field.vec_f, hpRingIndices, numIndices);
      else
        reorder_field(field.vec_d, hpRingIndices, numIndices);
    }
}

static int
define_reduced_grid(int gridID1, size_t ysize, std::vector<int> &hpReducedPoints)
{
  auto gridsize = gridInqSize(gridID1);
  auto gridID = gridCreate(GRID_GAUSSIAN_REDUCED, gridsize);
  gridDefYsize(gridID, ysize);
  gridDefReducedPoints(gridID, ysize, hpReducedPoints.data());
  return gridID;
}

class ZonStat
{
  int gridIDdestroy = -1, gridID1 = -1, gridID2 = -1;
  int zongridID = -1;
  int sourceGridIsRegular = true;

  int taxisID1;
  int taxisID2;

  CdoStreamID streamID1;
  CdoStreamID streamID2;

  Field field1, field2;
  std::vector<int> hpRingIndices;

  VarList varList1;

  int operfunc;
  int nlatmax;

  Varray2D<size_t> remapIndices;
  Varray2D<double> remapWeights;

  double pn = 0.0;

  void
  add_operators(void)
  {
    // clang-format off
    cdo_operator_add("zonmin",    FieldFunc_Min,    0, nullptr);
    cdo_operator_add("zonmax",    FieldFunc_Max,    0, nullptr);
    cdo_operator_add("zonrange",  FieldFunc_Range,  0, nullptr);
    cdo_operator_add("zonsum",    FieldFunc_Sum,    0, nullptr);
    cdo_operator_add("zonmean",   FieldFunc_Mean,   0, nullptr);
    cdo_operator_add("zonavg",    FieldFunc_Avg,    0, nullptr);
    cdo_operator_add("zonvar",    FieldFunc_Var,    0, nullptr);
    cdo_operator_add("zonvar1",   FieldFunc_Var1,   0, nullptr);
    cdo_operator_add("zonstd",    FieldFunc_Std,    0, nullptr);
    cdo_operator_add("zonstd1",   FieldFunc_Std1,   0, nullptr);
    cdo_operator_add("zonskew",   FieldFunc_Skew,   0, nullptr);
    cdo_operator_add("zonkurt",   FieldFunc_Kurt,   0, nullptr);
    cdo_operator_add("zonmedian", FieldFunc_Median, 0, nullptr);
    cdo_operator_add("zonpctl",   FieldFunc_Pctl,   0, nullptr);
    // clang-format on
  }

public:
  void
  init(void *process)
  {
    cdo_initialize(process);

    add_operators();

    auto operatorID = cdo_operator_id();
    operfunc = cdo_operator_f1(operatorID);

    if (operfunc == FieldFunc_Pctl)
      {
        operator_input_arg("percentile number");
        pn = parameter_to_double(cdo_operator_argv(0));
      }
    else if (cdo_operator_argc() == 1 && operfunc == FieldFunc_Mean)
      {
        sourceGridIsRegular = false;
        gridID2 = cdo_define_grid(cdo_operator_argv(0));
        auto gridtype = gridInqType(gridID2);
        if (gridtype != GRID_GAUSSIAN && gridtype != GRID_LONLAT) cdo_abort("Target grid type must be Gaussian or LonLat!");
        if (!gridInqYbounds(gridID2, NULL)) cdo_abort("Target grid cell bounds missing!");
        if (gridInqXsize(gridID2) > 1) cdo_abort("Target grid must be zonal!");
      }
    else
      {
        operator_check_argc(0);
      }

    streamID1 = cdo_open_read(0);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = vlistDuplicate(vlistID1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    auto ngrids = vlistNgrids(vlistID1);
    for (int index = 0; index < ngrids; ++index)
      {
        auto gridID = vlistGrid(vlistID1, index);
        auto xsize = gridInqXsize(gridID);
        auto ysize = gridInqYsize(gridID);
        auto gridtype = gridInqType(gridID);
        if (xsize > 1 || gridtype == GRID_GAUSSIAN_REDUCED || is_healpix_grid(gridID))
          {
            if (gridID1 == -1) gridID1 = gridID;
          }
        else
          {
            if (ysize > 1 && zongridID == -1) zongridID = gridID;
          }
      }

    int ndiffgrids = 0;
    for (int index = 0; index < ngrids; ++index)
      {
        auto gridID = vlistGrid(vlistID1, index);
        if (zongridID != -1 && zongridID == gridID) continue;
        if (gridID1 != gridID) ndiffgrids++;
      }
    if (zongridID == -1 && gridID1 == -1) cdo_abort("Unsupported grid type!");
    if (ndiffgrids) cdo_abort("Too many different grids!");

    if (gridID1 != -1)
      {
        auto gridtype = gridInqType(gridID1);
        if (sourceGridIsRegular)
          {
            if (gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN || gridtype == GRID_GAUSSIAN_REDUCED
                || gridtype == GRID_GENERIC || is_healpix_grid(gridID1))
              {
                gridID2 = (zongridID != -1 && gridInqYsize(zongridID) == gridInqYsize(gridID1)) ? zongridID : gridToZonal(gridID1);
              }
            else
              {
                if (operfunc == FieldFunc_Mean)
                  {
                    cdo_print("Add zonal grid description to calculate a zonal mean for data on non-rectangular grids.");
                    cdo_print("A predefined zonal description is zonal_<DY>. DY is the increment of the latitudes in degrees.");
                    cdo_print("Example for 2 degree latitude bins:  cdo zonmean,zonal_2 infile outfile");
                  }
                cdo_abort("Unsupported gridtype: %s", gridNamePtr(gridtype));
              }
          }
        else
          {
            auto gridID = generate_full_cell_grid(gridID1);
            if (gridID != gridID1) gridIDdestroy = gridID1 = gridID;
          }
      }
    else
      {
        gridID2 = zongridID;
        cdo_warning("Input stream already contains zonal data!");
      }

    if (gridID2 == -1) cdo_abort("Internal error, target grid undefined!");

    for (int index = 0; index < ngrids; ++index) vlistChangeGridIndex(vlistID2, index, gridID2);

    if (Options::cdoChunkType == CDI_UNDEFID)
      {
        auto nvars = vlistNvars(vlistID2);
        for (int varID = 0; varID < nvars; ++varID) cdiDefKeyInt(vlistID2, varID, CDI_KEY_CHUNKTYPE, CDI_CHUNK_AUTO);
      }

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);

    varListInit(varList1, vlistID1);

    nlatmax = gridInqYsize(gridID2);

    auto zonVar = varList1[0];
    zonVar.gridID = gridID2;
    zonVar.gridsize = nlatmax;
    zonVar.nlevels = 1;
    zonVar.memType = MemType::Double;
    field2.init(zonVar);

    if (!sourceGridIsRegular) remap_weights_zonal_mean(gridID1, gridID2, remapIndices, remapWeights);
    if (is_healpix_grid(gridID1))
      {
        auto [nside, order] = cdo::get_healpix_params(gridID1);
        std::vector<int> hpReducedPoints;
        hp_generate_ring_indices(order, nside, gridInqSize(gridID1), hpRingIndices, hpReducedPoints);
        gridIDdestroy = gridID1 = define_reduced_grid(gridID1, gridInqYsize(gridID2), hpReducedPoints);
      }
  }

  void
  step(int tsID, int nrecs)
  {
    cdo_taxis_copy_timestep(taxisID2, taxisID1);
    cdo_def_timestep(streamID2, tsID);

    for (int recID = 0; recID < nrecs; ++recID)
      {
        int varID, levelID;
        cdo_inq_record(streamID1, &varID, &levelID);
        field1.init(varList1[varID]);
        cdo_read_record(streamID1, field1);
        field1.grid = gridID1;

        field2.missval = field1.missval;

        if (zongridID != -1 && zongridID == field1.grid) { field_ncopy(nlatmax, field1, field2); }
        else if (sourceGridIsRegular)
          {
            if (is_healpix_grid(varList1[varID].gridID)) reorder_field(field1, hpRingIndices);
            (operfunc == FieldFunc_Pctl) ? zonal_pctl(field1, field2, pn) : zonal_function(field1, field2, operfunc);
          }
        else
          {
            remap_zonal_mean(remapIndices, remapWeights, field1, field2);
          }

        cdo_def_record(streamID2, varID, levelID);
        cdo_write_record(streamID2, field2);
      }
  }

  void
  run()
  {
    int tsID = 0;
    while (true)
      {
        auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
        if (nrecs == 0) return;
        step(tsID++, nrecs);
      }
  }

  void
  close()
  {

    cdo_stream_close(streamID2);
    cdo_stream_close(streamID1);

    if (gridIDdestroy != -1) gridDestroy(gridIDdestroy);

    cdo_finish();
  }
};

void *
Zonstat(void *process)
{
  ZonStat zonstat;
  zonstat.init(process);
  zonstat.run();
  zonstat.close();

  return nullptr;
}
