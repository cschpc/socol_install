/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

*/

#include <cdi.h>

#include "process_int.h"
#include <mpim_grid.h>
#include "cdi_lockedIO.h"
#include "field_functions.h"

static void
varrms(const Varray<double> &w, const FieldVector &field1, const FieldVector &field2, Field &field3)
{
  auto grid1 = field1[0].grid;
  auto grid2 = field2[0].grid;
  auto missval1 = field1[0].missval;
  auto missval2 = field2[0].missval;
  double rsum = 0.0, rsumw = 0.0;

  auto nlev = field1.size();
  auto len = gridInqSize(grid1);
  if (len != gridInqSize(grid2)) cdo_abort("fields have different size!");

  // if ( nmiss1 )
  {
    for (size_t k = 0; k < nlev; ++k)
      {
        auto array1 = field1[k].vec_d;
        auto array2 = field2[k].vec_d;
        for (size_t i = 0; i < len; ++i) /*	  if ( !DBL_IS_EQUAL(w[i], missval1) ) */
          {
            rsum = ADDMN(rsum, MULMN(w[i], MULMN(SUBMN(array2[i], array1[i]), SUBMN(array2[i], array1[i]))));
            rsumw = ADDMN(rsumw, w[i]);
          }
      }
  }
  /*
else
  {
    for ( i = 0; i < len; i++ )
      {
        rsum  += w[i] * array1[i];
        rsumw += w[i];
      }
  }
  */

  const double ravg = SQRTMN(DIVMN(rsum, rsumw));

  size_t rnmiss = 0;
  if (DBL_IS_EQUAL(ravg, missval1)) rnmiss++;

  field3.vec_d[0] = ravg;
  field3.nmiss = rnmiss;
}

class ModuleVarrms
{
private:
  CdoStreamID streamID1;
  CdoStreamID streamID2;
  CdoStreamID streamID3;

  int taxisID1;
  int taxisID3;

  int vlistID1;
  int vlistID3;

  FieldVector2D vars1;
  FieldVector2D vars2;

  int nvars;
  int lastgrid = -1;
  int oldcode = 0;

  bool needWeights = true;

  Field field3;
  Varray<double> weights;

public:
  void
  init(void *process)
  {

    cdo_initialize(process);

    operator_check_argc(0);

    streamID1 = cdo_open_read(0);
    streamID2 = cdo_open_read(1);

    vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = cdo_stream_inq_vlist(streamID2);

    double slon = 0.0, slat = 0.0;
    auto gridID3 = gridCreate(GRID_LONLAT, 1);
    gridDefXsize(gridID3, 1);
    gridDefYsize(gridID3, 1);
    gridDefXvals(gridID3, &slon);
    gridDefYvals(gridID3, &slat);

    vlistClearFlag(vlistID1);
    nvars = vlistNvars(vlistID1);
    for (int varID = 0; varID < nvars; ++varID) vlistDefFlag(vlistID1, varID, 0, true);

    vlistID3 = vlistCreate();
    cdo_vlist_copy_flag(vlistID3, vlistID1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID3 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID3, taxisID3);

    auto ngrids = vlistNgrids(vlistID1);
    int index = 0;
    auto gridID1 = vlistGrid(vlistID1, index);

    if (needWeights && gridInqType(gridID1) != GRID_LONLAT && gridInqType(gridID1) != GRID_GAUSSIAN)
      cdo_abort("Unsupported gridtype: %s", gridNamePtr(gridInqType(gridID1)));

    vlistChangeGridIndex(vlistID3, index, gridID3);
    if (ngrids > 1) cdo_abort("Too many different grids!");

    streamID3 = cdo_open_write(2);
    cdo_def_vlist(streamID3, vlistID3);

    fields_from_vlist(vlistID1, vars1, FIELD_VEC);
    fields_from_vlist(vlistID2, vars2, FIELD_VEC);

    auto gridsizemax = vlistGridsizeMax(vlistID1);
    if (needWeights) weights.resize(gridsizemax);

    field3.resize(1);
    field3.grid = gridID3;
  }

  void
  run()
  {
    int tsID = 0;
    while (true)
      {
        auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
        if (nrecs == 0) break;

        auto nrecs2 = cdo_stream_inq_timestep(streamID2, tsID);
        if (nrecs2 == 0) cdo_abort("Input streams have different number of timesteps!");

        cdo_taxis_copy_timestep(taxisID3, taxisID1);
        cdo_def_timestep(streamID3, tsID);

        for (int recID = 0; recID < nrecs; ++recID)
          {
            size_t nmiss;
            int varID, levelID;
            cdo_inq_record(streamID1, &varID, &levelID);
            cdo_read_record(streamID1, vars1[varID][levelID].vec_d.data(), &nmiss);
            if (nmiss) cdo_abort("Missing values unsupported for this operator!");

            cdo_inq_record(streamID2, &varID, &levelID);
            cdo_read_record(streamID2, vars2[varID][levelID].vec_d.data(), &nmiss);
            if (nmiss) cdo_abort("Missing values unsupported for this operator!");
          }

        for (int varID = 0; varID < nvars; ++varID)
          {
            auto wstatus = false;
            auto gridID = vars1[varID][0].grid;
            if (needWeights && gridID != lastgrid)
              {
                lastgrid = gridID;
                wstatus = gridcell_weights(gridID, weights);
              }
            auto code = vlistInqVarCode(vlistID1, varID);
            if (wstatus != 0 && tsID == 0 && code != oldcode)
              cdo_warning("Using constant area weights for code %d!", oldcode = code);

            field3.missval = vars1[varID][0].missval;
            varrms(weights, vars1[varID], vars2[varID], field3);

            cdo_def_record(streamID3, varID, 0);
            cdo_write_record(streamID3, field3.vec_d.data(), field3.nmiss);
          }

        tsID++;
      }
  }

  void
  close()
  {
    cdo_stream_close(streamID3);
    cdo_stream_close(streamID2);
    cdo_stream_close(streamID1);

    vlistDestroy(vlistID3);

    cdo_finish();
  }
};
void *
Varrms(void *process)
{
  ModuleVarrms varrms;
  varrms.init(process);
  varrms.run();
  varrms.close();

  return nullptr;
}
