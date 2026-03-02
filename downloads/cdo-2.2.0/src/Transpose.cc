/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Transpose  transxy         Transpose X/Y
*/

#include <cdi.h>

#include "process_int.h"
#include "matrix_view.h"

void
transxy(int gridID, const Varray<double> &v1, Varray<double> &v2)
{
  auto nx = gridInqXsize(gridID);
  auto ny = gridInqYsize(gridID);
  auto gridsize = gridInqSize(gridID);

  if (gridsize == (nx * ny))
    {
      MatrixView<const double> mV1(v1.data(), ny, nx);
      MatrixView<double> mV2(v2.data(), nx, ny);

      for (size_t j = 0; j < ny; ++j)
        for (size_t i = 0; i < nx; ++i) mV2[i][j] = mV1[j][i];
    }
  else
    {
      for (size_t i = 0; i < gridsize; ++i) v2[i] = v1[i];
    }
}

class ModuleTranspose
{
private:
  CdoStreamID streamID1;
  CdoStreamID streamID2;

  int taxisID1;
  int taxisID2;

  int vlistID1;

  Varray<double> array1;
  Varray<double> array2;

public:
  void
  init(void *process)
  {
    cdo_initialize(process);

    operator_check_argc(0);

    streamID1 = cdo_open_read(0);

    vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = vlistDuplicate(vlistID1);

    auto ngrids = vlistNgrids(vlistID1);
    for (int index = 0; index < ngrids; ++index)
      {
        auto gridID1 = vlistGrid(vlistID1, index);
        auto nx = gridInqXsize(gridID1);
        auto ny = gridInqYsize(gridID1);
        auto gridsize = gridInqSize(gridID1);
        if (gridsize == (nx * ny))
          {
            auto gridID2 = gridCreate(GRID_GENERIC, gridsize);
            gridDefXsize(gridID2, ny);
            gridDefYsize(gridID2, nx);
            vlistChangeGridIndex(vlistID2, index, gridID2);
          }
      }

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);

    auto gridsizemax = vlistGridsizeMax(vlistID1);
    array1 = Varray<double>(gridsizemax);
    array2 = Varray<double>(gridsizemax);
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
            size_t nmiss;
            int varID, levelID;
            cdo_inq_record(streamID1, &varID, &levelID);
            cdo_read_record(streamID1, array1.data(), &nmiss);

            auto gridID = vlistInqVarGrid(vlistID1, varID);
            transxy(gridID, array1, array2);

            cdo_def_record(streamID2, varID, levelID);
            cdo_write_record(streamID2, array2.data(), nmiss);
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
Transpose(void *process)
{
  ModuleTranspose transpose;
  transpose.init(process);
  transpose.run();
  transpose.close();
  return nullptr;
}
