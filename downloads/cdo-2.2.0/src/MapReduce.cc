/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Ralf Müller

*/

#include <cdi.h>

#include "process_int.h"
#include "cdi_lockedIO.h"
#include <mpim_grid.h>
#include "griddes.h"

// read only the first data variable from input filename into a given double pointer
static void
read_first_record(const std::string &filename, double *field)
{
  size_t nmiss;
  int varID, levelID;
  auto streamID = stream_open_read_locked(filename.c_str());
  streamInqTimestep(streamID, 0);
  streamInqRecord(streamID, &varID, &levelID);
  streamReadRecord(streamID, field, &nmiss);
  streamClose(streamID);
}

// count the number of locations, for which the mask is true
static int
countMask(const double *maskField, size_t gridSize, double falseVal)
{
  size_t counter = 0;

  for (size_t i = 0; i < gridSize; ++i)
    {
      if (!DBL_IS_EQUAL(maskField[i], falseVal)) counter++;
    }

  return counter;
}

class ModuleMapReduce
{
  CdoStreamID streamID1;
  CdoStreamID streamID2;

  int taxisID1;
  int taxisID2;

  int vlistID2;

  std::vector<bool> vars;
  size_t maskSize;

  Varray<double> arrayIn;
  Varray<double> arrayOut;

  std::vector<size_t> maskIndexList;

public:
  void
  init(void *process)
  {
    cdo_initialize(process);

    if (cdo_operator_argc() < 1) cdo_abort("Too few arguments!");

    // check input grid type and size - this will be used for selecting relevant variables from the input file
    const auto &maskfilename = cdo_operator_argv(0);
    auto inputGridID = cdo_define_grid(maskfilename);
    auto inputGridSize = gridInqSize(inputGridID);
    auto inputGridType = gridInqType(inputGridID);
    Debug(cdoDebug, "MapReduce: input gridSize: %zu", inputGridSize);

    // create an index list of the relevant locations
    {
      Varray<double> inputMaskField(inputGridSize);
      read_first_record(maskfilename, inputMaskField.data());

      // non-zero values mark the relevant points
      maskSize = countMask(inputMaskField.data(), inputGridSize, 0.0);
      Debug(cdoDebug, "MapReduce: maskSize = %zu", maskSize);

      maskIndexList.resize(maskSize, -1);

      size_t k = 0;
      for (size_t i = 0; i < inputGridSize; ++i)
        {
          if (!DBL_IS_EQUAL(inputMaskField[i], 0.0)) maskIndexList[k++] = i;
        }

      if (k == inputGridSize) cdo_warning("Number of mask values and gridsize is the same, no reduction!");
    }

    // check if coordinated bounds shound not be created
    bool nobounds = false, nocoords = false;
    if (2 <= cdo_operator_argc())
      {
        const auto &coordinatesLimitation = cdo_operator_argv(1);
        if (coordinatesLimitation == "nobounds") nobounds = true;
        if (coordinatesLimitation == "nocoords") nocoords = true;
      }
    // create unstructured output grid including bounds
    auto outputGridID = gridToUnstructuredSelecton(inputGridID, maskIndexList, nocoords, nobounds);

    /* create output vlist: Only variabes which have the same gridtype and
     * gridsize as the input mask should be proessed. Everything else is ignoreds
     * {{{ */
    streamID1 = cdo_open_read(0);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto nvars = vlistNvars(vlistID1);
    vars = std::vector<bool>(nvars, false);

    VarList varList1;
    varListInit(varList1, vlistID1);

    // use vlist flags for marking the corresponding variables
    vlistClearFlag(vlistID1);
    for (int varID = 0; varID < nvars; ++varID)
      {
        auto gridID = varList1[varID].gridID;
        if (inputGridType == gridInqType(gridID) && inputGridSize == gridInqSize(gridID))
          {
            vars[varID] = true;
            auto nlevels = varList1[varID].nlevels;
            for (int levID = 0; levID < nlevels; levID++) vlistDefFlag(vlistID1, varID, levID, true);
          }
        else { cdo_warning("Gridtype or gridsize differ, skipped variable %s!", varList1[varID].name); }
      }

    vlistID2 = vlistCreate();
    cdo_vlist_copy_flag(vlistID2, vlistID1);
    vlistDefNtsteps(vlistID2, vlistNtsteps(vlistID1));
    // }}}

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    // use the new selection grid for all output variables
    auto ngrids = vlistNgrids(vlistID2);
    for (int index = 0; index < ngrids; ++index) vlistChangeGridIndex(vlistID2, index, outputGridID);

    // loop over input fields and mask the data values {{{
    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);

    arrayIn = Varray<double>(inputGridSize);
    arrayOut = Varray<double>(maskSize);
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
                auto varID2 = vlistFindVar(vlistID2, varID);
                auto levelID2 = vlistFindLevel(vlistID2, varID, levelID);

                size_t nmiss;
                cdo_read_record(streamID1, arrayIn.data(), &nmiss);

                for (size_t i = 0; i < maskSize; ++i) arrayOut[i] = arrayIn[maskIndexList[i]];

                cdo_def_record(streamID2, varID2, levelID2);
                cdo_write_record(streamID2, arrayOut.data(), 0);
              }
          }
        tsID++;
      }
    // }}}
  }
  void
  close()
  {
    cdo_stream_close(streamID2);
    cdo_stream_close(streamID1);

    cdo_finish();
  }
};
/*
 * the operators argument has to be a single horizontal field,
 * non-zero values are used to mark the relevant locations
 */
void *
MapReduce(void *process)
{
  ModuleMapReduce mapreduce;
  mapreduce.init(process);
  mapreduce.run();
  mapreduce.close();

  return nullptr;
}
