/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Writerandom writerandom
*/

#include <cstdlib>
#include <cdi.h>

#include "process_int.h"

class ModuleWriteRandom
{
  CdoStreamID streamID1;
  CdoStreamID streamID2;

  int taxisID1;
  int taxisID2;

  VarList varList1;

public:
  void
  init(void *process)
  {

    cdo_initialize(process);

    operator_check_argc(0);

    streamID1 = cdo_open_read(0);

    const auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    const auto vlistID2 = vlistDuplicate(vlistID1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    streamID2 = cdo_open_write(1);

    cdo_def_vlist(streamID2, vlistID2);

    varListInit(varList1, vlistID1);
  }
  void
  run()
  {
    int tsID = 0;
    while (true)
      {
        const auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
        if (nrecs == 0) break;

        cdo_taxis_copy_timestep(taxisID2, taxisID1);
        cdo_def_timestep(streamID2, tsID);

        Varray2D<double> recdata(nrecs);
        std::vector<int> recvarID(nrecs), reclevelID(nrecs), recindex(nrecs);
        std::vector<size_t> recnmiss(nrecs);

        for (int recID = 0; recID < nrecs; ++recID)
          {
            int varID, levelID;
            cdo_inq_record(streamID1, &varID, &levelID);
            recvarID[recID] = varID;
            reclevelID[recID] = levelID;
            recdata[recID].resize(varList1[varID].gridsize);
            cdo_read_record(streamID1, recdata[recID].data(), &recnmiss[recID]);
          }

        for (int recID = 0; recID < nrecs; ++recID) recindex[recID] = -1;

        for (int rindex = nrecs - 1; rindex >= 0; rindex--)
          {
            const auto index = (int) (rindex * ((double) std::rand()) / ((double) RAND_MAX));
            //	printf("rindex %d %d\n", rindex, index);
            int ipos = -1;
            for (int recID = 0; recID < nrecs; ++recID)
              {
                if (recindex[recID] == -1) ipos++;
                if (recindex[recID] == -1 && ipos == index)
                  {
                    recindex[recID] = rindex;
                    break;
                  }
              }
          }

        // for ( int recID = 0; recID < nrecs; recID++ ) printf("recID %d %d\n", recID, recindex[recID]);

        for (int recID = 0; recID < nrecs; ++recID)
          if (recindex[recID] == -1) cdo_abort("Internal problem! Random initialize.");

        for (int recID = 0; recID < nrecs; ++recID)
          {
            const auto rindex = recindex[recID];
            const auto varID = recvarID[rindex];
            const auto levelID = reclevelID[rindex];
            cdo_def_record(streamID2, varID, levelID);
            cdo_write_record(streamID2, recdata[rindex].data(), recnmiss[rindex]);
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
Writerandom(void *process)
{
  ModuleWriteRandom writerandom;
  writerandom.init(process);
  writerandom.run();
  writerandom.close();

  return nullptr;
}
