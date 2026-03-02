/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include "process_int.h"
#include "varray.h"

void *
Test(void *process)
{
  cdo_initialize(process);

  /*
  const auto streamID1 = cdo_open_read(0);
  const auto streamID2 = cdo_open_write(1);

  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);
  */
  cdo_finish();

  return nullptr;
}

void *
Test2(void *process)
{
  cdo_initialize(process);

  /*
  const auto streamID1 = cdo_open_read(0);
  const auto streamID2 = cdo_open_read(1);
  const auto streamID3 = cdo_open_write(2);

  cdo_stream_close(streamID3);
  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);
  */
  cdo_finish();

  return nullptr;
}

void *
Testdata(void *process)
{
  cdo_initialize(process);

  operator_check_argc(0);

  const auto streamID1 = cdo_open_read(0);

  const auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  const auto taxisID1 = vlistInqTaxis(vlistID1);

  const auto streamID2 = cdo_open_write(1);

  const auto vlistID2 = vlistDuplicate(vlistID1);
  const auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  cdo_def_vlist(streamID2, vlistID2);

  auto gridsize = vlistGridsizeMax(vlistID1);
  Varray<double> array(gridsize);
  Varray<float> fval(gridsize);
  Varray<int> ival(gridsize);
  Varray<unsigned char> cval(gridsize * 4);
  Varray<unsigned char> cval2(gridsize * 4);

  auto fp = std::fopen("testdata", "w");

  int tsID2 = 0;
  int tsID1 = 0;
  while (true)
    {
      const auto nrecs = cdo_stream_inq_timestep(streamID1, tsID1);
      if (nrecs == 0) break;

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_def_timestep(streamID2, tsID2);

      for (int recID = 0; recID < nrecs; ++recID)
        {
          int varID, levelID;
          cdo_inq_record(streamID1, &varID, &levelID);
          cdo_def_record(streamID2, varID, levelID);

          size_t nmiss;
          cdo_read_record(streamID1, array.data(), &nmiss);

          gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
          for (size_t i = 0; i < gridsize; ++i)
            {
              fval[i] = (float) array[i];

              memcpy(&ival[i], &fval[i], 4);
              memcpy(&cval[i * 4], &fval[i], 4);

              cval2[i + gridsize * 0] = cval[i * 4 + 0];
              cval2[i + gridsize * 1] = cval[i * 4 + 1];
              cval2[i + gridsize * 2] = cval[i * 4 + 2];
              cval2[i + gridsize * 3] = cval[i * 4 + 3];

              if (tsID1 == 0 && recID == 0)
                printf("%4zu %3u %3u %3u %3u %d %g\n", i, (unsigned int) cval[4 * i + 0], (unsigned int) cval[4 * i + 1],
                       (unsigned int) cval[4 * i + 2], (unsigned int) cval[4 * i + 3], ival[i], fval[i]);
            }

          cdo_write_record(streamID2, array.data(), nmiss);

          fwrite(cval.data(), 4, gridsize, fp);
        }

      tsID1++;
      tsID2++;
    }

  std::fclose(fp);
  cdo_stream_close(streamID1);
  cdo_stream_close(streamID2);

  cdo_finish();

  return nullptr;
}
