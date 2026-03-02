/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Selrec     selrec          Select records
*/

#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "param_conversion.h"

void *
Selrec(void *process)
{
  cdo_initialize(process);

  if (process_self().m_ID != 0) cdo_abort("This operator can't be combined with other operators!");

  operator_input_arg("records");

  const auto intarr = cdo_argv_to_int(cdo_get_oper_argv());
  const int nsel = intarr.size();

  if (Options::cdoVerbose)
    {
      for (int i = 0; i < nsel; ++i) cdo_print("intarr entry: %d %d", i, intarr[i]);
    }

  const auto streamID1 = cdo_open_read(0);

  const auto filetype = cdo_inq_filetype(streamID1);

  if (filetype == CDI_FILETYPE_NC || filetype == CDI_FILETYPE_NC2 || filetype == CDI_FILETYPE_NC4 || filetype == CDI_FILETYPE_NC4C
      || filetype == CDI_FILETYPE_NC5 || filetype == CDI_FILETYPE_NCZARR)
    cdo_abort("This operator does not work on NetCDF data!");

  const auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  const auto vlistID2 = vlistDuplicate(vlistID1);

  const auto taxisID1 = vlistInqTaxis(vlistID1);
  const auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  const auto streamID2 = cdo_open_write(1);
  cdo_def_vlist(streamID2, vlistID2);

  int recordID = 0;
  int tsID = 0;
  while (true)
    {
      const auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
      if (nrecs == 0) break;

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_def_timestep(streamID2, tsID);

      for (int recID = 0; recID < nrecs; ++recID)
        {
          recordID++;
          int varID, levelID;
          cdo_inq_record(streamID1, &varID, &levelID);

          for (int i = 0; i < nsel; ++i)
            {
              if (recordID == intarr[i])
                {
                  cdo_def_record(streamID2, varID, levelID);
                  cdo_copy_record(streamID2, streamID1);

                  break;
                }
            }
        }

      tsID++;
    }

  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  cdo_finish();

  return nullptr;
}
