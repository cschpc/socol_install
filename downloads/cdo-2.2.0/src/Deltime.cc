/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include "process_int.h"
#include "util_string.h"
#include "cdo_options.h"

void *
Deltime(void *process)
{
  int dday, dmon;
  const char *cmons[] = { "", "jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec" };

  cdo_initialize(process);

  auto dataIsUnchanged = data_is_unchanged();

  // clang-format off
  auto DELDAY   = cdo_operator_add("delday",   0, 0, nullptr);
  auto DEL29FEB = cdo_operator_add("del29feb", 0, 0, nullptr);
  // clang-format on

  (void) (DELDAY);  // unused

  auto operatorID = cdo_operator_id();

  if (operatorID == DEL29FEB)
    {
      dday = 29;
      dmon = 2;
      operator_check_argc(0);
    }
  else
    {
      int nsel = cdo_operator_argc();
      if (nsel < 1) cdo_abort("Too few arguments!");
      if (nsel > 1) cdo_abort("Too many arguments!");
      const char *sarg = cdo_operator_argv(0).c_str();
      dday = atoi(sarg);
      dmon = 0;
      while (isdigit(*sarg)) sarg++;
      if (isalpha(*sarg))
        {
          char smon[32];
          strncpy(smon, sarg, sizeof(smon) - 1);
          smon[sizeof(smon) - 1] = 0;
          cstr_to_lower(smon);
          int im;
          for (im = 0; im < 12; ++im)
            if (memcmp(smon, cmons[im + 1], 3) == 0) break;

          if (im < 12) dmon = im + 1;
        }
    }

  if (Options::cdoVerbose) cdo_print("delete day %d%s", dday, cmons[dmon]);

  auto streamID1 = cdo_open_read(0);

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto vlistID2 = vlistDuplicate(vlistID1);

  auto taxisID1 = vlistInqTaxis(vlistID1);
  auto taxisID2 = taxisDuplicate(taxisID1);
  taxisDefCalendar(taxisID2, CALENDAR_365DAYS);
  vlistDefTaxis(vlistID2, taxisID2);

  auto streamID2 = cdo_open_write(1);
  cdo_def_vlist(streamID2, vlistID2);

  Field field;

  VarList varList1;
  varListInit(varList1, vlistID1);

  int nfound = 0;
  int tsID = 0;
  int tsID2 = 0;
  while (true)
    {
      auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
      if (nrecs == 0) break;

      auto vDateTime = taxisInqVdatetime(taxisID1);

      int year, month, day;
      cdiDate_decode(vDateTime.date, &year, &month, &day);

      auto copyTimestep = true;
      if (day == dday && (month == dmon || dmon == 0))
        {
          nfound++;
          copyTimestep = false;
          if (Options::cdoVerbose) cdo_print("Delete %4.4d-%2.2d-%2.2d at timestep %d", year, month, day, tsID + 1);
        }

      if (copyTimestep)
        {
          cdo_taxis_copy_timestep(taxisID2, taxisID1);
          cdo_def_timestep(streamID2, tsID2++);

          for (int recID = 0; recID < nrecs; ++recID)
            {
              int varID, levelID;
              cdo_inq_record(streamID1, &varID, &levelID);
              cdo_def_record(streamID2, varID, levelID);
              if (dataIsUnchanged) { cdo_copy_record(streamID2, streamID1); }
              else
                {
                  field.init(varList1[varID]);
                  cdo_read_record(streamID1, field);
                  cdo_write_record(streamID2, field);
                }
            }
        }

      tsID++;
    }

  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  if (nfound == 0) cdo_warning("Day %d%s not found!", dday, cmons[dmon]);

  vlistDestroy(vlistID2);

  cdo_finish();

  return nullptr;
}
