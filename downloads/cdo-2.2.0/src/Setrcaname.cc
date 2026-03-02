/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author:

*/

#include <cdi.h>

#include "process_int.h"
#include "readline.h"
#include "cdo_zaxis.h"

#define MAX_LINE_LEN 4096

void *
Setrcaname(void *process)
{
  const char *rcsnames;
  char line[MAX_LINE_LEN];
  char sname[CDI_MAX_NAME], sdescription[CDI_MAX_NAME], sunits[CDI_MAX_NAME];
  int scode, sltype, slevel;
  int zaxisID, ltype, code, nlev;
  int level;

  cdo_initialize(process);

  const auto dataIsUnchanged = data_is_unchanged();

  operator_input_arg("file name with RCA names");
  rcsnames = cdo_operator_argv(0).c_str();

  const auto streamID1 = cdo_open_read(0);

  const auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  const auto vlistID2 = vlistDuplicate(vlistID1);

  const auto nvars = vlistNvars(vlistID2);

  auto fp = std::fopen(rcsnames, "r");
  if (fp != nullptr)
    {
      while (cdo::readline(fp, line, MAX_LINE_LEN))
        {
          std::sscanf(line, "%d\t%d\t%d\t%s\t%s\t%s", &scode, &sltype, &slevel, sname, sdescription, sunits);
          /*
          printf("%s\n", line);
          printf("%d:%d:%d:%s:%s:%s\n", scode, sltype, slevel, sname,
          sdescription, sunits);
          */
          for (int varID = 0; varID < nvars; ++varID)
            {
              code = vlistInqVarCode(vlistID2, varID);
              zaxisID = vlistInqVarZaxis(vlistID2, varID);
              nlev = zaxisInqSize(zaxisID);

              ltype = zaxis_to_ltype(zaxisID);

              if (code == scode)
                {
                  if (ltype == 105)
                    {
                      if (nlev != 1)
                        {
                          cdo_warning("Number of levels should be 1 for level type 105!");
                          cdo_warning("Maybe environment variable SPLIT_LTYPE_105 is not set.");
                          continue;
                        }
                      level = (int) cdo_zaxis_inq_level(zaxisID, 0);
                      if (sltype == 105 && slevel == level)
                        {
                          cdiDefKeyString(vlistID2, varID, CDI_KEY_NAME, sname);
                          cdiDefKeyString(vlistID2, varID, CDI_KEY_LONGNAME, sdescription);
                          cdiDefKeyString(vlistID2, varID, CDI_KEY_UNITS, sunits);
                          break;
                        }
                    }
                  else if (sltype != 105)
                    {
                      cdiDefKeyString(vlistID2, varID, CDI_KEY_NAME, sname);
                      cdiDefKeyString(vlistID2, varID, CDI_KEY_LONGNAME, sdescription);
                      cdiDefKeyString(vlistID2, varID, CDI_KEY_UNITS, sunits);
                      break;
                    }
                }
            }
        }

      std::fclose(fp);
    }
  else { perror(rcsnames); }

  const auto taxisID1 = vlistInqTaxis(vlistID1);
  const auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  const auto streamID2 = cdo_open_write(1);

  cdo_def_vlist(streamID2, vlistID2);

  Field field;

  VarList varList1;
  varListInit(varList1, vlistID1);

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
          cdo_def_record(streamID2, varID, levelID);

          if (dataIsUnchanged) { cdo_copy_record(streamID2, streamID1); }
          else
            {
              field.init(varList1[varID]);
              cdo_read_record(streamID1, field);
              cdo_write_record(streamID2, field);
            }
        }

      tsID++;
    }

  cdo_stream_close(streamID1);
  cdo_stream_close(streamID2);

  vlistDestroy(vlistID2);

  cdo_finish();

  return nullptr;
}
