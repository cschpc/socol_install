/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include "cdo_options.h"
#include "cdo_query.h"
#include "process_int.h"
#include "param_conversion.h"
#include "pmlist.h"
#include "cdo_default_values.h"

void *
Query(void *process)
{
  cdo_initialize(process);

  cdo_operator_add("query", 0, 0, "query entries");

  auto dataIsUnchanged = data_is_unchanged();

  auto operatorID = cdo_operator_id();

  operator_input_arg(cdo_operator_enter(operatorID));

  auto natts = cdo_operator_argc();
  if (natts == 0) cdo_abort("Parameter missing!");

  if (cdo_assert_files_only() == false) cdo_abort("This operator can't be combined with other operators!");

  PMList pmlist;
  KVList kvlist;
  kvlist.name = cdo_module_name();
  if (kvlist.parse_arguments(natts, cdo_get_oper_argv()) != 0) cdo_abort("Parse error!");
  if (Options::cdoVerbose) kvlist.print();

  auto pkvlist = &kvlist;
  if (natts == 1)
    {
      KeyValues &kv = kvlist.front();
      if (kv.key == "FILE")
        {
          if (Options::cdoVerbose) cdo_print("Reading query from: %s", kv.values[0]);
          auto filename = parameter_to_word(kv.values[0].c_str());
          auto fp = std::fopen(filename, "r");
          if (fp == nullptr) cdo_abort("Open failed on: %s\n", filename);
          pmlist.read_namelist(fp, filename);
          pkvlist = &pmlist.front();
          std::fclose(fp);
          if (Options::cdoVerbose) pkvlist->print();
        }
    }

  auto query = cdiQueryCreate();
  set_query_parameter(*pkvlist, query);
  if (Options::cdoVerbose) cdiQueryPrint(query);

  auto streamID1 = streamOpenReadQuery(cdo_get_stream_name(0), query);
  if (streamID1 < 0) cdi_open_error(streamID1, "Open failed on >%s<", cdo_get_stream_name(0));

  cdiQueryPrintEntriesNotFound(query);
  cdiQueryDelete(query);

  auto filetype = streamInqFiletype(streamID1);
  if (CdoDefault::FileType == CDI_UNDEFID) CdoDefault::FileType = filetype;

  auto vlistID1 = streamInqVlist(streamID1);
  auto vlistID2 = vlistDuplicate(vlistID1);

  auto taxisID1 = vlistInqTaxis(vlistID1);
  auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  auto streamID2 = cdo_open_write(1);
  cdo_def_vlist(streamID2, vlistID2);

  Field field;

  VarList varList1;
  varListInit(varList1, vlistID1);

  int tsID = 0;
  while (true)
    {
      auto nrecs = streamInqTimestep(streamID1, tsID);
      if (nrecs == 0) break;

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_def_timestep(streamID2, tsID);

      for (int recID = 0; recID < nrecs; ++recID)
        {
          int varID, levelID;
          streamInqRecord(streamID1, &varID, &levelID);
          cdo_def_record(streamID2, varID, levelID);
          /*
          if (dataIsUnchanged)
            {
              streamCopyRecord(streamID2, streamID1);
            }
          else
          */
          {
            field.init(varList1[varID]);
            if (field.memType == MemType::Float)
              streamReadRecordF(streamID1, field.vec_f.data(), &field.nmiss);
            else
              streamReadRecord(streamID1, field.vec_d.data(), &field.nmiss);
            cdo_write_record(streamID2, field);
          }
        }

      tsID++;
    }

  streamClose(streamID1);
  cdo_stream_close(streamID2);

  cdo_finish();

  return nullptr;
}
