/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

     Timsort    timsort         Sort over the time
*/

#include <algorithm>  // sort

#include <cdi.h>

#include "process_int.h"
#include "cdo_vlist.h"
#include "cdo_options.h"
#include "cimdOmp.h"
#include "field_functions.h"

void *
Timsort(void *process)
{
  int varID, levelID;
  int nalloc = 0;

  cdo_initialize(process);

  operator_check_argc(0);

  auto streamID1 = cdo_open_read(0);

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto vlistID2 = vlistDuplicate(vlistID1);

  auto taxisID1 = vlistInqTaxis(vlistID1);
  auto taxisID2 = cdo_taxis_create(TAXIS_ABSOLUTE);
  vlistDefTaxis(vlistID2, taxisID2);

  auto streamID2 = cdo_open_write(1);
  cdo_def_vlist(streamID2, vlistID2);

  VarList varList;
  varListInit(varList, vlistID1);

  auto nvars = vlistNvars(vlistID1);
  FieldVector3D vars;
  std::vector<CdiDateTime> vDateTimes;

  int tsID = 0;
  while (true)
    {
      auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
      if (nrecs == 0) break;

      if (tsID >= nalloc)
        {
          constexpr int NALLOC_INC = 1024;
          nalloc += NALLOC_INC;
          vDateTimes.resize(nalloc);
          vars.resize(nalloc);
        }

      vDateTimes[tsID] = taxisInqVdatetime(taxisID1);

      fields_from_vlist(vlistID1, vars[tsID]);

      for (int recID = 0; recID < nrecs; ++recID)
        {
          cdo_inq_record(streamID1, &varID, &levelID);
          auto &field = vars[tsID][varID][levelID];
          field.init(varList[varID]);
          cdo_read_record(streamID1, field);
        }

      tsID++;
    }

  int nts = tsID;

  std::vector<Field> fields(Threading::ompNumThreads);

  for (varID = 0; varID < nvars; ++varID)
    {
      const auto &var = varList[varID];

      if (var.isConstant) continue;

      auto memType = var.memType;
      auto gridsize = var.gridsize;
      for (levelID = 0; levelID < var.nlevels; ++levelID)
        {
#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
          for (size_t i = 0; i < gridsize; ++i)
            {
              auto ompthID = cdo_omp_get_thread_num();

              if (memType == MemType::Float)
                {
                  auto &v = fields[ompthID].vec_f;
                  v.resize(nts);
                  for (int t = 0; t < nts; ++t) v[t] = vars[t][varID][levelID].vec_f[i];

                  std::sort(v.begin(), v.end());

                  for (int t = 0; t < nts; ++t) vars[t][varID][levelID].vec_f[i] = v[t];
                }
              else
                {
                  auto &v = fields[ompthID].vec_d;
                  v.resize(nts);
                  for (int t = 0; t < nts; ++t) v[t] = vars[t][varID][levelID].vec_d[i];

                  std::sort(v.begin(), v.end());

                  for (int t = 0; t < nts; ++t) vars[t][varID][levelID].vec_d[i] = v[t];
                }
            }
        }
    }

  for (tsID = 0; tsID < nts; ++tsID)
    {
      taxisDefVdatetime(taxisID2, vDateTimes[tsID]);
      cdo_def_timestep(streamID2, tsID);

      for (varID = 0; varID < nvars; ++varID)
        {
          const auto &var = varList[varID];
          for (levelID = 0; levelID < var.nlevels; ++levelID)
            {
              auto &field = vars[tsID][varID][levelID];
              if (field.hasData())
                {
                  cdo_def_record(streamID2, varID, levelID);
                  cdo_write_record(streamID2, field);
                }
            }
        }
    }

  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  cdo_finish();

  return nullptr;
}
