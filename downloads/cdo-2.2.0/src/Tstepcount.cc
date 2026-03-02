/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Tstepcount  tstepcount  Count number of timesteps
*/

#include <cdi.h>

#include "process_int.h"
#include "param_conversion.h"
#include "cdo_options.h"
#include "cimdOmp.h"
#include "field_functions.h"

template <typename T>
static T
tstepcount(long nts, T missval, const Varray<T> &v, T refval)
{
  if (DBL_IS_EQUAL(refval, missval)) return missval;

  long j;
  long n = 0;
  for (j = 0; j < nts; ++j)
    {
      n++;
      if (DBL_IS_EQUAL(v[j], refval)) break;
    }

  return (j == nts) ? missval : (T) n;
}

void *
Tstepcount(void *process)
{
  CdiDateTime vDateTime{};

  cdo_initialize(process);

  auto refval = (cdo_operator_argc() == 1) ? parameter_to_double(cdo_operator_argv(0)) : 0.0;

  const auto streamID1 = cdo_open_read(0);

  const auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  const auto vlistID2 = vlistDuplicate(vlistID1);

  vlistDefNtsteps(vlistID2, 1);

  const auto nvars = vlistNvars(vlistID1);
  for (int varID = 0; varID < nvars; ++varID) cdiDefKeyString(vlistID2, varID, CDI_KEY_UNITS, "steps");

  const auto taxisID1 = vlistInqTaxis(vlistID1);
  const auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  const auto streamID2 = cdo_open_write(1);
  cdo_def_vlist(streamID2, vlistID2);

  VarList varList1;
  varListInit(varList1, vlistID1);

  FieldVector3D vars;

  int tsID = 0;
  while (true)
    {
      const auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
      if (nrecs == 0) break;

      constexpr size_t NALLOC_INC = 1024;
      if ((size_t) tsID >= vars.size()) vars.resize(vars.size() + NALLOC_INC);

      vDateTime = taxisInqVdatetime(taxisID1);

      fields_from_vlist(vlistID1, vars[tsID]);

      for (int recID = 0; recID < nrecs; ++recID)
        {
          int varID, levelID;
          cdo_inq_record(streamID1, &varID, &levelID);
          auto &field = vars[tsID][varID][levelID];
          field.init(varList1[varID]);
          cdo_read_record(streamID1, field);
        }

      tsID++;
    }

  int nts = tsID;

  std::vector<Field> fields(Threading::ompNumThreads);

  for (int varID = 0; varID < nvars; ++varID)
    {
      auto memType = varList1[varID].memType;
      auto missval = varList1[varID].missval;
      auto gridsize = varList1[varID].gridsize;
      for (int levelID = 0; levelID < varList1[varID].nlevels; ++levelID)
        {
#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(dynamic, 1)
#endif
          for (size_t i = 0; i < gridsize; ++i)
            {
              const auto ompthID = cdo_omp_get_thread_num();

              if (memType == MemType::Float)
                {
                  auto &v = fields[ompthID].vec_f;
                  v.resize(nts);
                  for (int t = 0; t < nts; ++t) v[t] = vars[t][varID][levelID].vec_f[i];

                  const auto count = tstepcount(nts, (float) missval, v, (float) refval);

                  vars[0][varID][levelID].vec_f[i] = count;
                }
              else
                {
                  auto &v = fields[ompthID].vec_d;
                  v.resize(nts);
                  for (int t = 0; t < nts; ++t) v[t] = vars[t][varID][levelID].vec_d[i];

                  const auto count = tstepcount(nts, missval, v, refval);

                  vars[0][varID][levelID].vec_d[i] = count;
                }
            }
        }
    }

  taxisDefVdatetime(taxisID2, vDateTime);
  cdo_def_timestep(streamID2, 0);

  for (int varID = 0; varID < nvars; ++varID)
    {
      for (int levelID = 0; levelID < varList1[varID].nlevels; ++levelID)
        {
          cdo_def_record(streamID2, varID, levelID);
          auto &field1 = vars[0][varID][levelID];
          field_num_mv(field1);
          cdo_write_record(streamID2, field1);
        }
    }

  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  cdo_finish();

  return nullptr;
}
