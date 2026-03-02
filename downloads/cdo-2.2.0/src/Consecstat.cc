/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: 2010 Ralf Mueller

*/

/*
   This module contains the following operators:

      Consectstep  consecsum  For each timestep, the current number of
                              onsecutive timsteps is counted
      Consectstep  consects   For each period of consecutive timesteps, only
                              count its length + last contributing timesteps

   =============================================================================
   Created:  04/08/2010 11:58:01 AM
    Author:  Ralf Mueller (ram), ralf.mueller@mpimet.mpg.de
   Company:  Max-Planck-Institute for Meteorology
   =============================================================================
*/

#include <cdi.h>

#include "process_int.h"
#include "cdo_vlist.h"
#include "param_conversion.h"
#include "field_functions.h"

enum
{
  CONSECSUM,
  CONSECTS
};

#define SWITCHWARN "Hit default case! This should never happen (%s).\n"

static void
selEndOfPeriod(Field &periods, const Field &history, const Field &current, int isLastTimestep)
{
  auto pmissval = periods.missval;
  auto &parray = periods.vec_d;
  auto hmissval = history.missval;
  const auto &harray = history.vec_d;
  auto cmissval = current.missval;
  const auto &carray = current.vec_d;

  auto len = gridInqSize(periods.grid);
  if (len != gridInqSize(current.grid) || (gridInqSize(current.grid) != gridInqSize(history.grid)))
    cdo_abort("Fields have different gridsize (%s)", __func__);

  if (!isLastTimestep)
    {
      if (current.nmiss || history.nmiss)
        {
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
          for (size_t i = 0; i < len; ++i)
            {
              if (!DBL_IS_EQUAL(harray[i], hmissval))
                {
                  if (!DBL_IS_EQUAL(carray[i], cmissval))
                    parray[i] = (DBL_IS_EQUAL(carray[i], 0.0) && IS_NOT_EQUAL(harray[i], 0.0)) ? harray[i] : pmissval;
                  else  // DBL_IS_EQUAL(carray[i], cmissval)
                    parray[i] = (IS_NOT_EQUAL(harray[i], 0.0)) ? harray[i] : pmissval;
                }
              else /* DBL_IS_EQUAL(harray[i], hmissval) */ { parray[i] = pmissval; }
            }
        }
      else
        {
#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
          for (size_t i = 0; i < len; ++i)
            parray[i] = (DBL_IS_EQUAL(carray[i], 0.0) && IS_NOT_EQUAL(harray[i], 0.0)) ? harray[i] : pmissval;
        }
    }
  else
    {
      if (current.nmiss)
        {
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
          for (size_t i = 0; i < len; ++i)
            if (!DBL_IS_EQUAL(carray[i], cmissval))
              parray[i] = (DBL_IS_EQUAL(carray[i], 0.0)) ? pmissval : carray[i];
            else  // DBL_IS_EQUAL(carray[i], cmissval)
              parray[i] = pmissval;
        }
      else
        {
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
          for (size_t i = 0; i < len; ++i) parray[i] = DBL_IS_EQUAL(carray[i], 0.0) ? pmissval : carray[i];
        }
    }

  periods.nmiss = varray_num_mv(len, parray, pmissval);
}

void *
Consecstat(void *process)
{
  CdiDateTime vDateTime{};
  CdiDateTime histDateTime{};
  double refval = 0.0;

  cdo_initialize(process);
  cdo_operator_add("consecsum", CONSECSUM, 0, "refval");
  cdo_operator_add("consects", CONSECTS, 0, nullptr);
  const int operatorID = cdo_operator_id();

  if (operatorID == CONSECSUM)
    if (cdo_operator_argc() > 0) refval = parameter_to_double(cdo_operator_argv(0));

  const auto istreamID = cdo_open_read(0);

  const auto ivlistID = cdo_stream_inq_vlist(istreamID);
  const auto itaxisID = vlistInqTaxis(ivlistID);
  const auto ovlistID = vlistDuplicate(ivlistID);
  const auto otaxisID = taxisDuplicate(itaxisID);
  vlistDefTaxis(ovlistID, otaxisID);

  VarList varList;
  varListInit(varList, ovlistID);

  Field field;
  field.resize(vlistGridsizeMax(ovlistID));

  const auto nvars = vlistNvars(ivlistID);
  FieldVector2D vars, hist, periods;
  fields_from_vlist(ivlistID, vars, FIELD_VEC, 0);
  if (operatorID == CONSECTS)
    {
      fields_from_vlist(ivlistID, hist, FIELD_VEC);
      fields_from_vlist(ivlistID, periods, FIELD_VEC);
    }

  for (int varID = 0; varID < nvars; ++varID) cdiDefKeyString(ovlistID, varID, CDI_KEY_UNITS, "steps");  // TODO

  const auto ostreamID = cdo_open_write(1);
  cdo_def_vlist(ostreamID, ovlistID);

  int itsID = 0;
  int otsID = 0;
  while (true)
    {
      const auto nrecs = cdo_stream_inq_timestep(istreamID, itsID);
      if (nrecs == 0) break;

      vDateTime = taxisInqVdatetime(itaxisID);
      switch (operatorID)
        {
        case CONSECSUM:
          taxisDefVdatetime(otaxisID, vDateTime);
          cdo_def_timestep(ostreamID, otsID);
          break;
        case CONSECTS:
          if (itsID != 0)
            {
              taxisDefVdatetime(otaxisID, histDateTime);
              cdo_def_timestep(ostreamID, otsID - 1);
            }
          break;
        default: printf(SWITCHWARN, __func__); break;
        }

      for (int recID = 0; recID < nrecs; ++recID)
        {
          int varID, levelID;
          cdo_inq_record(istreamID, &varID, &levelID);
          cdo_read_record(istreamID, field.vec_d.data(), &field.nmiss);
          field.grid = varList[varID].gridID;
          field.size = varList[varID].gridsize;
          field.missval = varList[varID].missval;

          field2_sumtr(vars[varID][levelID], field, refval);

          switch (operatorID)
            {
            case CONSECSUM:
              cdo_def_record(ostreamID, varID, levelID);
              cdo_write_record(ostreamID, vars[varID][levelID].vec_d.data(), vars[varID][levelID].nmiss);
              break;
            case CONSECTS:
              if (itsID != 0)
                {
                  selEndOfPeriod(periods[varID][levelID], hist[varID][levelID], vars[varID][levelID], false);
                  cdo_def_record(ostreamID, varID, levelID);
                  cdo_write_record(ostreamID, periods[varID][levelID].vec_d.data(), periods[varID][levelID].nmiss);
                }
#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
              for (size_t i = 0; i < vars[varID][levelID].size; ++i) hist[varID][levelID].vec_d[i] = vars[varID][levelID].vec_d[i];
#else
              hist[varID][levelID].vec_d = vars[varID][levelID].vec_d;
#endif
              break;
            default: printf(SWITCHWARN, __func__); break;
            }
        }

      histDateTime = vDateTime;
      itsID++;
      otsID++;
    }

  if (operatorID == CONSECTS) /* Save the last timestep */
    {
      taxisDefVdatetime(otaxisID, vDateTime);
      cdo_def_timestep(ostreamID, otsID - 1);

      for (int varID = 0; varID < nvars; ++varID)
        {
          const auto nlevels = varList[varID].nlevels;
          for (int levelID = 0; levelID < nlevels; ++levelID)
            {
              selEndOfPeriod(periods[varID][levelID], hist[varID][levelID], vars[varID][levelID], true);
              cdo_def_record(ostreamID, varID, levelID);
              cdo_write_record(ostreamID, periods[varID][levelID].vec_d.data(), periods[varID][levelID].nmiss);
            }
        }
    }

  cdo_stream_close(istreamID);
  cdo_stream_close(ostreamID);

  cdo_finish();

  return nullptr;
}
