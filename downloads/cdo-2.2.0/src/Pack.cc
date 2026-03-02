/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Pack    pack         Pack
*/

#ifdef _OPENMP
#include <omp.h>
#endif

#include <climits>

#include <cdi.h>

#include "process_int.h"
#include "cdo_vlist.h"
#include "datetime.h"
#include "cdo_default_values.h"
#include "field_functions.h"

static int
get_type_values(const int datatype, double &tmin, double &tmax, double &tmv)
{
  int status = 0;

  // clang-format off
  switch (datatype)
    {
    case CDI_DATATYPE_INT8:   tmv = -SCHAR_MAX; tmin = -SCHAR_MAX + 1;  tmax = SCHAR_MAX;     break;
    case CDI_DATATYPE_UINT8:  tmv =  UCHAR_MAX; tmin = 0;               tmax = UCHAR_MAX - 1; break;
    case CDI_DATATYPE_INT16:  tmv =  -SHRT_MAX; tmin = -SHRT_MAX + 1;   tmax = SHRT_MAX;      break;
    case CDI_DATATYPE_UINT16: tmv =  USHRT_MAX; tmin = 0;               tmax = USHRT_MAX - 1; break;
    case CDI_DATATYPE_INT32:  tmv =   -INT_MAX; tmin = -INT_MAX + 1;    tmax = INT_MAX;       break;
    case CDI_DATATYPE_UINT32: tmv =   UINT_MAX; tmin = 0;               tmax = UINT_MAX - 1;  break;
    default: status = 1; break;
    }
  // clang-format on

  return status;
}

static int
compute_scale_and_offset(const int datatype, const double fmin, const double fmax, double &scaleFactor, double &addOffset)
{
  scaleFactor = 1.0;
  addOffset = 0.0;

  double tmin, tmax, tmv;
  if (get_type_values(datatype, tmin, tmax, tmv)) return 1;

  if (is_not_equal(fmin, fmax))
    {
      scaleFactor = (fmax - fmin) / (tmax - tmin);
      addOffset = ((fmax + fmin) - scaleFactor * (tmin + tmax)) / 2;
    }

  return 0;
}

static MinMax
field_min_max(Field &field)
{
  auto nmiss = field.nmiss;
  auto missval = field.missval;
  auto len = field.size;

  if (field.memType == MemType::Float)
    return nmiss ? varray_min_max_mv(len, field.vec_f, (float) missval) : varray_min_max(len, field.vec_f);
  else
    return nmiss ? varray_min_max_mv(len, field.vec_d, missval) : varray_min_max(len, field.vec_d);
}

static void
fieldChangeMissval(Field &field, double missval1, double missval2)
{
  auto len = field.size;

  if (field.memType == MemType::Float)
    {
      auto &v = field.vec_f;
      for (size_t i = 0; i < len; ++i)
        if (dbl_is_equal(v[i], (float) missval1)) v[i] = missval2;
    }
  else
    {
      auto &v = field.vec_d;
      for (size_t i = 0; i < len; ++i)
        if (dbl_is_equal(v[i], missval1)) v[i] = missval2;
    }
}

void *
Pack(void *process)
{
  int datatype = CDI_DATATYPE_INT16;
  DateTimeList dtlist;

  cdo_initialize(process);

  operator_check_argc(0);

  auto streamID1 = cdo_open_read(0);

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto vlistID2 = vlistDuplicate(vlistID1);

  auto taxisID1 = vlistInqTaxis(vlistID1);
  auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  VarList varList;
  varListInit(varList, vlistID1);

  auto nvars = vlistNvars(vlistID1);

  FieldVector3D vars;

  int tsID = 0;
  while (true)
    {
      auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
      if (nrecs == 0) break;

      constexpr size_t NALLOC_INC = 1024;
      if ((size_t) tsID >= vars.size()) vars.resize(vars.size() + NALLOC_INC);

      dtlist.taxis_inq_timestep(taxisID1, tsID);

      fields_from_vlist(vlistID1, vars[tsID]);

      for (int recID = 0; recID < nrecs; ++recID)
        {
          int varID, levelID;
          cdo_inq_record(streamID1, &varID, &levelID);
          auto &field = vars[tsID][varID][levelID];
          field.init(varList[varID]);
          cdo_read_record(streamID1, field);
        }

      tsID++;
    }

  auto nts = tsID;

  if (CdoDefault::DataType != CDI_UNDEFID)
    {
      if (CdoDefault::DataType == CDI_DATATYPE_FLT64 || CdoDefault::DataType == CDI_DATATYPE_FLT32)
        {
          cdo_warning("Changed default output datatype to int16");
          CdoDefault::DataType = datatype;
        }
      else { datatype = CdoDefault::DataType; }
    }

  CdoDefault::DataType = datatype;
  constexpr double undefValue = 1.0e300;

  for (int varID = 0; varID < nvars; ++varID)
    {
      const auto &var = varList[varID];

      double fmin = undefValue, fmax = -undefValue;
      size_t nmisspv = 0;

      for (int levelID = 0; levelID < var.nlevels; ++levelID)
        {
          for (int t = 0; t < nts; ++t)
            {
              if (t > 0 && var.isConstant) continue;

              auto &field = vars[t][varID][levelID];
              auto nmiss = field.nmiss;

              if (nmiss) nmisspv += nmiss;

              if (nmiss < var.gridsize)
                {
                  auto mm = field_min_max(field);
                  fmin = std::min(fmin, mm.min);
                  fmax = std::max(fmax, mm.max);
                }
            }
        }

      vlistDefVarDatatype(vlistID2, varID, datatype);

      auto hasValidData = (is_not_equal(fmin, undefValue) && is_not_equal(fmax, -undefValue));

      if (nmisspv)
        {
          double tmin, tmax, missval2;
          if (!get_type_values(datatype, tmin, tmax, missval2))
            {
              vlistDefVarMissval(vlistID2, varID, missval2);

              if (!(missval2 < tmin || missval2 > tmax))
                cdo_warning("new missing value %g is inside data range (%g - %g)!", missval2, tmin, tmax);

              for (int levelID = 0; levelID < var.nlevels; ++levelID)
                {
                  for (int t = 0; t < nts; ++t)
                    {
                      if (t > 0 && var.isConstant) continue;

                      auto &field = vars[t][varID][levelID];
                      if (field.nmiss) fieldChangeMissval(field, var.missval, missval2);
                    }
                }
            }
        }

      if (hasValidData)
        {
          auto memTypeIsFloat = (var.memType == MemType::Float);
          double scaleFactor, addOffset;
          if (!compute_scale_and_offset(datatype, fmin, fmax, scaleFactor, addOffset))
            {
              cdiDefKeyFloat(vlistID2, varID, CDI_KEY_ADDOFFSET, memTypeIsFloat ? (float) addOffset : addOffset);
              cdiDefKeyFloat(vlistID2, varID, CDI_KEY_SCALEFACTOR, memTypeIsFloat ? (float) scaleFactor : scaleFactor);
            }
        }
    }

  auto streamID2 = cdo_open_write(1);
  cdo_def_vlist(streamID2, vlistID2);

  for (tsID = 0; tsID < nts; ++tsID)
    {
      dtlist.taxis_def_timestep(taxisID2, tsID);
      cdo_def_timestep(streamID2, tsID);

      for (int varID = 0; varID < nvars; ++varID)
        {
          const auto &var = varList[varID];
          if (tsID > 0 && var.isConstant) continue;
          for (int levelID = 0; levelID < var.nlevels; ++levelID)
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
