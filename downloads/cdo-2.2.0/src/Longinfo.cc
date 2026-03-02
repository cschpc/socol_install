/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Longinfo       linfo            Long dataset information
*/

#include <cfloat>

#include <cdi.h>

#include "cdo_options.h"
#include "cdo_vlist.h"
#include "process_int.h"
#include "varray.h"
#include "datetime.h"
#include "printinfo.h"
#include "cdo_zaxis.h"
#include "field_functions.h"

void cdoPrintAttributes(FILE *fp, int cdiID, int varID, int nblanks);

struct LonginfoStat
{
  double min = DBL_MAX;
  double max = -DBL_MAX;
  double sum = 0.0;
  double sumi = 0.0;
  size_t nvals = 0;
};

static void
field_min_max_sum(const Field &field, double &min, double &max, double &sum)
{
  auto mms = MinMaxSum(min, max, sum);
  if (field.memType == MemType::Float)
    mms = varray_min_max_sum(field.size, field.vec_f, mms);
  else
    mms = varray_min_max_sum(field.size, field.vec_d, mms);

  min = mms.min;
  max = mms.max;
  sum = mms.sum;
}

size_t static field_min_max_sum_mv(const Field &field, double &min, double &max, double &sum)
{
  auto mms = MinMaxSum(min, max, sum);
  if (field.memType == MemType::Float)
    mms = varray_min_max_sum_mv(field.size, field.vec_f, static_cast<float>(field.missval), mms);
  else
    mms = varray_min_max_sum_mv(field.size, field.vec_d, field.missval, mms);

  min = mms.min;
  max = mms.max;
  sum = mms.sum;
  return mms.n;
}

static size_t
compute_stat_real(const Field &field, LonginfoStat &infostat, size_t gridsize)
{
  size_t imiss = 0;

  if (field.nmiss)
    {
      auto nvals = field_min_max_sum_mv(field, infostat.min, infostat.max, infostat.sum);
      imiss = gridsize - nvals;
      infostat.nvals += nvals;
    }
  else if (gridsize == 1)
    {
      double val = (field.memType == MemType::Float) ? field.vec_f[0] : field.vec_d[0];
      infostat.sum = (infostat.nvals == 0) ? val : infostat.sum + val;
      infostat.nvals += 1;
    }
  else
    {
      field_min_max_sum(field, infostat.min, infostat.max, infostat.sum);
      infostat.nvals += gridsize;
    }

  return imiss;
}

void *
Longinfo(void *process)
{
  char paramstr[32];

  cdo_initialize(process);

  // clang-format off
  cdo_operator_add("linfo",   0,  0, nullptr);
  // clang-format on

  operator_check_argc(0);

  auto streamID = cdo_open_read(0);
  auto vlistID = cdo_stream_inq_vlist(streamID);
  auto taxisID = vlistInqTaxis(vlistID);

  VarList varList;
  varListInit(varList, vlistID);

  Field field;

  int tsID = 0;
  while (true)
    {
      auto nrecs = cdo_stream_inq_timestep(streamID, tsID);
      if (nrecs == 0) break;

      fprintf(stdout, "timestep: %d\n", tsID + 1);

      auto vDateTime = taxisInqVdatetime(taxisID);
      fprintf(stdout, "\tdateTime: %s\n\n", datetime_to_string(vDateTime).c_str());

      for (int recID = 0; recID < nrecs; ++recID)
        {
          fprintf(stdout, "\tfield: %d of %d\n", recID + 1, nrecs);

          int varID, levelID;
          cdo_inq_record(streamID, &varID, &levelID);
          const auto &var = varList[varID];

          auto dig = (var.datatype == CDI_DATATYPE_FLT64) ? Options::CDO_dbl_digits : Options::CDO_flt_digits;

          fprintf(stdout, "\t\tvarIndex: %d\n", varID + 1);
          fprintf(stdout, "\t\tlevelIndex: %d\n", levelID + 1);
          fprintf(stdout, "\t\tlevel: %.*g\n", dig, cdo_zaxis_inq_level(var.zaxisID, levelID));
          fprintf(stdout, "\t\tname: %s\n", var.name.c_str());
          if (var.longname.size()) fprintf(stdout, "\t\tlongname: \"%s\"\n", var.longname.c_str());
          if (var.units.size()) fprintf(stdout, "\t\tunits: \"%s\"\n", var.units.c_str());
          cdiParamToString(var.param, paramstr, sizeof(paramstr));
          if (paramstr[0] && paramstr[0] != '-') fprintf(stdout, "\t\tparam: %s\n", paramstr);

          field.init(var);
          cdo_read_record(streamID, field);
          auto nmiss = field.nmiss;

          LonginfoStat infostat;

          fprintf(stdout, "\t\tdataType: %s\n", cdo::datatype_to_cstr(var.datatype));
          fprintf(stdout, "\t\tmemoryType: %s\n", (var.memType == MemType::Float) ? "float" : "double");
          fprintf(stdout, "\t\tgridsize: %zu\n", var.gridsize);
          fprintf(stdout, "\t\tnumMiss: %zu\n", nmiss);
          fprintf(stdout, "\t\tmissval: %.*g\n", dig, var.missval);

          double addoffset = 0.0, scalefactor = 1.0;
          auto haveAddoffset = (cdiInqKeyFloat(vlistID, varID, CDI_KEY_ADDOFFSET, &addoffset) == CDI_NOERR);
          auto haveScalefactor = (cdiInqKeyFloat(vlistID, varID, CDI_KEY_SCALEFACTOR, &scalefactor) == CDI_NOERR);
          if (haveAddoffset) fprintf(stdout, "\t\taddoffset: %.*g\n", dig, addoffset);
          if (haveScalefactor) fprintf(stdout, "\t\tscalefactor: %.*g\n", dig, scalefactor);

          auto imiss = compute_stat_real(field, infostat, var.gridsize);

          if (infostat.nvals > 1)
            {
              fprintf(stdout, "\t\trange: %.*g\n", dig, infostat.max - infostat.min);
              fprintf(stdout, "\t\tminimum: %.*g\n", dig, infostat.min);
              fprintf(stdout, "\t\tmaximum: %.*g\n", dig, infostat.max);
              fprintf(stdout, "\t\taverage: %.*g\n", dig, infostat.sum / infostat.nvals);
              fprintf(stdout, "\t\tstandardDev: %.*g\n", dig, field_std1(field));
              fprintf(stdout, "\t\tskewness: %.*g\n", dig, field_skew(field));
              fprintf(stdout, "\t\tkurtosis: %.*g\n", dig, field_kurt(field));
            }
          else if (infostat.nvals == 1) { fprintf(stdout, "\t\tvalue: %g\n", infostat.sum); }

          cdoPrintAttributes(stdout, vlistID, varID, 16);
        }

      // if (imiss != nmiss && nmiss) cdo_warning("Found %zu of %zu missing values!", imiss, nmiss);

      tsID++;
    }

  cdo_stream_close(streamID);

  cdo_finish();

  return nullptr;
}
