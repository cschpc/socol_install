/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Sinfo      sinfo           Short dataset information
*/

#include "cdi.h"
#include "julian_date.h"

#include <cstring>
#include "cdo_options.h"
#include "printinfo.h"
#include "mpmo_color.h"
#include "process_int.h"
#include "compare.h"
#include "util_string.h"
#include "datetime.h"
#include "cdo_default_values.h"

const char *get_steptype_name(const int tsteptype);

enum
{
  func_generic,
  func_param,
  func_name,
  func_code
};

static const char *
memtype_to_cstr(MemType memType)
{
  return (memType == MemType::Double) ? "F64" : "F32";
}

static const char *
num_values_to_byte_cstr(size_t numValues)
{
  static char cstring[128];
  cstring[0] = 0;

  size_t muindex = 0;
  const char *mu[] = { "Bytes", "Kbytes", "Mbytes", "Gbytes", "Tbytes", "Pbytes" };
  const size_t nmu = sizeof(mu) / sizeof(char *);
  while (numValues > 9999 && muindex < nmu - 1)
    {
      numValues /= 1024;
      muindex++;
    }
  std::snprintf(cstring, sizeof(cstring), "%zu %s", numValues, mu[muindex]);

  return cstring;
}

static size_t
get_num_input_bits(int datatype)
{
  // clang-format off
  if      (datatype == CDI_DATATYPE_INT8  ) return 8;
  else if (datatype == CDI_DATATYPE_UINT8 ) return 8;
  else if (datatype == CDI_DATATYPE_INT16 ) return 16;
  else if (datatype == CDI_DATATYPE_UINT16) return 16;
  else if (datatype == CDI_DATATYPE_INT32 ) return 32;
  else if (datatype == CDI_DATATYPE_UINT32) return 32;
  else if (datatype == CDI_DATATYPE_FLT32 ) return 32;
  else if (datatype == CDI_DATATYPE_FLT64 ) return 64;
  else if (datatype == CDI_DATATYPE_PACK8 ) return 8;
  else if (datatype == CDI_DATATYPE_PACK16) return 16;
  else if (datatype == CDI_DATATYPE_PACK32) return 24;
  else if (datatype == CDI_DATATYPE_PACK  ) return 8;        // unknown
  else if (datatype > 0 && datatype <= 32 ) return datatype; // Number of packed bits in GRIB
  else                                      return 64;
  // clang-format on
}

static size_t
get_num_output_bits(int datatype)
{
  if (CdoDefault::DataType == CDI_UNDEFID)
    {
      if (CdoDefault::FileType != CDI_UNDEFID) {}
    }
  else
    {
      datatype = CdoDefault::DataType;
    }

  return get_num_input_bits(datatype);
}

static void
limit_string_length(char *string, size_t maxlen)
{
  string[maxlen - 1] = 0;
  auto len = strlen(string);
  if (len > 10)
    {
      for (size_t i = 3; i < len; ++i)
        if (string[i] == ' ' || string[i] == ',' || (i > 10 && string[i] == '.'))
          {
            string[i] = 0;
            break;
          }
    }
}

static void
print_vars_info(int operfunc, bool ensembleInfo, int vlistID, bool xsInfo)
{
  char tmpname[CDI_MAX_NAME];
  char paramstr[32];

  VarList varList;
  varListInit(varList, vlistID);

  auto nvars = vlistNvars(vlistID);
  auto nsubtypes = vlistNsubtypes(vlistID);

  for (int varID = 0; varID < nvars; ++varID)
    {
      const auto &var = varList[varID];

      auto tabnum = tableInqNum(vlistInqVarTable(vlistID, varID));

      fprintf(stdout, "%6d", varID + 1);
      fprintf(stdout, " : ");

      set_text_color(stdout, BLUE);
      // institute info
      auto instptr = institutInqNamePtr(vlistInqVarInstitut(vlistID, varID));
      strcpy(tmpname, "unknown");
      if (instptr) strncpy(tmpname, instptr, CDI_MAX_NAME - 1);
      limit_string_length(tmpname, 32);
      fprintf(stdout, "%-8s ", tmpname);

      // source info
      auto modelptr = modelInqNamePtr(vlistInqVarModel(vlistID, varID));
      strcpy(tmpname, "unknown");
      if (modelptr) strncpy(tmpname, modelptr, CDI_MAX_NAME - 1);
      limit_string_length(tmpname, 32);
      fprintf(stdout, "%-8s ", tmpname);

      // timetype
      fprintf(stdout, "%c ", var.isConstant ? 'c' : 'v');

      // tsteptype
      fprintf(stdout, "%-8s ", get_steptype_name(var.tsteptype));

      // ensemble information
      if (ensembleInfo)
        {
          int perturbationNumber, numberOfForecastsInEnsemble;
          auto r1 = cdiInqKeyInt(vlistID, varID, CDI_KEY_PERTURBATIONNUMBER, &perturbationNumber);
          auto r2 = cdiInqKeyInt(vlistID, varID, CDI_KEY_NUMBEROFFORECASTSINENSEMBLE, &numberOfForecastsInEnsemble);
          if (r1 == 0 && r2 == 0)
            fprintf(stdout, "%2d/%-2d ", perturbationNumber, numberOfForecastsInEnsemble);
          else
            fprintf(stdout, "--/-- ");
        }

      if (nsubtypes > 1)
        {
          auto subtypeID = vlistInqVarSubtype(vlistID, varID);
          auto subtypesize = subtypeInqSize(subtypeID);
          fprintf(stdout, " %6d  ", subtypesize);
          fprintf(stdout, "%3d ", vlistSubtypeIndex(vlistID, subtypeID) + 1);
        }
      reset_text_color(stdout);

      // layer info
      set_text_color(stdout, GREEN);
      fprintf(stdout, "%6d ", var.nlevels);
      reset_text_color(stdout);
      fprintf(stdout, "%3d ", vlistZaxisIndex(vlistID, var.zaxisID) + 1);

      // grid info
      set_text_color(stdout, GREEN);
      fprintf(stdout, "%9zu ", var.gridsize);
      reset_text_color(stdout);
      fprintf(stdout, "%3d ", vlistGridIndex(vlistID, var.gridID) + 1);

      // datatype
      set_text_color(stdout, BLUE);
      fprintf(stdout, " %-3s", cdo::datatype_to_cstr(var.datatype));

      // memType
      if (xsInfo) fprintf(stdout, "   %-3s", memtype_to_cstr(var.memType));

      auto compType = vlistInqVarCompType(vlistID, varID);
      auto isCompressed = (compType != CDI_COMPRESS_NONE);
      fprintf(stdout, "%c ", isCompressed ? (int) comptype_to_name(compType)[0] : ' ');

      reset_text_color(stdout);

      fprintf(stdout, ": ");

      // parameter info
      cdiParamToString(var.param, paramstr, sizeof(paramstr));

      // set_text_color(stdout, GREEN);
      // clang-format off
      if      (operfunc == func_name) fprintf(stdout, "%-14s", var.name.c_str());
      else if (operfunc == func_code) fprintf(stdout, "%4d %4d", tabnum, var.code);
      else                            fprintf(stdout, "%-14s", paramstr);
      // clang-format on
      if (xsInfo && Options::cdoVerbose && operfunc == func_name && var.units.size()) fprintf(stdout, " [%s]", var.units.c_str());
      // reset_text_color(stdout);

      if (Options::cdoVerbose)
        {
          auto chunks = cdo::inq_key_string(vlistID, varID, CDI_KEY_CHUNKS);
          fprintf(stdout, " : %s", chunks.c_str());
        }

      fprintf(stdout, "\n");
    }
}

static void
print_time_info(int ntsteps, int taxisID)
{
  set_text_color(stdout, BRIGHT);
  fprintf(stdout, "   Time coordinate");
  reset_text_color(stdout);
  fprintf(stdout, " :\n");

  auto taxisName = taxisNamePtr(taxisID);
  auto tname = taxisName ? taxisName : "time";
  fprintf(stdout, "%33s : ", tname);

  set_text_color(stdout, GREEN);
  if (ntsteps == CDI_UNDEFID)
    fprintf(stdout, "unlimited steps\n");
  else
    fprintf(stdout, "%d step%s\n", ntsteps, (ntsteps == 1) ? "" : "s");
  reset_text_color(stdout);

  if (Options::cdoVerbose)
    {
      // int datatype;
      // cdiInqKeyInt(taxisID, CDI_GLOBAL, CDI_KEY_DATATYPE, &datatype);
      // fprintf(stdout, "%33s : %s\n", "datatype", cdo::datatype_to_cstr(datatype));
      fprintf(stdout, "%33s : %d\n", "taxisID", taxisID);
    }

  if (taxisID != CDI_UNDEFID)
    {
      if (taxisInqType(taxisID) != TAXIS_ABSOLUTE)
        {
          auto rDateTime = taxisInqRdatetime(taxisID);
          fprintf(stdout, "     RefTime = %s %s", date_to_string(rDateTime.date).c_str(), time_to_string(rDateTime.time).c_str());

          auto tunits = taxisInqTunit(taxisID);
          if (tunits != CDI_UNDEFID) fprintf(stdout, "  Units = %s", tunit_to_cstr(tunits));

          auto calendar = taxisInqCalendar(taxisID);
          if (calendar != CDI_UNDEFID) fprintf(stdout, "  Calendar = %s", calendar_to_cstr(calendar));

          if (taxisHasBounds(taxisID)) fprintf(stdout, "  Bounds = true");

          fprintf(stdout, "\n");

          if (taxisInqType(taxisID) == TAXIS_FORECAST)
            {
              auto fDateTime = taxisInqFdatetime(taxisID);
              fprintf(stdout, "     ForecastRefTime = %s\n", datetime_to_string(fDateTime).c_str());
            }
        }
    }
}

static int
print_time_info_xs(int ntsteps, int taxisID, CdoStreamID streamID)
{
  int numTimesteps = ntsteps;
  TimeIncrement timeIncrement, timeIncrement0;
  CdiDateTime vDateTimeFirst{};
  CdiDateTime vDateTimeLast{};

  set_text_color(stdout, BRIGHT);
  fprintf(stdout, "   Time coordinate");
  reset_text_color(stdout);
  fprintf(stdout, " :\n");

  if (taxisID != CDI_UNDEFID)
    {
      auto calendar = taxisInqCalendar(taxisID);

      int tsID = 0;
      while (true)
        {
          auto nrecs = cdo_stream_inq_timestep(streamID, tsID);
          if (nrecs == 0) break;

          auto vDateTime = taxisInqVdatetime(taxisID);

          if (tsID)
            {
              auto julianDate0 = julianDate_encode(calendar, vDateTimeLast);
              auto julianDate = julianDate_encode(calendar, vDateTime);
              auto jdelta = julianDate_to_seconds(julianDate_sub(julianDate, julianDate0));
              timeIncrement = get_time_increment(jdelta, vDateTimeLast.date, vDateTime.date);
            }
          else
            {
              vDateTimeFirst = vDateTime;
            }

          if (tsID == 1) { timeIncrement0 = timeIncrement; }
          else if (tsID > 1 && timeIncrement0 != timeIncrement)
            {
              timeIncrement0.period = 0;
            }

          vDateTimeLast = vDateTime;

          tsID++;
        }

      numTimesteps = tsID;
    }

  auto taxisName = taxisNamePtr(taxisID);
  auto tname = taxisName ? taxisName : "time";

  fprintf(stdout, "%33s : ", "steps");

  set_text_color(stdout, GREEN);
  if (ntsteps == CDI_UNDEFID)
    fprintf(stdout, "unlimited (%d currently)\n", numTimesteps);
  else
    fprintf(stdout, "%d\n", ntsteps);
  reset_text_color(stdout);

  if (taxisID != CDI_UNDEFID)
    {
      fprintf(stdout, "%33s : ", tname);

      fprintf(stdout, "%s", datetime_to_string(vDateTimeFirst).c_str());
      if (numTimesteps > 1)
        {
          fprintf(stdout, " to %s", datetime_to_string(vDateTimeLast).c_str());
          if (timeIncrement0.period)
            fprintf(stdout, " by %d %s%s", (int) timeIncrement0.period, time_units_cstr(timeIncrement0.units),
                    (timeIncrement0.period != 1) ? "s" : "");
        }
      fprintf(stdout, "\n");

      auto calendar = taxisInqCalendar(taxisID);

      if (taxisInqType(taxisID) != TAXIS_ABSOLUTE)
        {
          fprintf(stdout, "%33s : ", "units");

          auto tunits = taxisInqTunit(taxisID);
          fprintf(stdout, "%s", tunit_to_cstr(tunits));

          auto rDateTime = taxisInqRdatetime(taxisID);
          fprintf(stdout, " since %s\n", datetime_to_string(rDateTime).c_str());

          if (calendar != CDI_UNDEFID) fprintf(stdout, "%33s : %s\n", "calendar", calendar_to_cstr(calendar));

          if (taxisInqType(taxisID) == TAXIS_FORECAST)
            {
              auto fDateTime = taxisInqFdatetime(taxisID);
              fprintf(stdout, "%33s : %s\n", "forecastRefTime", datetime_to_string(fDateTime).c_str());
            }
        }

      if (taxisHasBounds(taxisID)) fprintf(stdout, "%33s : %s\n", "available", "bounds");
    }

  return numTimesteps;
}

static void
addOperators(void)
{
  // clang-format off
  cdo_operator_add("sinfo",   func_generic, 0, nullptr);
  cdo_operator_add("sinfop",  func_param,   0, nullptr);
  cdo_operator_add("sinfon",  func_name,    0, nullptr);
  cdo_operator_add("sinfoc",  func_code,    0, nullptr);
  cdo_operator_add("seinfo",  func_generic, 1, nullptr);
  cdo_operator_add("seinfop", func_param,   1, nullptr);
  cdo_operator_add("seinfon", func_name,    1, nullptr);
  cdo_operator_add("seinfoc", func_code,    1, nullptr);
  cdo_operator_add("xsinfo",  func_name,    2, nullptr);
  cdo_operator_add("xsinfop", func_param,   2, nullptr);
  cdo_operator_add("xsinfon", func_name,    2, nullptr);
  cdo_operator_add("xsinfoc", func_code,    2, nullptr);
  // clang-format on
}

class ModuleSinfo
{
private:
  int operfunc;
  int ensembleInfo;
  bool xsInfo;

public:
  void
  init(void *process)
  {
    cdo_initialize(process);

    addOperators();

    auto operatorID = cdo_operator_id();

    operfunc = cdo_operator_f1(operatorID);
    ensembleInfo = (cdo_operator_f2(operatorID) == 1);
    xsInfo = (cdo_operator_f2(operatorID) == 2);

    operator_check_argc(0);
  }

  void
  run()
  {
    for (int indf = 0; indf < cdo_stream_cnt(); indf++)
      {
        auto streamID = cdo_open_read(indf);
        auto vlistID = cdo_stream_inq_vlist(streamID);
        auto nsubtypes = vlistNsubtypes(vlistID);

        set_text_color(stdout, BRIGHT);
        fprintf(stdout, "   File format");
        reset_text_color(stdout);
        fprintf(stdout, " : ");
        print_filetype(streamID, vlistID);

        set_text_color(stdout, BRIGHT);
        fprintf(stdout, "%6d : Institut Source   T Steptype", -(indf + 1));
        if (ensembleInfo) fprintf(stdout, " Einfo");
        if (nsubtypes > 1) fprintf(stdout, " Subtypes");
        fprintf(stdout, " Levels Num    Points Num Dtype");
        if (xsInfo) fprintf(stdout, " Mtype");
        fprintf(stdout, " : %s",
                (operfunc == func_name) ? "Parameter name" : ((operfunc == func_code) ? "Table Code" : "Parameter ID"));

        if (Options::cdoVerbose) fprintf(stdout, " : Extra");
        reset_text_color(stdout);
        fprintf(stdout, "\n");

        int numVariables = vlistNvars(vlistID);
        int numRecordsConst = 0;
        int numRecordsVar = 0;
        int numValues = 0;
        size_t memorySize = 0;
        size_t inputSize = 0;
        size_t outputSize = 0;
        if (Options::test && xsInfo)
          {
            VarList varList;
            varListInit(varList, vlistID);
            for (int varID = 0; varID < numVariables; ++varID)
              {
                const auto &var = varList[varID];

                if (var.isConstant)
                  numRecordsConst += var.nlevels;
                else
                  numRecordsVar += var.nlevels;

                auto size = var.nlevels * var.gridsize * var.nwpv;
                numValues += size;

                auto numBytes = (var.memType == MemType::Double) ? 8 : 4;
                memorySize += size * numBytes;

                auto numBits = get_num_input_bits(var.datatype);
                inputSize += (size * numBits) / 8;

                numBits = get_num_output_bits(var.datatype);
                outputSize += (size * numBits) / 8;
              }
          }

        print_vars_info(operfunc, ensembleInfo, vlistID, xsInfo);

        set_text_color(stdout, BRIGHT);
        fprintf(stdout, "   Grid coordinates");
        reset_text_color(stdout);
        fprintf(stdout, " :\n");

        print_grid_info(vlistID);

        set_text_color(stdout, BRIGHT);
        fprintf(stdout, "   Vertical coordinates");
        reset_text_color(stdout);
        fprintf(stdout, " :\n");

        print_zaxis_info(vlistID);

        if (nsubtypes > 1)
          {
            fprintf(stdout, "   Subtypes");
            fprintf(stdout, " :\n");

            print_subtype_info(vlistID);
          }

        auto taxisID = vlistInqTaxis(vlistID);
        auto ntsteps = vlistNtsteps(vlistID);

        int numTimesteps = ntsteps;

        if (ntsteps != 0)
          {
            if (xsInfo) { numTimesteps = print_time_info_xs(ntsteps, taxisID, streamID); }
            else
              {
                print_time_info(ntsteps, taxisID);

                fprintf(stdout, "  YYYY-MM-DD hh:mm:ss  YYYY-MM-DD hh:mm:ss  YYYY-MM-DD hh:mm:ss  YYYY-MM-DD hh:mm:ss\n");

                set_text_color(stdout, MAGENTA);
                print_timesteps(streamID, taxisID, Options::cdoVerbose);
                reset_text_color(stdout);
                fprintf(stdout, "\n");
              }
          }

        if (Options::test && xsInfo)
          {
            set_text_color(stdout, BRIGHT);
            fprintf(stdout, "   Summary");
            reset_text_color(stdout);
            fprintf(stdout, " :\n");

            int numFields = numRecordsVar * numTimesteps + numRecordsConst;
            fprintf(stdout, "%33s : %d\n", "number of fields", numFields);
            fprintf(stdout, "%33s : %d\n", "number of variables", numVariables);
            if (numTimesteps) fprintf(stdout, "%33s : %d\n", "number of timesteps", numTimesteps);
            fprintf(stdout, "%33s : %d\n", "number of values", numTimesteps * numValues);
            // fprintf(stdout, "%33s : %s\n", "required memory", num_values_to_byte_cstr(memorySize));
            fprintf(stdout, "%33s : ~%s\n", "input size", num_values_to_byte_cstr(numTimesteps * inputSize));
            fprintf(stdout, "%33s : ~%s\n", "output size", num_values_to_byte_cstr(numTimesteps * outputSize));
          }

        cdo_stream_close(streamID);
      }
  }

  void
  close()
  {
    cdo_finish();
  }
};

void *
Sinfo(void *process)
{
  ModuleSinfo sinfo;
  sinfo.init(process);
  sinfo.run();
  sinfo.close();
  return nullptr;
}
