/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Exprf      expr            Evaluate expressions
      Exprf      exprf           Evaluate expressions from script file
      Exprf      aexpr           Append evaluated expressions
      Exprf      aexprf          Append evaluated expressions from script file
*/
/*
  Operatoren: +, -, *, \, ^, ==, !=, >, <, >=, <=, <=>, &&, ||, ?:
  Functions: sqrt, exp, log, log10, sin, cos, tan, asin, acos, atan
  Functions: min, max, avg, std, var
  Constansts: M_PI, M_E
*/

#include <fstream>
#include <algorithm>
#include <cassert>

#include "cdi.h"
#include "julian_date.h"

#include "cdo_options.h"
#include "cdo_vlist.h"
#include "dmemory.h"
#include "process_int.h"
#include <mpim_grid.h>
#include "expr.h"
#include "cdo_zaxis.h"
#include "cdi_lockedIO.h"
#include "util_string.h"

void gridcell_areas(int gridID, Varray<double> &array);
int get_surface_ID(int vlistID);  // from Vertstat.cc
struct yy_buffer_state *yy_scan_string(const char *str, void *scanner);

static std::string
exprs_from_argument(const std::vector<std::string> &exprArgv)
{
  std::string exprString{};

  if (exprArgv.size() > 0)
    {
      for (size_t i = 0; i < exprArgv.size(); ++i)
        {
          if (i > 0) exprString += ",";
          exprString += exprArgv[i];
        }
      if (exprString[exprString.size() - 1] != ';') exprString += ";";
    }
  else { operator_check_argc(1); }

  return exprString;
}

static std::string
exprs_from_file(const std::vector<std::string> &exprArgv)
{
  if (exprArgv.size() != 1) operator_check_argc(1);
  auto exprFile = exprArgv[0];
  std::ifstream stream(exprFile);
  if (!stream.is_open()) cdo_abort("Open failed on %s", exprFile);
  std::stringstream buffer;
  buffer << stream.rdbuf();
  return buffer.str();
}

#define MAX_PARAMS 4096

static std::size_t
replace_all(std::string &inout, const std::string &what, const std::string &with)
{
  std::size_t count{};
  for (std::string::size_type pos{}; inout.npos != (pos = inout.find(what.data(), pos, what.length()));
       pos += with.length(), ++count)
    {
      inout.replace(pos, what.length(), with.data(), with.length());
    }
  return count;
}

static std::string
exprs_expand(std::string &exprString, const VarList &varList)
{
  auto replaceTemplate = false;
  std::string templateName = "_ALL_";

  if (exprString.find(templateName) != std::string::npos)
    {
      replaceTemplate = true;
      for (size_t varID = 0; varID < varList.size(); ++varID)
        {
          if (templateName == varList[varID].name)
            {
              replaceTemplate = false;
              break;
            }
        }
    }

  if (replaceTemplate)
    {
      std::string exprStringNew{};
      auto exprStringArgv = split_with_seperator(exprString, ';');
      for (auto &string : exprStringArgv)
        {
          if (string.find(templateName) == std::string::npos) { exprStringNew += string + ";"; }
          else
            {
              for (size_t varID = 0; varID < varList.size(); ++varID)
                {
                  auto tmpString = string;
                  replace_all(tmpString, templateName, varList[varID].name);
                  exprStringNew += tmpString + ";";
                }
            }
        }

      return exprStringNew;
    }

  return exprString;
}

static void
params_init(std::vector<ParamEntry> &params, const VarList &varList, int vlistID)
{
  auto nvars = vlistNvars(vlistID);
  for (int varID = 0; varID < nvars; ++varID)
    {
      const auto &var = varList[varID];

      auto stdname = cdo::inq_key_string(vlistID, varID, CDI_KEY_STDNAME);

      params[varID].type = ParamType::VAR;
      params[varID].isValid = true;
      params[varID].hasMV = true;
      params[varID].gridID = var.gridID;
      params[varID].zaxisID = var.zaxisID;
      params[varID].datatype = var.datatype;
      params[varID].steptype = var.timetype;
      params[varID].nlat = gridInqYsize(var.gridID);
      params[varID].ngp = var.gridsize;
      params[varID].nlev = var.nlevels;
      params[varID].missval = var.missval;
      params[varID].name = strdup(var.name.c_str());
      if (var.longname.size()) params[varID].longname = strdup(var.longname.c_str());
      if (var.units.size()) params[varID].units = strdup(var.units.c_str());
      if (stdname.size()) params[varID].stdname = strdup(stdname.c_str());
    }
}

static void
params_delete(const std::vector<ParamEntry> &params)
{
  for (int varID = 0; varID < MAX_PARAMS; ++varID)
    {
      if (params[varID].data) Free(params[varID].data);
      if (params[varID].name) Free(params[varID].name);
      if (params[varID].longname) Free(params[varID].longname);
      if (params[varID].stdname) Free(params[varID].stdname);
      if (params[varID].units) Free(params[varID].units);
    }
}

static void
params_add_coord(ParseParamType &parseArg, int coord, int cdiID, size_t size, const char *units, const char *longname)
{
  auto ncoords = parseArg.ncoords;
  if (ncoords >= parseArg.maxCoords) cdo_abort("Too many coordinates (limit=%d)", parseArg.maxCoords);

  parseArg.coords[ncoords].needed = false;
  parseArg.coords[ncoords].coord = coord;
  parseArg.coords[ncoords].cdiID = cdiID;
  parseArg.coords[ncoords].size = size;
  if (units)
    {
      parseArg.coords[ncoords].units.resize(strlen(units) + 1);
      strcpy(parseArg.coords[ncoords].units.data(), units);
    }
  if (longname)
    {
      parseArg.coords[ncoords].longname.resize(strlen(longname) + 1);
      strcpy(parseArg.coords[ncoords].longname.data(), longname);
    }

  parseArg.ncoords++;
}

int
params_get_coord_ID(const ParseParamType &parseArg, int coord, int cdiID)
{
  auto ncoords = parseArg.ncoords;
  for (int coordID = 0; coordID < ncoords; ++coordID)
    {
      if (parseArg.coords[coordID].coord == coord && parseArg.coords[coordID].cdiID == cdiID) return coordID;
    }

  cdo_abort("%s: coordinate %c not found!", __func__, coord);

  return -1;
}

static void
params_add_coordinates(int vlistID, ParseParamType &parseArg)
{
  auto ngrids = vlistNgrids(vlistID);
  for (int index = 0; index < ngrids; ++index)
    {
      auto gridID = vlistGrid(vlistID, index);
      auto size = gridInqSize(gridID);
      auto xunits = cdo::inq_key_string(gridID, CDI_XAXIS, CDI_KEY_UNITS);
      params_add_coord(parseArg, 'x', gridID, size, xunits.c_str(), "longitude");
      auto yunits = cdo::inq_key_string(gridID, CDI_YAXIS, CDI_KEY_UNITS);
      params_add_coord(parseArg, 'y', gridID, size, yunits.c_str(), "latitude");

      params_add_coord(parseArg, 'a', gridID, size, "m^2", "grid cell area");
      params_add_coord(parseArg, 'w', gridID, size, nullptr, "grid cell area weights");
    }

  auto nzaxis = vlistNzaxis(vlistID);
  for (int index = 0; index < nzaxis; ++index)
    {
      auto zaxisID = vlistZaxis(vlistID, index);
      auto size = zaxisInqSize(zaxisID);
      auto zunits = cdo::inq_key_string(zaxisID, CDI_GLOBAL, CDI_KEY_UNITS);
      auto longname = cdo::inq_key_string(zaxisID, CDI_GLOBAL, CDI_KEY_LONGNAME);
      params_add_coord(parseArg, 'z', zaxisID, size, zunits.c_str(), longname.c_str());
      params_add_coord(parseArg, 'i', zaxisID, size, zunits.c_str(), "level index");
      params_add_coord(parseArg, 'd', zaxisID, size, zunits.c_str(), "delta z");
    }
}

static int
params_add_ts(ParseParamType &parseArg)
{
  auto &params = parseArg.params;

  auto varID = parseArg.nparams;
  if (varID >= parseArg.maxparams) cdo_abort("Too many parameter (limit=%d)", parseArg.maxparams);

  params[varID].name = strdup("_timestep_info");
  params[varID].gridID = parseArg.pointID;
  params[varID].zaxisID = parseArg.surfaceID;
  params[varID].steptype = TIME_VARYING;
  params[varID].ngp = CoordIndex::LEN;
  params[varID].nlev = 1;

  parseArg.nparams++;
  parseArg.cnparams++;

  return varID;
}

static void
parse_param_init(ParseParamType &parseArg, int vlistID, int pointID, int zonalID, int surfaceID)
{
  auto nvars = vlistNvars(vlistID);
  auto ngrids = vlistNgrids(vlistID);
  auto nzaxis = vlistNzaxis(vlistID);
  auto maxCoords = ngrids * 4 + nzaxis * 3;

  parseArg.maxparams = MAX_PARAMS;
  parseArg.params.resize(MAX_PARAMS);
  parseArg.nparams = nvars;
  parseArg.cnparams = nvars;
  parseArg.nvars1 = nvars;
  parseArg.init = true;
  parseArg.debug = (Options::cdoVerbose != 0);
  parseArg.pointID = pointID;
  parseArg.zonalID = zonalID;
  parseArg.surfaceID = surfaceID;
  parseArg.needed.resize(nvars);
  parseArg.coords.resize(maxCoords);
  parseArg.maxCoords = maxCoords;
  parseArg.ncoords = 0;
}

static int
genZonalID(int vlistID)
{
  int zonalID = -1;

  auto ngrids = vlistNgrids(vlistID);
  for (int index = 0; index < ngrids; ++index)
    {
      auto gridID = vlistGrid(vlistID, index);
      auto gridtype = gridInqType(gridID);
      if (gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN || gridtype == GRID_GENERIC)
        if (gridInqXsize(gridID) > 1 && gridInqYsize(gridID) >= 1)
          {
            zonalID = gridToZonal(gridID);
            break;
          }
    }

  return zonalID;
}

static void
set_date_and_time(ParamEntry &varts, int calendar, int tsID, const CdiDateTime &vDateTime0, const CdiDateTime &vDateTime)
{
  double jdelta = 0.0;

  if (tsID)
    {
      auto julianDate0 = julianDate_encode(calendar, vDateTime0);
      auto julianDate = julianDate_encode(calendar, vDateTime);
      jdelta = julianDate_to_seconds(julianDate_sub(julianDate, julianDate0));
    }

  varts.data[CoordIndex::TIMESTEP] = tsID + 1;
  varts.data[CoordIndex::DATE] = cdiDate_get(vDateTime.date);
  varts.data[CoordIndex::TIME] = cdiTime_get(vDateTime.time);
  varts.data[CoordIndex::DELTAT] = jdelta;

  int year, mon, day;
  int hour, minute, second, ms;
  cdiDate_decode(vDateTime.date, &year, &mon, &day);
  cdiTime_decode(vDateTime.time, &hour, &minute, &second, &ms);

  varts.data[CoordIndex::DAY] = day;
  varts.data[CoordIndex::MONTH] = mon;
  varts.data[CoordIndex::YEAR] = year;
  varts.data[CoordIndex::SECOND] = second;
  varts.data[CoordIndex::MINUTE] = minute;
  varts.data[CoordIndex::HOUR] = hour;
}

static void
addOperators(void)
{
  // clang-format off
  cdo_operator_add("expr",   1, 1, "expressions");
  cdo_operator_add("exprf",  1, 0, "expr script filename");
  cdo_operator_add("aexpr",  0, 1, "expressions");
  cdo_operator_add("aexprf", 0, 0, "expr script filename");
  // clang-format on
}

class ModuleExpr
{
  CdoStreamID streamID1;
  CdoStreamID streamID2;

  int vlistID2;

  int taxisID1;
  int taxisID2;

  int vartsID;

  int nvars1;
  int nvars2;

  ParseParamType parseArg;
  CdiDateTime vDateTime0{};

  std::vector<int> varIDmap;
  std::string exprString;

  void *scanner = nullptr;

  int pointID;
  int zonalID;
  int surfaceID;
  int calendar;

public:
  void
  init(void *process)
  {
    cdo_initialize(process);

    yylex_init(&scanner);

    yyset_extra(&parseArg, scanner);

    addOperators();

    auto operatorID = cdo_operator_id();
    bool replacesVariables = cdo_operator_f1(operatorID);
    bool readsCommandLine = cdo_operator_f2(operatorID);

    operator_input_arg(cdo_operator_enter(operatorID));

    const auto &exprArgv = cdo_get_oper_argv();

    exprString = readsCommandLine ? exprs_from_argument(exprArgv) : exprs_from_file(exprArgv);

    streamID1 = cdo_open_read(0);
    auto vlistID1 = cdo_stream_inq_vlist(streamID1);

    VarList varList1;
    varListInit(varList1, vlistID1);

    exprString = exprs_expand(exprString, varList1);
    if (Options::cdoVerbose) cdo_print(exprString);

    nvars1 = vlistNvars(vlistID1);

    pointID = gridCreate(GRID_GENERIC, 1);
    zonalID = genZonalID(vlistID1);
    surfaceID = get_surface_ID(vlistID1);

    parse_param_init(parseArg, vlistID1, pointID, zonalID, surfaceID);

    auto &params = parseArg.params;
    params_init(parseArg.params, varList1, vlistID1);

    // Set all input variables to 'needed' if replacing is switched off
    for (int varID = 0; varID < nvars1; ++varID) parseArg.needed[varID] = !replacesVariables;

    // init function rand()
    std::srand(Options::Random_Seed);

    vartsID = params_add_ts(parseArg);
    parseArg.tsID = vartsID;
    params_add_coordinates(vlistID1, parseArg);

    CDO_parser_errorno = 0;
    yy_scan_string(exprString.c_str(), scanner);
    yyparse(parseArg, scanner);
    if (CDO_parser_errorno != 0) cdo_abort("Syntax error!");

    parseArg.init = false;

    if (Options::cdoVerbose)
      for (int varID = 0; varID < nvars1; ++varID)
        if (parseArg.needed[varID]) cdo_print("Needed var: %d %s", varID, params[varID].name);

    if (Options::cdoVerbose)
      for (int varID = 0; varID < parseArg.nparams; ++varID)
        cdo_print("var: %d %s ngp=%zu nlev=%zu coord=%c", varID, params[varID].name, params[varID].ngp, params[varID].nlev,
                  (params[varID].coord == 0) ? ' ' : params[varID].coord);

    varIDmap = std::vector<int>(parseArg.nparams);

    vlistID2 = vlistCreate();
    vlistDefNtsteps(vlistID2, vlistNtsteps(vlistID1));
    vlistClearFlag(vlistID1);
    if (!replacesVariables)
      {
        int pidx = 0;
        for (int varID = 0; varID < nvars1; ++varID)
          {
            params[varID].select = false;
            if (!params[varID].remove)
              {
                varIDmap[pidx++] = varID;
                auto nlevels = varList1[varID].nlevels;
                // printf("Replace %d nlevs %d\n", varID, nlevels);
                for (int levID = 0; levID < nlevels; levID++) vlistDefFlag(vlistID1, varID, levID, true);
              }
          }
      }
    cdo_vlist_copy_flag(vlistID2, vlistID1);  // Copy global attributes

    //  printf("parseArg.nparams %d\n", parseArg.nparams);
    for (int pidx = 0; pidx < parseArg.nparams; pidx++)
      {
        const auto &param = params[pidx];
        if (pidx < nvars1 && !param.select) continue;
        if (pidx >= nvars1)
          {
            if (param.type == ParamType::CONST) continue;
            if (param.name[0] == '_') continue;
            if (param.remove) continue;
            if (param.coord) continue;
          }

        // printf("gridID %d zaxisID %d\n",  param.gridID, param.zaxisID);
        auto varID = vlistDefVar(vlistID2, param.gridID, param.zaxisID, param.steptype);
        cdiDefKeyString(vlistID2, varID, CDI_KEY_NAME, param.name);
        // printf("add: %d %s %d levs %d\n", pidx,  param.name, varID, zaxisInqSize(param.zaxisID));
        if (param.hasMV) vlistDefVarMissval(vlistID2, varID, param.missval);
        if (param.units) cdiDefKeyString(vlistID2, varID, CDI_KEY_UNITS, param.units);
        if (param.longname) cdiDefKeyString(vlistID2, varID, CDI_KEY_LONGNAME, param.longname);
        if (param.stdname) cdiDefKeyString(vlistID2, varID, CDI_KEY_STDNAME, param.stdname);
        auto len = strlen(param.name);
        if (len > 3 && memcmp(param.name, "var", 3) == 0)
          {
            if (isdigit(param.name[3]))
              {
                auto code = atoi(param.name + 3);
                vlistDefVarCode(vlistID2, varID, code);
              }
          }
        varIDmap[varID] = pidx;
      }

    if (Options::cdoVerbose)
      {
        for (int varID = 0; varID < nvars1; ++varID)
          if (parseArg.needed[varID]) printf("needed: %d %s\n", varID, parseArg.params[varID].name);
        cdo_print("vlistNvars(vlistID1)=%d, vlistNvars(vlistID2)=%d", vlistNvars(vlistID1), vlistNvars(vlistID2));
      }

    nvars2 = vlistNvars(vlistID2);
    if (nvars2 == 0) cdo_abort("No output variable found!");

    for (int varID = 0; varID < nvars1; ++varID)
      {
        if (parseArg.needed[varID])
          {
            auto nItems = std::max((size_t) 4, params[varID].ngp * params[varID].nlev);
            params[varID].data = (double *) Malloc(nItems * sizeof(double));
          }
      }

    for (int varID = parseArg.nvars1; varID < parseArg.nparams; ++varID)
      {
        auto nItems = std::max((size_t) 4, params[varID].ngp * params[varID].nlev);
        params[varID].data = (double *) Malloc(nItems * sizeof(double));
      }

    for (int i = 0; i < parseArg.ncoords; ++i)
      {
        if (parseArg.coords[i].needed)
          {
            auto &cdata = parseArg.coords[i].data;
            auto csize = parseArg.coords[i].size;
            auto cdiID = parseArg.coords[i].cdiID;
            auto coord = parseArg.coords[i].coord;
            if (coord == 'x' || coord == 'y' || coord == 'a' || coord == 'w')
              {
                auto gridID = cdiID;
                auto ngp = csize;
                cdata.resize(ngp);
                if (coord == 'x' || coord == 'y')
                  {
                    gridID = generate_full_point_grid(gridID);
                    if (!gridHasCoordinates(gridID)) cdo_abort("Cell center coordinates missing!");

                    if (coord == 'x') gridInqXvals(gridID, cdata.data());
                    if (coord == 'y') gridInqYvals(gridID, cdata.data());

                    if (gridID != parseArg.coords[i].cdiID) gridDestroy(gridID);
                  }
                else if (coord == 'a') { gridcell_areas(gridID, cdata); }
                else if (coord == 'w')
                  {
                    cdata[0] = 1;
                    if (ngp > 1)
                      {
                        auto wstatus = gridcell_weights(gridID, cdata);
                        if (wstatus) cdo_warning("Grid cell bounds not available, using constant grid cell area weights!");
                      }
                  }
              }
            else if (coord == 'z' || coord == 'i' || coord == 'd')
              {
                auto zaxisID = cdiID;
                auto nlev = csize;
                cdata.resize(nlev);
                if (coord == 'z') { cdo_zaxis_inq_levels(zaxisID, cdata.data()); }
                else if (coord == 'i')
                  {
                    for (size_t k = 0; k < nlev; ++k) cdata[k] = k + 1;
                    cdo_zaxis_inq_levels(zaxisID, cdata.data());
                  }
                else if (coord == 'd')
                  {
                    varray_fill(nlev, cdata, 1.0);
                    if (zaxisInqLbounds(zaxisID, nullptr) && zaxisInqUbounds(zaxisID, nullptr))
                      {
                        std::vector<double> lbounds(nlev), ubounds(nlev);
                        zaxisInqLbounds(zaxisID, lbounds.data());
                        zaxisInqUbounds(zaxisID, ubounds.data());
                        for (size_t k = 0; k < nlev; ++k) cdata[k] = ubounds[k] - lbounds[k];
                      }
                  }
              }
            else
              cdo_abort("Computation of coordinate %c not implemented!", coord);
          }
      }

    for (int varID = parseArg.nvars1; varID < parseArg.nparams; ++varID)
      {
        auto coord = params[varID].coord;
        if (coord)
          {
            if (coord == 'x' || coord == 'y' || coord == 'a' || coord == 'w')
              {
                auto coordID = params_get_coord_ID(parseArg, coord, params[varID].gridID);
                auto gridID = parseArg.coords[coordID].cdiID;
                auto ngp = parseArg.coords[coordID].size;
                auto &cdata = parseArg.coords[coordID].data;
                assert(gridID == params[varID].gridID);
                assert(!cdata.empty());

                array_copy(ngp, cdata.data(), params[varID].data);
              }
            else if (coord == 'z' || coord == 'i' || coord == 'd')
              {
                auto coordID = params_get_coord_ID(parseArg, coord, params[varID].zaxisID);
                auto zaxisID = parseArg.coords[coordID].cdiID;
                auto nlev = parseArg.coords[coordID].size;
                auto &cdata = parseArg.coords[coordID].data;
                assert(zaxisID == params[varID].zaxisID);
                assert(!cdata.empty());

                array_copy(nlev, cdata.data(), params[varID].data);
              }
            else
              cdo_abort("Computation of coordinate %c not implemented!", coord);
          }
      }

    if (Options::cdoVerbose) vlistPrint(vlistID2);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);

    calendar = taxisInqCalendar(taxisID1);
  }

  void
  run()
  {
    int tsID = 0;
    while (true)
      {

        auto &params = parseArg.params;
        auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
        if (nrecs == 0) break;

        auto vDateTime = taxisInqVdatetime(taxisID1);

        set_date_and_time(params[vartsID], calendar, tsID, vDateTime0, vDateTime);

        vDateTime0 = vDateTime;

        cdo_taxis_copy_timestep(taxisID2, taxisID1);

        cdo_def_timestep(streamID2, tsID);

        // for (int varID = 0; varID < nvars1; ++varID) printf(">>> %s %d\n", params[varID].name, params[varID].isValid);
        for (int varID = 0; varID < nvars1; ++varID) params[varID].isValid = true;
        for (int varID = 0; varID < nvars1; ++varID)
          if (tsID == 0 || params[varID].steptype != TIME_CONSTANT) params[varID].nmiss = 0;

        for (int recID = 0; recID < nrecs; ++recID)
          {
            int varID, levelID;
            cdo_inq_record(streamID1, &varID, &levelID);
            if (parseArg.needed[varID])
              {
                auto offset = params[varID].ngp * levelID;
                auto vardata = params[varID].data + offset;
                size_t nmiss;
                cdo_read_record(streamID1, vardata, &nmiss);
                params[varID].nmiss += nmiss;

                if (nmiss > 0) cdo_check_missval(params[varID].missval, params[varID].name);
              }
          }

        for (int varID = 0; varID < nvars2; ++varID)
          {
            auto pidx = varIDmap[varID];
            if (pidx < nvars1) continue;

            params[pidx].nmiss = 0;
            varray_fill(params[pidx].ngp * params[pidx].nlev, params[pidx].data, 0.0);
          }

        parseArg.cnparams = vartsID + 1;
        yy_scan_string(exprString.c_str(), scanner);
        yyparse(parseArg, scanner);

        for (int varID = 0; varID < nvars2; ++varID)
          {
            auto pidx = varIDmap[varID];

            if (tsID > 0 && params[pidx].steptype == TIME_CONSTANT) continue;

            auto missval = vlistInqVarMissval(vlistID2, varID);

            auto ngp = params[pidx].ngp;
            auto nlev = (int) params[pidx].nlev;
            for (int levelID = 0; levelID < nlev; ++levelID)
              {
                auto offset = ngp * levelID;
                double *vardata = params[pidx].data + offset;
                auto nmiss = array_num_mv(ngp, vardata, missval);
                cdo_def_record(streamID2, varID, levelID);
                cdo_write_record(streamID2, vardata, nmiss);
              }
          }

        tsID++;
      }
  }

  void
  close()
  {
    cdo_stream_close(streamID2);
    cdo_stream_close(streamID1);

    vlistDestroy(vlistID2);

    yylex_destroy(scanner);

    params_delete(parseArg.params);

    gridDestroy(pointID);
    if (zonalID != -1) gridDestroy(zonalID);

    cdo_finish();
  }
};

void *
Expr(void *process)
{
  ModuleExpr expr;
  expr.init(process);
  expr.run();
  expr.close();
  return nullptr;
}
