/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Selvar     selparam        Select parameters by identifier (format: code.tabnum  or  pnum.cat.dis)
      Selvar     delparam        Delete parameters by identifier (format: code.tabnum  or  pnum.cat.dis)
      Selvar     selcode         Select parameters by code number
      Selvar     delcode         Delete parameters by code number
      Selvar     selname         Select parameters by name
      Selvar     delname         Delete parameters by name
      Selvar     selstdname      Select parameters by CF standard name
      Selvar     sellevel        Select levels
      Selvar     sellevidx       Select levels by index
      Selvar     selgrid         Select grids
      Selvar     selzaxis        Select zaxis
      Selvar     seltabnum       Select parameter table number
      Selvar     selltype        Select GRIB level type
*/

#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "cdo_vlist.h"
#include "cdo_zaxis.h"
#include "util_wildcards.h"
#include "cdi_lockedIO.h"
#include "param_conversion.h"

void *
Selvar(void *process)
{
  int nsel = 0;
  char paramstr[32];
  char gridname[CDI_MAX_NAME];
  char zaxistypename[CDI_MAX_NAME];
  std::vector<int> intarr;
  std::vector<double> fltarr;

  cdo_initialize(process);

  auto dataIsUnchanged = data_is_unchanged();

#define INVERTS_SELECTION(id) (cdo_operator_f2(id) & 1)
#define TAKES_STRINGS(id) (cdo_operator_f2(id) & 2)
#define TAKES_INTEGERS(id) (cdo_operator_f2(id) & 4)
#define TAKES_FLOATS(id) (cdo_operator_f2(id) & 8)

  // clang-format off
  auto SELPARAM     = cdo_operator_add("selparam",     0, 2,   "parameters");
  auto SELCODE      = cdo_operator_add("selcode",      0, 4,   "code numbers");
  auto SELNAME      = cdo_operator_add("selname",      0, 2,   "variable names");
  auto SELSTDNAME   = cdo_operator_add("selstdname",   0, 2,   "standard names");
  auto SELLEVEL     = cdo_operator_add("sellevel",     0, 8,   "levels");
  auto SELLEVIDX    = cdo_operator_add("sellevidx",    0, 4,   "index of levels");
  auto SELGRID      = cdo_operator_add("selgrid",      0, 4|2, "list of grid names or numbers");
  auto SELZAXIS     = cdo_operator_add("selzaxis",     0, 4|2, "list of zaxis types or numbers");
  auto SELZAXISNAME = cdo_operator_add("selzaxisname", 0, 2,   "list of zaxis names");
  auto SELTABNUM    = cdo_operator_add("seltabnum",    0, 4,   "table numbers");
  auto DELPARAM     = cdo_operator_add("delparam",     1, 2|1, "parameter");
  auto DELCODE      = cdo_operator_add("delcode",      1, 1,   "code numbers");
  auto DELNAME      = cdo_operator_add("delname",      1, 2|1, "variable names");
  auto SELLTYPE     = cdo_operator_add("selltype",     0, 4,   "GRIB level types");
  // clang-format on

  auto operatorID = cdo_operator_id();

  operator_input_arg(cdo_operator_enter(operatorID));

  auto ldelete = (cdo_operator_f1(operatorID) == 1);

  int args_are_numeric = cdo_operator_argc() > 0 && isdigit(cdo_operator_argv(0)[0]);

  auto argnames = cdo_get_oper_argv();
  if (TAKES_STRINGS(operatorID) && !(TAKES_INTEGERS(operatorID) && args_are_numeric))
    {
      nsel = cdo_operator_argc();

      if (Options::cdoVerbose)
        for (int i = 0; i < nsel; ++i) cdo_print("name %d = %s", i + 1, argnames[i]);
    }
  else if (TAKES_FLOATS(operatorID))
    {
      fltarr = cdo_argv_to_flt(argnames);
      nsel = fltarr.size();

      if (Options::cdoVerbose)
        for (int i = 0; i < nsel; ++i) cdo_print("flt %d = %g", i + 1, fltarr[i]);
    }
  else
    {
      intarr = cdo_argv_to_int(argnames);
      nsel = intarr.size();

      if (Options::cdoVerbose)
        for (int i = 0; i < nsel; ++i) cdo_print("int %d = %d", i + 1, intarr[i]);
    }

  std::vector<bool> selfound;
  if (nsel)
    {
      selfound.resize(nsel);
      for (int i = 0; i < nsel; ++i) selfound[i] = false;
    }

  auto streamID1 = cdo_open_read(0);

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto nvars = vlistNvars(vlistID1);
  std::vector<bool> vars(nvars);

  VarList varList1;
  varListInit(varList1, vlistID1);

  if (operatorID == SELGRID && !args_are_numeric && nsel == 1 && argnames[0].compare(0, 4, "var=") == 0)
    {
      int gridnum = 0;
      const char *gridvarname = &argnames[0][4];
      if (*gridvarname == 0) cdo_abort("Variable name missing!");

      for (int varID = 0; varID < nvars; ++varID)
        {
          const auto &var = varList1[varID];
          if (var.name == gridvarname)
            {
              gridnum = 1 + vlistGridIndex(vlistID1, var.gridID);
              args_are_numeric = true;
              intarr.push_back(gridnum);
              break;
            }
        }

      if (!gridnum) cdo_abort("Variable %s not found!", gridvarname);
    }

  vlistClearFlag(vlistID1);
  for (int varID = 0; varID < nvars; ++varID)
    {
      vars[varID] = ldelete;

      const auto &var = varList1[varID];
      auto stdname = cdo::inq_key_string(vlistID1, varID, CDI_KEY_STDNAME);
      auto tabnum = tableInqNum(vlistInqVarTable(vlistID1, varID));
      auto grididx = vlistGridIndex(vlistID1, var.gridID);
      auto zaxisidx = vlistZaxisIndex(vlistID1, var.zaxisID);
      auto nlevs = zaxisInqSize(var.zaxisID);
      gridName(gridInqType(var.gridID), gridname);
      auto zaxisname = cdo::inq_key_string(var.zaxisID, CDI_GLOBAL, CDI_KEY_NAME);
      zaxisName(zaxisInqType(var.zaxisID), zaxistypename);

      cdiParamToString(var.param, paramstr, sizeof(paramstr));

      for (int levID = 0; levID < nlevs; levID++)
        {
          auto level = cdo_zaxis_inq_level(var.zaxisID, levID);

          if (ldelete) vlistDefFlag(vlistID1, varID, levID, true);

          for (int isel = 0; isel < nsel; isel++)
            {
              auto found = false;
              if (operatorID == SELCODE)
                found = (intarr[isel] == var.code);
              else if (operatorID == SELPARAM)
                found = (wildcardmatch(argnames[isel], paramstr) == 0);
              else if (operatorID == SELNAME)
                found = (wildcardmatch(argnames[isel], var.name) == 0);
              else if (operatorID == SELSTDNAME)
                found = (wildcardmatch(argnames[isel], stdname) == 0);
              else if (operatorID == SELLEVEL)
                found = (std::fabs(fltarr[isel] - level) < 0.0001);
              else if (operatorID == SELLEVIDX)
                found = (intarr[isel] == (levID + 1));
              else if (operatorID == SELGRID && args_are_numeric)
                found = (intarr[isel] == (grididx + 1));
              else if (operatorID == SELGRID && !args_are_numeric)
                found = (memcmp(argnames[isel].c_str(), gridname, argnames[isel].size()) == 0);
              else if (operatorID == SELZAXIS && args_are_numeric)
                found = (intarr[isel] == (zaxisidx + 1));
              else if (operatorID == SELZAXIS && !args_are_numeric)
                found = (memcmp(argnames[isel].c_str(), zaxistypename, argnames[isel].size()) == 0);
              else if (operatorID == SELZAXISNAME)
                found = (wildcardmatch(argnames[isel], zaxisname) == 0);
              else if (operatorID == SELTABNUM)
                found = (intarr[isel] == tabnum);
              else if (operatorID == DELCODE)
                found = (intarr[isel] == var.code);
              else if (operatorID == DELNAME)
                found = (wildcardmatch(argnames[isel], var.name) == 0);
              else if (operatorID == DELPARAM)
                found = (argnames[isel] == paramstr);
              else if (operatorID == SELLTYPE)
                found = (intarr[isel] == zaxis_to_ltype(var.zaxisID));

              if (found)
                {
                  vlistDefFlag(vlistID1, varID, levID, !INVERTS_SELECTION(operatorID));
                  selfound[isel] = true;
                  vars[varID] = !ldelete;
                }
            }
        }
    }

  int npar = 0;
  for (int varID = 0; varID < nvars; ++varID)
    if (vars[varID]) npar++;

  for (int varID = 0; varID < nvars; ++varID)
    {
      if (vars[varID])
        {
          const auto &var = varList1[varID];
          if (zaxisInqType(var.zaxisID) == ZAXIS_HYBRID)
            {
              auto psvarid = vlist_get_psvarid(vlistID1, var.zaxisID);
              if (psvarid != -1 && !vars[psvarid])
                {
                  vars[psvarid] = true;
                  vlistDefFlag(vlistID1, psvarid, 0, !INVERTS_SELECTION(operatorID));
                }
            }
        }
    }

  for (int isel = 0; isel < nsel; isel++)
    {
      if (selfound[isel] == false)
        {
          if (operatorID == SELCODE || operatorID == DELCODE)
            cdo_warning("Code number %d not found!", intarr[isel]);
          else if (operatorID == SELPARAM || operatorID == DELPARAM)
            cdo_warning("Parameter %s not found!", argnames[isel]);
          else if (operatorID == SELNAME || operatorID == DELNAME)
            cdo_warning("Variable name %s not found!", argnames[isel]);
          else if (operatorID == SELSTDNAME)
            cdo_warning("Variable with standard name %s not found!", argnames[isel]);
          else if (operatorID == SELLEVEL)
            cdo_warning("Level %g not found!", fltarr[isel]);
          else if (operatorID == SELLEVIDX)
            cdo_warning("Level index %d not found!", intarr[isel]);
          else if (operatorID == SELGRID && args_are_numeric)
            cdo_warning("Grid %d not found!", intarr[isel]);
          else if (operatorID == SELGRID && !args_are_numeric)
            cdo_warning("Grid name %s not found!", argnames[isel]);
          else if (operatorID == SELZAXIS && args_are_numeric)
            cdo_warning("Zaxis %d not found!", intarr[isel]);
          else if (operatorID == SELZAXIS && !args_are_numeric)
            cdo_warning("Zaxis type %s not found!", argnames[isel]);
          else if (operatorID == SELZAXISNAME)
            cdo_warning("Zaxis name %s not found!", argnames[isel]);
          else if (operatorID == SELTABNUM)
            cdo_warning("Table number %d not found!", intarr[isel]);
          else if (operatorID == SELLTYPE)
            cdo_warning("GRIB level type %d not found!", intarr[isel]);
        }
    }

  if (npar == 0) cdo_abort("No variables selected!");

  auto vlistID2 = vlistCreate();
  cdo_vlist_copy_flag(vlistID2, vlistID1);

  if (Options::cdoVerbose) vlistPrint(vlistID2);

  nvars = vlistNvars(vlistID2);
  {
    int varID;
    for (varID = 0; varID < nvars; ++varID)
      if (vlistInqVarTimetype(vlistID2, varID) != TIME_CONSTANT) break;
    if (varID == nvars) vlistDefNtsteps(vlistID2, 0);
  }

  auto taxisID1 = vlistInqTaxis(vlistID1);
  auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  auto streamID2 = cdo_open_write(1);
  cdo_def_vlist(streamID2, vlistID2);

  Field field;

  int tsID = 0;
  while (true)
    {
      auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
      if (nrecs == 0) break;

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_def_timestep(streamID2, tsID);

      for (int recID = 0; recID < nrecs; ++recID)
        {
          int varID, levelID;
          cdo_inq_record(streamID1, &varID, &levelID);
          if (vlistInqFlag(vlistID1, varID, levelID) == true)
            {
              auto varID2 = vlistFindVar(vlistID2, varID);
              auto levelID2 = vlistFindLevel(vlistID2, varID, levelID);
              cdo_def_record(streamID2, varID2, levelID2);

              if (dataIsUnchanged) { cdo_copy_record(streamID2, streamID1); }
              else
                {
                  const auto &var = varList1[varID];
                  field.init(var);
                  cdo_read_record(streamID1, field);
                  cdo_write_record(streamID2, field);
                }
            }
        }

      tsID++;
    }

  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  vlistDestroy(vlistID2);

  cdo_finish();

  return nullptr;
}
