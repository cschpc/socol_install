/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Setpartab  setpartab       Set parameter table
*/

#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "table.h"
#include "param_conversion.h"
#include "cdo_cmor.h"
#include "pmlist.h"
#include "mapping.h"
#include "convert_units.h"
#include "util_files.h"
#include "cdi_lockedIO.h"

enum pt_mode_t
{
  CODE_NUMBER,
  PARAMETER_ID,
  VARIABLE_NAME,
  STANDARD_NAME
};

static void
apply_parameterList(pt_mode_t ptmode, PMList &pmlist, int nvars, int vlistID2, std::vector<CmorVar> &vars)
{
  const std::vector<std::string> hentry = { "Header" };
  const std::vector<std::string> ventry = { "variable_entry", "parameter" };
  char valstr[CDI_MAX_NAME];
  char paramstr[32];
  int codenum = 0;

  // search for global missing value
  auto hasMissvals = false;
  double missval = 0.0;

  {
    auto kvlist = pmlist.getKVListVentry(hentry);
    if (kvlist)
      {
        auto kv = kvlist->search("missing_value");
        if (kv && kv->nvalues > 0)
          {
            hasMissvals = true;
            missval = parameter_to_double(kv->values[0]);
          }
      }
  }

  int numVarsFound = 0;
  for (int varID = 0; varID < nvars; ++varID)
    {
      auto &var = vars[varID];
      var.name = cdo::inq_var_name(vlistID2, varID);

      if (hasMissvals)
        {
          auto missval_old = vlistInqVarMissval(vlistID2, varID);
          if (!dbl_is_equal(missval, missval_old))
            {
              var.changemissval = true;
              var.missval_old = missval_old;
              vlistDefVarMissval(vlistID2, varID, missval);
            }
        }

      const KVList *kvlist = nullptr;
      if (ptmode == CODE_NUMBER)
        {
          codenum = vlistInqVarCode(vlistID2, varID);
          std::snprintf(valstr, sizeof(valstr), "%d", codenum);
          kvlist = pmlist.searchKVListVentry("code", valstr, ventry);
          if (kvlist)
            {
              auto tableID = vlistInqVarTable(vlistID2, varID);
              auto tabnum = tableInqNum(tableID);
              int levtype = 0;
              cdiInqKeyInt(vlistInqVarZaxis(vlistID2, varID), CDI_GLOBAL, CDI_KEY_TYPEOFFIRSTFIXEDSURFACE, &levtype);
              auto table = tabnum;
              auto ltype = levtype;
              {
                auto kv = kvlist->search("table");
                if (kv && kv->nvalues == 1) table = parameter_to_int(kv->values[0]);
              }
              {
                auto kv = kvlist->search("ltype");
                if (kv && kv->nvalues == 1) ltype = parameter_to_int(kv->values[0]);
              }
              if (!(tabnum == table && levtype == ltype)) kvlist = nullptr;
            }
        }
      else if (ptmode == PARAMETER_ID)
        {
          auto param = vlistInqVarParam(vlistID2, varID);
          cdiParamToString(param, paramstr, sizeof(paramstr));
          std::snprintf(valstr, sizeof(valstr), "%s", paramstr);
          kvlist = pmlist.searchKVListVentry("param", valstr, ventry);
          if (kvlist)
            {
              int levtype = 0;
              cdiInqKeyInt(vlistInqVarZaxis(vlistID2, varID), CDI_GLOBAL, CDI_KEY_TYPEOFFIRSTFIXEDSURFACE, &levtype);
              auto kv = kvlist->search("ltype");
              auto ltype = (kv && kv->nvalues == 1) ? parameter_to_int(kv->values[0]) : levtype;
              if (levtype != ltype) kvlist = nullptr;
            }
        }
      else if (ptmode == VARIABLE_NAME) { kvlist = pmlist.searchKVListVentry("name", var.name, ventry); }

      if (kvlist)
        {
          numVarsFound++;
          int pnum, ptab, pdum;
          cdiDecodeParam(vlistInqVarParam(vlistID2, varID), &pnum, &ptab, &pdum);
          auto hasValidMin = false, hasValidMax = false;

          for (const auto &kv : *kvlist)
            {
              const auto &key = kv.key;
              mapvar(vlistID2, varID, kv, key, &var, hasValidMin, hasValidMax, ptab, (ptmode != VARIABLE_NAME));
            }
          if (hasValidMin && hasValidMax) var.checkvalid = true;
        }
      else if (Options::cdoVerbose)
        {
          // clang-format off
          if      (ptmode == CODE_NUMBER)   cdo_print("Code number %d not found in parameter table!", codenum);
          else if (ptmode == PARAMETER_ID)  cdo_print("Parameter ID %s not found in parameter table!", paramstr);
          else if (ptmode == VARIABLE_NAME) cdo_print("Variable %s not found in parameter table!", var.name);
          // clang-format on
        }
    }

  if (numVarsFound == 0)
    {
      // clang-format off
      if      (ptmode == CODE_NUMBER)   cdo_warning("None of the input variables has a code number that matches the entries in the parameter table!");
      else if (ptmode == PARAMETER_ID)  cdo_warning("None of the input variables has a parameter ID that matches the entries in the parameter table!");
      else if (ptmode == VARIABLE_NAME) cdo_warning("None of the input variables has a name that matches the entries in the parameter table!");
      // clang-format on
    }
}

void *
Setpartab(void *process)
{
  int tableID = -1;
  int tableformat = 0;
  auto deleteVars = false;

  cdo_initialize(process);

  auto SETCODETAB = cdo_operator_add("setcodetab", 0, 0, "parameter code table name");
  auto SETPARTABC = cdo_operator_add("setpartabc", 0, 0, "parameter table name");
  auto SETPARTABP = cdo_operator_add("setpartabp", 0, 0, "parameter table name");
  auto SETPARTABN = cdo_operator_add("setpartabn", 0, 0, "parameter table name");

  auto operatorID = cdo_operator_id();

  operator_input_arg(cdo_operator_enter(operatorID));

  if (cdo_operator_argc() < 1) cdo_abort("Too few arguments!");

  auto convertData = false;
  if (cdo_operator_argc() == 2)
    {
      if (cdo_operator_argv(1) == "convert")
        convertData = true;
      else
        cdo_abort("Unknown parameter: >%s<", cdo_operator_argv(1));
    }

  if (cdo_operator_argc() > 2) cdo_abort("Too many arguments!");

  pt_mode_t ptmode = CODE_NUMBER;
  // clang-format off
  if      (operatorID == SETCODETAB) ptmode = CODE_NUMBER;
  else if (operatorID == SETPARTABC) ptmode = CODE_NUMBER;
  else if (operatorID == SETPARTABP) ptmode = PARAMETER_ID;
  else if (operatorID == SETPARTABN) ptmode = VARIABLE_NAME;
  // clang-format on

  if (ptmode == CODE_NUMBER)
    {
      const auto &partab = cdo_operator_argv(0);
      FILE *fp = FileUtils::file_exists(partab.c_str()) ? std::fopen(partab.c_str(), "r") : nullptr;
      if (fp != nullptr)
        {
          fseek(fp, 0L, SEEK_END);
          auto fsize = (size_t) ftell(fp);
          std::vector<char> parbuf(fsize + 1);
          fseek(fp, 0L, SEEK_SET);
          fread(parbuf.data(), fsize, 1, fp);
          parbuf[fsize] = 0;
          fseek(fp, 0L, SEEK_SET);

          if (atoi(parbuf.data()) == 0) tableformat = 1;

          std::fclose(fp);
        }

      if (tableformat == 0) tableID = cdo::define_table(partab);
    }
  else if (ptmode == PARAMETER_ID) { tableformat = 1; }
  else if (ptmode == VARIABLE_NAME) { tableformat = 1; }

  if (Options::cdoVerbose) cdo_print("Table format version %d", tableformat);

  auto streamID1 = cdo_open_read(0);

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto vlistID2 = vlistDuplicate(vlistID1);
  // vlistPrint(vlistID2);

  auto nvars = vlistNvars(vlistID2);
  std::vector<CmorVar> vars(nvars);

  if (convertData)
    for (int varID = 0; varID < nvars; ++varID) vars[varID].convert = true;

  if (tableformat == 0)
    {
      // for (int varID = 0; varID < nvars; ++varID) vlistDefVarTable(vlistID2, varID, tableID);
      char name[CDI_MAX_NAME], longname[CDI_MAX_NAME], units[CDI_MAX_NAME];
      for (int varID = 0; varID < nvars; ++varID)
        {
          auto param = vlistInqVarParam(vlistID2, varID);
          int pdis, pcat, pnum;
          cdiDecodeParam(param, &pnum, &pcat, &pdis);
          if (pdis == 255)
            {
              auto code = pnum;
              int ltype = 0;
              cdiInqKeyInt(vlistInqVarZaxis(vlistID2, varID), CDI_GLOBAL, CDI_KEY_TYPEOFFIRSTFIXEDSURFACE, &ltype);
              name[0] = 0;
              longname[0] = 0;
              units[0] = 0;
              tableInqEntry(tableID, code, ltype, name, longname, units);
              if (name[0])
                {
                  cdiDefKeyString(vlistID2, varID, CDI_KEY_NAME, name);
                  if (longname[0]) cdiDefKeyString(vlistID2, varID, CDI_KEY_LONGNAME, longname);
                  if (units[0]) cdiDefKeyString(vlistID2, varID, CDI_KEY_UNITS, units);
                }
            }
          vlistDefVarTable(vlistID2, varID, tableID);
        }
    }
  else
    {
      {
        auto filename = cdo_operator_argv(0).c_str();
        auto fp = std::fopen(filename, "r");
        if (fp == nullptr) cdo_abort("Open failed on: %s\n", filename);

        PMList pmlist;
        pmlist.read_namelist(fp, filename);
        std::fclose(fp);

        apply_parameterList(ptmode, pmlist, nvars, vlistID2, vars);
      }

      for (int varID = 0; varID < nvars; ++varID)
        if (vars[varID].remove)
          {
            deleteVars = true;
            break;
          }

      if (deleteVars)
        {
          vlistClearFlag(vlistID1);
          vlistClearFlag(vlistID2);

          for (int varID = 0; varID < nvars; ++varID)
            {
              auto zaxisID = vlistInqVarZaxis(vlistID2, varID);
              auto nlevs = zaxisInqSize(zaxisID);
              for (int levID = 0; levID < nlevs; levID++)
                {
                  vlistDefFlag(vlistID1, varID, levID, true);
                  vlistDefFlag(vlistID2, varID, levID, true);
                  if (vars[varID].remove)
                    {
                      vlistDefFlag(vlistID1, varID, levID, false);
                      vlistDefFlag(vlistID2, varID, levID, false);
                    }
                }
            }

          auto vlistIDx = vlistCreate();
          cdo_vlist_copy_flag(vlistIDx, vlistID2);

          vlistDestroy(vlistID2);

          vlistID2 = vlistIDx;
          if (vlistNvars(vlistID2) == 0) cdo_abort("No variable selected!");
        }

      for (int varID = 0; varID < nvars; ++varID)
        {
          auto &var = vars[varID];
          if (!var.convert) var.changeunits = false;
          if (var.changeunits)
            cdo_convert_units(&var.ut_converter, &var.changeunits, (char *) &var.units, (char *) &var.units_old, var.name.c_str());
        }
    }

  auto taxisID1 = vlistInqTaxis(vlistID1);
  auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  // vlistPrint(vlistID2);
  auto streamID2 = cdo_open_write(1);
  cdo_def_vlist(streamID2, vlistID2);

  auto gridsizemax = vlistGridsizeMax(vlistID1);
  if (vlistNumber(vlistID1) != CDI_REAL) gridsizemax *= 2;
  Varray<double> array(gridsizemax);

  VarList varList2;
  varListInit(varList2, vlistID2);

  int tsID = 0;
  while (true)
    {
      auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
      if (nrecs == 0) break;

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_def_timestep(streamID2, tsID);

      cmor_check_init(nvars, vars);

      for (int recID = 0; recID < nrecs; ++recID)
        {
          int varID, levelID;
          cdo_inq_record(streamID1, &varID, &levelID);

          auto &var = vars[varID];
          auto varID2 = varID;
          auto levelID2 = levelID;

          if (deleteVars)
            {
              if (var.remove) continue;

              if (vlistInqFlag(vlistID1, varID, levelID) == true)
                {
                  varID2 = vlistFindVar(vlistID2, varID);
                  levelID2 = vlistFindLevel(vlistID2, varID, levelID);
                }
            }

          cdo_def_record(streamID2, varID2, levelID2);

          size_t nmiss;
          cdo_read_record(streamID1, array.data(), &nmiss);

          auto missval = varList2[varID2].missval;
          auto gridsize = varList2[varID2].nwpv * varList2[varID2].gridsize;

          if (nmiss && var.changemissval)
            {
              for (size_t i = 0; i < gridsize; ++i)
                {
                  if (dbl_is_equal(array[i], var.missval_old)) array[i] = missval;
                }
            }

          if (var.lfactor)
            {
              for (size_t i = 0; i < gridsize; ++i)
                {
                  if (!dbl_is_equal(array[i], missval)) array[i] *= var.factor;
                }
            }

#ifdef HAVE_UDUNITS2
          if (var.changeunits)
            {
              int nerr = 0;
              for (size_t i = 0; i < gridsize; ++i)
                {
                  if (!dbl_is_equal(array[i], missval))
                    {
                      array[i] = cv_convert_double((const cv_converter *) var.ut_converter, array[i]);
                      if (ut_get_status() != UT_SUCCESS) nerr++;
                    }
                }
              if (nerr)
                {
                  cdo_warning("Udunits: Error converting units from [%s] to [%s], parameter: %s", var.units_old, var.units,
                              var.name);
                  var.changeunits = false;
                }
            }
#endif

          cdo_write_record(streamID2, array.data(), nmiss);

          cmor_check_prep(var, gridsize, missval, array.data());
        }

      cmor_check_eval(vlistID2, nvars, vars);

      tsID++;
    }

  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

#ifdef HAVE_UDUNITS2
  for (int varID = 0; varID < nvars; ++varID)
    if (vars[varID].changeunits) cdo_convert_free(vars[varID].ut_converter);

  cdo_convert_destroy();
#endif

  cdo_finish();

  return nullptr;
}
