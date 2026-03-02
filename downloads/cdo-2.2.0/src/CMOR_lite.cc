/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      CMORlite      cmorlite        CMOR lite
*/

#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "param_conversion.h"
#include "cdo_cmor.h"
#include "pmlist.h"
#include "convert_units.h"
#include "cdi_lockedIO.h"

static void
cdo_define_var_units(CmorVar &var, int vlistID2, int varID, const std::string &units)
{
  auto unitsOld = cdo::inq_var_units(vlistID2, varID);
  if (units != unitsOld)
    {
      if (unitsOld.size() > 0 && units.size() > 0)
        {
          var.changeunits = true;
          strcpy(var.units_old, unitsOld.c_str());
          strcpy(var.units, units.c_str());
        }

      cdiDefKeyString(vlistID2, varID, CDI_KEY_UNITS, units.c_str());
      cdiDefAttTxt(vlistID2, varID, "original_units", (int) unitsOld.size(), unitsOld.c_str());
    }
}

void
cmor_check_init(int nvars, std::vector<CmorVar> &vars)
{
  for (int varID = 0; varID < nvars; ++varID)
    {
      auto &var = vars[varID];
      if (var.checkvalid || var.check_min_mean_abs || var.check_max_mean_abs)
        {
          var.amean = 0;
          var.nvals = 0;
          var.n_lower_min = 0;
          var.n_greater_max = 0;
        }
    }
}

void
cmor_check_eval(int vlistID, int nvars, const std::vector<CmorVar> &vars)
{
  for (int varID = 0; varID < nvars; ++varID)
    {
      const auto &var = vars[varID];
      if (var.checkvalid || var.check_min_mean_abs || var.check_max_mean_abs)
        {
          auto amean = var.amean;
          auto nvals = var.nvals;

          if (nvals > 0) amean /= nvals;

          auto n_lower_min = var.n_lower_min;
          auto n_greater_max = var.n_greater_max;

          auto varname = cdo::inq_var_name(vlistID, varID);

          if (n_lower_min > 0)
            cdo_warning("Invalid value(s) detected for variable '%s': %ld values were lower than minimum valid value (%.4g).",
                        varname, n_lower_min, var.valid_min);
          if (n_greater_max > 0)
            cdo_warning("Invalid value(s) detected for variable '%s': %ld values were greater than maximum valid value (%.4g).",
                        varname, n_greater_max, var.valid_max);

          if (var.check_min_mean_abs)
            {
              if (amean < .1 * var.ok_min_mean_abs)
                cdo_warning("Invalid Absolute Mean for variable '%s' (%.5g) is lower by more than an order of magnitude than "
                            "minimum allowed: %.4g",
                            varname, amean, var.ok_min_mean_abs);

              if (amean < var.ok_min_mean_abs)
                cdo_warning("Invalid Absolute Mean for variable '%s' (%.5g) is lower than minimum allowed: %.4g", varname, amean,
                            var.ok_min_mean_abs);
            }

          if (var.check_max_mean_abs)
            {
              if (amean > 10. * var.ok_max_mean_abs)
                cdo_warning("Invalid Absolute Mean for variable '%s' (%.5g) is greater by more than an order of magnitude than "
                            "maximum allowed: %.4g",
                            varname, amean, var.ok_max_mean_abs);

              if (amean > var.ok_max_mean_abs)
                cdo_warning("Invalid Absolute Mean for variable '%s' (%.5g) is greater than maximum allowed: %.4g", varname, amean,
                            var.ok_max_mean_abs);
            }
        }
    }
}

void
cmor_check_prep(CmorVar &var, const long gridsize, const double missval, const double *const array)
{
  if (var.checkvalid || var.check_min_mean_abs || var.check_max_mean_abs)
    {
      double amean = 0;
      long nvals = 0;

      for (long i = 0; i < gridsize; ++i)
        {
          auto aval = array[i];
          if (!dbl_is_equal(aval, missval))
            {
              amean += std::fabs(aval);
              nvals++;
            }
        }

      var.amean += amean;
      var.nvals += nvals;

      long n_lower_min = 0, n_greater_max = 0;

      for (long i = 0; i < gridsize; ++i)
        {
          auto aval = array[i];
          if (!dbl_is_equal(aval, missval))
            {
              if (aval < var.valid_min) n_lower_min++;
              if (aval > var.valid_max) n_greater_max++;
            }
        }

      var.n_lower_min += n_lower_min;
      var.n_greater_max += n_greater_max;
    }
}

static void
apply_cmor_list(PMList &pmlist, int nvars, int vlistID2, std::vector<CmorVar> &vars)
{
  const std::vector<std::string> hentry = { "Header" };
  const std::vector<std::string> ventry = { "variable_entry", "parameter" };

  // search for global missing value
  auto hasMissvals = false;
  double missval = 0.0;

  {
    auto kvlist = pmlist.getKVListVentry(hentry);
    if (kvlist)
      {
        for (const auto &kv : *kvlist)
          {
            const auto &key = kv.key;
            const auto &value = kv.values[0];
            if (kv.nvalues != 1 || value.empty()) continue;

            if (key == "missing_value")
              {
                hasMissvals = true;
                missval = parameter_to_double(value);
              }
            else if (key == "table_id" || key == "modeling_realm" || key == "realm" || key == "project_id" || key == "frequency")
              {
                cdiDefAttTxt(vlistID2, CDI_GLOBAL, key.c_str(), (int) value.size(), value.c_str());
              }
          }
      }
  }

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

      auto kvlist = pmlist.searchKVListVentry("name", var.name, ventry);
      if (kvlist)
        {
          auto hasValidMin = false, hasValidMax = false;

          for (const auto &kv : *kvlist)
            {
              const auto &key = kv.key;
              const auto &value = kv.values[0];
              if (kv.nvalues != 1 || value.empty()) continue;
              auto value_cstr = value.c_str();

              // printf("key=%s  value=>%s<\n", key.c_str(), value.c_str());

              // clang-format off
              if      (key == "standard_name") cdiDefKeyString(vlistID2, varID, CDI_KEY_STDNAME, value_cstr);
              else if (key == "long_name")     cdiDefKeyString(vlistID2, varID, CDI_KEY_LONGNAME, value_cstr);
              else if (key == "units")         cdo_define_var_units(var, vlistID2, varID, value_cstr);
              else if (key == "name") ;     // cdiDefKeyString(vlistID2, varID, CDI_KEY_NAME, parameter_to_word(value));
              else if (key == "out_name")
                {
                  auto outname = parameter_to_word(value);
                  if (var.name != outname)
                    {
                      cdiDefKeyString(vlistID2, varID, CDI_KEY_NAME, outname.c_str());
                      cdiDefAttTxt(vlistID2, varID, "original_name", var.name.size(), var.name.c_str());
                    }
                }
              else if (key == "param")         vlistDefVarParam(vlistID2, varID, string_to_param(parameter_to_word(value)));
              else if (key == "out_param")     vlistDefVarParam(vlistID2, varID, string_to_param(parameter_to_word(value)));
              else if (key == "comment")       cdiDefAttTxt(vlistID2, varID, key.c_str(), (int) value.size(), value_cstr);
              else if (key == "cell_methods")  cdiDefAttTxt(vlistID2, varID, key.c_str(), (int) value.size(), value_cstr);
              else if (key == "cell_measures") cdiDefAttTxt(vlistID2, varID, key.c_str(), (int) value.size(), value_cstr);
              else if (key == "delete")        var.remove = parameter_to_bool(value);
              else if (key == "convert")       var.convert = parameter_to_bool(value);
              else if (key == "factor")
                {
                  var.lfactor = true;
                  var.factor = parameter_to_double(value);
                  if (Options::cdoVerbose) cdo_print("%s - scale factor %g", var.name, var.factor);
                }
              else if (key == "missval" || key == "missing_value")
                {
                  missval = parameter_to_double(value);
                  auto missval_old = vlistInqVarMissval(vlistID2, varID);
                  if (!dbl_is_equal(missval, missval_old))
                    {
                      if (Options::cdoVerbose) cdo_print("%s - change missval from %g to %g", var.name, missval_old, missval);
                      var.changemissval = true;
                      var.missval_old = missval_old;
                      vlistDefVarMissval(vlistID2, varID, missval);
                    }
                }
              else if (key == "valid_min")
                {
                  hasValidMin = true;
                  var.valid_min = parameter_to_double(value);
                }
              else if (key == "valid_max")
                {
                  hasValidMax = true;
                  var.valid_max = parameter_to_double(value);
                }
              else if (key == "ok_min_mean_abs")
                {
                  var.check_min_mean_abs = true;
                  var.ok_min_mean_abs = parameter_to_double(value);
                }
              else if (key == "ok_max_mean_abs")
                {
                  var.check_max_mean_abs = true;
                  var.ok_max_mean_abs = parameter_to_double(value);
                }
              else if (key == "datatype" || key == "type")
                {
                  auto datatype = cdo::str_to_datatype(parameter_to_word(value));
                  if (datatype != -1) vlistDefVarDatatype(vlistID2, varID, datatype);
                }
              else
                {
                  if (Options::cdoVerbose) cdo_print("Attribute %s:%s not supported!", var.name, key);
                }
              // clang-format on
            }

          if (hasValidMin && hasValidMax) var.checkvalid = true;
        }
      else { cdo_print("Variable %s not found in CMOR table!", var.name); }
    }
}

void *
CMOR_lite(void *process)
{
  auto deleteVars = false;

  cdo_initialize(process);

  Options::CMOR_Mode = 1;
  if (Options::CMOR_Mode) cdiDefGlobal("CMOR_MODE", Options::CMOR_Mode);

  cdo_operator_add("cmorlite", 0, 0, "parameter table name");

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

  auto streamID1 = cdo_open_read(0);

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto vlistID2 = vlistDuplicate(vlistID1);

  auto nvars = vlistNvars(vlistID2);
  std::vector<CmorVar> vars(nvars);

  if (convertData)
    for (int varID = 0; varID < nvars; ++varID) vars[varID].convert = true;

  auto filename = cdo_operator_argv(0).c_str();
  auto fp = std::fopen(filename, "r");
  if (fp == nullptr) cdo_abort("Open failed on: %s\n", filename);

  PMList pmlist;
  pmlist.read_cmor_table(fp, filename);
  std::fclose(fp);

  apply_cmor_list(pmlist, nvars, vlistID2, vars);

  VarList varList2;
  varListInit(varList2, vlistID2);

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
          for (int levID = 0; levID < varList2[varID].nlevels; levID++)
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

  auto taxisID1 = vlistInqTaxis(vlistID1);
  auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  // vlistPrint(vlistID2);
  auto streamID2 = cdo_open_write(1);
  cdo_def_vlist(streamID2, vlistID2);

  auto gridsizemax = vlistGridsizeMax(vlistID1);
  if (vlistNumber(vlistID1) != CDI_REAL) gridsizemax *= 2;
  Varray<double> array(gridsizemax);

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
          auto gridsize = varList2[varID2].gridsize;
          if (varList2[varID2].nwpv != CDI_REAL) gridsize *= 2;

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

  return 0;
}
