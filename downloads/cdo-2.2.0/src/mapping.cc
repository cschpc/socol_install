/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cstdlib>
#include <string>

#include "cdo_options.h"
#include "cdo_output.h"
#include "pmlist.h"
#include "convert_units.h"
#include "param_conversion.h"
#include "parse_literals.h"
#include "cdo_cmor.h"
#include "cdo_cdi_wrapper.h"
#include "compare.h"

int string_to_param(const char *paramstr);

void
mapvar(int vlistID, int varID, const KeyValues &kv, const std::string &key, CmorVar *var, bool &hasValidMin, bool &hasValidMax,
       int ptab, bool isnPtmodeName)
{
  const auto &value = kv.values[0];
  auto lv1 = (kv.nvalues == 1);

  // printf("key=%s  value=>%s<\n", key.c_str(), value.c_str());

  // clang-format off
  if (!var)
    {
      if (key == "cn")
        {
          auto name = cdo::inq_var_name(vlistID, varID);
          if (name[0] != 0) cdiDefAttTxt(vlistID, varID, "original_name", (int) name.size(), name.c_str());
          cdiDefKeyString(vlistID, varID, CDI_KEY_NAME, parameter_to_word(value.c_str()));
        }
      else if (key == "u")
        cdiDefKeyString(vlistID, varID, CDI_KEY_UNITS, value.c_str());
      else if (key == "cm")
        cdiDefAttTxt(vlistID, varID, "cell_methods", (int) value.size(), value.c_str());
      else if (key == "ca")
        cdiDefAttTxt(vlistID, varID, "character_axis", (int) value.size(), value.c_str());
      else if (key == "za")
        cdiDefAttTxt(vlistID, varID, "z_axis", (int) value.size(), value.c_str());
      else if (key == "vc")
        cdiDefAttTxt(vlistID, varID, "variable_comment", (int) value.size(), value.c_str());
      else if (key == "p")
        {
          if (!isspace(value[0])) cdiDefAttTxt(vlistID, varID, "positive", (int) value.size(), value.c_str());
        }
      else
        {
          if (Options::cdoVerbose) cdo_print("In applying the mapping table:\n          Key: '%s' is ignored.", key);
        }
    }
  else
    {
      if      (lv1 && key == "standard_name") cdiDefKeyString(vlistID, varID, CDI_KEY_STDNAME, value.c_str());
      else if (lv1 && key == "long_name") cdiDefKeyString(vlistID, varID, CDI_KEY_LONGNAME, value.c_str());
      else if (lv1 && key == "units") cdiDefKeyString(vlistID, varID, CDI_KEY_UNITS, value.c_str());
      else if (lv1 && key == "name")
        {
          if (isnPtmodeName) cdiDefKeyString(vlistID, varID, CDI_KEY_NAME, parameter_to_word(value.c_str()));
        }
      else if (lv1 && key == "out_name")
        {
          auto outname = parameter_to_word(value);
          if (var->name != outname)
            {
              cdiDefKeyString(vlistID, varID, CDI_KEY_NAME, outname.c_str());
              cdiDefAttTxt(vlistID, varID, "original_name", (int) var->name.size(), var->name.c_str());
            }
        }
      else if (lv1 && key == "param")
        vlistDefVarParam(vlistID, varID, string_to_param(parameter_to_word(value)));
      else if (lv1 && key == "out_param")
        vlistDefVarParam(vlistID, varID, string_to_param(parameter_to_word(value)));
      else if (lv1 && key == "code")
        vlistDefVarParam(vlistID, varID, cdiEncodeParam(parameter_to_int(value), ptab, 255));
      else if (lv1 && key == "out_code")
        vlistDefVarParam(vlistID, varID, cdiEncodeParam(parameter_to_int(value), ptab, 255));
      else if (lv1 && key == "uvRelativeToGrid")
        cdiDefKeyInt(vlistID, varID, CDI_KEY_UVRELATIVETOGRID, parameter_to_bool(value));
      else if (lv1 && key == "comment")
        cdiDefAttTxt(vlistID, varID, key.c_str(), (int) value.size(), value.c_str());
      else if (lv1 && key == "chunktype") ;
      else if (lv1 && key == "cell_methods")
        cdiDefAttTxt(vlistID, varID, key.c_str(), (int) value.size(), value.c_str());
      else if (lv1 && key == "cell_measures")
        cdiDefAttTxt(vlistID, varID, key.c_str(), (int) value.size(), value.c_str());
      else if (lv1 && key == "delete") var->remove = parameter_to_bool(value);
      else if (lv1 && key == "convert") var->convert = parameter_to_bool(value);
      else if (lv1 && key == "factor")
        {
          var->lfactor = true;
          var->factor = parameter_to_double(value);
          if (Options::cdoVerbose) cdo_print("%s - scale factor %g", var->name, var->factor);
        }
      else if (lv1 && (key == "missval" || key == "missing_value"))
        {
          auto missval = parameter_to_double(value);
          auto missval_old = vlistInqVarMissval(vlistID, varID);
          if (!DBL_IS_EQUAL(missval, missval_old))
            {
              if (Options::cdoVerbose) cdo_print("%s - change missval from %g to %g", var->name, missval_old, missval);
              var->changemissval = true;
              var->missval_old = missval_old;
              vlistDefVarMissval(vlistID, varID, missval);
            }
        }
      else if (lv1 && key == "valid_min")
        {
          hasValidMin = true;
          var->valid_min = parameter_to_double(value);
        }
      else if (lv1 && key == "valid_max")
        {
          hasValidMax = true;
          var->valid_max = parameter_to_double(value);
        }
      else if (lv1 && key == "ok_min_mean_abs")
        {
          var->check_min_mean_abs = true;
          var->ok_min_mean_abs = parameter_to_double(value);
        }
      else if (lv1 && key == "ok_max_mean_abs")
        {
          var->check_max_mean_abs = true;
          var->ok_max_mean_abs = parameter_to_double(value);
        }
      else if (lv1 && (key == "datatype" || key == "type"))
        {
          auto datatype = cdo::str_to_datatype(parameter_to_word(value));
          if (datatype != -1) vlistDefVarDatatype(vlistID, varID, datatype);
        }
      else if (lv1 && key == "dimensions")
        {
        }
      else
        {
          const auto &values = kv.values;
          const auto &rvalue = kv.values[0];
          int nvalues = kv.nvalues;
          if (nvalues == 1 && rvalue.empty()) nvalues = 0;

          int dtype = literals_find_datatype(nvalues, values);

          if (dtype == CDI_DATATYPE_INT8 || dtype == CDI_DATATYPE_INT16 || dtype == CDI_DATATYPE_INT32)
            {
              std::vector<int> ivals(nvalues);
              for (int i = 0; i < nvalues; ++i) ivals[i] = literal_to_int(values[i]);
              cdiDefAttInt(vlistID, varID, key.c_str(), dtype, nvalues, ivals.data());
            }
          else if (dtype == CDI_DATATYPE_FLT32 || dtype == CDI_DATATYPE_FLT64)
            {
              std::vector<double> dvals(nvalues);
              for (int i = 0; i < nvalues; ++i) dvals[i] = literal_to_double(values[i]);
              cdiDefAttFlt(vlistID, varID, key.c_str(), dtype, nvalues, dvals.data());
            }
          else
            {
              cdiDefAttTxt(vlistID, varID, key.c_str(), (int)rvalue.size(), rvalue.c_str());
            }
        }
    }
  // clang-format on
}
