/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "param_conversion.h"
#include "parse_literals.h"
#include "pmlist.h"
#include "util_wildcards.h"

constexpr int Delimiter = '@';

static int
get_datatype(char *buffer)
{
  int dtype = -1;

  auto slen = strlen(buffer);
  if (slen >= 3 && buffer[slen - 2] == ':')
    {
      if (slen >= 4 && buffer[slen - 3] == '\\')
        {
          for (int i = 2; i >= 0; i--) buffer[slen - i - 1] = buffer[slen - i];
        }
      else
        {
          auto type = buffer[slen - 1];
          // clang-format off
          if      (type == 's') dtype = 999;
          else if (type == 'd') dtype = CDI_DATATYPE_FLT64;
          else if (type == 'i') dtype = CDI_DATATYPE_INT32;
          else cdo_abort("Attribute type '%c' not supported!", type);
          // clang-format on
          buffer[slen - 2] = 0;
        }
    }

  return dtype;
}

static std::pair<std::string, std::string>
split_var_attr(const std::string &key, int delimiter)
{
  std::string varName, attName;

  auto sz = key.find_first_of(delimiter);

  if (sz == std::string::npos) { attName = key; }
  else
    {
      varName = key.substr(0, sz);
      attName = key.substr(sz + 1);
    }

  return std::make_pair(varName, attName);
}

static std::vector<int>
find_variables(const std::string &varName, int vlistID, std::vector<std::string> &wnames, int &cdiID)
{
  std::vector<int> varIDs;

  constexpr int Undefined = -99;
  cdiID = Undefined;
  auto nvars = vlistNvars(vlistID);
  // auto ngrids = vlistNgrids(vlistID);
  auto nzaxis = vlistNzaxis(vlistID);
  if (varName.size())
    {
      for (int idx = 0; idx < nvars; idx++)
        {
          auto name = cdo::inq_var_name(vlistID, idx);
          if (wildcardmatch(varName, name) == 0)
            {
              cdiID = vlistID;
              varIDs.push_back(idx);
            }
        }

      if (cdiID == Undefined)
        {
          /*
          for ( int idx = 0; idx < ngrids; idx++ )
            {
              int gridID = vlistGrid(vlistID, idx);
              auto xname = cdo::inq_key_string(gridID, CDI_XAXIS, CDI_KEY_NAME);
              if (wildcardmatch(varName, xname) == 0)
                {
                  cdiID = gridID;
                  varIDs.push_back(CDI_GLOBAL);
                }
              auto yname = cdo::inq_key_string(gridID, CDI_YAXIS, CDI_KEY_NAME);
              if (wildcardmatch(varName, yname) == 0)
                {
                  cdiID = gridID;
                  varIDs.push_back(CDI_GLOBAL);
                }
              }
          */
          for (int idx = 0; idx < nzaxis; idx++)
            {
              auto zaxisID = vlistZaxis(vlistID, idx);
              auto key = cdo::inq_key_string(zaxisID, CDI_GLOBAL, CDI_KEY_NAME);
              if (wildcardmatch(varName, key) == 0)
                {
                  cdiID = zaxisID;
                  varIDs.push_back(CDI_GLOBAL);
                }
            }
        }

      if (cdiID == Undefined)
        {
          auto printWarning = true;
          for (size_t i = 0; i < wnames.size(); ++i)
            {
              if (wnames[i].empty())
                {
                  wnames[i] = varName;
                  break;
                }
              if (wnames[i] == varName)
                {
                  printWarning = false;
                  break;
                }
            }
          if (printWarning) cdo_warning("Variable >%s< not found!", varName);
        }
    }
  else
    {
      cdiID = vlistID;
      varIDs.push_back(CDI_GLOBAL);
    }

  return varIDs;
}

static std::vector<std::string>
find_attribute(int cdiID, int varID, const std::string &attrName, int &dtype)
{
  std::vector<std::string> attrValues;

  int natts = 0;
  cdiInqNatts(cdiID, varID, &natts);
  for (int ia = 0; ia < natts; ia++)
    {
      char attname[CDI_MAX_NAME];
      int atttype, attlen;
      cdiInqAtt(cdiID, varID, ia, attname, &atttype, &attlen);

      if (attrName == attname)
        {
          if (atttype == CDI_DATATYPE_TXT)
            {
              std::vector<char> atttxt(attlen + 1);
              cdiInqAttTxt(cdiID, varID, attname, attlen, atttxt.data());
              atttxt[attlen] = '\0';
              attrValues.push_back(atttxt.data());
            }
          else if (atttype == CDI_DATATYPE_INT32)
            {
              std::vector<int> attint(attlen);
              cdiInqAttInt(cdiID, varID, attname, attlen, attint.data());
              for (int i = 0; i < attlen; ++i) attrValues.push_back(std::to_string(attint[i]));
              if (dtype == -1) dtype = atttype;
            }
          else if (atttype == CDI_DATATYPE_FLT32 || atttype == CDI_DATATYPE_FLT64)
            {
              std::vector<double> attflt(attlen);
              cdiInqAttFlt(cdiID, varID, attname, attlen, attflt.data());
              for (int i = 0; i < attlen; ++i) attrValues.push_back(std::to_string(attflt[i]));
              if (dtype == -1) dtype = atttype;
            }
          else { cdo_warning("Unsupported type %d name %s", atttype, attname); }
        }
    }

  if (varID != CDI_GLOBAL && attrValues.empty())
    {
      auto stdname = cdo::inq_key_string(cdiID, varID, CDI_KEY_STDNAME);
      auto longname = cdo::inq_var_longname(cdiID, varID);
      auto units = cdo::inq_var_units(cdiID, varID);

      auto param = vlistInqVarParam(cdiID, varID);
      char paramstr[32];
      param_to_string(param, paramstr, sizeof(paramstr));

      auto code = vlistInqVarCode(cdiID, varID);
      auto table = vlistInqVarTable(cdiID, varID);
      auto missval = vlistInqVarMissval(cdiID, varID);

      double addoffset = 0.0, scalefactor = 1.0;
      cdiInqKeyFloat(cdiID, varID, CDI_KEY_ADDOFFSET, &addoffset);
      cdiInqKeyFloat(cdiID, varID, CDI_KEY_SCALEFACTOR, &scalefactor);

      // clang-format off
      if      (attrName == "long_name")     attrValues.push_back(longname);
      else if (attrName == "standard_name") attrValues.push_back(stdname);
      else if (attrName == "units")         attrValues.push_back(units);
      else if (attrName == "param")         attrValues.push_back(paramstr);
      else if (attrName == "code")          attrValues.push_back(std::to_string(code));
      else if (attrName == "table")         attrValues.push_back(std::to_string(table));
      else if (attrName == "missing_value") attrValues.push_back(std::to_string(missval));
      else if (attrName == "add_offset")    attrValues.push_back(std::to_string(addoffset));
      else if (attrName == "scale_factor")  attrValues.push_back(std::to_string(scalefactor));
      // clang-format off
    }

  return attrValues;
}

static std::vector<std::string>
get_attribute(int vlistID, const std::string &varAttr, int &dtype)
{
  auto [varName, attName] = split_var_attr(varAttr, Delimiter);
  if (attName.empty()) cdo_abort("Attribute name missing in >%s<!", varAttr);

  int cdiID = CDI_UNDEFID;
  int varID = CDI_UNDEFID;
  if (varName.size())
    {
      auto nvars = vlistNvars(vlistID);
      for (int idx = 0; idx < nvars; idx++)
        {
          auto name = cdo::inq_var_name(vlistID, idx);
          if (varName == name)
            {
              cdiID = vlistID;
              varID = idx;
              break;
            }
        }
      if (varID == CDI_UNDEFID) cdo_abort("Variable %s not found!", varName);
    }
  else
    {
      cdiID = vlistID;
      varID = CDI_GLOBAL;
    }

  auto attrValues = find_attribute(cdiID, varID, attName, dtype);
  if (attrValues.empty()) cdo_abort("Attribute %s not found!", varAttr);

  return attrValues;
}

static void
delete_attribute(int cdiID, int varID, const std::string &attName)
{
  auto status = cdiDelAtt(cdiID, varID, attName.c_str());
  if (status != CDI_NOERR)  // check CDI keys
    {
      // clang-format off
      if      (attName == "long_name") cdiDeleteKey(cdiID, varID, CDI_KEY_LONGNAME);
      else if (attName == "units")     cdiDeleteKey(cdiID, varID, CDI_KEY_UNITS);
      // clang-format off
    }
}

static void
set_attribute(int cdiID, int varID, const std::string &attName, int dtype, const std::vector<std::string> &values)
{
  int nvalues = values.size();

  if (dtype == CDI_DATATYPE_INT8 || dtype == CDI_DATATYPE_INT16 || dtype == CDI_DATATYPE_INT32)
    {
      std::vector<int> ivals(nvalues);
      for (int i = 0; i < nvalues; ++i) ivals[i] = literal_to_int(values[i]);
      cdiDefAttInt(cdiID, varID, attName.c_str(), dtype, nvalues, ivals.data());
    }
  else if (dtype == CDI_DATATYPE_FLT32 || dtype == CDI_DATATYPE_FLT64)
    {
      Varray<double> dvals(nvalues);
      for (int i = 0; i < nvalues; ++i) dvals[i] = literal_to_double(values[i]);
      cdiDefAttFlt(cdiID, varID, attName.c_str(), dtype, nvalues, dvals.data());
    }
  else
    {
      if (nvalues > 1) cdo_abort("Multidimensional string attributes not supported! %s=\"%s\"", attName, values[1]);
      const auto &value = values[0];
      auto len = (int) value.size();
      int outlen = 0;
      std::vector<char> outvalue(len);
      for (int i = 0; i < len; ++i)
        {
          if (i > 0 && value[i - 1] == '\\' && value[i] == 'n')
            outvalue[outlen - 1] = '\n';
          else
            outvalue[outlen++] = value[i];
        }
      cdiDefAttTxt(cdiID, varID, attName.c_str(), outlen, outvalue.data());
    }
}

static void
set_attributes(const KVList &kvlist, int vlistID)
{
  int kvn = kvlist.size();
  std::vector<std::string> wnames(kvn);

  char buffer[CDI_MAX_NAME];
  for (const auto &kv : kvlist)
    {
      strcpy(buffer, kv.key.c_str());
      auto dtype = get_datatype(buffer);

      auto [varName, attName] = split_var_attr(buffer, Delimiter);
      if (attName.empty()) cdo_abort("Attribute name missing in >%s<!", kv.key);

      int cdiID;
      auto varIDs = find_variables(varName, vlistID, wnames, cdiID);
      int numVars = varIDs.size();

      if (cdiID >= -1 && numVars > 0)
        {
          if (kv.nvalues == 0 || (kv.nvalues == 1 && kv.values[0].empty()))
            {
              for (int idx = 0; idx < numVars; ++idx)
                delete_attribute(cdiID, varIDs[idx], attName);
            }
          else
            {
              std::vector<std::string> attrValues(1);
              auto useAttrValues = false;
              if (kv.nvalues == 1)
                {
                  const auto &value = kv.values[0];
                  if (value.size() > 2 && value[0] == '{' && value[value.size()-1] == '}')
                    {
                      attrValues = get_attribute(vlistID, value.substr(1,value.size()-2), dtype);
                      useAttrValues = true;
                    }
                }

              const auto &values = useAttrValues ? attrValues : kv.values;

              if (dtype == -1) dtype = literals_find_datatype(values.size(), values);

              for (int idx = 0; idx < numVars; ++idx)
                set_attribute(cdiID, varIDs[idx], attName, dtype, values);
            }
        }
    }
}

void *
Setattribute(void *process)
{
  cdo_initialize(process);

  cdo_operator_add("setattribute", 0, 0, "attributes");

  auto dataIsUnchanged = data_is_unchanged();

  auto operatorID = cdo_operator_id();

  operator_input_arg(cdo_operator_enter(operatorID));

  auto natts = cdo_operator_argc();
  if (natts == 0) cdo_abort("Parameter missing!");

  PMList pmlist;
  KVList kvlist;
  kvlist.name = cdo_module_name();
  if (kvlist.parse_arguments(natts, cdo_get_oper_argv()) != 0) cdo_abort("Parse error!");
  if (Options::cdoVerbose) kvlist.print();

  auto pkvlist = &kvlist;
  if (natts == 1)
    {
      KeyValues &kv = kvlist.front();
      if (kv.key == "FILE")
        {
          if (Options::cdoVerbose) cdo_print("Reading attributes from: %s", kv.values[0]);
          auto filename = parameter_to_word(kv.values[0].c_str());
          auto fp = std::fopen(filename, "r");
          if (fp == nullptr) cdo_abort("Open failed on: %s\n", filename);
          pmlist.read_namelist(fp, filename);
          pkvlist = &pmlist.front();
          std::fclose(fp);
          if (Options::cdoVerbose) pkvlist->print();
        }
    }

  auto streamID1 = cdo_open_read(0);

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto vlistID2 = vlistDuplicate(vlistID1);

  set_attributes(*pkvlist, vlistID2);

  auto taxisID1 = vlistInqTaxis(vlistID1);
  auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  auto streamID2 = cdo_open_write(1);
  cdo_def_vlist(streamID2, vlistID2);

  Field field;

  VarList varList1;
  varListInit(varList1, vlistID1);

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
          cdo_def_record(streamID2, varID, levelID);

          if (dataIsUnchanged)
            {
              cdo_copy_record(streamID2, streamID1);
            }
          else
            {
              field.init(varList1[varID]);
              cdo_read_record(streamID1, field);
              cdo_write_record(streamID2, field);
            }
        }

      tsID++;
    }

  cdo_stream_close(streamID1);
  cdo_stream_close(streamID2);

  cdo_finish();

  return nullptr;
}
