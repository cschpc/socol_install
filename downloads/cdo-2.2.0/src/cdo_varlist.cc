/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include "cdo_varlist.h"
#include "cdo_cdi_wrapper.h"
#include "cdo_output.h"
#include "util_string.h"
#include "compare.h"
#include "stdnametable.h"

static bool
isIntType(int dataType)
{
  return (dataType == CDI_DATATYPE_UINT8 || dataType == CDI_DATATYPE_UINT16 || dataType == CDI_DATATYPE_INT16);
}

static bool
isFloatType(int dataType)
{
  return (dataType == CDI_DATATYPE_FLT32 || dataType == CDI_DATATYPE_CPX32);
}

void
varListInit(VarList &varList, int vlistID)
{
  auto numVars = vlistNvars(vlistID);
  varList.resize(numVars);

  for (int varID = 0; varID < numVars; ++varID)
    {
      auto &var = varList[varID];
      var.name = cdo::inq_var_name(vlistID, varID);
      var.longname = cdo::inq_var_longname(vlistID, varID);
      var.units = cdo::inq_var_units(vlistID, varID);
      var.gridID = vlistInqVarGrid(vlistID, varID);
      var.zaxisID = vlistInqVarZaxis(vlistID, varID);
      var.timetype = vlistInqVarTimetype(vlistID, varID);
      var.tsteptype = vlistInqVarTsteptype(vlistID, varID);
      var.gridsize = gridInqSize(var.gridID);
      var.nlevels = zaxisInqSize(var.zaxisID);
      var.datatype = vlistInqVarDatatype(vlistID, varID);
      var.missval = vlistInqVarMissval(vlistID, varID);
      var.code = vlistInqVarCode(vlistID, varID);
      var.param = vlistInqVarParam(vlistID, varID);
      var.nwpv = vlistInqVarNumber(vlistID, varID);
      var.isConstant = (var.timetype == TIME_CONSTANT);

      if (Options::CDO_Memtype == MemType::Native)
        {
          double addoffset = 0.0, scalefactor = 1.0;
          auto haveAddoffset = (cdiInqKeyFloat(vlistID, varID, CDI_KEY_ADDOFFSET, &addoffset) == CDI_NOERR);
          auto haveScalefactor = (cdiInqKeyFloat(vlistID, varID, CDI_KEY_SCALEFACTOR, &scalefactor) == CDI_NOERR);
          auto isPacked = (haveAddoffset || haveScalefactor);
          auto useFloatType = isFloatType(var.datatype) || (isIntType(var.datatype) && !isPacked);
          var.memType = useFloatType ? MemType::Float : MemType::Double;
        }
      else { var.memType = Options::CDO_Memtype; }
    }
}

void
varListSetMemtype(VarList &varList, MemType memType)
{
  for (auto &var : varList) var.memType = memType;
}

void
varListSetUniqueMemtype(VarList &varList)
{
  int numVars = varList.size();
  if (numVars)
    {
      auto memtype = varList[0].memType;
      int varID;
      for (varID = 1; varID < numVars; ++varID)
        {
          if (varList[varID].memType != memtype) break;
        }
      if (varID < numVars) varListSetMemtype(varList, MemType::Double);
    }
}

int
varList_numConstVars(const VarList &varList)
{
  int numConstVars = 0;
  int numVars = varList.size();
  for (int varID = 0; varID < numVars; ++varID)
    {
      const auto &var = varList[varID];
      if (var.timetype == TIME_CONSTANT) numConstVars++;
    }
  return numConstVars;
}

int
varList_numVaryingVars(const VarList &varList)
{
  int numVaryingVars = 0;
  int numVars = varList.size();
  for (int varID = 0; varID < numVars; ++varID)
    {
      const auto &var = varList[varID];
      if (var.timetype == TIME_VARYING) numVaryingVars++;
    }
  return numVaryingVars;
}

VarIDs
search_varIDs(const VarList &varList, int vlistID, int numFullLevels)
{
  VarIDs varIDs;

  auto numVars = vlistNvars(vlistID);

  auto useTable = false;
  for (int varID = 0; varID < numVars; ++varID)
    {
      auto tableNum = tableInqNum(vlistInqVarTable(vlistID, varID));
      if (tableNum > 0 && tableNum < 255)
        {
          useTable = true;
          break;
        }
    }

  if (Options::cdoVerbose && useTable) cdo_print("Using code tables!");

  char paramstr[32];
  gribcode_t gribcodes;

  for (int varID = 0; varID < numVars; ++varID)
    {
      auto &var = varList[varID];
      auto nlevels = var.nlevels;
      auto instNum = institutInqCenter(vlistInqVarInstitut(vlistID, varID));
      auto tableNum = tableInqNum(vlistInqVarTable(vlistID, varID));

      auto code = var.code;

      cdiParamToString(var.param, paramstr, sizeof(paramstr));
      int pnum, pcat, pdis;
      cdiDecodeParam(var.param, &pnum, &pcat, &pdis);
      if (pdis >= 0 && pdis < 255) code = -1;

      if (useTable)
        {
          if (tableNum == 2) { wmo_gribcodes(&gribcodes); }
          else if (tableNum == 128 || tableNum == 0 || tableNum == 255) { echam_gribcodes(&gribcodes); }
          //  KNMI: HIRLAM model version 7.2 uses tableNum=1    (LAMH_D11*)
          //  KNMI: HARMONIE model version 36 uses tableNum=1   (grib*) (opreational NWP version)
          //  KNMI: HARMONIE model version 38 uses tableNum=253 (grib,grib_md) and tableNum=1 (grib_sfx) (research version)
          else if (tableNum == 1 || tableNum == 253) { hirlam_harmonie_gribcodes(&gribcodes); }
        }
      else { echam_gribcodes(&gribcodes); }

      if (Options::cdoVerbose)
        cdo_print("Center=%d  TableNum=%d  Code=%d  Param=%s  Varname=%s  varID=%d", instNum, tableNum, code, paramstr, var.name,
                  varID);

      if (code <= 0 || code == 255)
        {
          auto varname = string_to_lower(cdo::inq_var_name(vlistID, varID));
          auto stdname = string_to_lower(cdo::inq_key_string(vlistID, varID, CDI_KEY_STDNAME));

          code = stdname_to_echamcode(stdname);
          if (code == -1)
            {
              //                                  ECHAM                 ECMWF
              // clang-format off
              if      (-1 == varIDs.sgeopotID && (varname == "geosp" || varname == "z")) code = gribcodes.geopot;
              else if (-1 == varIDs.tempID    && (varname == "st"    || varname == "t")) code = gribcodes.temp;
              else if (-1 == varIDs.psID      && (varname == "aps"   || varname == "sp")) code = gribcodes.ps;
              else if (-1 == varIDs.psID      &&  varname == "ps") code = gribcodes.ps;
              else if (-1 == varIDs.lnpsID    && (varname == "lsp"   || varname == "lnsp")) code = gribcodes.lsp;
              else if (-1 == varIDs.lnpsID2   &&  varname == "lnps") code = 777;
              else if (-1 == varIDs.geopotID  &&  stdname == "geopotential_full") code = gribcodes.geopot;
              else if (-1 == varIDs.tempID    &&  varname == "t") code = gribcodes.temp;
              else if (-1 == varIDs.humID     &&  varname == "q") code = gribcodes.hum;
              // else if (varname == "clwc") code = 246;
              // else if (varname == "ciwc") code = 247;
              // clang-format on
            }
        }

      // clang-format off
      if      (code == gribcodes.geopot  && nlevels == 1)             varIDs.sgeopotID = varID;
      else if (code == gribcodes.geopot  && nlevels == numFullLevels) varIDs.geopotID = varID;
      else if (code == gribcodes.temp    && nlevels == numFullLevels) varIDs.tempID = varID;
      else if (code == gribcodes.ps      && nlevels == 1)             varIDs.psID = varID;
      else if (code == gribcodes.lsp     && nlevels == 1)             varIDs.lnpsID = varID;
      else if (code == 777               && nlevels == 1)             varIDs.lnpsID2 = varID;
      else if (code == gribcodes.gheight && nlevels == numFullLevels) varIDs.gheightID = varID;
      else if (code == gribcodes.hum     && nlevels == numFullLevels) varIDs.humID = varID;
      // else if (code == 246 && nlevels == numFullLevels) varIDs.clwcID = varID;
      // else if (code == 247 && nlevels == numFullLevels) varIDs.ciwcID = varID;
      // clang-format on
    }

  return varIDs;
}
