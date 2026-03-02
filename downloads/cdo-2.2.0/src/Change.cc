/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Change     chcode          Change code number
      Change     chtabnum        Change GRIB1 parameter table number
      Change     chparam         Change parameter identifier
      Change     chname          Change variable or coordinate name
      Change     chlevel         Change level
      Change     chlevelc        Change level of one code
      Change     chlevelv        Change level of one variable
      Change     chltype         Change GRIB level type
*/

#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "param_conversion.h"

static void
changeCode(const VarList &varList1, int vlistID2, int nch, const std::vector<int> &chints)
{
  auto nvars = vlistNvars(vlistID2);
  for (int varID = 0; varID < nvars; ++varID)
    {
      auto code = varList1[varID].code;
      for (int i = 0; i < nch; i += 2)
        if (code == chints[i]) vlistDefVarCode(vlistID2, varID, chints[i + 1]);
    }
}

static void
changeTabnum(int vlistID2, int nch, const std::vector<int> &chints)
{
  auto nvars = vlistNvars(vlistID2);
  for (int varID = 0; varID < nvars; ++varID)
    {
      auto tabnum = tableInqNum(vlistInqVarTable(vlistID2, varID));
      for (int i = 0; i < nch; i += 2)
        if (tabnum == chints[i])
          {
            auto tableID = tableDef(-1, chints[i + 1], nullptr);
            vlistDefVarTable(vlistID2, varID, tableID);
          }
    }
}

static void
changeParam(const VarList &varList1, int vlistID2, int nch, const std::vector<const char *> &chnames)
{
  auto nvars = vlistNvars(vlistID2);
  for (int varID = 0; varID < nvars; ++varID)
    {
      auto param = varList1[varID].param;
      if (Options::cdoVerbose)
        {
          int pnum, pcat, pdis;
          cdiDecodeParam(param, &pnum, &pcat, &pdis);
          cdo_print("pnum, pcat, pdis: %d.%d.%d", pnum, pcat, pdis);
        }
      for (int i = 0; i < nch; i += 2)
        if (param == string_to_param(chnames[i])) vlistDefVarParam(vlistID2, varID, string_to_param(chnames[i + 1]));
    }
}

static void
changeName(const VarList &varList1, int vlistID2, int nch, const std::vector<const char *> &chnames)
{
  auto npairs = nch / 2;
  std::vector<std::pair<const char *, const char *>> vpairs(npairs);
  for (int i = 0; i < npairs; ++i) vpairs[i].first = chnames[i * 2];
  for (int i = 0; i < npairs; ++i) vpairs[i].second = chnames[i * 2 + 1];

  auto nvars = vlistNvars(vlistID2);
  std::vector<bool> namefound(npairs, false);
  for (int varID = 0; varID < nvars; ++varID)
    {
      for (int i = 0; i < npairs; ++i)
        if (varList1[varID].name == vpairs[i].first)
          {
            namefound[i] = true;
            cdiDefKeyString(vlistID2, varID, CDI_KEY_NAME, vpairs[i].second);
            break;
          }
    }

  auto searchForGridName = false;
  for (int i = 0; i < npairs; ++i)
    if (!namefound[i])
      {
        searchForGridName = true;
        break;
      }

  if (searchForGridName)
    {
      auto ngrids = vlistNgrids(vlistID2);
      for (int index = 0; index < ngrids; ++index)
        {
          int gridID2 = -1;
          auto gridID1 = vlistGrid(vlistID2, index);
          auto xname = cdo::inq_key_string(gridID1, CDI_XAXIS, CDI_KEY_NAME);
          auto yname = cdo::inq_key_string(gridID1, CDI_YAXIS, CDI_KEY_NAME);
          auto xfound = false, yfound = false;
          for (int i = 0; i < npairs; ++i)
            {
              if (!namefound[i])
                {
                  if (xname == vpairs[i].first)
                    {
                      xfound = true;
                      namefound[i] = true;
                      if (gridID2 == -1) gridID2 = gridDuplicate(gridID1);
                      cdiDefKeyString(gridID2, CDI_XAXIS, CDI_KEY_NAME, vpairs[i].second);
                    }
                }
              if (!namefound[i])
                {
                  if (yname == vpairs[i].first)
                    {
                      yfound = true;
                      namefound[i] = true;
                      if (gridID2 == -1) gridID2 = gridDuplicate(gridID1);
                      cdiDefKeyString(gridID2, CDI_YAXIS, CDI_KEY_NAME, vpairs[i].second);
                    }
                }

              if (xfound && yfound) break;
            }

          if (gridID2 != -1) vlistChangeGrid(vlistID2, gridID1, gridID2);
        }
    }

  auto searchForZaxisName = false;
  for (int i = 0; i < npairs; ++i)
    if (!namefound[i])
      {
        searchForZaxisName = true;
        break;
      }

  if (searchForZaxisName)
    {
      auto nzaxis = vlistNzaxis(vlistID2);
      for (int index = 0; index < nzaxis; ++index)
        {
          auto zaxisID1 = vlistZaxis(vlistID2, index);
          auto varname = cdo::inq_key_string(zaxisID1, CDI_GLOBAL, CDI_KEY_NAME);
          for (int i = 0; i < npairs; ++i)
            {
              if (!namefound[i])
                {
                  if (varname == vpairs[i].first)
                    {
                      namefound[i] = true;
                      auto zaxisID2 = zaxisDuplicate(zaxisID1);
                      cdiDefKeyString(zaxisID2, CDI_GLOBAL, CDI_KEY_NAME, vpairs[i].second);
                      vlistChangeZaxis(vlistID2, zaxisID1, zaxisID2);
                      break;
                    }
                }
            }
        }
    }

  for (int i = 0; i < npairs; ++i)
    if (!namefound[i]) cdo_warning("Variable name %s not found!", vpairs[i].first);
}

static void
changeUnit(const VarList &varList1, int vlistID2, int nch, const std::vector<const char *> &chnames)
{
  auto nvars = vlistNvars(vlistID2);
  for (int varID = 0; varID < nvars; ++varID)
    {
      for (int i = 0; i < nch; i += 2)
        if (varList1[varID].units == chnames[i]) cdiDefKeyString(vlistID2, varID, CDI_KEY_UNITS, chnames[i + 1]);
    }
}

static void
changeLevel(int vlistID2, int nch, const std::vector<double> &chlevels)
{
  auto nzaxis = vlistNzaxis(vlistID2);
  for (int index = 0; index < nzaxis; ++index)
    {
      auto zaxisID1 = vlistZaxis(vlistID2, index);
      if (zaxisInqLevels(zaxisID1, nullptr))
        {
          auto nlevs = zaxisInqSize(zaxisID1);
          Varray<double> levels(nlevs);
          zaxisInqLevels(zaxisID1, &levels[0]);

          int nfound = 0;
          for (int i = 0; i < nch; i += 2)
            for (int k = 0; k < nlevs; ++k)
              if (std::fabs(levels[k] - chlevels[i]) < 0.0001) nfound++;

          if (nfound)
            {
              Varray<double> newlevels = levels;
              auto zaxisID2 = zaxisDuplicate(zaxisID1);
              for (int i = 0; i < nch; i += 2)
                for (int k = 0; k < nlevs; ++k)
                  if (std::fabs(levels[k] - chlevels[i]) < 0.001) newlevels[k] = chlevels[i + 1];

              zaxisDefLevels(zaxisID2, &newlevels[0]);
              vlistChangeZaxis(vlistID2, zaxisID1, zaxisID2);
            }
        }
    }
}

static void
changeVarLevel(int varID, int vlistID2, const std::vector<double> &chlevels)
{
  auto zaxisID1 = vlistInqVarZaxis(vlistID2, varID);
  if (zaxisInqLevels(zaxisID1, nullptr))
    {
      auto nlevs = zaxisInqSize(zaxisID1);
      Varray<double> levels(nlevs);
      zaxisInqLevels(zaxisID1, &levels[0]);

      int nfound = 0;
      for (int k = 0; k < nlevs; ++k)
        if (std::fabs(levels[k] - chlevels[0]) < 0.0001) nfound++;

      if (nfound)
        {
          auto zaxisID2 = zaxisDuplicate(zaxisID1);
          for (int k = 0; k < nlevs; ++k)
            if (std::fabs(levels[k] - chlevels[0]) < 0.001) levels[k] = chlevels[1];

          zaxisDefLevels(zaxisID2, &levels[0]);
          vlistChangeVarZaxis(vlistID2, varID, zaxisID2);
        }
      else
        cdo_abort("Level %g not found!", chlevels[0]);
    }
}

static void
changeLevelByCode(int chcode, int vlistID2, const std::vector<double> &chlevels)
{
  int varID;
  auto nvars = vlistNvars(vlistID2);
  for (varID = 0; varID < nvars; ++varID)
    {
      auto code = vlistInqVarCode(vlistID2, varID);
      if (code == chcode) break;
    }
  if (varID == nvars) cdo_abort("Code %d not found!", chcode);

  changeVarLevel(varID, vlistID2, chlevels);
}

static void
changeLevelByName(const char *chname, const VarList &varList1, int vlistID2, const std::vector<double> &chlevels)
{
  int varID;
  auto nvars = vlistNvars(vlistID2);
  for (varID = 0; varID < nvars; ++varID)
    {
      if (varList1[varID].name == chname) break;
    }
  if (varID == nvars) cdo_abort("Variable name %s not found!", chname);

  changeVarLevel(varID, vlistID2, chlevels);
}

static void
changeLtype(int vlistID2, int nch, const std::vector<int> &chltypes)
{
  auto nzaxis = vlistNzaxis(vlistID2);
  for (int index = 0; index < nzaxis; ++index)
    {
      auto zaxisID1 = vlistZaxis(vlistID2, index);
      auto zaxisID2 = zaxisDuplicate(zaxisID1);
      int ltype = 0;
      cdiInqKeyInt(zaxisID1, CDI_GLOBAL, CDI_KEY_TYPEOFFIRSTFIXEDSURFACE, &ltype);

      for (int i = 0; i < nch; i += 2)
        {
          auto ltype1 = chltypes[i];
          auto ltype2 = chltypes[i + 1];
          if (ltype1 == ltype)
            {
              zaxisChangeType(zaxisID2, ZAXIS_GENERIC);
              cdiDefKeyInt(zaxisID2, CDI_GLOBAL, CDI_KEY_TYPEOFFIRSTFIXEDSURFACE, ltype2);
              vlistChangeZaxis(vlistID2, zaxisID1, zaxisID2);
            }
        }
    }
}

void *
Change(void *process)
{
  const char *chname = nullptr;
  int chcode = 0;
  std::vector<const char *> chnames;
  std::vector<int> chints, chltypes;
  std::vector<double> chlevels;

  cdo_initialize(process);

  // clang-format off
  auto CHCODE   = cdo_operator_add("chcode",   0, 0, "pairs of old and new code numbers");
  auto CHTABNUM = cdo_operator_add("chtabnum", 0, 0, "pairs of old and new GRIB1 table numbers");
  auto CHPARAM  = cdo_operator_add("chparam",  0, 0, "pairs of old and new parameter identifiers");
  auto CHNAME   = cdo_operator_add("chname",   0, 0, "pairs of old and new variable names");
  auto CHUNIT   = cdo_operator_add("chunit",   0, 0, "pairs of old and new variable units");
  auto CHLEVEL  = cdo_operator_add("chlevel",  0, 0, "pairs of old and new levels");
  auto CHLEVELC = cdo_operator_add("chlevelc", 0, 0, "code number, old and new level");
  auto CHLEVELV = cdo_operator_add("chlevelv", 0, 0, "variable name, old and new level");
  auto CHLTYPE  = cdo_operator_add("chltype",  0, 0, "pairs of old and new type");
  // clang-format on

  auto operatorID = cdo_operator_id();

  operator_input_arg(cdo_operator_enter(operatorID));

  auto nch = cdo_operator_argc();

  if (operatorID == CHCODE || operatorID == CHTABNUM)
    {
      if (nch % 2) cdo_abort("Odd number of input arguments!");
      chints.resize(nch);
      for (int i = 0; i < nch; ++i) chints[i] = parameter_to_int(cdo_operator_argv(i));
    }
  else if (operatorID == CHPARAM || operatorID == CHNAME || operatorID == CHUNIT)
    {
      if (nch % 2) cdo_abort("Odd number of input arguments!");
      chnames.resize(nch);
      for (int i = 0; i < nch; ++i) chnames[i] = &cdo_operator_argv(i)[0];
    }
  else if (operatorID == CHLEVEL)
    {
      if (nch % 2) cdo_abort("Odd number of input arguments!");
      chlevels.resize(nch);
      for (int i = 0; i < nch; ++i) chlevels[i] = parameter_to_double(cdo_operator_argv(i));
    }
  else if (operatorID == CHLEVELC)
    {
      operator_check_argc(3);

      chcode = parameter_to_int(cdo_operator_argv(0));
      chlevels.resize(2);
      chlevels[0] = parameter_to_double(cdo_operator_argv(1));
      chlevels[1] = parameter_to_double(cdo_operator_argv(2));
    }
  else if (operatorID == CHLEVELV)
    {
      operator_check_argc(3);

      chname = cdo_operator_argv(0).c_str();
      chlevels.resize(2);
      chlevels[0] = parameter_to_double(cdo_operator_argv(1));
      chlevels[1] = parameter_to_double(cdo_operator_argv(2));
    }
  else if (operatorID == CHLTYPE)
    {
      if (nch % 2) cdo_abort("Odd number of input arguments!");
      chltypes.resize(nch);
      for (int i = 0; i < nch; ++i) chltypes[i] = parameter_to_int(cdo_operator_argv(i));
    }

  auto streamID1 = cdo_open_read(0);

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto vlistID2 = vlistDuplicate(vlistID1);

  auto taxisID1 = vlistInqTaxis(vlistID1);
  auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  VarList varList1;
  varListInit(varList1, vlistID1);

  // clang-format off
  if      (operatorID == CHCODE)   changeCode(varList1, vlistID2, nch, chints);
  else if (operatorID == CHTABNUM) changeTabnum(vlistID2, nch, chints);
  else if (operatorID == CHPARAM)  changeParam(varList1, vlistID2, nch, chnames);
  else if (operatorID == CHNAME)   changeName(varList1, vlistID2, nch, chnames);
  else if (operatorID == CHUNIT)   changeUnit(varList1, vlistID2, nch, chnames);
  else if (operatorID == CHLEVEL)  changeLevel(vlistID2, nch, chlevels);
  else if (operatorID == CHLEVELC) changeLevelByCode(chcode, vlistID2, chlevels);
  else if (operatorID == CHLEVELV) changeLevelByName(chname, varList1, vlistID2, chlevels);
  else if (operatorID == CHLTYPE)  changeLtype(vlistID2, nch, chltypes);
  // clang-format on

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
          cdo_def_record(streamID2, varID, levelID);

          field.init(varList1[varID]);
          cdo_read_record(streamID1, field);
          cdo_write_record(streamID2, field);
        }

      tsID++;
    }

  cdo_stream_close(streamID1);
  cdo_stream_close(streamID2);

  cdo_finish();

  return nullptr;
}
