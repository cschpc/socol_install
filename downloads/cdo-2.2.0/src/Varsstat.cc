/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include "process_int.h"
#include "field_functions.h"

static void
checkUniqueZaxis(int vlistID)
{
  auto nzaxis = vlistNzaxis(vlistID);
  auto zaxisID = vlistZaxis(vlistID, 0);
  auto nlevels = zaxisInqSize(zaxisID);
  for (int index = 1; index < nzaxis; ++index)
    {
      if (nlevels != zaxisInqSize(vlistZaxis(vlistID, index))) cdo_abort("Number of level differ!");
    }
}

static void
checkUniqueGridsize(int vlistID)
{
  auto ngrids = vlistNgrids(vlistID);
  auto gridID = vlistGrid(vlistID, 0);
  auto gridsize = gridInqSize(gridID);
  for (int index = 0; index < ngrids; ++index)
    {
      if (gridsize != gridInqSize(vlistGrid(vlistID, index))) cdo_abort("Horizontal gridsize differ!");
    }
}

static void
setAttributes(const VarList &varList1, int vlistID2, int varID2, int operatorID)
{
  const auto &var0 = varList1[0];
  auto paramIsEqual = true;
  auto name = var0.name;
  auto param = var0.param;
  const int nvars = varList1.size();
  for (int varID = 1; varID < nvars; ++varID)
    {
      if (param != varList1[varID].param || name != varList1[varID].name)
        {
          paramIsEqual = false;
          break;
        }
    }

  if (!paramIsEqual) name = cdo_operator_name(operatorID);
  cdiDefKeyString(vlistID2, varID2, CDI_KEY_NAME, name.c_str());
  if (paramIsEqual)
    {
      if (param >= 0) vlistDefVarParam(vlistID2, varID2, param);
      if (var0.longname.size()) cdiDefKeyString(vlistID2, varID2, CDI_KEY_LONGNAME, var0.longname.c_str());
      if (var0.units.size()) cdiDefKeyString(vlistID2, varID2, CDI_KEY_UNITS, var0.units.c_str());
    }
}

static void
addOperators(void)
{
  // clang-format off
  cdo_operator_add("varsrange", FieldFunc_Range, 0, nullptr);
  cdo_operator_add("varsmin",   FieldFunc_Min,   0, nullptr);
  cdo_operator_add("varsmax",   FieldFunc_Max,   0, nullptr);
  cdo_operator_add("varssum",   FieldFunc_Sum,   0, nullptr);
  cdo_operator_add("varsmean",  FieldFunc_Mean,  0, nullptr);
  cdo_operator_add("varsavg",   FieldFunc_Avg,   0, nullptr);
  cdo_operator_add("varsvar",   FieldFunc_Var,   0, nullptr);
  cdo_operator_add("varsvar1",  FieldFunc_Var1,  0, nullptr);
  cdo_operator_add("varsstd",   FieldFunc_Std,   0, nullptr);
  cdo_operator_add("varsstd1",  FieldFunc_Std1,  0, nullptr);
  // clang-format on
}

void *
Varsstat(void *process)
{
  cdo_initialize(process);

  addOperators();

  auto operatorID = cdo_operator_id();
  auto operfunc = cdo_operator_f1(operatorID);

  operator_check_argc(0);

  auto lrange = (operfunc == FieldFunc_Range);
  auto lmean = (operfunc == FieldFunc_Mean || operfunc == FieldFunc_Avg);
  auto lstd = (operfunc == FieldFunc_Std || operfunc == FieldFunc_Std1);
  auto lvarstd = (lstd || operfunc == FieldFunc_Var || operfunc == FieldFunc_Var1);
  auto lvars2 = (lvarstd || lrange);
  const int divisor = (operfunc == FieldFunc_Std1 || operfunc == FieldFunc_Var1);

  auto field2_stdvar_func = lstd ? field2_std : field2_var;
  auto fieldc_stdvar_func = lstd ? fieldc_std : fieldc_var;

  auto streamID1 = cdo_open_read(0);
  auto vlistID1 = cdo_stream_inq_vlist(streamID1);

  VarList varList1;
  varListInit(varList1, vlistID1);

  checkUniqueZaxis(vlistID1);
  auto zaxisID = vlistZaxis(vlistID1, 0);
  auto nlevels = zaxisInqSize(zaxisID);

  checkUniqueGridsize(vlistID1);
  auto gridID = vlistGrid(vlistID1, 0);
  auto gridsize = gridInqSize(gridID);

  auto timetype = varList1[0].timetype;
  auto nvars = vlistNvars(vlistID1);
  for (int varID = 1; varID < nvars; ++varID)
    {
      if (timetype != varList1[varID].timetype) cdo_abort("Number of timesteps differ!");
    }

  auto vlistID2 = vlistCreate();
  vlistDefNtsteps(vlistID2, vlistNtsteps(vlistID1));

  auto varID2 = vlistDefVar(vlistID2, gridID, zaxisID, timetype);
  setAttributes(varList1, vlistID2, varID2, operatorID);

  auto taxisID1 = vlistInqTaxis(vlistID1);
  auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  auto streamID2 = cdo_open_write(1);
  cdo_def_vlist(streamID2, vlistID2);

  Field field;

  FieldVector vars1(nlevels), samp1(nlevels), vars2;
  if (lvars2) vars2.resize(nlevels);

  for (int levelID = 0; levelID < nlevels; ++levelID)
    {
      auto missval = varList1[0].missval;

      samp1[levelID].grid = gridID;
      samp1[levelID].missval = missval;
      samp1[levelID].memType = MemType::Double;
      vars1[levelID].grid = gridID;
      vars1[levelID].missval = missval;
      vars1[levelID].memType = MemType::Double;
      vars1[levelID].resize(gridsize);
      if (lvars2)
        {
          vars2[levelID].grid = gridID;
          vars2[levelID].missval = missval;
          vars2[levelID].memType = MemType::Double;
          vars2[levelID].resize(gridsize);
        }
    }

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

          auto &rsamp1 = samp1[levelID];
          auto &rvars1 = vars1[levelID];

          rvars1.nsamp++;
          if (lrange) vars2[levelID].nsamp++;

          if (varID == 0)
            {
              cdo_read_record(streamID1, rvars1);
              if (lrange)
                {
                  vars2[levelID].nmiss = rvars1.nmiss;
                  vars2[levelID].vec_d = rvars1.vec_d;
                }

              if (lvarstd) field2_moq(vars2[levelID], rvars1);

              if (rvars1.nmiss || !rsamp1.empty())
                {
                  if (rsamp1.empty()) rsamp1.resize(rvars1.size);
                  field2_vinit(rsamp1, rvars1);
                }
            }
          else
            {
              field.init(varList1[varID]);
              cdo_read_record(streamID1, field);

              if (field.nmiss || !rsamp1.empty())
                {
                  if (rsamp1.empty()) rsamp1.resize(rvars1.size, rvars1.nsamp);
                  field2_vincr(rsamp1, field);
                }

              // clang-format off
              if      (lvarstd) field2_sumsumq(rvars1, vars2[levelID], field);
              else if (lrange)  field2_maxmin(rvars1, vars2[levelID], field);
              else              field2_function(rvars1, field, operfunc);
              // clang-format on
            }
        }

      for (int levelID = 0; levelID < nlevels; ++levelID)
        {
          const auto &rsamp1 = samp1[levelID];
          auto &rvars1 = vars1[levelID];

          if (rvars1.nsamp)
            {
              if (lmean)
                {
                  if (!rsamp1.empty())
                    field2_div(rvars1, rsamp1);
                  else
                    fieldc_div(rvars1, (double) rvars1.nsamp);
                }
              else if (lvarstd)
                {
                  if (!rsamp1.empty())
                    field2_stdvar_func(rvars1, vars2[levelID], rsamp1, divisor);
                  else
                    fieldc_stdvar_func(rvars1, vars2[levelID], rsamp1.nsamp, divisor);
                }
              else if (lrange) { field2_sub(rvars1, vars2[levelID]); }

              cdo_def_record(streamID2, 0, levelID);
              cdo_write_record(streamID2, rvars1);
              rvars1.nsamp = 0;
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
