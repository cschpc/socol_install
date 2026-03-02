/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Vertstat   vertrange       Vertical range
      Vertstat   vertmin         Vertical minimum
      Vertstat   vertmax         Vertical maximum
      Vertstat   vertsum         Vertical sum
      Vertstat   vertint         Vertical integral
      Vertstat   vertmean        Vertical mean
      Vertstat   vertavg         Vertical average
      Vertstat   vertvar         Vertical variance
      Vertstat   vertvar1        Vertical variance [Normalize by (n-1)]
      Vertstat   vertstd         Vertical standard deviation
      Vertstat   vertstd1        Vertical standard deviation [Normalize by (n-1)]
*/

#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "cdo_vlist.h"
#include "cdo_zaxis.h"
#include "param_conversion.h"
#include "pmlist.h"
#include "cdi_lockedIO.h"
#include "field_functions.h"

#define IS_SURFACE_LEVEL(zaxisID) (zaxisInqType(zaxisID) == ZAXIS_SURFACE && zaxisInqSize(zaxisID) == 1)

int
get_surface_ID(int vlistID)
{
  int surfID = -1;

  auto nzaxis = vlistNzaxis(vlistID);
  for (int index = 0; index < nzaxis; ++index)
    {
      auto zaxisID = vlistZaxis(vlistID, index);
      if (IS_SURFACE_LEVEL(zaxisID))
        {
          surfID = vlistZaxis(vlistID, index);
          break;
        }
    }

  if (surfID == -1) surfID = zaxis_from_name("surface");

  return surfID;
}

static void
set_surface_ID(const int vlistID, const int surfID)
{
  auto nzaxis = vlistNzaxis(vlistID);
  for (int index = 0; index < nzaxis; ++index)
    {
      auto zaxisID = vlistZaxis(vlistID, index);
      if (zaxisID != surfID || !IS_SURFACE_LEVEL(zaxisID)) vlistChangeZaxisIndex(vlistID, index, surfID);
    }
}

static void
vertstat_get_parameter(bool &weights, bool &genbounds)
{
  auto pargc = cdo_operator_argc();
  if (pargc)
    {
      const auto &pargv = cdo_get_oper_argv();

      KVList kvlist;
      kvlist.name = cdo_module_name();
      if (kvlist.parse_arguments(pargc, pargv) != 0) cdo_abort("Parse error!");
      if (Options::cdoVerbose) kvlist.print();

      for (const auto &kv : kvlist)
        {
          const auto &key = kv.key;
          if (kv.nvalues > 1) cdo_abort("Too many values for parameter key >%s<!", key);
          if (kv.nvalues < 1) cdo_abort("Missing value for parameter key >%s<!", key);
          const auto &value = kv.values[0];

          // clang-format off
          if      (key == "weights")   weights = parameter_to_bool(value);
          else if (key == "genbounds") genbounds = parameter_to_bool(value);
          else cdo_abort("Invalid parameter key >%s<!", key);
          // clang-format on
        }
    }
}

class ModuleVertstat
{
  struct VertInfo
  {
    int zaxisID = -1;
    int status = -1;
    int numLevels = 0;
    Varray<double> thickness;
    Varray<double> weights;
  };

  int operatorID;
  int operfunc;

  CdoStreamID streamID1;
  CdoStreamID streamID2;

  int vlistID2;

  int taxisID1;
  int taxisID2;

  int nzaxis;
  int nvars;

  bool needWeights;
  bool lrange;
  bool lmean;
  bool lstd;
  bool lvarstd;
  int divisor;

  VarList varList;

  FieldVector vars1;
  FieldVector vars2;
  FieldVector samp1;

  Field field;

  std::vector<VertInfo> vert;
  int VERTINT;

public:
  void
  init(void *process)
  {
    cdo_initialize(process);

    // clang-format off
               cdo_operator_add("vertrange", FieldFunc_Range, 0, nullptr);
               cdo_operator_add("vertmin",   FieldFunc_Min,   0, nullptr);
               cdo_operator_add("vertmax",   FieldFunc_Max,   0, nullptr);
               cdo_operator_add("vertsum",   FieldFunc_Sum,   0, nullptr);
    VERTINT  = cdo_operator_add("vertint",   FieldFunc_Sum,   1, nullptr);
               cdo_operator_add("vertmean",  FieldFunc_Mean,  1, nullptr);
               cdo_operator_add("vertavg",   FieldFunc_Avg,   1, nullptr);
               cdo_operator_add("vertvar",   FieldFunc_Var,   1, nullptr);
               cdo_operator_add("vertvar1",  FieldFunc_Var1,  1, nullptr);
               cdo_operator_add("vertstd",   FieldFunc_Std,   1, nullptr);
               cdo_operator_add("vertstd1",  FieldFunc_Std1,  1, nullptr);
    // clang-format on

    operatorID = cdo_operator_id();
    operfunc = cdo_operator_f1(operatorID);
    needWeights = cdo_operator_f2(operatorID);

    lrange = (operfunc == FieldFunc_Range);
    lmean = (operfunc == FieldFunc_Mean || operfunc == FieldFunc_Avg);
    lstd = (operfunc == FieldFunc_Std || operfunc == FieldFunc_Std1);
    lvarstd = (lstd || operfunc == FieldFunc_Var || operfunc == FieldFunc_Var1);
    divisor = (operfunc == FieldFunc_Std1 || operfunc == FieldFunc_Var1);

    // int applyWeights = lmean;

    streamID1 = cdo_open_read(0);
    auto vlistID1 = cdo_stream_inq_vlist(streamID1);

    vlistClearFlag(vlistID1);
    nvars = vlistNvars(vlistID1);
    for (int varID = 0; varID < nvars; ++varID) vlistDefFlag(vlistID1, varID, 0, true);

    vlistID2 = vlistCreate();
    cdo_vlist_copy_flag(vlistID2, vlistID1);
    vlistDefNtsteps(vlistID2, vlistNtsteps(vlistID1));

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    auto surfID = get_surface_ID(vlistID1);
    set_surface_ID(vlistID2, surfID);

    nzaxis = vlistNzaxis(vlistID1);
    vert = std::vector<VertInfo>(nzaxis);
    if (needWeights)
      {
        auto useweights = true;
        auto genbounds = false;
        vertstat_get_parameter(useweights, genbounds);

        if (!useweights)
          {
            genbounds = false;
            cdo_print("Using constant vertical weights!");
          }

        for (int index = 0; index < nzaxis; ++index)
          {
            auto zaxisID = vlistZaxis(vlistID1, index);
            auto nlev = zaxisInqSize(zaxisID);
            vert[index].numLevels = 0;
            vert[index].status = 0;
            vert[index].zaxisID = zaxisID;
            // if (nlev > 1)
            {
              vert[index].numLevels = nlev;
              vert[index].thickness.resize(nlev);
              vert[index].weights.resize(nlev);
              vert[index].status
                  = get_layer_thickness(useweights, genbounds, index, zaxisID, nlev, vert[index].thickness, vert[index].weights);
            }
            if (!useweights) vert[index].status = 3;
          }
      }

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);

    varListInit(varList, vlistID1);

    vars1 = FieldVector(nvars);
    vars2 = FieldVector(nvars);
    samp1 = FieldVector(nvars);
    if (lvarstd || lrange) vars2.resize(nvars);

    for (int varID = 0; varID < nvars; ++varID)
      {
        const auto &var = varList[varID];
        samp1[varID].grid = var.gridID;
        samp1[varID].missval = var.missval;
        samp1[varID].memType = MemType::Double;
        vars1[varID].grid = var.gridID;
        vars1[varID].missval = var.missval;
        vars1[varID].memType = MemType::Double;
        vars1[varID].resize(var.gridsize);
        if (lvarstd || lrange)
          {
            vars2[varID].grid = var.gridID;
            vars2[varID].missval = var.missval;
            vars2[varID].memType = MemType::Double;
            vars2[varID].resize(var.gridsize);
          }
      }
  }

  void
  run()
  {

    auto field2_stdvar_func = lstd ? field2_std : field2_var;
    auto fieldc_stdvar_func = lstd ? fieldc_std : fieldc_var;
    int tsID = 0;
    while (true)
      {
        auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
        if (nrecs == 0) break;

        cdo_taxis_copy_timestep(taxisID2, taxisID1);
        cdo_def_timestep(streamID2, tsID);

        std::vector<bool> varsLevelInit(nvars, false);

        for (int recID = 0; recID < nrecs; ++recID)
          {
            int varID, levelID;
            cdo_inq_record(streamID1, &varID, &levelID);

            const auto &var = varList[varID];

            auto &rsamp1 = samp1[varID];
            auto &rvars1 = vars1[varID];

            rvars1.nsamp++;
            if (lrange) vars2[varID].nsamp++;

            auto gridsize = var.gridsize;

            auto layerWeight = 1.0;
            auto layerThickness = 1.0;
            if (needWeights)
              {
                for (int index = 0; index < nzaxis; ++index)
                  if (vert[index].zaxisID == var.zaxisID)
                    {
                      if (vert[index].status == 0 && tsID == 0 && levelID == 0 && var.nlevels > 1)
                        {
                          cdo_warning("Layer bounds not available, using constant vertical weights for variable %s!", var.name);
                        }
                      else
                        {
                          layerWeight = vert[index].weights[levelID];
                          layerThickness = vert[index].thickness[levelID];
                        }

                      break;
                    }
              }

            if (varsLevelInit[varID] == false)
              {
                varsLevelInit[varID] = true;
                cdo_read_record(streamID1, rvars1);
                if (lrange)
                  {
                    vars2[varID].nmiss = rvars1.nmiss;
                    vars2[varID].vec_d = rvars1.vec_d;
                  }

                if (operatorID == VERTINT && is_not_equal(layerThickness, 1.0)) fieldc_mul(rvars1, layerThickness);
                if (lmean && is_not_equal(layerWeight, 1.0)) fieldc_mul(rvars1, layerWeight);

                if (lvarstd)
                  {
                    if (is_not_equal(layerWeight, 1.0))
                      {
                        field2_moqw(vars2[varID], rvars1, layerWeight);
                        fieldc_mul(rvars1, layerWeight);
                      }
                    else { field2_moq(vars2[varID], rvars1); }
                  }

                if (rvars1.nmiss || !rsamp1.empty() || needWeights)
                  {
                    if (rsamp1.empty()) rsamp1.resize(gridsize);

                    for (size_t i = 0; i < gridsize; ++i)
                      rsamp1.vec_d[i] = (dbl_is_equal(rvars1.vec_d[i], rvars1.missval)) ? 0.0 : layerWeight;
                  }
              }
            else
              {
                field.init(var);
                cdo_read_record(streamID1, field);

                if (operatorID == VERTINT && is_not_equal(layerThickness, 1.0)) fieldc_mul(field, layerThickness);
                if (lmean && is_not_equal(layerWeight, 1.0)) fieldc_mul(field, layerWeight);

                if (field.nmiss || !rsamp1.empty())
                  {
                    if (rsamp1.empty()) rsamp1.resize(gridsize, rvars1.nsamp);

                    if (field.memType == MemType::Float)
                      {
                        for (size_t i = 0; i < gridsize; ++i)
                          if (!dbl_is_equal(field.vec_f[i], (float) rvars1.missval)) rsamp1.vec_d[i] += layerWeight;
                      }
                    else
                      {
                        for (size_t i = 0; i < gridsize; ++i)
                          if (!dbl_is_equal(field.vec_d[i], rvars1.missval)) rsamp1.vec_d[i] += layerWeight;
                      }
                  }

                if (lvarstd)
                  {
                    if (is_not_equal(layerWeight, 1.0))
                      {
                        field2_sumqw(vars2[varID], field, layerWeight);
                        field2_sumw(rvars1, field, layerWeight);
                      }
                    else { field2_sumsumq(rvars1, vars2[varID], field); }
                  }
                else if (lrange) { field2_maxmin(rvars1, vars2[varID], field); }
                else { field2_function(rvars1, field, operfunc); }
              }
          }

        for (int varID = 0; varID < nvars; ++varID)
          {
            const auto &rsamp1 = samp1[varID];
            auto &rvars1 = vars1[varID];

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
                      field2_stdvar_func(rvars1, vars2[varID], rsamp1, divisor);
                    else
                      fieldc_stdvar_func(rvars1, vars2[varID], rvars1.nsamp, divisor);
                  }
                else if (lrange) { field2_sub(rvars1, vars2[varID]); }

                cdo_def_record(streamID2, varID, 0);
                cdo_write_record(streamID2, rvars1);
                rvars1.nsamp = 0;
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

    cdo_finish();

  }
};

void *
Vertstat(void *process)
{
  ModuleVertstat vertstat;
  vertstat.init(process);
  vertstat.run();
  vertstat.close();

  return nullptr;
}
