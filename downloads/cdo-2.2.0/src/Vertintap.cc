/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Vertint    ap2pl           Model air pressure level to pressure level interpolation
*/

#include "cdo_options.h"
#include "process_int.h"
#include "cdo_vlist.h"
#include "field_vinterp.h"
#include "stdnametable.h"
#include "util_string.h"
#include "const.h"
#include "cdo_zaxis.h"
#include "param_conversion.h"
#include "vertint_util.h"

static void
check_range_ps(int stepNum, const Field &psProg)
{
  auto mm = field_min_max(psProg);
  if (mm.min < MIN_PS || mm.max > MAX_PS)
    cdo_warning("Surface pressure out of range (min=%g max=%g) [timestep:%d]!", mm.min, mm.max, stepNum);
}

static bool
is_height_axis(int zaxisID)
{
  auto isHeight = false;
  if (zaxisInqType(zaxisID) == ZAXIS_REFERENCE)
    {
      // auto units = cdo::inq_key_string(zaxisID, CDI_GLOBAL, CDI_KEY_UNITS);
      auto stdname = cdo::inq_key_string(zaxisID, CDI_GLOBAL, CDI_KEY_STDNAME);
      // if (stdname == "height" && units.empty()) isHeight = true;
      if (stdname == "height") isHeight = true;
    }
  return isHeight;
}

static void
change_height_zaxis(int nhlev, int vlistID1, int vlistID2, int zaxisID2)
{
  auto nzaxis = vlistNzaxis(vlistID1);
  for (int iz = 0; iz < nzaxis; ++iz)
    {
      auto zaxisID = vlistZaxis(vlistID1, iz);
      auto nlevel = zaxisInqSize(zaxisID);
      if ((nlevel == nhlev || nlevel == (nhlev + 1)) && is_height_axis(zaxisID)) vlistChangeZaxisIndex(vlistID2, iz, zaxisID2);
    }
}

template <typename T>
static void
calc_half_press(size_t gridsize, size_t nhlevf, const Varray<T> &fullPress, size_t nhlevh, Varray<T> &halfPress)
{
  for (size_t i = 0; i < gridsize; ++i) halfPress[i] = 0;
#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
  for (size_t k = 1; k < nhlevf; ++k)
    {
      auto fullPress_km1 = &fullPress[(k - 1) * gridsize];
      auto fullPress_k = &fullPress[k * gridsize];
      auto halfPress_k = &halfPress[k * gridsize];
      for (size_t i = 0; i < gridsize; ++i) halfPress_k[i] = 0.5 * (fullPress_km1[i] + fullPress_k[i]);
    }
  for (size_t i = 0; i < gridsize; ++i) halfPress[(nhlevh - 1) * gridsize + i] = fullPress[(nhlevf - 1) * gridsize + i];
}

static void
calc_half_press(const Field3D &fullPress, Field3D &halfPress)
{
  if (fullPress.memType == MemType::Float)
    calc_half_press(fullPress.gridsize, fullPress.nlevels, fullPress.vec_f, halfPress.nlevels, halfPress.vec_f);
  else
    calc_half_press(fullPress.gridsize, fullPress.nlevels, fullPress.vec_d, halfPress.nlevels, halfPress.vec_d);
}

void *
Vertintap(void *process)
{
  enum
  {
    func_pl,
    func_hl
  };
  int nhlev = 0, nhlevf = 0, nhlevh = 0;
  int apressID = -1, dpressID = -1;
  int psID = -1;

  cdo_initialize(process);

  // clang-format off
                   cdo_operator_add("ap2pl",     func_pl, 0, "pressure levels in pascal");
  auto AP2PLX    = cdo_operator_add("ap2plx",    func_pl, 0, "pressure levels in pascal");
                   cdo_operator_add("ap2hl",     func_hl, 0, "height levels in meter");
  auto AP2HLX    = cdo_operator_add("ap2hlx",    func_hl, 0, "height levels in meter");
  // clang-format on

  auto operatorID = cdo_operator_id();
  auto useHeightLevel = (cdo_operator_f1(operatorID) == func_hl);

  auto extrapolate = (operatorID == AP2PLX || operatorID == AP2HLX);
  if (extrapolate == false) extrapolate = getenv_extrapolate();

  operator_input_arg(cdo_operator_enter(operatorID));

  Varray<double> plev;
  if (cdo_operator_argc() == 1 && cdo_operator_argv(0) == "default")
    {
      if (useHeightLevel)
        plev = { 10, 50, 100, 500, 1000, 5000, 10000, 15000, 20000, 25000, 30000 };
      else
        plev
            = { 100000, 92500, 85000, 70000, 60000, 50000, 40000, 30000, 25000, 20000, 15000, 10000, 7000, 5000, 3000, 2000, 1000 };
    }
  else { plev = cdo_argv_to_flt(cdo_get_oper_argv()); }

  int nplev = plev.size();

  auto streamID1 = cdo_open_read(0);

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto vlistID2 = vlistDuplicate(vlistID1);

  auto taxisID1 = vlistInqTaxis(vlistID1);
  auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  auto gridsize = vlist_check_gridsize(vlistID1);

  auto zaxistype = useHeightLevel ? ZAXIS_HEIGHT : ZAXIS_PRESSURE;
  auto zaxisIDp = zaxisCreate(zaxistype, nplev);
  zaxisDefLevels(zaxisIDp, plev.data());

  VarList varList1;
  varListInit(varList1, vlistID1);
  varListSetUniqueMemtype(varList1);
  auto memtype = varList1[0].memType;

  auto nvars = vlistNvars(vlistID1);

  for (int varID = 0; varID < nvars; ++varID)
    {
      auto stdname = string_to_lower(cdo::inq_key_string(vlistID1, varID, CDI_KEY_STDNAME));

      // clang-format off
      if      (stdname == var_stdname(surface_air_pressure)) psID = varID;
      else if (stdname == var_stdname(air_pressure))         apressID = varID;
      else if (stdname == var_stdname(pressure_thickness))   dpressID = varID;
      // clang-format on
    }

  if (Options::cdoVerbose)
    {
      cdo_print("Found:");
      // clang-format off
      if (-1 != psID)     cdo_print("  %s -> %s", var_stdname(surface_air_pressure), varList1[psID].name);
      if (-1 != apressID) cdo_print("  %s -> %s", var_stdname(air_pressure), varList1[apressID].name);
      if (-1 != dpressID) cdo_print("  %s -> %s", var_stdname(pressure_thickness), varList1[dpressID].name);
      // clang-format on
    }

  if (apressID == -1) cdo_abort("%s not found!", var_stdname(air_pressure));

  int zaxisIDh = -1;
  auto nzaxis = vlistNzaxis(vlistID1);
  for (int i = 0; i < nzaxis; ++i)
    {
      auto zaxisID = vlistZaxis(vlistID1, i);
      if (zaxisID == varList1[apressID].zaxisID)
        {
          auto mono_level = true;
          auto nlevels = zaxisInqSize(zaxisID);

          if (nlevels > 1 && is_height_axis(zaxisID))
            {
              Varray<double> level(nlevels);
              cdo_zaxis_inq_levels(zaxisID, &level[0]);
              int l;
              for (l = 0; l < nlevels; ++l)
                {
                  if ((l + 1) != (int) (level[l] + 0.5)) break;
                }
              if (l == nlevels) mono_level = true;
            }

          if (nlevels > 1 && is_height_axis(zaxisID) && mono_level)
            {
              zaxisIDh = zaxisID;
              nhlev = nlevels;
              nhlevf = nhlev;
              nhlevh = nhlevf + 1;

              break;
            }
        }
    }

  change_height_zaxis(nhlev, vlistID1, vlistID2, zaxisIDp);

  VarList varList2;
  varListInit(varList2, vlistID2);
  varListSetMemtype(varList2, memtype);

  std::vector<bool> processVars(nvars), interpVars(nvars);
  Varray2D<size_t> varnmiss(nvars);
  Field3DVector vardata1(nvars), vardata2(nvars);

  auto maxlev = (nhlevh > nplev) ? nhlevh : nplev;

  Varray<size_t> pnmiss;
  if (!extrapolate) pnmiss.resize(nplev);

  // check levels
  if (zaxisIDh != -1)
    {
      auto nlev = zaxisInqSize(zaxisIDh);
      if (nlev != nhlev) cdo_abort("Internal error, wrong number of height level!");
    }

  std::vector<int> vertIndex;
  Field psProg;
  Field3D fullPress, halfPress;
  if (zaxisIDh != -1 && gridsize > 0)
    {
      vertIndex.resize(gridsize * nplev);

      CdoVar var3Dfull, var3Dhalf;
      var3Dfull.gridsize = gridsize;
      var3Dfull.nlevels = nhlevf;
      var3Dfull.memType = memtype;
      fullPress.init(var3Dfull);

      var3Dhalf.gridsize = gridsize;
      var3Dhalf.nlevels = nhlevh;
      var3Dhalf.memType = memtype;
      halfPress.init(var3Dhalf);
    }
  else
    {
      cdo_warning("No 3D variable with generalized height levels found!");
      cdo_print(
          "Generalized height levels are defined by the attributes: standard_name=\"height\" and long_name=\"generalized height\"");
    }

  if (useHeightLevel)
    {
      Varray<double> phlev(nplev);
      height_to_pressure(phlev.data(), plev.data(), nplev);

      if (Options::cdoVerbose)
        for (int i = 0; i < nplev; ++i) cdo_print("level = %d   height = %g   pressure = %g", i + 1, plev[i], phlev[i]);

      plev = phlev;
    }

  for (int varID = 0; varID < nvars; ++varID)
    {
      const auto &var = varList1[varID];
      auto isHeightAxis = is_height_axis(var.zaxisID);

      if (gridInqType(var.gridID) == GRID_SPECTRAL) cdo_abort("Spectral data unsupported!");

      vardata1[varID].init(var);

      interpVars[varID]
          = (var.zaxisID == zaxisIDh || (isHeightAxis && zaxisIDh != -1 && (var.nlevels == nhlevh || var.nlevels == nhlevf)));

      if (interpVars[varID])
        {
          varnmiss[varID].resize(maxlev, 0);
          vardata2[varID].init(varList2[varID]);
        }
      else
        {
          if (isHeightAxis && zaxisIDh != -1 && var.nlevels > 1)
            cdo_warning("Parameter %d has wrong number of levels, skipped! (name=%s nlevel=%d)", varID + 1, var.name,
                        var.nlevels);

          varnmiss[varID].resize(var.nlevels);
        }
    }

  if (zaxisIDh != -1 && psID == -1)
    {
      if (dpressID != -1)
        cdo_warning("Surface pressure not found - set to vertical sum of %s!", var_stdname(pressure_thickness));
      else
        cdo_warning("Surface pressure not found - set to lower bound of %s!", var_stdname(air_pressure));
    }

  for (int varID = 0; varID < nvars; ++varID)
    {
      if (interpVars[varID] && varList1[varID].isConstant) vlistDefVarTimetype(vlistID2, varID, TIME_VARYING);
    }

  auto streamID2 = cdo_open_write(1);

  cdo_def_vlist(streamID2, vlistID2);

  int tsID = 0;
  while (true)
    {
      auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
      if (nrecs == 0) break;

      for (int varID = 0; varID < nvars; ++varID)
        {
          processVars[varID] = false;
          const auto &var = varList1[varID];
          for (int levelID = 0; levelID < var.nlevels; ++levelID) varnmiss[varID][levelID] = 0;
        }

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_def_timestep(streamID2, tsID);

      for (int recID = 0; recID < nrecs; ++recID)
        {
          int varID, levelID;
          cdo_inq_record(streamID1, &varID, &levelID);
          cdo_read_record(streamID1, vardata1[varID], levelID, &varnmiss[varID][levelID]);
          processVars[varID] = true;
        }

      for (int varID = 0; varID < nvars; ++varID)
        if (interpVars[varID]) processVars[varID] = true;

      if (zaxisIDh != -1)
        {
          if (tsID == 1 && varList1[apressID].timetype == TIME_CONSTANT)
            cdo_warning("%s does not vary in time!", var_stdname(air_pressure));

          if (psID != -1)
            {
              psProg.init(varList1[psID]);
              field_copy(vardata1[psID], psProg);
            }
          else if (dpressID != -1)
            {
              psProg.init(varList1[dpressID]);
              field_fill(psProg, 0);
              for (int k = 0; k < nhlevf; ++k) field_add(psProg, vardata1[dpressID], k);
            }
          else
            {
              psProg.init(varList1[apressID]);
              field_copy(vardata1[apressID], nhlevf - 1, psProg);
            }

          // check range of psProg
          check_range_ps(tsID + 1, psProg);

          field_copy(vardata1[apressID], fullPress);

          calc_half_press(fullPress, halfPress);

          gen_vert_index(vertIndex, plev, fullPress, gridsize);
          if (!extrapolate) gen_vert_index_mv(vertIndex, plev, gridsize, psProg, pnmiss);
        }

      for (int varID = 0; varID < nvars; ++varID)
        {
          if (processVars[varID])
            {
              const auto &var = varList1[varID];

              if (tsID > 0 && !interpVars[varID] && var.isConstant) continue;

              if (interpVars[varID])
                {
                  if (var.nlevels != nhlevf && var.nlevels != nhlevh)
                    cdo_abort("Number of generalized height level differ from full/half level (param=%s)!", var.name);

                  for (int levelID = 0; levelID < var.nlevels; ++levelID)
                    {
                      if (varnmiss[varID][levelID]) cdo_abort("Missing values unsupported for this operator!");
                    }

                  const auto &levels3D = (var.nlevels == nhlevf) ? fullPress : halfPress;
                  // vertIndex on half levels missing; do we need halfPress???
                  vertical_interp_X(levels3D, vardata1[varID], vardata2[varID], vertIndex, plev, gridsize);

                  if (!extrapolate) varray_copy(nplev, pnmiss, varnmiss[varID]);
                }

              for (int levelID = 0; levelID < varList2[varID].nlevels; ++levelID)
                {
                  cdo_def_record(streamID2, varID, levelID);
                  cdo_write_record(streamID2, interpVars[varID] ? vardata2[varID] : vardata1[varID], levelID,
                                   varnmiss[varID][levelID]);
                }
            }
        }

      tsID++;
    }

  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  cdo_finish();

  return nullptr;
}
