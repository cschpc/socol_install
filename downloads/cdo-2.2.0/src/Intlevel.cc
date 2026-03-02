/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Intlevel   intlevel        Linear level interpolation
*/

#include <cdi.h>

#include "cimdOmp.h"
#include "cdo_options.h"
#include "process_int.h"
#include "cdo_vlist.h"
#include "cdo_zaxis.h"
#include "parse_literals.h"
#include "pmlist.h"
#include "param_conversion.h"

template <typename T>
static double
vert_interp_lev_kernel(float w1, float w2, const T var1L1, const T var1L2, const T missval)
{
  if (dbl_is_equal(var1L1, missval)) w1 = 0.0f;
  if (dbl_is_equal(var1L2, missval)) w2 = 0.0f;

  // clang-format off
  if      (is_equal(w1, 0.0f) && is_equal(w2, 0.0f)) return missval;
  else if (is_equal(w1, 0.0f)) return (w2 >= 0.5f) ? var1L2 : missval;
  else if (is_equal(w2, 0.0f)) return (w1 >= 0.5f) ? var1L1 : missval;
  else                         return var1L1 * (double)w1 + var1L2 * (double)w2;
  // clang-format on
}

constexpr int BottomLevel = 32000;
constexpr int TopLevel = 32001;

static void
restore_index_and_weights(int nlev1, int idx, float wgt, int &idx1, int &idx2, float &wgt1, float &wgt2)
{
  if (idx == BottomLevel)
    {
      idx1 = 0;
      idx2 = 0;
      wgt1 = 0.0f;
      wgt2 = wgt;
    }
  else if (idx == TopLevel)
    {
      idx1 = nlev1 - 1;
      idx2 = nlev1 - 1;
      wgt1 = wgt;
      wgt2 = 0.0f;
    }
  else
    {
      idx1 = (idx < 0) ? -idx : idx;
      idx2 = (idx < 0) ? idx1 - 1 : idx1 + 1;
      wgt1 = wgt;
      wgt2 = 1.0f - wgt;
    }
  // printf("%d %d %g %g\n", idx1, idx2, wgt1, wgt2);
}

//  1D vertical interpolation
template <typename T>
static void
vert_interp_lev(size_t gridsize, int nlev1, T missval, const Varray<T> &vardata1, Varray<T> &vardata2, int nlev2,
                const Varray<int> &lev_idx, const Varray<float> &lev_wgt)
{
  for (int ilev = 0; ilev < nlev2; ++ilev)
    {
      auto idx = lev_idx[ilev];
      auto wgt = lev_wgt[ilev];
      int idx1, idx2;
      float wgt1, wgt2;
      restore_index_and_weights(nlev1, idx, wgt, idx1, idx2, wgt1, wgt2);

      // upper/lower values from input field
      auto var1L1 = &vardata1[gridsize * idx1];
      auto var1L2 = &vardata1[gridsize * idx2];

      auto var2 = &vardata2[gridsize * ilev];

#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
      for (size_t i = 0; i < gridsize; ++i) { var2[i] = vert_interp_lev_kernel(wgt1, wgt2, var1L1[i], var1L2[i], missval); }
    }
}

static void
vert_interp_lev(size_t gridsize, int nlev1, double missval, const Field3D &field1, Field3D &field2, int nlev2,
                const Varray<int> &lev_idx, const Varray<float> &lev_wgt)
{
  if (field1.memType != field2.memType) cdo_abort("Interal error, memType of field1 and field2 differ!");

  if (field1.memType == MemType::Float)
    vert_interp_lev(gridsize, nlev1, (float) missval, field1.vec_f, field2.vec_f, nlev2, lev_idx, lev_wgt);
  else
    vert_interp_lev(gridsize, nlev1, missval, field1.vec_d, field2.vec_d, nlev2, lev_idx, lev_wgt);
}

//  3D vertical interpolation
template <typename T>
void
vert_interp_lev3d(size_t gridsize, int nlev1, T missval, const Varray<T> &vardata1, Varray<T> &vardata2, int nlev2,
                  const Varray<int> &lev_idx, const Varray<float> &lev_wgt)
{
  for (int ilev = 0; ilev < nlev2; ilev++)
    {
      auto offset = ilev * gridsize;
      auto var2 = &vardata2[offset];

#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
      for (size_t i = 0; i < gridsize; ++i)
        {
          auto idx = lev_idx[offset + i];
          auto wgt = lev_wgt[offset + i];
          int idx1, idx2;
          float wgt1, wgt2;
          restore_index_and_weights(nlev1, idx, wgt, idx1, idx2, wgt1, wgt2);

          // upper/lower values from input field
          auto var1L1 = vardata1[idx1 * gridsize + i];
          auto var1L2 = vardata1[idx2 * gridsize + i];

          var2[i] = vert_interp_lev_kernel(wgt1, wgt2, var1L1, var1L2, missval);
        }
    }
}

// Explicit instantiation
template void vert_interp_lev3d(size_t gridsize, int nlev1, float missval, const Varray<float> &vardata1, Varray<float> &vardata2,
                                int nlev2, const Varray<int> &lev_idx, const Varray<float> &lev_wgt);
template void vert_interp_lev3d(size_t gridsize, int nlev1, double missval, const Varray<double> &vardata1,
                                Varray<double> &vardata2, int nlev2, const Varray<int> &lev_idx, const Varray<float> &lev_wgt);

void
vert_interp_lev3d(size_t gridsize, int nlev1, double missval, const Field3D &field1, Field3D &field2, int nlev2,
                  const Varray<int> &lev_idx, const Varray<float> &lev_wgt)
{
  if (field1.memType != field2.memType) cdo_abort("Interal error, memType of field1 and field2 differ!");

  if (field1.memType == MemType::Float)
    vert_interp_lev3d(gridsize, nlev1, (float) missval, field1.vec_f, field2.vec_f, nlev2, lev_idx, lev_wgt);
  else
    vert_interp_lev3d(gridsize, nlev1, missval, field1.vec_d, field2.vec_d, nlev2, lev_idx, lev_wgt);
}

void
vert_gen_weights(int expol, int nlev1, const Varray<double> &lev1, int nlev2, const Varray<double> &lev2, Varray<int> &lev_idx,
                 Varray<float> &lev_wgt)
{
  for (int i2 = 0; i2 < nlev2; ++i2)
    {
      int idx1 = 0, idx2 = 0;
      double val1, val2 = 0.0;

      // Because 2 levels were added to the source vertical coordinate (one on top, one at the bottom), its loop starts at 1
      int i1;
      for (i1 = 1; i1 < nlev1; ++i1)
        {
          auto lev1_isUp = (lev1[i1 - 1] < lev1[i1]);
          idx1 = lev1_isUp ? i1 - 1 : i1;
          idx2 = lev1_isUp ? i1 : i1 - 1;
          val1 = lev1[idx1];
          val2 = lev1[idx2];
          if (lev2[i2] > val1 && lev2[i2] <= val2) break;
        }

      if (i1 == nlev1) cdo_abort("Level %g not found!", lev2[i2]);

      if (i1 - 1 == 0)  // destination levels is not covert by the first two input z levels
        {
          lev_idx[i2] = BottomLevel;
          lev_wgt[i2] = static_cast<float>(expol || is_equal(lev2[i2], val2));
        }
      else if (i1 == nlev1 - 1)  // destination level is beyond the last value of the input z field
        {
          lev_idx[i2] = TopLevel;
          lev_wgt[i2] = static_cast<float>(expol || is_equal(lev2[i2], val2));
        }
      else  // target z values has two bounday values in input z field
        {
          lev_idx[i2] = idx1 - 1;
          if (idx1 > idx2) lev_idx[i2] = -lev_idx[i2];
          lev_wgt[i2] = (lev1[idx2] - lev2[i2]) / (lev1[idx2] - lev1[idx1]);
        }
      // printf("%d %g %d %d %g %g %d %g\n", i2, lev2[i2], idx1, idx2, lev1[idx1], lev1[idx2], lev_idx[i2], lev_wgt[i2]);
    }
}

bool
levelDirUp(const int nlev, const double *const lev)
{
  auto lup = (nlev > 1 && lev[1] > lev[0]);
  for (int k = 1; k < nlev - 1; ++k)
    if (lup && lev[k + 1] <= lev[k]) return false;

  return lup;
}

bool
levelDirDown(const int nlev, const double *const lev)
{
  auto ldown = (nlev > 1 && lev[1] < lev[0]);
  for (int k = 1; k < nlev - 1; ++k)
    if (ldown && lev[k + 1] >= lev[k]) return false;

  return ldown;
}

template <typename T>
static void
vert_gen_weights3d1d(bool expol, size_t gridsize, int nlev1, const Varray<T> &xlev1, int nlev2, const std::vector<double> &lev2,
                     Varray<int> &xlev_idx, Varray<float> &xlev_wgt)
{
  auto nthreads = Threading::ompNumThreads;
  Varray2D<double> lev1p2(nthreads, Varray<double>(nlev1 + 2));
  Varray2D<float> lev_wgt(nthreads, Varray<float>(nlev2));
  Varray2D<int> lev_idx(nthreads, Varray<int>(nlev2));

  // Check monotony of vertical levels
  for (int k = 0; k < nlev1; ++k) lev1p2[0][k] = xlev1[k * gridsize];
  auto lup = levelDirUp(nlev1, lev1p2[0].data());
  auto ldown = levelDirDown(nlev1, lev1p2[0].data());
  if (!lup && !ldown) cdo_abort("Non monotonic zaxis!");
  double level_0 = lup ? -1.e33 : 1.e33;
  double level_N = lup ? 1.e33 : -1.e33;

#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
  for (size_t i = 0; i < gridsize; ++i)
    {
      auto ompthID = cdo_omp_get_thread_num();

      lev1p2[ompthID][0] = level_0;
      lev1p2[ompthID][nlev1 + 1] = level_N;
      for (int k = 0; k < nlev1; ++k) lev1p2[ompthID][k + 1] = xlev1[k * gridsize + i];

      vert_gen_weights(expol, nlev1 + 2, lev1p2[ompthID], nlev2, lev2, lev_idx[ompthID], lev_wgt[ompthID]);

      for (int k = 0; k < nlev2; ++k) xlev_idx[k * gridsize + i] = lev_idx[ompthID][k];
      for (int k = 0; k < nlev2; ++k) xlev_wgt[k * gridsize + i] = lev_wgt[ompthID][k];
    }
}

static void
vert_gen_weights3d1d(bool expol, size_t gridsize, int nlev1, Field3D &field1, int nlev2, const std::vector<double> &lev2,
                     Varray<int> &lev_idx, Varray<float> &lev_wgt)
{
  if (field1.memType == MemType::Float)
    vert_gen_weights3d1d(expol, gridsize, nlev1, field1.vec_f, nlev2, lev2, lev_idx, lev_wgt);
  else
    vert_gen_weights3d1d(expol, gridsize, nlev1, field1.vec_d, nlev2, lev2, lev_idx, lev_wgt);
}

static int
create_zaxis_from_zvar(const std::vector<double> &levels, int vlistID, int varID)
{
  auto nlevels = levels.size();
  auto zaxisID = zaxisCreate(ZAXIS_GENERIC, nlevels);

  std::string name = "zlev";
  cdiDefKeyString(zaxisID, CDI_GLOBAL, CDI_KEY_NAME, name.c_str());
  auto longname = cdo::inq_var_longname(vlistID, varID);
  if (longname.size()) cdiDefKeyString(zaxisID, CDI_GLOBAL, CDI_KEY_LONGNAME, longname.c_str());
  auto units = cdo::inq_var_units(vlistID, varID);
  if (units.size()) cdiDefKeyString(zaxisID, CDI_GLOBAL, CDI_KEY_UNITS, units.c_str());

  zaxisDefLevels(zaxisID, levels.data());

  return zaxisID;
}

static int
create_zaxis_from_zaxis(const std::vector<double> &levels, int zaxisID1)
{
  auto nlevels = levels.size();
  auto zaxisID2 = zaxisCreate(zaxisInqType(zaxisID1), nlevels);

  cdiCopyKey(zaxisID1, CDI_GLOBAL, CDI_KEY_NAME, zaxisID2);
  cdiCopyKey(zaxisID1, CDI_GLOBAL, CDI_KEY_LONGNAME, zaxisID2);
  cdiCopyKey(zaxisID1, CDI_GLOBAL, CDI_KEY_UNITS, zaxisID2);
  cdiCopyKey(zaxisID1, CDI_GLOBAL, CDI_KEY_DATATYPE, zaxisID2);

  zaxisDefLevels(zaxisID2, levels.data());

  return zaxisID2;
}

static std::vector<double>
intlevel_get_parameter(std::string &zdescription, std::string &zvarname)
{
  std::vector<double> levels;

  auto pargc = cdo_operator_argc();
  if (pargc)
    {
      const auto &pargv = cdo_get_oper_argv();

      KVList kvlist;
      kvlist.name = cdo_module_name();
      if (kvlist.parse_arguments(pargc, pargv) != 0) cdo_abort("Parse error!");
      if (Options::cdoVerbose) kvlist.print();

      auto def_level = false;
      auto def_zdescription = false;
      for (const auto &kv : kvlist)
        {
          const auto &key = kv.key;
          // if (kv.nvalues > 1) cdo_abort("Too many values for parameter key >%s<!", key);
          if (kv.nvalues < 1) cdo_abort("Missing value for parameter key >%s<!", key);
          const auto &values = kv.values;
          const auto &value = kv.values[0];
          int nvalues = kv.nvalues;
          if (nvalues == 1 && value.empty()) nvalues = 0;

          // clang-format off
         if (key == "level")
            {
              if (def_zdescription) cdo_abort("Parameter level and zdescription can't be mixed!");
              levels.resize(nvalues);
              for (int i = 0; i < nvalues; ++i) levels[i] = literal_to_double(values[i]);
              def_level = true;
            }
          else if (key == "zdescription")
            {
              if (def_level) cdo_abort("Parameter zdescription and level can't be mixed!");
              zdescription = value;
              def_zdescription = true;
            }
          else if (key == "zvarname")
            {
              zvarname = value;
            }
          else cdo_abort("Invalid parameter key >%s<!", key);
          // clang-format on
        }
    }

  return levels;
}

static void
handle_zvar(size_t &zvarGridsize, size_t &wisize, int nvars, const std::string &zvarname, const VarList &varList1, int &zvarID,
            bool &zvarIsVarying, int &nlevel, const std::vector<double> &lev2, int &nlev1, const int &nlev2, const int &vlistID1,
            int &zaxisID2, int &zaxisID1)
{
  for (int varID = 0; varID < nvars; ++varID)
    {
      if (zvarname == varList1[varID].name)
        {
          zvarID = varID;
          break;
        }
    }

  if (zvarID == CDI_UNDEFID) cdo_abort("Variable %s not found!", zvarname);
  zvarIsVarying = (varList1[zvarID].timetype == TIME_VARYING);
  zvarGridsize = varList1[zvarID].gridsize;
  nlev1 = varList1[zvarID].nlevels;

  if (zaxisID2 == CDI_UNDEFID) zaxisID2 = create_zaxis_from_zvar(lev2, vlistID1, zvarID);

  wisize = zvarGridsize * nlev2;

  int i;
  auto nzaxis = vlistNzaxis(vlistID1);
  for (i = 0; i < nzaxis; ++i)
    {
      auto zaxisID = vlistZaxis(vlistID1, i);
      nlevel = zaxisInqSize(zaxisID);
      if (nlevel == nlev1)
        {
          zaxisID1 = zaxisID;
          break;
        }
    }
  if (i == nzaxis) cdo_abort("No processable variable found!");
}

static void
handle_empty_zvar(size_t &wisize, int &nlevel, std::vector<double> &lev1, const std::vector<double> &lev2, int &nlev1,
                  const int &nlev2, const int &vlistID1, int &zaxisID2, int &zaxisID1)
{
  int i;
  auto nzaxis = vlistNzaxis(vlistID1);
  for (i = 0; i < nzaxis; ++i)
    {
      auto zaxisID = vlistZaxis(vlistID1, i);
      nlevel = zaxisInqSize(zaxisID);
      // if (zaxisInqType(zaxisID) != ZAXIS_HYBRID && zaxisInqType(zaxisID) != ZAXIS_HYBRID_HALF)
      if (nlevel > 1)
        {
          zaxisID1 = zaxisID;
          break;
        }
    }
  if (i == nzaxis) cdo_abort("No processable variable found!");

  if (zaxisID2 == CDI_UNDEFID) zaxisID2 = create_zaxis_from_zaxis(lev2, zaxisID1);

  nlev1 = nlevel;
  lev1.resize(nlev1 + 2);
  cdo_zaxis_inq_levels(zaxisID1, &lev1[1]);

  auto lup = levelDirUp(nlev1, &lev1[1]);
  auto ldown = levelDirDown(nlev1, &lev1[1]);
  if (!lup && !ldown) cdo_abort("Non monotonic zaxis!");
  lev1[0] = lup ? -1.e33 : 1.e33;
  lev1[nlev1 + 1] = lup ? 1.e33 : -1.e33;

  if (Options::cdoVerbose)
    for (i = 0; i < nlev1 + 2; ++i) cdo_print("lev1 %d: %g", i, lev1[i]);

  wisize = nlev2;
}

void *
Intlevel(void *process)
{
  int zaxisID1 = -1;
  int nlevel = 0;

  cdo_initialize(process);

  // clang-format off
  auto INTLEVEL  = cdo_operator_add("intlevel",  0, 0, nullptr);
  auto INTLEVELX = cdo_operator_add("intlevelx", 0, 0, nullptr);
  // clang-format on

  (void) (INTLEVEL);  // unused

  auto operatorID = cdo_operator_id();

  auto expol = (operatorID == INTLEVELX);

  operator_input_arg("level|zdescription, [zvarname]");

  std::vector<double> lev2;
  int zaxisID2 = CDI_UNDEFID;
  std::string zvarname;
  auto argv = cdo_get_oper_argv();

  if (isdigit((int) argv[0][0])) { lev2 = cdo_argv_to_flt(argv); }
  else
    {
      std::string zdescription;
      lev2 = intlevel_get_parameter(zdescription, zvarname);
      if (!zdescription.empty())
        {
          auto zfilename = zdescription.c_str();
          auto zfp = std::fopen(zfilename, "r");
          if (!zfp) cdo_abort("Open failed on %s", zfilename);
          zaxisID2 = zaxis_from_file(zfp, zfilename);
          std::fclose(zfp);
          if (zaxisID2 == CDI_UNDEFID) cdo_abort("Invalid zaxis description file %s!", zfilename);
          auto nlevels = zaxisInqSize(zaxisID2);
          lev2.resize(nlevels);
          zaxisInqLevels(zaxisID2, lev2.data());
        }
    }

  const int nlev2 = lev2.size();

  if (Options::cdoVerbose)
    for (int i = 0; i < nlev2; ++i) cdo_print("lev2 %d: %g", i, lev2[i]);

  auto streamID1 = cdo_open_read(0);

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto vlistID2 = vlistDuplicate(vlistID1);

  auto taxisID1 = vlistInqTaxis(vlistID1);
  auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  VarList varList1;
  varListInit(varList1, vlistID1);
  varListSetUniqueMemtype(varList1);
  auto memType = varList1[0].memType;

  auto nvars = vlistNvars(vlistID1);

  // Find z-variable
  int nlev1 = 0;
  Varray<double> lev1;
  int zvarID = CDI_UNDEFID;
  bool zvarIsVarying = false;
  size_t zvarGridsize = 0;
  size_t wisize = 0;

  if (zvarname.empty()) { handle_empty_zvar(wisize, nlevel, lev1, lev2, nlev1, nlev2, vlistID1, zaxisID2, zaxisID1); }
  else
    {
      handle_zvar(zvarGridsize, wisize, nvars, zvarname, varList1, zvarID, zvarIsVarying, nlevel, lev2, nlev1, nlev2, vlistID1,
                  zaxisID2, zaxisID1);
    }

  auto nzaxis = vlistNzaxis(vlistID1);
  for (int index = 0; index < nzaxis; ++index)
    {
      auto zaxisID = vlistZaxis(vlistID1, index);
      auto nlevels = zaxisInqSize(zaxisID);
      if (zaxisID == zaxisID1 || nlevels == nlev1) vlistChangeZaxisIndex(vlistID2, index, zaxisID2);
    }

  VarList varList2;
  varListInit(varList2, vlistID2);
  varListSetMemtype(varList2, memType);

  Varray<int> lev_idx(wisize);
  Varray<float> lev_wgt(wisize);

  auto streamID2 = cdo_open_write(1);
  cdo_def_vlist(streamID2, vlistID2);

  std::vector<bool> processVars(nvars);
  std::vector<bool> interpVars(nvars);
  std::vector<std::vector<size_t>> varnmiss(nvars);
  Field3DVector vardata1(nvars), vardata2(nvars);

  const int maxlev = (nlev1 > nlev2) ? nlev1 : nlev2;

  for (int varID = 0; varID < nvars; ++varID)
    {
      auto nlevels = varList1[varID].nlevels;

      vardata1[varID].init(varList1[varID]);

      interpVars[varID] = (varList1[varID].zaxisID == zaxisID1 || nlevels == nlev1);

      if (interpVars[varID])
        {
          varnmiss[varID].resize(maxlev, 0);
          vardata2[varID].init(varList2[varID]);
        }
      else { varnmiss[varID].resize(nlevels); }
    }

  int tsID = 0;
  while (true)
    {
      auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
      if (nrecs == 0) break;

      for (int varID = 0; varID < nvars; ++varID) processVars[varID] = false;

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_def_timestep(streamID2, tsID);

      for (int recID = 0; recID < nrecs; ++recID)
        {
          int varID, levelID;
          cdo_inq_record(streamID1, &varID, &levelID);
          cdo_read_record(streamID1, vardata1[varID], levelID, &varnmiss[varID][levelID]);
          processVars[varID] = true;
        }

      if (tsID == 0 || zvarIsVarying)
        {
          if (!zvarname.empty())
            vert_gen_weights3d1d(expol, zvarGridsize, nlev1, vardata1[zvarID], nlev2, lev2, lev_idx, lev_wgt);
          else
            vert_gen_weights(expol, nlev1 + 2, lev1, nlev2, lev2, lev_idx, lev_wgt);
        }

      for (int varID = 0; varID < nvars; ++varID)
        {
          if (processVars[varID] && interpVars[varID])
            {
              auto missval = varList1[varID].missval;
              auto gridsize = varList1[varID].gridsize;

              if (!zvarname.empty())
                vert_interp_lev3d(gridsize, nlev1, missval, vardata1[varID], vardata2[varID], nlev2, lev_idx, lev_wgt);
              else
                vert_interp_lev(gridsize, nlev1, missval, vardata1[varID], vardata2[varID], nlev2, lev_idx, lev_wgt);

              for (int levelID = 0; levelID < nlev2; ++levelID)
                {
                  auto offset = gridsize * levelID;
                  if (memType == MemType::Float)
                    varnmiss[varID][levelID] = array_num_mv(gridsize, &vardata2[varID].vec_f[offset], (float) missval);
                  else
                    varnmiss[varID][levelID] = array_num_mv(gridsize, &vardata2[varID].vec_d[offset], missval);
                }
            }
        }

      for (int varID = 0; varID < nvars; ++varID)
        {
          if (processVars[varID])
            {
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
