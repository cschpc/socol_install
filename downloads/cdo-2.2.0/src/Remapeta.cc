/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Remapeta     remapeta          Model to model level interpolation
*/

#include <cdi.h>

#include "process_int.h"
#include "cdo_vlist.h"
#include "readline.h"
#include "hetaeta.h"
#include "vertical_interp.h"
#include "stdnametable.h"
#include "util_string.h"
#include "timer.h"
#include "const.h"
#include "cdo_options.h"
#include "cdo_zaxis.h"
#include "cdi_lockedIO.h"

template <typename T>
static void
setmissval(long nvals, const Varray<int> &imiss, double missval, T *array)
{
  if (!imiss.empty())
    for (long i = 0; i < nvals; ++i)
      if (imiss[i]) array[i] = missval;
}

template <typename T>
static void
corr_hum(long gridsize, T *q, double q_min)
{
  for (long i = 0; i < gridsize; ++i)
    {
      if (q[i] < q_min) q[i] = q_min;
    }
}

static long
ncctop(double cptop, long nlev, long nlevp1, const double *vct_a, const double *vct_b)
{
  /*
    Description:
    Defines highest level *ncctop* where condensation is allowed.

    Author:

    E. Roeckner, MPI, October 2001
  */
  long nctop = 0;
  Varray<double> zph(nlevp1), zp(nlev);
  // double    cptop  =  1000.;   /* min. pressure level for cond. */

  // half level pressure values, assuming 101320. Pa surface pressure

  for (long jk = 0; jk < nlevp1; ++jk)
    {
      auto za = vct_a[jk];
      auto zb = vct_b[jk];
      zph[jk] = za + zb * 101320.;
    }

  // full level pressure

  for (long jk = 0; jk < nlev; ++jk) zp[jk] = (zph[jk] + zph[jk + 1]) * 0.5;

  // search for pressure level cptop (Pa)

  for (long jk = 0; jk < nlev; ++jk)
    {
      nctop = jk;
      if (zp[jk] >= cptop) break;
    }

  return nctop;
}

static Varray<double>
vctFromFile(const char *filename)
{
  char line[1024], *pline;
  int i = 0;
  constexpr int maxvct = 8192;

  auto fp = std::fopen(filename, "r");
  if (fp == nullptr)
    {
      perror(filename);
      exit(EXIT_FAILURE);
    }

  Varray<double> vct;
  vct.resize(maxvct);

  while (cdo::readline(fp, line, 1024))
    {
      if (line[0] == '#' || line[0] == '\0') continue;

      pline = line;
      auto num = (int) strtod(pline, &pline);
      if (pline == nullptr) cdo_abort("Format error in VCT file %s!", filename);
      if (num != i) cdo_warning("Inconsistent VCT file, entry %d is %d.", i, num);

      if (i + maxvct / 2 >= maxvct - 1) cdo_abort("Too many values in VCT file!");

      vct[i] = strtod(pline, &pline);
      if (pline == nullptr) cdo_abort("Format error in VCT file %s!", filename);

      vct[i + maxvct / 2] = strtod(pline, &pline);

      i++;
    }

  std::fclose(fp);

  auto nvct = 2 * i;
  auto nlevh = i - 1;

  for (i = 0; i < nlevh + 1; ++i) vct[i + nvct / 2] = vct[i + maxvct / 2];

  vct.resize(nvct);

  return vct;
}

template <typename T>
static void
vertSum(Varray<double> &sum, const Varray<T> &var3d, size_t gridsize, size_t nlevels)
{
  for (size_t i = 0; i < gridsize; ++i) sum[i] = 0;

  for (size_t k = 0; k < nlevels; ++k)
    for (size_t i = 0; i < gridsize; ++i) { sum[i] += var3d[k * gridsize + i]; }
}

template <typename T>
static void
vertSumw(Varray<double> &sum, const Varray<T> &var3d, size_t gridsize, size_t nlevels, const Varray<double> &deltap)
{
  for (size_t i = 0; i < gridsize; ++i) sum[i] = 0;

  for (size_t k = 0; k < nlevels; ++k)
    for (size_t i = 0; i < gridsize; ++i) { sum[i] += var3d[k * gridsize + i] * deltap[k * gridsize + i]; }
}

template <typename T>
void
field_copy_array(size_t len, const Field &field, T *array)
{
  if (field.memType == MemType::Float)
    for (size_t i = 0; i < len; ++i) array[i] = field.vec_f[i];
  else
    for (size_t i = 0; i < len; ++i) array[i] = field.vec_d[i];
}

#define MAX_VARS3D 1024

template <typename T>
void
remapeta(MemType memType)
{
  constexpr double cconst = 1.0E-6;
  size_t nfis2gp = 0;
  int nvars3D = 0;
  Varray<double> fis2;
  size_t nmissout = 0;
  bool lfis2 = false;
  int varids[MAX_VARS3D];
  Varray<int> imiss;
  double missval = 0.0;
  double cptop = 0.0;  // min. pressure level for cond.

  const int timer_hetaeta = (Options::Timer) ? timer_new("Remapeta_hetaeta") : 0;

  // clang-format off
                   cdo_operator_add("remapeta",   0, 0, "VCT file name");
  auto REMAPETAS = cdo_operator_add("remapeta_s", 0, 0, "VCT file name");
  auto REMAPETAZ = cdo_operator_add("remapeta_z", 0, 0, "VCT file name");
  // clang-format on

  auto operatorID = cdo_operator_id();

  operator_input_arg(cdo_operator_enter(operatorID));

  auto envstr = getenv("REMAPETA_PTOP");
  if (envstr)
    {
      auto fval = atof(envstr);
      if (fval > 0.0)
        {
          cptop = fval;
          cdo_print("Set REMAPETA_PTOP to %g", cptop);
        }
    }

  auto vct2 = vctFromFile(cdo_operator_argv(0).c_str());
  const int nvct2 = vct2.size();
  auto numFullLevels2 = nvct2 / 2 - 1;

  auto a2 = vct2.data();
  auto b2 = vct2.data() + nvct2 / 2;

  if (Options::cdoVerbose)
    for (int i = 0; i < numFullLevels2 + 1; ++i) cdo_print("vct2: %5d %25.17f %25.17f", i, vct2[i], vct2[nvct2 / 2 + i]);

  auto streamID1 = cdo_open_read(0);

  if (cdo_operator_argc() == 2)
    {
      lfis2 = true;

      const char *fname = cdo_operator_argv(1).c_str();
      auto streamID = stream_open_read_locked(fname);
      auto vlistID1 = streamInqVlist(streamID);

      int varID, levelID;
      streamInqRecord(streamID, &varID, &levelID);
      auto gridID = vlistInqVarGrid(vlistID1, varID);
      nfis2gp = gridInqSize(gridID);

      fis2.resize(nfis2gp);

      size_t nmiss;
      streamReadRecord(streamID, fis2.data(), &nmiss);

      if (nmiss)
        {
          missval = vlistInqVarMissval(vlistID1, varID);
          imiss.resize(nfis2gp);
          for (size_t i = 0; i < nfis2gp; ++i) imiss[i] = dbl_is_equal(fis2[i], missval);

          nmissout = nmiss;
        }

      // check range of surface_geopotential
      auto mm = array_min_max_mask(nfis2gp, fis2.data(), imiss);
      if (mm.min < MIN_FIS || mm.max > MAX_FIS)
        cdo_warning("%s out of range (min=%g max=%g)!", var_stdname(surface_geopotential), mm.min, mm.max);

      if (mm.min < -1.e10 || mm.max > 1.e10) cdo_abort("%s out of range!", var_stdname(surface_geopotential));

      streamClose(streamID);
    }

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto vlistID2 = vlistDuplicate(vlistID1);

  auto taxisID1 = vlistInqTaxis(vlistID1);
  auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  VarList varList1;
  varListInit(varList1, vlistID1);

  auto gridID0 = vlistGrid(vlistID1, 0);
  if (gridInqType(gridID0) == GRID_SPECTRAL) cdo_abort("Spectral data unsupported!");

  auto gridsize = vlist_check_gridsize(vlistID1);

  auto zaxisID2 = zaxisCreate(ZAXIS_HYBRID, numFullLevels2);

  {
    Varray<double> lev2(numFullLevels2);
    for (int i = 0; i < numFullLevels2; ++i) lev2[i] = i + 1;
    zaxisDefLevels(zaxisID2, lev2.data());
  }

  if (nvct2 == 0) cdo_abort("Internal problem, vct2 undefined!");
  zaxisDefVct(zaxisID2, nvct2, vct2.data());

  auto surfaceID = zaxis_from_name("surface");

  int zaxisID_ML = -1;
  int numHybridLevels = 0, numFullLevels1 = 0, numHalfLevels1 = 0;
  auto vct1 = vlist_read_vct(vlistID1, zaxisID_ML, numHybridLevels, numFullLevels1, numHalfLevels1);
  const int nvct1 = vct1.size();

  vlist_change_hybrid_zaxis(vlistID1, vlistID2, zaxisID_ML, zaxisID2);

  auto nzaxis = vlistNzaxis(vlistID1);
  for (int i = 0; i < nzaxis; ++i)
    {
      auto zaxisID = vlistZaxis(vlistID1, i);
      auto nlevels = zaxisInqSize(zaxisID);
      if (zaxisInqType(zaxisID) == ZAXIS_HYBRID && nlevels == 1) vlistChangeZaxisIndex(vlistID2, i, surfaceID);
    }

  const auto *a1 = vct1.data();
  const auto *b1 = vct1.data() + nvct1 / 2;
  if (Options::cdoVerbose)
    for (int i = 0; i < nvct1 / 2; ++i) cdo_print("vct1: %5d %25.17f %25.17f", i, vct1[i], vct1[nvct1 / 2 + i]);

  auto streamID2 = cdo_open_write(1);
  cdo_def_vlist(streamID2, vlistID2);

  VarList varList2;
  varListInit(varList2, vlistID2);

  if (zaxisID_ML == -1) cdo_warning("No 3D variable with hybrid sigma pressure coordinate found!");

  auto nvars = vlistNvars(vlistID1);

  VarIDs varIDs = search_varIDs(varList1, vlistID1, numFullLevels1);

  if (Options::cdoVerbose)
    {
      cdo_print("Found:");
      // clang-format off
      if (-1 != varIDs.tempID)    cdo_print("  %s", var_stdname(air_temperature));
      if (-1 != varIDs.psID)      cdo_print("  %s", var_stdname(surface_air_pressure));
      if (-1 != varIDs.lnpsID)    cdo_print("  LOG(%s)", var_stdname(surface_air_pressure));
      if (-1 != varIDs.sgeopotID) cdo_print("  %s", var_stdname(surface_geopotential));
      if (-1 != varIDs.humID)     cdo_print("  %s", var_stdname(specific_humidity));
      // clang-format on
    }

  for (int varID = 0; varID < nvars; ++varID)
    {
      auto gridID = varList1[varID].gridID;
      auto zaxisID = varList1[varID].zaxisID;
      auto nlevels = varList1[varID].nlevels;

      if (gridInqType(gridID) == GRID_SPECTRAL && zaxisInqType(zaxisID) == ZAXIS_HYBRID)
        cdo_abort("Spectral data on model level unsupported!");

      if (gridInqType(gridID) == GRID_SPECTRAL) cdo_abort("Spectral data unsupported!");

      if (zaxisInqType(zaxisID) == ZAXIS_HYBRID && zaxisID_ML != -1 && nlevels == numFullLevels1)
        {
          if (!(varID == varIDs.tempID || varID == varIDs.humID)) varids[nvars3D++] = varID;
        }
      else
        {
          if (varID == varIDs.tempID) varIDs.tempID = -1;
          if (varID == varIDs.humID) varIDs.humID = -1;
        }
    }

  auto ltq = (varIDs.tempID != -1 && varIDs.humID != -1);

  if (!ltq)
    {
      if (varIDs.tempID != -1) cdo_abort("Temperature without humidity unsupported!");
      if (varIDs.humID != -1) cdo_abort("Humidity without temperature unsupported!");
    }

  Varray<double> sum1, sum2;
  if (operatorID == REMAPETAS || operatorID == REMAPETAZ)
    {
      sum1.resize(gridsize);
      sum2.resize(gridsize);
    }

  Varray<double> deltap1, deltap2;
  Varray<double> halfPress1, halfPress2;
  if (operatorID == REMAPETAZ)
    {
      deltap1.resize(gridsize * numFullLevels1);
      deltap2.resize(gridsize * numFullLevels2);
      halfPress1.resize(gridsize * (numFullLevels1 + 1));
      halfPress2.resize(gridsize * (numFullLevels2 + 1));
    }

  Field field;
  if (memType == MemType::Float)
    field.resizef(gridsize);
  else
    field.resize(gridsize);

  Varray<double> fis1(gridsize);
  Varray<double> ps1(gridsize);

  if (!lfis2) fis2.resize(gridsize);
  if (lfis2 && gridsize != nfis2gp) cdo_abort("Orographies have different grid size!");

  Varray<double> ps2(gridsize);

  Varray<T> t1, t2;
  Varray<T> q1, q2;
  Varray<double> tscor, pscor, secor;
  if (ltq)
    {
      tscor.resize(gridsize);
      pscor.resize(gridsize);
      secor.resize(gridsize);

      t1.resize(gridsize * numFullLevels1);
      q1.resize(gridsize * numFullLevels1);

      t2.resize(gridsize * numFullLevels2);
      q2.resize(gridsize * numFullLevels2);
    }

  Varray2D<T> vars1, vars2;
  if (nvars3D)
    {
      vars1.resize(nvars);
      vars2.resize(nvars);
      for (int varID = 0; varID < nvars3D; ++varID) vars1[varID].resize(gridsize * numFullLevels1);
      for (int varID = 0; varID < nvars3D; ++varID) vars2[varID].resize(gridsize * numFullLevels2);
    }

  if (zaxisID_ML != -1 && varIDs.sgeopotID == -1)
    {
      varray_fill(fis1, 0.0);
      if (ltq) cdo_warning("%s not found - set to zero!", var_stdname(surface_geopotential));
    }

  int presID = varIDs.lnpsID;
  if (zaxisID_ML != -1 && varIDs.lnpsID == -1)
    {
      if (varIDs.psID == -1)
        cdo_abort("%s not found!", var_stdname(surface_air_pressure));
      else
        presID = varIDs.psID;
    }

  if (Options::cdoVerbose)
    {
      if (presID == varIDs.lnpsID)
        cdo_print("using LOG(%s)", var_stdname(surface_air_pressure));
      else
        cdo_print("using %s", var_stdname(surface_air_pressure));
    }

  if (Options::cdoVerbose) cdo_print("nvars3D = %d   ltq = %d", nvars3D, (int) ltq);

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
          cdo_read_record(streamID1, field);

          if (zaxisID_ML != -1)
            {
              auto zaxisID = varList1[varID].zaxisID;
              auto nlevels = varList1[varID].nlevels;
              auto offset = gridsize * levelID;

              if (varID == varIDs.sgeopotID)
                field_copy_array(gridsize, field, fis1.data());
              else if (varID == presID)
                {
                  if (varIDs.lnpsID != -1)
                    {
                      if (field.memType == MemType::Float)
                        for (size_t i = 0; i < gridsize; ++i) ps1[i] = std::exp((double) field.vec_f[i]);
                      else
                        for (size_t i = 0; i < gridsize; ++i) ps1[i] = std::exp(field.vec_d[i]);
                    }
                  else if (varIDs.psID != -1)
                    field_copy_array(gridsize, field, ps1.data());
                }
              else if (ltq && varID == varIDs.tempID)
                field_copy_array(gridsize, field, &t1[offset]);
              else if (ltq && varID == varIDs.humID)
                field_copy_array(gridsize, field, &q1[offset]);
              // else if ( zaxisID == zaxisID_ML )
              else if (zaxisInqType(zaxisID) == ZAXIS_HYBRID && nlevels == numFullLevels1)
                {
                  int i;
                  for (i = 0; i < nvars3D; ++i)
                    if (varID == varids[i]) break;

                  if (i == nvars3D) cdo_abort("Internal error, 3D variable not found!");

                  field_copy_array(gridsize, field, &vars1[i][offset]);
                }
              else
                {
                  cdo_def_record(streamID2, varID, levelID);
                  cdo_write_record(streamID2, field);
                }
            }
          else
            {
              cdo_def_record(streamID2, varID, levelID);
              cdo_write_record(streamID2, field);
            }
        }

      if (zaxisID_ML != -1)
        {
          // check range of psProg
          auto mm = array_min_max_mask(gridsize, ps1.data(), imiss);
          if (mm.min < MIN_PS || mm.max > MAX_PS) cdo_warning("Surface pressure out of range (min=%g max=%g)!", mm.min, mm.max);

          // check range of geop
          mm = array_min_max_mask(gridsize, fis1.data(), imiss);
          if (mm.min < MIN_FIS || mm.max > MAX_FIS) cdo_warning("Orography out of range (min=%g max=%g)!", mm.min, mm.max);
        }

      if (!lfis2)
        for (size_t i = 0; i < gridsize; ++i) fis2[i] = fis1[i];

      if (ltq)
        {
          int varID = varIDs.tempID;
          for (int levelID = 0; levelID < varList1[varID].nlevels; ++levelID)
            {
              auto offset = gridsize * levelID;
              auto mm = array_min_max_mask(gridsize, &t1[offset], imiss);
              if (mm.min < MIN_T || mm.max > MAX_T)
                cdo_warning("Input temperature at level %d out of range (min=%g max=%g)!", levelID + 1, mm.min, mm.max);
            }

          varID = varIDs.humID;
          for (int levelID = 0; levelID < varList1[varID].nlevels; ++levelID)
            {
              auto offset = gridsize * levelID;
              corr_hum(gridsize, &q1[offset], MIN_Q);

              auto mm = array_min_max_mask(gridsize, &q1[offset], imiss);
              if (mm.min < MIN_Q || mm.max > MAX_Q)
                cdo_warning("Input humidity at level %d out of range (min=%g max=%g)!", levelID + 1, mm.min, mm.max);
            }
        }

      if (nvars3D || ltq)
        {
          if (Options::Timer) timer_start(timer_hetaeta);
          hetaeta(ltq, gridsize, imiss.data(), numFullLevels1, a1, b1, fis1, ps1, t1, q1, numFullLevels2, a2, b2, fis2, ps2, t2, q2,
                  nvars3D, vars1, vars2, tscor, pscor, secor);
          if (Options::Timer) timer_stop(timer_hetaeta);
        }

      const long nctop = (cptop > 0) ? ncctop(cptop, (long) numFullLevels2, (long) numFullLevels2 + 1, a2, b2) : 0;

      if (zaxisID_ML != -1 && varIDs.sgeopotID != -1)
        {
          int varID = varIDs.sgeopotID;
          int levelID = 0;
          setmissval(gridsize, imiss, missval, fis2.data());
          cdo_def_record(streamID2, varID, levelID);
          cdo_write_record(streamID2, fis2.data(), nmissout);
        }

      if (zaxisID_ML != -1 && varIDs.lnpsID != -1)
        for (size_t i = 0; i < gridsize; ++i) ps2[i] = std::log(ps2[i]);

      if (zaxisID_ML != -1 && presID != -1)
        {
          int varID = presID;
          int levelID = 0;
          setmissval(gridsize, imiss, missval, ps2.data());
          cdo_def_record(streamID2, varID, levelID);
          cdo_write_record(streamID2, ps2.data(), nmissout);
        }

      if (ltq)
        {
          int varID = varIDs.tempID;
          for (int levelID = 0; levelID < varList2[varID].nlevels; ++levelID)
            {
              auto offset = gridsize * levelID;
              auto single2 = &t2[offset];

              auto mm = array_min_max_mask(gridsize, single2, imiss);
              if (mm.min < MIN_T || mm.max > MAX_T)
                cdo_warning("Output temperature at level %d out of range (min=%g max=%g)!", levelID + 1, mm.min, mm.max);

              setmissval(gridsize, imiss, missval, single2);
              cdo_def_record(streamID2, varID, levelID);

              if (memType == MemType::Float)
                cdo_write_record_f(streamID2, (float *) single2, nmissout);
              else
                cdo_write_record(streamID2, (double *) single2, nmissout);
            }

          varID = varIDs.humID;
          for (int levelID = 0; levelID < varList2[varID].nlevels; ++levelID)
            {
              auto offset = gridsize * levelID;
              auto single2 = &q2[offset];

              corr_hum(gridsize, single2, MIN_Q);

              if (levelID < nctop)
                for (size_t i = 0; i < gridsize; ++i) single2[i] = cconst;

              auto mm = array_min_max_mask(gridsize, single2, imiss);
              if (mm.min < MIN_Q || mm.max > MAX_Q)
                cdo_warning("Output humidity at level %d out of range (min=%g max=%g)!", levelID + 1, mm.min, mm.max);

              setmissval(gridsize, imiss, missval, single2);
              cdo_def_record(streamID2, varID, levelID);

              if (memType == MemType::Float)
                cdo_write_record_f(streamID2, (float *) single2, nmissout);
              else
                cdo_write_record(streamID2, (double *) single2, nmissout);
            }
        }

      for (int iv = 0; iv < nvars3D; ++iv)
        {
          int varID = varids[iv];

          auto nlevels = varList2[varID].nlevels;

          if (operatorID == REMAPETAS)
            {
              vertSum(sum1, vars1[iv], gridsize, numFullLevels1);
              vertSum(sum2, vars2[iv], gridsize, numFullLevels2);
            }
          else if (operatorID == REMAPETAZ)
            {
              vct_to_hybrid_pressure((double *) nullptr, halfPress1.data(), vct1.data(), ps1.data(), numFullLevels1, gridsize);
              for (int k = 0; k < numFullLevels1; ++k)
                for (size_t i = 0; i < gridsize; ++i)
                  {
                    deltap1[k * gridsize + i] = halfPress1[(k + 1) * gridsize + i] - halfPress1[k * gridsize + i];
                    deltap1[k * gridsize + i] = std::log(deltap1[k * gridsize + i]);
                  }
              vertSumw(sum1, vars1[iv], gridsize, numFullLevels1, deltap1);

              vct_to_hybrid_pressure((double *) nullptr, halfPress2.data(), vct2.data(), ps1.data(), numFullLevels2, gridsize);
              for (int k = 0; k < numFullLevels2; ++k)
                for (size_t i = 0; i < gridsize; ++i)
                  {
                    deltap2[k * gridsize + i] = halfPress2[(k + 1) * gridsize + i] - halfPress2[k * gridsize + i];
                    deltap2[k * gridsize + i] = std::log(deltap2[k * gridsize + i]);
                  }
              vertSumw(sum2, vars2[iv], gridsize, numFullLevels2, deltap2);
            }

          for (int levelID = 0; levelID < nlevels; ++levelID)
            {
              auto offset = gridsize * levelID;
              auto single2 = &vars2[iv][offset];

              if (operatorID == REMAPETAS || operatorID == REMAPETAZ)
                for (size_t i = 0; i < gridsize; ++i) single2[i] = single2[i] * sum1[i] / sum2[i];

              setmissval(gridsize, imiss, missval, single2);
              cdo_def_record(streamID2, varID, levelID);

              if (memType == MemType::Float)
                cdo_write_record_f(streamID2, (float *) single2, nmissout);
              else
                cdo_write_record(streamID2, (double *) single2, nmissout);
            }
        }

      tsID++;
    }

  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);
}

void *
Remapeta(void *process)
{
  cdo_initialize(process);

  auto memType = Options::CDO_Memtype;
  if (memType == MemType::Float)
    remapeta<float>(memType);
  else
    remapeta<double>(memType);

  cdo_finish();

  return nullptr;
}
