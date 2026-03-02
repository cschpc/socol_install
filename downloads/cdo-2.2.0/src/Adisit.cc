/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Adisit      adisit          compute insitu from potential temperature
      Adisit      adipot          compute potential from insitu temperature
*/

#include <cdi.h>

#include "cdo_options.h"
#include "cdo_zaxis.h"
#include "process_int.h"
#include "cdo_vlist.h"
#include "param_conversion.h"
#include "util_string.h"
#include "field_functions.h"

/*
  Transformation from potential to in situ temperature according to Bryden, 1973,
  "New polynomials for thermal expansion, adiabatic temperature gradient and potential temperature of sea water".
  Deep Sea Research and Oceanographic Abstracts. 20, 401-408 (GILL P.602), which gives the inverse
  transformation for an approximate value, all terms linear in t are taken after that one newton step.
  For the check value 8.4678516 the accuracy is 0.2 mikrokelvin.
*/

// compute insitu temperature from potential temperature
static inline double
adisit(double tpot, double sal, double p)
{
  constexpr double a_a1 = 3.6504E-4, a_a2 = 8.3198E-5, a_a3 = 5.4065E-7, a_a4 = 4.0274E-9, a_b1 = 1.7439E-5, a_b2 = 2.9778E-7,
                   a_c1 = 8.9309E-7, a_c2 = 3.1628E-8, a_c3 = 2.1987E-10, a_d = 4.1057E-9, a_e1 = 1.6056E-10, a_e2 = 5.0484E-12;

  double qc = p * (a_a1 + p * (a_c1 - a_e1 * p));
  double qv = p * (a_b1 - a_d * p);
  double dc = 1.0 + p * (-a_a2 + p * (a_c2 - a_e2 * p));
  double dv = a_b2 * p;
  double qnq = -p * (-a_a3 + p * a_c3);
  double qn3 = -p * a_a4;

  double tpo = tpot;
  double qvs = qv * (sal - 35.0) + qc;
  double dvs = dv * (sal - 35.0) + dc;
  double t = (tpo + qvs) / dvs;
  double fne = -qvs + t * (dvs + t * (qnq + t * qn3)) - tpo;
  double fst = dvs + t * (2.0 * qnq + 3.0 * qn3 * t);
  t = t - fne / fst;

  return t;
}

// compute potential temperature from insitu temperature
// Ref: Gill, p. 602, Section A3.5:Potential Temperature
static inline double
adipot(double t, double s, double p)
{
  constexpr double a_a1 = 3.6504E-4, a_a2 = 8.3198E-5, a_a3 = 5.4065E-7, a_a4 = 4.0274E-9, a_b1 = 1.7439E-5, a_b2 = 2.9778E-7,
                   a_c1 = 8.9309E-7, a_c2 = 3.1628E-8, a_c3 = 2.1987E-10, a_d = 4.1057E-9, a_e1 = 1.6056E-10, a_e2 = 5.0484E-12;

  double s_rel = s - 35.0;

  double aa = (a_a1 + t * (a_a2 - t * (a_a3 - a_a4 * t)));
  double bb = s_rel * (a_b1 - a_b2 * t);
  double cc = (a_c1 + t * (-a_c2 + a_c3 * t));
  double cc1 = a_d * s_rel;
  double dd = (-a_e1 + a_e2 * t);

  double tpot = t - p * (aa + bb + p * (cc - cc1 + p * dd));

  return tpot;
}

static void
calc_adisit(size_t gridsize, size_t nlevel, const Varray<double> &pressure, const FieldVector &tho, const FieldVector &sao,
            FieldVector &tis)
{
  // pressure units: hPa
  // tho units:      Celsius
  // sao units:      psu

  for (size_t levelID = 0; levelID < nlevel; ++levelID)
    {
      const auto &thovec = tho[levelID].vec_d;
      const auto &saovec = sao[levelID].vec_d;
      auto &tisvec = tis[levelID].vec_d;
      auto thoMissval = tho[levelID].missval;
      auto saoMissval = sao[levelID].missval;
      auto tisMissval = tis[levelID].missval;
      for (size_t i = 0; i < gridsize; ++i)
        {
          auto isMissing = (DBL_IS_EQUAL(thovec[i], thoMissval) || DBL_IS_EQUAL(saovec[i], saoMissval));
          tisvec[i] = isMissing ? tisMissval : adisit(thovec[i], saovec[i], pressure[levelID]);
        }
    }
}

static void
calc_adipot(size_t gridsize, size_t nlevel, const Varray<double> &pressure, const FieldVector &t, const FieldVector &s,
            FieldVector &tpot)
{
  // pressure units: hPa
  // t units:        Celsius
  // s units:        psu

  for (size_t levelID = 0; levelID < nlevel; ++levelID)
    {
      const auto &tvec = t[levelID].vec_d;
      const auto &svec = s[levelID].vec_d;
      auto &tpotvec = tpot[levelID].vec_d;
      auto tMissval = t[levelID].missval;
      auto sMissval = s[levelID].missval;
      auto tpotMissval = tpot[levelID].missval;
      for (size_t i = 0; i < gridsize; ++i)
        {
          auto isMissing = (DBL_IS_EQUAL(tvec[i], tMissval) || DBL_IS_EQUAL(svec[i], sMissval));
          tpotvec[i] = isMissing ? tpotMissval : adipot(tvec[i], svec[i], pressure[levelID]);
        }
    }
}

int
get_code(int vlistID1, int varID, const std::string &cname)
{
  auto code = vlistInqVarCode(vlistID1, varID);
  if (code <= 0)
    {
      auto varname = string_to_lower(cdo::inq_var_name(vlistID1, varID));
      auto stdname = string_to_lower(cdo::inq_key_string(vlistID1, varID, CDI_KEY_STDNAME));

      if (varname == "s" || stdname == "sea_water_salinity") { code = 5; }
      else if (varname == "t")
        {
          code = 2;
        }

      if (stdname == cname) code = 2;
    }

  return code;
}

struct IOSettings
{
  CdoStreamID streamID2;
  int vlistID2;
  size_t gridsize;
  int nlevels;
  int taxisID1;
  int taxisID2;
  int tisID2;
  int saoID2;
};

IOSettings
configureOutput(const std::function<void(const int, const int)> &outputSettingFunc, int vlistID1, int thoID, int saoID,
                FieldVector &tho, FieldVector &sao, FieldVector &tis, Varray<double> &pressure)
{
  double pin = (cdo_operator_argc() == 1) ? parameter_to_double(cdo_operator_argv(0)) : -1.0;

  VarList varList1;
  varListInit(varList1, vlistID1);

  auto units = varList1[thoID].units;
  if (units.empty()) units = "Celcius";

  auto gridID = vlistGrid(vlistID1, 0);
  auto gridsize = vlist_check_gridsize(vlistID1);

  auto nlevels1 = varList1[saoID].nlevels;
  auto nlevels2 = varList1[thoID].nlevels;
  auto zaxisID = varList1[thoID].zaxisID;

  if (nlevels1 != nlevels2) cdo_abort("temperature and salinity have different number of levels!");
  auto nlevels = nlevels1;

  pressure.resize(nlevels);
  cdo_zaxis_inq_levels(zaxisID, pressure.data());

  if (pin >= 0)
    for (int i = 0; i < nlevels; ++i) pressure[i] = pin;
  else
    for (int i = 0; i < nlevels; ++i) pressure[i] /= 10;

  if (Options::cdoVerbose)
    {
      cdo_print("Level Pressure");
      for (int i = 0; i < nlevels; ++i) cdo_print("%5d  %g", i + 1, pressure[i]);
    }

  tho.resize(nlevels, Field{});
  sao.resize(nlevels, Field{});
  tis.resize(nlevels, Field{});
  for (int levelID = 0; levelID < nlevels; ++levelID)
    {
      tho[levelID].resize(gridsize);
      sao[levelID].resize(gridsize);
      tis[levelID].resize(gridsize);
      tho[levelID].missval = varList1[thoID].missval;
      sao[levelID].missval = varList1[saoID].missval;
      tis[levelID].missval = tho[levelID].missval;
    }

  int datatype = CDI_DATATYPE_FLT32;
  if (varList1[thoID].datatype == CDI_DATATYPE_FLT64 && varList1[saoID].datatype == CDI_DATATYPE_FLT64)
    datatype = CDI_DATATYPE_FLT64;

  auto vlistID2 = vlistCreate();
  vlistDefNtsteps(vlistID2, vlistNtsteps(vlistID1));

  auto tisID2 = vlistDefVar(vlistID2, gridID, zaxisID, TIME_VARYING);

  outputSettingFunc(vlistID2, tisID2);
  cdiDefKeyString(vlistID2, tisID2, CDI_KEY_UNITS, units.c_str());
  vlistDefVarMissval(vlistID2, tisID2, tis[0].missval);
  vlistDefVarDatatype(vlistID2, tisID2, datatype);

  auto saoID2 = vlistDefVar(vlistID2, gridID, zaxisID, TIME_VARYING);
  vlistDefVarParam(vlistID2, saoID2, cdiEncodeParam(5, 255, 255));
  cdiDefKeyString(vlistID2, saoID2, CDI_KEY_NAME, "s");
  cdiDefKeyString(vlistID2, saoID2, CDI_KEY_LONGNAME, "Sea water salinity");
  cdiDefKeyString(vlistID2, saoID2, CDI_KEY_STDNAME, "sea_water_salinity");
  cdiDefKeyString(vlistID2, saoID2, CDI_KEY_UNITS, "psu");
  vlistDefVarMissval(vlistID2, saoID2, sao[0].missval);
  vlistDefVarDatatype(vlistID2, saoID2, datatype);

  auto taxisID1 = vlistInqTaxis(vlistID1);
  auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  auto streamID2 = cdo_open_write(1);
  cdo_def_vlist(streamID2, vlistID2);

  return IOSettings{ streamID2, vlistID2, gridsize, nlevels, taxisID1, taxisID2, tisID2, saoID2 };
}

void *
Adisit(void *process)
{
  int thoID = -1, saoID = -1;

  cdo_initialize(process);

  auto ADISIT = cdo_operator_add("adisit", 0, 0, "");
  auto ADIPOT = cdo_operator_add("adipot", 0, 0, "");

  auto operatorID = cdo_operator_id();

  auto streamID1 = cdo_open_read(0);
  auto vlistID1 = cdo_stream_inq_vlist(streamID1);

  auto nvars = vlistNvars(vlistID1);

  std::string cname_ADISIT{ "sea_water_potential_temperature" };
  std::string cname_ADIPOT{ "sea_water_temperature" };
  std::string cname = std::string{ (operatorID == ADISIT) ? cname_ADISIT : cname_ADIPOT };

  for (int varID = 0; varID < nvars; ++varID)
    {
      auto code = get_code(vlistID1, varID, cname);

      if (code == 2) { thoID = varID; }
      else if (code == 20 && operatorID == ADIPOT)
        {
          thoID = varID;
        }
      else if (code == 5)
        {
          saoID = varID;
        }
    }

  if (saoID == -1) cdo_abort("Sea water salinity not found!");
  if (thoID == -1) cdo_abort("%s temperature not found!", (operatorID == ADISIT) ? "Potential" : "Insitu");

  auto outputSetting_ADISIT = [](int vlistID2, int tisID2) {
    vlistDefVarParam(vlistID2, tisID2, cdiEncodeParam(20, 255, 255));
    cdiDefKeyString(vlistID2, tisID2, CDI_KEY_NAME, "to");
    cdiDefKeyString(vlistID2, tisID2, CDI_KEY_LONGNAME, "Sea water temperature");
    cdiDefKeyString(vlistID2, tisID2, CDI_KEY_STDNAME, "sea_water_temperature");
  };

  auto outputSetting_ADIPOT = [](int vlistID2, int tisID2) {
    vlistDefVarParam(vlistID2, tisID2, cdiEncodeParam(2, 255, 255));
    cdiDefKeyString(vlistID2, tisID2, CDI_KEY_NAME, "tho");
    cdiDefKeyString(vlistID2, tisID2, CDI_KEY_LONGNAME, "Sea water potential temperature");
    cdiDefKeyString(vlistID2, tisID2, CDI_KEY_STDNAME, "sea_water_potential_temperature");
  };

  auto outputSettingFunc = (operatorID == ADISIT) ? outputSetting_ADISIT : outputSetting_ADIPOT;

  FieldVector tho, sao, tis;
  Varray<double> pressure;
  const auto &configResults = configureOutput(outputSettingFunc, vlistID1, thoID, saoID, tho, sao, tis, pressure);

  auto streamID2 = configResults.streamID2;
  auto vlistID2 = configResults.vlistID2;
  auto gridsize = configResults.gridsize;
  auto nlevels = configResults.nlevels;
  auto taxisID1 = configResults.taxisID1;
  auto taxisID2 = configResults.taxisID2;
  auto tisID2 = configResults.tisID2;
  auto saoID2 = configResults.saoID2;

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
          if (varID == thoID) cdo_read_record(streamID1, tho[levelID].vec_d.data(), &tho[levelID].nmiss);
          if (varID == saoID) cdo_read_record(streamID1, sao[levelID].vec_d.data(), &sao[levelID].nmiss);

          if (varID == thoID)
            {
              constexpr double MIN_T = -10.0;
              constexpr double MAX_T = 40.0;
              auto mm = field_min_max(tho[levelID]);
              if (mm.min < MIN_T || mm.max > MAX_T)
                cdo_warning("Temperature in degree Celsius out of range (min=%g max=%g) [timestep:%d levelIndex:%d]!", mm.min,
                            mm.max, tsID + 1, levelID + 1);
            }
        }

      if (operatorID == ADISIT)
        calc_adisit(gridsize, nlevels, pressure, tho, sao, tis);
      else
        calc_adipot(gridsize, nlevels, pressure, tho, sao, tis);

      for (int levelID = 0; levelID < nlevels; ++levelID)
        {
          cdo_def_record(streamID2, tisID2, levelID);
          cdo_write_record(streamID2, tis[levelID].vec_d.data(), field_num_miss(tis[levelID]));
        }

      for (int levelID = 0; levelID < nlevels; ++levelID)
        {
          cdo_def_record(streamID2, saoID2, levelID);
          cdo_write_record(streamID2, sao[levelID].vec_d.data(), field_num_miss(sao[levelID]));
        }

      tsID++;
    }

  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  vlistDestroy(vlistID2);

  cdo_finish();

  return nullptr;
}
