/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Vertwind    vertwind      Convert the vertical velocity to [m/s]
*/

#include <cdi.h>

#include "process_int.h"
#include "cdo_vlist.h"
#include "cdi_lockedIO.h"
#include "vertical_interp.h"
#include "util_string.h"
#include "cdo_zaxis.h"

constexpr double R = 287.07;   // spezielle Gaskonstante fuer Luft
constexpr double G = 9.80665;  // Erdbeschleunigung

class ModuleVertwind
{
  CdoStreamID streamID1;
  CdoStreamID streamID2;

  int vlistID2;

  int taxisID1 = -1;
  int taxisID2 = -1;

  int zaxisID;

  size_t gridsize;
  int nlevels;

  int tempID = -1, sqID = -1, psID = -1, omegaID = -1;
  Varray<double> vct;
  Varray<double> hpress, psProg;
  Varray<double> temp;
  Varray<double> sq;
  Varray<double> omega;
  Varray<double> wms;
  Varray<double> fpress;

  double missval_t;
  double missval_sq;
  double missval_wap;
  double missval_out;

public:
  void
  init(void *process)
  {

    cdo_initialize(process);

    operator_check_argc(0);

    streamID1 = cdo_open_read(0);
    auto vlistID1 = cdo_stream_inq_vlist(streamID1);

    vlist_check_gridsize(vlistID1);

    VarList varList1;
    varListInit(varList1, vlistID1);

    int temp_code = 130;
    int sq_code = 133;
    int ps_code = 134;
    int omega_code = 135;

    int nvars = vlistNvars(vlistID1);
    for (int varID = 0; varID < nvars; ++varID)
      {
        int code = varList1[varID].code;

        if (code <= 0)
          {
            auto varname = string_to_lower(varList1[varID].name);

            if (varname == "st") { code = temp_code; }
            else if (varname == "sq")
              {
                code = sq_code;
              }
            else if (varname == "aps")
              {
                code = ps_code;
              }
            else if (varname == "omega")
              {
                code = omega_code;
              }
          }

        if (code == temp_code) { tempID = varID; }
        else if (code == sq_code)
          {
            sqID = varID;
          }
        else if (code == ps_code)
          {
            psID = varID;
          }
        else if (code == omega_code)
          {
            omegaID = varID;
          }
      }

    if (tempID == -1 || sqID == -1 || omegaID == -1)
      {
        if (tempID == -1) { cdo_warning("Temperature (code 130) not found!"); }
        if (sqID == -1) { cdo_warning("Specific humidity (code 133) not found!"); }
        if (omegaID == -1) { cdo_warning("Vertical velocity (code 135) not found!"); }

        cdo_abort("Parameter not found!");
      }

    // Get missing values
    missval_t = varList1[tempID].missval;
    missval_sq = varList1[sqID].missval;
    missval_wap = varList1[omegaID].missval;
    missval_out = missval_wap;

    auto gridID = varList1[omegaID].gridID;
    zaxisID = varList1[omegaID].zaxisID;

    if (psID == -1 && zaxisInqType(zaxisID) == ZAXIS_HYBRID) cdo_abort("Surface pressure (code 134) not found!");

    gridsize = gridInqSize(gridID);
    nlevels = zaxisInqSize(zaxisID);
    Varray<double> levels(nlevels);
    cdo_zaxis_inq_levels(zaxisID, levels.data());

    temp = Varray<double>(gridsize * nlevels);
    sq = Varray<double>(gridsize * nlevels);
    omega = Varray<double>(gridsize * nlevels);
    wms = Varray<double>(gridsize * nlevels);
    fpress = Varray<double>(gridsize * nlevels);

    if (zaxisInqType(zaxisID) == ZAXIS_PRESSURE)
      {
        for (int levelID = 0; levelID < nlevels; ++levelID)
          {
            auto offset = (size_t) levelID * gridsize;
            for (size_t i = 0; i < gridsize; ++i) fpress[offset + i] = levels[levelID];
          }
      }
    else if (zaxisInqType(zaxisID) == ZAXIS_HYBRID)
      {
        psProg.resize(gridsize);
        hpress.resize(gridsize * (nlevels + 1));

        auto nvct = zaxisInqVctSize(zaxisID);
        if (nlevels == (nvct / 2 - 1))
          {
            vct.resize(nvct);
            zaxisInqVct(zaxisID, vct.data());
          }
        else
          cdo_abort("Unsupported vertical coordinate table format!");
      }
    else
      cdo_abort("Unsupported Z-Axis type!");

    vlistClearFlag(vlistID1);
    for (int levelID = 0; levelID < nlevels; ++levelID) vlistDefFlag(vlistID1, omegaID, levelID, true);

    vlistID2 = vlistCreate();
    cdo_vlist_copy_flag(vlistID2, vlistID1);
    vlistDefNtsteps(vlistID2, vlistNtsteps(vlistID1));

    vlistDefVarCode(vlistID2, 0, 40);
    cdiDefKeyString(vlistID2, 0, CDI_KEY_NAME, "W");
    cdiDefKeyString(vlistID2, 0, CDI_KEY_LONGNAME, "Vertical velocity");
    cdiDefKeyString(vlistID2, 0, CDI_KEY_UNITS, "m/s");
    vlistDefVarMissval(vlistID2, 0, missval_out);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    streamID2 = cdo_open_write(1);

    cdo_def_vlist(streamID2, vlistID2);
  }

  void
  run()
  {
    int tsID = 0;
    while (true)
      {
        auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
        if (nrecs == 0) break;

        cdo_taxis_copy_timestep(taxisID2, taxisID1);
        cdo_def_timestep(streamID2, tsID);

        for (int recID = 0; recID < nrecs; ++recID)
          {
            size_t nmiss;
            int varID, levelID;
            cdo_inq_record(streamID1, &varID, &levelID);

            auto offset = (size_t) levelID * gridsize;

            if (varID == tempID)
              cdo_read_record(streamID1, &temp[offset], &nmiss);
            else if (varID == sqID)
              cdo_read_record(streamID1, &sq[offset], &nmiss);
            else if (varID == omegaID)
              cdo_read_record(streamID1, &omega[offset], &nmiss);
            else if (varID == psID && zaxisInqType(zaxisID) == ZAXIS_HYBRID)
              cdo_read_record(streamID1, psProg.data(), &nmiss);
          }

        if (zaxisInqType(zaxisID) == ZAXIS_HYBRID)
          vct_to_hybrid_pressure(fpress.data(), hpress.data(), vct.data(), psProg.data(), nlevels, gridsize);

        for (int levelID = 0; levelID < nlevels; ++levelID)
          {
            auto offset = (size_t) levelID * gridsize;

            for (size_t i = 0; i < gridsize; ++i)
              {
                if (dbl_is_equal(temp[offset + i], missval_t) || dbl_is_equal(omega[offset + i], missval_wap)
                    || dbl_is_equal(sq[offset + i], missval_sq))
                  {
                    wms[offset + i] = missval_out;
                  }
                else
                  {
                    // Virtuelle Temperatur bringt die Feuchteabhaengigkeit hinein
                    auto tv = temp[offset + i] * (1. + 0.608 * sq[offset + i]);

                    // Die Dichte erhaelt man nun mit der Gasgleichung rho=p/(R*tv) Level in Pa!
                    auto rho = fpress[offset + i] / (R * tv);
                    /*
                      Nun daraus die Vertikalgeschwindigkeit im m/s, indem man die Vertikalgeschwindigkeit
                      in Pa/s durch die Erdbeschleunigung und die Dichte teilt
                    */
                    wms[offset + i] = omega[offset + i] / (G * rho);
                  }
              }
          }

        for (int levelID = 0; levelID < nlevels; ++levelID)
          {
            auto offset = (size_t) levelID * gridsize;

            size_t nmiss_out = 0;
            for (size_t i = 0; i < gridsize; ++i)
              if (dbl_is_equal(wms[offset + i], missval_out)) nmiss_out++;

            cdo_def_record(streamID2, 0, levelID);
            cdo_write_record(streamID2, &wms[offset], nmiss_out);
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
Vertwind(void *process)
{
  ModuleVertwind vertwind;
  vertwind.init(process);
  vertwind.run();
  vertwind.close();

  return nullptr;
}
