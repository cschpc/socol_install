/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Wind       uv2dv           U and V wind to divergence and vorticity
      Wind       dv2uv           Divergence and vorticity to U and V wind
      Wind       dv2ps           Divergence and vorticity to velocity potential and stream function
*/

#include <cdi.h>

#include "cdo_vlist.h"
#include "process_int.h"
#include "param_conversion.h"
#include <mpim_grid.h>
#include "griddes.h"
#include "specspace.h"
#include "util_string.h"

static void
defineAttributesDV(int vlistID2, int gridID2, int varID1, int varID2)
{
  vlistChangeVarGrid(vlistID2, varID1, gridID2);
  vlistChangeVarGrid(vlistID2, varID2, gridID2);
  vlistDefVarParam(vlistID2, varID1, cdiEncodeParam(155, 128, 255));
  vlistDefVarParam(vlistID2, varID2, cdiEncodeParam(138, 128, 255));
  cdiDefKeyString(vlistID2, varID1, CDI_KEY_NAME, "sd");
  cdiDefKeyString(vlistID2, varID2, CDI_KEY_NAME, "svo");
  cdiDefKeyString(vlistID2, varID1, CDI_KEY_LONGNAME, "divergence");
  cdiDefKeyString(vlistID2, varID2, CDI_KEY_LONGNAME, "vorticity");
  cdiDefKeyString(vlistID2, varID1, CDI_KEY_UNITS, "1/s");
  cdiDefKeyString(vlistID2, varID2, CDI_KEY_UNITS, "1/s");
}

static void
defineAttributesUV(int vlistID2, int gridID2, int varID1, int varID2)
{
  vlistChangeVarGrid(vlistID2, varID1, gridID2);
  vlistChangeVarGrid(vlistID2, varID2, gridID2);
  vlistDefVarParam(vlistID2, varID1, cdiEncodeParam(131, 128, 255));
  vlistDefVarParam(vlistID2, varID2, cdiEncodeParam(132, 128, 255));
  cdiDefKeyString(vlistID2, varID1, CDI_KEY_NAME, "u");
  cdiDefKeyString(vlistID2, varID2, CDI_KEY_NAME, "v");
  cdiDefKeyString(vlistID2, varID1, CDI_KEY_LONGNAME, "u-velocity");
  cdiDefKeyString(vlistID2, varID2, CDI_KEY_LONGNAME, "v-velocity");
  cdiDefKeyString(vlistID2, varID1, CDI_KEY_UNITS, "m/s");
  cdiDefKeyString(vlistID2, varID2, CDI_KEY_UNITS, "m/s");
}

static void
defineAttributesPS(int vlistID2, int varID1, int varID2)
{
  vlistDefVarParam(vlistID2, varID1, cdiEncodeParam(149, 128, 255));
  vlistDefVarParam(vlistID2, varID2, cdiEncodeParam(148, 128, 255));
  cdiDefKeyString(vlistID2, varID1, CDI_KEY_NAME, "velopot");
  cdiDefKeyString(vlistID2, varID2, CDI_KEY_NAME, "stream");
  cdiDefKeyString(vlistID2, varID1, CDI_KEY_LONGNAME, "velocity potential");
  cdiDefKeyString(vlistID2, varID2, CDI_KEY_LONGNAME, "streamfunction");
  cdiDefKeyString(vlistID2, varID1, CDI_KEY_UNITS, "m^2/s");
  cdiDefKeyString(vlistID2, varID2, CDI_KEY_UNITS, "m^2/s");
}

class ModuleWind
{
  CdoStreamID streamID1;
  CdoStreamID streamID2;

  int taxisID1;
  int taxisID2;

  Varray<double> array1;

  Varray<double> ivar1, ivar2, ovar1, ovar2;

  bool dataIsUnchanged = true;

  bool luv2dv;
  bool ldv2uv;
  bool linear;

  int operatorID;

  int UV2DV;
  int UV2DVL;
  int DV2UV;
  int DV2UVL;
  int DV2PS;

  size_t nlev = 0;
  int gridID1 = -1, gridID2 = -1;
  size_t nmiss;
  long ntr = -1;
  int varID1 = -1, varID2 = -1;
  SP_Transformation spTrans;
  DV_Transformation dvTrans;

public:
  void
  init(void *process)
  {
    cdo_initialize(process);

    dataIsUnchanged = data_is_unchanged();

    // clang-format off
    UV2DV  = cdo_operator_add("uv2dv",  0, 0, nullptr);
    UV2DVL = cdo_operator_add("uv2dvl", 0, 0, nullptr);
    DV2UV  = cdo_operator_add("dv2uv",  0, 0, nullptr);
    DV2UVL = cdo_operator_add("dv2uvl", 0, 0, nullptr);
    DV2PS  = cdo_operator_add("dv2ps",  0, 0, nullptr);

    operatorID = cdo_operator_id();

     luv2dv = (operatorID == UV2DV || operatorID == UV2DVL);
     ldv2uv = (operatorID == DV2UV || operatorID == DV2UVL);
     linear = (operatorID == UV2DVL || operatorID == DV2UVL);

    int (*nlat2ntr)(int) = linear ? nlat_to_ntr_linear : nlat_to_ntr;
    const char *ctype = linear ? "l" : "";

    if ((luv2dv || ldv2uv) && cdo_operator_argc() == 1)
      {
        std::string type = parameter_to_word(cdo_operator_argv(0));
        if      (type == "linear")    { nlat2ntr = nlat_to_ntr_linear; ctype = "l"; }
        else if (type == "cubic")     { nlat2ntr = nlat_to_ntr_cubic; ctype = "c"; }
        else if (type == "quadratic") { nlat2ntr = nlat_to_ntr; }
        else cdo_abort("Unsupported type: %s\n", type);
      }
    // clang-format on

    streamID1 = cdo_open_read(0);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = vlistDuplicate(vlistID1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    // find variables
    auto nvars = vlistNvars(vlistID2);
    for (int varID = 0; varID < nvars; ++varID)
      {
        auto varname = string_to_lower(cdo::inq_var_name(vlistID1, varID));
        auto param = vlistInqVarParam(vlistID2, varID);
        int pnum, pcat, pdis;
        cdiDecodeParam(param, &pnum, &pcat, &pdis);
        int code = pnum;
        if (operatorID == UV2DV || operatorID == UV2DVL)
          {
            // search for u and v wind
            if (pdis != 255 || code <= 0)
              {
                if (varname == "u") code = 131;
                if (varname == "v") code = 132;
              }

            if (code == 131) varID1 = varID;
            if (code == 132) varID2 = varID;
          }
        else if (operatorID == DV2UV || operatorID == DV2UVL || operatorID == DV2PS)
          {
            // search for divergence and vorticity
            if (pdis != 255)  // GRIB2
              {
                if (varname == "d") code = 155;
                if (varname == "vo") code = 138;
              }
            else if (code <= 0)
              {
                if (varname == "sd") code = 155;
                if (varname == "svo") code = 138;
              }

            if (code == 155) varID1 = varID;
            if (code == 138) varID2 = varID;
          }
        else
          cdo_abort("Unexpected operatorID %d", operatorID);
      }

    auto gridIDsp = vlist_get_first_spectral_grid(vlistID1);
    auto gridIDgp = vlist_get_first_gaussian_grid(vlistID1);

    // define output grid
    if (luv2dv)
      {
        if (varID1 == -1) cdo_warning("U-wind not found!");
        if (varID2 == -1) cdo_warning("V-wind not found!");

        if (varID1 != -1 && varID2 != -1)
          {
            gridID1 = vlistInqVarGrid(vlistID1, varID1);

            if (gridInqType(gridID1) != GRID_GAUSSIAN) cdo_abort("U-wind is not on Gaussian grid!");
            if (gridID1 != vlistInqVarGrid(vlistID1, varID2)) cdo_abort("U and V wind must have the same grid represention!");

            auto numLPE = gridInqNP(gridID1);
            const long nlon = gridInqXsize(gridID1);
            const long nlat = gridInqYsize(gridID1);
            const long ntr1 = nlat2ntr(nlat);

            if (numLPE > 0 && nlat != (numLPE * 2)) cdo_abort("U and V fields on Gaussian grid are not global!");

            if (gridIDsp != -1)
              if (ntr1 != gridInqTrunc(gridIDsp)) gridIDsp = -1;

            if (gridIDsp == -1)
              {
                gridIDsp = gridCreate(GRID_SPECTRAL, (ntr1 + 1) * (ntr1 + 2));
                gridDefTrunc(gridIDsp, ntr1);
                gridDefComplexPacking(gridIDsp, 1);
              }

            if (gridIDsp == -1) cdo_abort("No Gaussian grid data found!");

            gridID2 = gridIDsp;

            defineAttributesDV(vlistID2, gridID2, varID1, varID2);

            ntr = gridInqTrunc(gridID2);
            nlev = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID1));

            spTrans.init(nlon, nlat, ntr, PolFlag::UV2DV, nlev);
          }
      }
    else if (ldv2uv)
      {
        if (varID1 == -1) cdo_warning("Divergence not found!");
        if (varID2 == -1) cdo_warning("Vorticity not found!");

        if (varID1 != -1 && varID2 != -1)
          {
            gridID1 = vlistInqVarGrid(vlistID1, varID2);

            if (gridInqType(gridID1) != GRID_SPECTRAL) cdo_abort("Vorticity is not on spectral grid!");

            if (gridID1 != vlistInqVarGrid(vlistID1, varID1))
              cdo_abort("Divergence and vorticity must have the same grid represention!");

            if (gridIDgp != -1)
              {
                const long nlat = gridInqYsize(gridIDgp);
                const long ntr1 = nlat2ntr(nlat);
                if (gridInqTrunc(gridIDsp) != ntr1) gridIDgp = -1;
              }

            if (gridIDgp == -1)
              {
                char gridname[20];
                std::snprintf(gridname, sizeof(gridname), "t%s%dgrid", ctype, gridInqTrunc(gridIDsp));
                gridIDgp = grid_from_name(gridname);
              }

            gridID2 = gridIDgp;

            defineAttributesUV(vlistID2, gridID2, varID1, varID2);

            const long nlon = gridInqXsize(gridID2);
            const long nlat = gridInqYsize(gridID2);
            ntr = gridInqTrunc(gridID1);
            nlev = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID1));

            spTrans.init(nlon, nlat, ntr, PolFlag::SP2FC, nlev);
            dvTrans.init(ntr);
          }
      }
    else if (operatorID == DV2PS)
      {
        if (varID1 == -1) cdo_warning("Divergence not found!");
        if (varID2 == -1) cdo_warning("Vorticity not found!");

        if (varID1 != -1 && varID2 != -1)
          {
            gridID1 = vlistInqVarGrid(vlistID1, varID2);

            if (gridInqType(gridID1) != GRID_SPECTRAL) cdo_abort("Vorticity is not on spectral grid!");

            if (gridID1 != vlistInqVarGrid(vlistID1, varID1))
              cdo_abort("Divergence and vorticity must have the same grid represention!");

            defineAttributesPS(vlistID2, varID1, varID2);

            ntr = gridInqTrunc(gridID1);
            gridID2 = gridID1;
          }
      }

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);

    auto gridsizemax = vlistGridsizeMax(vlistID1);
    array1 = Varray<double>(gridsizemax);

    if (varID1 != -1 && varID2 != -1)
      {
        nlev = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID1));

        auto gridsize = gridInqSize(gridID1);
        ivar1.resize(nlev * gridsize);
        ivar2.resize(nlev * gridsize);

        gridsize = gridInqSize(gridID2);
        ovar1.resize(nlev * gridsize);
        ovar2.resize(nlev * gridsize);
      }
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
            int varID, levelID;
            cdo_inq_record(streamID1, &varID, &levelID);

            if ((varID1 != -1 && varID2 != -1) && (varID == varID1 || varID == varID2))
              {
                cdo_read_record(streamID1, array1.data(), &nmiss);
                if (nmiss) cdo_abort("Missing values unsupported for spectral data!");

                auto gridsize = gridInqSize(gridID1);
                auto offset = gridsize * levelID;

                if (varID == varID1)
                  array_copy(gridsize, array1.data(), &ivar1[offset]);
                else if (varID == varID2)
                  array_copy(gridsize, array1.data(), &ivar2[offset]);
              }
            else
              {
                cdo_def_record(streamID2, varID, levelID);
                if (dataIsUnchanged) { cdo_copy_record(streamID2, streamID1); }
                else
                  {
                    cdo_read_record(streamID1, array1.data(), &nmiss);
                    cdo_write_record(streamID2, array1.data(), nmiss);
                  }
              }
          }

        if (varID1 != -1 && varID2 != -1)
          {
            if (luv2dv)
              trans_uv2dv(spTrans, nlev, gridID1, ivar1.data(), ivar2.data(), gridID2, ovar1.data(), ovar2.data());
            else if (ldv2uv)
              trans_dv2uv(spTrans, dvTrans, nlev, gridID1, ivar1.data(), ivar2.data(), gridID2, ovar1.data(), ovar2.data());
            else if (operatorID == DV2PS)
              {
                dv2ps(ivar1.data(), ovar1.data(), nlev, ntr);
                dv2ps(ivar2.data(), ovar2.data(), nlev, ntr);
              }

            auto gridsize = gridInqSize(gridID2);
            if (luv2dv || operatorID == DV2PS)
              {
                for (size_t levelID = 0; levelID < nlev; ++levelID)
                  {
                    auto offset = gridsize * levelID;
                    cdo_def_record(streamID2, varID2, levelID);
                    cdo_write_record(streamID2, &ovar2[offset], 0);
                  }
                for (size_t levelID = 0; levelID < nlev; ++levelID)
                  {
                    auto offset = gridsize * levelID;
                    cdo_def_record(streamID2, varID1, levelID);
                    cdo_write_record(streamID2, &ovar1[offset], 0);
                  }
              }
            else if (ldv2uv)
              {
                for (size_t levelID = 0; levelID < nlev; ++levelID)
                  {
                    auto offset = gridsize * levelID;
                    cdo_def_record(streamID2, varID1, levelID);
                    cdo_write_record(streamID2, &ovar1[offset], 0);
                  }
                for (size_t levelID = 0; levelID < nlev; ++levelID)
                  {
                    auto offset = gridsize * levelID;
                    cdo_def_record(streamID2, varID2, levelID);
                    cdo_write_record(streamID2, &ovar2[offset], 0);
                  }
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

    cdo_finish();
  }
};
void *
Wind(void *process)
{
  ModuleWind wind;
  wind.init(process);
  wind.run();
  wind.close();

  return nullptr;
}
