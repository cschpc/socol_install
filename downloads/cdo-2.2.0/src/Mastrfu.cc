/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Mastrfu    mastrfu         Mass stream function
*/

#include <cdi.h>

#include "varray.h"
#include "process_int.h"
#include <mpim_grid.h>
#include "cdo_zaxis.h"

static void
mastrfu(int gridID, int zaxisID, const Varray2D<double> &field1, Varray2D<double> &field2, size_t nmiss, double missval)
{
  auto fact = 4.0 * std::atan(1.0) * 6371000.0 / 9.81;

  auto nlat = gridInqSize(gridID);
  auto nlev = zaxisInqSize(zaxisID);
  Varray<double> phi(nlat), cosphi(nlat), plevel(nlev);

  cdo_zaxis_inq_levels(zaxisID, plevel.data());

  auto plevel1 = plevel[0];
  auto plevelN = plevel[nlev - 1];
  if (plevel1 < plevelN) cdo_abort("The 3d pressure level data is upside down! Use the operator invertlev to invert the levels.");

  gridInqYvals(gridID, phi.data());

  auto units = cdo::inq_key_string(gridID, CDI_YAXIS, CDI_KEY_UNITS);

  if (units.rfind("degree", 0) == 0)
    for (size_t ilat = 0; ilat < nlat; ilat++) phi[ilat] *= DEG2RAD;

  for (size_t ilat = 0; ilat < nlat; ilat++) phi[ilat] = std::sin(phi[ilat]);

  for (size_t ilat = 0; ilat < nlat; ilat++) cosphi[ilat] = std::sqrt(1.0 - phi[ilat] * phi[ilat]);

  for (int ilev = 0; ilev < nlev; ilev++)
    for (size_t ilat = 0; ilat < nlat; ilat++) field2[ilev][ilat] = 0.0;

  if (nmiss == 0)
    {
      for (int ilev = nlev - 1; ilev >= 0; ilev--)
        for (int n = ilev; n < nlev - 1; ++n)
          for (size_t ilat = 0; ilat < nlat; ilat++)
            {
              field2[ilev][ilat] += fact * (field1[n][ilat] + field1[n + 1][ilat]) * cosphi[ilat] * (plevel[n] - plevel[n + 1]);
            }
    }
  else
    {
      for (size_t ilat = 0; ilat < nlat; ilat++)
        for (int ilev = nlev - 1; ilev >= 0; ilev--)
          for (int n = ilev; n < nlev - 1; ++n)
            {
              if (dbl_is_equal(field1[n][ilat], missval) || dbl_is_equal(field1[n + 1][ilat], missval))
                {
                  field2[ilev][ilat] = missval;
                  break;
                }
              else
                field2[ilev][ilat] += fact * (field1[n][ilat] + field1[n + 1][ilat]) * cosphi[ilat] * (plevel[n] - plevel[n + 1]);
            }
    }
}

void *
Mastrfu(void *process)
{
  cdo_initialize(process);

  auto streamID1 = cdo_open_read(0);

  operator_check_argc(0);

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);

  auto nvars = vlistNvars(vlistID1);
  if (nvars != 1) cdo_abort("This operator works only with one variable!");

  auto code = vlistInqVarCode(vlistID1, 0);
  if (code > 0 && code != 132) cdo_warning("Unexpected code %d!", code);

  auto missval = vlistInqVarMissval(vlistID1, 0);

  auto zaxisID = vlistInqVarZaxis(vlistID1, 0);
  if (zaxisInqType(zaxisID) != ZAXIS_PRESSURE && zaxisInqType(zaxisID) != ZAXIS_GENERIC)
    {
      cdo_warning("Unexpected vertical grid %s!", cdo::inq_key_string(zaxisID, CDI_GLOBAL, CDI_KEY_LONGNAME));
    }

  auto gridID = vlistInqVarGrid(vlistID1, 0);
  if (gridInqXsize(gridID) > 1) cdo_abort("Grid must be a zonal mean!");

  auto nlat = gridInqSize(gridID);
  auto nlev = zaxisInqSize(zaxisID);

  auto vlistID2 = vlistDuplicate(vlistID1);

  auto taxisID1 = vlistInqTaxis(vlistID1);
  auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  vlistDefVarCode(vlistID2, 0, 272);
  cdiDefKeyString(vlistID2, 0, CDI_KEY_NAME, "mastrfu");
  cdiDefKeyString(vlistID2, 0, CDI_KEY_LONGNAME, "mass stream function");
  cdiDefKeyString(vlistID2, 0, CDI_KEY_UNITS, "kg/s");
  vlistDefVarDatatype(vlistID2, 0, CDI_DATATYPE_FLT32);

  auto streamID2 = cdo_open_write(1);
  cdo_def_vlist(streamID2, vlistID2);

  Varray2D<double> array1(nlev, Varray<double>(nlat));
  Varray2D<double> array2(nlev, Varray<double>(nlat));

  int tsID = 0;
  while (true)
    {
      auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
      if (nrecs == 0) break;

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_def_timestep(streamID2, tsID);

      size_t nmiss = 0;
      for (int recID = 0; recID < nrecs; ++recID)
        {
          int varID, levelID;
          cdo_inq_record(streamID1, &varID, &levelID);
          size_t nmiss1;
          cdo_read_record(streamID1, array1[levelID].data(), &nmiss1);
          nmiss += nmiss1;
        }

      mastrfu(gridID, zaxisID, array1, array2, nmiss, missval);

      for (int recID = 0; recID < nrecs; ++recID)
        {
          int varID = 0;
          int levelID = recID;
          cdo_def_record(streamID2, varID, levelID);
          nmiss = array_num_mv(nlat, array2[levelID].data(), missval);
          cdo_write_record(streamID2, array2[levelID].data(), nmiss);
        }

      tsID++;
    }

  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  cdo_finish();

  return nullptr;
}
