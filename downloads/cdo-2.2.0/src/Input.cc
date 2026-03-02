/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Input     input          ASCII input
      Input     inputsrv       SERVICE input
      Input     inputext       EXTRA input
*/

#include <cdi.h>

#include <mpim_grid.h>
#include "process_int.h"
#include "griddes.h"
#include "cdo_zaxis.h"
#include "grid_read_pingo.h"
#include "cdo_default_values.h"

static size_t
input_iarray(size_t nval, int *array)
{
  size_t ival = 0;

  for (size_t i = 0; i < nval; ++i)
    {
      auto n = scanf("%d", &array[i]);
      if (n != 1) break;

      ival++;
    }

  return ival;
}

void *
Input(void *process)
{
  int varID = 0;
  size_t gridsize0 = 0, gridsize = 0;
  int taxisID = 0;
  CdoStreamID streamID = CDO_STREAM_UNDEF;
  int vlistID = -1;
  int code = 0, level = 0, date = 0, time = 0, nlon = 0, nlat = 0;
  int output_filetype = CDI_FILETYPE_GRB;
  int ihead[8];
  double missval = 0;
  std::vector<double> array;

  cdo_initialize(process);

  // clang-format off
  auto INPUT    = cdo_operator_add("input",    0, 0, nullptr);
  auto INPUTSRV = cdo_operator_add("inputsrv", 0, 0, nullptr);
  auto INPUTEXT = cdo_operator_add("inputext", 0, 0, nullptr);
  // clang-format on

  auto operatorID = cdo_operator_id();

  int gridID = -1;
  int zaxisID = -1;
  if (operatorID == INPUT)
    {
      operator_input_arg("grid description file or name");
      operator_check_argc((cdo_operator_argc() == 1) ? 1 : 2);

      gridID = cdo_define_grid(cdo_operator_argv(0));
      if (cdo_operator_argc() == 2) zaxisID = cdo_define_zaxis(cdo_operator_argv(1));
    }
  else { operator_check_argc(0); }

  const int nlevs = (zaxisID == -1) ? 1 : zaxisInqSize(zaxisID);

  double dlevel = 0;
  int nrecs = 0;

  int tsID = 0;
  while (true)
    {
      if (operatorID == INPUT)
        {
          output_filetype = cdo_filetype();

          code = -1;
          gridsize = gridInqSize(gridID);
          date = 0;
          time = 0;

          if (nrecs == 0) array.resize(gridsize * nlevs);

          cdo_print("Enter all %zu elements of timestep %d!", gridsize * nlevs, nrecs + 1);

          auto rval = input_darray(stdin, gridsize * nlevs, array);

          if (nrecs > 0 && rval == 0) break;

          if (rval != gridsize * nlevs) cdo_abort("Too few input elements (%zu of %zu)!", rval, gridsize * nlevs);

          if (feof(stdin)) break;
        }
      else if (operatorID == INPUTEXT)
        {
          output_filetype = CdoDefault::FileType;
          if (output_filetype == CDI_UNDEFID) output_filetype = CDI_FILETYPE_EXT;

          cdo_print("Enter header (date,code,level,gridsize) of record %d (or EOF(=^D))!", nrecs + 1);

          auto rval = input_iarray(4, ihead);
          if (feof(stdin) && nrecs == 0) cdo_abort("Too few header elements (%d of %d)!", rval, 4);
          if (feof(stdin)) break;
          if (rval != 4) cdo_abort("Invalid header input!");

          date = ihead[0];
          code = ihead[1];
          level = ihead[2];
          gridsize = ihead[3];

          time = 0;

          if (nrecs == 0)
            {
              dlevel = level;
              gridsize0 = gridsize;
              array.resize(gridsize);
              gridID = gridCreate(GRID_GENERIC, gridsize);
            }
          else
            {
              if (gridsize != gridsize0) cdo_abort("Gridsize must not change!");
            }

          cdo_print("Enter all %zu elements of record %d!", gridsize, nrecs + 1);

          rval = input_darray(stdin, gridsize, array);
          if (rval != gridsize) cdo_abort("Invalid data input!");
        }
      else if (operatorID == INPUTSRV)
        {
          output_filetype = CdoDefault::FileType;
          if (output_filetype == CDI_UNDEFID) output_filetype = CDI_FILETYPE_SRV;

          cdo_print("Enter header (code,level,date,time,nlon,nlat,dispo1,dispo2) of record %d (or EOF(=^D))!", nrecs + 1);

          auto rval = input_iarray(8, ihead);
          if (feof(stdin) && nrecs == 0) cdo_abort("Too few header elements (%d of %d)!", rval, 8);
          if (feof(stdin)) break;
          if (rval != 8) cdo_abort("Invalid header input!");

          code = ihead[0];
          level = ihead[1];
          date = ihead[2];
          time = ihead[3];
          nlon = ihead[4];
          nlat = ihead[5];

          gridsize = nlon * nlat;

          if (nrecs == 0)
            {
              dlevel = level;
              gridsize0 = gridsize;
              array.resize(gridsize);
              gridID = gridCreate(GRID_GENERIC, gridsize);
              gridDefXsize(gridID, nlon);
              gridDefYsize(gridID, nlat);
            }
          else
            {
              if (gridsize != gridsize0) cdo_abort("Gridsize must not change!");
            }

          cdo_print("Enter all %zu elements of record %d!", gridsize, nrecs + 1);

          rval = input_darray(stdin, gridsize, array);
          if (rval != gridsize) cdo_abort("Invalid data input!");
        }

      if (nrecs == 0)
        {
          if (zaxisID == -1)
            {
              zaxisID = zaxisCreate(ZAXIS_SURFACE, 1);
              zaxisDefLevels(zaxisID, &dlevel);
            }

          vlistID = vlistCreate();
          varID = vlistDefVar(vlistID, gridID, zaxisID, TIME_VARYING);
          vlistDefVarParam(vlistID, varID, cdiEncodeParam(code, 255, 255));

          missval = vlistInqVarMissval(vlistID, varID);

          taxisID = cdo_taxis_create(TAXIS_RELATIVE);
          vlistDefTaxis(vlistID, taxisID);

          streamID = cdo_open_write(0, output_filetype);
          cdo_def_vlist(streamID, vlistID);
        }

      CdiDateTime vDateTime{};
      vDateTime.date = cdiDate_set(date);
      vDateTime.time = cdiTime_set(time);
      taxisDefVdatetime(taxisID, vDateTime);
      cdo_def_timestep(streamID, tsID);

      for (int levelID = 0; levelID < nlevs; ++levelID)
        {
          auto offset = gridsize * levelID;
          auto nmiss = array_num_mv(gridsize, &array[offset], missval);
          cdo_def_record(streamID, varID, levelID);
          cdo_write_record(streamID, &array[offset], nmiss);
        }

      nrecs++;
      tsID++;
    }

  if (streamID)
    {
      cdo_stream_close(streamID);
      vlistDestroy(vlistID);
    }

  cdo_finish();

  return nullptr;
}
