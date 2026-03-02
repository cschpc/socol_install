/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Ninfo      npar            Number of parameters
      Ninfo      nlevel          Number of levels
      Ninfo      nyear           Number of years
      Ninfo      nmon            Number of months
      Ninfo      ndate           Number of dates
      Ninfo      ntime           Number of timesteps
      Ninfo      ngridpoints     Number of gridpoints
      Ninfo      ngrids          Number of grids
*/

#include <cdi.h>

#include "process_int.h"

void *
Ninfo(void *process)
{
  enum
  {
    NYEAR,
    NMON,
    NDATE,
    NTIME,
    NPAR,
    NLEVEL,
    NGRIDPOINTS,
    NGRIDS
  };

  cdo_initialize(process);

  // clang-format off
  cdo_operator_add("nyear"      , NYEAR      , 0, nullptr);
  cdo_operator_add("nmon"       , NMON       , 0, nullptr);
  cdo_operator_add("ndate"      , NDATE      , 0, nullptr);
  cdo_operator_add("ntime"      , NTIME      , 0, nullptr);
  cdo_operator_add("npar"       , NPAR       , 0, nullptr);
  cdo_operator_add("nlevel"     , NLEVEL     , 0, nullptr);
  cdo_operator_add("ngridpoints", NGRIDPOINTS, 0, nullptr);
  cdo_operator_add("ngrids"     , NGRIDS     , 0, nullptr);
  // clang-format on

  auto operatorID = cdo_operator_id();
  auto operfunc = cdo_operator_f1(operatorID);

  operator_check_argc(0);

  auto streamID = cdo_open_read(0);
  auto vlistID = cdo_stream_inq_vlist(streamID);

  auto nvars = vlistNvars(vlistID);
  auto taxisID = vlistInqTaxis(vlistID);
  auto ntsteps = vlistNtsteps(vlistID);
  auto ngrids = vlistNgrids(vlistID);

  VarList varList;
  varListInit(varList, vlistID);

  switch (operfunc)
    {
    case NYEAR:
      {
        int year0 = 0;
        int nyear = 0;
        int tsID = 0;
        if (ntsteps != 0)
          while (cdo_stream_inq_timestep(streamID, tsID))
            {
              int year = taxisInqVdatetime(taxisID).date.year;

              if (tsID == 0 || year0 != year)
                {
                  year0 = year;
                  nyear++;
                }

              tsID++;
            }
        fprintf(stdout, "%d\n", nyear);
        break;
      }
    case NMON:
      {
        int month0 = 0;
        int nmonth = 0;
        int tsID = 0;
        if (ntsteps != 0)
          while (cdo_stream_inq_timestep(streamID, tsID))
            {
              int month = taxisInqVdatetime(taxisID).date.month;
              if (tsID == 0 || month0 != month)
                {
                  month0 = month;
                  nmonth++;
                }

              tsID++;
            }
        fprintf(stdout, "%d\n", nmonth);
        break;
      }
    case NDATE:
      {
        CdiDate date0 = {};
        int ndate = 0;
        int tsID = 0;
        if (ntsteps != 0)
          while (cdo_stream_inq_timestep(streamID, tsID))
            {
              auto vDate = taxisInqVdatetime(taxisID).date;
              if (tsID == 0 || !cdiDate_isEQ(date0, vDate))
                {
                  date0 = vDate;
                  ndate++;
                }

              tsID++;
            }
        fprintf(stdout, "%d\n", ndate);
        break;
      }
    case NTIME:
      {
        int tsID = (ntsteps > 0) ? ntsteps : 0;
        if (tsID == 0)
          while (cdo_stream_inq_timestep(streamID, tsID)) tsID++;
        fprintf(stdout, "%d\n", tsID);
        break;
      }
    case NPAR: fprintf(stdout, "%d\n", nvars); break;
    case NLEVEL:
      for (int varID = 0; varID < nvars; ++varID) { fprintf(stdout, "%d\n", varList[varID].nlevels); }
      break;
    case NGRIDPOINTS:
      for (int varID = 0; varID < nvars; ++varID) { fprintf(stdout, "%zu\n", varList[varID].gridsize); }
      break;
    case NGRIDS: fprintf(stdout, "%d\n", ngrids); break;
    default: cdo_abort("operator not implemented!"); break;
    }

  cdo_stream_close(streamID);

  cdo_finish();

  return nullptr;
}
