/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Set        setcode         Set code number
      Set        setparam        Set parameter identifier
      Set        setname         Set variable name
      Set        setlevel        Set level
      Set        setltype        Set GRIB level type
*/

#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "param_conversion.h"
#include "cdo_zaxis.h"

static void
set_level(int vlistID2, double newlevel)
{
  const auto nzaxis = vlistNzaxis(vlistID2);
  for (int index = 0; index < nzaxis; ++index)
    {
      const auto zaxisID1 = vlistZaxis(vlistID2, index);
      const auto zaxisID2 = zaxisDuplicate(zaxisID1);
      const auto nlevs = zaxisInqSize(zaxisID2);
      Varray<double> levels(nlevs);
      cdo_zaxis_inq_levels(zaxisID2, levels.data());
      levels[0] = newlevel;
      zaxisDefLevels(zaxisID2, levels.data());
      vlistChangeZaxis(vlistID2, zaxisID1, zaxisID2);
    }
}

static void
set_ltype(int vlistID2, double newval)
{
  const auto nzaxis = vlistNzaxis(vlistID2);
  for (int index = 0; index < nzaxis; ++index)
    {
      const auto zaxisID1 = vlistZaxis(vlistID2, index);
      const auto zaxisID2 = zaxisDuplicate(zaxisID1);
      const auto zaxistype = ZAXIS_GENERIC;
      zaxisChangeType(zaxisID2, zaxistype);
      cdiDefKeyInt(zaxisID2, CDI_GLOBAL, CDI_KEY_TYPEOFFIRSTFIXEDSURFACE, newval);
      vlistChangeZaxis(vlistID2, zaxisID1, zaxisID2);
    }
}

void *
Set(void *process)
{
  int maxSteps = -1;
  int newval = -1, tabnum = 0;
  int newparam = 0;
  const char *newname = nullptr, *newunit = nullptr;
  double newlevel = 0;

  cdo_initialize(process);

  // clang-format off
  const auto SETCODE     = cdo_operator_add("setcode",     0, 0, "code number");
  const auto SETPARAM    = cdo_operator_add("setparam",    0, 0, "parameter identifier (format: code[.tabnum] or num[.cat[.dis]])");
  const auto SETNAME     = cdo_operator_add("setname",     0, 0, "variable name");
  const auto SETUNIT     = cdo_operator_add("setunit",     0, 0, "variable unit");
  const auto SETLEVEL    = cdo_operator_add("setlevel",    0, 0, "level");
  const auto SETLTYPE    = cdo_operator_add("setltype",    0, 0, "GRIB level type");
  const auto SETTABNUM   = cdo_operator_add("settabnum",   0, 0, "GRIB table number");
  const auto SETMAXSTEPS = cdo_operator_add("setmaxsteps", 0, 0, "max. number of timesteps");
  // clang-format on

  const auto operatorID = cdo_operator_id();

  operator_input_arg(cdo_operator_enter(operatorID));
  if (operatorID == SETCODE || operatorID == SETLTYPE) { newval = parameter_to_int(cdo_operator_argv(0)); }
  else if (operatorID == SETPARAM) { newparam = string_to_param(cdo_operator_argv(0)); }
  else if (operatorID == SETNAME) { newname = cdo_operator_argv(0).c_str(); }
  else if (operatorID == SETUNIT) { newunit = cdo_operator_argv(0).c_str(); }
  else if (operatorID == SETTABNUM) { tabnum = parameter_to_int(cdo_operator_argv(0)); }
  else if (operatorID == SETLEVEL) { newlevel = parameter_to_double(cdo_operator_argv(0)); }
  else if (operatorID == SETMAXSTEPS) { maxSteps = parameter_to_int(cdo_operator_argv(0)); }

  const auto streamID1 = cdo_open_read(0);

  const auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  const auto vlistID2 = vlistDuplicate(vlistID1);
  // vlistPrint(vlistID2);

  const auto taxisID1 = vlistInqTaxis(vlistID1);
  const auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  if (operatorID == SETCODE)
    {
      const auto nvars = vlistNvars(vlistID2);
      for (int varID = 0; varID < nvars; ++varID) vlistDefVarCode(vlistID2, varID, newval);
    }
  else if (operatorID == SETPARAM) { vlistDefVarParam(vlistID2, 0, newparam); }
  else if (operatorID == SETNAME) { cdiDefKeyString(vlistID2, 0, CDI_KEY_NAME, newname); }
  else if (operatorID == SETUNIT) { cdiDefKeyString(vlistID2, 0, CDI_KEY_UNITS, newunit); }
  else if (operatorID == SETTABNUM)
    {
      const auto tableID = tableDef(-1, tabnum, nullptr);
      const auto nvars = vlistNvars(vlistID2);
      for (int varID = 0; varID < nvars; ++varID) vlistDefVarTable(vlistID2, varID, tableID);
    }
  else if (operatorID == SETLEVEL) { set_level(vlistID2, newlevel); }
  else if (operatorID == SETLTYPE) { set_ltype(vlistID2, newval); }
  else if (operatorID == SETMAXSTEPS) { vlistDefNtsteps(vlistID2, maxSteps); }

  VarList varList1;
  varListInit(varList1, vlistID1);

  Field field;

  const auto streamID2 = cdo_open_write(1);
  cdo_def_vlist(streamID2, vlistID2);

  int tsID = 0;
  while (true)
    {
      const auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
      if (nrecs == 0) break;

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_def_timestep(streamID2, tsID);

      for (int recID = 0; recID < nrecs; ++recID)
        {
          int varID, levelID;
          cdo_inq_record(streamID1, &varID, &levelID);
          cdo_def_record(streamID2, varID, levelID);

          field.init(varList1[varID]);
          cdo_read_record(streamID1, field);
          cdo_write_record(streamID2, field);
        }

      tsID++;
    }

  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  cdo_finish();

  return nullptr;
}
