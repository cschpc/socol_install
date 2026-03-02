/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Arithc     addc            Add by constant
      Arithc     subc            Subtract by constant
      Arithc     mulc            Multiply by constant
      Arithc     divc            Divide by constant
      Arithc     mod             Modulo operator
*/

#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "cdo_vlist.h"
#include "param_conversion.h"
#include "field_functions.h"

static void
fill_vars(const VarList &varList, std::vector<bool> &vars)
{
  auto nvars = vars.size();

  if (Options::cdo_num_varnames() > 0)
    {
      auto found = false;
      for (size_t varID = 0; varID < nvars; ++varID)
        {
          vars[varID] = false;
          for (size_t i = 0; i < Options::cdo_num_varnames(); ++i)
            if (varList[varID].name == Options::cdoVarnames[i])
              {
                vars[varID] = true;
                found = true;
                break;
              }
        }

      if (!found) cdo_abort("Variable %s%s not found!", Options::cdoVarnames[0], (Options::cdo_num_varnames() > 1) ? ",..." : "");
    }
  else
    {
      for (size_t varID = 0; varID < nvars; ++varID) vars[varID] = true;
    }
}

static void
addOperators(void)
{
  // clang-format off
  cdo_operator_add("addc", FieldFunc_Add, 1, "constant value");
  cdo_operator_add("subc", FieldFunc_Sub, 1, "constant value");
  cdo_operator_add("mulc", FieldFunc_Mul, 1, "constant value");
  cdo_operator_add("divc", FieldFunc_Div, 1, "constant value");
  cdo_operator_add("minc", FieldFunc_Min, 0, "constant value");
  cdo_operator_add("maxc", FieldFunc_Max, 0, "constant value");
  cdo_operator_add("mod",  FieldFunc_Mod, 0, "divisor");
  // clang-format on
}

void *
Arithc(void *process)
{
  cdo_initialize(process);

  addOperators();

  auto operatorID = cdo_operator_id();
  auto operfunc = cdo_operator_f1(operatorID);
  const bool opercplx = cdo_operator_f2(operatorID);

  operator_input_arg(cdo_operator_enter(operatorID));
  if (cdo_operator_argc() < 1) cdo_abort("Too few arguments!");
  if (cdo_operator_argc() > 2) cdo_abort("Too many arguments!");
  auto rconst = parameter_to_double(cdo_operator_argv(0));

  double rconstcplx[2] = { rconst, 0.0 };
  if (cdo_operator_argc() == 2) rconstcplx[1] = parameter_to_double(cdo_operator_argv(1));

  auto streamID1 = cdo_open_read(0);

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto vlistID2 = vlistDuplicate(vlistID1);

  VarList varList1;
  varListInit(varList1, vlistID1);

  auto nvars = vlistNvars(vlistID1);
  // for (int varID = 0; varID < nvars; ++varID) varList1[varID].memType = MemType::Double;

  std::vector<bool> vars(nvars);
  fill_vars(varList1, vars);

  auto taxisID1 = vlistInqTaxis(vlistID1);
  auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  auto streamID2 = cdo_open_write(1);
  cdo_def_vlist(streamID2, vlistID2);

  auto nwpv = (vlistNumber(vlistID1) == CDI_COMP) ? 2 : 1;
  if (nwpv == 2 && !opercplx) cdo_abort("Fields with complex numbers are not supported by this operator!");

  Field field;

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
          field.init(varList1[varID]);
          cdo_read_record(streamID1, field);

          if (vars[varID])
            {
              if (field.nwpv == 2)
                fieldc_function_complex(field, rconstcplx, operfunc);
              else
                fieldc_function(field, rconst, operfunc);

              // recalculate number of missing values
              field_num_mv(field);
            }

          cdo_def_record(streamID2, varID, levelID);
          cdo_write_record(streamID2, field);
        }

      tsID++;
    }

  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  vlistDestroy(vlistID2);

  cdo_finish();

  return nullptr;
}
