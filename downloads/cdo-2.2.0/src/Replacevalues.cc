/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Setvals     setvals       Set list of old values to new values
      Setrtoc     setrtoc       Set range to new value
      Setrtoc2    setrtoc2      Set range to new value others to value2
*/

#include <cdi.h>

#include "process_int.h"
#include "param_conversion.h"

void *
Replacevalues(void *process)
{
  int nvals = 0;
  std::vector<double> fltarr;
  double rmin = 0, rmax = 0;
  double newval = 0, newval2 = 0;

  cdo_initialize(process);

  // clang-format off
  const auto SETVALS  = cdo_operator_add("setvals",  0, 0, "I1,O1,...,In,On");
  const auto SETRTOC  = cdo_operator_add("setrtoc",  0, 0, "range (min, max), value");
  const auto SETRTOC2 = cdo_operator_add("setrtoc2", 0, 0, "range (min, max), value1, value2");
  // clang-format on

  const auto operatorID = cdo_operator_id();

  operator_input_arg(cdo_operator_enter(operatorID));

  if (operatorID == SETVALS)
    {
      fltarr = cdo_argv_to_flt(cdo_get_oper_argv());
      nvals = fltarr.size();
      if (nvals < 2) cdo_abort("Too few arguments!");
      if (nvals % 2 != 0) cdo_abort("Need pairs of arguments!");
      nvals = nvals / 2;
    }
  else if (operatorID == SETRTOC)
    {
      operator_check_argc(3);
      rmin = parameter_to_double(cdo_operator_argv(0));
      rmax = parameter_to_double(cdo_operator_argv(1));
      newval = parameter_to_double(cdo_operator_argv(2));
    }
  else if (operatorID == SETRTOC2)
    {
      operator_check_argc(4);
      rmin = parameter_to_double(cdo_operator_argv(0));
      rmax = parameter_to_double(cdo_operator_argv(1));
      newval = parameter_to_double(cdo_operator_argv(2));
      newval2 = parameter_to_double(cdo_operator_argv(3));
    }

  const auto streamID1 = cdo_open_read(0);

  const auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  const auto vlistID2 = vlistDuplicate(vlistID1);

  const auto taxisID1 = vlistInqTaxis(vlistID1);
  const auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  const auto streamID2 = cdo_open_write(1);

  cdo_def_vlist(streamID2, vlistID2);

  const auto gridsizemax = vlistGridsizeMax(vlistID1);
  Varray<double> array(gridsizemax);

  VarList varList1;
  varListInit(varList1, vlistID1);

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
          size_t nmiss;
          cdo_read_record(streamID1, array.data(), &nmiss);

          const auto gridsize = varList1[varID].gridsize;
          const auto missval = varList1[varID].missval;

          if (operatorID == SETVALS)
            {
              for (size_t i = 0; i < gridsize; ++i)
                if (!DBL_IS_EQUAL(array[i], missval))
                  {
                    for (int j = 0; j < nvals; ++j)
                      {
                        if (DBL_IS_EQUAL(array[i], fltarr[j * 2]))
                          {
                            array[i] = fltarr[j * 2 + 1];
                            break;
                          }
                      }
                  }
            }
          else if (operatorID == SETRTOC)
            {
              for (size_t i = 0; i < gridsize; ++i)
                if (!DBL_IS_EQUAL(array[i], missval))
                  {
                    if (array[i] >= rmin && array[i] <= rmax) array[i] = newval;
                  }
            }
          else if (operatorID == SETRTOC2)
            {
              for (size_t i = 0; i < gridsize; ++i)
                if (!DBL_IS_EQUAL(array[i], missval)) { array[i] = (array[i] >= rmin && array[i] <= rmax) ? newval : newval2; }
            }

          cdo_def_record(streamID2, varID, levelID);
          cdo_write_record(streamID2, array.data(), nmiss);
        }

      tsID++;
    }

  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  cdo_finish();

  return nullptr;
}
