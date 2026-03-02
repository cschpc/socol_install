/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Setzaxis   setzaxis        Set zaxis
*/

#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "param_conversion.h"
#include "cdo_zaxis.h"

int
getkeyval_dp(const std::string &p_keyval, const char *p_key, double *val)
{
  int status = 0;
  auto keyval = p_keyval.c_str();
  auto keylen = strlen(p_key);

  if (strncmp(keyval, p_key, keylen) == 0)
    {
      auto pkv = keyval + keylen;
      if (pkv[0] == '=' && pkv[1] != 0)
        {
          *val = parameter_to_double(&pkv[1]);
          status = 1;
        }
      else { cdo_abort("Syntax error for parameter %s!", keyval); }
    }

  return status;
}

void *
Setzaxis(void *process)
{
  int zaxisID1, zaxisID2 = -1;
  int nzaxis, index;
  bool hasZtop = false, hasZbot = false;
  double ztop = 0, zbot = 0;

  cdo_initialize(process);

  // clang-format off
  auto SETZAXIS       = cdo_operator_add("setzaxis",        0, 0, "zaxis description file");
  auto GENLEVELBOUNDS = cdo_operator_add("genlevelbounds",  0, 0, nullptr);
  // clang-format on

  int operatorID = cdo_operator_id();

  if (operatorID == SETZAXIS)
    {
      operator_input_arg(cdo_operator_enter(operatorID));
      operator_check_argc(1);
      zaxisID2 = cdo_define_zaxis(cdo_operator_argv(0));
    }
  else if (operatorID == GENLEVELBOUNDS)
    {
      auto npar = cdo_operator_argc();
      auto parnames = cdo_get_oper_argv();

      for (int i = 0; i < npar; ++i)
        {
          if (Options::cdoVerbose) cdo_print("keyval[%d]: %s", i + 1, parnames[i]);

          if (!hasZbot && getkeyval_dp(parnames[i], "zbot", &zbot))
            hasZbot = true;
          else if (!hasZtop && getkeyval_dp(parnames[i], "ztop", &ztop))
            hasZtop = true;
          else
            cdo_abort("Parameter >%s< unsupported! Supported parameter are: zbot, ztop", parnames[i]);
        }
    }

  auto streamID1 = cdo_open_read(0);

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto vlistID2 = vlistDuplicate(vlistID1);

  auto taxisID1 = vlistInqTaxis(vlistID1);
  auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  auto streamID2 = cdo_open_write(1);

  if (operatorID == SETZAXIS)
    {
      int found = 0;
      nzaxis = vlistNzaxis(vlistID1);
      for (index = 0; index < nzaxis; ++index)
        {
          zaxisID1 = vlistZaxis(vlistID1, index);

          if (zaxisInqSize(zaxisID1) == zaxisInqSize(zaxisID2))
            {
              vlistChangeZaxisIndex(vlistID2, index, zaxisID2);
              found++;
            }
        }
      if (!found) cdo_warning("No zaxis with %d levels found!", zaxisInqSize(zaxisID2));
    }
  else if (operatorID == GENLEVELBOUNDS)
    {
      nzaxis = vlistNzaxis(vlistID1);
      for (index = 0; index < nzaxis; ++index)
        {
          zaxisID1 = vlistZaxis(vlistID1, index);
          auto nlevels = zaxisInqSize(zaxisID1);
          if (nlevels > 1)
            {
              Varray<double> levels(nlevels), lbounds(nlevels), ubounds(nlevels);

              cdo_zaxis_inq_levels(zaxisID1, levels.data());
              zaxisID2 = zaxisDuplicate(zaxisID1);
              if (!zaxisInqLevels(zaxisID1, nullptr)) zaxisDefLevels(zaxisID2, levels.data());

              gen_layer_bounds(nlevels, levels, lbounds, ubounds);

              auto isPositive = !(levels[0] < 0.0 && levels[nlevels - 1] < 0.0);
              auto isReverse = (levels[0] > levels[nlevels - 1]);
              auto positveIsDown = (isPositive && !isReverse && positive_is_down(zaxisID1));

              if (hasZbot)
                {
                  if (positveIsDown)
                    ubounds[nlevels - 1] = zbot;
                  else
                    lbounds[0] = zbot;
                }
              if (hasZtop)
                {
                  if (positveIsDown)
                    lbounds[0] = ztop;
                  else
                    ubounds[nlevels - 1] = ztop;
                }
              zaxisDefLbounds(zaxisID2, lbounds.data());
              zaxisDefUbounds(zaxisID2, ubounds.data());
              vlistChangeZaxisIndex(vlistID2, index, zaxisID2);
            }
        }
    }

  cdo_def_vlist(streamID2, vlistID2);

  auto gridsize = vlistGridsizeMax(vlistID1);
  if (vlistNumber(vlistID1) != CDI_REAL) gridsize *= 2;
  Varray<double> array(gridsize);

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
          cdo_def_record(streamID2, varID, levelID);

          size_t nmiss;
          cdo_read_record(streamID1, array.data(), &nmiss);
          cdo_write_record(streamID2, array.data(), nmiss);
        }

      tsID++;
    }

  cdo_stream_close(streamID1);
  cdo_stream_close(streamID2);

  cdo_finish();

  return nullptr;
}
