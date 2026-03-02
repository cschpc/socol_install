/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

        Timstat3        varquot2test
        Timstat3        meandiff2test
*/

#include <cdi.h>

#include "process_int.h"
#include "cdo_vlist.h"
#include "param_conversion.h"
#include "statistic.h"

#define NIN 2
#define NOUT 1
#define NFWORK 4
#define NIWORK 2

void *
Timstat3(void *process)
{
  CdoStreamID streamID[NIN];
  CdiDateTime vDateTime{};
  int vlistID[NIN], vlistID2 = -1;
  int is;
  Varray3D<int> iwork[NIWORK];
  FieldVector2D fwork[NFWORK];
  int reached_eof[NIN];
  constexpr int n_in = NIN;

  cdo_initialize(process);

  // clang-format off
  auto VARQUOT2TEST  = cdo_operator_add("varquot2test",  0, 0, nullptr);
  auto MEANDIFF2TEST = cdo_operator_add("meandiff2test", 0, 0, nullptr);
  // clang-format on

  auto operatorID = cdo_operator_id();

  operator_input_arg("constant and risk (e.g. 0.05)");
  operator_check_argc(2);
  auto rconst = parameter_to_double(cdo_operator_argv(0));
  auto risk = parameter_to_double(cdo_operator_argv(1));

  if (operatorID == VARQUOT2TEST)
    {
      if (rconst <= 0) cdo_abort("Constant must be positive!");
      if (risk <= 0 || risk >= 1) cdo_abort("Risk must be greater than 0 and lower than 1!");
    }

  for (is = 0; is < NIN; ++is)
    {
      streamID[is] = cdo_open_read(is);
      vlistID[is] = cdo_stream_inq_vlist(streamID[is]);
      if (is > 0)
        {
          vlistID2 = cdo_stream_inq_vlist(streamID[is]);
          vlist_compare(vlistID[0], vlistID2, CMP_ALL);
        }
    }

  auto vlistID3 = vlistDuplicate(vlistID[0]);

  auto gridsizemax = vlistGridsizeMax(vlistID[0]);
  auto nvars = vlistNvars(vlistID[0]);

  VarList varList0;
  varListInit(varList0, vlistID[0]);

  auto maxrecs = vlistNrecs(vlistID[0]);
  std::vector<RecordInfo> recList(maxrecs);

  auto taxisID1 = vlistInqTaxis(vlistID[0]);
  auto taxisID3 = taxisDuplicate(taxisID1);

  vlistDefTaxis(vlistID3, taxisID3);
  auto streamID3 = cdo_open_write(2);
  cdo_def_vlist(streamID3, vlistID3);

  for (int i = 0; i < NIN; ++i) reached_eof[i] = 0;

  Field in[NIN], out[NOUT];
  for (int i = 0; i < NIN; ++i) in[i].resize(gridsizemax);
  for (int i = 0; i < NOUT; ++i) out[i].resize(gridsizemax);

  for (int iw = 0; iw < NFWORK; ++iw) fwork[iw].resize(nvars);
  for (int iw = 0; iw < NIWORK; ++iw) iwork[iw].resize(nvars);

  for (int varID = 0; varID < nvars; ++varID)
    {
      const auto &var = varList0[varID];

      auto gridsize = vlistGridsizeMax(vlistID[0]);

      for (int iw = 0; iw < NFWORK; ++iw) fwork[iw][varID].resize(var.nlevels);
      for (int iw = 0; iw < NIWORK; ++iw) iwork[iw][varID].resize(var.nlevels);

      for (int levelID = 0; levelID < var.nlevels; ++levelID)
        {
          for (int iw = 0; iw < NFWORK; ++iw)
            {
              fwork[iw][varID][levelID].grid = var.gridID;
              fwork[iw][varID][levelID].missval = var.missval;
              fwork[iw][varID][levelID].resize(gridsize);
            }

          for (int iw = 0; iw < NIWORK; ++iw) iwork[iw][varID][levelID].resize(gridsize, 0);
        }
    }

  int tsID = 0;
  while (true)
    {
      for (is = 0; is < NIN; ++is)
        {
          if (reached_eof[is]) continue;

          auto nrecs = cdo_stream_inq_timestep(streamID[is], tsID);
          if (nrecs == 0)
            {
              reached_eof[is] = 1;
              continue;
            }

          vDateTime = taxisInqVdatetime(taxisID1);

          for (int recID = 0; recID < nrecs; ++recID)
            {
              int varID, levelID;
              cdo_inq_record(streamID[is], &varID, &levelID);

              auto gridsize = gridInqSize(vlistInqVarGrid(vlistID[is], varID));

              in[is].missval = vlistInqVarMissval(vlistID[is], varID);

              if (tsID == 0 && is == 0) recList[recID].set(varID, levelID);

              cdo_read_record(streamID[is], in[is].vec_d.data(), &in[is].nmiss);

              for (size_t i = 0; i < gridsize; ++i)
                {
                  // if ( ( ! dbl_is_equal(array1[i], missval1) ) && ( ! dbl_is_equal(array2[i], missval2) ) )
                  {
                    fwork[NIN * is + 0][varID][levelID].vec_d[i] += in[is].vec_d[i];
                    fwork[NIN * is + 1][varID][levelID].vec_d[i] += in[is].vec_d[i] * in[is].vec_d[i];
                    iwork[is][varID][levelID][i]++;
                  }
                }
            }
        }

      for (is = 0; is < NIN; ++is)
        if (!reached_eof[is]) break;

      if (is == NIN) break;

      tsID++;
    }

  taxisDefVdatetime(taxisID3, vDateTime);
  cdo_def_timestep(streamID3, 0);

  for (int recID = 0; recID < maxrecs; ++recID)
    {
      auto [varID, levelID] = recList[recID].get();

      auto missval1 = fwork[0][varID][levelID].missval;
      auto missval2 = missval1;

      if (operatorID == VARQUOT2TEST)
        {
          auto fwork0 = fwork[0][varID][levelID];
          auto fwork1 = fwork[1][varID][levelID];
          auto fwork2 = fwork[2][varID][levelID];
          auto fwork3 = fwork[3][varID][levelID];
          for (size_t i = 0; i < gridsizemax; ++i)
            {
              auto fnvals0 = iwork[0][varID][levelID][i];
              auto fnvals1 = iwork[1][varID][levelID][i];

              auto temp0 = DIVMN(MULMN(fwork0.vec_d[i], fwork0.vec_d[i]), fnvals0);
              auto temp1 = DIVMN(MULMN(fwork2.vec_d[i], fwork2.vec_d[i]), fnvals1);
              auto temp2 = SUBMN(fwork1.vec_d[i], temp0);
              auto temp3 = SUBMN(fwork3.vec_d[i], temp1);
              auto statistic = DIVMN(temp2, ADDMN(temp2, MULMN(rconst, temp3)));

              double fractil_1 = missval1, fractil_2 = missval1;
              if (fnvals0 > 1 && fnvals1 > 1)
                cdo::beta_distr_constants((fnvals0 - 1) / 2, (fnvals1 - 1) / 2, 1 - risk, &fractil_1, &fractil_2);

              out[0].vec_d[i] = dbl_is_equal(statistic, missval1)                    ? missval1
                                : (statistic <= fractil_1 || statistic >= fractil_2) ? 1.0
                                                                                     : 0.0;
            }
        }
      else if (operatorID == MEANDIFF2TEST)
        {
          double mean_factor[NIN], var_factor[NIN];

          mean_factor[0] = 1.0;
          mean_factor[1] = -1.0;
          var_factor[0] = var_factor[1] = 1.0;

          for (size_t i = 0; i < gridsizemax; ++i)
            {
              double temp0 = 0.0;
              double deg_of_freedom = -n_in;
              for (int j = 0; j < n_in; ++j)
                {
                  auto fnvals = iwork[j][varID][levelID][i];
                  auto tmp = DIVMN(MULMN(fwork[2 * j][varID][levelID].vec_d[i], fwork[2 * j][varID][levelID].vec_d[i]), fnvals);
                  temp0 = ADDMN(temp0, DIVMN(SUBMN(fwork[2 * j + 1][varID][levelID].vec_d[i], tmp), var_factor[j]));
                  deg_of_freedom = ADDMN(deg_of_freedom, fnvals);
                }

              if (!dbl_is_equal(temp0, missval1) && temp0 < 0) temp0 = 0;  // This is possible because of rounding errors

              auto stddev_estimator = SQRTMN(DIVMN(temp0, deg_of_freedom));
              auto mean_estimator = -rconst;
              for (int j = 0; j < n_in; ++j)
                {
                  auto fnvals = iwork[j][varID][levelID][i];
                  mean_estimator
                      = ADDMN(mean_estimator, MULMN(mean_factor[j], DIVMN(fwork[2 * j][varID][levelID].vec_d[i], fnvals)));
                }

              double temp1 = 0.0;
              for (int j = 0; j < n_in; ++j)
                {
                  auto fnvals = iwork[j][varID][levelID][i];
                  temp1 = ADDMN(temp1, DIVMN(MUL(MUL(mean_factor[j], mean_factor[j]), var_factor[j]), fnvals));
                }

              auto norm = SQRTMN(temp1);

              auto temp2 = DIVMN(DIVMN(mean_estimator, norm), stddev_estimator);
              auto fractil = (deg_of_freedom < 1) ? missval1 : cdo::student_t_inv(deg_of_freedom, 1 - risk / 2);

              out[0].vec_d[i]
                  = dbl_is_equal(temp2, missval1) || dbl_is_equal(fractil, missval1) ? missval1 : std::fabs(temp2) >= fractil;
            }
        }

      cdo_def_record(streamID3, varID, levelID);
      cdo_write_record(streamID3, out[0].vec_d.data(), field_num_miss(out[0]));
    }

  cdo_stream_close(streamID3);
  for (is = 0; is < NIN; ++is) cdo_stream_close(streamID[is]);

  cdo_finish();

  return nullptr;
}
