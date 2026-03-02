/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Harmonic   harmonic        Harmonic
*/

#include <cdi.h>

#include "process_int.h"
#include "cdo_vlist.h"
#include "param_conversion.h"
#include "util_files.h"
#include "varray.h"

void *
Harmonic(void *process)
{
  CdiDateTime vDateTime{};

  cdo_initialize(process);

  operator_input_arg("wave number and wave length of first harmonic in number of timesteps");

  operator_check_argc(2);

  auto n_out = parameter_to_int(cdo_operator_argv(0));
  auto n = parameter_to_int(cdo_operator_argv(1));

  if (n_out > 9) cdo_abort("Maximum number of wave numbers is 9!");

  if (n < 1 || n < 2 * n_out)
    cdo_abort("The wave length must be positive and smaller than 2 times the number of requested harmonics (=%d)!", n_out);

  const auto streamID1 = cdo_open_read(0);

  const auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  const auto vlistID2 = vlistDuplicate(vlistID1);

  VarList varList1, varList2;
  varListInit(varList1, vlistID1);
  varListInit(varList2, vlistID2);

  const auto taxisID1 = vlistInqTaxis(vlistID1);
  const auto taxisID2 = cdo_taxis_create(TAXIS_ABSOLUTE);
  vlistDefTaxis(vlistID2, taxisID2);

  auto refname = cdo_get_stream_name(0);
  char filesuffix[32] = { 0 };
  FileUtils::gen_suffix(filesuffix, sizeof(filesuffix), cdo_inq_filetype(streamID1), vlistID1, refname);

  std::vector<CdoStreamID> streamIDs(n_out);

  char filename[8192];
  strcpy(filename, cdo_get_stream_name(1));
  int nchars = strlen(filename);

  for (int j = 0; j < n_out; ++j)
    {
      sprintf(filename + nchars, "%1d", j + 1);
      if (filesuffix[0]) sprintf(filename + nchars + 1, "%s", filesuffix);

      const auto streamID2 = cdo_open_write(filename);
      cdo_def_vlist(streamID2, vlistID2);
      streamIDs[j] = streamID2;
    }

  const auto nvars = vlistNvars(vlistID1);

  Varray3D<double> out(n_out);
  Varray3D<double> work(2 * n_out);

  for (int j = 0; j < n_out; ++j)
    {
      out[j].resize(nvars);
      for (int varID = 0; varID < nvars; ++varID) out[j][varID].resize(varList1[varID].nlevels * varList1[varID].gridsize);
    }

  for (int j = 0; j < n_out * 2; ++j)
    {
      work[j].resize(nvars);
      for (int varID = 0; varID < nvars; ++varID)
        {
          const auto gridsize = varList1[varID].gridsize;
          const auto nlevels = varList1[varID].nlevels;
          work[j][varID].resize(gridsize * nlevels);
          varray_fill(gridsize * nlevels, work[j][varID], 0.0);
        }
    }

  const auto gridsizemax = vlistGridsizeMax(vlistID1);
  Varray<double> array(gridsizemax);

  int tsID = 0;
  while (true)
    {
      const auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
      if (nrecs == 0) break;

      if (tsID == 0) vDateTime = taxisInqVdatetime(taxisID1);

      for (int recID = 0; recID < nrecs; ++recID)
        {
          int varID, levelID;
          cdo_inq_record(streamID1, &varID, &levelID);
          size_t nmiss;
          cdo_read_record(streamID1, array.data(), &nmiss);

          if (nmiss) cdo_abort("Missing values are not allowed!");

          const auto gridsize = varList1[varID].gridsize;
          const auto offset = gridsize * levelID;

          for (int j = 0; j < n_out; ++j)
            {
              const auto scarg = 2 * M_PI * (((j + 1) * (tsID + 1)) % n) / n;
              const auto sine = std::sin(scarg);
              const auto cosine = std::cos(scarg);
              for (size_t i = 0; i < gridsize; ++i)
                {
                  work[j][varID][i + offset] += array[i] * sine;
                  work[n_out + j][varID][i + offset] += array[i] * cosine;
                }
            }
        }

      tsID++;
    }

  const auto nts = tsID;

  cdo_stream_close(streamID1);

  if (nts % n) { cdo_abort("The length of first harmonic (=%d) does not divide the number of timesteps (=%d)!", n, nts); }

  for (int j = 0; j < n_out && 2 * (j + 1) < n; ++j)
    {
      for (int varID = 0; varID < nvars; ++varID)
        {
          const auto gridsize = varList2[varID].gridsize;
          const auto nlevels = varList2[varID].nlevels;
          for (int levelID = 0; levelID < nlevels; ++levelID)
            {
              const auto offset = gridsize * levelID;
              for (size_t i = 0; i < gridsize; ++i)
                {
                  const auto sqrwork1 = work[j][varID][i + offset] * work[j][varID][i + offset];
                  const auto sqrwork2 = work[n_out + j][varID][i + offset] * work[n_out + j][varID][i + offset];
                  out[j][varID][i + offset] = std::sqrt(sqrwork1 + sqrwork2) * 2 / nts;
                }
            }
        }
    }

  if (2 * n_out == n)
    {
      for (int varID = 0; varID < nvars; ++varID)
        {
          const auto gridsize = varList2[varID].gridsize;
          const auto nlevels = varList2[varID].nlevels;
          for (int levelID = 0; levelID < nlevels; ++levelID)
            {
              const auto offset = gridsize * levelID;
              for (size_t i = 0; i < gridsize; ++i)
                out[n_out - 1][varID][i + offset] = work[2 * n_out - 1][varID][i + offset] / nts;
            }
        }
    }

  auto nout = n_out;

  taxisDefVdatetime(taxisID2, vDateTime);
  for (int j = 0; j < nout; ++j)
    {
      const auto streamID2 = streamIDs[j];
      cdo_def_timestep(streamID2, 0);

      for (int varID = 0; varID < nvars; ++varID)
        {
          const auto gridsize = varList2[varID].gridsize;
          const auto nlevels = varList2[varID].nlevels;
          for (int levelID = 0; levelID < nlevels; ++levelID)
            {
              const auto offset = gridsize * levelID;
              cdo_def_record(streamID2, varID, levelID);
              cdo_write_record(streamID2, &out[j][varID][offset], 0);
            }
        }
    }

  for (int j = 0; j < n_out && 2 * (j + 1) < n; ++j)
    {
      for (int varID = 0; varID < nvars; ++varID)
        {
          const auto gridsize = varList2[varID].gridsize;
          const auto nlevels = varList2[varID].nlevels;
          const auto missval = vlistInqVarMissval(vlistID2, varID);
          for (int levelID = 0; levelID < nlevels; ++levelID)
            {
              const auto offset = gridsize * levelID;
              for (size_t i = 0; i < gridsize; ++i)
                {
                  const auto work1 = work[j][varID][i + offset];
                  const auto work2 = work[n_out + j][varID][i + offset];
                  const auto tmpatan2 = std::atan2(work1, work2) * n / (j + 1) / 2 / M_PI;
                  out[j][varID][i + offset] = (work1 > 0.0 || work2 > 0.0) ? tmpatan2 : missval;
                  if (out[j][varID][i + offset] < 0) out[j][varID][i + offset] += n / (j + 1.);
                }
            }
        }
    }

  nout = n_out;
  if (2 * n_out == n) nout -= 1;

  taxisDefVdatetime(taxisID2, vDateTime);
  for (int j = 0; j < nout; ++j)
    {
      const auto streamID2 = streamIDs[j];
      cdo_def_timestep(streamID2, 1);

      for (int varID = 0; varID < nvars; ++varID)
        {
          const auto gridsize = varList2[varID].gridsize;
          const auto nlevels = varList2[varID].nlevels;
          const auto missval = vlistInqVarMissval(vlistID2, varID);
          for (int levelID = 0; levelID < nlevels; ++levelID)
            {
              const auto offset = gridsize * levelID;
              const auto nmiss = array_num_mv(gridsize, &out[j][varID][offset], missval);
              cdo_def_record(streamID2, varID, levelID);
              cdo_write_record(streamID2, &out[j][varID][offset], nmiss);
            }
        }
    }

  for (int j = 0; j < n_out; ++j) cdo_stream_close(streamIDs[j]);

  cdo_finish();

  return nullptr;
}
