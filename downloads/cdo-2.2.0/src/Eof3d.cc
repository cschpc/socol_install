/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Cedrick Ansorge
          Uwe Schulzweida

*/

/*
   This module contains the following operators:

     EOF3d        eof3d             3D-EOF in spatial or time space
     EOF3d        eof3dspatial      3D-EOF in spatial space
     EOF3d        eof3dtime         3D-EOF in time space
*/
/*
 * TODO:
 * Role of the weights for eofs. Should not be mixed up with division with
 * number of contributing values during summation.
 */

#ifdef _OPENMP
#include <omp.h>
#endif

#include "cdi.h"
#include "julian_date.h"

#include "cdo_options.h"
#include "cdo_vlist.h"
#include "process_int.h"
#include "param_conversion.h"
#include <mpim_grid.h>
#include "eigen_solution.h"
#include "timer.h"
#include "datetime.h"
#include "eof_mode.h"

// NO MISSING VALUE SUPPORT ADDED SO FAR

void *
EOF3d(void *process)
{
  enum
  {
    EOF3D_,
    EOF3D_TIME,
    EOF3D_SPATIAL
  };

  size_t temp_size = 0, npack = 0;
  int varID, levelID;
  bool missval_warning = false;
  size_t nmiss;
  int ngrids;
  int timer_cov = 0, timer_eig = 0;

  int calendar = CALENDAR_STANDARD;

  double sum_w;

  if (Options::Timer)
    {
      timer_cov = timer_new("Timeof cov");
      timer_eig = timer_new("Timeof eig");
    }

  cdo_initialize(process);

  // clang-format off
  cdo_operator_add("eof3d",        EOF3D_,        0, nullptr);
  cdo_operator_add("eof3dtime",    EOF3D_TIME,    0, nullptr);
  cdo_operator_add("eof3dspatial", EOF3D_SPATIAL, 0, nullptr);
  // clang-format on

  const auto operatorID = cdo_operator_id();
  const auto operfunc = cdo_operator_f1(operatorID);

  operator_input_arg("Number of eigen functions to write out");
  auto n_eig = parameter_to_int(cdo_operator_argv(0));

  auto eigen_mode = get_eigenmode();
  auto weight_mode = get_weightmode();

  // eigenvalues

  if (operfunc == EOF3D_SPATIAL) cdo_abort("Operator not Implemented - use eof3d or eof3dtime instead");

  auto streamID1 = cdo_open_read(0);
  auto vlistID1 = cdo_stream_inq_vlist(streamID1);

  VarList varList1;
  varListInit(varList1, vlistID1);

  // COUNT NUMBER OF TIMESTEPS if EOF3D_ or EOF3D_TIME
  auto nts = vlistNtsteps(vlistID1);
  if (nts == -1)
    {
      nts = 0;
      while (cdo_stream_inq_timestep(streamID1, nts)) nts++;

      if (Options::cdoVerbose) cdo_print("Counted %d timeSteps", nts);

      cdo_stream_close(streamID1);

      streamID1 = cdo_open_read(0);
      vlistID1 = cdo_stream_inq_vlist(streamID1);
    }
  else if (Options::cdoVerbose)
    cdo_print("Found %d timeSteps", nts);

  auto taxisID1 = vlistInqTaxis(vlistID1);

  // reset the requested number of eigen-function to the maximum if neccessary
  if (n_eig > nts)
    {
      cdo_warning("Solving in time-space:");
      cdo_warning("Number of eigen-functions to write out is bigger than number of time-steps.");
      cdo_warning("Setting n_eig to %d.", nts);
      n_eig = nts;
    }

  size_t n = nts;

  if (Options::cdoVerbose) cdo_print("counted %d timesteps", n);

  auto nvars = vlistNvars(vlistID1);
  int nrecs;

  auto gridID1 = varList1[0].gridID;
  auto gridsizemax = vlistGridsizeMax(vlistID1);

  // allocation of temporary fields and output structures
  Varray<double> in(gridsizemax);
  Varray2D<int> datacounts(nvars);
  Varray3D<double> datafields(nvars);
  Varray3D<double> eigenvectors(nvars);
  Varray3D<double> eigenvalues(nvars);

  size_t maxlevs = 0;
  for (varID = 0; varID < nvars; ++varID)
    {
      const auto gridsize = vlistGridsizeMax(vlistID1);
      const size_t nlevs = varList1[varID].nlevels;
      temp_size = gridsize * nlevs;
      const auto missval = varList1[varID].missval;

      if (nlevs > maxlevs) maxlevs = nlevs;

      datafields[varID].resize(nts);
      for (int tsID = 0; tsID < nts; ++tsID) datafields[varID][tsID].resize(temp_size, 0.0);

      datacounts[varID].resize(temp_size, 0);

      eigenvectors[varID].resize(n_eig);
      eigenvalues[varID].resize(nts);

      for (size_t i = 0; i < n; ++i)
        {
          if (i < (size_t) n_eig) eigenvectors[varID][i].resize(temp_size, missval);

          eigenvalues[varID][i].resize(1, missval);
        }
    }

  if (Options::cdoVerbose)
    cdo_print("Allocated eigenvalue/eigenvector with nts=%d, n=%d, gridsize=%zu for processing in %s", nts, n, gridsizemax,
              "time_space");

  Varray<double> weights(maxlevs * gridsizemax, 1.0);

  if (weight_mode == WEIGHT_ON)
    {
      const auto wstatus = gridcell_weights(gridID1, weights);
      if (wstatus != 0)
        {
          weight_mode = WEIGHT_OFF;
          cdo_warning("Using constant grid cell area weights!");
        }
      else
        {
          for (size_t k = 1; k < maxlevs; ++k)
            for (size_t i = 0; i < gridsizemax; ++i) weights[k * gridsizemax + i] = weights[i];
        }
    }

  int tsID = 0;

  // read the data and create covariance matrices for each var & level
  while (true)
    {
      nrecs = cdo_stream_inq_timestep(streamID1, tsID);
      if (nrecs == 0) break;

      for (int recID = 0; recID < nrecs; ++recID)
        {
          cdo_inq_record(streamID1, &varID, &levelID);

          const auto gridsize = varList1[varID].gridsize;
          const auto missval = varList1[varID].missval;

          cdo_read_record(streamID1, in.data(), &nmiss);

          const auto offset = gridsize * levelID;
          for (size_t i = 0; i < gridsize; ++i)
            {
              if (!DBL_IS_EQUAL(in[i], missval))
                {
                  datafields[varID][tsID][offset + i] = in[i];
                  datacounts[varID][offset + i]++;
                }
              else
                {
                  if (datacounts[varID][offset + i] != 0) cdo_abort("Missing values unsupported!");
                  if (!missval_warning)
                    {
                      // cdo_warning("Missing Value Support not checked for this Operator!");
                      // cdo_warning("Does not work with changing locations of missing values in time.");
                      missval_warning = true;
                    }
                  datafields[varID][tsID][i + offset] = 0;
                }
            }
        }
      tsID++;
    }

  if (Options::cdoVerbose) cdo_print("Read data for %d variables", nvars);

  Varray<size_t> pack(temp_size);  // TODO

  for (varID = 0; varID < nvars; ++varID)
    {
      const auto gridsize = varList1[varID].gridsize;
      const auto nlevs = varList1[varID].nlevels;
      auto missval = varList1[varID].missval;
      temp_size = gridsize * nlevs;

      if (Options::cdoVerbose)
        {
          cdo_print("============================================================================");
          cdo_print("Calculating covariance matrix and SVD for var%d (%s)", varID + 1, varList1[varID].name);
        }

      npack = 0;  // TODO already set to 0

      if (Options::Timer) timer_start(timer_cov);

      for (size_t i = 0; i < temp_size; ++i)
        {
          if (datacounts[varID][i] > 1)
            {
              pack[npack] = i;
              npack++;
            }
        }

      sum_w = 1;
      if (weight_mode == WEIGHT_ON)
        {
          sum_w = 0;
          for (size_t i = 0; i < npack; ++i) sum_w += weights[pack[i]];
        }

      if (npack < 1)
        {
          cdo_warning("Refusing to calculate EOF from a single time step for var%d (%s)", varID + 1, varList1[varID].name);
          continue;
        }

      Varray<double> eigv(n);
      Varray2D<double> covar(nts);
      for (int j1 = 0; j1 < nts; ++j1) covar[j1].resize(nts);

      if (Options::cdoVerbose)
        {
          cdo_print("varID %d allocated eigv and cov with nts=%d and n=%d", varID + 1, nts, n);
          cdo_print("   npack=%zu, nts=%d temp_size=%zu", npack, nts, temp_size);
        }

      const auto &data = datafields[varID];
#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
      for (int j1 = 0; j1 < nts; ++j1)
        {
          const auto &df1p = data[j1];
          for (int j2 = j1; j2 < nts; ++j2)
            {
              const auto &df2p = data[j2];
              double sum = 0.0;
              for (size_t i = 0; i < npack; ++i) sum += weights[pack[i] % gridsizemax] * df1p[pack[i]] * df2p[pack[i]];
              covar[j2][j1] = covar[j1][j2] = sum / sum_w / nts;
            }
        }

      if (Options::cdoVerbose) cdo_print("calculated cov-matrix");

      // SOLVE THE EIGEN PROBLEM
      if (Options::Timer) timer_stop(timer_cov);

      if (Options::Timer) timer_start(timer_eig);

      if (Options::cdoVerbose) cdo_print("Processed correlation matrix for var %d | npack: %zu", varID + 1, n);

      if (eigen_mode == JACOBI)
        parallel_eigen_solution_of_symmetric_matrix(covar, eigv, n, __func__);
      else
        eigen_solution_of_symmetric_matrix(covar, eigv, n, __func__);
      // NOW: covar contains the eigenvectors, eigv the eigenvalues

      if (Options::cdoVerbose) cdo_print("Processed SVD decomposition for var %d from %dx%d matrix", varID + 1, n, n);

      for (size_t eofID = 0; eofID < n; eofID++) eigenvalues[varID][eofID][0] = eigv[eofID];

      if (Options::Timer) timer_stop(timer_eig);

      for (int eofID = 0; eofID < n_eig; eofID++)
        {
          double *eigenvec = eigenvectors[varID][eofID].data();

#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
          for (size_t i = 0; i < npack; ++i)
            {
              double sum = 0.0;
              for (int j = 0; j < nts; ++j) sum += datafields[varID][j][pack[i]] * covar[eofID][j];

              eigenvec[pack[i]] = sum;
            }

          // NORMALIZING
          double sum = 0.0;

#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static) reduction(+ : sum)
#endif
          for (size_t i = 0; i < npack; ++i) sum += weights[pack[i] % gridsizemax] * eigenvec[pack[i]] * eigenvec[pack[i]];

          if (sum > 0)
            {
              sum = std::sqrt(sum);
#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
              for (size_t i = 0; i < npack; ++i) eigenvec[pack[i]] /= sum;
            }
          else
            {
#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
              for (size_t i = 0; i < npack; ++i) eigenvec[pack[i]] = missval;
            }
        }  // for ( eofID = 0; eofID < n_eig; eofID++ )
    }      // for ( varID = 0; varID < nvars; varID++ )

  // write files with eigenvalues (ID3) and eigenvectors (ID2)

  // eigenvalues
  const auto streamID2 = cdo_open_write(1);

  const auto vlistID2 = vlistDuplicate(vlistID1);
  const auto taxisID2 = taxisDuplicate(taxisID1);
  taxisDefRdate(taxisID2, 0);
  taxisDefRtime(taxisID2, 0);
  vlistDefTaxis(vlistID2, taxisID2);

  const auto gridID2 = gridCreate(GRID_LONLAT, 1);
  gridDefXsize(gridID2, 1);
  gridDefYsize(gridID2, 1);
  double xvals = 0.0, yvals = 0.0;
  gridDefXvals(gridID2, &xvals);
  gridDefYvals(gridID2, &yvals);

  ngrids = vlistNgrids(vlistID2);
  for (int i = 0; i < ngrids; ++i) vlistChangeGridIndex(vlistID2, i, gridID2);

  const auto zaxisID2 = zaxisCreate(ZAXIS_GENERIC, 1);
  double zvals = 0;
  zaxisDefLevels(zaxisID2, &zvals);
  cdiDefKeyString(zaxisID2, CDI_GLOBAL, CDI_KEY_NAME, "zaxis_Reduced");
  cdiDefKeyString(zaxisID2, CDI_GLOBAL, CDI_KEY_LONGNAME, "Reduced zaxis from EOF3D - only one eigen value per 3D eigen vector");

  const auto nzaxis = vlistNzaxis(vlistID2);
  for (int i = 0; i < nzaxis; ++i) vlistChangeZaxisIndex(vlistID2, i, zaxisID2);

  // eigenvectors
  const auto streamID3 = cdo_open_write(2);

  const auto vlistID3 = vlistDuplicate(vlistID1);
  const auto taxisID3 = taxisDuplicate(taxisID1);
  taxisDefRdate(taxisID3, 0);
  taxisDefRtime(taxisID3, 0);
  vlistDefTaxis(vlistID3, taxisID3);

  cdo_def_vlist(streamID2, vlistID2);
  cdo_def_vlist(streamID3, vlistID3);

  auto julianDate = julianDate_encode(calendar, cdiDateTime_set(10101, 0));

  for (tsID = 0; tsID < (int) n; ++tsID)
    {
      julianDate = julianDate_add_seconds(julianDate, 60);
      const auto vDateTime = julianDate_decode(calendar, julianDate);

      taxisDefVdatetime(taxisID2, vDateTime);
      cdo_def_timestep(streamID2, tsID);

      if (tsID < n_eig)
        {
          taxisDefVdatetime(taxisID3, vDateTime);
          cdo_def_timestep(streamID3, tsID);
        }

      for (varID = 0; varID < nvars; ++varID)
        {
          const auto missval = varList1[varID].missval;
          const auto nlevs = varList1[varID].nlevels;
          for (levelID = 0; levelID < (int) nlevs; ++levelID)
            {
              const auto offset = levelID * gridsizemax;
              if (tsID < n_eig)
                {
                  nmiss = array_num_mv(gridsizemax, &eigenvectors[varID][tsID][offset], missval);
                  cdo_def_record(streamID3, varID, levelID);
                  cdo_write_record(streamID3, &eigenvectors[varID][tsID][offset], nmiss);
                }
            }

          nmiss = (DBL_IS_EQUAL(eigenvalues[varID][tsID][0], missval)) ? 1 : 0;

          cdo_def_record(streamID2, varID, 0);
          cdo_write_record(streamID2, eigenvalues[varID][tsID].data(), nmiss);
        }  // for ( varID = 0; ... )
    }      // for ( tsID = 0; ... )

  cdo_stream_close(streamID3);
  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  cdo_finish();

  return nullptr;
}
