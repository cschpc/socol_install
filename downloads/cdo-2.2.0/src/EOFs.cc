/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Cedrick Ansorge
          Uwe Schulzweida

*/

/*
   This module contains the following operators:

     Timeof        eof             EOF in spatial or time space
     Timeof        eofspatial      EOF in spatial space
     Timeof        eoftime         EOF in time space
*/
/*
 * TODO:
 * Role of the weights for eofs. Should not be mixed up with division with number of contributing values during summation.
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

// No missing value support added so far!

static void
scale_eigvec_grid(Varray<double> &out, int tsID, size_t npack, const std::vector<size_t> &pack, const Varray<double> &weight,
                  const Varray2D<double> &covar, double sum_w)
{
  for (size_t i = 0; i < npack; ++i) out[pack[i]] = covar[tsID][i] / std::sqrt(weight[pack[i]] / sum_w);
}

static void
scale_eigvec_time(Varray<double> &out, int tsID, int nts, size_t npack, const std::vector<size_t> &pack,
                  const Varray<double> &weight, const Varray2D<double> &covar, const Varray2D<double> &data, double missval,
                  double sum_w)
{
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(npack, nts, tsID, pack, data, covar, out)
#endif
  for (size_t i = 0; i < npack; ++i)
    {
      double sum = 0.0;
      for (int j = 0; j < nts; ++j) sum += data[j][i] * covar[tsID][j];

      out[pack[i]] = sum;
    }
  /*
  for ( size_t j = 0; j < nts; ++j )
    {
      for ( size_t i = 0; i < npack; ++i )
        out[pack[i]] += data[j][i] * covar[tsID][j];
    }
  */

  // Normalizing
  double sum = 0.0;

#ifdef _OPENMP
#pragma omp parallel for default(none) reduction(+ : sum) shared(out, weight, pack, npack)
#endif
  for (size_t i = 0; i < npack; ++i)
    {
      // do not need to account for weights as eigenvectors are non-weighted
      sum += weight[pack[i]] * out[pack[i]] * out[pack[i]];
    }

  if (sum > 0.0)
    {
      sum = std::sqrt(sum / sum_w);
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(npack, pack, sum, out)
#endif
      for (size_t i = 0; i < npack; ++i) out[pack[i]] /= sum;
    }
  else
    {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(npack, pack, out, missval)
#endif
      for (size_t i = 0; i < npack; ++i) out[pack[i]] = missval;
    }
}

void *
EOFs(void *process)
{
  enum
  {
    EOF_,
    EOF_TIME,
    EOF_SPATIAL
  };

  size_t nmiss;
  int varID, levelID;
  int nts = 1;
  int grid_space = 0, time_space = 0;
  int timer_cov = 0, timer_eig = 0;

  auto calendar = CALENDAR_STANDARD;

  struct eofdata_t
  {
    bool init = false;
    bool first_call = true;
    Varray<double> eig_val;
    Varray2D<double> covar;
    Varray2D<double> data;
  };

  if (Options::Timer)
    {
      timer_cov = timer_new("Timeof cov");
      timer_eig = timer_new("Timeof eig");
    }

  cdo_initialize(process);

  // clang-format off
  cdo_operator_add("eof",        EOF_,        0, nullptr);
  cdo_operator_add("eoftime",    EOF_TIME,    0, nullptr);
  cdo_operator_add("eofspatial", EOF_SPATIAL, 0, nullptr);
  // clang-format on

  auto operatorID = cdo_operator_id();
  auto operfunc = cdo_operator_f1(operatorID);

  operator_input_arg("Number of eigen functions to write out");
  auto n_eig = parameter_to_int(cdo_operator_argv(0));

  auto eigen_mode = get_eigenmode();
  auto weight_mode = get_weightmode();

  auto streamID1 = cdo_open_read(0);
  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto taxisID1 = vlistInqTaxis(vlistID1);

  VarList varList1;
  varListInit(varList1, vlistID1);

  auto gridID1 = varList1[0].gridID;
  auto gridsizemax = vlistGridsizeMax(vlistID1);
  auto nvars = vlistNvars(vlistID1);

  auto ngrids = vlistNgrids(vlistID1);
  for (int index = 1; index < ngrids; ++index)
    if (vlistGrid(vlistID1, 0) != vlistGrid(vlistID1, index)) cdo_abort("Too many different grids!");

  // eigenvalues

  // Count number of timesteps if EOF_ or EOF_TIME
  if (operfunc == EOF_ || operfunc == EOF_TIME)
    {
      if (Options::cdoVerbose) cdo_print("Counting timesteps in ifile");

      nts = vlistNtsteps(vlistID1);
      if (nts == 0) nts = 1;
      if (nts == -1)
        {
          nts = 0;
          while (cdo_stream_inq_timestep(streamID1, nts)) nts++;

          if (Options::cdoVerbose) cdo_print("Counted %d timeSteps", nts);

          cdo_stream_close(streamID1);

          streamID1 = cdo_open_read(0);
          vlistID1 = cdo_stream_inq_vlist(streamID1);
          taxisID1 = vlistInqTaxis(vlistID1);
        }
      else if (Options::cdoVerbose)
        cdo_print("Found %d timeSteps", nts);

      if ((size_t) nts < gridsizemax || operfunc == EOF_TIME)
        {
          time_space = 1;
          grid_space = 0;
        }
      else
        {
          time_space = 0;
          grid_space = 1;
        }
    }
  else if (operfunc == EOF_SPATIAL)
    {
      time_space = 0;
      grid_space = 1;
    }

  // reset the requested number of eigen-function to the maximum if neccessary
  size_t n = 0;
  if (time_space)
    {
      if (n_eig > nts)
        {
          cdo_warning("Solving in time-space:");
          cdo_warning("Number of eigen-functions to write out is bigger than number of time-steps.");
          cdo_warning("Setting n_eig to %d.", nts);
          cdo_warning("If You want to force a solution in grid-space use operator eofspatial");
          n_eig = nts;
        }
      n = nts;
    }
  else if (grid_space)
    {
      if (((double) gridsizemax) * gridsizemax > (double) SIZE_MAX) cdo_abort("Grid space too large!");

      if ((size_t) n_eig > gridsizemax)
        {
          cdo_warning("Solving in spatial space");
          cdo_warning("Number of eigen-functions to write out is bigger than grid size");
          cdo_warning("Setting n_eig to %zu", gridsizemax);
          cdo_warning("If You want to force a solution in time-space use operator eoftime");
          n_eig = gridsizemax;
        }
      n = gridsizemax;
    }

  if (Options::cdoVerbose)
    cdo_print("Calculating %d eigenvectors and %zu eigenvalues in %s", n_eig, n, (grid_space == 1) ? "grid_space" : "time_space");

  Varray<double> weights(gridsizemax, 1.0);

  if (weight_mode == WEIGHT_ON)
    {
      auto wstatus = gridcell_weights(gridID1, weights);
      if (wstatus != 0)
        {
          weight_mode = WEIGHT_OFF;
          cdo_warning("Using constant grid cell area weights!");
        }
    }

  // allocation of temporary fields and output structures
  size_t npack = SIZE_MAX;
  std::vector<size_t> pack(gridsizemax);
  Varray<double> in(gridsizemax);
  std::vector<std::vector<eofdata_t>> eofdata(nvars);

  for (varID = 0; varID < nvars; ++varID)
    {
      auto nlevs = varList1[varID].nlevels;
      eofdata[varID].resize(nlevs);

      if (time_space)
        for (levelID = 0; levelID < nlevs; ++levelID) eofdata[varID][levelID].data.resize(nts);
    }

  if (Options::cdoVerbose) cdo_print("Allocated eigenvalue/eigenvector structures with nts=%d gridsize=%zu", nts, gridsizemax);

  double sum_w = 1.0;

  int tsID = 0;

  // read the data and create covariance matrices for each var & level
  while (true)
    {
      auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
      if (nrecs == 0) break;

      for (int recID = 0; recID < nrecs; ++recID)
        {
          cdo_inq_record(streamID1, &varID, &levelID);
          cdo_read_record(streamID1, in.data(), &nmiss);

          auto gridsize = varList1[varID].gridsize;
          auto missval = varList1[varID].missval;

          if (npack == SIZE_MAX)
            {
              npack = 0;
              for (size_t i = 0; i < gridsize; ++i)
                {
                  if (!DBL_IS_EQUAL(weights[i], 0.0) && !DBL_IS_EQUAL(weights[i], missval) && !DBL_IS_EQUAL(in[i], missval))
                    {
                      pack[npack] = i;
                      npack++;
                    }
                }

              if (weight_mode == WEIGHT_ON)
                {
                  sum_w = 0.0;
                  for (size_t i = 0; i < npack; ++i) sum_w += weights[pack[i]];
                }
            }

          {
            size_t ipack = 0;
            for (size_t i = 0; i < gridsize; ++i)
              {
                if (!DBL_IS_EQUAL(weights[i], 0.0) && !DBL_IS_EQUAL(weights[i], missval) && !DBL_IS_EQUAL(in[i], missval))
                  {
                    if (pack[ipack] != i) cdo_abort("Missing values unsupported!");
                    ipack++;
                  }
              }
            if (ipack != npack) cdo_abort("Missing values unsupported!");
          }

          if (grid_space)
            {
              if (!eofdata[varID][levelID].init)
                {
                  eofdata[varID][levelID].covar.resize(npack);
                  for (size_t i = 0; i < npack; ++i) eofdata[varID][levelID].covar[i].resize(npack, 0.0);
                }

              auto &covar = eofdata[varID][levelID].covar;
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(npack, covar, in, pack)
#endif
              for (size_t ipack = 0; ipack < npack; ++ipack)
                {
                  auto &covar_i = covar[ipack];
                  auto in_i = in[pack[ipack]];
                  for (size_t jpack = ipack; jpack < npack; ++jpack) covar_i[jpack] += in_i * in[pack[jpack]];
                }
            }
          else if (time_space)
            {
              eofdata[varID][levelID].data[tsID].resize(npack);
              auto &data = eofdata[varID][levelID].data[tsID];

              for (size_t ipack = 0; ipack < npack; ipack++) data[ipack] = in[pack[ipack]];
            }

          eofdata[varID][levelID].init = true;
        }

      tsID++;
    }

  if (grid_space) nts = tsID;

  if (tsID == 1) cdo_abort("File consists of only one timestep!");

  // write files with eigenvalues (ID3) and eigenvectors (ID2)

  // eigenvalues
  auto streamID2 = cdo_open_write(1);

  auto vlistID2 = vlistDuplicate(vlistID1);
  auto taxisID2 = taxisDuplicate(taxisID1);
  taxisDefRdate(taxisID2, 0);
  taxisDefRtime(taxisID2, 0);
  vlistDefTaxis(vlistID2, taxisID2);

  auto gridID2 = gridCreate(GRID_LONLAT, 1);
  gridDefXsize(gridID2, 1);
  gridDefYsize(gridID2, 1);
  double xvals = 0.0, yvals = 0.0;
  gridDefXvals(gridID2, &xvals);
  gridDefYvals(gridID2, &yvals);
  for (int i = 0; i < ngrids; ++i) vlistChangeGridIndex(vlistID2, i, gridID2);

  // eigenvectors
  auto streamID3 = cdo_open_write(2);

  auto vlistID3 = vlistDuplicate(vlistID1);
  auto taxisID3 = taxisDuplicate(taxisID1);
  taxisDefRdate(taxisID3, 0);
  taxisDefRtime(taxisID3, 0);
  vlistDefTaxis(vlistID3, taxisID3);

  cdo_def_vlist(streamID2, vlistID2);
  cdo_def_vlist(streamID3, vlistID3);

  int64_t vdate = 10101;
  int vtime = 0;
  CdiDateTime vDateTime = cdiDateTime_set(vdate, vtime);
  auto julianDate = julianDate_encode(calendar, vDateTime);

  Varray<double> &out = in;

  int nts_out = (npack < (size_t) nts) ? npack : nts;

  for (tsID = 0; tsID < nts_out; ++tsID)
    {
      julianDate = julianDate_add_seconds(julianDate, 60);
      vDateTime = julianDate_decode(calendar, julianDate);

      taxisDefVdatetime(taxisID2, vDateTime);
      cdo_def_timestep(streamID2, tsID);

      if (tsID < n_eig)
        {
          taxisDefVdatetime(taxisID3, vDateTime);
          cdo_def_timestep(streamID3, tsID);
        }

      for (varID = 0; varID < nvars; ++varID)
        {
          const auto &vname = varList1[varID].name;
          auto gridsize = varList1[varID].gridsize;
          auto nlevs = varList1[varID].nlevels;
          auto missval = varList1[varID].missval;

          for (levelID = 0; levelID < nlevs; ++levelID)
            {
              const auto &data = eofdata[varID][levelID].data;
              auto &covar = eofdata[varID][levelID].covar;

              if (eofdata[varID][levelID].first_call)
                {
                  eofdata[varID][levelID].first_call = false;

                  if (Options::cdoVerbose)
                    cdo_print("Calculating covar matrices for %d levels of var%i (%s)", nlevs, varID + 1, vname);

                  if (Options::Timer) timer_start(timer_cov);

                  if (Options::cdoVerbose) cdo_print("processing level %d", levelID + 1);

                  if (grid_space)
                    {
                      if (npack) eofdata[varID][levelID].eig_val.resize(npack);

                      for (size_t ipack = 0; ipack < npack; ++ipack)
                        {
                          size_t i = pack[ipack];
                          for (size_t jpack = 0; jpack < npack; ++jpack)
                            {
                              if (jpack < ipack) { covar[ipack][jpack] = covar[jpack][ipack]; }
                              else
                                {
                                  auto j = pack[jpack];
                                  covar[ipack][jpack] = covar[ipack][jpack] *                                    // covariance
                                                        std::sqrt(weights[i]) * std::sqrt(weights[j]) / sum_w /  // weights
                                                        nts;  // number of data contributing
                                }
                            }
                        }
                    }
                  else if (time_space)
                    {
                      if (Options::cdoVerbose) cdo_print("allocating covar with %dx%d elements | npack=%zu", nts, nts, npack);

                      covar.resize(nts);
                      for (int i = 0; i < nts; ++i) covar[i].resize(nts);

                      eofdata[varID][levelID].eig_val.resize(nts);

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(nts, data, covar, weights, npack, pack, sum_w) schedule(static)
#endif
                      for (int j1 = 0; j1 < nts; ++j1)
                        {
                          const auto &df1p = data[j1];
                          for (int j2 = j1; j2 < nts; ++j2)
                            {
                              const auto &df2p = data[j2];
                              double sum = 0.0;
                              for (size_t i = 0; i < npack; ++i) sum += weights[pack[i]] * df1p[i] * df2p[i];
                              covar[j2][j1] = covar[j1][j2] = sum / sum_w / nts;
                            }
                        }

                      if (Options::cdoVerbose) cdo_print("finished calculation of covar-matrix for var %s", vname);
                    }

                  if (Options::Timer) timer_stop(timer_cov);

                  // Solve the eigen problem
                  if (Options::Timer) timer_start(timer_eig);

                  auto &eig_val = eofdata[varID][levelID].eig_val;
                  if (eigen_mode == JACOBI)
                    // TODO: use return status (>0 okay, -1 did not converge at all)
                    parallel_eigen_solution_of_symmetric_matrix(covar, eig_val, n, __func__);
                  else
                    eigen_solution_of_symmetric_matrix(covar, eig_val, n, __func__);

                  if (Options::Timer) timer_stop(timer_eig);
                  // NOW: covar contains the eigenvectors, eig_val the eigenvalues

                  for (size_t i = 0; i < gridsize; ++i) out[i] = missval;

                  // for ( int i = 0; i < n; i++ ) eig_val[i] *= sum_w;
                }  // first_call

              if (tsID < n_eig)
                {
                  if (grid_space)
                    scale_eigvec_grid(out, tsID, npack, pack, weights, covar, sum_w);
                  else if (time_space)
                    scale_eigvec_time(out, tsID, nts, npack, pack, weights, covar, data, missval, sum_w);

                  nmiss = varray_num_mv(gridsize, out, missval);
                  cdo_def_record(streamID3, varID, levelID);
                  cdo_write_record(streamID3, out.data(), nmiss);
                }  // loop n_eig

              auto eig_val = eofdata[varID][levelID].eig_val.data();

              nmiss = (DBL_IS_EQUAL(eig_val[tsID], missval)) ? 1 : 0;

              cdo_def_record(streamID2, varID, levelID);
              cdo_write_record(streamID2, &eig_val[tsID], nmiss);
            }  // loop nlevs
        }      // loop nvars
    }

  cdo_stream_close(streamID3);
  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  vlistDestroy(vlistID2);
  vlistDestroy(vlistID3);

  gridDestroy(gridID2);

  //  taxisDestroy(taxisID2);
  //  taxisDestroy(taxisID3);

  cdo_finish();

  return nullptr;
}
