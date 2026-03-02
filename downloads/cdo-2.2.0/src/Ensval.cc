/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Cedrick Ansorge

*/

/*
   This module contains the following operators:
   Ensval       enscrps          Ensemble cumulative ranked probability score & decomposition
   Ensval       ensbrs           Ensemble Brier score & decomposition

   The implementation of the decomposition and score calculation as carried out in this routine follows the paper
     Hans Hersbach (2000): Decomposition of the Continuous Ranked Probability
     Score for Ensemble Prediction Systems, in: Weather and Forecasting (15) pp. 559-570
*/

#include <algorithm>  // sort

#include <cdi.h>

#include "process_int.h"
#include "cdo_vlist.h"
#include "cdo_options.h"
#include "param_conversion.h"
#include "util_files.h"

enum OPERTYPE
{
  CRPS,
  BRS
};

enum RESTYPE_BRS
{
  BRS_RES,
  BRS_RELI,
  BRS_RESOL,
  BRS_UNCTY
};

enum RESTYPE_CRPS
{
  CRPS_RES,
  CRPS_RELI,
  CRPS_POT
};

void *
Ensval(void *process)
{
  int k;
  int nrecs = 0, nostreams = 0, ngrids;
  size_t nmiss;
  int levelID = -1, varID = -1;
  size_t gridsize = 0;
  int vlistID;
  int have_miss = 0;
  int stream = 0;
  CdoStreamID streamID = 0;
  double missval = 0;
  double sum_weights = 0;
  double crps_reli = 0, crps_pot = 0, crps = 0;
  double heavyside0, heavysideN;
  double brs_reli, brs_resol, brs_uncty, brs_thresh = 0;

  int fileID;
  char type_suffix[10];

  struct ens_file_t
  {
    CdoStreamID streamID;
    int vlistID;
    Varray<double> array;
  };

  Varray<double> alpha, beta, alpha_weights, beta_weights;
  Varray<double> brs_g, brs_o;
  Varray<double> weights;
  // int vlistCheck, gridsizeCheck;

  cdo_initialize(process);

  cdo_operator_add("enscrps", CRPS, 0, nullptr);
  cdo_operator_add("ensbrs", BRS, 0, nullptr);

  const auto operatorID = cdo_operator_id();
  const auto operfunc = cdo_operator_f1(operatorID);

  const auto nfiles = cdo_stream_cnt();
  const auto nens = nfiles - 1;

  if (operfunc == CRPS) { nostreams = 3; }
  else if (operfunc == BRS)
    {
      operator_input_arg("Threshold for Brier score?");
      operator_check_argc(1);
      brs_thresh = parameter_to_double(cdo_operator_argv(0));
      nostreams = 4;

      fprintf(stderr, "brs_thres %10.6f\n", brs_thresh);
    }

  // allocate array to hold results
  Varray<double> r(nostreams);

  // one stream for each value of the decomposition
  std::vector<CdoStreamID> streamID2(nostreams);
  std::vector<int> vlistID2(nostreams), taxisID2(nostreams), zaxisID2(nostreams);

  Varray<double> val(nens, 0);

  if (operfunc == CRPS)
    {
      alpha.resize(nens + 1, 0);
      beta.resize(nens + 1, 0);
      alpha_weights.resize(nens + 1, 0);
      beta_weights.resize(nens + 1, 0);
    }
  else if (operfunc == BRS)
    {
      brs_g.resize(nens + 1, 0);
      brs_o.resize(nens + 1, 0);
    }

  if (Options::cdoVerbose) cdo_print("Ensemble over %d files.", nens);

  std::vector<ens_file_t> ef(nfiles);

  for (fileID = 0; fileID < nfiles; ++fileID)
    {
      streamID = cdo_open_read(fileID);
      vlistID = cdo_stream_inq_vlist(streamID);

      ef[fileID].streamID = streamID;
      ef[fileID].vlistID = vlistID;
    }

  const auto streamID1 = ef[0].streamID;

  if (Options::cdoVerbose) cdo_print("Opened %i Input Files for Ensemble Operator", nfiles);

  // check for identical contents of all ensemble members
  const auto nvars = vlistNvars(ef[0].vlistID);
  if (Options::cdoVerbose) cdo_print("nvars %i", nvars);

  for (fileID = 1; fileID < nfiles; ++fileID) vlist_compare(ef[0].vlistID, ef[fileID].vlistID, CMP_ALL);

  const auto vlistID1 = ef[0].vlistID;
  const auto taxisID1 = vlistInqTaxis(vlistID1);
  const auto zaxisID1 = vlistInqVarZaxis(vlistID1, 0);

  const auto gridID2 = gridCreate(GRID_LONLAT, 1);
  gridDefXsize(gridID2, 1);
  gridDefYsize(gridID2, 1);
  double xval = 0, yval = 0;
  gridDefXvals(gridID2, &xval);
  gridDefYvals(gridID2, &yval);

  auto ofilebase = cdo_get_obase();

  auto refname = cdo_get_stream_name(0);
  char filesuffix[32] = { 0 };
  FileUtils::gen_suffix(filesuffix, sizeof(filesuffix), cdo_inq_filetype(streamID1), vlistID1, refname);

  for (stream = 0; stream < nostreams; stream++)
    {
      int namelen = ofilebase.size() + 9  // type_suffix
                    + 32                  // filesuffix
                    + 3;                  // separating dots and EOS

      switch (operfunc)
        {
        case CRPS:
          switch (stream)
            {
            case 0: sprintf(type_suffix, "crps"); break;
            case 1: sprintf(type_suffix, "crps_reli"); break;
            case 2: sprintf(type_suffix, "crps_pot"); break;
            }
          break;
        case BRS:
          switch (stream)
            {
            case 0: sprintf(type_suffix, "brs"); break;
            case 1: sprintf(type_suffix, "brs_reli"); break;
            case 2: sprintf(type_suffix, "brs_reso"); break;
            case 3: sprintf(type_suffix, "brs_unct"); break;
            }
          break;
        }

      std::vector<char> ofilename(namelen, 0);

      std::snprintf(ofilename.data(), namelen, "%s.%s%s", ofilebase.c_str(), type_suffix, filesuffix);
      // fprintf(stderr, "StreamID %i: %s\n", stream, ofilename);

      if (!Options::cdoOverwriteMode && FileUtils::file_exists(ofilename.data()) && !FileUtils::user_file_overwrite(ofilename.data()))
        cdo_abort("Outputfile %s already exists!", ofilename.data());

      streamID2[stream] = cdo_open_write(ofilename.data());

      zaxisID2[stream] = zaxisDuplicate(zaxisID1);
      taxisID2[stream] = taxisDuplicate(taxisID1);
      vlistID2[stream] = vlistDuplicate(vlistID1);

      ngrids = vlistNgrids(vlistID2[stream]);
      for (int i = 0; i < ngrids; ++i) vlistChangeGridIndex(vlistID2[stream], i, gridID2);

      vlistDefTaxis(vlistID2[stream], taxisID2[stream]);
      cdo_def_vlist(streamID2[stream], vlistID2[stream]);
    }

  if (Options::cdoVerbose) cdo_print(" sum_weights %10.6f", sum_weights);

  VarList varList1;
  varListInit(varList1, vlistID1);

  int tsID = 0;
  do {
      const auto nrecs0 = cdo_stream_inq_timestep(ef[0].streamID, tsID);
      for (fileID = 1; fileID < nfiles; ++fileID)
        {
          streamID = ef[fileID].streamID;
          nrecs = cdo_stream_inq_timestep(streamID, tsID);
          if (nrecs != nrecs0)
            cdo_abort("Number of records at time step %d of %s and %s differ!", tsID + 1, cdo_get_stream_name(0),
                      cdo_get_stream_name(fileID));
        }

      for (stream = 0; stream < nostreams; stream++)
        {
          cdo_taxis_copy_timestep(taxisID2[stream], taxisID1);
          if (nrecs0 > 0) cdo_def_timestep(streamID2[stream], tsID);
        }

      for (int recID = 0; recID < nrecs0; ++recID)
        {
          for (fileID = 0; fileID < nfiles; ++fileID)
            {
              streamID = ef[fileID].streamID;
              cdo_inq_record(streamID, &varID, &levelID);

              if (fileID == 0)
                {
                  gridsize = varList1[varID].gridsize;
                  missval = varList1[varID].missval;
                  weights.resize(gridsize);
                }

              ef[fileID].array.resize(gridsize);

              cdo_read_record(streamID, ef[fileID].array.data(), &nmiss);
            }

          // xsize = gridInqXsize(gridID);
          // ysize = gridInqYsize(gridID);

          /*
          if (xsize > 1 && ysize > 1)
          {
            gridcell_weights(gridID, weights);
            sum_weights = varray_sum(gridsize, weights);
          }
          else
          */
          {
            varray_fill(weights, 1.0 / gridsize);
            sum_weights = 1.0;
          }

          nmiss = 0;
          heavyside0 = 0;
          heavysideN = 0;

          for (size_t i = 0; i < gridsize; ++i)
            {
              double xa = 0.0;
              have_miss = 0;
              for (fileID = 0; fileID < nfiles; ++fileID)
                {
                  auto value = ef[fileID].array[i];
                  if (fileID == 0)
                    xa = value;  // 1st file contains reference
                  else
                    val[fileID - 1] = value;  // Ensembles start at 2nd file

                  if (dbl_is_equal(value, missval))
                    {
                      have_miss = 1;
                      break;
                    }
                }

              auto &x = val;
              std::sort(x.begin(), x.end());  // Sort The Ensemble Array to ascending order

              // only process if no missing value in ensemble
              if (!have_miss && operfunc == CRPS)
                {
                  if (xa < x[0])
                    { /* Consider outliers            */
                      beta[0] += (x[0] - xa) * weights[i];
                      heavyside0 += 1.;
                    }
                  if (xa > x[nens - 1])
                    {
                      alpha[nens] += (xa - x[nens - 1]) * weights[i];
                      alpha_weights[nens] += weights[i];
                      heavysideN += 1.;
                    }

                  // Loop start at zero ==> 1st ensemble (c-indexing)
                  for (k = 0; k < nens - 1; ++k)
                    {                     // Cumulate alpha and beta
                      if (xa > x[k + 1])  // left of heavyside
                        alpha[k + 1] += (x[k + 1] - x[k]) * weights[i];
                      else if (xa < x[k])  // right of heavyside
                        beta[k + 1] += (x[k + 1] - x[k]) * weights[i];
                      else if (x[k + 1] >= xa && xa >= x[k])  // hitting jump pf heavyside (occurs exactly once!)
                        beta[k + 1] += (x[k + 1] - xa) * weights[i];
                    }
                }
              else if (operfunc == BRS)
                {
                  // int occ = xa > brs_thresh? 1 : 0;

                  // brs_g[i] - number of enemble members with rank i that
                  // forecast event
                  //          - event: value > brs_thresh
                  //
                  if (x[0] > brs_thresh)
                    brs_g[0] += weights[i];
                  else if (x[nens - 1] < brs_thresh)
                    brs_g[nens] += weights[i];
                  else
                    for (k = 0; k < nens - 1; ++k)
                      {
                        if (x[k + 1] >= brs_thresh && brs_thresh >= x[k])
                          {
                            brs_g[k + 1] += weights[i];
                            break;
                          }
                      }

                  // brs_o[i] - number of times that the obs is between Ensemble i-1 and i
                  if (1)
                    {
                      if (x[0] > xa)
                        brs_o[0] += weights[i];
                      else if (x[nens - 1] < xa)
                        brs_o[nens] += weights[i];
                      else
                        for (k = 0; k < nens - 1; ++k)
                          {
                            if (x[k + 1] >= xa && xa >= x[k])
                              {
                                brs_o[k + 1] += weights[i];
                                break;
                              }
                          }
                    }
                }
            }  // for ( i=0; i<gridsize; i++ )

          if (operfunc == CRPS)
            {
              // First Bin
              double p = 0.0, g = 0.0;
              double o = heavyside0 / gridsize;
              if (o > 0.) g = beta[0] / o;
              crps_reli = g * (o - p) * (o - p);
              crps_pot = g * o * (1. - o);
              crps = g * ((1. - o) * p * p + o * (1. - p) * (1. - p));

              // Middle Bins
              for (k = 1; k < nens; ++k)
                {
                  p = (double) k / (double) nens;

                  if (!dbl_is_equal(sum_weights, 1.0))
                    {
                      alpha[k] /= sum_weights;
                      beta[k] /= sum_weights;
                    }

                  g = alpha[k] + beta[k];
                  o = beta[k] / (alpha[k] + beta[k]);

                  crps_reli += g * (o - p) * (o - p);
                  crps_pot += g * o * (1. - o);
                  crps += g * ((1. - o) * p * p + o * (1. - p) * (1. - p));
                }

              // Last Bin
              p = 1.0;
              g = 0.0;
              o = 1.0 - heavysideN / gridsize;
              if (is_not_equal(o, 1.))
                {
                  g = alpha[nens] / (1 - o);

                  crps_reli += g * (o - p) * (o - p);
                  crps_pot += g * o * (1 - o);
                  crps += g * ((1 - o) * p * p + o * (1 - p) * (1 - p));
                }
              r[CRPS_RES] = crps;
              r[CRPS_RELI] = crps_reli;
              r[CRPS_POT] = crps_pot;
            }
          else if (operfunc == BRS)
            {
              brs_reli = 0;
              brs_resol = 0;
              brs_uncty = 0;

              double gsum = 0, obar = 0, osum = 0;
              for (k = 0; k <= nens; ++k)
                {
                  obar += brs_g[k] * brs_o[k];
                  gsum += brs_g[k];
                  osum += brs_o[k];
                }

              if (std::fabs(osum - 1) > 1.e-06 || std::fabs(gsum - 1) > 1.e-06)
                {
                  cdo_abort("Internal error - normalization constraint of problem not fulfilled");
                  cdo_abort("This is likely due to missing values");
                }

              double o = 0, p = 0, g = 0;
              brs_uncty = obar * (1 - obar);

              for (k = 0; k <= nens; ++k)
                {

                  g = brs_g[k];
                  o = brs_o[k];
                  p = 1.0 - k / (float) nens;
                  // need p = 1 - k/nens here as k=0 if all members forecast
                  // event and k=nens if none does so.

                  brs_reli += g * (o - p) * (o - p);
                  brs_resol += g * (o - obar) * (o - obar);
                }

              r[BRS_RES] = brs_reli - brs_resol + brs_uncty;
              r[BRS_RELI] = brs_reli;
              r[BRS_RESOL] = brs_resol;
              r[BRS_UNCTY] = brs_uncty;

              if (Options::cdoVerbose)
                {
                  // cdo_print("Brier score for var %i level %i calculated",varID, levelID);
                  cdo_print("BRS: obar %12.6g brs  %12.6g reli %12.6g resol %12.6g u %12.6g", obar,
                            brs_reli - brs_resol + brs_uncty, brs_reli, brs_resol, brs_uncty);
                }
            }

          if (Options::cdoVerbose && operfunc == CRPS)
            cdo_print("CRPS:%12.6g reli:%12.6g crps_pot:%12.6g crps:%12.6g", crps, crps_reli, crps_pot, crps_reli + crps_pot);

          for (stream = 0; stream < nostreams; stream++)
            {
              cdo_def_record(streamID2[stream], varID, levelID);
              if (std::isnan(r[stream]))
                {
                  r[stream] = missval;
                  have_miss = 1;
                }
              cdo_write_record(streamID2[stream], &r[stream], have_miss);
            }

          switch (operfunc)
            {
            case (CRPS):
              varray_fill(alpha, 0.0);
              varray_fill(beta, 0.0);
              heavyside0 = 0;
              heavysideN = 0;
              break;
            case (BRS):
              varray_fill(brs_o, 0.0);
              varray_fill(brs_g, 0.0);
              break;
            }
        }  // for ( int recID = 0; recID < nrecs; recID++ )
      tsID++;
    }
  while (nrecs);

  for (fileID = 0; fileID < nfiles; ++fileID)
    {
      streamID = ef[fileID].streamID;
      cdo_stream_close(streamID);
    }

  for (stream = 0; stream < nostreams; stream++) cdo_stream_close(streamID2[stream]);

  for (stream = 0; stream < nostreams; stream++)
    {
      vlistDestroy(vlistID2[stream]);
      taxisDestroy(taxisID2[stream]);
      zaxisDestroy(zaxisID2[stream]);
    }

  gridDestroy(gridID2);

  cdo_finish();

  return nullptr;
}
