/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Cedrick Ansorge

*/

/*
   This module contains the following operators:
   Ensstat3       ensrkhistspace   Ensemble ranked histogram averaged over time
   Ensstat3       ensrkhisttime    Ensemble ranked histogram averaged over space
   Ensstat3       ensroccurve      Ensamble Receiver Operating Characteristics
*/

#include <cdi.h>

#include "process_int.h"
#include "cdo_vlist.h"
#include "param_conversion.h"
#include "cdo_options.h"
#include "util_files.h"
#include "cimdOmp.h"
#include "field_functions.h"

// Defines for rank histogram
enum TDATA_TYPE
{
  TIME,
  SPACE
};

#define time_data TIME
#define space_data SPACE

// Defines for Receiver Operating Characteristics (ROC)
#define DEBUG_ROC 0
enum CONTINGENCY_TYPE
{
  TP,  // TP - True positive  ( event     forecast and     occured)  HIT
  FP,  // FP - False positive ( event     forecast and not occured)  false ALARM
  FN,  // FN - False negative ( event not forecast and     ocurred)  MISSED
  TN   // TN - True negative  ( event not forecast and not ocurred)  CORRECT REJECTION
};

enum ROC_ENUM_TYPE
{
  TPR,  // TPR = True Positive Rate = TP / ( TP + FN )
  FPR   // FNR = False Negtive Rate = FN / ( FP + TN )
};

static double
roc_curve_integrate(const Varray2D<double> &roc, int n)
{
  double area = 0.0;

  for (int i = 1; i <= n; ++i)
    {
      auto x1 = roc[i][FPR];
      auto x0 = roc[i - 1][FPR];
      auto y1 = roc[i][TPR];
      auto y0 = roc[i - 1][TPR];
      auto dx = x1 - x0;
      auto dy = y1 - y0;

      auto step_area = -0.5 * dx * dy - dx * y0;
      area += step_area;
    }

  return area - 0.5;
}

void *
Ensstat3(void *process)
{
  int j;
  int nrecs = 0, nrecs0;
  size_t nmiss = 0;
  int chksum = 0;  // for check of histogram population
  int levelID, varID, binID = 0;
  int gridID, gridID2;
  int have_miss = 0;
  CdoStreamID streamID2 = 0;

  struct ens_file_t
  {
    CdoStreamID streamID;
    int vlistID;
    Varray<double> array;
  };

  cdo_initialize(process);

  enum
  {
    func_roc,
    func_rank
  };

  // clang-format off
  cdo_operator_add("ensroc",         func_roc,  0,          nullptr);
  cdo_operator_add("ensrkhistspace", func_rank, space_data, nullptr);
  cdo_operator_add("ensrkhisttime",  func_rank, time_data,  nullptr);
  // clang-format on

  auto operatorID = cdo_operator_id();
  auto operfunc = cdo_operator_f1(operatorID);
  auto datafunc = cdo_operator_f2(operatorID);

  int nbins = 0;
  if (operfunc == func_roc)
    {
      operator_input_arg("Number of eigen functions to write out");
      nbins = parameter_to_int(cdo_operator_argv(0));
    }

  auto nfiles = cdo_stream_cnt() - 1;
  auto nens = nfiles - 1;

  if (Options::cdoVerbose) cdo_print("Ensemble over %d files.", nfiles);

  auto ofilename = cdo_get_stream_name(nfiles);

  if (!Options::cdoOverwriteMode && FileUtils::file_exists(ofilename) && !FileUtils::user_file_overwrite(ofilename))
    cdo_abort("Outputfile %s already exists!", ofilename);

  std::vector<ens_file_t> ef(nfiles);

  /* *************************************************** */
  /* should each thread be allocating memory locally???? */
  /* ("first touch strategy")                            */
  /* --> #pragma omp parallel for ...                    */
  /* *************************************************** */
  FieldVector field(Threading::ompNumThreads);
  for (int i = 0; i < Threading::ompNumThreads; ++i)
    {
      field[i].resize(nfiles);
      field[i].weightv.resize(nfiles);
      for (int fileID = 0; fileID < nfiles; ++fileID) field[i].weightv[fileID] = 1;
    }

  for (int fileID = 0; fileID < nfiles; ++fileID)
    {
      auto streamID = cdo_open_read(fileID);
      ef[fileID].streamID = streamID;
      ef[fileID].vlistID = cdo_stream_inq_vlist(streamID);
    }

  // check for identical contents of all ensemble members
  for (int fileID = 1; fileID < nfiles; ++fileID) vlist_compare(ef[0].vlistID, ef[fileID].vlistID, CMP_ALL);

  auto vlistID1 = ef[0].vlistID;
  auto vlistID2 = vlistCreate();
  vlistDefNtsteps(vlistID2, vlistNtsteps(vlistID1));

  auto nvars = vlistNvars(vlistID1);
  std::vector<int> varIDs2(nvars);

  auto zaxisID2 = zaxisCreate(ZAXIS_GENERIC, nfiles);
  {
    Varray<double> levs(nfiles, 0);
    for (int i = 0; i < nfiles; ++i) levs[i] = i;
    zaxisDefLevels(zaxisID2, levs.data());
    cdiDefKeyString(zaxisID2, CDI_GLOBAL, CDI_KEY_NAME, "histogram_binID");
  }

  int time_mode = (datafunc == TIME) ? TIME_VARYING : TIME_CONSTANT;

  for (varID = 0; varID < nvars; ++varID)
    {
      /* **************************************************************** */
      /* nfiles includes the observation, so there are nfiles-1 ensembles */
      /* and exactly nfiles bins, in which the observation could fall     */
      /* **************************************************************** */
      if (datafunc == TIME)
        {
          double val = 0.0;
          gridID2 = gridCreate(GRID_LONLAT, 1);
          gridDefXsize(gridID2, 1);
          gridDefYsize(gridID2, 1);
          gridDefXvals(gridID2, &val);
          gridDefYvals(gridID2, &val);
        }
      else  // datafunc == SPACE
        gridID2 = vlistInqVarGrid(vlistID1, varID);

      varIDs2[varID] = vlistDefVar(vlistID2, gridID2, zaxisID2, time_mode);
    }

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  for (varID = 0; varID < nvars; ++varID)
    {
      if (zaxisInqSize(vlistInqVarZaxis(vlistID1, varID)) > 1)
        {
          cdo_warning("More than one level not supported when processing ranked histograms.");
          cdo_warning("Try to use `cdo splitlevel` to split the dataset into levels and apply");
          cdo_warning("the operator seperately to each level.");
          cdo_abort("Exit due to unsupported file structure");
        }
    }

  if (operfunc != func_roc)
    {
      streamID2 = cdo_open_write(nfiles);
      cdo_def_vlist(streamID2, vlistID2);
    }

  auto gridsize = vlistGridsizeMax(vlistID1);

  for (int fileID = 0; fileID < nfiles; ++fileID) ef[fileID].array.resize(gridsize);

  std::vector<std::vector<int>> array2(nfiles + 1);
  if (operfunc == func_rank && datafunc == SPACE)
    {  // need to memorize data for entire grid before writing
      for (binID = 0; binID < nfiles; binID++) array2[binID].resize(gridsize, 0);
    }
  else if (operfunc == func_rank)
    {  // can process data separately for each timestep and only need to cumulate values over the grid
      for (binID = 0; binID < nfiles; binID++) array2[binID].resize(1, 0);
    }

  std::vector<int> hist;
  Varray2D<double> roc;             // receiver operating characteristics table
  Varray2D<int> ctg_tab;            // contingency table and histogram
  Varray<double> uThresh, lThresh;  // thresholds for histograms
  if (operfunc == func_roc)
    {
      hist.resize(nbins);

      roc.resize(nbins + 1);
      ctg_tab.resize(nbins + 1);
      for (int i = 0; i <= nbins; ++i)
        {
          roc[i].resize(2, 0);
          ctg_tab[i].resize(4, 0);
        }

      uThresh.resize(nbins);
      lThresh.resize(nbins);
      for (int i = 0; i < nbins; ++i)
        {
          uThresh[i] = ((double) i + 1) / nbins;
          lThresh[i] = (double) i / nbins;
        }
    }

  int tsID = 0;
  do {
      nrecs0 = cdo_stream_inq_timestep(ef[0].streamID, tsID);
      for (int fileID = 1; fileID < nfiles; ++fileID)
        {
          auto streamID = ef[fileID].streamID;
          nrecs = cdo_stream_inq_timestep(streamID, tsID);
          if (nrecs != nrecs0)
            {
              if (nrecs == 0)
                cdo_abort("Inconsistent ensemble file, too few time steps in %s!", cdo_get_stream_name(fileID));
              else
                cdo_abort("Inconsistent ensemble file, number of records at time step %d of %s and %s differ!", tsID + 1,
                          cdo_get_stream_name(0), cdo_get_stream_name(fileID));
            }
        }

      if (operfunc == func_rank && (datafunc == TIME || tsID == 0))
        {
          cdo_taxis_copy_timestep(taxisID2, taxisID1);
          if (nrecs0 > 0) cdo_def_timestep(streamID2, tsID);
        }

      //      fprintf(stderr,"TIMESTEP %i varID %i rec %i\n",tsID,varID,recID);

      for (int recID = 0; recID < nrecs0; ++recID)
        {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(ef, nfiles) private(nmiss) lastprivate(varID, levelID)
#endif
          for (int fileID = 0; fileID < nfiles; ++fileID)
            {
              auto streamID = ef[fileID].streamID;
              cdo_inq_record(streamID, &varID, &levelID);
              cdo_read_record(streamID, ef[fileID].array.data(), &nmiss);
            }

          gridID = vlistInqVarGrid(vlistID1, varID);
          gridsize = gridInqSize(gridID);
          auto missval = vlistInqVarMissval(vlistID1, varID);

          nmiss = 0;
          if (datafunc == TIME && operfunc == func_rank)
            for (binID = 0; binID < nfiles; binID++) array2[binID][0] = 0;

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(binID)
#endif
          for (size_t i = 0; i < gridsize; ++i)
            {
              auto ompthID = cdo_omp_get_thread_num();

              field[ompthID].missval = missval;
              field[ompthID].nmiss = 0;
              have_miss = 0;
              for (int fileID = 0; fileID < nfiles; ++fileID)
                {
                  field[ompthID].vec_d[fileID] = ef[fileID].array[i];
                  if (DBL_IS_EQUAL(field[ompthID].vec_d[fileID], missval))
                    {
                      have_miss = 1;
                      break;
                    }
                }

              // need to ignore all data for a gridpoint if a single ensemble
              // has a missing value at that gridpoint.
              if (!have_miss)  // only process if no missing value in ensemble
                {
                  switch (operfunc)
                    {
                    case (func_rank):
                      /* ****************/
                      /* RANK HISTOGRAM */
                      /* ************** */
                      // for ( j=0; j<nfiles; j++ )
                      //   fprintf(stderr,"%5.2g ",field[ompthID].vec_d[j]);
#ifdef _OPENMP
#pragma omp critical
#endif
                      binID = (int) field_rank(field[ompthID]);
                      // fprintf(stderr,"-->%i\n",binID);

                      if (datafunc == SPACE && !have_miss)
                        array2[binID][i]++;
                      else if (!have_miss)
                        array2[binID][0]++;

                      break;

                    case (func_roc):
                      /* ********************************** */
                      /* RECEIVER OPERATING CHARACTERISTICS */
                      /* ********************************** */
                      auto dat = &field[ompthID].vec_d[1];
                      auto ival = (field[ompthID].vec_d[0] > 0.5) ? 1 : 0;

                      for (binID = 0; binID < nbins; binID++) hist[binID] = 0;

                      for (j = 0; j < nens; ++j)
                        for (binID = 0; binID < nbins; binID++)
                          if (dat[j] >= lThresh[binID] && dat[j] < uThresh[binID]) hist[binID]++;

                      chksum = 0;
                      for (binID = 0; binID < nbins; binID++) chksum += hist[binID];

                      if (chksum != nens) exit(1);

                      int cum = 0;
                      if (ival == 1)
                        {
                          // all true positives in first bin
                          ctg_tab[0][TP] += nens;

                          cum += hist[0];
                          for (binID = 1; binID < nbins; binID++)
                            {
                              ctg_tab[binID][TP] += nens - cum;
                              ctg_tab[binID][FN] += cum;
                              cum += hist[binID];
                            }
                          ctg_tab[binID][TP] += nens - cum;
                          ctg_tab[binID][FN] += cum;
                        }
                      else if (ival == 0)
                        {
                          // all false positives in first bin
                          ctg_tab[0][FP] += nens;
                          cum += hist[0];
                          for (binID = 1; binID < nbins; binID++)
                            {
                              ctg_tab[binID][FP] += nens - cum;
                              ctg_tab[binID][TN] += cum;
                              cum += hist[binID];
                            }
                          ctg_tab[binID][FP] += nens - cum;
                          ctg_tab[binID][TN] += cum;
                        }
                      break;

                    }  // switch ( operfunc )
                }      // if ( ! have_miss )
            }          // for ( i=0; i<gridsize; i++ )

          if (datafunc == TIME && operfunc == func_rank)
            {
              for (binID = 0; binID < nfiles; binID++)
                {
                  double val = (double) array2[binID][0];
                  //		fprintf(stderr,"%i ",(int)val);
                  cdo_def_record(streamID2, varIDs2[varID], binID);
                  cdo_write_record(streamID2, &val, nmiss);
                }
              // fprintf(stderr,"\n");
            }
          else if (operfunc == func_roc)
            {
              if (DEBUG_ROC)
                {
                  fprintf(stderr, "#             :     TP     FP     FN     TN         TPR        FPR\n");

                  for (binID = 0; binID <= nbins; binID++)
                    {
                      int p = ctg_tab[binID][TP] + ctg_tab[binID][FN];
                      int n = ctg_tab[binID][FP] + ctg_tab[binID][TN];
                      double tpr = ctg_tab[binID][TP] / (double) p;
                      double fpr = ctg_tab[binID][FP] / (double) n;
                      chksum += ctg_tab[binID][0] + ctg_tab[binID][1] + ctg_tab[binID][2] + ctg_tab[binID][3];

                      roc[binID][TPR] = tpr;
                      roc[binID][FPR] = fpr;

                      fprintf(stderr, "%3i %10.4g: %6i %6i %6i %6i: %10.4g %10.4g\n", binID, (binID < nbins) ? lThresh[binID] : 1,
                              ctg_tab[binID][0], ctg_tab[binID][1], ctg_tab[binID][2], ctg_tab[binID][3], tpr, fpr);
                    }
                  fprintf(stderr, "nbins %10i\n", nbins);
                  fprintf(stderr, "#ROC CurveArea: %10.6f\n", roc_curve_integrate(roc, nbins));
                }  // if ( DEBUG_ROC )
            }      // else if (operfunc == func_roc )
        }          // for ( recID=0; recID<nrecs; recID++ )
      tsID++;
    }  // do [...]
  while (nrecs0 > 0);

  if (operfunc == func_rank)
    {
      int osize = (datafunc == TIME) ? 1 : gridsize;
      Varray<double> tmpdoub(osize);

      for (binID = 0; binID < nfiles; binID++)
        {
          for (int i = 0; i < osize; ++i) tmpdoub[i] = (double) array2[binID][i];

          cdo_def_record(streamID2, varIDs2[varID], binID);
          cdo_write_record(streamID2, tmpdoub.data(), nmiss);
        }
    }
  else if (operfunc == func_roc)
    {
      fprintf(stdout, "#             :     TP     FP     FN     TN         TPR        FPR\n");

      for (int i = 0; i <= nbins; ++i)
        {
          int p = ctg_tab[i][TP] + ctg_tab[i][FN];
          int n = ctg_tab[i][FP] + ctg_tab[i][TN];
          double tpr = ctg_tab[i][TP] / (double) p;
          double fpr = ctg_tab[i][FP] / (double) n;
          //chksum += ctg_tab[i][0] + ctg_tab[i][1] + ctg_tab[i][2] + ctg_tab[i][3];

          roc[i][TPR] = tpr;
          roc[i][FPR] = fpr;

          int sum = ctg_tab[i][TP] + ctg_tab[i][TN] + ctg_tab[i][FP] + ctg_tab[i][FN];

          fprintf(stdout, "%3i %10.4g: %6i %6i %6i %6i (%6i): %10.4g %10.4g\n", i, (i < nbins) ? lThresh[i] : 1, ctg_tab[i][0],
                  ctg_tab[i][1], ctg_tab[i][2], ctg_tab[i][3], sum, tpr, fpr);
        }

      fprintf(stdout, "#ROC CurveArea: %10.6f\n", roc_curve_integrate(roc, nbins));
    }

  for (int fileID = 0; fileID < nfiles; ++fileID) cdo_stream_close(ef[fileID].streamID);

  if (operfunc != func_roc) cdo_stream_close(streamID2);

  cdo_finish();

  return nullptr;
}
