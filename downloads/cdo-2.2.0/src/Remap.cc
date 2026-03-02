/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Interpolate remapcon        First order conservative remapping
      Interpolate remapcon2       Second order conservative remapping
      Interpolate remapbil        Bilinear interpolation
      Interpolate remapbic        Bicubic interpolation
      Interpolate remapdis        Distance-weighted averaging
      Interpolate remapnn         Nearest neighbor remapping
      Interpolate remaplaf        Largest area fraction remapping
      Genweights  gencon          Generate first order conservative remap
      Genweights  gencon2         Generate second order conservative remap
      Genweights  genbil          Generate bilinear interpolation weights
      Genweights  genbic          Generate bicubic interpolation weights
      Genweights  gendis          Generate distance-weighted averaging weights
      Genweights  gennn           Generate nearest neighbor weights
      Genweights  genlaf          Generate largest area fraction weights
      Remap       remap           SCRIP grid remapping
*/

#include <algorithm>

#include <cdi.h>

#include "process_int.h"
#include "param_conversion.h"
#include "remap.h"
#include <mpim_grid.h>
#include "griddes.h"
#include "cdo_options.h"
#include "timer.h"

enum
{
  REMAPSCON,
  REMAPCON2,
  REMAPBIL,
  REMAPBIC,
  REMAPDIS,
  REMAPNN,
  REMAPLAF,
  REMAPAVG,
  GENSCON,
  GENCON2,
  GENBIL,
  GENBIC,
  GENDIS,
  GENNN,
  GENLAF,
  REMAP,
  REMAPCON,
  REMAPYCON2,
  GENCON,
  GENYCON2
};

enum
{
  HEAP_SORT,
  MERGE_SORT
};

static RemapSwitches
operfunc_to_maptype(int operfunc)
{
  RemapSwitches remapSwitches;
  remapSwitches.remapOrder = 1;

  switch (operfunc)
    {
    case REMAPCON:
    case GENCON: remapSwitches.mapType = RemapMethod::CONSERV; break;
    case REMAPYCON2:
    case GENYCON2:
      remapSwitches.mapType = RemapMethod::CONSERV;
      remapSwitches.remapOrder = 2;
      break;
    case REMAPSCON:
    case GENSCON: remapSwitches.mapType = RemapMethod::CONSERV_SCRIP; break;
    case REMAPCON2:
    case GENCON2:
      remapSwitches.mapType = RemapMethod::CONSERV_SCRIP;
      remapSwitches.remapOrder = 2;
      break;
    case REMAPLAF:
    case GENLAF:
      remapSwitches.mapType = RemapMethod::CONSERV;
      remapSwitches.submapType = SubmapType::LAF;
      break;
    case REMAPAVG:
      remapSwitches.mapType = RemapMethod::CONSERV;
      remapSwitches.submapType = SubmapType::AVG;
      break;
    case REMAPBIL:
    case GENBIL: remapSwitches.mapType = RemapMethod::BILINEAR; break;
    case REMAPBIC:
    case GENBIC: remapSwitches.mapType = RemapMethod::BICUBIC; break;
    case REMAPDIS:
    case GENDIS:
      remapSwitches.mapType = RemapMethod::DISTWGT;
      if (remapSwitches.numNeighbors == 0) remapSwitches.numNeighbors = 4;
      break;
    case REMAPNN:
    case GENNN:
      remapSwitches.mapType = RemapMethod::DISTWGT;
      remapSwitches.numNeighbors = 1;
      break;
    default: cdo_abort("Unknown mapping method"); break;
    }

  return remapSwitches;
}

static int
maptype_to_operfunc(const RemapSwitches &remapSwitches)
{
  int operfunc = -1;

  if (remapSwitches.mapType == RemapMethod::CONSERV_SCRIP)
    operfunc = (remapSwitches.remapOrder == 2) ? REMAPCON2 : REMAPSCON;
  else if (remapSwitches.mapType == RemapMethod::CONSERV)
    operfunc = (remapSwitches.submapType == SubmapType::LAF) ? REMAPLAF : ((remapSwitches.remapOrder == 2) ? REMAPYCON2 : REMAPCON);
  else if (remapSwitches.mapType == RemapMethod::BILINEAR)
    operfunc = REMAPBIL;
  else if (remapSwitches.mapType == RemapMethod::BICUBIC)
    operfunc = REMAPBIC;
  else if (remapSwitches.mapType == RemapMethod::DISTWGT)
    operfunc = (remapSwitches.numNeighbors == 1) ? REMAPNN : REMAPDIS;
  else
    cdo_abort("Unsupported mapping method (mapType = %d)", remapSwitches.mapType);

  return operfunc;
}

static void
remap_print_info(int operfunc, bool remap_genweights, const RemapGrid &srcGrid, const RemapGrid &tgtGrid, size_t nmiss)
{
  auto srcGridName = is_healpix_grid(srcGrid.gridID) ? "healpix" : gridNamePtr(gridInqType(srcGrid.gridID));
  auto tgtGridName = is_healpix_grid(tgtGrid.gridID) ? "healpix" : gridNamePtr(gridInqType(tgtGrid.gridID));

  std::string outStr;
  // clang-format off
  if      (operfunc == REMAPBIL   || operfunc == GENBIL)   outStr = "Bilinear";
  else if (operfunc == REMAPBIC   || operfunc == GENBIC)   outStr = "Bicubic";
  else if (operfunc == REMAPNN    || operfunc == GENNN)    outStr = "Nearest neighbor";
  else if (operfunc == REMAPDIS   || operfunc == GENDIS)   outStr = "Distance-weighted average";
  else if (operfunc == REMAPSCON  || operfunc == GENSCON)  outStr = "SCRIP first order conservative";
  else if (operfunc == REMAPCON2  || operfunc == GENCON2)  outStr = "SCRIP second order conservative";
  else if (operfunc == REMAPLAF   || operfunc == GENLAF)   outStr = "YAC largest area fraction";
  else if (operfunc == REMAPCON   || operfunc == GENCON)   outStr = "YAC first order conservative";
  else if (operfunc == REMAPYCON2 || operfunc == GENYCON2) outStr = "YAC second order conservative";
  else if (operfunc == REMAPAVG)                           outStr = "Average";
  else                                                     outStr = "Unknown";
  // clang-format on

  outStr += remap_genweights ? " weights from " : " remapping from ";
  outStr += srcGridName;
  outStr += " (" + std::to_string(srcGrid.dims[0]);
  if (srcGrid.rank == 2) outStr += "x" + std::to_string(srcGrid.dims[1]);
  outStr += ") to ";
  outStr += tgtGridName;
  outStr += " (" + std::to_string(tgtGrid.dims[0]);
  if (tgtGrid.rank == 2) outStr += "x" + std::to_string(tgtGrid.dims[1]);
  outStr += ") grid";

  if (nmiss) outStr += ", with source mask (" + std::to_string(gridInqSize(srcGrid.gridID) - nmiss) + ")";

  cdo_print(outStr);
}

static void
remap_print_warning(const std::string &remapWeightsFile, int operfunc, RemapGrid &srcGrid, size_t nmiss)
{
  (void) operfunc;

  std::string outStr = "Remap weights from ";
  outStr += remapWeightsFile;
  outStr += " not used, ";
  outStr += gridNamePtr(gridInqType(srcGrid.gridID));
  outStr += " (" + std::to_string(srcGrid.dims[0]);
  if (srcGrid.rank == 2) outStr += "x" + std::to_string(srcGrid.dims[1]);
  outStr += ") grid";

  if (nmiss) outStr += " with mask (" + std::to_string(gridInqSize(srcGrid.gridID) - nmiss) + ")";

  outStr += " not found!";

  cdo_warning(outStr);
}

double pointSearchRadius = 180.0;

static double remap_threshhold = 2.0;
static int remapTest = 0;
static int remap_non_global = false;
static int remap_num_srch_bins = 180;
static bool lremap_num_srch_bins = false;
static bool remap_extrapolate = false;
static bool doExtrapolate = false;
static int MaxRemaps = -1;
static int sort_mode = HEAP_SORT;
static double remapFracMin = 0.0;

std::string
getenv_string(const std::string &envVar)
{
  std::string envString;
  auto envCstring = getenv(envVar.c_str());
  if (envCstring) envString = envCstring;
  return envString;
}

static void
remapGetenv(void)
{
  char *envstr;

  envstr = getenv("MAX_REMAPS");
  if (envstr && *envstr)
    {
      auto ival = atoi(envstr);
      if (ival > 0)
        {
          MaxRemaps = ival;
          if (Options::cdoVerbose) cdo_print("Set MAX_REMAPS to %d", MaxRemaps);
        }
    }

  envstr = getenv("REMAP_MAX_ITER");
  if (envstr && *envstr)
    {
      auto ival = atoi(envstr);
      if (ival > 0)
        {
          remap_set_int(REMAP_MAX_ITER, ival);
          if (Options::cdoVerbose) cdo_print("Set REMAP_MAX_ITER to %d", ival);
        }
    }

  envstr = getenv("REMAP_TEST");
  if (envstr && *envstr)
    {
      auto ival = atoi(envstr);
      if (ival > 0)
        {
          remapTest = ival;
          if (Options::cdoVerbose) cdo_print("Set REMAP_TEST to %d", remapTest);
        }
    }

#ifdef _OPENMP
  sort_mode = (Threading::ompNumThreads == 1) ? HEAP_SORT : MERGE_SORT;
#endif

  {
    auto envString = getenv_string("REMAP_SORT_MODE");
    if (envString.size())
      {
        if (envString == "heap")
          sort_mode = HEAP_SORT;
        else if (envString == "merge")
          sort_mode = MERGE_SORT;

        if (Options::cdoVerbose) cdo_print("Set sort_mode to %s", (sort_mode == HEAP_SORT) ? "HEAP_SORT" : "MERGE_SORT");
      }
  }

  envstr = getenv("REMAP_THRESHHOLD");
  if (envstr && *envstr)
    {
      auto fval = atof(envstr);
      if (fval > 0.0)
        {
          remap_threshhold = fval;
          if (Options::cdoVerbose) cdo_print("Set REMAP_THRESHHOLD to %g", remap_threshhold);
        }
    }

  remap_set_threshhold(remap_threshhold);

  {
    auto envString = getenv_string("CDO_REMAP_RADIUS");
    if (envString.size())
      {
        auto fval = radius_str_to_deg(envString);
        if (fval < 0.0 || fval > 180.0) cdo_abort("%s=%g out of bounds (0-180 deg)!", "CDO_REMAP_RADIUS", fval);
        pointSearchRadius = fval;
        if (Options::cdoVerbose) cdo_print("Set CDO_REMAP_RADIUS to %g", pointSearchRadius);
      }
  }

  {
    auto envString = getenv_string("CDO_GRIDSEARCH_RADIUS");
    if (envString.size())
      {
        auto fval = radius_str_to_deg(envString);
        if (fval < 0.0 || fval > 180.0) cdo_abort("%s=%g out of bounds (0-180 deg)!", "CDO_GRIDSEARCH_RADIUS", fval);
        pointSearchRadius = fval;
        if (Options::cdoVerbose) cdo_print("Set CDO_GRIDSEARCH_RADIUS to %g", pointSearchRadius);
      }
  }

  if (Options::cdoVerbose) cdo_print("Point search radius = %g deg", pointSearchRadius);

  envstr = getenv("REMAP_AREA_MIN");
  if (envstr && *envstr)
    {
      auto fval = atof(envstr);
      if (fval > 0.0)
        {
          remapFracMin = fval;
          if (Options::cdoVerbose) cdo_print("Set REMAP_AREA_MIN to %g", remapFracMin);
        }
    }

  envstr = getenv("REMAP_NUM_SRCH_BINS");
  if (envstr && *envstr)
    {
      auto ival = atoi(envstr);
      if (ival > 0)
        {
          remap_num_srch_bins = ival;
          lremap_num_srch_bins = true;
          if (Options::cdoVerbose) cdo_print("Set REMAP_NUM_SRCH_BINS to %d", remap_num_srch_bins);
        }
    }

  envstr = getenv("REMAP_NON_GLOBAL");
  if (envstr && *envstr)
    {
      auto ival = atoi(envstr);
      if (ival >= 0)
        {
          remap_non_global = ival;
          if (Options::cdoVerbose) cdo_print("Set REMAP_NON_GLOBAL to %d", remap_non_global);
        }
    }

  {
    auto envString = getenv_string("REMAP_EXTRAPOLATE");
    if (envString.size())
      {
        if (envString == "ON" || envString == "on")
          {
            doExtrapolate = true;
            remap_extrapolate = true;
          }
        else if (envString == "OFF" || envString == "off")
          {
            doExtrapolate = true;
            remap_extrapolate = false;
          }
        else
          cdo_warning("Environment variable REMAP_EXTRAPOLATE has wrong value!");

        if (Options::cdoVerbose) cdo_print("Extrapolation %s!", remap_extrapolate ? "enabled" : "disabled");
      }
  }

  {
    auto envString = getenv_string("CDO_REMAP_GENWEIGHTS");
    if (envString.size())
      {
        if (envString == "ON" || envString == "on")
          Options::REMAP_genweights = true;
        else if (envString == "OFF" || envString == "off")
          Options::REMAP_genweights = false;
        else
          cdo_warning("Environment variable CDO_REMAP_GENWEIGHTS has wrong value!");

        if (Options::cdoVerbose) cdo_print("Generation of weights %s!", Options::REMAP_genweights ? "enabled" : "disabled");
      }
  }
}

static bool
gridIsGlobal(int gridID)
{
  auto isGlobalGrid = true;
  auto isNonGlobal = (remap_non_global || !gridIsCircular(gridID));

  auto gridtype = gridInqType(gridID);
  if (gridProjIsSupported(gridID) || (gridtype == GRID_LONLAT && isNonGlobal) || (gridtype == GRID_CURVILINEAR && isNonGlobal))
    isGlobalGrid = false;

  return isGlobalGrid;
}

template <typename T>
static void
scale_gridbox_area(size_t gridsize1, const Varray<T> &array1, size_t gridsize2, Varray<T> &array2, const Varray<double> &grid2_area)
{
  auto array1sum = varray_sum(gridsize1, array1);
  auto array2sum = varray_sum(gridsize2, grid2_area);
  for (size_t i = 0; i < gridsize2; ++i) array2[i] = grid2_area[i] / array2sum * array1sum;

  static auto hasGridboxInfo = true;
  if (hasGridboxInfo)
    {
      cdo_print("gridbox_area replaced and scaled to %g", array1sum);
      hasGridboxInfo = false;
    }
}

static void
scale_gridbox_area(const Field &field1, Field &field2, const Varray<double> &grid2_area)
{
  scale_gridbox_area(field1.gridsize, field1.vec_f, field2.gridsize, field2.vec_f, grid2_area);
  scale_gridbox_area(field1.gridsize, field1.vec_d, field2.gridsize, field2.vec_d, grid2_area);
}

template <typename T>
void gme_grid_restore(T *p, int ni, int nd);

template <typename T>
static void
restore_gme_grid(int gridID, size_t srvGridSize, Varray<T> &v, size_t tgtGridSize, const Varray<int> &vgpm)
{
  int nd, ni, ni2, ni3;
  gridInqParamGME(gridID, &nd, &ni, &ni2, &ni3);

  auto n = tgtGridSize;
  for (size_t i = srvGridSize; i > 0; i--)
    if (vgpm[i - 1]) v[i - 1] = v[--n];

  gme_grid_restore(v.data(), ni, nd);
}

static void
restore_gme_grid(Field &field, const RemapGrid &tgtGrid)
{
  if (field.memType == MemType::Float)
    restore_gme_grid(field.grid, field.gridsize, field.vec_f, tgtGrid.size, tgtGrid.vgpm);
  else
    restore_gme_grid(field.grid, field.gridsize, field.vec_d, tgtGrid.size, tgtGrid.vgpm);
}

template <typename T>
static void
store_gme_grid(size_t gridsize, Varray<T> &v, const Varray<int> &vgpm)
{
  for (size_t i = 0, j = 0; i < gridsize; ++i)
    if (vgpm[i]) v[j++] = v[i];
}

static void
store_gme_grid(Field &field, const Varray<int> &vgpm)
{
  if (field.memType == MemType::Float)
    store_gme_grid(field.gridsize, field.vec_f, vgpm);
  else
    store_gme_grid(field.gridsize, field.vec_d, vgpm);
}

static int
set_remapgrids(int filetype, int vlistID, const VarList &varList, int ngrids, std::vector<bool> &remapgrids)
{
  (void) filetype;
  for (int index = 0; index < ngrids; ++index)
    {
      remapgrids[index] = true;

      auto gridID = vlistGrid(vlistID, index);
      auto gridtype = gridInqType(gridID);
      auto hasProjParams = ((gridtype == GRID_PROJECTION) && grid_has_proj_params(gridID));

      if (!gridProjIsSupported(gridID) && !hasProjParams && gridtype != GRID_LONLAT && gridtype != GRID_GAUSSIAN
          && gridtype != GRID_GME && gridtype != GRID_CURVILINEAR && gridtype != GRID_UNSTRUCTURED
          && gridtype != GRID_GAUSSIAN_REDUCED)
        {
          if (gridtype == GRID_GENERIC && gridInqSize(gridID) <= 2) { remapgrids[index] = false; }
          else
            {
              auto nvars = vlistNvars(vlistID);
              for (int varID = 0; varID < nvars; ++varID)
                if (gridID == varList[varID].gridID)
                  {
                    cdo_abort("Unsupported %s coordinates (Variable: %s)!", gridNamePtr(gridtype), varList[varID].name);
                    break;
                  }
            }
        }
    }

  int index;
  for (index = 0; index < ngrids; ++index)
    if (remapgrids[index]) break;
  if (index == ngrids) cdo_abort("No remappable grid found!");

  return index;
}

static int
set_max_remaps(int vlistID)
{
  int maxRemaps = 0;

  auto nzaxis = vlistNzaxis(vlistID);
  for (int index = 0; index < nzaxis; ++index)
    {
      auto zaxisID = vlistZaxis(vlistID, index);
      auto zaxisSize = zaxisInqSize(zaxisID);
      if (zaxisSize > maxRemaps) maxRemaps = zaxisSize;
    }

  auto nvars = vlistNvars(vlistID);
  if (nvars > maxRemaps) maxRemaps = nvars;

  maxRemaps++;

  if (Options::cdoVerbose) cdo_print("Set maxRemaps to %d", maxRemaps);

  return maxRemaps;
}

static NormOpt
get_normOpt(void)
{
  // clang-format off
  NormOpt normOpt(NormOpt::FRACAREA);

  auto envString = getenv_string("CDO_REMAP_NORM");
  if (envString.size())
    {
      if      (envString == "frac") normOpt = NormOpt::FRACAREA;
      else if (envString == "dest") normOpt = NormOpt::DESTAREA;
      else if (envString == "none") normOpt = NormOpt::NONE;
      else cdo_warning("CDO_REMAP_NORM=%s unsupported!", envString);
    }

  if (Options::cdoVerbose)
    {
      const char *outStr = "none";
      if      (normOpt == NormOpt::FRACAREA) outStr = "frac";
      else if (normOpt == NormOpt::DESTAREA) outStr = "dest";
      cdo_print("Normalization option: %s", outStr);
    }
  // clang-format on

  return normOpt;
}

template <typename T>
static void
remap_normalize_field(NormOpt normOpt, size_t gridsize, Varray<T> &array, T missval, const RemapGrid &tgtGrid)
{
  // used to check the result of remapcon

  if (normOpt == NormOpt::NONE)
    {
      for (size_t i = 0; i < gridsize; ++i)
        {
          if (!dbl_is_equal(array[i], missval))
            {
              auto gridError = tgtGrid.cell_frac[i] * tgtGrid.cell_area[i];
              array[i] = (std::fabs(gridError) > 0.0) ? array[i] / gridError : missval;
            }
        }
    }
  else if (normOpt == NormOpt::DESTAREA)
    {
      for (size_t i = 0; i < gridsize; ++i)
        {
          if (!dbl_is_equal(array[i], missval))
            {
              array[i] = (std::fabs(tgtGrid.cell_frac[i]) > 0.0) ? array[i] / tgtGrid.cell_frac[i] : missval;
            }
        }
    }
}

static void
remap_normalize_field(NormOpt normOpt, Field &field, const RemapGrid &tgtGrid)
{
  if (field.memType == MemType::Float)
    remap_normalize_field(normOpt, field.gridsize, field.vec_f, (float) field.missval, tgtGrid);
  else
    remap_normalize_field(normOpt, field.gridsize, field.vec_d, field.missval, tgtGrid);
}

template <typename T>
static void
remap_set_frac_min(size_t gridsize, Varray<T> &array, T missval, RemapGrid *tgtGrid)
{
  if (remapFracMin > 0.0)
    {
      for (size_t i = 0; i < gridsize; ++i)
        if (tgtGrid->cell_frac[i] < remapFracMin) array[i] = missval;
    }
}

static void
remap_set_frac_min(Field &field, RemapGrid *tgtGrid)
{
  if (field.memType == MemType::Float)
    remap_set_frac_min(field.gridsize, field.vec_f, (float) field.missval, tgtGrid);
  else
    remap_set_frac_min(field.gridsize, field.vec_d, field.missval, tgtGrid);
}

static void
remap_set_mask(size_t gridsize, const Field &field1, double missval, Varray<short> &imask)
{
  if (field1.memType == MemType::Float)
    remap_set_mask(gridsize, field1.vec_f, (float) missval, imask);
  else
    remap_set_mask(gridsize, field1.vec_d, missval, imask);
}

int timer_remap, timer_remap_init, timer_remap_sort;

static void
remapTimerInit(void)
{
  timer_remap = timer_new("remap");
  timer_remap_init = timer_new("remap init");
  timer_remap_sort = timer_new("remap sort");
}

static void
remapLinksPerValue(RemapVars &rv)
{
  long lpv = -1;

  rv.linksPerValue = lpv;
  return;
  /*
  size_t numLinks = rv.numLinks;

  if ( numLinks > 0 )
    {
      lpv = 1;
      auto &tgt_add = rv.tgtCellIndices;
      size_t n = 0;
      size_t ival = tgt_add[n];
      for ( n = 1; n < numLinks; ++n )
        if ( tgt_add[n] == ival ) lpv++;
        else break;

    printf("lpv %zu\n", lpv);

      if ( numLinks%lpv != 0 ) lpv = 0;

      n++;
      if ( n < numLinks )
        {
          size_t lpv2 = 0;
          size_t ival2 = tgt_add[n];
          for ( ; n < numLinks; ++n )
            if ( tgt_add[n] == ival2 ) lpv2++;
            else if ( lpv == lpv2 )
              {
                lpv2 = 0;
                ival2 = tgt_add[n];
              }
            else
              {
                lpv = 0;
                break;
              }
        }

  printf("lpv %zu\n", lpv);

      if ( lpv == 1 )
        {
          for ( size_t n = 1; n < numLinks; ++n )
            {
              if ( tgt_add[n] == tgt_add[n-1] )
                {
                  lpv = 0;
                  break;
                }
            }
        }
      else if ( lpv > 1 )
        {
          for ( size_t n = 1; n < numLinks/lpv; ++n )
            {
              ival = tgt_add[n*lpv];
              for ( size_t k = 1; k < lpv; ++k )
                {
                  if ( tgt_add[n*lpv+k] != ival )
                    {
                      lpv = 0;
                      break;
                    }
                }
              if ( lpv == 0 ) break;
            }
        }
    }

  printf("lpv %zu\n", lpv);

  rv.linksPerValue = lpv;
  */
}

static void
remapSortAddr(RemapVars &rv)
{
  if (Options::Timer) timer_start(timer_remap_sort);
  if (sort_mode == MERGE_SORT)
    { /*
      ** use a combination of the old sort_add and a split and merge approach.
      ** The chunk size is determined by MERGE_SORT_LIMIT_SIZE in remaplib.c.
      ** OpenMP parallelism is supported
      */
      sort_iter(rv.numLinks, rv.num_wts, &rv.tgtCellIndices[0], &rv.srcCellIndices[0], &rv.wts[0], Threading::ompNumThreads);
    }
  else
    { /* use a pure heap sort without any support of parallelism */
      sort_add(rv.numLinks, rv.num_wts, &rv.tgtCellIndices[0], &rv.srcCellIndices[0], &rv.wts[0]);
    }
  if (Options::Timer) timer_stop(timer_remap_sort);
}

static int
remap_gen_num_bins(int ysize)
{
  constexpr int maxbins = 720;
  int num_srch_bins = ysize / 2 + ysize % 2;
  if (num_srch_bins > maxbins) num_srch_bins = maxbins;
  if (num_srch_bins < 1) num_srch_bins = 1;
  return num_srch_bins;
}

void
remap_init(RemapType &remap)
{
  remap.nused = 0;
  remap.gridID = -1;
  remap.gridsize = 0;
  remap.nmiss = 0;
}

static void
remap_gen_weights(const RemapSwitches &remapSwitches, RemapType &remap)
{
  auto mapType = remapSwitches.mapType;
  auto numNeighbors = remapSwitches.numNeighbors;
  // clang-format off
  if      (mapType == RemapMethod::CONSERV_SCRIP) remap_conserv_weights_scrip(remap.search, remap.vars);
  else if (mapType == RemapMethod::BILINEAR)      remap_bilinear_weights(remap.search, remap.vars);
  else if (mapType == RemapMethod::BICUBIC)       remap_bicubic_weights(remap.search, remap.vars);
  else if (mapType == RemapMethod::DISTWGT)       remap_distwgt_weights(numNeighbors, remap.search, remap.vars);
  else if (mapType == RemapMethod::CONSERV)       remap_conserv_weights(remap.search, remap.vars);
  // clang-format on

  if (remap.vars.sort_add) remapSortAddr(remap.vars);
  if (remap.vars.linksPerValue == -1) remapLinksPerValue(remap.vars);
}

static void
remap_field(const RemapSwitches &remapSwitches, RemapType &remap, const Field &field1, Field &field2)
{
  auto mapType = remapSwitches.mapType;
  auto numNeighbors = remapSwitches.numNeighbors;
  // clang-format off
  if      (mapType == RemapMethod::BILINEAR) remap_bilinear(remap.search, field1, field2);
  else if (mapType == RemapMethod::BICUBIC)  remap_bicubic(remap.search, field1, field2);
  else if (mapType == RemapMethod::DISTWGT)  remap_dist_wgt(numNeighbors, remap.search, field1, field2);
  else if (mapType == RemapMethod::CONSERV)  remap_conserv(remap.vars.normOpt, remap.search, field1, field2);
  // clang-format on
}

static RemapSwitches
remap_read_weights(const std::string &remapWeightsFile, int gridID1, int gridID2, RemapType &remap0)
{
  auto remapSwitches = remap_read_data_scrip(remapWeightsFile, gridID1, gridID2, remap0.srcGrid, remap0.tgtGrid, remap0.vars);

  if (remap0.vars.linksPerValue == 0) remapLinksPerValue(remap0.vars);

  auto gridsize = remap0.srcGrid.size;
  remap0.gridID = gridID1;
  remap0.gridsize = gridInqSize(gridID1);

  if (remapSwitches.mapType == RemapMethod::DISTWGT && !doExtrapolate) remap_extrapolate = true;
  if (gridIsCircular(gridID1) && !doExtrapolate) remap_extrapolate = true;

  if (gridInqType(gridID1) == GRID_GME) gridsize = remap0.srcGrid.nvgp;

  if (gridsize != remap0.gridsize) cdo_abort("Size of source grid and weights from %s differ!", remapWeightsFile);

  if (gridInqType(gridID1) == GRID_GME) gridsize = remap0.srcGrid.size;

  for (size_t i = 0; i < gridsize; ++i)
    if (remap0.srcGrid.mask[i] == false) remap0.nmiss++;

  auto gridsize2 = gridInqSize(gridID2);
  if (gridInqType(gridID2) == GRID_GME)
    {
      remap0.tgtGrid.nvgp = gridInqSize(gridID2);
      remap0.tgtGrid.vgpm.resize(gridInqSize(gridID2));
      auto gridID2_gme = gridToUnstructured(gridID2, NeedCorners::Yes);
      gridInqMaskGME(gridID2_gme, &remap0.tgtGrid.vgpm[0]);
      gridDestroy(gridID2_gme);
      size_t isize = 0;
      for (size_t i = 0; i < gridsize2; ++i)
        if (remap0.tgtGrid.vgpm[i]) isize++;
      gridsize2 = isize;
    }
  // printf("grid2 %zu %d %zu\n", gridsize2, remap0.tgtGrid.nvgp, remap0.tgtGrid.size);
  if (remap0.tgtGrid.size != gridsize2) cdo_abort("Size of target grid and weights from %s differ!", remapWeightsFile);

  return remapSwitches;
}

static void
add_operators(void)
{
  // clang-format off
  cdo_operator_add("remap",        REMAP,        0, nullptr);
  cdo_operator_add("remapcon",     REMAPCON,     0, nullptr);
  cdo_operator_add("remapycon2test",REMAPYCON2,   0, nullptr);
  cdo_operator_add("remapscon",    REMAPSCON,    0, nullptr);
  cdo_operator_add("remapcon2",    REMAPCON2,    0, nullptr);
  cdo_operator_add("remapbil",     REMAPBIL,     0, nullptr);
  cdo_operator_add("remapbic",     REMAPBIC,     0, nullptr);
  cdo_operator_add("remapdis",     REMAPDIS,     0, nullptr);
  cdo_operator_add("remapnn",      REMAPNN,      0, nullptr);
  cdo_operator_add("remaplaf",     REMAPLAF,     0, nullptr);
  cdo_operator_add("remapavgtest", REMAPAVG,     0, nullptr);
  cdo_operator_add("gencon",       GENCON,       1, nullptr);
  cdo_operator_add("genycon2test", GENYCON2,     1, nullptr);
  cdo_operator_add("genscon",      GENSCON,      1, nullptr);
  cdo_operator_add("gencon2",      GENCON2,      1, nullptr);
  cdo_operator_add("genbil",       GENBIL,       1, nullptr);
  cdo_operator_add("genbic",       GENBIC,       1, nullptr);
  cdo_operator_add("gendis",       GENDIS,       1, nullptr);
  cdo_operator_add("gennn",        GENNN,        1, nullptr);
  cdo_operator_add("genlaf",       GENLAF,       1, nullptr);
  // clang-format on
}

void *
Remap(void *argument)
{
  RemapSwitches remapSwitches;
  bool remap_genweights = Options::REMAP_genweights;
  int remapIndex = -1;
  int numRemaps = 0;
  int numNeighbors = 0;
  std::string remapWeightsFile;

  if (Options::Timer) remapTimerInit();

  cdo_initialize(argument);

  add_operators();

  auto operatorID = cdo_operator_id();
  auto operfunc = cdo_operator_f1(operatorID);
  bool writeRemapWeightsOnly = cdo_operator_f2(operatorID);
  auto doRemap = (operfunc == REMAP);

  remap_set_int(REMAP_WRITE_REMAP, writeRemapWeightsOnly);

  if (operfunc == REMAPDIS || operfunc == GENDIS || operfunc == REMAPNN || operfunc == GENNN) remap_extrapolate = true;

  remapGetenv();

  if (Options::cdoVerbose) cdo_print("Extrapolation %s!", remap_extrapolate ? "enabled" : "disabled");

  if (doRemap)
    {
      operator_input_arg("grid description file or name, remap weights file (SCRIP NetCDF)");
      operator_check_argc(2);
      remapWeightsFile = cdo_operator_argv(1);
    }
  else
    {
      operator_input_arg("grid description file or name");
      if (operfunc == REMAPDIS && cdo_operator_argc() == 2)
        {
          auto inum = parameter_to_int(cdo_operator_argv(1));
          if (inum < 1) cdo_abort("Number of nearest neighbors out of range (>0)!");
          numNeighbors = inum;
        }
      else
        {
          operator_check_argc(1);
        }
    }

  auto gridID2 = cdo_define_grid(cdo_operator_argv(0));
  if (gridInqType(gridID2) == GRID_GENERIC) cdo_abort("Unsupported target grid type (generic)!");

  auto streamID1 = cdo_open_read(0);

  auto filetype = cdo_inq_filetype(streamID1);

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto vlistID2 = vlistDuplicate(vlistID1);

  auto taxisID1 = vlistInqTaxis(vlistID1);
  auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  VarList varList1;
  varListInit(varList1, vlistID1);

  auto ngrids = vlistNgrids(vlistID1);
  std::vector<bool> remapgrids(ngrids);
  auto gridID1 = vlistGrid(vlistID1, set_remapgrids(filetype, vlistID1, varList1, ngrids, remapgrids));

  for (int index = 0; index < ngrids; ++index)
    if (remapgrids[index]) vlistChangeGridIndex(vlistID2, index, gridID2);

  if (MaxRemaps == -1) MaxRemaps = set_max_remaps(vlistID1);
  if (MaxRemaps < 1) cdo_abort("max_remaps out of range (>0)!");

  std::vector<RemapType> remaps(MaxRemaps);

  if (writeRemapWeightsOnly || doRemap) remap_genweights = true;

  if (doRemap)
    {
      numRemaps = 1;
      remapSwitches = remap_read_weights(remapWeightsFile, gridID1, gridID2, remaps[0]);
      operfunc = maptype_to_operfunc(remapSwitches);

      if (remapTest) remap_vars_reorder(remaps[0].vars);
    }
  else
    {
      remapSwitches = operfunc_to_maptype(operfunc);
      if (numNeighbors) remapSwitches.numNeighbors = numNeighbors;
    }

  auto mapType = remapSwitches.mapType;
  auto remapOrder = remapSwitches.remapOrder;

  VarList varList2;
  varListInit(varList2, vlistID2);

  // if (!remap_genweights && (mapType == RemapMethod::CONSERV_SCRIP || mapType == RemapMethod::CONSERV)) remap_genweights = true;
  if (!remap_genweights && mapType == RemapMethod::CONSERV_SCRIP) remap_genweights = true;

  bool useMask = !(!remap_genweights
                   && (mapType == RemapMethod::BILINEAR || mapType == RemapMethod::BICUBIC || mapType == RemapMethod::DISTWGT
                       || mapType == RemapMethod::CONSERV));

  remap_set_int(REMAP_GENWEIGHTS, (int) remap_genweights);

  NormOpt normOpt(NormOpt::NONE);
  if (mapType == RemapMethod::CONSERV_SCRIP || mapType == RemapMethod::CONSERV) normOpt = get_normOpt();

  auto grid1sizemax = vlistGridsizeMax(vlistID1);

  auto needGradients = (mapType == RemapMethod::BICUBIC);
  if ((mapType == RemapMethod::CONSERV_SCRIP || mapType == RemapMethod::CONSERV) && remapOrder == 2)
    {
      if (Options::cdoVerbose) cdo_print("Second order remapping");
      needGradients = true;
    }

  RemapGradients gradients;
  if (needGradients) gradients.init(grid1sizemax);

  if (remap_genweights)
    {
      // remap() gives rounding errors on target arrays with single precision floats
      auto nvars = vlistNvars(vlistID2);
      for (int varID = 0; varID < nvars; ++varID) varList2[varID].memType = MemType::Double;
    }

  Field field1, field2;

  Varray<short> imask;

  CdoStreamID streamID2 = CDO_STREAM_UNDEF;
  if (!writeRemapWeightsOnly)
    {
      streamID2 = cdo_open_write(1);
      cdo_def_vlist(streamID2, vlistID2);
    }

  int tsID = 0;
  while (true)
    {
      auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
      if (nrecs == 0) break;

      cdo_taxis_copy_timestep(taxisID2, taxisID1);

      if (!writeRemapWeightsOnly) cdo_def_timestep(streamID2, tsID);

      for (int recID = 0; recID < nrecs; ++recID)
        {
          int varID, levelID;
          cdo_inq_record(streamID1, &varID, &levelID);
          field1.init(varList1[varID]);
          cdo_read_record(streamID1, field1);
          auto nmiss1 = useMask ? field1.nmiss : 0;

          field2.init(varList2[varID]);

          auto gridID = varList1[varID].gridID;
          auto missval = varList1[varID].missval;
          auto gridsize = varList1[varID].gridsize;

          auto skipVar = false;
          if (!remapgrids[vlistGridIndex(vlistID1, gridID)])
            {
              if (writeRemapWeightsOnly) continue;

              field_copy(field1, field2);
              skipVar = true;
            }

          if (!skipVar)
            {
              if (mapType != RemapMethod::CONSERV_SCRIP && mapType != RemapMethod::CONSERV && gridInqType(gridID) == GRID_GME)
                cdo_abort("Only conservative remapping is available to remap between GME grids!");

              if (gridIsCircular(gridID) && !doExtrapolate) remap_extrapolate = true;

              imask.resize(gridsize, 1);

              if (nmiss1) remap_set_mask(gridsize, field1, missval, imask);

              for (remapIndex = numRemaps - 1; remapIndex >= 0; remapIndex--)
                {
                  auto &remap = remaps[remapIndex];
                  if (gridID == remap.gridID)
                    {
                      if ((useMask && nmiss1 == remap.nmiss && imask == remap.srcGrid.mask) || !useMask)
                        {
                          remap.nused++;
                          break;
                        }
                    }
                }

              if (Options::cdoVerbose && remapIndex >= 0) cdo_print("Using remap %d", remapIndex);

              if (remapIndex < 0)
                {
                  if (numRemaps < MaxRemaps)
                    {
                      remapIndex = numRemaps;
                      numRemaps++;
                    }
                  else
                    {
                      int remapIndex0 = (MaxRemaps > 1 && remaps[0].nused > remaps[1].nused);
                      remap_vars_free(remaps[remapIndex0].vars);
                      remap_grid_free(remaps[remapIndex0].srcGrid);
                      remap_grid_free(remaps[remapIndex0].tgtGrid);
                      remap_search_free(remaps[remapIndex0].search);
                      remapIndex = remapIndex0;
                      remap_init(remaps[remapIndex]);
                    }

                  auto &remap = remaps[remapIndex];
                  if (remap.gridID != gridID)
                    {
                      if (gridIsCircular(gridID) && !doExtrapolate) remap_extrapolate = true;

                      //  remap.srcGrid.luse_cell_area = false;
                      //  remap.tgtGrid.luse_cell_area = false;

                      if (gridInqType(gridID) != GRID_UNSTRUCTURED && !lremap_num_srch_bins)
                        {
                          remap_num_srch_bins = (!remap_extrapolate && mapType == RemapMethod::DISTWGT)
                                                    ? 1
                                                    : remap_gen_num_bins(gridInqYsize(gridID));
                        }

                      remap_set_int(REMAP_NUM_SRCH_BINS, remap_num_srch_bins);

                      remap.vars.normOpt = normOpt;
                      remap.vars.pinit = false;

                      if ((mapType == RemapMethod::BILINEAR || mapType == RemapMethod::BICUBIC)
                          && (gridInqType(gridID) == GRID_GME || gridInqType(gridID) == GRID_UNSTRUCTURED))
                        cdo_abort("Bilinear/bicubic interpolation doesn't support unstructured source grids!");

                      // Initialize grid information for both grids
                      if (Options::Timer) timer_start(timer_remap_init);
                      remap_init_grids(mapType, remap_extrapolate, gridID, remap.srcGrid, gridID2, remap.tgtGrid);
                      remap_search_init(mapType, remap.search, remap.srcGrid, remap.tgtGrid);
                      if (Options::Timer) timer_stop(timer_remap_init);
                    }

                  remap.gridID = gridID;
                  remap.nmiss = nmiss1;

                  if (gridInqType(gridID) == GRID_GME)
                    {
                      for (size_t i = 0, j = 0; i < gridsize; ++i)
                        if (remap.srcGrid.vgpm[i]) imask[j++] = imask[i];
                    }

                  varray_copy(remap.srcGrid.size, imask, remap.srcGrid.mask);

                  if (mapType == RemapMethod::CONSERV_SCRIP || mapType == RemapMethod::CONSERV)
                    {
                      varray_fill(remap.srcGrid.cell_area, 0.0);
                      varray_fill(remap.srcGrid.cell_frac, 0.0);
                      varray_fill(remap.tgtGrid.cell_area, 0.0);
                    }
                  varray_fill(remap.tgtGrid.cell_frac, 0.0);

                  // initialize some remapping variables
                  remap_vars_init(mapType, remapOrder, remap.vars);

                  remap_print_info(operfunc, remap_genweights, remap.srcGrid, remap.tgtGrid, nmiss1);

                  if (remap_genweights)
                    {
                      if (needGradients)
                        {
                          if (remap.srcGrid.rank != 2 && remapOrder == 2)
                            cdo_abort("Second order remapping is not available for unstructured grids!");
                        }

                      remap_gen_weights(remapSwitches, remap);

                      if (writeRemapWeightsOnly) goto WRITE_REMAP;

                      if (remapTest) remap_vars_reorder(remap.vars);
                    }
                }

              auto &remap = remaps[remapIndex];
              if (gridInqType(gridID) == GRID_GME) store_gme_grid(field1, remaps[remapIndex].srcGrid.vgpm);

              auto gridsize2 = gridInqSize(gridID2);

              if (remap_genweights)
                {
                  remap.nused++;

                  if (needGradients)
                    {
                      if (remap.srcGrid.rank != 2 && remapOrder == 2)
                        cdo_abort("Second order remapping is not available for unstructured grids!");

                      remap_gradients(remap.srcGrid, field1, gradients);
                    }

                  if (operfunc == REMAPLAF)
                    remap_laf(field2, missval, gridsize2, remap.vars, field1);
                  else if (operfunc == REMAPAVG)
                    remap_avg(field2, missval, gridsize2, remap.vars, field1);
                  else
                    remap_field(field2, missval, gridsize2, remap.vars, field1, gradients);
                }
              else
                {
                  remap_field(remapSwitches, remap, field1, field2);
                }

              if (operfunc == REMAPSCON || operfunc == REMAPCON2 || operfunc == REMAPCON || operfunc == REMAPYCON2)
                {
                  // used only to check the result of remapcon
                  if (0) remap_normalize_field(remap.vars.normOpt, field2, remap.tgtGrid);

                  remap_set_frac_min(field2, &remap.tgtGrid);
                }

              if (operfunc == REMAPSCON || operfunc == REMAPCON2 || operfunc == REMAPCON || operfunc == REMAPYCON2)
                {
                  if (varList1[varID].name == "gridbox_area") scale_gridbox_area(field1, field2, remap.tgtGrid.cell_area);
                }

              // calculate some statistics
              if (Options::cdoVerbose) remap_stat(remapOrder, remap.srcGrid, remap.tgtGrid, remap.vars, field1, field2);

              if (gridInqType(gridID2) == GRID_GME) restore_gme_grid(field2, remap.tgtGrid);

              field2.nmiss = field_num_mv(field2);
            }

          cdo_def_record(streamID2, varID, levelID);
          cdo_write_record(streamID2, field2);
        }

      tsID++;
    }

  cdo_stream_close(streamID2);

WRITE_REMAP:

  if (writeRemapWeightsOnly)
    remap_write_data_scrip(cdo_get_stream_name(1), remapSwitches, remaps[remapIndex].srcGrid, remaps[remapIndex].tgtGrid,
                           remaps[remapIndex].vars);

  cdo_stream_close(streamID1);

  if (doRemap && remap_genweights && remaps[0].nused == 0)
    remap_print_warning(remapWeightsFile, operfunc, remaps[0].srcGrid, remaps[0].nmiss);

  for (remapIndex = 0; remapIndex < numRemaps; remapIndex++)
    {
      remap_vars_free(remaps[remapIndex].vars);
      remap_grid_free(remaps[remapIndex].srcGrid);
      remap_grid_free(remaps[remapIndex].tgtGrid);
      remap_search_free(remaps[remapIndex].search);
    }

  cdo_finish();

  return nullptr;
}
