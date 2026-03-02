/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/
#ifndef REMAP_H
#define REMAP_H

#include <cstdint>
#include <cmath>

#include "varray.h"
#include "remap_vars.h"
#include "remap_grid_cell_search.h"
#include "grid_point_search.h"
#include "mpim_grid/grid_healpix.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

constexpr double PI = M_PI;
constexpr double PI2 = (2.0 * PI);
constexpr double PIH = (0.5 * PI);
constexpr float PI_f = PI;
constexpr float PI2_f = PI2;
constexpr float PIH_f = PIH;

constexpr double TINY = 1.e-14;

enum class RemapGridType
{
  Undefined,
  HealPix,
  Reg2D,
  Unstruct
};

#define REMAP_GRID_BASIS_SRC 1
#define REMAP_GRID_BASIS_TGT 2

struct LonLatPoint
{
  double lon = 0.0, lat = 0.0;
  LonLatPoint(){};
  LonLatPoint(double _lon, double _lat) : lon(_lon), lat(_lat){};
};

// clang-format off
struct  // RemapGrid
#ifdef WARN_UNUSED
[[gnu::warn_unused]]
#endif
RemapGrid
// clang-format on
{
  RemapGridType type = RemapGridType::Undefined;
  int gridID = -1;
  int tmpgridID = -1;
  int rank = 0;                 // rank of the grid
  size_t size = 0;              // total points on the grid
  size_t num_cell_corners = 0;  // number of corners for each grid cell

  int nside = 0;
  HpOrder order = HpOrder::Undef;

  bool needCellCorners = false;
  bool useCellCorners = false;  // use corners for bounding boxes

  bool doExtrapolate = false;
  bool isCyclic = false;

  size_t dims[2] = { 0, 0 };  // size of grid dimension

  int nvgp = 0;        // size of vgpm
  Varray<int> vgpm;    // flag which cells are valid
  Varray<short> mask;  // flag which cells participate

  Varray<double> reg2d_center_lon;  // reg2d lon/lat coordinates for
  Varray<double> reg2d_center_lat;  // each grid center in radians
  Varray<double> reg2d_corner_lon;  // reg2d lon/lat coordinates for
  Varray<double> reg2d_corner_lat;  // each grid corner in radians

  Varray<double> cell_center_lon;  // lon/lat coordinates for
  Varray<double> cell_center_lat;  // each grid center in radians
  Varray<double> cell_corner_lon;  // lon/lat coordinates for
  Varray<double> cell_corner_lat;  // each grid corner in radians

  Varray<double> cell_area;  // total area of each grid cell
  Varray<double> cell_frac;  // fractional area of grid cells participating in remapping
};

// clang-format off
struct  // GridSearchBins
#ifdef WARN_UNUSED
[[gnu::warn_unused]]
#endif
GridSearchBins
// clang-format on
{
  unsigned nbins;                // num of bins for restricted search
  size_t ncells;                 // total number of grid cells (cell_bound_box)
  Varray<size_t> bin_addr;       // min,max adds for grid cells in this lat bin
  Varray<float> bin_lats;        // min,max latitude for each search bin
  Varray<float> cell_bound_box;  // lon/lat bounding box for use
};

// clang-format off
struct  // RemapSearch
#ifdef WARN_UNUSED
[[gnu::warn_unused]]
#endif
RemapSearch
// clang-format on
{
  RemapGrid *srcGrid;
  RemapGrid *tgtGrid;

  GridSearchBins srcBins;
  GridSearchBins tgtBins;

  GridPointSearch gps;
  GridCellSearch gcs;
};

// clang-format off
struct  // RemapType
#ifdef WARN_UNUSED
[[gnu::warn_unused]]
#endif
RemapType
// clang-format on
{
  int nused = 0;
  int gridID = -1;
  size_t gridsize = 0;
  size_t nmiss = 0;
  RemapGrid srcGrid;
  RemapGrid tgtGrid;
  RemapVars vars;
  RemapSearch search;
};

#define REMAP_WRITE_REMAP 2
#define REMAP_MAX_ITER 3
#define REMAP_NUM_SRCH_BINS 4
#define REMAP_GENWEIGHTS 5

int remap_check_mask_indices(const size_t (&indices)[4], const Varray<short> &mask);

void remap_set_threshhold(double threshhold);
void remap_set_int(int remapvar, int value);

void remap_init_grids(RemapMethod mapType, bool doExtrapolate, int gridID1, RemapGrid &srcGrid, int gridID2, RemapGrid &tgtGrid);

void remap_grid_free(RemapGrid &grid);
void remap_grid_alloc(RemapMethod mapType, RemapGrid &grid);
void remap_search_init(RemapMethod mapType, RemapSearch &search, RemapGrid &srcGrid, RemapGrid &tgtGrid);
void remap_search_free(RemapSearch &search);

void remap_search_points(RemapSearch &rsearch, const LonLatPoint &llpoint, knnWeightsType &knnWeights);
int remap_search_square(RemapSearch &rsearch, const LonLatPoint &llpoint, size_t (&src_add)[4], double (&srcLats)[4],
                        double (&srcLons)[4]);
size_t remap_search_cells(RemapSearch &rsearch, bool isReg2dCell, GridCell &gridCell, Varray<size_t> &srchAddr);

void remap_bilinear_weights(RemapSearch &rsearch, RemapVars &rv);
void remap_bicubic_weights(RemapSearch &rsearch, RemapVars &rv);
void remap_distwgt_weights(size_t numNeighbors, RemapSearch &rsearch, RemapVars &rv);
void remap_conserv_weights(RemapSearch &rsearch, RemapVars &rv);
void remap_conserv_weights_scrip(RemapSearch &rsearch, RemapVars &rv);

void remap_bilinear(RemapSearch &rsearch, const Field &field1, Field &field2);
void remap_bicubic(RemapSearch &rsearch, const Field &field1, Field &field2);
void remap_dist_wgt(size_t numNeighbors, RemapSearch &rsearch, const Field &field1, Field &field2);
void remap_conserv(NormOpt normOpt, RemapSearch &rsearch, const Field &field1, Field &field2);

void remap_stat(int remapOrder, RemapGrid &srcGrid, RemapGrid &tgtGrid, RemapVars &rv, const Field &field1, const Field &field2);

template <typename T>
void remap_gradients(RemapGrid &grid, const Varray<short> &mask, const Varray<T> &array, RemapGradients &gradients);
void remap_gradients(RemapGrid &grid, const Field &field, RemapGradients &gradients);

void sort_add(size_t numLinks, size_t num_wts, size_t *add1, size_t *add2, double *weights);
void sort_iter(size_t numLinks, size_t num_wts, size_t *add1, size_t *add2, double *weights, int parent);

void remap_write_data_scrip(const char *weightsfile, const RemapSwitches &remapSwitches, RemapGrid &srcGrid, RemapGrid &tgtGrid,
                            RemapVars &rv);
RemapSwitches remap_read_data_scrip(const std::string &weightsfile, int gridID1, int gridID2, RemapGrid &srcGrid,
                                    RemapGrid &tgtGrid, RemapVars &rv);

void calc_lat_bins(GridSearchBins &searchBins);
size_t get_srch_cells(size_t tgtCellIndex, GridSearchBins &tgtBins, GridSearchBins &srcBins, float *tgt_cell_bound_box,
                      Varray<size_t> &srch_add);

int grid_search_square_reg_2d_NN(size_t nx, size_t ny, size_t *nbr_add, double *nbr_dist, double plat, double plon,
                                 const Varray<double> &src_center_lat, const Varray<double> &src_center_lon);

int grid_search_square_reg_2d(RemapGrid *srcGrid, size_t (&src_add)[4], double (&srcLats)[4], double (&srcLons)[4], double plat,
                              double plon);

bool point_in_quad(bool isCyclic, size_t nx, size_t ny, size_t i, size_t j, size_t adds[4], double lons[4], double lats[4],
                   double plon, double plat, const double *centerLon, const double *centerLat);

int grid_search_square_curv_2d_scrip(RemapGrid *srcGrid, size_t (&src_add)[4], double (&srcLats)[4], double (&srcLons)[4],
                                     double plat, double plon, GridSearchBins &srcBins);

std::pair<double, double> remap_find_weights(const LonLatPoint &llpoint, const double (&srcLons)[4], const double (&srcLats)[4]);

int rect_grid_search(size_t &ii, size_t &jj, double x, double y, size_t nxm, size_t nym, const Varray<double> &xm,
                     const Varray<double> &ym);

LonLatPoint remapgrid_get_lonlat(RemapGrid *grid, size_t index);

void remap_check_area(size_t grid_size, const Varray<double> &cell_area, const char *name);

template <typename T>
void remap_set_mask(size_t gridsize, const Varray<T> &v, T missval, Varray<short> &mask);

#endif /* REMAP_H */
