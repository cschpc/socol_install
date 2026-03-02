/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#ifndef MPIM_GRID_H
#define MPIM_GRID_H

#include <cstdio>
#include <vector>
#include <string>

#include <cdi.h>

#include "grid_proj.h"
#include "grid_rot.h"
#include "grid_convert.h"
#include "varray.h"

enum class NeedCorners
{
  No,
  Yes,
  IfAvail
};

extern "C" int qu2reg3_double(double *pfield, int *kpoint, int klat, int klon, double msval, int *kret, int omisng, int operio,
                              int oveggy);

extern bool gridVerbose;

static inline bool
gridProjIsSupported(int gridID)
{
  auto gridtype = gridInqType(gridID);
  auto pt = (gridtype == GRID_PROJECTION) ? gridInqProjType(gridID) : -1;
  return (pt == CDI_PROJ_RLL || pt == CDI_PROJ_LCC || pt == CDI_PROJ_LAEA || pt == CDI_PROJ_STERE || pt == CDI_PROJ_SINU
          || pt == CDI_PROJ_HEALPIX);
}

static inline bool
gridHasCoordinates(int gridID)
{
  return (gridInqXvals(gridID, nullptr) && gridInqYvals(gridID, nullptr));
}

static inline bool
gridHasBounds(int gridID)
{
  return (gridInqXbounds(gridID, nullptr) && gridInqYbounds(gridID, nullptr));
}

void gridEnableVerbose(bool enable);

int nfc_to_nlat(int nfc, int ntr);
int nlat_to_ntr(int nlat);
int nlat_to_ntr_linear(int nlat);
int nlat_to_ntr_cubic(int nlat);
int ntr_to_nlat(int ntr);
int ntr_to_nlat_linear(int ntr);
int ntr_to_nlat_cubic(int ntr);
int nlat_to_nlon(int nlat);
int nlat_to_nlon_cubic(int nlat);

void grid_copy_names(int gridID1, int gridID2);
void grid_copy_mapping(int gridID1, int gridID2);

bool grid_is_distance_generic(int gridID);

void grid_to_radian(const std::string &units, size_t nvals, double *values, const std::string &description);
void cdo_grid_to_radian(int gridID, int varID, size_t nvals, double *values, const std::string &description);
void cdo_grid_to_degree(int gridID, int varID, size_t nvals, double *values, const std::string &description);

void grid_gen_corners(size_t n, const double *vals, double *corners);
void grid_gen_bounds(size_t n, const std::vector<double> &vals, std::vector<double> &bounds);
void grid_check_lat_borders(int n, double *ybounds);

void grid_gen_xbounds2D(size_t nx, size_t ny, const std::vector<double> &xbounds, std::vector<double> &xbounds2D);
void grid_gen_ybounds2D(size_t nx, size_t ny, const std::vector<double> &ybounds, std::vector<double> &ybounds2D);

int gridcell_weights(int gridID, Varray<double> &weights);
int gridGenArea(int gridID, double *area);
int gridGenAreaReg2Dweights(int gridID, double *area);

int gridToZonal(int gridID);
int gridToMeridional(int gridID);
int gridToUnstructured(int gridID, NeedCorners needCorners = NeedCorners::No);
int gridToUnstructuredSelecton(int gridID1, const std::vector<size_t> &selectionIndexList, int nocoords, int nobounds);
int gridToCurvilinear(int gridID, NeedCorners needCorners = NeedCorners::No);
int gridCurvilinearToRegular(int gridID);
int gridToRegular(int gridID);
void field2regular(int gridID1, int gridID2, double missval, double *array, size_t nmiss, int lnearest);

// GME grid
void gme_factorni(int kni, int *kni2, int *kni3);
void gme_grid(int withBounds, size_t gridsize, double *rlon, double *rlat, double *blon, double *blat, int *imask, int ni, int nd,
              int ni2, int ni3);

void cdo_print_griddes(int gridID, int opt);

bool grid_has_proj_params(int gridID);
std::vector<char> grid_get_proj_params(int gridID);

bool is_point_grid(int gridID);
int generate_full_point_grid(int gridID);
int generate_full_cell_grid(int gridID);
int generate_full_grid(int gridID);

static inline bool
is_reg2d_grid(int gridID)
{
  return (gridInqType(gridID) == GRID_LONLAT || gridInqType(gridID) == GRID_GAUSSIAN);
}

static inline bool
is_healpix_grid(int gridID)
{
  return (gridInqType(gridID) == GRID_PROJECTION && gridInqProjType(gridID) == CDI_PROJ_HEALPIX);
}

#endif /* MPIM_GRID_H */
