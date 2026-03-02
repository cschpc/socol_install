#ifndef CDO_VLIST_H
#define CDO_VLIST_H

#include <string>
#include <map>

#include <cdi.h>

#include "cdo_varlist.h"
#include "cdo_options.h"
#include "varray.h"

enum cmp_flag
{
  CMP_NAME = 1,
  CMP_GRID = 2,
  CMP_NLEVEL = 4,
  CMP_GRIDSIZE = 8,
  CMP_HRD = CMP_NAME | CMP_GRIDSIZE,
  CMP_DIM = CMP_GRIDSIZE | CMP_NLEVEL | CMP_GRID,
  CMP_ALL = CMP_NAME | CMP_GRIDSIZE | CMP_NLEVEL | CMP_GRID
};

void vlist_define_timestep_type(int vlistID, int operfunc);

void vlist_map(int vlistID1, int vlistID2, int flag, int mapflag, std::map<int, int> &mapOfVarIDs);
void vlist_compare(int vlistID1, int vlistID2, int flag);
int vlist_compare_x(int vlistID1, int vlistID2, int flag);
bool vlist_is_szipped(int vlistID);

int vlist_inq_nwpv(int vlistID, int varID);
size_t vlist_check_gridsize(int vlistID);
int vlist_get_psvarid(int vlistID, int zaxisID);
Varray<double> vlist_read_vct(int vlistID, int &zaxisID_ML, int &numHybridLevels, int &numFullLevels, int &numHalfLevels);
void vlist_change_hybrid_zaxis(int vlistID1, int vlistID2, int zaxisID1, int zaxisID2);

void cdo_compare_grids(int gridID1, int gridID2);  // TODO: Check if this belongs here or if it should be in griddes.*

int vlist_get_first_gaussian_grid(int vlistID);
int vlist_get_first_spectral_grid(int vlistID);
int vlist_get_first_fourier_grid(int vlistID);

void cdo_check_missval(double missval, const std::string &varname);

#endif
