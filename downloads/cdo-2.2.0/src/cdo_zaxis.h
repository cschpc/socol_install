#ifndef CDO_ZAXIS_H
#define CDO_ZAXIS_H

#include <string>
#include "varray.h"

int cdo_define_zaxis(const std::string &zaxisfile);
void define_zaxis(const char *zaxisarg);
int zaxis_from_name(const char *zaxisname);
int zaxis_from_file(FILE *zfp, const char *filename);
int zaxis_to_ltype(int zaxisID);
double cdo_zaxis_inq_level(int zaxisID, int levelID);
int cdo_zaxis_inq_levels(int zaxisID, double *levels);

void gen_layer_bounds(int nlev, const Varray<double> &levels, Varray<double> &lbounds, Varray<double> &ubounds);
int get_layer_thickness(bool useWeights, bool genBounds, int index, int zaxisID, int nlev, Varray<double> &thickness,
                        Varray<double> &weights);

static inline bool
positive_is_down(int zaxisID)
{
  return (zaxisInqPositive(zaxisID) == 2);
}

#endif
