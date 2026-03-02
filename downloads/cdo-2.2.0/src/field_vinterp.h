/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/
#ifndef FIELD_VINTERP_H
#define FIELD_VINTERP_H

#include "field.h"
#include "vertical_interp.h"

void gen_vert_index(std::vector<int> &vertIndex, Varray<double> &plev, Field3D &full_level, size_t gridsize, bool lreverse = false);

void gen_vert_index_mv(std::vector<int> &vertIndex, Varray<double> &plev, size_t gridsize, Field &level0, Varray<size_t> &pnmiss,
                       bool lreverse = false);

void vertical_interp_T(size_t nlevels, Field3D &full_level, Field3D &half_level, Field3D &field1, Field3D &field2, Field &sgeopot,
                       std::vector<int> &vertIndex, Varray<double> &plev, size_t gridsize);

void vertical_interp_Z(size_t nlevels, Field3D &full_level, Field3D &half_level, Field3D &field1, Field3D &field2, Field3D &temp,
                       Field &sgeopot, std::vector<int> &vertIndex, Varray<double> &plev, size_t gridsize);

void vertical_interp_X(const Field3D &levels3D, const Field3D &field1, Field3D &field2, const std::vector<int> &vertIndex3D,
                       const Varray<double> &levels2, size_t gridsize);

#endif
