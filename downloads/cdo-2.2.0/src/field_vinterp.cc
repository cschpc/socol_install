/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include "field_vinterp.h"

void
gen_vert_index(std::vector<int> &vertIndex, Varray<double> &plev, Field3D &full_level, size_t gridsize, bool lreverse)
{
  const auto nplev = plev.size();
  const auto nhlevf = full_level.nlevels;
  if (full_level.memType == MemType::Float)
    gen_vert_index(vertIndex.data(), plev.data(), full_level.vec_f.data(), gridsize, nplev, nhlevf, lreverse);
  else
    gen_vert_index(vertIndex.data(), plev.data(), full_level.vec_d.data(), gridsize, nplev, nhlevf, lreverse);
}

void
gen_vert_index_mv(std::vector<int> &vertIndex, Varray<double> &plev, size_t gridsize, Field &level0, Varray<size_t> &pnmiss,
                  bool lreverse)
{
  const auto nplev = plev.size();
  if (level0.memType == MemType::Float)
    gen_vert_index_mv(vertIndex.data(), plev.data(), gridsize, nplev, level0.vec_f.data(), pnmiss.data(), lreverse);
  else
    gen_vert_index_mv(vertIndex.data(), plev.data(), gridsize, nplev, level0.vec_d.data(), pnmiss.data(), lreverse);
}

void
vertical_interp_T(size_t nlevels, Field3D &full_level, Field3D &half_level, Field3D &field1, Field3D &field2, Field &sgeopot,
                  std::vector<int> &vertIndex, Varray<double> &plev, size_t gridsize)
{
  const auto nplev = plev.size();
  const auto missval = field1.missval;
  if (field1.memType == MemType::Float)
    vertical_interp_T(sgeopot.vec_f.data(), field1.vec_f.data(), field2.vec_f.data(), full_level.vec_f.data(),
                      half_level.vec_f.data(), &vertIndex[0], plev.data(), nplev, gridsize, nlevels, missval);
  else
    vertical_interp_T(sgeopot.vec_d.data(), field1.vec_d.data(), field2.vec_d.data(), full_level.vec_d.data(),
                      half_level.vec_d.data(), &vertIndex[0], plev.data(), nplev, gridsize, nlevels, missval);
}

void
vertical_interp_Z(size_t nlevels, Field3D &full_level, Field3D &half_level, Field3D &field1, Field3D &field2, Field3D &temp,
                  Field &sgeopot, std::vector<int> &vertIndex, Varray<double> &plev, size_t gridsize)
{
  const auto nplev = plev.size();
  const auto missval = field1.missval;
  if (field1.memType == MemType::Float)
    vertical_interp_Z(sgeopot.vec_f.data(), field1.vec_f.data(), field2.vec_f.data(), full_level.vec_f.data(),
                      half_level.vec_f.data(), &vertIndex[0], temp.vec_f.data(), plev.data(), nplev, gridsize, nlevels, missval);
  else
    vertical_interp_Z(sgeopot.vec_d.data(), field1.vec_d.data(), field2.vec_d.data(), full_level.vec_d.data(),
                      half_level.vec_d.data(), &vertIndex[0], temp.vec_d.data(), plev.data(), nplev, gridsize, nlevels, missval);
}

void
vertical_interp_X(const Field3D &levels3D, const Field3D &field1, Field3D &field2, const std::vector<int> &vertIndex3D,
                  const Varray<double> &levels2, size_t gridsize)
{
  const auto numLevels2 = levels2.size();
  const auto missval = field1.missval;
  if (field1.memType == MemType::Float)
    vertical_interp_X(field1.vec_f.data(), field2.vec_f.data(), levels3D.vec_f.data(), vertIndex3D.data(), levels2.data(),
                      numLevels2, gridsize, levels3D.nlevels, missval);
  else
    vertical_interp_X(field1.vec_d.data(), field2.vec_d.data(), levels3D.vec_d.data(), vertIndex3D.data(), levels2.data(),
                      numLevels2, gridsize, levels3D.nlevels, missval);
}
