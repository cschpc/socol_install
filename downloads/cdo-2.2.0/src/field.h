/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/
#ifndef FIELD_H
#define FIELD_H

#include <cstdio>
#include <utility>
#include "varray.h"
#include "cdo_options.h"
#include "cdo_vlist.h"
#include "compare.h"

// clang-format off
const auto memtype_is_float_float   = [](auto a, auto b) noexcept { return (a == MemType::Float  && b == MemType::Float); };
const auto memtype_is_double_float  = [](auto a, auto b) noexcept { return (a == MemType::Double && b == MemType::Float); };
const auto memtype_is_float_double  = [](auto a, auto b) noexcept { return (a == MemType::Float  && b == MemType::Double); };
const auto memtype_is_double_double = [](auto a, auto b) noexcept { return (a == MemType::Double && b == MemType::Double); };
// clang-format on

enum field_flag
{
  FIELD_VEC = 2,   // allocated memory
  FIELD_FLT = 4,   // 32-bit float
  FIELD_DBL = 8,   // 64-bit float
  FIELD_NAT = 16,  // native: 32-bit float for 32-bit float data, otherweise 64-bit float
};

// clang-format off
class  // Field
#ifdef WARN_UNUSED
[[gnu::warn_unused]]
#endif
Field
// clang-format on
{
public:
  int fpeRaised = 0;
  int nwpv = 1;  // number of words per value; real:1  complex:2
  int grid = -1;
  MemType memType = MemType::Native;  // MemType::Float or MemType::Double

  size_t gridsize = 0;
  size_t size = 0;
  size_t nsamp = 0;

  size_t nmiss = 0;
  double missval = 0;

  Varray<float> vec_f;
  Varray<double> vec_d;
  Varray<double> weightv;

  Field() {}
  void init(const CdoVar &var);
  void resize(size_t count);
  void resize(size_t count, double value);
  void resizef(size_t count);
  void resizef(size_t count, float value);
  bool empty() const;
  void check_gridsize() const;

  bool
  hasData() const
  {
    return (memType == MemType::Float) ? !vec_f.empty() : !vec_d.empty();
  }

private:
  size_t m_count = 0;
};

// clang-format off
class  // Field3D
#ifdef WARN_UNUSED
[[gnu::warn_unused]]
#endif
Field3D : public Field
// clang-format on
{
public:
  size_t nlevels = 0;

  Field3D() {}
  void init(const CdoVar &var);
};

struct RecordInfo
{
  int varID = 0;
  int levelID = 0;
  void set(int _varID, int _levelID) { this->varID = _varID; this->levelID = _levelID; }
  std::pair<int, int> get() const { return std::make_pair(varID, levelID); }
};

using FieldVector = std::vector<Field>;
using FieldVector2D = std::vector<FieldVector>;
using FieldVector3D = std::vector<FieldVector2D>;

using Field3DVector = std::vector<Field3D>;

void field_fill(Field &field, double value);
void field_ncopy(size_t n, const Field &fieldIn, Field &fieldOut);
void field_copy(const Field &fieldIn, Field &fieldOut);
void field_copy(const Field3D &fieldIn, Field3D &fieldOut);
void field_copy(const Field3D &fieldIn, int levelID, Field &fieldOut);
void field_add(Field &field1, const Field3D &field2, int levelID);
size_t field_num_miss(const Field &field);
size_t field_num_mv(Field &field);
MinMax field_min_max(const Field &field);

template <class UnaryOperation>
void
field_transform(const Field &fieldIn, Field &fieldOut, UnaryOperation unary_op)
{
  if (fieldIn.memType == MemType::Float && fieldOut.memType == MemType::Float)
    varray_transform(fieldIn.vec_f, fieldOut.vec_f, unary_op);
  else
    varray_transform(fieldIn.vec_d, fieldOut.vec_d, unary_op);
}

#endif /* FIELD_H */
