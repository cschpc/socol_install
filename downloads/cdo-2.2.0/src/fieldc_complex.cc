/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include "process_int.h"
#include "field_functions.h"

template <typename T>
static void
fieldc_mul_complex(size_t len, size_t nmiss, Varray<T> &v, const double missval, const double rconst[2])
{
  (void) nmiss;
  // z1 x z2 = (x1x2 - y1y2) + i(x1y2 + x2y1)
  const auto missval1 = missval;
  const auto missval2 = missval;
  for (size_t i = 0; i < len; ++i)
    {
      const double v1r = v[2 * i];
      const double v1i = v[2 * i + 1];
      v[2 * i] = SUBMN(MULMN(v1r, rconst[0]), MULMN(v1i, rconst[1]));
      v[2 * i + 1] = ADDMN(MULMN(v1r, rconst[1]), MULMN(v1i, rconst[0]));
    }
}

static void
fieldc_mul_complex(Field &field, const double rconst[2])
{
  if (field.memType == MemType::Float)
    fieldc_mul_complex(field.size, field.nmiss, field.vec_f, field.missval, rconst);
  else
    fieldc_mul_complex(field.size, field.nmiss, field.vec_d, field.missval, rconst);
}

template <typename T>
static void
fieldc_div_complex(size_t len, size_t nmiss, Varray<T> &v, const double missval, const double rconst[2])
{
  (void) nmiss;
  // z1 / z2 = (x1x2 + y1y2) / (x2x2 + y2y2) + i (y1x2 - x1y2) / (x2x2 + y2y2)
  const auto missval1 = missval;
  const auto missval2 = missval;
  for (size_t i = 0; i < len; ++i)
    {
      const double v1r = v[2 * i];
      const double v1i = v[2 * i + 1];
      const double denominator = ADDMN(MULMN(rconst[0], rconst[0]), MULMN(rconst[1], rconst[1]));
      v[2 * i] = DIVMN(ADDMN(MULMN(v1r, rconst[0]), MULMN(v1i, rconst[1])), denominator);
      v[2 * i + 1] = DIVMN(SUBMN(MULMN(v1i, rconst[0]), MULMN(v1r, rconst[1])), denominator);
    }
}

static void
fieldc_div_complex(Field &field, const double rconst[2])
{
  if (field.memType == MemType::Float)
    fieldc_div_complex(field.size, field.nmiss, field.vec_f, field.missval, rconst);
  else
    fieldc_div_complex(field.size, field.nmiss, field.vec_d, field.missval, rconst);
}

template <typename T>
static void
fieldc_add_complex(size_t len, size_t nmiss, Varray<T> &v, const double missval, const double rconst[2])
{
  (void) nmiss;
  const auto missval1 = missval;
  const auto missval2 = missval;
  for (size_t i = 0; i < len; ++i)
    {
      v[2 * i] = ADDMN(v[2 * i], rconst[0]);
      v[2 * i + 1] = ADDMN(v[2 * i + 1], rconst[1]);
    }
}

static void
fieldc_add_complex(Field &field, const double rconst[2])
{
  if (field.memType == MemType::Float)
    fieldc_add_complex(field.size, field.nmiss, field.vec_f, field.missval, rconst);
  else
    fieldc_add_complex(field.size, field.nmiss, field.vec_d, field.missval, rconst);
}

template <typename T>
static void
fieldc_sub_complex(size_t len, size_t nmiss, Varray<T> &v, const double missval, const double rconst[2])
{
  (void) nmiss;
  const auto missval1 = missval;
  const auto missval2 = missval;
  for (size_t i = 0; i < len; ++i)
    {
      v[2 * i] = SUBMN(v[2 * i], rconst[0]);
      v[2 * i + 1] = SUBMN(v[2 * i + 1], rconst[1]);
    }
}

static void
fieldc_sub_complex(Field &field, const double rconst[2])
{
  if (field.memType == MemType::Float)
    fieldc_sub_complex(field.size, field.nmiss, field.vec_f, field.missval, rconst);
  else
    fieldc_sub_complex(field.size, field.nmiss, field.vec_d, field.missval, rconst);
}

void
fieldc_function_complex(Field &field, const double rconst[2], int function)
{
  switch (function)
    {
    case FieldFunc_Add: fieldc_add_complex(field, rconst); break;
    case FieldFunc_Sub: fieldc_sub_complex(field, rconst); break;
    case FieldFunc_Mul: fieldc_mul_complex(field, rconst); break;
    case FieldFunc_Div: fieldc_div_complex(field, rconst); break;
    default: cdo_abort("%s: function %d not implemented!", __func__, function);
    }
}
