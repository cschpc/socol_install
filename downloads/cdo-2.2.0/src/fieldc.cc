/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include "cdo_output.h"
#include "field_functions.h"

template <typename T>
static void
fieldc_mul(size_t len, size_t nmiss, Varray<T> &v, const T missval, const double rconst)
{
  const auto missval1 = missval;
  const auto missval2 = missval;

  if (nmiss)
    {
      for (size_t i = 0; i < len; ++i) v[i] = MULMN(v[i], rconst);
    }
  else
    {
      for (size_t i = 0; i < len; ++i) v[i] *= rconst;
    }
}

void
fieldc_mul(Field &field, const double rconst)
{
  if (field.memType == MemType::Float)
    fieldc_mul(field.size, field.nmiss, field.vec_f, (float) field.missval, rconst);
  else
    fieldc_mul(field.size, field.nmiss, field.vec_d, field.missval, rconst);
}

template <typename T>
static size_t
fieldc_div(size_t len, size_t nmiss, Varray<T> &v, const T missval, const double rconst)
{
  const auto missval1 = missval;
  const auto missval2 = missval;

  if (nmiss || IS_EQUAL(rconst, 0))
    {
      for (size_t i = 0; i < len; ++i) v[i] = DIVMN(v[i], rconst);

      if (IS_EQUAL(rconst, 0)) nmiss = len;
    }
  else
    {
      for (size_t i = 0; i < len; ++i) v[i] /= rconst;
    }

  return nmiss;
}

void
fieldc_div(Field &field, const double rconst)
{
  if (field.memType == MemType::Float)
    field.nmiss = fieldc_div(field.size, field.nmiss, field.vec_f, (float) field.missval, rconst);
  else
    field.nmiss = fieldc_div(field.size, field.nmiss, field.vec_d, field.missval, rconst);
}

template <typename T>
static void
fieldc_add(size_t len, size_t nmiss, Varray<T> &v, const T missval, const double rconst)
{
  const auto missval1 = missval;
  const auto missval2 = missval;

  if (nmiss)
    {
      for (size_t i = 0; i < len; ++i) v[i] = ADDMN(v[i], rconst);
    }
  else
    {
      for (size_t i = 0; i < len; ++i) v[i] += rconst;
    }
}

void
fieldc_add(Field &field, const double rconst)
{
  if (field.memType == MemType::Float)
    fieldc_add(field.size, field.nmiss, field.vec_f, (float) field.missval, rconst);
  else
    fieldc_add(field.size, field.nmiss, field.vec_d, field.missval, rconst);
}

void
fieldc_sub(Field &field, const double rconst)
{
  fieldc_add(field, -rconst);
}

template <typename T>
static void
fieldc_min(size_t len, size_t nmiss, Varray<T> &v, const T missval, const T rconst)
{
  const auto missval1 = missval;

  if (nmiss)
    {
      for (size_t i = 0; i < len; ++i)
        if (DBL_IS_EQUAL(v[i], missval1) || v[i] > rconst) v[i] = rconst;
    }
  else
    {
      for (size_t i = 0; i < len; ++i)
        if (v[i] > rconst) v[i] = rconst;
    }
}

void
fieldc_min(Field &field, const double rconst)
{
  if (field.memType == MemType::Float)
    fieldc_min(field.size, field.nmiss, field.vec_f, (float) field.missval, (float) rconst);
  else
    fieldc_min(field.size, field.nmiss, field.vec_d, field.missval, rconst);
}

template <typename T>
static void
fieldc_max(size_t len, size_t nmiss, Varray<T> &v, const T missval, const T rconst)
{
  const auto missval1 = missval;

  if (nmiss)
    {
      for (size_t i = 0; i < len; ++i)
        if (DBL_IS_EQUAL(v[i], missval1) || v[i] < rconst) v[i] = rconst;
    }
  else
    {
      for (size_t i = 0; i < len; ++i)
        if (v[i] < rconst) v[i] = rconst;
    }
}

void
fieldc_max(Field &field, const double rconst)
{
  if (field.memType == MemType::Float)
    fieldc_max(field.size, field.nmiss, field.vec_f, (float) field.missval, (float) rconst);
  else
    fieldc_max(field.size, field.nmiss, field.vec_d, field.missval, rconst);
}

template <typename T>
static void
fieldc_mod(size_t len, size_t nmiss, Varray<T> &v, const T missval, const double divisor)
{
  (void) nmiss;
  const auto missval1 = missval;

  for (size_t i = 0; i < len; ++i) { v[i] = DBL_IS_EQUAL(v[i], missval1) ? missval1 : std::fmod((double) v[i], divisor); }
}

void
fieldc_mod(Field &field, const double divisor)
{
  if (field.memType == MemType::Float)
    fieldc_mod(field.size, field.nmiss, field.vec_f, (float) field.missval, divisor);
  else
    fieldc_mod(field.size, field.nmiss, field.vec_d, field.missval, divisor);
}

void
fieldc_function(Field &field, const double rconst, const int function)
{
  switch (function)
    {
    case FieldFunc_Add: fieldc_add(field, rconst); break;
    case FieldFunc_Sub: fieldc_sub(field, rconst); break;
    case FieldFunc_Mul: fieldc_mul(field, rconst); break;
    case FieldFunc_Div: fieldc_div(field, rconst); break;
    case FieldFunc_Min: fieldc_min(field, rconst); break;
    case FieldFunc_Max: fieldc_max(field, rconst); break;
    case FieldFunc_Mod: fieldc_mod(field, rconst); break;
    default: cdo_abort("%s: function %d not implemented!", __func__, function);
    }
}
