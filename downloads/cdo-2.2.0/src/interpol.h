/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/
#ifndef INTERPOL_H
#define INTERPOL_H

#include <stdio.h>

class Field;

void intgridbil(Field &field1, Field &field2);
void intgridcon(Field &field1, Field &field2);
void interpolate(Field &field1, Field &field2);
void intgriddis(Field &field1, Field &field2, size_t nnn);
void intgridnn(Field &field1, Field &field2);

constexpr double
intlin(const double x, const double y1, const double x1, const double y2, const double x2)
{
  // intlin - lineare interpolation

  // Uwe Schulzweida  04/05/1995

  return (y2 * (x - x1) + y1 * (x2 - x)) / (x2 - x1);
}

#endif
