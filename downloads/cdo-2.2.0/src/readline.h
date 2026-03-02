/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/
#ifndef READLINE_H
#define READLINE_H

#include <cstdio>

namespace cdo
{

int readline(FILE *fp, char *line, int len);

}

#endif
