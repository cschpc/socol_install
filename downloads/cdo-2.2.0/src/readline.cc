/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include "readline.h"

namespace cdo
{

int
readline(FILE *fp, char *line, int len)
{
  int ichar, ipos = 0;

  while ((ichar = fgetc(fp)) != EOF)
    {
      if (ichar == '\r')
        {
          if ((ichar = fgetc(fp)) != EOF)
            if (ichar != '\n') ungetc(ichar, fp);
          break;
        }
      if (ichar == '\n') break;
      line[ipos++] = ichar;
      if (ipos >= len)
        {
          fprintf(stderr, "readline Warning: end of line not found (maxlen = %d)!\n", len);
          break;
        }
    }
  line[ipos] = 0;

  if (feof(fp) && ipos == 0) return 0;

  return 1;
}

}  // namespace cdo
