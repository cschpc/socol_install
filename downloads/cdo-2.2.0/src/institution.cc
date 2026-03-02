/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include "process_int.h"
#include "readline.h"
#include "cdo_default_values.h"

static int
readInstitution(const char *instfile)
{
  char line[1024];
  int lnr = 0;
  int nvar = 0, maxvar = 4;
  char name[1024], longname[1024];
  int center = CDI_UNDEFID, subcenter = CDI_UNDEFID;

  auto instfp = std::fopen(instfile, "r");
  if (instfp == nullptr) return CDI_UNDEFID;

  while (cdo::readline(instfp, line, 1024))
    {
      lnr++;
      if (line[0] == '#') continue;
      if (nvar == maxvar) break;
      nvar++;

      char *pline = line;
      while (isspace((int) *pline)) pline++;

      if (nvar == 1) maxvar = isdigit((int) pline[0]) ? 4 : 2;

      if (nvar == 1 && maxvar == 4) center = atoi(pline);

      if (nvar == 2 && maxvar == 4)
        {
          if (!isdigit((int) pline[0])) cdo_abort("wrong format in line %d. Missing subcenter!", lnr);

          subcenter = atoi(pline);
        }

      if ((nvar == 3 && maxvar == 4) || (nvar == 1 && maxvar == 2)) strcpy(name, pline);

      if ((nvar == 4 && maxvar == 4) || (nvar == 2 && maxvar == 2)) strcpy(longname, pline);
    }

  std::fclose(instfp);

  auto instID = institutInq(center, subcenter, name, longname);
  if (instID == CDI_UNDEFID) instID = institutDef(center, subcenter, name, longname);

  return instID;
}

void
define_institution(const char *instarg)
{
  const char *instname = instarg;
  int instID = readInstitution(instname);

  if (instID == CDI_UNDEFID) instID = institutInq(0, 0, instname, nullptr);

  if (instID == CDI_UNDEFID) cdo_abort("institution <%s> not found", instname);

  CdoDefault::InstID = instID;
}
