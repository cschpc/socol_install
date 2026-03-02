#ifndef CDO_VARLIST_H
#define CDO_VARLIST_H

#include <string>
#include <vector>

#include <cdi.h>

#include "cdo_varlist.h"
#include "cdo_options.h"

struct CdoVar
{
  std::string name;
  std::string longname;
  std::string units;
  MemType memType = MemType::Native;
  int gridID = -1;
  int zaxisID = -1;
  int timetype = -1;
  int tsteptype = -1;
  size_t gridsize = 0;
  int nlevels = 0;
  int datatype = -1;
  double missval = 0;
  int code = 0;
  int param = 0;
  int nwpv = 1;  // number of words per value; real:1  complex:2
  bool isConstant = false;
};

using VarList = std::vector<CdoVar>;
void varListInit(VarList &varList, int vlistID);
void varListSetUniqueMemtype(VarList &varList);
void varListSetMemtype(VarList &varList, MemType memType);
int varList_numConstVars(const VarList &varList);
int varList_numVaryingVars(const VarList &varList);

struct VarIDs
{
  int sgeopotID = CDI_UNDEFID;
  int geopotID = CDI_UNDEFID;
  int tempID = CDI_UNDEFID;
  int psID = CDI_UNDEFID;
  int lnpsID = CDI_UNDEFID;
  int lnpsID2 = CDI_UNDEFID;
  int gheightID = CDI_UNDEFID;
  int humID = CDI_UNDEFID;
  int clwcID = CDI_UNDEFID;
  int ciwcID = CDI_UNDEFID;
};

VarIDs search_varIDs(const VarList &varList, int vlistID, int numFullLevels);

#endif
