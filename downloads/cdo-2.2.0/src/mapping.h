#ifndef MAPVAR_H
#define MAPVAR_H

#include <string>
#include "cdo_cmor.h"
#include "pmlist.h"

void mapvar(int vlistID, int varID, const KeyValues &kv, const std::string &key, CmorVar *var, bool &hasValidMin, bool &hasValidMax,
            int ptab, bool isnPtmodeName);

#endif
