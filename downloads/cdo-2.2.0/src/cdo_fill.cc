/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include "cdo_fill.h"
#include "cdo_vlist.h"
#include "process_int.h"

void
cdo_fill_ts(int p_vlistID, Varray2D<double> &p_varData)
{
  const auto nvars = vlistNvars(p_vlistID);
  p_varData.resize(nvars);
  for (int varID = 0; varID < nvars; ++varID)
    {
      const auto gridsize = gridInqSize(vlistInqVarGrid(p_vlistID, varID));
      const auto nlev = zaxisInqSize(vlistInqVarZaxis(p_vlistID, varID));
      p_varData[varID].resize(nlev * gridsize);
    }
}

void
cdo_fill_ts(int p_vlistID, Varray2D<double> &p_varData, Varray2D<size_t> &p_varNmiss)
{
  const auto nvars = vlistNvars(p_vlistID);
  p_varData.resize(nvars);
  p_varNmiss.resize(nvars);
  for (int varID = 0; varID < nvars; ++varID)
    {
      const auto gridsize = gridInqSize(vlistInqVarGrid(p_vlistID, varID));
      const auto nlev = zaxisInqSize(vlistInqVarZaxis(p_vlistID, varID));
      p_varData[varID].resize(nlev * gridsize);
      p_varNmiss[varID].resize(nlev);
    }
}
