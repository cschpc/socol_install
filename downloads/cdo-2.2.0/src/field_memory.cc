/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include "field.h"
#include "cdo_vlist.h"
#include "cdo_options.h"

static void
fields_from_vlist_kernel(const int vlistID, FieldVector2D &field2D, int ptype, bool lfill, double fillValue)
{
  const auto allocateData = (ptype & FIELD_VEC);
  const auto nvars = vlistNvars(vlistID);
  field2D.resize(nvars);

  for (int varID = 0; varID < nvars; ++varID)
    {
      const auto gridID = vlistInqVarGrid(vlistID, varID);
      const auto gridsize = gridInqSize(gridID);
      const size_t nwpv = vlist_inq_nwpv(vlistID, varID);  // number of words per value; real:1  complex:2
      const auto size = gridsize * nwpv;
      const auto zaxisID = vlistInqVarZaxis(vlistID, varID);
      const auto nlevels = zaxisInqSize(zaxisID);
      const auto missval = vlistInqVarMissval(vlistID, varID);
      const auto datatype = vlistInqVarDatatype(vlistID, varID);
      auto memType = (ptype & FIELD_FLT) ? MemType::Float : MemType::Double;
      if (ptype & FIELD_NAT)
        {
          if (Options::CDO_Memtype == MemType::Native)
            memType = (datatype == CDI_DATATYPE_FLT32 || datatype == CDI_DATATYPE_CPX32) ? MemType::Float : MemType::Double;
          else
            memType = Options::CDO_Memtype;
        }

      field2D[varID].resize(nlevels);

      for (int levelID = 0; levelID < nlevels; ++levelID)
        {
          auto &field = field2D[varID][levelID];
          field.nwpv = nwpv;
          field.grid = gridID;
          field.size = size;
          field.memType = memType;
          field.missval = missval;

          if (allocateData)
            {
              if (memType == MemType::Float)
                {
                  if (lfill)
                    field.resizef(size, (float) fillValue);
                  else
                    field.resizef(size);
                }
              else
                {
                  if (lfill)
                    field.resize(size, fillValue);
                  else
                    field.resize(size);
                }
            }
        }
    }
}

void
fields_from_vlist(const int vlistID, FieldVector2D &field2D)
{
  fields_from_vlist_kernel(vlistID, field2D, 0, false, 0);
}

void
fields_from_vlist(const int vlistID, FieldVector2D &field2D, int ptype)
{
  fields_from_vlist_kernel(vlistID, field2D, ptype, false, 0);
}

void
fields_from_vlist(const int vlistID, FieldVector2D &field2D, int ptype, double fillValue)
{
  fields_from_vlist_kernel(vlistID, field2D, ptype, true, fillValue);
}
