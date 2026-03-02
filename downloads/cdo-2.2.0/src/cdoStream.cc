/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <string>
#include <cassert>

#include "cdoStream.h"
#include "cdo_options.h"

static int nextCdoStreamID = 1;

CdoStream::~CdoStream() {}
CdoStream::CdoStream() : m_cdoStreamID(nextCdoStreamID++) {}

int
CdoStream::get_id()
{
  return m_cdoStreamID;
}

int
CdoStream::getTsID()
{
  return m_tsID;
}

int
CdoStream::getVarID()
{
  return m_varID;
}

int
CdoStream::getVlistID()
{
  return m_vlistID;
}

void
CdoStream::defDatarangeList(int p_vlistID)
{
  auto filetype = m_filetype;

  if (m_vlistID != -1) cdo_abort("Internal problem, vlist already defined!");

  if (m_datarangelist.size() != 0) cdo_abort("Internal problem, datarangelist already allocated!");

  auto nvars = vlistNvars(p_vlistID);
  assert(nvars > 0);

  m_datarangelist.resize(nvars);

  for (int varID = 0; varID < nvars; ++varID)
    {
      m_datarangelist[varID].gridsize = gridInqSize(vlistInqVarGrid(p_vlistID, varID));
      m_datarangelist[varID].datatype = vlistInqVarDatatype(p_vlistID, varID);
      m_datarangelist[varID].missval = vlistInqVarMissval(p_vlistID, varID);

      double addoffset = 0.0, scalefactor = 1.0;
      auto haveAddoffset = (cdiInqKeyFloat(p_vlistID, varID, CDI_KEY_ADDOFFSET, &addoffset) == CDI_NOERR);
      auto haveScalefactor = (cdiInqKeyFloat(p_vlistID, varID, CDI_KEY_SCALEFACTOR, &scalefactor) == CDI_NOERR);
      if (haveAddoffset) m_datarangelist[varID].addoffset = addoffset;
      if (haveScalefactor) m_datarangelist[varID].scalefactor = scalefactor;

      m_datarangelist[varID].checkDatarange = false;

      auto datatype = m_datarangelist[varID].datatype;

      if (filetype == CDI_FILETYPE_NC || filetype == CDI_FILETYPE_NC2 || filetype == CDI_FILETYPE_NC4
          || filetype == CDI_FILETYPE_NC4C || filetype == CDI_FILETYPE_NC5 || filetype == CDI_FILETYPE_NCZARR)
        {
          if (datatype == CDI_DATATYPE_UINT8
              && (filetype == CDI_FILETYPE_NC || filetype == CDI_FILETYPE_NC2 || filetype == CDI_FILETYPE_NC5))
            {
              datatype = CDI_DATATYPE_INT16;
              m_datarangelist[varID].datatype = datatype;
            }

          if (datatype == CDI_DATATYPE_UINT16
              && (filetype == CDI_FILETYPE_NC || filetype == CDI_FILETYPE_NC2 || filetype == CDI_FILETYPE_NC5))
            {
              datatype = CDI_DATATYPE_INT32;
              m_datarangelist[varID].datatype = datatype;
            }

          if (haveAddoffset || haveScalefactor)
            {
              if (datatype == CDI_DATATYPE_INT8 || datatype == CDI_DATATYPE_UINT8 || datatype == CDI_DATATYPE_INT16
                  || datatype == CDI_DATATYPE_UINT16)
                m_datarangelist[varID].checkDatarange = true;
            }
          else if (Options::CheckDatarange)
            {
              m_datarangelist[varID].checkDatarange = true;
            }
        }
    }

  m_vlistID = p_vlistID;  // used for -r/-a
}
