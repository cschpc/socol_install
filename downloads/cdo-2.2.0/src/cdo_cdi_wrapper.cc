#include <cdi.h>
#include <cstring>
#include <cstdlib>
#include <vector>

#include "cdi_int.h"
#include "cdo_default_values.h"
#include "cdo_cdi_wrapper.h"
#include "cdo_output.h"

namespace cdo
{

const char *
filetype_to_cstr(int filetype)
{
  switch (filetype)
    {
    // clang-format off
    case CDI_FILETYPE_GRB:    return "GRIB";
    case CDI_FILETYPE_GRB2:   return "GRIB2";
    case CDI_FILETYPE_NC:     return "NetCDF";
    case CDI_FILETYPE_NC2:    return "NetCDF2";
    case CDI_FILETYPE_NC4:    return "NetCDF4";
    case CDI_FILETYPE_NC4C:   return "NetCDF4 classic";
    case CDI_FILETYPE_NC5:    return "NetCDF5";
    case CDI_FILETYPE_NCZARR: return "NCZarr";
    case CDI_FILETYPE_SRV:    return "SERVICE";
    case CDI_FILETYPE_EXT:    return "EXTRA";
    case CDI_FILETYPE_IEG:    return "IEG";
    default:                  return "";
    // clang-format on
    }
}

const char *
datatype_to_cstr(int datatype)
{
  static char cstr[20] = { 0 };
  if (datatype > 0 && datatype <= 32) std::snprintf(cstr, sizeof(cstr), "P%d", datatype);

  // clang-format off
  if      (datatype == CDI_DATATYPE_PACK  ) return "P0";
  else if (datatype > 0 && datatype <= 32 ) return cstr;
  else if (datatype == CDI_DATATYPE_CPX32 ) return "C32";
  else if (datatype == CDI_DATATYPE_CPX64 ) return "C64";
  else if (datatype == CDI_DATATYPE_FLT32 ) return "F32";
  else if (datatype == CDI_DATATYPE_FLT64 ) return "F64";
  else if (datatype == CDI_DATATYPE_INT8  ) return "I8";
  else if (datatype == CDI_DATATYPE_INT16 ) return "I16";
  else if (datatype == CDI_DATATYPE_INT32 ) return "I32";
  else if (datatype == CDI_DATATYPE_UINT8 ) return "U8";
  else if (datatype == CDI_DATATYPE_UINT16) return "U16";
  else if (datatype == CDI_DATATYPE_UINT32) return "U32";
  else                                      return "";
  // clang-format on
}

int
str_to_datatype(const std::string &datatypeStr)
{
  if (datatypeStr.size() > 1)
    {
      auto ilen = atoi(datatypeStr.c_str() + 1);
      // clang-format off
      if      (datatypeStr == "P0")     return CDI_DATATYPE_PACK;
      else if (datatypeStr[0] == 'P' && ilen > 0 && ilen <= 32) return ilen;
      else if (datatypeStr == "C32")    return CDI_DATATYPE_CPX32;
      else if (datatypeStr == "C64")    return CDI_DATATYPE_CPX64;
      else if (datatypeStr == "F32")    return CDI_DATATYPE_FLT32;
      else if (datatypeStr == "F64")    return CDI_DATATYPE_FLT64;
      else if (datatypeStr == "I8")     return CDI_DATATYPE_INT8;
      else if (datatypeStr == "I16")    return CDI_DATATYPE_INT16;
      else if (datatypeStr == "I32")    return CDI_DATATYPE_INT32;
      else if (datatypeStr == "U8")     return CDI_DATATYPE_UINT8;
      else if (datatypeStr == "U16")    return CDI_DATATYPE_UINT16;
      else if (datatypeStr == "U32")    return CDI_DATATYPE_UINT32;
      else if (datatypeStr == "real")   return CDI_DATATYPE_FLT32;
      else if (datatypeStr == "double") return CDI_DATATYPE_FLT64;
      // clang-format on
    }

  return -1;
}

}  // namespace cdo

int
cdo_taxis_create(int taxisType)
{
  if (CdoDefault::TaxisType != CDI_UNDEFID) taxisType = CdoDefault::TaxisType;
  return taxisCreate(taxisType);
}

void
cdo_taxis_copy_timestep(int taxisIDdes, int taxisIDsrc)
{
  taxisCopyTimestep(taxisIDdes, taxisIDsrc);
}

void
cdo_def_table_id(int tableID)
{
  cdiDefTableID(tableID);
}

void
grid_gen_xvals(int xsize, double xfirst, double xlast, double xinc, double *xvals)
{
  gridGenXvals(xsize, xfirst, xlast, xinc, xvals);
}

void
grid_gen_yvals(int gridtype, int ysize, double yfirst, double ylast, double yinc, double *yvals)
{
  gridGenYvals(gridtype, ysize, yfirst, ylast, yinc, yvals);
}

namespace cdo
{

int
inq_att_int(int cdiID, int varID, const std::string &attname)
{
  int attint = -1;
  cdiInqAttInt(cdiID, varID, attname.c_str(), 1, &attint);
  return attint;
}

std::string
inq_att_string(int cdiID, int varID, const std::string &attname)
{
  int attlen = cdiInqAttLen(cdiID, varID, attname.c_str());
  std::vector<char> atttxt(1, 0);
  if (attlen > 0)
    {
      atttxt.resize(attlen + 1);
      cdiInqAttTxt(cdiID, varID, attname.c_str(), attlen, atttxt.data());
      atttxt[attlen] = 0;
    }

  return std::string(atttxt.data());
}

std::string
inq_key_string(int cdiID, int varID, int key)
{
  char cstr[CDI_MAX_NAME] = { 0 };
  int length = CDI_MAX_NAME;
  cdiInqKeyString(cdiID, varID, key, cstr, &length);

  return std::string(cstr);
}

std::string
inq_var_name(int vlistID, int varID)
{
  char cstr[CDI_MAX_NAME] = { 0 };
  vlistInqVarName(vlistID, varID, cstr);
  return std::string(cstr);
}

std::string
inq_var_longname(int vlistID, int varID)
{
  char cstr[CDI_MAX_NAME] = { 0 };
  vlistInqVarLongname(vlistID, varID, cstr);
  return std::string(cstr);
}

std::string
inq_var_units(int vlistID, int varID)
{
  char cstr[CDI_MAX_NAME] = { 0 };
  vlistInqVarUnits(vlistID, varID, cstr);
  return std::string(cstr);
}

std::pair<int, HpOrder>
get_healpix_params(int gridID)
{
  auto projection = "healpix";
  auto nside = cdo::inq_att_int(gridID, CDI_GLOBAL, "healpix_nside");
  auto order = cdo::inq_att_string(gridID, CDI_GLOBAL, "healpix_order");
  if (nside == -1 || order.empty())
    {
      if (order.empty()) cdo_warning("%s mapping parameter %s missing!", projection, "healpix_order");
      if (nside == -1) cdo_warning("%s mapping parameter %s missing!", projection, "healpix_nside");
      cdo_abort("%s mapping parameter missing!", "healpix");
    }
  if (nside < 1) cdo_abort("%s mapping parameter %s < 1!", projection, "healpix_nside");
  auto hpOrder = hp_get_order(order);
  if (hpOrder == HpOrder::Undef) cdo_abort("%s mapping parameter healpix_order=%s unsupported!", projection, order);

  return std::make_pair(nside, hpOrder);
}

}
