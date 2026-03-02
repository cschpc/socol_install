/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/
#include <cdi.h>

#include "cdi_uuid.h"
#include "cdo_options.h"
#include "process_int.h"

void cdoPrintAttributes(FILE *fp, int cdiID, int varID, int nblanks);

static void
printDblsPrefixAutoBrk(FILE *fp, int dig, const char *prefix, size_t n, const double vals[], size_t extbreak)
{
  int nbyte0 = strlen(prefix);
  fputs(prefix, fp);
  int nbyte = nbyte0;
  for (size_t i = 0; i < n; ++i)
    {
      if (nbyte > 80 || (i && i == extbreak))
        {
          fprintf(fp, "\n%*s", nbyte0, "");
          nbyte = nbyte0;
        }
      nbyte += fprintf(fp, "%.*g ", dig, vals[i]);
    }
  fputs("\n", fp);
}

static void
zaxisPrintKernel(int zaxisID, FILE *fp)
{
  auto type = zaxisInqType(zaxisID);
  auto nlevels = zaxisInqSize(zaxisID);
  int datatype = CDI_UNDEFID;
  cdiInqKeyInt(zaxisID, CDI_GLOBAL, CDI_KEY_DATATYPE, &datatype);
  size_t nvals = (size_t) zaxisInqLevels(zaxisID, nullptr);

  int dig = (datatype == CDI_DATATYPE_FLT64) ? Options::CDO_dbl_digits : Options::CDO_flt_digits;

  fprintf(fp, "zaxistype = %s\n", zaxisNamePtr(type));
  fprintf(fp, "size      = %d\n", nlevels);
  // clang-format off
  if      (datatype == CDI_DATATYPE_FLT32) fprintf(fp, "datatype  = float\n");
  else if (datatype == CDI_DATATYPE_INT32) fprintf(fp, "datatype  = int\n");
  else if (datatype == CDI_DATATYPE_INT16) fprintf(fp, "datatype  = short\n");
  // clang-format on

  if (nlevels == 1 && zaxisInqScalar(zaxisID)) fprintf(fp, "scalar    = true\n");

  auto zname = cdo::inq_key_string(zaxisID, CDI_GLOBAL, CDI_KEY_NAME);
  auto zlongname = cdo::inq_key_string(zaxisID, CDI_GLOBAL, CDI_KEY_LONGNAME);
  auto zunits = cdo::inq_key_string(zaxisID, CDI_GLOBAL, CDI_KEY_UNITS);
  if (zname.size()) fprintf(fp, "name      = %s\n", zname.c_str());
  if (zlongname.size()) fprintf(fp, "longname  = \"%s\"\n", zlongname.c_str());
  if (zunits.size()) fprintf(fp, "units     = \"%s\"\n", zunits.c_str());

  std::vector<double> vals;
  if (nvals) vals.resize(nvals);

  if (nvals)
    {
      zaxisInqLevels(zaxisID, vals.data());
      printDblsPrefixAutoBrk(fp, dig, "levels    = ", nvals, vals.data(), 0);
    }
  else if (type == ZAXIS_CHAR)
    {
      int clen = zaxisInqCLen(zaxisID);
      char **cvals = nullptr;
      zaxisInqCVals(zaxisID, &cvals);
      fprintf(fp, "levels    = \n");
      for (int i = 0; i < nlevels; ++i)
        {
          fprintf(fp, "     [%2d] = %.*s\n", i, clen, cvals[i]);
          free(cvals[i]);
        }
      if (cvals) free(cvals);
    }

  if (zaxisInqLbounds(zaxisID, nullptr) && zaxisInqUbounds(zaxisID, nullptr))
    {
      zaxisInqLbounds(zaxisID, vals.data());
      printDblsPrefixAutoBrk(fp, dig, "lbounds   = ", nvals, vals.data(), 0);

      zaxisInqUbounds(zaxisID, vals.data());
      printDblsPrefixAutoBrk(fp, dig, "ubounds   = ", nvals, vals.data(), 0);
    }

  if (type == ZAXIS_HYBRID || type == ZAXIS_HYBRID_HALF)
    {
      int vctsize = zaxisInqVctSize(zaxisID);
      if (vctsize)
        {
          fprintf(fp, "vctsize   = %d\n", vctsize);
          std::vector<double> vct(vctsize);
          zaxisInqVct(zaxisID, vct.data());
          printDblsPrefixAutoBrk(fp, dig, "vct       = ", vctsize, vct.data(), vctsize / 2);
        }
    }

  if (type == ZAXIS_REFERENCE)
    {
      unsigned char uuid[CDI_UUID_SIZE] = { 0 };
      int length = CDI_UUID_SIZE;
      cdiInqKeyBytes(zaxisID, CDI_GLOBAL, CDI_KEY_UUID, uuid, &length);
      if (!cdiUUIDIsNull(uuid))
        {
          char uuidStr[uuidNumHexChars + 1] = { 0 };
          if (cdiUUID2Str(uuid, uuidStr) == uuidNumHexChars) fprintf(fp, "uuid      = %s\n", uuidStr);
        }
    }

  cdoPrintAttributes(fp, zaxisID, CDI_GLOBAL, 0);
}

void
cdoPrintZaxis(int zaxisID)
{
  zaxisPrintKernel(zaxisID, stdout);
}
