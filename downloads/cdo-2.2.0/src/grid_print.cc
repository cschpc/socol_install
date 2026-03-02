/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>
#include "cdi_uuid.h"

#include "cdo_options.h"
#include "dmemory.h"
#include "process_int.h"
#include "util_string.h"

static void
printDblsPrefixAutoBrk(FILE *fp, int dig, const char *prefix, size_t n, const double vals[])
{
  int nbyte0 = strlen(prefix);
  fputs(prefix, fp);
  auto nbyte = nbyte0;
  for (size_t i = 0; i < n; ++i)
    {
      if (nbyte > 80)
        {
          fprintf(fp, "\n%*s", nbyte0, "");
          nbyte = nbyte0;
        }
      nbyte += fprintf(fp, "%.*g ", dig, vals[i]);
    }
  fputs("\n", fp);
}

static void
printIntsPrefixAutoBrk(FILE *fp, const char *prefix, size_t n, const int vals[])
{
  int nbyte0 = strlen(prefix);
  fputs(prefix, fp);
  auto nbyte = nbyte0;
  for (size_t i = 0; i < n; ++i)
    {
      if (nbyte > 80)
        {
          fprintf(fp, "\n%*s", nbyte0, "");
          nbyte = nbyte0;
        }
      nbyte += fprintf(fp, "%d ", vals[i]);
    }
  fputs("\n", fp);
}

static void
print_bounds(FILE *fp, int dig, const char *prefix, size_t n, size_t nvertex, const double bounds[])
{
  int nbyte0 = strlen(prefix);
  fputs(prefix, fp);
  for (size_t i = 0; i < n; ++i)
    {
      if (i > 0) fprintf(fp, "\n%*s", nbyte0, "");
      for (size_t iv = 0; iv < nvertex; iv++) fprintf(fp, "%.*g ", dig, bounds[i * nvertex + iv]);
    }
  fputs("\n", fp);
}

static void
printMask(FILE *fp, const char *prefix, size_t n, const int mask[])
{
  int nbyte0 = strlen(prefix);
  fputs(prefix, fp);
  auto nbyte = nbyte0;
  for (size_t i = 0; i < n; ++i)
    {
      if (nbyte > 80)
        {
          fprintf(fp, "\n%*s", (int) nbyte0, "");
          nbyte = nbyte0;
        }
      nbyte += (size_t) fprintf(fp, "%d ", mask[i]);
    }
  fputs("\n", fp);
}

static void
grid_print_attributes(FILE *fp, int gridID)
{
  int cdiID = gridID;
  int varID = CDI_GLOBAL;
  int atttype, attlen;
  char attname[CDI_MAX_NAME + 1];
  char fltstr[128];

  int natts;
  cdiInqNatts(cdiID, varID, &natts);

  for (int iatt = 0; iatt < natts; ++iatt)
    {
      cdiInqAtt(cdiID, varID, iatt, attname, &atttype, &attlen);
      if (attlen == 0) continue;

      if (cdo_cmpstr(attname, "grid_mapping_name")) continue;

      if (atttype == CDI_DATATYPE_TXT)
        {
          std::vector<char> atttxt(attlen + 1);
          cdiInqAttTxt(cdiID, varID, attname, attlen, atttxt.data());
          atttxt[attlen] = 0;
          if (strchr(atttxt.data(), '"'))
            fprintf(fp, "%s = '%s'\n", attname, atttxt.data());
          else
            fprintf(fp, "%s = \"%s\"\n", attname, atttxt.data());
        }
      else if (atttype == CDI_DATATYPE_INT8 || atttype == CDI_DATATYPE_UINT8 || atttype == CDI_DATATYPE_INT16
               || atttype == CDI_DATATYPE_UINT16 || atttype == CDI_DATATYPE_INT32 || atttype == CDI_DATATYPE_UINT32)
        {
          std::vector<int> attint(attlen);
          cdiInqAttInt(cdiID, varID, attname, attlen, attint.data());
          fprintf(fp, "%s =", attname);
          for (int i = 0; i < attlen; ++i) fprintf(fp, " %d", attint[i]);
          fprintf(fp, "\n");
        }
      else if (atttype == CDI_DATATYPE_FLT32 || atttype == CDI_DATATYPE_FLT64)
        {
          std::vector<double> attflt(attlen);
          cdiInqAttFlt(cdiID, varID, attname, attlen, attflt.data());
          fprintf(fp, "%s =", attname);
          if (atttype == CDI_DATATYPE_FLT32)
            for (int i = 0; i < attlen; ++i)
              fprintf(fp, " %sf", double_to_att_str(Options::CDO_flt_digits, fltstr, sizeof(fltstr), attflt[i]));
          else
            for (int i = 0; i < attlen; ++i)
              fprintf(fp, " %s", double_to_att_str(Options::CDO_dbl_digits, fltstr, sizeof(fltstr), attflt[i]));
          fprintf(fp, "\n");
        }
    }
}

static void
grid_print_kernel(int gridID, int opt, FILE *fp)
{
  size_t xdim, ydim;
  auto nxvals = gridInqXvals(gridID, nullptr);
  auto nyvals = gridInqYvals(gridID, nullptr);
  auto nxbounds = gridInqXbounds(gridID, nullptr);
  auto nybounds = gridInqYbounds(gridID, nullptr);

  auto gridsize = gridInqSize(gridID);
  auto xsize = gridInqXsize(gridID);
  auto ysize = gridInqYsize(gridID);
  auto type = gridInqType(gridID);
  auto nvertex = gridInqNvertex(gridID);
  auto xstrlen = gridInqXIsc(gridID);
  auto ystrlen = gridInqYIsc(gridID);
  int datatype = CDI_UNDEFID;
  cdiInqKeyInt(gridID, CDI_GLOBAL, CDI_KEY_DATATYPE, &datatype);

  char **xcvals = nullptr;
  if (xstrlen)
    {
      xcvals = (char **) Malloc(xsize * sizeof(char *));
      for (size_t i = 0; i < xsize; ++i) xcvals[i] = (char *) Malloc((xstrlen + 1) * sizeof(char));
      gridInqXCvals(gridID, xcvals);
      for (size_t i = 0; i < xsize; ++i) xcvals[i][xstrlen] = 0;
      for (size_t i = 0; i < xsize; ++i)
        for (int k = xstrlen - 1; k; k--)
          {
            if (xcvals[i][k] == ' ')
              xcvals[i][k] = 0;
            else
              break;
          }
    }

  char **ycvals = nullptr;
  if (ystrlen)
    {
      ycvals = (char **) Malloc(ysize * sizeof(char *));
      for (size_t i = 0; i < ysize; ++i) ycvals[i] = (char *) Malloc((ystrlen + 1) * sizeof(char));
      gridInqYCvals(gridID, ycvals);
      for (size_t i = 0; i < ysize; ++i) ycvals[i][ystrlen] = 0;
      for (size_t i = 0; i < ysize; ++i)
        for (int k = ystrlen - 1; k; k--)
          {
            if (ycvals[i][k] == ' ')
              ycvals[i][k] = 0;
            else
              break;
          }
    }

  int dig = (datatype == CDI_DATATYPE_FLT64) ? Options::CDO_dbl_digits : Options::CDO_flt_digits;

  fprintf(fp, "gridtype  = %s\n", gridNamePtr(type));
  fprintf(fp, "gridsize  = %zu\n", gridsize);
  if (datatype == CDI_DATATYPE_FLT32) fprintf(fp, "datatype  = float\n");

  if (type != GRID_GME)
    {
      if (type != GRID_UNSTRUCTURED && type != GRID_SPECTRAL && type != GRID_FOURIER)
        {
          if (xsize > 0) fprintf(fp, "xsize     = %zu\n", xsize);
          if (ysize > 0) fprintf(fp, "ysize     = %zu\n", ysize);
        }

      if (nxvals > 0 || xcvals)
        {
          auto name = cdo::inq_key_string(gridID, CDI_XAXIS, CDI_KEY_NAME);
          auto dimname = cdo::inq_key_string(gridID, CDI_XAXIS, CDI_KEY_DIMNAME);
          auto longname = cdo::inq_key_string(gridID, CDI_XAXIS, CDI_KEY_LONGNAME);
          auto units = cdo::inq_key_string(gridID, CDI_XAXIS, CDI_KEY_UNITS);
          if (name.size()) fprintf(fp, "xname     = %s\n", name.c_str());
          if (dimname.size() && dimname != name) fprintf(fp, "xdimname  = %s\n", dimname.c_str());
          if (longname.size()) fprintf(fp, "xlongname = \"%s\"\n", longname.c_str());
          if (units.size()) fprintf(fp, "xunits    = \"%s\"\n", units.c_str());
        }

      if (nyvals > 0 || ycvals)
        {
          auto name = cdo::inq_key_string(gridID, CDI_YAXIS, CDI_KEY_NAME);
          auto dimname = cdo::inq_key_string(gridID, CDI_YAXIS, CDI_KEY_DIMNAME);
          auto longname = cdo::inq_key_string(gridID, CDI_YAXIS, CDI_KEY_LONGNAME);
          auto units = cdo::inq_key_string(gridID, CDI_YAXIS, CDI_KEY_UNITS);
          if (name.size()) fprintf(fp, "yname     = %s\n", name.c_str());
          if (dimname.size() && dimname != name) fprintf(fp, "ydimname  = %s\n", dimname.c_str());
          if (longname.size()) fprintf(fp, "ylongname = \"%s\"\n", longname.c_str());
          if (units.size()) fprintf(fp, "yunits    = \"%s\"\n", units.c_str());
        }

      if (type == GRID_UNSTRUCTURED || type == GRID_CURVILINEAR)
        {
          auto vdimName = cdo::inq_key_string(gridID, CDI_GLOBAL, CDI_KEY_VDIMNAME);
          if (vdimName.size()) fprintf(fp, "vdimname  = %s\n", vdimName.c_str());
        }
      if (type == GRID_UNSTRUCTURED && nvertex > 0) fprintf(fp, "nvertex   = %d\n", nvertex);
    }

  switch (type)
    {
    case GRID_LONLAT:
    case GRID_GAUSSIAN:
    case GRID_GAUSSIAN_REDUCED:
    case GRID_GENERIC:
    case GRID_PROJECTION:
    case GRID_CURVILINEAR:
    case GRID_UNSTRUCTURED:
    case GRID_TRAJECTORY:
    case GRID_CHARXY:
      {
        if (type == GRID_GAUSSIAN || type == GRID_GAUSSIAN_REDUCED) fprintf(fp, "numLPE    = %d\n", gridInqNP(gridID));

        if (type == GRID_CURVILINEAR || type == GRID_UNSTRUCTURED)
          {
            xdim = gridsize;
            ydim = gridsize;
          }
        else if (type == GRID_GAUSSIAN_REDUCED)
          {
            xdim = 2;
            ydim = ysize;
          }
        else
          {
            xdim = xsize;
            ydim = ysize;
          }

        if (type == GRID_UNSTRUCTURED)
          {
            int number = 0;
            cdiInqKeyInt(gridID, CDI_GLOBAL, CDI_KEY_NUMBEROFGRIDUSED, &number);
            if (number > 0)
              {
                fprintf(fp, "number    = %d\n", number);
                int position = 0;
                cdiInqKeyInt(gridID, CDI_GLOBAL, CDI_KEY_NUMBEROFGRIDINREFERENCE, &position);
                if (position >= 0) fprintf(fp, "position  = %d\n", position);
              }

            int length = 0;
            if (CDI_NOERR == cdiInqKeyLen(gridID, CDI_GLOBAL, CDI_KEY_REFERENCEURI, &length))
              {
                char referenceLink[8192];
                cdiInqKeyString(gridID, CDI_GLOBAL, CDI_KEY_REFERENCEURI, referenceLink, &length);
                fprintf(fp, "uri       = %s\n", referenceLink);
              }
          }

        if (nxvals > 0)
          {
            double xfirst = 0.0, xinc = 0.0;

            if (type == GRID_LONLAT || type == GRID_GAUSSIAN || type == GRID_PROJECTION || type == GRID_GENERIC)
              {
                xfirst = gridInqXval(gridID, 0);
                xinc = gridInqXinc(gridID);
              }

            if (is_not_equal(xinc, 0.0) && opt)
              {
                fprintf(fp, "xfirst    = %.*g\n", dig, xfirst);
                fprintf(fp, "xinc      = %.*g\n", dig, xinc);
              }
            else
              {
                std::vector<double> xvals(nxvals);
                gridInqXvals(gridID, xvals.data());
                printDblsPrefixAutoBrk(fp, dig, "xvals     = ", nxvals, xvals.data());
              }
          }

        if (nxbounds)
          {
            std::vector<double> xbounds(nxbounds);
            gridInqXbounds(gridID, xbounds.data());
            print_bounds(fp, dig, "xbounds   = ", xdim, nvertex, xbounds.data());
          }

        if (xcvals)
          {
            fprintf(fp, "xcvals    = \"%.*s\"", xstrlen, xcvals[0]);
            for (size_t i = 1; i < xsize; ++i) fprintf(fp, ", \"%.*s\"", xstrlen, xcvals[i]);
            fprintf(fp, "\n");

            for (size_t i = 0; i < xsize; ++i) Free(xcvals[i]);
            Free(xcvals);
          }

        if (nyvals > 0)
          {
            double yfirst = 0.0, yinc = 0.0;

            if (type == GRID_LONLAT || type == GRID_GENERIC || type == GRID_PROJECTION || type == GRID_GENERIC)
              {
                yfirst = gridInqYval(gridID, 0);
                yinc = gridInqYinc(gridID);
              }

            if (is_not_equal(yinc, 0.0) && opt)
              {
                fprintf(fp, "yfirst    = %.*g\n", dig, yfirst);
                fprintf(fp, "yinc      = %.*g\n", dig, yinc);
              }
            else
              {
                std::vector<double> yvals(nyvals);
                gridInqYvals(gridID, yvals.data());
                printDblsPrefixAutoBrk(fp, dig, "yvals     = ", nyvals, yvals.data());
              }
          }

        if (nybounds)
          {
            std::vector<double> ybounds(nybounds);
            gridInqYbounds(gridID, ybounds.data());
            print_bounds(fp, dig, "ybounds   = ", ydim, nvertex, ybounds.data());
          }

        if (ycvals)
          {
            fprintf(fp, "ycvals    = \"%.*s\"", ystrlen, ycvals[0]);
            for (size_t i = 1; i < ysize; ++i) fprintf(fp, ", \"%.*s\"", ystrlen, ycvals[i]);
            fprintf(fp, "\n");

            for (size_t i = 0; i < ysize; ++i) Free(ycvals[i]);
            Free(ycvals);
          }

        if (gridHasArea(gridID))
          {
            std::vector<double> area(gridsize);
            gridInqArea(gridID, area.data());
            printDblsPrefixAutoBrk(fp, dig, "area      = ", gridsize, area.data());
          }

        if (type == GRID_GAUSSIAN_REDUCED)
          {
            std::vector<int> reducedPoints(ysize);
            gridInqReducedPoints(gridID, reducedPoints.data());
            printIntsPrefixAutoBrk(fp, "reducedPoints = ", (size_t) (ysize > 0 ? ysize : 0), reducedPoints.data());
          }

#ifdef HIRLAM_EXTENSIONS
        {
          int scanningMode = 0;
          cdiInqKeyInt(gridID, CDI_GLOBAL, CDI_KEY_SCANNINGMODE, &scanningMode);
          fprintf(fp, "scanningMode = %d\n", scanningMode);
        }
#endif

        if (type == GRID_PROJECTION)
          {
            auto gridMapping = cdo::inq_key_string(gridID, CDI_GLOBAL, CDI_KEY_GRIDMAP_VARNAME);
            auto gridMappingName = cdo::inq_key_string(gridID, CDI_GLOBAL, CDI_KEY_GRIDMAP_NAME);
            if (gridMapping.size()) fprintf(fp, "grid_mapping = %s\n", gridMapping.c_str());
            if (gridMappingName.size()) fprintf(fp, "grid_mapping_name = %s\n", gridMappingName.c_str());
            grid_print_attributes(fp, gridID);
          }

        break;
      }
    case GRID_SPECTRAL:
      {
        fprintf(fp, "truncation = %d\n", gridInqTrunc(gridID));
        fprintf(fp, "complexpacking = %d\n", gridInqComplexPacking(gridID));
        break;
      }
    case GRID_FOURIER:
      {
        fprintf(fp, "truncation = %d\n", gridInqTrunc(gridID));
        break;
      }
    case GRID_GME:
      {
        int nd, ni, ni2, ni3;
        gridInqParamGME(gridID, &nd, &ni, &ni2, &ni3);
        fprintf(fp, "ni        = %d\n", ni);
        break;
      }
    default:
      {
        fprintf(stderr, "Unsupported grid type: %s\n", gridNamePtr(type));
        break;
      }
    }

  unsigned char uuid[CDI_UUID_SIZE] = { 0 };
  int length = CDI_UUID_SIZE;
  cdiInqKeyBytes(gridID, CDI_GLOBAL, CDI_KEY_UUID, uuid, &length);
  if (!cdiUUIDIsNull(uuid))
    {
      char uuidStr[uuidNumHexChars + 1] = { 0 };
      if (cdiUUID2Str(uuid, uuidStr) == uuidNumHexChars) fprintf(fp, "uuid      = %s\n", uuidStr);
    }

  if (gridInqMask(gridID, nullptr))
    {
      std::vector<int> mask(gridsize);
      gridInqMask(gridID, mask.data());
      printMask(fp, "mask      = ", (size_t) (gridsize > 0 ? gridsize : 0), mask.data());
    }

  auto projID = gridInqProj(gridID);
  if (projID != CDI_UNDEFID && gridInqType(projID) == GRID_PROJECTION) grid_print_kernel(projID, opt, fp);
}

void
cdo_print_griddes(int gridID, int opt)
{
  grid_print_kernel(gridID, opt, stdout);
}
