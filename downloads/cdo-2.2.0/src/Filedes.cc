/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Filedes    codetab         Parameter code table
      Filedes    griddes         Grid description
      Filedes    vct             Vertical coordinate table
*/

#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include <mpim_grid.h>
#include "util_string.h"

void cdoPrintZaxis(int zaxisID);

void
cdoPrintAttributes(FILE *fp, int cdiID, int varID, int nblanks)
{
  int natts;
  cdiInqNatts(cdiID, varID, &natts);

  for (int ia = 0; ia < natts; ++ia)
    {
      char attname[CDI_MAX_NAME];
      int atttype, attlen;
      cdiInqAtt(cdiID, varID, ia, attname, &atttype, &attlen);

      if (atttype == CDI_DATATYPE_INT8 || atttype == CDI_DATATYPE_UINT8 || atttype == CDI_DATATYPE_INT16
          || atttype == CDI_DATATYPE_UINT16 || atttype == CDI_DATATYPE_INT32 || atttype == CDI_DATATYPE_UINT32)
        {
          std::vector<int> attint(attlen);
          cdiInqAttInt(cdiID, varID, attname, attlen, attint.data());
          fprintf(fp, "%*s", nblanks, "");
          fprintf(fp, "%s = ", attname);
          for (int i = 0; i < attlen; ++i)
            {
              if (i) fprintf(fp, ", ");
              fprintf(fp, "%d", attint[i]);
            }
          fprintf(fp, "\n");
        }
      else if (atttype == CDI_DATATYPE_FLT32 || atttype == CDI_DATATYPE_FLT64)
        {
          char fltstr[128];
          std::vector<double> attflt(attlen);
          cdiInqAttFlt(cdiID, varID, attname, attlen, attflt.data());
          fprintf(fp, "%*s", nblanks, "");
          fprintf(fp, "%s = ", attname);
          for (int i = 0; i < attlen; ++i)
            {
              if (i) fprintf(fp, ", ");
              if (atttype == CDI_DATATYPE_FLT32)
                fprintf(fp, "%sf", double_to_att_str(Options::CDO_flt_digits, fltstr, sizeof(fltstr), attflt[i]));
              else
                fprintf(fp, "%s", double_to_att_str(Options::CDO_dbl_digits, fltstr, sizeof(fltstr), attflt[i]));
            }
          fprintf(fp, "\n");
        }
      else if (atttype == CDI_DATATYPE_TXT)
        {
          std::vector<char> atttxt(attlen + 1);
          cdiInqAttTxt(cdiID, varID, attname, attlen, atttxt.data());
          atttxt[attlen] = 0;
          fprintf(fp, "%*s", nblanks, "");
          fprintf(fp, "%s = \"%s\"\n", attname, atttxt.data());
        }
    }
}

static void
printSource(FILE *fp, int vlistID, int varID)
{
  // institute info
  auto instptr = institutInqLongnamePtr(vlistInqVarInstitut(vlistID, varID));
  if (instptr) fprintf(fp, "  institution = \"%s\"\n", instptr);

  // source info
  auto modelptr = modelInqNamePtr(vlistInqVarModel(vlistID, varID));
  if (modelptr) fprintf(fp, "  source = \"%s\"\n", modelptr);
}

static void
printVCT(int vlistID, bool lvct)
{
  auto nzaxis = vlistNzaxis(vlistID);
  for (int index = 0; index < nzaxis; ++index)
    {
      auto zaxisID = vlistZaxis(vlistID, index);
      auto type = zaxisInqType(zaxisID);
      if (type == ZAXIS_HYBRID || type == ZAXIS_HYBRID_HALF)
        {
          auto vctsize = zaxisInqVctSize(zaxisID);
          auto vct = zaxisInqVctPtr(zaxisID);

          if (vctsize % 2 == 0)
            {
              if (lvct)
                {
                  fprintf(stdout, "#   k         vct_a(k) [Pa]             vct_b(k) []\n");
                  for (int i = 0; i < vctsize / 2; ++i) fprintf(stdout, "%5d %25.17f %25.17f\n", i, vct[i], vct[vctsize / 2 + i]);
                }
              else
                {
                  fprintf(stdout, "vctsize   = %d\n", vctsize);
                  int nbyte0 = fprintf(stdout, "vct       = ");
                  int nbyte = nbyte0;
                  for (int i = 0; i < vctsize; ++i)
                    {
                      if (nbyte > 70 || i == vctsize / 2)
                        {
                          fprintf(stdout, "\n%*s", nbyte0, "");
                          nbyte = nbyte0;
                        }
                      nbyte += fprintf(stdout, "%.9g ", vct[i]);
                    }
                  fprintf(stdout, "\n");
                }
            }
          else
            for (int i = 0; i < vctsize; ++i) fprintf(stdout, "%5d %25.17f\n", i, vct[i]);

          break;
        }
    }
}

static void
printCodeTable(const VarList &varList)
{
  int nvars = varList.size();
  for (int varID = 0; varID < nvars; ++varID)
    {
      fprintf(stdout, "%4d  %-12s", varList[varID].code, varList[varID].name.c_str());
      if (varList[varID].longname.size())
        {
          fprintf(stdout, "  %s", varList[varID].longname.c_str());
          if (varList[varID].units.size()) fprintf(stdout, " [%s]", varList[varID].units.c_str());
        }
      fprintf(stdout, "\n");
    }
}

static void
partab(FILE *fp, int vlistID, const VarList &varList, int option)
{
  int varID, datatype = -1;
  char paramstr[32];

  int nvars = varList.size();
  auto linebreak = (option != 4);

  if (option == 2)
    {
      int natts;
      cdiInqNatts(vlistID, CDI_GLOBAL, &natts);
      if (natts > 0)
        {
          fprintf(fp, "&parameter\n");
          fprintf(fp, "  name = _GLOBAL_\n");
          printSource(fp, vlistID, 0);
          cdoPrintAttributes(fp, vlistID, CDI_GLOBAL, 2);
          fprintf(fp, "/\n");
        }
    }

  if (nvars > 1)
    {
      datatype = varList[0].datatype;
      for (varID = 1; varID < nvars; ++varID)
        {
          if (datatype != varList[varID].datatype)
            {
              datatype = -1;
              break;
            }
        }

      if (datatype != -1)
        {
          fprintf(fp, "&parameter");
          if (linebreak) fprintf(fp, "\n");
          fprintf(fp, "  name = _default_");
          if (linebreak) fprintf(fp, "\n");
          auto datatypestr = cdo::datatype_to_cstr(datatype);
          if (*datatypestr)
            {
              fprintf(fp, "  datatype = %s", datatypestr);
              if (linebreak) fprintf(fp, "\n");
            }
          fprintf(fp, " /\n");
        }
    }

  for (varID = 0; varID < nvars; ++varID)
    {
      const auto &var = varList[varID];

      fprintf(fp, "&parameter");
      if (linebreak) fprintf(fp, "\n");

      auto stdname = cdo::inq_key_string(vlistID, varID, CDI_KEY_STDNAME);

      fprintf(fp, "  name = %s", varList[varID].name.c_str());
      if (linebreak) fprintf(fp, "\n");

      if (var.param >= 0)
        {
          cdiParamToString(var.param, paramstr, sizeof(paramstr));
          fprintf(fp, "  param = %s", paramstr);
          if (linebreak) fprintf(fp, "\n");
        }
      if (stdname.size())
        {
          fprintf(fp, "  standard_name = %s", stdname.c_str());
          if (linebreak) fprintf(fp, "\n");
        }
      if (var.longname.size())
        {
          fprintf(fp, "  long_name = \"%s\"", var.longname.c_str());
          if (linebreak) fprintf(fp, "\n");
        }
      if (var.units.size())
        {
          fprintf(fp, "  units = \"%s\"", var.units.c_str());
          if (linebreak) fprintf(fp, "\n");
        }

      if (datatype == -1)
        {
          auto datatypestr = cdo::datatype_to_cstr(varList[varID].datatype);
          if (*datatypestr)
            {
              fprintf(fp, "  datatype = %s", datatypestr);
              if (linebreak) fprintf(fp, "\n");
            }
        }

      int uvRelativeToGrid = 0;
      if (cdiInqKeyInt(vlistID, varID, CDI_KEY_UVRELATIVETOGRID, &uvRelativeToGrid) == CDI_NOERR)
        {
          fprintf(fp, "  uvRelativeToGrid = %d", uvRelativeToGrid);
          if (linebreak) fprintf(fp, "\n");
        }

      int chunkType = -1;
      cdiInqKeyInt(vlistID, varID, CDI_KEY_CHUNKTYPE, &chunkType);
      const char *chunkName = (chunkType == CDI_CHUNK_AUTO)
                                  ? "auto"
                                  : ((chunkType == CDI_CHUNK_GRID) ? "grid" : ((chunkType == CDI_CHUNK_LINES) ? "lines" : nullptr));

      if (chunkName)
        {
          fprintf(fp, "  chunkType = %s", chunkName);
          if (linebreak) fprintf(fp, "\n");
        }

      if (option == 2)
        {
          fprintf(fp, "  missing_value = %g\n", var.missval);
          cdoPrintAttributes(fp, vlistID, varID, 2);
        }

      if (!linebreak) fprintf(fp, "  ");
      fprintf(fp, "/\n");
    }
}

static void
filedes(CdoStreamID streamID)
{
  printf("\n");
  auto filetype = cdo_inq_filetype(streamID);

  auto filetypestr = cdo::filetype_to_cstr(filetype);
  if (filetypestr == nullptr || *filetypestr == 0)
    printf("  unsupported filetype %d\n", filetype);
  else
    printf("  %s data\n", filetypestr);

  switch (filetype)
    {
    case CDI_FILETYPE_SRV:
    case CDI_FILETYPE_EXT:
    case CDI_FILETYPE_IEG:
      {
        auto byteorder = cdo_inq_byteorder(streamID);
        switch (byteorder)
          {
          case CDI_BIGENDIAN: printf("  byteorder is BIGENDIAN\n"); break;
          case CDI_LITTLEENDIAN: printf("  byteorder is LITTLEENDIAN\n"); break;
          default: printf("  byteorder %d undefined\n", byteorder); break;
          }
      }
    }

  printf("\n");
}

void *
Filedes(void *process)
{
  cdo_initialize(process);

  // clang-format off
  auto GRIDDES  = cdo_operator_add("griddes",   0, 0, nullptr);
  auto GRIDDES2 = cdo_operator_add("griddes2",  0, 0, nullptr);
  auto ZAXISDES = cdo_operator_add("zaxisdes",  0, 0, nullptr);
  auto VCT      = cdo_operator_add("vct",       0, 0, nullptr);
  auto VCT2     = cdo_operator_add("vct2",      0, 0, nullptr);
  auto CODETAB  = cdo_operator_add("codetab",   0, 0, nullptr);
  auto FILEDES  = cdo_operator_add("filedes",   0, 0, nullptr);
  auto VLIST    = cdo_operator_add("vlist",     0, 0, nullptr);
  auto SPARTAB  = cdo_operator_add("spartab",   0, 0, nullptr);
  auto PARTAB   = cdo_operator_add("partab",    0, 0, nullptr);
  auto PARTAB2  = cdo_operator_add("partab2",   0, 0, nullptr);
  // clang-format on

  auto operatorID = cdo_operator_id();

  operator_check_argc(0);

  auto streamID = cdo_open_read(0);
  auto vlistID = cdo_stream_inq_vlist(streamID);

  VarList varList;
  varListInit(varList, vlistID);

  if (operatorID == GRIDDES || operatorID == GRIDDES2)
    {
      auto opt = (operatorID == GRIDDES) ? 1 : 0;
      auto ngrids = vlistNgrids(vlistID);
      for (int index = 0; index < ngrids; ++index)
        {
          printf("#\n# gridID %d\n#\n", index + 1);
          cdo_print_griddes(vlistGrid(vlistID, index), opt);
          auto nsubtypes = vlistNsubtypes(vlistID);
          for (int i = 0; i < nsubtypes; ++i) subtypePrint(vlistSubtype(vlistID, i));
        }
    }
  else if (operatorID == ZAXISDES)
    {
      auto nzaxis = vlistNzaxis(vlistID);
      for (int index = 0; index < nzaxis; ++index)
        {
          printf("#\n# zaxisID %d\n#\n", index + 1);
          cdoPrintZaxis(vlistZaxis(vlistID, index));
        }
    }
  else if (operatorID == VCT || operatorID == VCT2) { printVCT(vlistID, operatorID == VCT); }
  else if (operatorID == VLIST) { vlistPrint(vlistID); }
  else if (operatorID == CODETAB) { printCodeTable(varList); }
  else if (operatorID == PARTAB || operatorID == SPARTAB || operatorID == PARTAB2)
    {
      auto option = (operatorID == SPARTAB) ? 4 : ((operatorID == PARTAB2) ? 2 : 1);
      partab(stdout, vlistID, varList, option);
    }
  else if (operatorID == FILEDES) { filedes(streamID); }

  cdo_stream_close(streamID);

  if (!cdo::stdoutIsTerminal) Options::silentMode = true;

  cdo_finish();

  return nullptr;
}
