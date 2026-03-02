/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "util_wildcards.h"
#include "util_string.h"

void
print_attributes(const VarList &varList, int vlistID, int varOrGlobal, int natts, char *argument)
{
  if (varOrGlobal != CDI_GLOBAL)
    {
      const auto &var = varList[varOrGlobal];
      auto stdname = cdo::inq_key_string(vlistID, varOrGlobal, CDI_KEY_STDNAME);

      double addoffset = 0.0, scalefactor = 1.0;
      auto haveAddoffset = (cdiInqKeyFloat(vlistID, varOrGlobal, CDI_KEY_ADDOFFSET, &addoffset) == CDI_NOERR);
      auto haveScalefactor = (cdiInqKeyFloat(vlistID, varOrGlobal, CDI_KEY_SCALEFACTOR, &scalefactor) == CDI_NOERR);

      if (argument)
        {
          if (stdname.size() && wildcardmatch(argument, "standard_name") == 0)
            fprintf(stdout, "   standard_name = \"%s\"\n", stdname.c_str());
          if (var.longname.size() && wildcardmatch(argument, "long_name") == 0)
            fprintf(stdout, "   long_name = \"%s\"\n", var.longname.c_str());
          if (var.units.size() && wildcardmatch(argument, "units") == 0) fprintf(stdout, "   units = \"%s\"\n", var.units.c_str());
          if (wildcardmatch(argument, "missing_value") == 0) fprintf(stdout, "   missing_value = %g\n", var.missval);
          if (haveAddoffset && wildcardmatch(argument, "add_offset") == 0) fprintf(stdout, "   add_offset = %g\n", addoffset);
          if (haveScalefactor && wildcardmatch(argument, "scale_factor") == 0)
            fprintf(stdout, "   scale_factor = %g\n", scalefactor);
        }
      else
        {
          if (stdname.size()) fprintf(stdout, "   standard_name = \"%s\"\n", stdname.c_str());
          if (var.longname.size()) fprintf(stdout, "   long_name = \"%s\"\n", var.longname.c_str());
          if (var.units.size()) fprintf(stdout, "   units = \"%s\"\n", var.units.c_str());
          fprintf(stdout, "   missing_value = %g\n", var.missval);
          if (haveAddoffset) fprintf(stdout, "   add_offset = %g\n", addoffset);
          if (haveScalefactor) fprintf(stdout, "   scale_factor = %g\n", scalefactor);
        }
    }

  for (int ia = 0; ia < natts; ia++)
    {
      char attname[CDI_MAX_NAME];
      int atttype, attlen;
      cdiInqAtt(vlistID, varOrGlobal, ia, attname, &atttype, &attlen);

      if (argument && wildcardmatch(argument, attname) != 0) continue;

      if (atttype == CDI_DATATYPE_TXT)
        {
          std::vector<char> atttxt(attlen + 1);
          cdiInqAttTxt(vlistID, varOrGlobal, attname, attlen, atttxt.data());
          atttxt[attlen] = 0;
          fprintf(stdout, "   %s = \"", attname);
          for (int i = 0; i < attlen; ++i)
            {
              if (atttxt[i] == '\n')
                {
                  printf("\\n");
                  if (atttxt[i + 1] != 0)
                    {
                      printf("\"\n");
                      printf("             \"");
                    }
                }
              else
                printf("%c", atttxt[i]);
            }
          printf("\"\n");
        }
      else if (atttype == CDI_DATATYPE_INT32)
        {
          std::vector<int> attint(attlen);
          cdiInqAttInt(vlistID, varOrGlobal, attname, attlen, attint.data());
          fprintf(stdout, "   %s = ", attname);
          for (int i = 0; i < attlen; ++i)
            {
              if (i) printf(", ");
              printf("%d", attint[i]);
            }
          printf("\n");
        }
      else if (atttype == CDI_DATATYPE_FLT32 || atttype == CDI_DATATYPE_FLT64)
        {
          char fltstr[128];
          std::vector<double> attflt(attlen);
          cdiInqAttFlt(vlistID, varOrGlobal, attname, attlen, attflt.data());
          fprintf(stdout, "   %s = ", attname);
          for (int i = 0; i < attlen; ++i)
            {
              if (i) printf(", ");
              if (atttype == CDI_DATATYPE_FLT32)
                printf("%sf", double_to_att_str(Options::CDO_flt_digits, fltstr, sizeof(fltstr), attflt[i]));
              else
                printf("%s", double_to_att_str(Options::CDO_dbl_digits, fltstr, sizeof(fltstr), attflt[i]));
            }
          printf("\n");
        }
      else { cdo_warning("Unsupported type %i name %s", atttype, attname); }
    }
}

void
check_varname_and_print(const VarList &varList, int vlistID, int nvars, char *checkvarname, char *attname)
{
  auto lfound = false;
  for (int varID = 0; varID < nvars; ++varID)
    {
      const auto &var = varList[varID];
      if (!checkvarname || (wildcardmatch(checkvarname, var.name) == 0))
        {
          lfound = true;
          fprintf(stdout, "%s:\n", var.name.c_str());
          int natts;
          cdiInqNatts(vlistID, varID, &natts);
          print_attributes(varList, vlistID, varID, natts, attname);
          if (!checkvarname) break;
        }
    }
  if (!lfound && checkvarname) cdo_abort("Could not find variable %s!", checkvarname);
}

void *
Showattribute(void *process)
{
  constexpr int delim = '@';
  cdo_initialize(process);

  auto SHOWATTRIBUTE = cdo_operator_add("showattribute", 0, 0, nullptr);
  auto SHOWATTSVAR = cdo_operator_add("showattsvar", 0, 0, nullptr);

  auto operatorID = cdo_operator_id();

  auto streamID = cdo_open_read(0);
  auto vlistID = cdo_stream_inq_vlist(streamID);

  auto nvars = vlistNvars(vlistID);

  VarList varList;
  varListInit(varList, vlistID);

  auto nargs = cdo_operator_argc();
  if (nargs == 0)
    {
      if (operatorID == SHOWATTSVAR)
        check_varname_and_print(varList, vlistID, nvars, nullptr, nullptr);
      else
        {
          for (int varID = 0; varID < nvars; ++varID)
            {
              const auto &var = varList[varID];
              fprintf(stdout, "%s:\n", var.name.c_str());

              int nattsvar;
              cdiInqNatts(vlistID, varID, &nattsvar);
              print_attributes(varList, vlistID, varID, nattsvar, nullptr);
            }

          int natts;
          cdiInqNatts(vlistID, CDI_GLOBAL, &natts);
          if (natts) fprintf(stdout, "Global:\n");
          print_attributes(varList, vlistID, CDI_GLOBAL, natts, nullptr);
        }
    }
  else
    {
      auto params = cdo_get_oper_argv();
      char buffer[CDI_MAX_NAME];
      for (int i = 0; i < nargs; ++i)
        {
          strcpy(buffer, params[i].c_str());
          char *result = strrchr(buffer, delim);
          char *input = buffer;
          if (result == nullptr)
            {
              if (operatorID == SHOWATTRIBUTE)
                {
                  int natts;
                  cdiInqNatts(vlistID, CDI_GLOBAL, &natts);
                  if (natts) fprintf(stdout, "Global:\n");
                  print_attributes(varList, vlistID, CDI_GLOBAL, natts, input);
                }
              else if (operatorID == SHOWATTSVAR)
                check_varname_and_print(varList, vlistID, nvars, input, nullptr);
            }
          else
            {
              if (operatorID == SHOWATTRIBUTE)
                {
                  input = result + 1;
                  if (*input == 0) input = nullptr;
                  *result = 0;
                  char *varname = buffer;
                  if (*varname == 0) cdo_abort("Variable name not specified!");
                  check_varname_and_print(varList, vlistID, nvars, varname, input);
                }
              else if (operatorID == SHOWATTSVAR)
                check_varname_and_print(varList, vlistID, nvars, input, nullptr);
            }
        }
    }

  cdo_stream_close(streamID);

  cdo_finish();

  return nullptr;
}
