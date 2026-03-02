/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Fabian Wachsmann

*/

#include <cdi.h>
#include "julian_date.h"
#include <signal.h>

#include "process_int.h"
#include "param_conversion.h"
#include "util_files.h"
#include "cdi_lockedIO.h"
#include "cdo_options.h"
#include "varray.h"
#include "dmemory.h"

#ifdef HAVE_LIBCMOR
#include <cassert>
#include <unistd.h>

extern "C" {
#include "cmor.h"
}

#include "netcdf.h"
#include "util_string.h"
#include "pmlist.h"
#include "datetime.h"
#include "mapping.h"
#include "merge_axis.h"
#include "listbuffer.h"

#define CMOR_UNDEFID (CMOR_MAX_AXES + 1)

static void
get_stringcode(int vlistID, int varID, char *varcodestring)
{
  int varcode = vlistInqVarCode(vlistID, varID);
  sprintf(varcodestring, "%03d", varcode);
}

static int
get_ifilevalue_code(int vlistID, const char *value, int nvars)
{
  int code = atol(value);
  if ( code > 0 && code < 1000 )
    {
      char newcode[4];
      sprintf(newcode, "%03d", code);
      for (int varID = 0; varID < nvars; ++varID)
        {
          char codestring[CDI_MAX_NAME];
          get_stringcode(vlistID, varID, codestring);
          if (strcmp(codestring, newcode) == 0)
            return varID;
        }
      return CDI_UNDEFID;
    }
  else
    {
      cdo_print("'%s' could not be transformed into the code format (three digit integer). Code wont be used.", value);
      return CDI_UNDEFID;
    }
}

static int
get_ifilevalue_name(int vlistID, const char *value, int nvars)
{
  char ifilevalue[CDI_MAX_NAME];
  for (int varID = 0; varID < nvars; ++varID)
    {
      vlistInqVarName(vlistID, varID, ifilevalue);
      if (strcmp(ifilevalue, value) == 0)
        return varID;
    }
  return CDI_UNDEFID;
}

static int
getVarIDToMap(int vlistID, int nvars, const char *key, const char *value)
{
  if ( strcmp(key, "code") == 0 )
    return get_ifilevalue_code(vlistID, value, nvars);
  else
    return get_ifilevalue_name(vlistID, value, nvars);
}

#if (CMOR_VERSION_MAJOR == 3)
static void
removeDataset()
{
  char cwd[1024];
  getcwd(cwd, sizeof(cwd));
  cwd[strlen(cwd)] = '\0';
  int procID = getpid();
  char *dataset_path = (char *) Malloc((strlen(cwd) +  
                                        strlen("/dataset.json") + 
                                        floor (log10 (abs (procID))) + 1 +
                                       1) * sizeof(char));
  sprintf(dataset_path, "%s/dataset%d.json", cwd, procID);
  remove(dataset_path);
  Free(dataset_path);
}
#endif

void
sigfunc(int sig)
{
  if (sig == SIGTERM)
    {
#if (CMOR_VERSION_MAJOR == 3)
      removeDataset();
#endif
      cdo_abort("ERROR (infile: '%s')! Program terminated by CMOR. A temporary ofile can outlive which needs to be deleted manually.", cdo_get_stream_name(0));
    }
}

static char *
kv_get_a_val(KVList *kvl, const char *key, const char *replacer)
{
    return kvl->get_first_value(key,replacer);
}

KVList 
maptab_search_miptab(PMList pmlist, const char *cmorname, const char *miptab, const char *key)
{
  KVList listlatest;
  if (pmlist.size() && cmorname && miptab)
    {
      for (auto &node : pmlist)
        {
          const KeyValues *kvcn = node.search(key);
          if (kvcn && kvcn->nvalues > 0 && *(kvcn->values[0].c_str()) == *cmorname && kvcn->values[0] == cmorname)
            {
              listlatest = node;
              const KeyValues *kvmt = node.search("project_mip_table");
              if ((kvmt && kvmt->nvalues > 0 && *(kvmt->values[0].c_str()) == *miptab && kvmt->values[0] == miptab) || !kvmt)
                break;
            }
        }
    }

  return listlatest;
}

static void
handleError(char *filename, int errnum, const char *argument)
{
  if ( !filename )
    cdo_abort("ERROR (infile: '%s')! In parsing the command line:\n          More than 150 values for a key are not supported.", cdo_get_stream_name(0));
  switch (errnum)
    {
    case (1):
      cdo_abort("ERROR (infile: '%s')! In parsing file '%s':\n          "
               "Unexpected blank in line:\n          '%s'\n          Check syntax.",
               cdo_get_stream_name(0), filename, argument); break;
    case (2):
      cdo_abort("ERROR (infile: '%s')! In parsing file '%s':\n          "
               "Unexpected separator sign ',' in a key of line:\n          '%s'\n          Check syntax.",
               cdo_get_stream_name(0), filename, argument); break;
    case (3):
      cdo_abort("ERROR (infile: '%s')! In parsing file '%s':\n          "
               "More than 150 values for a key are not supported.", cdo_get_stream_name(0), filename); break;
    case (4):
      cdo_abort("ERROR (infile: '%s')! In parsing file '%s':\n          "
               "No values found for a keyword in line:\n          '%s'.\n          Check syntax.",
               cdo_get_stream_name(0), filename, argument); break;
    case (5):
      cdo_abort("ERROR (infile: '%s')! In parsing file '%s':\n          "
               "A value for a keyword begins with ',' in line:\n          '%s'.\n          Check syntax.",
               cdo_get_stream_name(0), filename, argument); break;
    case (6):
      cdo_abort("ERROR (infile: '%s')! In parsing file '%s':\n          "
               "A Value for a keyword has a start quote sign but no end quote sign in line:\n          '%s'.\n"
               "          Check syntax.",
               cdo_get_stream_name(0), filename, argument); break;
    case (7):
      cdo_abort("ERROR (infile: '%s')! In parsing file '%s':\n          "
               "Unexpected separator sign '=' or ':' is found in values in line:\n          '%s'.",
               cdo_get_stream_name(0), filename, argument); break;
    case (8):
      cdo_abort("ERROR (infile: '%s')! In parsing file '%s':\n          "
               "Connected lines for one keyvalue contain more than the allowed 4096 characters.",
               cdo_get_stream_name(0), filename); break;
    case (9):
      cdo_abort("ERROR (infile: '%s')! In parsing file '%s':\n          "
               "A ',' is found at end of a line without information in next line.",
               cdo_get_stream_name(0), filename); break;
    case (10): cdo_abort("ERROR (infile: '%s')! In parsing file '%s':\n          A ',' is found at end of file.", cdo_get_stream_name(0), filename); break;
    case (11): cdo_abort("ERROR (infile: '%s')! In parsing file '%s':\n          A character is found after the end of a string terminated with \".", cdo_get_stream_name(0), filename); break;
    case (12): cdo_abort("ERROR (infile: '%s')! In parsing file '%s':\n          A string value is larger than the maximum size allowed by CMOR (1024 signs).", cdo_get_stream_name(0), filename); break;
    }
}

static int
copy_value(char *value, char **values, int *nvalues)
{
  if (*nvalues > 150) return 3;
  values[*nvalues] = strdup(value);
  (*nvalues)++;
  values[*nvalues] = nullptr;
  return 0;
}

static void
free_array(char **tofree)
{
  int i = 0;
  while (tofree[i])
    {
      free(tofree[i]);
      i++;
    }
  Free(tofree);
}

static void
quote_replace(char **values, int nvalues, int i)
{
  char *source = values[nvalues];
  char *useful;
  source++;
  source[i - 2] = 0;
  useful = strdup(source);
  free(values[nvalues]);
  values[nvalues] = strdup(useful);
  free(useful);
}

static char **
parse_string_to_values(char *workfile, char *pline, int *nvalues, const char *keyword)
{
  char **values = (char **) Malloc( 150 * sizeof(char*));
  while (isspace((int) *pline)) pline++;
  size_t len = strlen(pline);
  while (isspace((int) *(pline + (int) len))) len--;
  *(pline + len) = 0;
  if ((int) len == 0)
    {
      handleError(workfile, 4, keyword);
    }
  *nvalues = 0;
  size_t i = 0;
  int errh;
  while (i < len && len)
    {
      if (*(pline + i) == ',')
        {
          if (i == 0)
            {
              handleError(workfile, 5, pline);
            }
          errh = copy_value(pline, values, nvalues);
          if ( errh )
            handleError(workfile, errh, pline);
          if (*(values[*nvalues - 1] + i - 1) == '"' || *(values[*nvalues - 1] + i - 1) == '\'')
            quote_replace(values, *nvalues - 1, i);
          else
            *(values[*nvalues - 1] + i) = 0;
          
          if ( strlen(values[*nvalues-1]) > CMOR_MAX_STRING )
            {
              handleError(workfile, 12, pline);
            }

          i++;
          pline += i;
          len -= i;
          i = 0;
        }
      else if (*(pline + i) == '"')
        {
          i++;
          while (*(pline + i) != '"')
            {
              i++;
              if (*(pline + i) == 0)
                {
                  handleError(workfile, 6, pline);
                }
            }
          i++;
          if ( i != len && *(pline + i) != ',' && !isspace((int) *(pline + i)) && *(pline + i) != 0 && *(pline + i) != '\n' )
            {
              handleError(workfile, 11, pline);
            }
        }
      else if (isspace((int) *(pline + i)))
        break;
      else if (*(pline + i) == '=' )
        {
          cdo_warning("Value to parse '%s' contains '='."
                     " This will not be used as assignment.", pline);
          i++;
        }
      else
        i++;
    }
  errh = copy_value(pline, values, nvalues);
  if ( errh )
    handleError(workfile, errh, pline);
  if (*(values[*nvalues - 1] + i - 1) == '"')
    quote_replace(values, *nvalues - 1, i);
  else
    *(values[*nvalues - 1] + i) = 0;
  pline += i;
  if ( strlen(values[*nvalues-1]) > CMOR_MAX_STRING )
    {
      handleError(workfile, 12, pline);
    }
  return values;
}


static void
kv_insert_vals(KVList *kvl, const char *key, char *string, bool lreplace, bool lparse)
{
/* For internal model and institution check if string is not null */
  if ( string != nullptr )
    {
      std::string cppkey(key);
      int nvalues = 0;

      if (lreplace)
        {
          kvl->remove(cppkey);
        }
      if ( lparse )
        {
          const char **values = (const char **)parse_string_to_values(kv_get_a_val(kvl, "workfile4err", nullptr), string, &nvalues, key);
          kvl->append(key, values, nvalues);
        }
      else
        {
          char *values[1];
          values[0] = strdup(string);
          kvl->append(key, values, 1);
        }
    }
}

static std::vector<std::string> 
kv_get_vals(KVList *kvl, const char *key, int *numvals)
{
  std::vector<std::string> result;
  const KeyValues *kv = kvl->search(key);
  if (kv)
    {
      result = kv->values;
      *numvals = kv->nvalues;
    }
  return result;
}

PMList 
cdo_parse_cmor_file(const char *filename, bool lismap)
{
  PMList pml;

  auto fp = std::fopen(filename, "r");

  if (fp == nullptr)
    {
      if ( lismap )
        fprintf(stderr, "In reading the mapping table:\n          Open failed on %s: %s.", filename, strerror(errno));
      else
        fprintf(stderr, "In reading info table:\n          Open failed on %s: %s.", filename, strerror(errno));
      return pml;
    }

  ListBuffer listBuffer;
  auto status = listBuffer.read(fp, filename);
  if ( status && lismap) cdo_abort("Read error on mapping table %s!", filename);
  else if ( status ) cdo_abort("Read error on info table %s!", filename);

  std::fclose(fp);

  NamelistParser p;
  status = parse_list_buffer(p, listBuffer);
  if (status && lismap) cdo_abort("Mapping table not found!");
  else if (status) cdo_abort("Info table not found!");

  parse_namelist(pml, p, listBuffer.buffer.data(), true);

  return pml;
}

static const char *
check_short_key(char *key)
{
  const char *short_keys[] = { "cn", "n",   "c",  "u",  "cm", "vc", "p", "i",  "ca", "za",
                               "gi", "rtu", "mt", "om", "ms", "dr", "d", "lc", "dj", "kaa",
                               "member", "project_id", "vd", "dl",  "ta", "ci", "di", "tp",
                               "sc", "ml",
                               nullptr };
  const char *long_keys[]
      = { "cmor_name",     "name",        "code",           "units",    "cell_methods", "variable_comment",
          "positive",      "info",        "character_axis", "z_axis",   "grid_info",    "required_time_units",
          "mapping_table", "output_mode", "max_size",       "drs_root", "drs",          "last_chunk",
          "dataset_json",  "keep_all_attributes", "member", "project_id", "version_date", "deflate_level",
          "t_axis",  "climatology_interval", "decadal_interval", "tracking_prefix",
          "save_chunk",  "move_longitudes", nullptr };

  for (int i = 0; short_keys[i]; ++i)
    if (strcmp(key, short_keys[i]) == 0 || strcmp(key, long_keys[i]) == 0) return short_keys[i];
  return nullptr;
}

static void
map_it(KVList *kvl, int vlistID, int varID, const char *var2map)
{
  for (const auto &kv : *kvl)
    {
      const char *key = (check_short_key((char *) kv.key.c_str())) ? check_short_key((char *) kv.key.c_str()) : nullptr;
      if (!key)
        {
          if (kv.key != "project_mip_table")
            cdo_warning("In variable mapping:\n           you try to assign '%s' to variable '%s'.\n          This "
                       "mapping table keyword is skipped. Check allowed mapping table keywords.",
                       kv.key.c_str(), var2map);
          continue;
        }
      std::string cppkey = key;
      const char *value = kv.values[0].c_str();
      if (!value) continue;

      bool hasValidMin = false, hasValidMax = false;
      mapvar(vlistID, varID, kv, cppkey, nullptr, hasValidMin, hasValidMax, 0, false);
    }
}

static int
change_name_via_name(int vlistID, char *map_name, const char *cmor_name)
{
  char name[CDI_MAX_NAME];
  for (int varID = 0; varID < vlistNvars(vlistID); ++varID)
    {
      vlistInqVarName(vlistID, varID, name);
      if (strcmp(name, map_name) == 0)
        {
          cdiDefKeyString(vlistID, varID, CDI_KEY_NAME, parameter_to_word(cmor_name));
          return 1;
        }
    }
  return 0;
}

static int
change_name_via_code(int vlistID, char *map_code, const char *cmor_name)
{
  int codeproof = atol(map_code);
  if ( !codeproof || codeproof > 1000 )
    {
      cdo_print("'%s' could not be transformed into the code format (three digit integer). Code wont be used.", map_code);
      return 0;
    }

  int code;
  for (int varID = 0; varID < vlistNvars(vlistID); ++varID)
    {
      code = vlistInqVarCode(vlistID, varID);
      if ( codeproof == code )
        {
          cdiDefKeyString(vlistID, varID, CDI_KEY_NAME, parameter_to_word((const char *) cmor_name));
          return 1;
        }
    }
  return 0;
}

struct mapping
{
  int help_var;
  int cdi_varID;
  int cmor_varID;
  int zfactor_id;
  int charvars;
  char datatype;
  void *data;
};

static struct mapping *
construct_var_mapping(int vlistID)
{
  int nvars_max = vlistNvars(vlistID);
  struct mapping *vars = (struct mapping *) Malloc((nvars_max + 1) * sizeof(struct mapping));
  vars[0].cdi_varID = CDI_UNDEFID;
  vars[0].cmor_varID = CMOR_UNDEFID;
  vars[0].data = nullptr;
  vars[0].charvars = 0;
  return vars;
}

static void
destruct_var_mapping(struct mapping vars[])
{
  for (int i = 0; vars[i].cdi_varID != CDI_UNDEFID; ++i) Free(vars[i].data);
  Free(vars);
}

static struct mapping *
map_var(int cdi_varID, struct mapping vars[])
{
  for (int i = 0; vars[i].cdi_varID != CDI_UNDEFID; ++i)
    if (cdi_varID == vars[i].cdi_varID) return &vars[i];
  return nullptr;
}

static struct mapping *
new_var_mapping(struct mapping vars[])
{
  int i;
  for (i = 0; vars[i].cdi_varID != CDI_UNDEFID; ++i)
    ;
  vars[i + 1].cdi_varID = CDI_UNDEFID;
  vars[i + 1].cmor_varID = CMOR_UNDEFID;
  vars[i + 1].data = nullptr;
  vars[i + 1].charvars = 0;
  return &vars[i];
}

static int *
new_axis_id(int *axis_ids)
{
  int i;
  for (i = 0; axis_ids[i] != CMOR_UNDEFID; ++i)
    ;
  axis_ids[i + 1] = CMOR_UNDEFID;
  return &axis_ids[i];
}

static int
count_axis_ids(int *axis_ids)
{
  int i;
  for (i = 0; axis_ids[i] != CMOR_UNDEFID; ++i)
    ;
  return i;
}


static const KVList *
check_for_charvars(KVList *cmorVarLine, const char */*cn*/, char */*miptabfreq*/, const char *key)
{
  /***/
  /* If a infile variable selector (name or code) has more than one value, it must be a character coordinate*/
  /* If it is given as a string and the string contains a ',', */
  /* it must be divided into several values and is a variable with character coordinate */
  /***/
  if (Options::cdoVerbose) cdo_print("Start to check for variables to merge as a character coordinate.");
  const KeyValues *kvn = nullptr;
  bool isCode = false;
  if (key)
    kvn = cmorVarLine->search(key);
  else
    {
      kvn = cmorVarLine->search("name");
      if (!kvn)
        {
          kvn = cmorVarLine->search("code");
          if ( kvn )
            isCode = true;
        }
    }
  if (kvn && kvn->nvalues > 1 ) return cmorVarLine;

  if (kvn && strstr(kvn->values[0].c_str(), ",") && kvn->nvalues == 1)
    {
      if (Options::cdoVerbose) cdo_print("Start to replace identifier string with its values.");
      char *unfilteredComma = strdup(kvn->values[0].c_str());
      if ( isCode )
        kv_insert_vals(cmorVarLine, "code", unfilteredComma, true, true);
      else
        kv_insert_vals(cmorVarLine, "name", unfilteredComma, true, true);
      Free(unfilteredComma);
      if (Options::cdoVerbose) cdo_print("Successfully replaced identifier string with its values.");
      return cmorVarLine;

    }
  if (Options::cdoVerbose) cdo_print("Successfully checked for variables to merge as a character coordinate.");
  return nullptr;
}

static void
addcharvar(const KeyValues *charvars, int vlistID, const char *key, struct mapping vars[])
{
  MergeVarsOnAxis withnewcharaxis;
  withnewcharaxis.inputNames = *charvars;

  if (Options::cdoVerbose) cdo_print("Start to merge variables to one character coordinate.");
  int nvars = vlistNvars(vlistID);
  for (int i = 0; i < charvars->nvalues; ++i)
    {
      MergeVarKeys temp;
      temp.vlistID = vlistID;
      temp.varID = getVarIDToMap(vlistID, nvars, key, withnewcharaxis.inputNames.values[i].c_str());
      if (temp.varID == CDI_UNDEFID)
        cdo_abort("ERROR (infile: '%s')! In merging variables to a variable with a character coordinate:\n          Could not find "
                 "'%s' in infile '%s' to build a variable with character coordinate.",
                 cdo_get_stream_name(0), withnewcharaxis.inputNames.values[i].c_str(), cdo_get_stream_name(0));
      temp.gridID = vlistInqVarGrid(vlistID, temp.varID);
      temp.projID = gridInqProj(temp.gridID);
      temp.zaxisID = vlistInqVarZaxis(vlistID, temp.varID);
      withnewcharaxis.inputKeys.push_back(temp);
    }

  std::vector<int> axissize(3);
  axissize[0] = gridInqXsize(withnewcharaxis.inputKeys[0].gridID);
  axissize[1] = gridInqYsize(withnewcharaxis.inputKeys[0].gridID);
  axissize[2] = zaxisInqSize(withnewcharaxis.inputKeys[0].zaxisID);

  int oldgridsize = gridInqSize(withnewcharaxis.inputKeys[0].gridID);
  if ( axissize[0] * axissize[1] == 0 )
    oldgridsize = 1;

  if (axissize[0] != 1 && axissize[1] != 1 && axissize[2] != 1)
    {
      if (Options::cdoVerbose)
        cdo_print("In merging variables to an axis:\n          Spatial dimensions are already set. A fourth axis is created.");
      char ids[CDI_MAX_NAME];
      sprintf(ids, "%d", withnewcharaxis.inputKeys[0].varID);
      for (int i = 1; i < charvars->nvalues; ++i)
        {
          char tempint[sizeof(int)];
          sprintf(tempint,"%d",withnewcharaxis.inputKeys[i].varID);
          strcat(ids,",");
          strcat(ids,tempint);
        }
      cdiDefAttTxt(vlistID, withnewcharaxis.inputKeys[0].varID, "merge_axis", (int) strlen(ids), ids);
      return;
    }
  int ntsteps = vlistNtsteps(vlistID);

  if (cdo_assert_files_only() == false)
    cdo_abort("ERROR (infile: '%s')! Cannot merge variables to a new axis because the input file needs to be opened twice but you piped operators so the input file is not clear.", cdo_get_stream_name(0));

  const auto streamID2 = streamOpenRead(cdo_get_stream_name(0));
  if (ntsteps == -1)
    {
      ntsteps = 0;
      int dummy;
      while ((dummy = streamInqTimestep(streamID2, ntsteps++)))
        ;
    }

  axissize = withnewcharaxis.define_new_axes(axissize);
  withnewcharaxis.define_var_structure(vlistID, ntsteps, axissize);

  struct mapping *var = new_var_mapping(vars);
  var->cdi_varID = withnewcharaxis.output.varID;
  var->datatype = withnewcharaxis.output.datatype;

  withnewcharaxis.read_cmor_charvar(axissize, streamID2, oldgridsize);

  var->data = Malloc(ntsteps * axissize[0] * axissize[1] * axissize[2] * sizeof(double));
  for ( int i = 0; i < ntsteps * axissize[0] * axissize[1] * axissize[2]; i++ )
    {
      if ( withnewcharaxis.output.datatype == 'd' )
        ((double *) var->data)[i] = withnewcharaxis.data[i];
      else
        ((float *) var->data)[i] = (float)withnewcharaxis.data[i];

    }
  var->charvars = 1;

  streamClose(streamID2);

  if (Options::cdoVerbose)
    cdo_print(
        "Successfully merged variables into one character axis. The final variable is called '%s' and has the ID: '%d'",
        charvars->values[0].c_str(), var->cdi_varID);
}

static int
maptab_via_key(char *tablename, PMList pml, int vlistID, int varID, const char *key, char *miptabfreq)
{
  char ifilevalue[CDI_MAX_NAME];
  if ((strcmp(key, "cmor_name") == 0) || (strcmp(key, "name") == 0) )
    vlistInqVarName(vlistID, varID, ifilevalue);
  else
    get_stringcode(vlistID, varID, ifilevalue);

  if (ifilevalue[0])
    {
      KVList kvl = maptab_search_miptab(pml, ifilevalue, miptabfreq, key);
      if (kvl.size())
        {
          if (Options::cdoVerbose) cdo_print("Start to map via '%s'.", key);
          map_it(&kvl, vlistID, varID, ifilevalue);
          return 1;
        }
      if (Options::cdoVerbose)
        cdo_print("In variable mapping with table '%s':\n          "
                 "Variable named '%s' with varID '%d' could not be mapped via '%s' "
                 "because no corresponding key '%s' was found in mapping table file.",
                 tablename, ifilevalue, varID, key, key);
      return 0;
    }
  else
    {
      if (Options::cdoVerbose)
        cdo_print("In variable mapping with table '%s':\n          "
                 "Variable with varID '%d' could not be mapped via '%s' because it "
                 "does not possess a '%s' in infile.",
                 tablename, varID, key, key);
      return 0;
    }
}

static int
maptab_via_cn_and_key(KVList *kvl_oname, int vlistID, int nvars, const char *key)
{
  const KeyValues *kv = kvl_oname->search(key);
  const KeyValues *kvcn = kvl_oname->search("cmor_name");
  if (kv)
    {
      int varID = (strcmp(key, "cmor_name") == 0) ? getVarIDToMap(vlistID, nvars, "name", kv->values[0].c_str())
                                                  : getVarIDToMap(vlistID, nvars, key, kv->values[0].c_str());
      int newvar = getVarIDToMap(vlistID, nvars, "name", kvcn->values[0].c_str());
      if (varID != CDI_UNDEFID)
        {
          if ( newvar != CDI_UNDEFID && newvar != varID)
            {
              cdo_warning("In variable mapping : \n          "
                  "You try to define a variable '%s' that is already in the input stream.\n"
                 "The already existing infile variable will be renamed.\n"
                 "This may lead to errors. Choose the respective var with '-selname' before applying cmor.",
                    kvcn->values[0].c_str());
              cdiDefKeyString(vlistID, newvar, CDI_KEY_NAME, "overwritten");          
            }
          if (Options::cdoVerbose) cdo_print("Started mapping of variable via '%s'.", key);
          map_it(kvl_oname, vlistID, varID, kv->values[0].c_str());
          return 1;
        }
      cdo_print("In variable mapping:\n          Variable '%s' configured via cmor_name\n          could not be mapped "
               "via key '%s' because no infile variable '%s' equals '%s'.",
               kvcn->values[0].c_str(), key, key, kv->values[0].c_str());
    }
  else if (Options::cdoVerbose)
    cdo_print("In variable mapping:\n          Variable '%s' configured via cmor_name\n          could not be mapped "
             "via key '%s' because it possesses no corresponding key '%s' in mapping file.",
             kvcn->values[0].c_str(), key, key);
  return 0;
}

static void
maptab_via_cmd(char *tablename, PMList pml, const char *origValue, int vlistID, const char *key, char *cmorName,
               char *miptabfreq, int filetype, char *maptab)
{
  KVList cmorVarLine;

  int nvars = vlistNvars(vlistID);
  int varIDToMap = getVarIDToMap(vlistID, nvars, key, origValue);
  if (varIDToMap == CDI_UNDEFID)
    cdo_abort("ERROR (infile: '%s')! In variable mapping with table '%s':\n          "
                   "Variable with '%s': '%s' configured via cmdline could not be "
                   "found in infile '%s'.",
                    cdo_get_stream_name(0), tablename, key, origValue, cdo_get_stream_name(0));

  cmorVarLine = maptab_search_miptab(pml, cmorName, miptabfreq, "cmor_name");
  const KeyValues *kvcn = cmorVarLine.search("cmor_name");
  int newvar = getVarIDToMap(vlistID, nvars, "name", kvcn->values[0].c_str());
  if ( newvar != CDI_UNDEFID && newvar != varIDToMap )
    {
      cdo_warning("In variable mapping : \n          "
                  "You try to define a variable '%s' that is already in the input stream.\n"
                 "The already existing infile variable will be renamed.\n"
                 "This may lead to errors. Choose the respective var with '-selname' before applying cmor.",
                    kvcn->values[0].c_str());
      cdiDefKeyString(vlistID, newvar, CDI_KEY_NAME, "overwritten");
    }
  if ( !cmorVarLine.size() )
    {
      cdo_warning("In variable mapping with table '%s':\n          "
                 "The registered cmor_name '%s' via cmdline could not be found in "
                 "mapping table.\n          No mapping table is applied.",
                 tablename, cmorName);
      cdiDefKeyString(vlistID, varIDToMap, CDI_KEY_NAME, parameter_to_word((const char *) cmorName));
    }
  else
    {
      if ( (filetype == CDI_FILETYPE_GRB || filetype == CDI_FILETYPE_GRB2) && strcmp(key, "name") == 0 )
        {
          if ( Options::cdoVerbose ) cdo_print("5.1. In applying the mapping table:\n          Note that you use 'name' as selector keyword "
                     "allthough the type of infile is GRB.");
        }
      if (Options::cdoVerbose) cdo_print("Started mapping of variable via '%s'.", key);
        map_it(&cmorVarLine, vlistID, varIDToMap, cmorName);

      if (Options::cdoVerbose) cdo_print("5. Successfully found, read and applied mapping table '%s'.", maptab);
    }
}

static void
maptab_via_cn(char *tablename, PMList pml, std::vector<std::string> request, int vlistID, int numvals, 
              char *miptabfreq, int filetype, struct mapping vars[], bool isWarn)
{
  for (int j = 0; j < numvals; ++j)
    {
      KVList kvl_oname = maptab_search_miptab(pml, request[j].c_str(), miptabfreq, "cmor_name");
      if (kvl_oname.size())
        {
          if ( vars && check_for_charvars(&kvl_oname, request[j].c_str(), miptabfreq, nullptr) )
            {
              const KeyValues *charcode = kvl_oname.search("name");
              if ( !charcode )
                {
                  charcode = kvl_oname.search("code");
                  addcharvar(charcode, vlistID, "code", vars);
                }
              else
                addcharvar(charcode, vlistID, "name", vars);
            }
          int nvars = vlistNvars(vlistID);
          if (filetype == CDI_FILETYPE_GRB || filetype == CDI_FILETYPE_GRB2)
            {
              if (maptab_via_cn_and_key(&kvl_oname, vlistID, nvars, "code"))
                {
                  if (Options::cdoVerbose) cdo_print("Successfully mapped variable via code to cmor_name '%s'.", request[j]);
                  continue;
                }
            }
          if (maptab_via_cn_and_key(&kvl_oname, vlistID, nvars, "name"))
            {
              if (Options::cdoVerbose) cdo_print("Successfully mapped variable via name to cmor_name: '%s'.", request[j]);
              continue;
            }
          else
            {
              if (Options::cdoVerbose)
                cdo_print("In variable mapping with table '%s':\n          "
                         "Try to use cmor_name for selecting the infile variable.",
                         tablename, request[j]);
              if (maptab_via_cn_and_key(&kvl_oname, vlistID, nvars, "cmor_name"))
                {
                  cdo_print("Successfully mapped variable via cmor_name to cmor_name '%s'.", request[j]);
                  continue;
                }
              if ( isWarn)
                cdo_warning("In variable mapping with table '%s':\n          "
                         "Mapping table line of cmor_name '%s' could neither be mapped "
                         "via 'name', 'code' nor 'cmor_name'.\n          No mapping for cmor_name: '%s'.",
                         tablename, request[j], request[j]);
              continue;
            }
        }
      else
        {
          if ( isWarn)
            cdo_warning("In variable mapping with table '%s':\n          "
                     "Requested cmor_name: '%s' is not found in mapping table.\n       "
                     "   No mapping for cmor_name: '%s'",
                     tablename, request[j], request[j]);
          continue;
        }
    }
}

/*
static char *
trim(char *s)
{
  if (s == nullptr) return s;
  while (*s != '\0' && (isspace(*s) || *s == '"')) s++;
  int n = strlen(s);
  while (n > 0 && (isspace(s[n - 1]) || s[n - 1] == '"')) n--;
  s[n] = '\0';
  return s;
}
*/

static bool
file_exist(const char *tfilename, bool force, const char *fileart, bool print)
{
  assert(tfilename != nullptr);
  size_t filesize = FileUtils::size(tfilename);
  if (filesize == 0 && force)
    cdo_abort("ERROR (infile: '%s')!\n          Empty '%s' file: '%s'.", cdo_get_stream_name(0), fileart, tfilename);
  else if (filesize == 0 && !force)
    {
      if (print && Options::cdoVerbose) cdo_print("Empty '%s' file: '%s'.", fileart, tfilename);
      return false;
    }
  if (strstr(tfilename, ".nc") || strstr(tfilename, ".grb")) return 1;
  FILE *fp = std::fopen(tfilename, "r");
  if (fp == nullptr && force)
    cdo_abort("ERROR (infile: '%s')!\n           Open failed on '%s' file: '%s'.", cdo_get_stream_name(0), fileart, tfilename);
  else if (fp == nullptr && !force)
    {
      if (print && Options::cdoVerbose )
        cdo_print("Open failed on '%s' file: '%s'.", fileart, tfilename);
      return false;
    }

  std::fclose(fp);
  return true;
}

static int
parse_kv_file(KVList *kvl, const char *filename)
{
  PMList pmkv = cdo_parse_cmor_file(filename, false);
  if ( !pmkv.size() )
    return 0;
  for ( auto &kv : pmkv.front() )
    {
      char *keyname = strdup(kv.key.c_str());
      const KeyValues *kvfromlist = kvl->search(kv.key.c_str());
      if ( kvfromlist )
        continue;

      char **values = (char **) Malloc((kv.nvalues+1) * sizeof(char *));
      int k = 0;
      for ( k = 0; k<kv.nvalues; k++ )
        values[k] = strdup(kv.values[k].c_str());

      kvl->append(keyname, values, kv.nvalues);
    }
  return 1;
}

static void
check_compare_set(char **finalset, char *attribute, const char *attname, const char *defaultstr)
{
  if (!(*finalset))
    {
      if (!attribute)
        {
          if (defaultstr)
            *finalset = strdup(defaultstr);
          else
            cdo_abort("ERROR (infile: '%s')! Function 'attribute check, compare and set':\n          Required value for attribute '%s' "
                     "is neither found in input file nor in the configuration.",
                     cdo_get_stream_name(0), attname);
        }
      else
        *finalset = strdup(attribute);
    }
  else if (attribute)
    {
      if (strcmp(attribute, *finalset) != 0 && strcmp(attribute, "") != 0 )
        {
          if (Options::cdoVerbose)
            cdo_print("Function 'attribute check, compare and set':\n          Be aware of differences between infile and "
                   "user specification.\n          Attribute '%s' in input file: '%s' does not agree with user "
                   "specification '%s'.",
                   attname, *finalset, attribute);
          cdo_print("Attribute '%s' = '%s'", attname, attribute);
          Free(*finalset);
          *finalset = strdup(attribute);
        }
    }
}

static char *
get_infile_attvalue(int vlistID, int varID, char *name, int type, int len)
{
  char *infile_attvalue = (char *) Malloc(CMOR_MAX_STRING * sizeof(char));
  if (type == CDI_DATATYPE_INT32)
    {
      int *values = (int *) Malloc(len * sizeof(int));
      cdiInqAttInt(vlistID, varID, name, len, &values[0]);
      sprintf(infile_attvalue, "%i", values[0]);
      for (int l = 1; l < len; ++l)
        {
          char tempint[sizeof(values[l])];
          sprintf(tempint,"%i",values[l]);
          strcat(infile_attvalue, tempint);
        }
      Free(values);
    }
  else if (type == CDI_DATATYPE_FLT32 || type == CDI_DATATYPE_FLT64)
    {
      char fltstr[128];
      std::vector<double> attflt(len);
      cdiInqAttFlt(vlistID, varID, name, len, attflt.data());
      for (int i = 0; i < len; ++i)
        {
          if (i) strcat(infile_attvalue, ", ");
          /*strcat(infile_attvalue, i ? "," : " ");*/
          char tempflt[64];
          if (type == CDI_DATATYPE_FLT32)
            {
              sprintf(tempflt,
                      "%sf",
                      double_to_att_str(Options::CDO_flt_digits, fltstr, sizeof(fltstr), attflt[i]));
            }
          else
            {
              sprintf(tempflt,
                      "%s",
                      double_to_att_str(Options::CDO_dbl_digits, fltstr, sizeof(fltstr), attflt[i]));
            }
          if (i)
            strcat(infile_attvalue, tempflt);
          else
            sprintf(infile_attvalue, "%s", tempflt);
        }
      infile_attvalue[CMOR_MAX_STRING] = '\0';
      cdo_print("%s", infile_attvalue);
    }
  else
    {
      cdiInqAttTxt(vlistID, varID, name, len, infile_attvalue);
      infile_attvalue[len] = '\0';
    }
  return infile_attvalue;
}

static char *
get_infile_attname(int vlistID, int varID, int natt, int *type, int *len)
{
  char *infile_attname = (char *) Malloc(CDI_MAX_NAME * sizeof(char));
  cdiInqAtt(vlistID, varID, natt, infile_attname, type, len);
  return infile_attname;
}

static char *
get_txtatt(int vlistID, int varID, const char *key)
{
  int natts;
  char *returnvalue = nullptr;
  cdiInqNatts(vlistID, varID, &natts);
  for (int i = 0; i < natts; ++i)
    {
      int type, len;
      char *infile_attname = get_infile_attname(vlistID, varID, i, &type, &len);
      if (strcmp(infile_attname, key) == 0)
        {
          returnvalue = get_infile_attvalue(vlistID, varID, infile_attname, type, len);
        }
      Free(infile_attname);
    }
  return returnvalue;
}

static void
get_all_atts(KVList *kvl, int vlistID, const char **infileAttSpecLong, const char **infileAttSpecShort)
{
  int natts, type = 0, len = 0;
  cdiInqNatts(vlistID, CDI_GLOBAL, &natts);
  for (int i = 0; i < natts; ++i)
    {
      char *infile_attname = get_infile_attname(vlistID, CDI_GLOBAL, i, &type, &len);

      if ( strcmp(infile_attname, "history") == 0 )
        continue;
      int j = 0;
      while (infileAttSpecLong[j] != nullptr)
        {
          if (strcmp(infile_attname, infileAttSpecLong[j]) == 0)
            {
              kv_insert_vals(kvl, infileAttSpecShort[j], get_infile_attvalue(vlistID, CDI_GLOBAL, infile_attname, type, len),
                             false, true);
              break;
            }
          j++;
        }
      if (infileAttSpecLong[j] == nullptr)
        {
          kv_insert_vals(kvl, infile_attname, get_infile_attvalue(vlistID, CDI_GLOBAL, infile_attname, type, len), false, false);
        }
      Free(infile_attname);

    }

  kv_insert_vals(kvl, "institution", (char *)institutInqLongnamePtr(vlistInqInstitut(vlistID)), false, false);
/*  if ( strcmp(project_id, "CMIP6") != 0 ) */
    kv_insert_vals(kvl, "source", (char *) modelInqNamePtr(vlistInqModel(vlistID)), false, false);
}

static void
add_globalhybrids(KVList *kvl, int vlistID)
{
/* Do not define history, source or institution as one of those: */
  const char *infileAttSpecLong[] = { "branch_dates", "climatology_interval", "decadal_interval", nullptr };
  const char *infileAttSpecShort[] = { "branch_dates", "ci", "di", nullptr };
  int i = 0;
  if (strcmp(kv_get_a_val(kvl, "kaa", "n"), "y") == 0 || strcmp(kv_get_a_val(kvl, "keep_all_attributes", "n"), "y") == 0)
    get_all_atts(kvl, vlistID, infileAttSpecLong, infileAttSpecShort);
  else
    {
      while (infileAttSpecLong[i] != nullptr)
        {
          char *att = get_txtatt(vlistID, CDI_GLOBAL, infileAttSpecLong[i]);
          if (att) kv_insert_vals(kvl, infileAttSpecShort[i], att, false, false);
          if (att) Free(att);
          i++;
        }

      const char *infileAtt[] = { "member", "variant_label", nullptr };
      i = 0;
      while (infileAtt[i] != nullptr)
        {
          char *att = get_txtatt(vlistID, CDI_GLOBAL, infileAtt[i]);
          if (att) kv_insert_vals(kvl, infileAtt[i], att, false, false);
          i++;
          if (att) Free(att);
        }
    }
  const char *longAtt[] = { "required_time_units", "grid_info"     , "mapping_table", "keep_all_attributes",
                            "drs_root"           , "drs"           , "t_axis"  , "version_date",
                            "deflate_level"      , "character_axis", "z_axis"       , "output_mode",
                            "max_size"           , "last_chunk"    , "save_chunk",  "move_longitudes", nullptr };
  const char *shortAtt[] = { "rtu", "gi", "mt", "kaa", "dr", "d", "ta", "vd", "dl", "ca", "za", "om", "ms", "lc", "sc", "ml", nullptr };
  i = 0;
  while (longAtt[i] != nullptr)
    {
      for (auto &kv_latt : *kvl)
        {
          if (kv_latt.key == longAtt[i])
            {
              kv_insert_vals(kvl, shortAtt[i], (char *) kv_latt.values[0].c_str(), false, false);
              kvl->remove(kv_latt.key);
              break;
            }
        }
      i++;
    }
}

static int
check_attarray(KVList *kvl, const char **reqAtt, int vlistID)
{
  int i = 0;
  while (reqAtt[i] != nullptr)
    {
      const KeyValues *kv_reqatt = kvl->search(reqAtt[i]);
      if (!kv_reqatt || kv_reqatt->nvalues == 0 )
        {
          if (strcmp(kv_get_a_val(kvl, "kaa", "y"), "y") != 0)
            {
              char *infileatt = get_txtatt(vlistID, CDI_GLOBAL, reqAtt[i]);
              if (infileatt && strcmp(infileatt, "notSet") != 0)
                {
                  if (Options::cdoVerbose) cdo_print("Attribute (from infile) '%s' is '%s'. ", reqAtt[i], infileatt);
                  const char *infileatts[] = { infileatt };
                  kvl->append(reqAtt[i], infileatts, 1);
                  Free(infileatt);
                }
              else if ( infileatt )
                {
                  Free(infileatt);
                }
              return i;
            }
          else
            return i;
        }
      else if (kv_reqatt->values[0] == "notSet")
        return i;
      else if (Options::cdoVerbose)
        cdo_print("Attribute (from meta data file) '%s' is '%s'. ", reqAtt[i], kv_reqatt->values[0]);
      i++;
    }
  return -1;
}

static void
attErr(const char **reqAtt, int errnum)
{
  char errStr[CMOR_MAX_STRING];
  int i = 1;

  sprintf(errStr, "ERROR! Attribute '%s' is required. Either it is missing, 'notSet', or the value is invalid.\n       "
                  "   Make sure that you have configured all following attributes:\n   "
                  "       %s",
          reqAtt[errnum], reqAtt[0]);
  while (reqAtt[i])
    {
      strcat(errStr, ", ");
      strcat(errStr,reqAtt[i]);
      i++;
    }
  cdo_abort(errStr);
}

static void
check_attr(KVList *kvl, char *project_id, int vlistID)
{
  /* Project id moved to main void fct */
  const char *reqAtt[] = { "institution", "source", "experiment_id", "rtu", nullptr };
  const char *reqAttCMIP5[] = { "institute_id", "product", "member", nullptr };
  const char *reqAttCORDEX[]
      = { "institute_id", "product", "member", "CORDEX_domain", "driving_model_id", "rcm_version_id", nullptr };
  /*  const char *reqAttCMIP6CMOR3[] = {"outpath", "output_path_template", "output_file_template", "tracking_prefix",
  nullptr};
  "further_info_url", "variant_label",*/
  const char *reqAttCMIP6[] = { "activity_id", "experiment",         "grid",      "grid_label",  "institution_id",
                                "license",     "nominal_resolution", "source_id", "source_type", nullptr };
  const char *expdepAttCMIP6[]
      = { "parent_experiment_id", "parent_activity_id", "parent_mip_era",    "parent_source_id", "parent_variant_label",
          "parent_time_units",    "sub_experiment",     "sub_experiment_id", "branch_method",    nullptr };
  /* In all Projects needed Attributes are tested first */

  int errnum = 0;
  if ((errnum = check_attarray(kvl, reqAtt, vlistID)) != -1) attErr(reqAtt, errnum);

#if (CMOR_VERSION_MAJOR == 2)
  /* Set default attributes */
  const char *reqAttCMOR2[] = { "contact", "model_id", nullptr };
  if ((errnum = check_attarray(kvl, reqAttCMOR2, vlistID)) != -1) attErr(reqAttCMOR2, errnum);

  const KeyValues *kv = kvl->search("references");
  if (!kv || kv->nvalues == 0 || kv->values[0] == "notSet")
    {
      const KeyValues *kv_model_id = kvl->search("model_id");
      char *references = (char *) Malloc(strlen(kv_model_id->values[0].c_str()) + 28);
      strcpy(references, "No references available for ");
      strcat(references, kv_model_id->values[0].c_str());
      cdo_print("Attribute 'references' is set to '%s' ", references);
      kv_insert_vals(kvl, "references", references, false, false);
      Free(references);
    }
#endif
  /* Special check for CMIP or CORDEX projects */
  if (strcmp(project_id, "CMIP5") == 0)
    {
      if (Options::cdoVerbose) cdo_print("Since the project id is CMIP5 further attributes are tested. ");
      if ((errnum = check_attarray(kvl, reqAttCMIP5, vlistID)) != -1) attErr(reqAttCMIP5, errnum);
#if (CMOR_VERSION_MAJOR == 3)
      const char *reqAttCMIP5CMOR3[] = { "modeling_realm", nullptr };
      /**************/
      /* Add additional attributes for CMIP5 */
      /* allthough using CMOR 3 */
      /**************/
      if ((errnum = check_attarray(kvl, reqAttCMIP5CMOR3, vlistID)) != -1) attErr(reqAttCMIP5CMOR3, errnum);
#endif
    }
  else if (strcmp(project_id, "CORDEX") == 0)
    {
      if (Options::cdoVerbose) cdo_print("Since the project id is CORDEX further attributes are tested.");
      if ((errnum = check_attarray(kvl, reqAttCORDEX, vlistID)) != -1) attErr(reqAttCORDEX, errnum);
    }
  else if (strcmp(project_id, "CMIP6") == 0)
    {
      kv_insert_vals(kvl, "_cmip6_option", (char *) "CMIP6", true, false);
      if (Options::cdoVerbose) cdo_print("Since the project_id is CMIP6 further attributes are tested.");
      if ((errnum = check_attarray(kvl, reqAttCMIP6, vlistID)) != -1) attErr(reqAttCMIP6, errnum);
      int j = 0;
      char *pei = nullptr;
      if ( strcmp(kv_get_a_val(kvl, "parent_experiment_id", "no parent"), "no parent") != 0 )
        {
          pei = kv_get_a_val(kvl, "parent_experiment_id", nullptr);
          cdo_print("Since you set attribute 'parent_experiment_id'='%s', further attributes are checked.", pei);
        }
      while (expdepAttCMIP6[j] != nullptr)
        {
          const KeyValues *kv_reqatt = kvl->search(expdepAttCMIP6[j]);
          if (!kv_reqatt || kv_reqatt->nvalues == 0 || kv_reqatt->values[0] == "notSet")
            {
              if (Options::cdoVerbose || pei)
                cdo_print("Depending on the experiment, attribute '%s' may be required. Either it is missing or notSet",
                         expdepAttCMIP6[j]);
            }
          else if (Options::cdoVerbose)
            {
              cdo_print("Attribute (from meta data file) '%s' is '%s'. ", expdepAttCMIP6[j], kv_reqatt->values[0]);
            }
          j++;
        }
    }
#if (CMOR_VERSION_MAJOR == 3)
  else
    {
      const char *customProjectCMOR3[]= {"_controlled_vocabulary_file",
                                         "_FORMULA_VAR_FILE",
                                         "_AXIS_ENTRY_FILE",
                                         "_history_template",
                                         "_further_info_url_tmpl",
                                         "output_path_template", 
                                         "output_file_template"};
      const char *customProjectCMOR3default[]= {"CMIP6_CV.json",
                                         FORMULA_VAR_FILENAME,
                                         AXIS_ENTRY_FILENAME,
                                         CMOR_DEFAULT_HISTORY_TEMPLATE,
                                         CMOR_DEFAULT_FURTHERURL_TEMPLATE,
                                         CMOR_DEFAULT_PATH_TEMPLATE,
                                         CMOR_DEFAULT_FILE_TEMPLATE};
      const KeyValues *kv = kvl->search("_controlled_vocabulary_file");
      if (!kv || kv->nvalues == 0 || kv->values[0] == "notSet")
        {
          if (Options::cdoVerbose )
                cdo_print("Since you have specified a project different than CMIP and you use CMOR version 3, we recommend that you also specify the following variables: '%s' '%s' '%s'. The default values are: '%s '%s' '%s' and need to be available in the MIP-tables directory.",
                         customProjectCMOR3[0],
                         customProjectCMOR3[1],
                         customProjectCMOR3[2],
                         customProjectCMOR3default[0],
                         customProjectCMOR3default[1],
                         customProjectCMOR3default[2]);
        }
    }
#endif
}

static void
check_mem(KVList *kvl, char *project_id)
{
/*Check if both is registered */
  std::vector<std::string> ensindexCMIP5 = { "realization", "initialization_method", "physics_version" };
  std::vector<std::string> ensindexCMIP6 = { "realization_index", "initialization_index", "physics_index", "forcing_index" };
  std::string ripf = "ripf";
  std::string rip = "rip";

  if (strcmp(project_id, "CMIP5") == 0 ||
      strcmp(project_id, "CORDEX") == 0)
    {
      char *cm = kv_get_a_val(kvl, "cm", " ");

      if ( cm[0] == 'n' && (strcmp(project_id, "CMIP5") == 0 || strcmp(project_id, "CORDEX") == 0 ) )
        {
          if ( strcmp(project_id, "CORDEX") == 0 )
            kv_insert_vals(kvl, "driving_model_ensemble_member", (char *)"r0i0p0", true, false);
          kv_insert_vals(kvl, "member", (char *)"r0i0p0", true, false);
        }
      else
        {
           char *vlabel = kv_get_a_val(kvl, "member", nullptr);
           bool lanyindex = false, lallindices = true;
           long int indexvalues[3] = { 0 };
           int indexint = 0;

           for (std::vector<std::string>::iterator index=ensindexCMIP5.begin(); index!=ensindexCMIP5.end(); ++index)
             {
               char *indexstring = kv_get_a_val(kvl, index->c_str(), nullptr);
               if ( indexstring )
                 {
                   indexvalues[indexint] = atol((const char *) indexstring);
                   if ( !indexvalues[indexint] )
                     {
                       cdo_warning("Could not parse value '%s' of attribute '%s' to integer.", indexstring, index->c_str());
                       lallindices = false;
                     }
                   else
                     lanyindex = true;
                 }
               else
                 lallindices = false;
               indexint++;
             }
           if ( lanyindex && vlabel )
             cdo_warning("You specified both variant_label and all indices attributes.\n"
                       "Note that indices attributes have higher priority.");
           else if ( !lallindices && vlabel == nullptr )
             cdo_abort("ERROR (infile: '%s')!\n"
                 "Could not find all index values required for describing the ensemble member (rip).\n"
                 "Make sure you have specified either 'variant_label' or all 3 'rip' indexes.", cdo_get_stream_name(0));

           if ( !lanyindex )
             {
               int scanres = std::sscanf((const char *)vlabel, "r%ldi%ldp%ld",
                               &indexvalues[0], &indexvalues[1], &indexvalues[2]);
               if ( !scanres )
                 cdo_abort("ERROR (infile: '%s')!\n"
                 "Could not scan all integers from attribute 'member'. \n"
                     "Make sure it has the format 'rINTiINTpINT'", cdo_get_stream_name(0));
               indexint = 0;
               for (std::vector<std::string>::iterator index=ensindexCMIP5.begin(); index!=ensindexCMIP5.end(); ++index)
                 {
                   std::string index2string = std::to_string(indexvalues[indexint]);
                   kv_insert_vals(kvl, index->c_str(), (char *)index2string.c_str(), true, false);
                   indexint++;
                 }
             }
           else if ( !lallindices )
             {
               long int vlvalues[3] = { 0 };
               int scanres = std::sscanf((const char *)vlabel, "r%ldi%ldp%ld",
                               &vlvalues[0], &vlvalues[1], &vlvalues[2]);
              if ( !scanres )
                cdo_warning("Could not scan all integers from attribute 'member'. \n"
                     "Make sure it has the format 'rINTiINTpINT'");
               indexint = 0;
               while ( indexint < 3 )
                 {
                   if ( !indexvalues[indexint] )
                     {
                       indexvalues[indexint] = vlvalues[indexint];
                       if ( !indexvalues[indexint] )
                         cdo_abort("ERROR (infile: '%s')!\n"
                           "You did not provide a value for attribute '%s'.\n"
                          "Make sure you have specified either 'variant_label' or all 3 'rip' indexes.",
                          cdo_get_stream_name(0) , ensindexCMIP5[indexint]);
                       std::string index2string = std::to_string(indexvalues[indexint]);
                       kv_insert_vals(kvl, ensindexCMIP5[indexint].c_str(), (char *)index2string.c_str(), true, false);
                     }
                   indexint++;
                 }
             }
           char member[CMOR_MAX_STRING];
           sprintf(member, "r%ldi%ldp%ld", indexvalues[0], indexvalues[1], indexvalues[2]);
           kv_insert_vals(kvl, "member", member, true, false);
           kv_insert_vals(kvl, "variant_label", member, true, false);
           indexint = 0;
        }
    }
  else
    {
      char *vlabel = kv_get_a_val(kvl, "variant_label", nullptr);
      bool lanyindex = false, lallindices = true;
      long int indexvalues[4] = { 0 };
      int indexint = 0;

      for (std::vector<std::string>::iterator index=ensindexCMIP6.begin(); index!=ensindexCMIP6.end(); ++index)
        {
          char *indexstring = kv_get_a_val(kvl, index->c_str(), nullptr);
          if ( indexstring )
            {
              indexvalues[indexint] = atol((const char *) indexstring);
              if ( !indexvalues[indexint] )
                {
                  cdo_warning("Could not parse value '%s' of attribute '%s' to integer.", indexstring, index->c_str());
                  lallindices = false;
                }
              else
                lanyindex = true;
            }
          else
            lallindices = false;
          indexint++;
        }
      if ( lanyindex && vlabel )
        cdo_warning("You specified both variant_label and all indices attributes.\n"
                       "Note that indices attributes have higher priority.");
      else if ( !lallindices && vlabel == nullptr )
        cdo_abort("ERROR (infile: '%s')!\n"
                 "Could not find all index values required for describing the ensemble member (ripf).\n"
                 "Make sure you have specified either 'variant_label' or all 4 'ripf' indexes.",
                  cdo_get_stream_name(0));

      if ( !lanyindex )
        {
          int scanres = std::sscanf((const char *)vlabel, "r%ldi%ldp%ldf%ld",
                               &indexvalues[0], &indexvalues[1], &indexvalues[2], &indexvalues[3]);
          if ( !scanres )
            cdo_abort("ERROR (infile: '%s')!\n"
                 "Could not scan all integers from attribute 'variant_label'. \n"
                     "Make sure it has the format 'rINTiINTpINTfINT'", cdo_get_stream_name(0));
          indexint = 0;
          for (std::vector<std::string>::iterator index=ensindexCMIP6.begin(); index!=ensindexCMIP6.end(); ++index)
            {
              std::string index2string = std::to_string(indexvalues[indexint]);
              kv_insert_vals(kvl, index->c_str(), (char *)index2string.c_str(), true, false);
              indexint++;
            }
        }
      else if ( !lallindices )
        {
          long int vlvalues[4] = { 0 };
          int scanres = std::sscanf((const char *)vlabel, "r%ldi%ldp%ldf%ld",
                               &vlvalues[0], &vlvalues[1], &vlvalues[2], &vlvalues[3]);
          if ( !scanres )
            cdo_warning("Could not scan all integers from attribute 'variant_label'. \n"
                     "Make sure it has the format 'rINTiINTpINTfINT'");
          indexint = 0;
          while ( indexint < 4 )
            {
              if ( !indexvalues[indexint] )
                {
                  indexvalues[indexint] = vlvalues[indexint];
                  if ( !indexvalues[indexint] )
                    cdo_abort("ERROR (infile: '%s')!\n"
                       "You did not provide a value for attribute '%s'.\n"
                       "Make sure you have specified either 'variant_label' or all 4 'ripf' indexes.",
                       cdo_get_stream_name(0), ensindexCMIP6[indexint]);
                  std::string index2string = std::to_string(indexvalues[indexint]);
                  kv_insert_vals(kvl, ensindexCMIP6[indexint].c_str(), (char *)index2string.c_str(), true, false);
                }
              indexint++;
            }
        }
      char member[CMOR_MAX_STRING];
      sprintf(member, "r%ldi%ldp%ldf%ld", indexvalues[0], indexvalues[1], indexvalues[2], indexvalues[3]);
      kv_insert_vals(kvl, "member", member, true, false);
      kv_insert_vals(kvl, "variant_label", member, true, false);
    }
}

static void
read_config_files(KVList *kvl)
{
  if (Options::cdoVerbose) cdo_print("1. Start to read configuration files.");
  /* Files from info key in command line. */
  const KeyValues *info = kvl->search("i");
  int i = 0;
  if (info)
    while (i < info->nvalues)
      {
        if (Options::cdoVerbose) cdo_print("1.1. Try to parse file: '%s' configured with key 'info'.", info->values[i]);
        if (parse_kv_file(kvl, info->values[i].c_str()) == 0) cdo_abort("ERROR (infile: '%s')! File '%s' does not exist.", cdo_get_stream_name(0), info->values[i]);
        if (Options::cdoVerbose) cdo_print("1.1. Successfully parsed file: '%s' configured with key 'info'.", info->values[i]);
        i++;
      }
  else
    if (Options::cdoVerbose) cdo_print("1.1. No info file was passed to the operator.");
  /* Config file in user's $cwd directory. */
  char cwd[1024];
  getcwd(cwd, sizeof(cwd));
  cwd[strlen(cwd)] = '\0';
  const char *dotconfig = ".cdocmorinfo";
  char *workfile = (char *) Malloc((strlen(cwd) + strlen(dotconfig) + 2) * sizeof(char));
  sprintf(workfile, "%s/%s", cwd, dotconfig);
  if (Options::cdoVerbose) cdo_print("1.2. Try to parse default file: '%s'.", workfile);
  if (parse_kv_file(kvl, workfile) == 0 && Options::cdoVerbose)
    cdo_warning("Default file for keyword 'info': '%s' does not exist.", workfile);
  else if (Options::cdoVerbose)
    cdo_print("1.2. Successfully parsed default file: '%s'.", workfile);
  Free(workfile);

  if (i == 0)
    {
      const KeyValues *info2 = kvl->search("i");
      if (info2)
        while (i < info2->nvalues)
          {
            if (Options::cdoVerbose)
              cdo_print("1.3. Try to parse file: '%s' configured with key 'info' in file '.cdocmorinfo'.", info2->values[i]);
            if (parse_kv_file(kvl, info2->values[i].c_str()) == 0)
              cdo_abort("ERROR (infile: '%s')! File '%s' does not exist.", cdo_get_stream_name(0), info2->values[i]);
            if (Options::cdoVerbose)
              cdo_print("1.3. Successfully parsed file: '%s' configured with key 'info' in file '.cdocmorinfo'.", info2->values[i]);
            i++;
          }
    }
  if (Options::cdoVerbose) cdo_print("1. Successfully read configuration files.");
}

static int
in_list(std::vector<std::string> list, const char *needle, int num)
{
  for (int i = 0; i < num; ++i)
    if (list[i] == needle) return 1;
  return 0;
}

static int
get_netcdf_file_action(KVList *kvl, char *proj)
{
  char *chunk = kv_get_a_val(kvl, "om", "a");
  if (chunk[0] == 'r')
    {
#if (CMOR_VERSION_MAJOR == 3)
      return CMOR_REPLACE;
#else
      return CMOR_REPLACE_4;
#endif
    }
  else if (chunk[0] == 'a')
    {
#if (CMOR_VERSION_MAJOR == 3)
      return CMOR_APPEND;
#else
      return CMOR_APPEND_4;
#endif
    }
  else if (chunk[0] == 'p')
    {
#if (CMOR_VERSION_MAJOR == 3)
      return CMOR_PRESERVE;
#else
      return CMOR_PRESERVE_4;
#endif
    }
  else
    {
      cdo_warning("You set output_mode = '%s', but valid are 'a' for append ,'r' for replace or 'p' for preserve.\n "
                     "         CMOR output mode is set to: replace.",
                     chunk);
      return CMOR_REPLACE;
    }
}

static int
get_cmor_verbosity(KVList *kvl)
{
  char *verbos = kv_get_a_val(kvl, "set_verbosity", nullptr);
  if (!verbos) return CMOR_NORMAL;
  if (strcmp(verbos, "CMOR_QUIET") == 0)
    return CMOR_QUIET;
  else
    return CMOR_NORMAL;
}

static int
get_cmor_exit_control(KVList *kvl)
{
  char *exit = kv_get_a_val(kvl, "exit_control", nullptr);
  if (!exit) return CMOR_NORMAL;
  if (strcasecmp(exit, "CMOR_EXIT_ON_MAJOR") == 0)
    return CMOR_EXIT_ON_MAJOR;
  else if (strcasecmp(exit, "CMOR_EXIT_ON_WARNING") == 0)
    return CMOR_EXIT_ON_WARNING;
  else
    return CMOR_NORMAL;
}

static char *
get_calendar_ptr(int calendar)
{
  char *calendar_ptr = (char *) Malloc(CMOR_MAX_STRING * sizeof(char));
  switch (calendar)
    {
    case CALENDAR_STANDARD: strcpy(calendar_ptr, "standard"); break;
    case CALENDAR_GREGORIAN: strcpy(calendar_ptr, "gregorian"); break;
    case CALENDAR_PROLEPTIC: strcpy(calendar_ptr, "proleptic_gregorian"); break;
    case CALENDAR_360DAYS: strcpy(calendar_ptr, "360_day"); break;
    case CALENDAR_365DAYS: strcpy(calendar_ptr, "noleap"); break;
    case CALENDAR_366DAYS: strcpy(calendar_ptr, "all_leap"); break;
    default: Free(calendar_ptr); return nullptr;
    }
  return calendar_ptr;
}

static int
get_calendar_int(char *calendar)
{
  if (!calendar)
    return -1;
  else if (strcmp(calendar, "standard") == 0)
    return CALENDAR_STANDARD;
  else if (strcmp(calendar, "gregorian") == 0)
    return CALENDAR_GREGORIAN;
  else if (strcmp(calendar, "proleptic_gregorian") == 0)
    return CALENDAR_PROLEPTIC;
  else if (strcmp(calendar, "360_day") == 0)
    return CALENDAR_360DAYS;
  else if (strcmp(calendar, "noleap") == 0)
    return CALENDAR_365DAYS;
  else if (strcmp(calendar, "all_leap") == 0)
    return CALENDAR_366DAYS;
  else
    {
      cdo_warning("You set calendar type = '%s' which is not supported by CMOR.", calendar);
      return -1;
    }
}

/***********************************************/
/*Time related functions************************/
/***********************************************/

static char *
get_time_units(int taxisID)
{
  char *units = (char *) Malloc(CMOR_MAX_STRING * sizeof(char));
  int timeunit = taxisInqTunit(taxisID);
  int year, month, day, hour, minute, second, ms;
  const auto rDateTime = taxisInqRdatetime(taxisID);
  cdiDate_decode(rDateTime.date, &year, &month, &day);
  cdiTime_decode(rDateTime.time, &hour, &minute, &second, &ms);
  if (timeunit == TUNIT_QUARTER || timeunit == TUNIT_30MINUTES) timeunit = TUNIT_MINUTE;
  if (timeunit == TUNIT_3HOURS || timeunit == TUNIT_6HOURS || timeunit == TUNIT_12HOURS) timeunit = TUNIT_HOUR;

  sprintf(units, "%s since %d-%d-%d %02d:%02d:%02d", tunitNamePtr(timeunit), year, month, day, hour, minute, second);
  return units;
}

static int
get_time_step_int(char *time_step)
{
  if (strcmp(time_step, "seconds") == 0)
    return TUNIT_SECOND;
  else if (strcmp(time_step, "minutes") == 0)
    return TUNIT_MINUTE;
  else if (strcmp(time_step, "hours") == 0)
    return TUNIT_HOUR;
  else if (strcmp(time_step, "days") == 0)
    return TUNIT_DAY;
  else if (strcmp(time_step, "months") == 0)
    return TUNIT_MONTH;
  else if (strcmp(time_step, "years") == 0)
    return TUNIT_YEAR;
  else
    {
      cdo_warning(
          "You set required_time_units = '%s since...'.\n          This time step is not yet implemented in cmor.",
          time_step);
      return 0;
    }
}

static int
check_time_units(char *time_units)
{
  /* Required attribute in check_att */
  int attyear, attmonth, attday, atthour, attminute, attsecond;
  char time_step[CMOR_MAX_STRING];
  if (sscanf(time_units, "%s since %d-%d-%d%*1s%02d:%02d:%02d%*1s", time_step, &attyear, &attmonth, &attday, &atthour,
             &attminute, &attsecond)
      != 7)
    {
      cdo_warning("You set required_time_units = '%s'\n          but it requires the form 'timestep since "
                 "year-month-day hour:minute:second.\n          Could not read all 7 required time unit values.",
                 time_units);
      return 0;
    }
  if (!get_time_step_int(time_step)) return 0;
  return 1;
}

static void
get_time_method(KVList *kvl, int vlistID, int varID, char *cmor_time_name, char *project_id, int miptab_freq,
                int *time_axis)
{
  if ((strcmp(project_id, "CMIP5") == 0 || strcmp(project_id, "CMIP6") == 0) && miptab_freq) switch (miptab_freq)
      {
      case 1:
        strcpy(cmor_time_name, "time2");
        *time_axis = 2;
        break;
      case 2:
        strcpy(cmor_time_name, "time");
        *time_axis = 0;
        break;
      case 4:
        strcpy(cmor_time_name, "time");
        *time_axis = 0;
        break;
      case 5:
        strcpy(cmor_time_name, "time1");
        *time_axis = 1;
        break;
      case 6:
        strcpy(cmor_time_name, "time1");
        *time_axis = 1;
        break;
      case 7:
        strcpy(cmor_time_name, "time3");
        *time_axis = 3;
        break;
      }
  if (cmor_time_name[0] != 't')
    {
      char *time_method = get_txtatt(vlistID, varID, "cell_methods");
      char *att_time_method = kv_get_a_val(kvl, "cm", nullptr);
      check_compare_set(&time_method, att_time_method, "cell_methods", " ");
      if (time_method[0] == 'm')
        {
          strcpy(cmor_time_name, "time \0");
          *time_axis = 0;
        }
      else if (time_method[0] == 'p')
        {
          strcpy(cmor_time_name, "time1\0");
          *time_axis = 1;
        }
      else if (time_method[0] == 'c')
        {
          strcpy(cmor_time_name, "time2\0");
          *time_axis = 2;
        }
      else if (time_method[0] == 'd')
        {
          strcpy(cmor_time_name, "time3\0");
          *time_axis = 3;
        }
      else if (time_method[0] == 'n')
        {
          strcpy(cmor_time_name, "none\0");
          *time_axis = 4;
        }
      else
        {
          if (time_method[0] == ' ')
            {
              if (Options::cdoVerbose)
                cdo_print("No value found for attribute 'cell_method' of variable with ID '%d'.\n          Depending on "
                         "the aggregation method, you can specifiy one of: \n          'n', 'm', 'p', 'c', 'd'.\n      "
                         "    The default ('m') is used.",
                         varID);
            }
          else
            cdo_warning("The value for attribute cell method = '%s' of variable with ID '%d' is not valid.\n          "
                       "Depending on the aggregation method, you can specifiy one of: \n          'n', 'm', 'p', 'c', "
                       "'d'.\n          The default ('m') is used.",
                       time_method, varID);
          strcpy(cmor_time_name, "time \0");
        }
      Free(time_method);
    }
  if (Options::cdoVerbose) cdo_print("Successfully determined time_axis = '%d'.", *time_axis);
}

static CdiDateTime
get_taxis(char *required_time_units, int *timeunit)
{
  int attyear, attmonth, attday, atthour, attminute, attsecond;
  char atttimeunit[CMOR_MAX_STRING];

  std::sscanf(required_time_units, "%s since %d-%d-%d%*1s%02d:%02d:%02d%*1s", atttimeunit, &attyear, &attmonth, &attday,
         &atthour, &attminute, &attsecond);
  *timeunit = get_time_step_int(atttimeunit);
  CdiDateTime sDateTime{};
  sDateTime.date = cdiDate_encode(attyear, attmonth, attday);
  sDateTime.time = cdiTime_encode(atthour, attminute, attsecond, 0);
  return sDateTime;
}

static double *
get_branch_times(KVList *kvl, int calendar, char *time_units, char *project_id)
{
  if (Options::cdoVerbose) cdo_print("6.1.2. Start to compute attribute 'branch_time'.");
  char *btip = kv_get_a_val(kvl, "branch_time_in_parent", nullptr);
  char *btic = kv_get_a_val(kvl, "branch_time_in_child", nullptr);
  char *ptu = kv_get_a_val(kvl, "parent_time_units", nullptr);
  double *branch_time = (double *) Malloc(2 * sizeof(double));
  branch_time[0] = 0.0;
  branch_time[1] = 0.0;

  if ( strcmp(kv_get_a_val(kvl, "parent_experiment_id", "no parent"), "no parent") != 0 && ( !btip || ( strcmp(project_id, "CMIP6") == 0 && !btic ) ) )
    {

  int numdates = 0;
  std::vector<std::string> branch_dates_p = kv_get_vals(kvl, "branch_dates", &numdates);

  if (numdates == 2 && ptu)
    {
      int branchdates[2], branchyears[2], branchmonths[2], branchdays[2];
      CdiDateTime branchDateTimes[2]{};
      for (int i = 0; i < 2; ++i)
        {
          branchdates[i] = atol(branch_dates_p[i].c_str());
          branchyears[i] = branchdates[i] / 100 / 100;
          branchmonths[i] = (branchdates[i] - branchyears[i] * 100 * 100) / 100;
          branchdays[i] = branchdates[i] - branchyears[i] * 100 * 100 - branchmonths[i] * 100;
          branchDateTimes[i].date = cdiDate_encode(branchyears[i], branchmonths[i], branchdays[i]);
        }

      int parenttimeunit;
      const auto parentsDateTime = get_taxis(ptu, &parenttimeunit);
      const auto parentstartdate = julianDate_encode(calendar, parentsDateTime);
      const auto parentbranchdate = julianDate_encode(calendar, branchDateTimes[0]);

      int childtimeunit;
      const auto childsDateTime = get_taxis(time_units, &childtimeunit);
      const auto childstartdate = julianDate_encode(calendar, childsDateTime);
      const auto childbranchdate = julianDate_encode(calendar, branchDateTimes[1]);

      /* If time unit is always "days since.." */
      branch_time[0] = julianDate_to_seconds(julianDate_sub(parentbranchdate, parentstartdate)) / 86400;
      branch_time[1] = julianDate_to_seconds(julianDate_sub(childbranchdate, childstartdate)) / 86400;
    }
  else
    {
      cdo_warning("Since you specified 'parent_experiment_id' you have to set 'parent_time_units' and either 'branch_dates' or both 'branch_time_in_parent' and 'branch_time_in_child'. Could not find a correct configuration.");
    }
    }
  else
    {
      if ( btip )
        branch_time[0] = atol(btip);
      if ( btic )
        branch_time[1] = atol(btic);
    }
  if (Options::cdoVerbose) cdo_print("6.1.2. Successfully computed 'branch_time_in_parent'='%f' and 'branch_time_in_child'='%f'.", branch_time[0], branch_time[1]);
  return branch_time;
}

static char *
check_required_time_units(KVList *kvl, int taxisID)
{
  if (Options::cdoVerbose) cdo_print("4.1. Start to check attribute 'required_time_units'.");
  char *time_units = get_time_units(taxisID);
  char *required_time_units = kv_get_a_val(kvl, "rtu", nullptr);
  if (required_time_units && check_time_units(required_time_units))
    check_compare_set(&time_units, required_time_units, "time_units", nullptr);
  else
    cdo_warning("Required Attribute 'required_time_units' from configuration is invalid!\n          Continue with infile time units instead.");
  kv_insert_vals(kvl, "rtu", time_units, true, false);
  if (Options::cdoVerbose) cdo_print("4.1. Successfully checked attribute 'required_time_units'.");
  return time_units;
}

static char *
check_calendar(KVList *kvl, int taxisID, int *calendar)
{
  if (Options::cdoVerbose) cdo_print("6.1.1. Start to check attribute 'calendar'.");
  char *attcalendar = kv_get_a_val(kvl, "calendar", nullptr);
  char *calendarptr = get_calendar_ptr(taxisInqCalendar(taxisID));
  if ((*calendar = get_calendar_int(attcalendar)) > -1)
    check_compare_set(&calendarptr, attcalendar, "calendar", nullptr);
  else if ( get_calendar_int(calendarptr) < 0 )
    {
      if ( attcalendar )
        cdo_print("6.1.1. Cannot use calendar '%s' from configuration.", attcalendar);
      if ( calendarptr )
        cdo_print("6.1.1. Cannot use calendar '%s' from infile.", calendarptr);
      if ( Options::cdoVerbose )
        cdo_print("6.1.1. You did not provide a valid calendar. Default 'standard' is used. Valid calendars are:\n          "
                 "'standard', 'gregorian', 'proleptic_gregorian', '360_day', 'noleap' and 'all_leap'.");
      strcpy(calendarptr, "standard");
    }
  else if ( get_calendar_int(calendarptr) > -1 && Options::cdoVerbose )
    cdo_print("6.1.1. Use calendar '%s' from infile.", calendarptr);
  if (Options::cdoVerbose) cdo_print("6.1.1. Successfully retrived calendar: '%s'.", calendarptr);
  return calendarptr;
}

/*********/
/* main: */
/*********/

static bool keep_this_attribute(KeyValues *kv, const char **array)
{
  if (kv->key.c_str() && kv->nvalues == 1)
    {
      int j = 0;
      while (array[j])
        {
          if (kv->key == array[j]) break;
          j++;
        }
      if (!array[j] && strncmp(kv->values[0].c_str(), "notSet", 6) != 0)
        {
          return true;
        }
      else
        return false;
    }
  else
    return false;
}

#if (CMOR_VERSION_MAJOR == 3)

static const char*
copyCV(char *directory)
{
  const char *cvwithout = "CMIP6_CV_without_prefix.json";
  char *cvname = (char *) Malloc((strlen(directory)+strlen(cvwithout)+2)*sizeof(char));
  sprintf(cvname, "%s/%s", directory, cvwithout);
  if ( Options::cdoVerbose ) 
    cdo_print("Check whether a CV without tracking prefix check exists.");
  if ( file_exist((const char *)cvname, false, "CV", false) )
    return cvwithout;
  else
    {
      if ( Options::cdoVerbose ) 
        cdo_print("Try to create a CV without tracking prefix check with program 'sed'.");

      char command[CDI_MAX_NAME];
      sprintf(command, "sed 's/\"hdl:21.14100\\/\\.\\*\"/\"\\^\\.*\"/' %s/CMIP6_CV.json >%s", directory, cvname);
      int dir_err = system(command);
      if (dir_err != 0)
        {
          if ( Options::cdoVerbose ) 
            cdo_print("Creation of a CV without tracking prefix check failed. "
                 "This can be due to missing 'sed' program or missing write permissions. "
                 "Continue with tracking prefix. "
                 "Consider that the output tracking id needs to be registered as a PID." );
          sprintf(cvname, "CMIP6_CV.json");
          return cvname;
        }
      if ( Options::cdoVerbose ) 
        cdo_print("Successfully created a CV without tracking prefix check.");
      return cvwithout;
    }
}
#endif

static void
setup_dataset(KVList *kvl, CdoStreamID streamID, int *calendar, char *project_id)
{
  if (Options::cdoVerbose) cdo_print("6. Start to process cmor_setup and cmor_dataset.");
  int netcdf_file_action = get_netcdf_file_action(kvl, project_id);
  int set_verbosity = get_cmor_verbosity(kvl);
  int exit_control = get_cmor_exit_control(kvl);
  int creat_subs = 1;
  char *drs = kv_get_a_val(kvl, "d", "y");
  if (drs[0] == 'n')
    creat_subs = 0;
  else if (drs[0] != 'y')
    {
      cdo_warning("In preparing cmor_setup:\n          You set 'd' = '%s' which is not valid.\n          Allowed are: "
                 "'n' or 'y'. d is set to 'y'.",
                 drs);
      kv_insert_vals(kvl, "d", (char *) "y", true, false);
    }
  if (strcmp(project_id, "CORDEX") == 0 && creat_subs)
    {
      if (Options::cdoVerbose)
        cdo_print("Since you produce CORDEX compliant output, a path is constructed with DRS elements according to the "
                 "project defined template.");

      char *miptabfreq = kv_get_a_val(kvl, "miptab_freq", nullptr);
      char freq[CDI_MAX_NAME];
      if (strcmp(miptabfreq, "6h") == 0)
        strcpy(freq, "6hr");
      else if (strcmp(miptabfreq, "3h") == 0)
        strcpy(freq, "3hr");
      else if (strcmp(miptabfreq, "1h") == 0)
        strcpy(freq, "1hr");
      else
        strcpy(freq, miptabfreq);

      char cordexDir[CDI_MAX_NAME];
      char cordexFileTem[CDI_MAX_NAME];
      sprintf(cordexDir, "%s/%s/%s/%s/%s/%s/%s/%s/%s/%s/%s", kv_get_a_val(kvl, "dr", "./"), project_id,
              kv_get_a_val(kvl, "product", nullptr), kv_get_a_val(kvl, "CORDEX_domain", nullptr),
              kv_get_a_val(kvl, "institute_id", nullptr), kv_get_a_val(kvl, "driving_model_id", nullptr),
              kv_get_a_val(kvl, "experiment_id", nullptr), kv_get_a_val(kvl, "member", nullptr),
              kv_get_a_val(kvl, "model_id", nullptr), kv_get_a_val(kvl, "rcm_version_id", nullptr), freq);
      sprintf(cordexFileTem, "%s_%s_%s_%s_%s_%s_%s", kv_get_a_val(kvl, "CORDEX_domain", nullptr),
              kv_get_a_val(kvl, "driving_model_id", nullptr), kv_get_a_val(kvl, "experiment_id", nullptr),
              kv_get_a_val(kvl, "member", nullptr), kv_get_a_val(kvl, "model_id", nullptr),
              kv_get_a_val(kvl, "rcm_version_id", nullptr), freq);

      kv_insert_vals(kvl, "dr", (char *) "./", true, false);
      kv_insert_vals(kvl, "cordexDir", cordexDir, true, false);
      kv_insert_vals(kvl, "cordexFileTem", cordexFileTem, true, false);
      creat_subs = 0;
    }

  int vlistID = cdo_stream_inq_vlist(streamID);

  int cmf = cmor_setup(kv_get_a_val(kvl, "inpath", "/usr/share/cmor/"), &netcdf_file_action, &set_verbosity,
                       &exit_control, kv_get_a_val(kvl, "logfile", nullptr), &creat_subs);
  if (cmf != 0) cdo_abort("ERROR (infile: '%s')! Function cmor_setup failed!", cdo_get_stream_name(0));

  signal(SIGTERM, sigfunc);
  int taxisID = vlistInqTaxis(vlistID);

  /*
    char *attcomment = kv_get_a_val(kvl, "comment", nullptr);
    char *comment = get_txtatt(vlistID, CDI_GLOBAL, "comment");
  */

  /* First compare file calendar and config calendar and retrieve pointer and integer
     Then check the required time units from config and retrieve
     Then compute branch_time_in_parent and branch_time_in_child */

  if (Options::cdoVerbose)
    cdo_print("6.1. Start to check model calendar as well as 'required_time_units' and 'branch_times' attributes.");
  char *calendarptr = check_calendar(kvl, taxisID, calendar);
  char *time_units = kv_get_a_val(kvl, "rtu", nullptr);
  double *branch_times = get_branch_times(kvl, *calendar, time_units, project_id);

  if (Options::cdoVerbose) cdo_print("6.1. Successfully found valid calendar, 'required_time_units' and 'branch_times'.");

/* if keep_all_atts is set: */
  const char *kaa_notneeded_general[] = { "cn"           , "n"                , "c"           , "u", "cm",
                                    "vc"           , "p"                , "i"           , "ca", "za",
                                    "gi"           , "rtu"              , "mt"          , "om", "ms",
                                    "dr"           , "d"                , "lc"          , "dj", "workfile4err",
                                    "kaa"          , "mtproof"          , "miptab_freq" , "mip_table_dir",
                                    "grid_info_dir", "mapping_table_dir", "branch_dates"          , "member",
                                    "ta"           , "firsttimeval"     , "calendar"    , "cordexDir",
                                    "cordexFileTem", "branch_time"      , "dl"          , "vd" , 
/*Following attributes are set by CMOR: */
                                    "tracking_id"  , "creation_date"    , "table_id"    , "di",
                                    "ci"           , "sc"               , "ml",
                                    nullptr };
#if defined(CMOR_VERSION_MAJOR)
#if (CMOR_VERSION_MAJOR == 2)
  const char *kaa_cmor2[] = { "parent_experiment", "modeling_realm", nullptr};
  const char *datasetvals[] = { "dr"                  , "experiment_id",
                                "institution"         , "source",
                                "realization"         , "contact",
                                "history"             , "comment",
                                "references"          , "leap_year",
                                "leap_month"          , "model_id",
                                "forcing"             , "initialization_method",
                                "physics_version"     , "institute_id",
                                "parent_experiment_id", "parent_experiment_rip", 
                                "product", nullptr};
  cmf = cmor_dataset(kv_get_a_val(kvl, (char *) datasetvals[0], "./"), kv_get_a_val(kvl, (char *) datasetvals[1], ""),
                       kv_get_a_val(kvl, (char *) datasetvals[2], ""), kv_get_a_val(kvl, (char *) datasetvals[3], ""), calendarptr,
                       atol(kv_get_a_val(kvl, (char *) datasetvals[4], "")), kv_get_a_val(kvl, (char *) datasetvals[5], ""),
                       kv_get_a_val(kvl,(char *) datasetvals[6], ""), kv_get_a_val(kvl, (char *) datasetvals[7], ""),
                       kv_get_a_val(kvl,(char *) datasetvals[8], ""), atol(kv_get_a_val(kvl, (char *) datasetvals[9], "")),
                       atol(kv_get_a_val(kvl, (char *) datasetvals[10], "")), nullptr, kv_get_a_val(kvl, (char *) datasetvals[11], ""),
                       kv_get_a_val(kvl, (char *) datasetvals[12], ""), atol(kv_get_a_val(kvl, (char *) datasetvals[13], "")),
                       atol(kv_get_a_val(kvl, (char *) datasetvals[14], "")), kv_get_a_val(kvl, (char *) datasetvals[15], ""),
                       kv_get_a_val(kvl, (char *) datasetvals[16], ""), &(branch_times[0]),
                       kv_get_a_val(kvl, (char *) datasetvals[17], ""));
  if (cmf != 0) cdo_abort("ERROR (infile: '%s')! Function cmor_dataset failed!", cdo_get_stream_name(0));
  const char *allneeded2[] = { "CORDEX_domain",
                               "driving_experiment",
                               "driving_model_id",
                               "driving_model_ensemble_member",
                               "driving_experiment_name",
                               "rcm_version_id",
                               nullptr };
  int ind = 0;
  if (strcmp(project_id, "CORDEX") == 0)
    while (allneeded2[ind])
      {
        char *tmp = kv_get_a_val(kvl, allneeded2[ind], nullptr);
        if (tmp) cmf = cmor_set_cur_dataset_attribute((char *) allneeded2[ind], tmp, 1);
        if (cmf != 0) cdo_abort("ERROR (infile: '%s')! Function cmor_set_cur_dataset_attribute failed!", cdo_get_stream_name(0));
        ind++;
      }
  if (strcmp(kv_get_a_val(kvl, "kaa", "n"), "y") == 0)
    {
      char notincluded[2048];
      strcpy(notincluded, "The following attributes are not included in the global attributes list.\n          Reasons can be: 1. Attribute is an internal keyword 2. No valaue is available 3. CMOR creates the attribute itself:\n          ");
      size_t inilen = strlen(notincluded); 
      size_t strlens = inilen;
      for (auto &kv : *kvl)
        {
          if ( keep_this_attribute(&kv, kaa_notneeded_general) && keep_this_attribute(&kv, allneeded2) &&
               keep_this_attribute(&kv, datasetvals) && keep_this_attribute(&kv, kaa_cmor2) )
            {
              cmf = cmor_set_cur_dataset_attribute((char *) kv.key.c_str(), (char *)kv.values[0].c_str(), 1);
              if (cmf != 0) cdo_abort("ERROR (infile: '%s')! Function cmor_set_cur_dataset_attribute failed!", cdo_get_stream_name(0));
            }
          else if ( !keep_this_attribute(&kv, kaa_notneeded_general) || !keep_this_attribute(&kv, kaa_cmor2) )
            {
              strlens+=(strlen(kv.key.c_str())+2);
              if ( strlens < 2048 )
                {
                  strcat(notincluded, kv.key.c_str());
                  strcat(notincluded, ", ");
                }
            }
        }
      if ( strlens > inilen )
        if (Options::cdoVerbose) cdo_print("%s", notincluded);
    }
#elif (CMOR_VERSION_MAJOR == 3)
  {
    /***/
    /* Could not give CMOR all attributes separately because some are required to be in a json file (outpath,...). */
    /* Better collect them in this file. */
    /* todo this **/
    /* If a Json file is denoted, read this file and check attributes */
    /***/

    /*
          char *filename = kv_get_a_val(kvl, "dj", nullptr); */

    size_t cwdsize=1024 ;
    char *cwd = (char *) Malloc((cwdsize+1)*sizeof(char));
    getcwd(cwd, cwdsize);
    cwd[strlen(cwd)] = '\0';
    int procID = getpid();
    FILE *dataset_json;
    char *dataset_path = (char *) Malloc((strlen(cwd) + 
                                          strlen("/dataset.json") + 
                                          floor (log10 (abs (procID))) + 1+1) * sizeof(char));

    sprintf(dataset_path, "%s/dataset%d.json", cwd, procID);
    dataset_json = std::fopen(dataset_path, "w+");
    if (!dataset_json)
      cdo_abort("ERROR (infile: '%s')! In preparing cmor_dataset:\n          Could not open a dataset file '%s' for cmor_dataset.",
               cdo_get_stream_name(0), dataset_path);
    fputs("{\n", dataset_json);

    int i = 0;
    if (strcmp(kv_get_a_val(kvl, "kaa", "n"), "y") == 0)
      {
        for (auto &kv : *kvl)
          {
            if ( keep_this_attribute(&kv, kaa_notneeded_general) )
              {
                int linelen = strlen(kv.key.c_str()) + strlen(kv.values[0].c_str()) + 10;
                std::vector<char> line(linelen);
                sprintf(line.data(), "\"%s\" : \"%s\",\n", kv.key.c_str(), kv.values[0].c_str());
                fputs((const char *) line.data(), dataset_json);
              }
          }
      }
    else
      {
        const char *allneeded[] = /*CMIP5*/ {
          "project_id", "experiment_id", "institution", "source", "realization", "contact", "history", "comment",
          "references", "leap_year", "leap_month", "source_id", "model_id", "forcing", "initialization_method",
          "modeling_realm", "physics_version", "institute_id", "parent_experiment_rip",
          /*CORDEX */
          "CORDEX_domain", "driving_experiment", "driving_model_id", "driving_model_ensemble_member",
          "driving_experiment_name", "rcm_version_id",
          /* CMIP6: */
          /* Glob Atts */
          "_cmip6_option", "Conventions", "activity_id", "branch_method", "experiment", "experiment_id",
          "forcing_index", "further_info_url", "grid", "grid_label", "initialization_index", "institution",
          "institution_id", "license", "mip_era", "nominal_resolution", "physics_index", "product", "realization_index",
          "source", "source_id", "source_type", "sub_experiment", "sub_experiment_id", "table_id", "variant_label",
          "parent_experiment_id", "parent_activity_id", "parent_mip_era", "parent_source_id", "parent_variant_label",
          "parent_time_units", "variant_info", "title",
          /* Others */
          "_controlled_vocabulary_file", "_FORMULA_VAR_FILE", "_AXIS_ENTRY_FILE", nullptr
        };
        while (allneeded[i])
          {
            char *tmp = kv_get_a_val(kvl, allneeded[i], "notSet");
            if (strncmp(tmp, "notSet", 6) != 0)
              {
                int linelen = strlen(allneeded[i]) + strlen(tmp) + 10;
                std::vector<char> line(linelen);
                sprintf(line.data(), "\"%s\" : \"%s\",\n", allneeded[i], tmp);
                fputs((const char *) line.data(), dataset_json);
              }
            i++;
          }
      }

    char branch_time_in_parent[CMOR_MAX_STRING];
    char branch_time_in_child[CMOR_MAX_STRING];
    sprintf(branch_time_in_parent, "%.12f", branch_times[0]);
    sprintf(branch_time_in_child, "%.12f", branch_times[1]);

    /* CMOR internal */
    fputs("\"outpath\" : \"", dataset_json);
    fputs(kv_get_a_val(kvl, "dr", "./"), dataset_json);
    fputs("\",\n", dataset_json);
    fputs("\"output_path_template\" : \"", dataset_json);
    if ( strcmp(project_id, "PalMod2") == 0 )
      {
        fputs(kv_get_a_val(kvl, "output_path_template", "<activity_id><institution_id><source_id><experiment_id><"
                                                    "member_id><table><variable_id><grid_label>"),
        dataset_json);
        fputs("\",\n", dataset_json);
        fputs("\"_history_template\" : \"", dataset_json);
        fputs(kv_get_a_val(kvl, "_history_template", 
                           "%s ; CMOR rewrote data to be consistent with <activity_id>, <Conventions> and CF standards."),
             dataset_json);
      }
    else
      fputs(kv_get_a_val(kvl, "output_path_template", "<mip_era><activity_id><institution_id><source_id><experiment_id><"
                                                    "member_id><table><variable_id><grid_label>"),

          dataset_json);
    fputs(kv_get_a_val(kvl, "vd", "<version>"), dataset_json);
    fputs("\",\n", dataset_json);
    fputs("\"output_file_template\" : \"", dataset_json);
    fputs(kv_get_a_val(kvl, "output_file_template",
                       "<variable_id><table><source_id><experiment_id><member_id><grid_label>"),
          dataset_json);
    fputs("\",\n", dataset_json);
    /*          fputs("\"frequency\" : \"", dataset_json);
              fputs(freq, dataset_json);
              fputs("\",\n", dataset_json);  */

    /* cdo cmor preprocessed: */
    fputs("\"calendar\" : \"", dataset_json);
    fputs(calendarptr, dataset_json);
    fputs("\",\n", dataset_json);
    if ( strcmp(kv_get_a_val(kvl, "parent_experiment_id", "no parent"), "no parent") != 0 )
      {
        fputs("\"branch_time_in_parent\" : \"", dataset_json);
        fputs(branch_time_in_parent, dataset_json);
        fputs("\",\n", dataset_json);
        fputs("\"branch_time_in_child\" : \"", dataset_json);
        fputs(branch_time_in_child, dataset_json);
        fputs("\",\n", dataset_json);
      }

    if ( strcmp(kv_get_a_val(kvl, "tp", "y"), "y") == 0 )
      {
        if ( strcmp(project_id, "CMIP6") == 0 )
          {
            fputs("\"tracking_prefix\" : \"hdl:21.14100\"", dataset_json);
            fputs(",\n", dataset_json);
          }
        else if ( strcmp(project_id, "PalMod2") == 0 )
          {
            fputs("\"tracking_prefix\" : \"hdl:21.14105\"", dataset_json);
            fputs(",\n", dataset_json);
          }
        else
          cdo_warning("Don't know the tracking_prefix for project '%s'.\n          "
                     "Continue without setting a tracking prefix.", project_id);
      }
    else
      {
        char *cvname = (char *)copyCV(kv_get_a_val(kvl, "mip_table_dir", nullptr));
        fputs("\"_controlled_vocabulary_file\" : \"", dataset_json);
        fputs(cvname, dataset_json);
        fputs("\",\n", dataset_json);
      }
    if ( strcmp(kv_get_a_val(kvl, "cvn", "y"), "y") == 0 )
      {
        fputs("\"CDO\" : \"", dataset_json);
        fputs(cdo_comment(), dataset_json);
        fputs("\",\n", dataset_json);
      }
    if ( strcmp(kv_get_a_val(kvl, "cdi_grid_type", "n"), "unstructured") == 0 )
      {
        fputs("\"CDI_grid_type\" : \"unstructured\"", dataset_json);
        fputs(",\n", dataset_json);
        fputs("\"grid_type\" : \"unstructured\"", dataset_json);
        fputs(",\n", dataset_json);
      }

    fputs("}\n", dataset_json);
    std::fclose(dataset_json);
    cmf = cmor_dataset_json(dataset_path);
    if (cmf != 0) cdo_abort("ERROR (infile: '%s')! Function cmor_dataset_json failed!", cdo_get_stream_name(0));

    Free(dataset_path);
    removeDataset();
    /*      Free(freq); */
  }

#else
  cdo_abort("ERROR (infile: '%s')! Cmor version %d not yet enabled!", cdo_get_stream_name(0), (int) CMOR_VERSION_MAJOR);
#endif
#else
  cdo_abort("ERROR (infile: '%s')! It is not clear which CMOR version is installed since\n          Makros CMOR_VERSION_MAJOR and "
           "CMOR_VERSION_MINOR are not available.", cdo_get_stream_name(0));
#endif
  Free(calendarptr);
  Free(branch_times);
  if (Options::cdoVerbose) cdo_print("6. Successfully finished cmor_setup and cmor_dataset.");
}

static void
gen_bounds(int n, double *vals, double *bounds)
{
  for (int i = 0; i < n - 1; ++i)
    {
      bounds[2 * i + 1] = 0.5 * (vals[i] + vals[i + 1]);
      bounds[2 * (i + 1)] = 0.5 * (vals[i] + vals[i + 1]);
    }

  bounds[0] = 2 * vals[0] - bounds[1];
  bounds[2 * n - 1] = 2 * vals[n - 1] - bounds[2 * (n - 1)];
}

static bool
get_zcell_bounds(int zaxisID, double *zcell_bounds, double *levels, int zsize)
{
  bool selfGenerated = false;
  double *lbounds;
  lbounds = (double *) Malloc(zsize * sizeof(double));
  zaxisInqLbounds(zaxisID, lbounds);
  double *ubounds;
  ubounds = (double *) Malloc(zsize * sizeof(double));
  zaxisInqUbounds(zaxisID, ubounds);
  if (!lbounds || !ubounds || std::pow((ubounds[1] - ubounds[0]), 2) < 0.001 || std::pow((lbounds[1] - lbounds[0]), 2) < 0.001)
    {
      gen_bounds(zsize, levels, zcell_bounds);
      selfGenerated = true;
    }
  else
    {
      if (lbounds)
        zcell_bounds[0] = lbounds[0];
      else
        zcell_bounds[0] = 0;
      for (int i = 0; i < zsize - 1; ++i)
        {
          zcell_bounds[2 * i + 1] = ubounds[i];
          zcell_bounds[2 * (i + 1)] = lbounds[i + 1];
        }
      if (ubounds)
        zcell_bounds[2 * zsize - 1] = ubounds[zsize - 1];
      else
        zcell_bounds[2 * zsize - 1] = levels[zsize - 1] + (levels[zsize - 1] - zcell_bounds[2 * zsize - 2]);
      selfGenerated = false;
    }
  Free(lbounds);
  Free(ubounds);
  return selfGenerated;
}

static void
get_zhybrid_half(int zaxisID, double *p0, double *alev_val, double *b_val, double *ap_val)
{
  int zsize = zaxisInqSize(zaxisID);
  int vctsize = zaxisInqVctSize(zaxisID);
  if ( vctsize < 1 )
    cdo_abort("ERROR (infile: '%s')! Missing z-axis description. Please provide all parameters"
             " to calculate the formula of a hybrid sigma pressure axis:"
             "\n          p = ap + b *ps", cdo_get_stream_name(0));
  double *vct = (double *) Malloc(vctsize * sizeof(double));
  zaxisInqVct(zaxisID, vct);
  for (int i = 0; i < zsize; ++i)
    {
      ap_val[i] = vct[i];
      b_val[i] = vct[zsize + i];
    }
  for (int i = 0; i < zsize; ++i)
    {
      alev_val[i] = ap_val[i] / p0[0] + b_val[i];
    }
  Free(vct);
}

static void
get_zhybrid(int zaxisID, double *p0, double *alev_val, double *alev_bnds, double *b_val, double *b_bnds, double *ap_val,
            double *ap_bnds)
{
  int zsize = zaxisInqSize(zaxisID);
  int vctsize = zaxisInqVctSize(zaxisID);
  if ( vctsize < 1 )
    cdo_abort("ERROR (infile: '%s')! Missing z-axis description. Please provide all parameters"
             " to calculate the formula of a hybrid sigma pressure axis:"
             "\n          p = ap + b *ps", cdo_get_stream_name(0));
  double *vct = (double *) Malloc(vctsize * sizeof(double));
  zaxisInqVct(zaxisID, vct);
  for (int i = 0; i < (zsize + 1); ++i)
    {
      ap_bnds[i] = vct[i];
      b_bnds[i] = vct[zsize + 1 + i];
    }
  for (int i = 0; i < zsize; ++i)
    {
      ap_val[i] = (ap_bnds[i] + ap_bnds[i + 1]) / 2.0;
      b_val[i] = (b_bnds[i] + b_bnds[i + 1]) / 2.0;
      alev_val[i] = ap_val[i] / p0[0] + b_val[i];
      alev_bnds[i] = ap_bnds[i] / p0[0] + b_bnds[i];
    }
  alev_bnds[zsize] = ap_bnds[zsize] / p0[0] + b_bnds[zsize];
  Free(vct);
}

static size_t
get_strmaxlen(std::vector<std::string> array, size_t len)
{
  size_t result = 0, i;
  for (i = 0; i < len; ++i)
    if (result < strlen(array[i].c_str())) result = strlen(array[i].c_str());
  return result;
}

static void
get_charvals_and_bnds(KVList *kvl, char *chardim, std::vector<std::string> &fvalss, std::vector<std::string> &fbndss, char **funits, int *nofvals, int *nofbnds, char *cmor_name)
{
  bool fivedim = true;
  char *charvalstring = (char *) Malloc(CMOR_MAX_STRING * sizeof(char));
  sprintf(charvalstring, "char_axis_%s_%s", chardim, cmor_name);
  fvalss = kv_get_vals(kvl, charvalstring, nofvals);
  if ( !fvalss.size() )
    {
      if (Options::cdoVerbose)
        cdo_print("Start to register char_axis_%s", chardim);
      sprintf(charvalstring, "char_axis_%s", chardim);
      fvalss = kv_get_vals(kvl, charvalstring, nofvals);
      fivedim = false;
    }
  else
    {
      if (Options::cdoVerbose)
        cdo_print("Start to register char_axis_%s_%s", chardim, cmor_name);
    }
  if ( !fvalss.size() )
    cdo_warning("You specify variables to merge to an axis and the axis name but its values are missing!"
             "\n          Specify its values via attribute 'char_axis_$name' in infofile.");


  if ( fivedim )
    sprintf(charvalstring, "char_axis_%s_%s_bounds", chardim, cmor_name);
  else
    sprintf(charvalstring, "char_axis_%s_bounds", chardim);
  fbndss = kv_get_vals(kvl, charvalstring, nofbnds);
  if ( !fbndss.size() )
    if (Options::cdoVerbose)
      cdo_print("You specified variables to merge to an axis, the axis name and its values."
             "\n          Note that sometimes bounds are required for these axes."
             "\n          Specify its bounds via attribute 'char_axis_$name_bounds' in infofile.");

  if ( fivedim )
    sprintf(charvalstring, "char_axis_%s_%s_units", chardim, cmor_name);
  else
    sprintf(charvalstring, "char_axis_%s_units", chardim);
  *funits = kv_get_a_val(kvl, charvalstring, "");
  if ( !*funits )
    if (Options::cdoVerbose)
      cdo_print("You specified variables to merge to an axis, the axis name and its values."
             "\n          Note that units are required if the axis has digital values."
             "\n          Specify its units via attribute 'char_axis_$name_units' in infofile.");

  Free(charvalstring);
}

static void
register_char_axis(int numchar, std::vector<std::string> charvals, int *axis_ids, char *chardim)
{
  if (Options::cdoVerbose) cdo_print("Start to register a character axis.");
  size_t maxlen = get_strmaxlen(charvals, numchar);
  void *charcmor = (void *) Malloc((numchar * maxlen + 1)* sizeof(char));
  sprintf((char *) charcmor, "%.*s", (int) strlen(charvals[0].c_str()), charvals[0].c_str());
  std::vector<char> blanks(maxlen);
  for (size_t i = 0; i < maxlen; ++i) blanks[i] = ' ';
  char tempb[CMOR_MAX_STRING];  
  sprintf(tempb, "%.*s", (int) (maxlen - strlen(charvals[0].c_str())), blanks.data());
  strcat((char *) charcmor,tempb);
  for (int i = 1; i < numchar; ++i)
    {
      strcat((char *)charcmor,charvals[i].c_str());
      char tempblanks[CMOR_MAX_STRING];
      sprintf(tempblanks,"%.*s",(int)(maxlen - strlen(charvals[i].c_str())),blanks.data());
      strcat((char *)charcmor,tempblanks);
    }
  int cmf = cmor_axis(new_axis_id(axis_ids), chardim, (char *) "", numchar, (void *) charcmor, 'c', NULL, maxlen, NULL);
  if (cmf != 0) cdo_abort("ERROR (infile: '%s')! Function cmor_axis failed!", cdo_get_stream_name(0));
  Free(charcmor);
  if (Options::cdoVerbose) cdo_print("Successfully registered a character axis.");
}


static void
register_fourth_axis(KVList *kvl, int vlistID, int varID, char *varname, int *axis_ids, char *project_id, int miptab_freq, int *mergeIDs)
{
  char *fa = get_txtatt(vlistID, varID, "merge_axis");
  char *chardimatt = kv_get_a_val(kvl, "ca", NULL);
  char *chardim = get_txtatt(vlistID, varID, "character_axis");
  check_compare_set(&chardim, chardimatt, "character_axis", "notSet");
  if ( strcmp(chardim, "notSet") == 0 && fa )
    cdo_abort("ERROR (infile: '%s')! You specify variables to merge to an axis but the axis name and its values are missing!"
             "\n          Specify its name via variable attribute 'character_axis' and its values via 'char_axis_$name' in infofile.", cdo_get_stream_name(0));
  if ( strcmp(chardim, "notSet") != 0 && strcmp(kv_get_a_val(kvl, "capr", "n"), "n") == 0 )
    {
      int nofvals = 0, nofbnds = 0;
      std::vector<std::string> fvalss, fbndss;
      char *funits;
      get_charvals_and_bnds(kvl, chardim, fvalss, fbndss, &funits, &nofvals, &nofbnds, varname);
      if ( nofvals )
        {
          if ( fa )
            {
              if (Options::cdoVerbose) cdo_print("Try to merge variables to a fourth axis.");
              int nvalues = 0;
              char **idss = parse_string_to_values(kv_get_a_val(kvl, "workfile4err", nullptr), fa, &nvalues, "fourth");

              if ( nvalues > 150 )
                cdo_abort("ERROR (infile: '%s')! Only supported for a maximum of 150 variables.", cdo_get_stream_name(0));
              for ( int i = 0; i < nvalues; ++i)
                {
                  if ( std::sscanf( idss[i], "%i", &mergeIDs[i]) == EOF )
                    cdo_abort("ERROR (infile: '%s')! Internal error.", cdo_get_stream_name(0));
                }  
              mergeIDs[nvalues] = CMOR_UNDEFID;
              if ( nofvals != nvalues )
                cdo_abort("ERROR (infile: '%s')! Number of variables '%d' to merge to a fourth axis '%s' is not equal to specified values for this axis: '%d'.", cdo_get_stream_name(0), nvalues, chardim, nofvals);
              free_array(idss);
              Free(fa);
            }
          else if ( nofvals == zaxisInqSize(vlistInqVarZaxis(vlistID, varID)) )
            {
              if (Options::cdoVerbose) cdo_print("Try to replace the vertical axis with a character_axis.");
            }
          else if ( nofvals == 1 )
            {
              if (Options::cdoVerbose) cdo_print("Try to register a forth axis with one value.");
            }
          else
            {
              cdo_abort("ERROR (infile: '%s')! No data to substitute found in infile for specified number of values '%d' for axis '%s'.", cdo_get_stream_name(0), nofvals, chardim);
            }
          double *fbnds = NULL;
          bool daxis = false;
          int i = 0;
          if ( fbndss.size() )
            {
              fbnds = (double *) Malloc(nofbnds * sizeof(double));
              for ( i = 0; i < nofbnds; ++i)
                {
                  int scanreturn = std::sscanf( fbndss[i].c_str(), "%lf", &fbnds[i]);
                  if ( scanreturn == EOF || scanreturn == 0 )
                    {
                      cdo_warning("Could not convert the '%d'th bnds value of the fourth spatial axis which is to merge into a double format.",i);
                      break;
                    }
                }     
            } 
          if ( i == nofbnds && i > 0 )
            daxis = true;
          double *fvals = NULL;
          i = 0;
          if ( nofvals > 0 )
            {
              fvals = (double *) Malloc(nofvals * sizeof(double));
              for ( i = 0; i < nofvals; ++i)
                {
                  int scanreturn = std::sscanf( fvalss[i].c_str(), "%lf", &fvals[i]);
                  if ( scanreturn == EOF || scanreturn == 0 )
                    {
                      if (Options::cdoVerbose) cdo_print("Could not convert the '%d'th value of the fourth spatial axis which is to merge into a double format.", i);
                      break;
                    }
                }
            }
          if ( i != nofvals && daxis )
            cdo_abort("ERROR (infile: '%s')! Could convert bnds of '%s' to double values but failed to convert values of '%s'."
                     "Check axis attribute 'char_axis_%s'.", cdo_get_stream_name(0), chardim, chardim, chardim);
          if ( i == nofvals )
            {  
              if ( cmor_axis(new_axis_id(axis_ids), (char *) chardim, (char *) funits, nofvals,
                              (void *) fvals, 'd', fbnds, 2, NULL) != 0 )
                cdo_abort("ERROR (infile: '%s')! Could not register axis '%s' with double values.", cdo_get_stream_name(0), chardim);
            }
          else
            {
              register_char_axis(nofvals, fvalss, axis_ids, chardim);
            }
          if ( fvals ) Free(fvals);
          if ( fbnds ) Free(fbnds);
          Free(chardim);
        }
    }
  else
    if (Options::cdoVerbose) cdo_print("Fourth axis is '%s'.", chardim);
}

static int
getRegisteredPsid(struct mapping vars[], int ps_index)
{
  int psID = CDI_UNDEFID;
  int i = 0;
  for ( i = 0; vars[i].cdi_varID != CDI_UNDEFID; ++i)
    {
       if ( vars[i].cdi_varID == ps_index )
         psID = i;
    }
  return psID;
}

static int
registerPsid(struct mapping vars[], int psindex, int vlistID)
{
  struct mapping *var = new_var_mapping(vars);
  size_t gridsize = vlistGridsizeMax(vlistID);
  var->cdi_varID = psindex;
  var->cmor_varID = -5;
  cdiDefKeyString(vlistID, psindex, CDI_KEY_NAME, parameter_to_word("ps"));
  if (vlistInqVarDatatype(vlistID, psindex) == CDI_DATATYPE_FLT64)
     {
       var->datatype = 'd';
       var->data = Malloc(gridsize * sizeof(double));
     }
  else
     {
       var->datatype = 'f';
       var->data = Malloc(gridsize * sizeof(float));
     }
  int psID = getRegisteredPsid(vars, psindex);
  if (Options::cdoVerbose) cdo_print("9. Successfully registered surface pressure '%d'.", psID);

  return psID;
}

static void
register_z_axis(KVList *kvl, int vlistID, int varID, int zaxisID, char *varname, int *axis_ids, int *zfactor_id,
                char *project_id, int miptab_freq, int *mergeIDs, struct mapping vars[], CdoStreamID streamID)
{
  char zaxisunits[CDI_MAX_NAME];
  int length = CDI_MAX_NAME;
  cdiInqKeyString(zaxisID, CDI_GLOBAL, CDI_KEY_UNITS, zaxisunits, &length);
  *zfactor_id = 0;
  int cmf = 0;
  int zsize = zaxisInqSize(zaxisID);
  double *levels;

  char *chardimatt = kv_get_a_val(kvl, "ca", nullptr);
  char *chardim = get_txtatt(vlistID, varID, "character_axis");
  check_compare_set(&chardim, chardimatt, "character_axis", "notSet");
  char *capr = kv_get_a_val(kvl, "capr", "n");
  if ( strcmp(capr, "n") != 0 )
    strcpy(chardim, "notSet");

  char *zaxis = get_txtatt(vlistID, varID, "z_axis");
  char *attzaxis = kv_get_a_val(kvl, "za", nullptr);
  check_compare_set(&zaxis, attzaxis, "z_axis", "notSet");

/*
  if (strcmp(chardim, "vegtype") == 0 || strcmp(chardim, "landUse") == 0 || strcmp(chardim, "soilpools") == 0 ) */

  if (zsize > 1)
    {
      levels = (double *) Malloc(zsize * sizeof(double));
      zaxisInqLevels(zaxisID, levels);
      double *zcell_bounds;
      zcell_bounds = (double *) Malloc(2 * zsize * sizeof(double));
      bool selfGenerated = get_zcell_bounds(zaxisID, zcell_bounds, levels, zsize);
      char *zaxisname = (char *) Malloc(CDI_MAX_NAME * sizeof(char));
      length = CDI_MAX_NAME;
      cdiInqKeyString(zaxisID, CDI_GLOBAL, CDI_KEY_NAME, zaxisname, &length);
      if (zaxisInqType(zaxisID) == ZAXIS_PRESSURE)
        {
          if (strcmp(zaxisunits, "") == 0 || zaxisunits[0] == ' ') strcpy(zaxisunits, "Pa");
          if (strcmp(zaxis, "alevel") == 0 )
            cmf = cmor_axis(new_axis_id(axis_ids), (char *) "plevs", (char *) zaxisunits, zsize, (void *) levels, 'd', nullptr, 0,
                            nullptr);
          else if (strcmp(zaxis, "notSet") != 0)
            {
              if ( strncmp(zaxis, "plev7", 5) == 0 )
                cmf = cmor_axis(new_axis_id(axis_ids), zaxis, (char *) zaxisunits, zsize, (void *) levels, 'd', zcell_bounds, 2,
                            nullptr);
              else
                cmf = cmor_axis(new_axis_id(axis_ids), zaxis, (char *) zaxisunits, zsize, (void *) levels, 'd', nullptr, 0,
                            nullptr);
            }
          else if (strcmp(project_id, "CMIP5") != 0 && strcmp(project_id, "CMIP6") != 0)
            cmf = cmor_axis(new_axis_id(axis_ids), (char *) "plevs", (char *) zaxisunits, zsize, (void *) levels,
                            'd', nullptr, 0, nullptr);
          else
            {
              if (strcmp(project_id, "CMIP6") == 0)
                {
                  switch (miptab_freq)
                    {
                    case 4:
                      cmf = cmor_axis(new_axis_id(axis_ids), (char *) "plev8", (char *) zaxisunits, zsize,
                                      (void *) levels, 'd', nullptr, 0, nullptr);
                      break;
                    default:
                      cdo_warning("CMIP6 requires a zaxis name for zaxis type PRESSURE.\n          The operator "
                                 "tries to use a zaxis name matching the number of levels found in infile.");
                      char zaxisname2[CDI_MAX_NAME];
                      std::snprintf(zaxisname2, CDI_MAX_NAME, "plev%d", zsize);
                      cmf = cmor_axis(new_axis_id(axis_ids), zaxisname2, (char *) zaxisunits, zsize, (void *) levels,
                                      'd', nullptr, 0, nullptr);
                      break;
                    }
                }
              else
                {
                  switch (miptab_freq)
                    {
                    case 3:
                      cmf = cmor_axis(new_axis_id(axis_ids), (char *) "plev7", (char *) zaxisunits, zsize,
                                      (void *) levels, 'd', nullptr, 0, nullptr);
                      break;
                    case 4:
                      cmf = cmor_axis(new_axis_id(axis_ids), (char *) "plev8", (char *) zaxisunits, zsize,
                                      (void *) levels, 'd', nullptr, 0, nullptr);
                      break;
                    case 5:
                      cmf = cmor_axis(new_axis_id(axis_ids), (char *) "plev3", (char *) zaxisunits, zsize,
                                      (void *) levels, 'd', nullptr, 0, nullptr);
                      break;
                    default:
                      cmf = cmor_axis(new_axis_id(axis_ids), (char *) "plevs", (char *) zaxisunits, zsize,
                                      (void *) levels, 'd', nullptr, 0, nullptr);
                      break;
                    }
                }
            }
        }
          else if (zaxisInqType(zaxisID) == ZAXIS_HYBRID || zaxisInqType(zaxisID) == ZAXIS_HYBRID_HALF)
            {
              char *time_method = get_txtatt(vlistID, varID, "cell_methods");
              char *att_time_method = kv_get_a_val(kvl, "cm", nullptr);
              check_compare_set(&time_method, att_time_method, "cell_methods", " ");
              int zaxistype = zaxisInqType(zaxisID);
              int vctsize = zaxisInqVctSize(zaxisID);
              if ( 2 * zsize == vctsize )
                zaxistype = ZAXIS_HYBRID_HALF;
              double *alev_val, *alev_bnds = nullptr, *ap_val, *ap_bnds = nullptr, *b_val, *b_bnds = nullptr;
              double *p0 = (double *) Malloc(sizeof(double));
              p0[0] = 101325.0;

              if ( zaxistype == ZAXIS_HYBRID )
                {
                  alev_val = (double *) Malloc(zsize * sizeof(double));
                  alev_bnds = (double *) Malloc((zsize + 1) * sizeof(double));
                  ap_val = (double *) Malloc(zsize * sizeof(double));
                  ap_bnds = (double *) Malloc((zsize + 1) * sizeof(double));
                  b_val = (double *) Malloc(zsize * sizeof(double));
                  b_bnds = (double *) Malloc((zsize + 1) * sizeof(double));
                }
              else
                {
                  alev_val = (double *) Malloc(zsize * sizeof(double));
                  ap_val = (double *) Malloc(zsize * sizeof(double));
                  b_val = (double *) Malloc(zsize * sizeof(double));
                }

              char *mtproof = kv_get_a_val(kvl, "mtproof", nullptr);
              if (mtproof)
                {
                  if (Options::cdoVerbose) cdo_print("Mapping table: '%s' is applied for ps. ", mtproof);
/*                  kv_insert_vals(kvl, "capr", (char *)"y", 1); */
                  int filetype = cdo_inq_filetype(streamID);
                  kv_insert_vals(kvl, "workfile4err", mtproof, true, false);
                  PMList pml = cdo_parse_cmor_file(mtproof, true);
                  std::vector<std::string> tempo = { "ps" };
                  maptab_via_cn(mtproof, pml, tempo, vlistID, 1,
                                kv_get_a_val(kvl, "miptab_freq", nullptr), filetype, nullptr, false);
                }
              int psindex = getVarIDToMap(vlistID, vlistNvars(vlistID), "name", "ps");
              if (psindex == CDI_UNDEFID)
                psindex = getVarIDToMap(vlistID, vlistNvars(vlistID), "code", "134");
              if (psindex == CDI_UNDEFID)
                cdo_abort("ERROR (infile: '%s')! In registration of a vertical axis:\n          Could not find a surface pressure "
                         "variable in infile. Cannot register a hybrid zaxis without surface pressure.", cdo_get_stream_name(0));
              int psID = getRegisteredPsid(vars, psindex);
              if ( psID == CDI_UNDEFID )
                {
                  if (Options::cdoVerbose) cdo_print("Start to register auxiliary ps. ");
                  psID = registerPsid(vars, psindex, vlistID);
                  if (Options::cdoVerbose) cdo_print("Successfully registered auxiliary ps. ");
                }

              if ( zaxistype == ZAXIS_HYBRID )
                {
                  get_zhybrid(zaxisID, p0, alev_val, alev_bnds, b_val, b_bnds, ap_val, ap_bnds);
                  cmf = cmor_axis(new_axis_id(axis_ids), (char *) "alternate_hybrid_sigma", (char *) "", zsize,
                              (void *) alev_val, 'd', alev_bnds, 1, nullptr);
                }
              else
                {
                  get_zhybrid_half(zaxisID, p0, alev_val, b_val, ap_val);
                  cmf = cmor_axis(new_axis_id(axis_ids), (char *) "alternate_hybrid_sigma_half", (char *) "", zsize,
                              (void *) alev_val, 'd', nullptr, 1, nullptr);
                }
              /*cmor_zfactor (int *zfactor_id,int zaxis_id, char *zfactor_name, char *units, int ndims, int axis_ids[],
               * char type, void *zfactor_values, void *zfactor_bounds)*/



              int lev_id = axis_ids[count_axis_ids(axis_ids) - 1];
              int lev_id_array[2];
              lev_id_array[0] = lev_id;
              cmf = cmor_zfactor(zfactor_id, lev_id, (char *) "p0", (char *) "Pa", 0, 0, 'd', (void *) p0, nullptr);
              if ( zaxistype == ZAXIS_HYBRID )
                {
                  cmf = cmor_zfactor(zfactor_id, lev_id, (char *) "b", (char *) "", 1, &lev_id_array[0], 'd',
                                 (void *) b_val, (void *) b_bnds);
                  cmf = cmor_zfactor(zfactor_id, lev_id, (char *) "ap", (char *) "Pa", 1, &lev_id_array[0], 'd',
                                 (void *) ap_val, (void *) ap_bnds);
                  if ( time_method[0] == 'p' )
                    {
                      if ( count_axis_ids(axis_ids) == 3 )
                        cmf = cmor_zfactor(zfactor_id, lev_id, (char *) "ps3", (char *) "Pa", count_axis_ids(axis_ids) - 1,
                                 axis_ids, vars[psID].datatype, nullptr, nullptr);
                      else
                        cmf = cmor_zfactor(zfactor_id, lev_id, (char *) "ps1", (char *) "Pa", count_axis_ids(axis_ids) - 1,
                                 axis_ids, vars[psID].datatype, nullptr, nullptr);
                    }
                  else if ( time_method[0] == 'c' )
                    cmf = cmor_zfactor(zfactor_id, lev_id, (char *) "ps2", (char *) "Pa", count_axis_ids(axis_ids) - 1,
                                 axis_ids, vars[psID].datatype, nullptr, nullptr);
                  else
                    cmf = cmor_zfactor(zfactor_id, lev_id, (char *) "ps", (char *) "Pa", count_axis_ids(axis_ids) - 1,
                                 axis_ids, vars[psID].datatype, nullptr, nullptr);
                }
              else
                {
                  cmf = cmor_zfactor(zfactor_id, lev_id, (char *) "b_half", (char *) "", 1, &lev_id_array[0], 'd',
                                 (void *) b_val, nullptr);
                  cmf = cmor_zfactor(zfactor_id, lev_id, (char *) "ap_half", (char *) "Pa", 1, &lev_id_array[0], 'd',
                                 (void *) ap_val, nullptr);
                  if ( time_method[0] == 'p' )
                    {
                      if ( count_axis_ids(axis_ids) == 3 )
                        cmf = cmor_zfactor(zfactor_id, lev_id, (char *) "ps3", (char *) "Pa", count_axis_ids(axis_ids) - 1,
                                 axis_ids, vars[psID].datatype, nullptr, nullptr);
                      else
                        cmf = cmor_zfactor(zfactor_id, lev_id, (char *) "ps1", (char *) "Pa", count_axis_ids(axis_ids) - 1,
                                 axis_ids, vars[psID].datatype, nullptr, nullptr);
                    }
                  else if ( time_method[0] == 'c' )
                    cmf = cmor_zfactor(zfactor_id, lev_id, (char *) "ps2", (char *) "Pa", count_axis_ids(axis_ids) - 1,
                                 axis_ids, vars[psID].datatype, nullptr, nullptr);
                  else
                    cmf = cmor_zfactor(zfactor_id, lev_id, (char *) "ps", (char *) "Pa", count_axis_ids(axis_ids) - 1,
                                 axis_ids, vars[psID].datatype, nullptr, nullptr);
                }
              Free(alev_val);
              Free(ap_val);
              Free(b_val);
              if ( zaxistype == ZAXIS_HYBRID )
                {
                  Free(ap_bnds);
                  Free(alev_bnds);
                  Free(b_bnds);
                }
            }
      else if (zaxisInqType(zaxisID) == ZAXIS_DEPTH_BELOW_SEA || zaxisInqType(zaxisID) == ZAXIS_DEPTH_BELOW_LAND)
        {
          if (!zaxisunits[0] || zaxisunits[0] == ' ' || zaxisunits[0] == '\0')
            {
              if ( zaxisInqType(zaxisID) == ZAXIS_DEPTH_BELOW_SEA )
                strcpy(zaxisunits, "m");
              else
                strcpy(zaxisunits, "cm");
            }
          if ( selfGenerated )
            {
              zcell_bounds[0] = (double) 0;
              zcell_bounds[2*zsize-1] = levels[zsize-1];
            }
          char longname[CDI_MAX_NAME];
          length = CDI_MAX_NAME;
          cdiInqKeyString(zaxisID, CDI_GLOBAL, CDI_KEY_LONGNAME, longname, &length);
          if ( strcmp(longname, "depth_below_sea") == 0 || strcmp(zaxis, "olevel") == 0 || strcmp(longname, "ocean depth coordinate") == 0 )
            cmf = cmor_axis(new_axis_id(axis_ids), (char *) "depth_coord", (char *) zaxisunits, zsize,
                          (void *) levels, 'd', zcell_bounds, 2, nullptr);
          else
            cmf = cmor_axis(new_axis_id(axis_ids), (char *) "sdepth", (char *) zaxisunits, zsize, (void *) levels,
                          'd', zcell_bounds, 2, nullptr);
        }
      else if (zaxisInqType(zaxisID) == ZAXIS_ALTITUDE )
        {
          if (strcmp(zaxis, "alevel") == 0 )
              cmf = cmor_axis(new_axis_id(axis_ids), (char *) "alt", (char *) zaxisunits, zsize, (void *) levels,
                          'd', zcell_bounds, 2, nullptr);
          else if ( strncmp(zaxis, "alt", 3) == 0 )
            {
              cmf = cmor_axis(new_axis_id(axis_ids), zaxis, (char *) zaxisunits, zsize, (void *) levels,
                          'd', zcell_bounds, 2, nullptr);
            }
          else
            {
              char zaxisname2[CDI_MAX_NAME];
              std::snprintf(zaxisname2, CDI_MAX_NAME, "alt%d", zsize);
              cmf = cmor_axis(new_axis_id(axis_ids), zaxisname2, (char *) zaxisunits, zsize, (void *) levels,
                          'd', zcell_bounds, 2, nullptr);
            }
        }
      else if (zaxisInqType(zaxisID) == ZAXIS_GENERIC && strcmp(chardim, "notSet") != 0 )
        {
          register_fourth_axis(kvl, vlistID, varID, varname, axis_ids, project_id, miptab_freq, mergeIDs);
          kv_insert_vals(kvl, "capr", (char *)"y", true, false);
        }
      else if (zaxisInqType(zaxisID) == ZAXIS_GENERIC || zaxisInqType(zaxisID) == ZAXIS_HEIGHT)
        {
          if (strcmp(zaxisname, "rho") == 0)
            {
              if (strcmp(zaxisunits, "kg m-3") != 0)
                {
                  cdo_abort("ERROR (infile: '%s')! For zaxis with name 'rho' the units must be kg m-3 but are: '%s'", cdo_get_stream_name(0), zaxisunits);
                }
              else
                {
                  cmf = cmor_axis(new_axis_id(axis_ids), (char *) "rho", (char *) "kg m-3", zsize, (void *) levels,
                                  'd', zcell_bounds, 2, nullptr);
                }
              Free(zaxisunits);
            }
          else
            cdo_abort("ERROR (infile: '%s')! In registration of a vertical axis:\n          Z-axis type %d with name '%s' not yet "
                     "enabled.",
                     cdo_get_stream_name(0), zaxisInqType(zaxisID), zaxisname);
          Free(zaxisname);
        }
      else
        cdo_abort("ERROR (infile: '%s')! In registration of a vertical axis:\n          Invalid Z-axis type %d . ",
                 cdo_get_stream_name(0), zaxisInqType(zaxisID));
      Free(zcell_bounds);
      Free(levels);
    }

  else if (zsize == 1 && strcmp(zaxis, "notSet") != 0)
    {
      /*      if ( strcmp(zaxis, "notSet") == 0 )
              {
                if (
              }
            else
              {*/
      char *szc_value = kv_get_a_val(kvl, (const char *) zaxis, nullptr);
      if (szc_value)
        {
          levels = (double *) Malloc(sizeof(double));
          levels[0] = (double) atof(szc_value);

          char *szc_key = (char *) Malloc(CMOR_MAX_STRING * sizeof(char));
          sprintf(szc_key, "%s_bounds", zaxis);

          int numchar = 0;
          std::vector<std::string> szc_bndss = kv_get_vals(kvl, szc_key, &numchar);
          if ( numchar != 0 && numchar != 2 )
            cdo_warning("Scalar z coordinate bounds need to be exactly two values! You configured '%d'.", numchar);

          double szc_bnds[2];
          if ( !szc_bndss.size() )
            {
              if ( std::sscanf( szc_bndss[0].c_str(), "%lf", &szc_bnds[0]) == EOF )
                cdo_abort("ERROR (infile: '%s')! Internal error.", cdo_get_stream_name(0));
              if ( std::sscanf( szc_bndss[1].c_str(), "%lf", &szc_bnds[1]) == EOF )
                cdo_abort("ERROR (infile: '%s')! Internal error.", cdo_get_stream_name(0));
            }

          sprintf(szc_key, "%s_units", zaxis);
          char *szcunits = kv_get_a_val(kvl, (const char *) szc_key, "m");

          if (Options::cdoVerbose)
            cdo_print("Scalar z coordinate name is: '%s'\n          Scalar z coordinate value is: '%f' meter", zaxis,
                     levels[0]);
          if ( !szc_bndss.size() )
            cmf = cmor_axis(new_axis_id(axis_ids), zaxis, (char *) szcunits, zsize, (void *) levels, 'd', szc_bnds, 2, nullptr);
          else
            cmf = cmor_axis(new_axis_id(axis_ids), zaxis, (char *) szcunits, zsize, (void *) levels, 'd', nullptr, 0, nullptr);
          Free(levels);
        }
      else
        cdo_print("You specified z_axis='%s'.\n          No value has been specified, axis will be created with the "
                 "default value.",
                 zaxis);
    }
  else
    cdo_print("Vertical axis is either default and scalar or not available.");

  if (zaxis) Free(zaxis);
  Free(chardim);
  if (cmf != 0) cdo_abort("ERROR (infile: '%s')! Function cmor_axis failed!", cdo_get_stream_name(0));
}

/*
static void register_character_dimension(int *axis_ids, char *filename)
{
  printf("The grid type is generic and a dimension 'basin' is found.\nTherefore, it is tried to read the character
dimension.\n");
  int nc_file_id, nfiledims, nvars, ngatts, unlimdimid;
  nc_type xtypep;
  int varndims, varnattsp;
  int *vardimids;

  char *varname = Malloc(36 * sizeof(char));
  char *dimname = Malloc(36 * sizeof(char));

  size_t dimlength, dimstrlength;

  nc_open(filename, NC_NOWRITE, &nc_file_id);
  nc_inq(nc_file_id, &nfiledims, &nvars, &ngatts, &unlimdimid);
  vardimids = Malloc(nfiledims * sizeof(int));
  void *final_chardim;
  for ( int i = 0; i < nvars; i++ )
    {
      nc_inq_var(nc_file_id, i, varname, &xtypep, &varndims, vardimids, &varnattsp);
      if ( strcmp(varname, "region") == 0 )
        {
          nc_inq_dim(nc_file_id, vardimids[1], dimname, &dimstrlength);
          nc_inq_dim(nc_file_id, vardimids[0], dimname, &dimlength);

          final_chardim = (void *)Malloc(dimstrlength * dimlength *sizeof(char));
          nc_get_var(nc_file_id, i, final_chardim);
        }
    }
  nc_close(nc_file_id);
  cmor_axis(new_axis_id(axis_ids), dimname, "", dimlength, final_chardim, 'c',  nullptr, dimstrlength, nullptr);
  Free(varname);
  Free(dimname);
  Free(vardimids);
}
*/
static void
change_zaxis(KVList *kvl, int vlistID, int zaxisID, int zaxisID2, char *grid_file)
{
  int a, b;
  a = zaxisInqSize(zaxisID);
  b = zaxisInqSize(zaxisID2);

  bool noZaxis = false;
  if (zaxisInqType(zaxisID2) == ZAXIS_SURFACE && b == 1)
    {
      double level[1];
      zaxisInqLevels(zaxisID2, level);
      if (level[0] < 1) noZaxis = true;
    }

  if (!noZaxis)
    {
      if (a != b)
        {
          cdo_warning("Could not use zaxis from file '%s' configured via attribute 'ginfo'\n          because total "
                 "size of infile: '%d' is not identical to total size of ginfo file: '%d'.",
                 grid_file, a, b);
        }
      else if (zaxisID != zaxisID2)
        {
          vlistChangeZaxis(vlistID, zaxisID, zaxisID2);
          char zaxisname[CDI_MAX_NAME];
          int length = CDI_MAX_NAME;
          cdiInqKeyString(zaxisID, CDI_GLOBAL, CDI_KEY_NAME, zaxisname, &length);
          kv_insert_vals(kvl, "za", zaxisname, false, false);
          cdo_print("Successfully substituted zaxis.");
        }
      else
        cdo_print("Zaxis from grid info file %s is equal to infile Zaxis.", grid_file);
    }
}

static void
change_grid(int vlistID, int gridID, int gridID2, char *grid_file)
{
  if (!gridID2)
    cdo_abort("ERROR (infile: '%s')! Could not use grid from file '%s' configured via attribute 'ginfo'\n          because of internal "
             "problems.",
             cdo_get_stream_name(0), grid_file);

  int a, b;
  bool lswitch = true;
  a = gridInqSize(gridID);
  b = gridInqSize(gridID2);

  if (b == 1)
    cdo_print("Grid is not substituted because the size of the grid found in the grid info file is one.");
  else
    {
      if (a != b)
        {
          lswitch = false;
          cdo_warning("Could not use grid from file '%s' configured via attribute 'ginfo'\n          because total "
                 "size of infile: '%d' is not identical to total size of ginfo file: '%d'.",
                 grid_file, a, b);
        }

      a = gridInqYsize(gridID);
      b = gridInqYsize(gridID2);
      if (a != b)
        {
          if ( gridInqType(gridID) != GRID_UNSTRUCTURED )
            lswitch = false;
          cdo_warning("Could not use grid from file '%s' configured via attribute 'ginfo'\n          because ysize of grid info file '%d' does not match ysize of infile '%d'.", grid_file, b, a);
        }

      a = gridInqXsize(gridID);
      b = gridInqXsize(gridID2);
      if (a != b)
        {
          if ( gridInqType(gridID) != GRID_UNSTRUCTURED )
            lswitch = false;
          lswitch = false;
          cdo_warning("Could not use grid from file '%s' configured via attribute 'ginfo'\n          because xsize of grid info file '%d' does not match xsize of infile '%d'.", grid_file, b, a);
        }

      if ( lswitch && gridID != gridID2)
        {
          vlistChangeGrid(vlistID, gridID, gridID2);
          cdo_print("Successfully substituted grid.");
        }
      else if ( lswitch )
        cdo_print("Grid from grid info file %s is equal to infile grid.", grid_file);
    }
}

static void
move_lons(double *xcoord_vals, double *xcell_bounds, int xsize, int xboundsize, int xnbounds)
{
  int testbool = 0;
  for (int i = 0; i < xsize; ++i)
    if (xcoord_vals[i] < 0.0)
      {
        testbool = 1;
        break;
      }
  if (testbool > 0)
    for (int i = 0; i < xsize; ++i)
      if (xcoord_vals[i] < 0) xcoord_vals[i] += 360.0;
  if (xnbounds > 1 && testbool > 0)
    for (int j = 0; j < xboundsize; ++j)
      if (xcell_bounds[j] < 0) xcell_bounds[j] += 360.0;
}

static void
inquire_vals_and_bounds(int gridID, int *xnbounds, int *ynbounds, double *xcoord_vals, double *ycoord_vals,
                        double *xcell_bounds, double *ycell_bounds)
{
  char unitstring[CDI_MAX_NAME];
  gridInqYvals(gridID, ycoord_vals);
  *ynbounds = gridInqYbounds(gridID, ycell_bounds);
  int length = CDI_MAX_NAME;
  cdiInqKeyString(gridID, CDI_YAXIS, CDI_KEY_UNITS, unitstring, &length);
  if ( strcmp(unitstring,"radian") == 0 )
    {
      for ( size_t i = 0; i < gridInqXsize(gridID); i++ )
        ycoord_vals[i] = (180. / M_PI) * ycoord_vals[i];
      for ( size_t i = 0; i < gridInqNvertex(gridID)*gridInqSize(gridID); i++ )
        ycell_bounds[i] = (180. / M_PI) * ycell_bounds[i];
    }
  gridInqXvals(gridID, xcoord_vals);
  *xnbounds = gridInqXbounds(gridID, xcell_bounds);
  length = CDI_MAX_NAME;
  cdiInqKeyString(gridID, CDI_XAXIS, CDI_KEY_UNITS, unitstring, &length);
  if ( strcmp(unitstring,"radian") == 0 )
    {
      for ( size_t i = 0; i < gridInqXsize(gridID); i++ )
        xcoord_vals[i] = (180. / M_PI) * xcoord_vals[i];
      for ( size_t i = 0; i < gridInqNvertex(gridID)*gridInqSize(gridID); i++ )
        xcell_bounds[i] = (180. / M_PI) * xcell_bounds[i];
    }

}

static void
get_cmor_table(KVList *kvl, char *project_id)
{
  int gridtable_id;
  int cmf = 0;
  char gridtable[CMOR_MAX_STRING];
  char *mip_table_dir = kv_get_a_val(kvl, "mip_table_dir", nullptr);
#if (CMOR_VERSION_MAJOR == 2)
  if (mip_table_dir[strlen(mip_table_dir) - 1] == '/')
    sprintf(gridtable, "%s%s_grids", mip_table_dir, project_id);
  else
    sprintf(gridtable, "%s/%s_grids", mip_table_dir, project_id);
#elif (CMOR_VERSION_MAJOR == 3)
  sprintf(gridtable, "%s/%s_grids.json", mip_table_dir, project_id);
#endif
  if (file_exist(gridtable, false, "Cmor-grid_table", false))
    {
      cmf = cmor_load_table(gridtable, &gridtable_id);
      cmf = cmor_set_table(gridtable_id);
    }
  else
    cdo_abort("ERROR (infile: '%s')! In grid registration:\n          File '%s' not found!\n          "
             "A project grid table is required for this type of grid but not "
             "found in the mip table directory '%s'.",
             cdo_get_stream_name(0), gridtable, mip_table_dir);
  if (cmf != 0) cdo_abort("ERROR (infile: '%s')! Function cmor_load_table or cmor_set_table failed!", cdo_get_stream_name(0));
}

static void
check_and_gen_bounds(int gridID, int nbounds, int length, double *coord_vals, double *cell_bounds, int x)
{
  if (nbounds != 2 * length)
    {
      gen_bounds(length, coord_vals, cell_bounds);
      if (x)
        {
          gridDefNvertex(gridID, 2);
          gridDefXbounds(gridID, cell_bounds);
        }
      else
        gridDefYbounds(gridID, cell_bounds);
    }
}

static double
lonbnds_mids_trans_check(double value1, double value2)
{
  if (std::fabs(value1 - value2) < 180.0)
    return (value1 + value2) * 0.5;
  else
    {
      if (value1 + value2 < 360.0)
        return (value1 + value2 + 360.0) * 0.5;
      else
        return (value1 + value2 + 360.0) * 0.5 - 360.0;
    }
}

static double
lonbnds_bnds_trans_check(double value1, double value2)
{
  if (std::fabs(value1 - value2) < 180)
    {
      if (2 * value1 < value2)
        return (2 * value1 - value2 + 360.0);
      else if (2 * value1 > value2 + 360.0)
        return (2 * value1 - value2 - 360.0);
      else
        return (2 * value1 - value2);
    }
  else if (value1 - value2 > 180)
    return (2 * value1 - value2 - 360.0);
  else
    return (2 * value1 - value2 + 360.0);
}

static void
check_and_gen_bounds_curv(int gridID, int totalsize, int xnbounds, int xlength, double *xcoord_vals,
                          double *xcell_bounds, int ynbounds, int ylength, double *ycoord_vals, double *ycell_bounds)
{
  if (xnbounds != 4 * totalsize || ynbounds != 4 * totalsize
      || (IS_EQUAL(xcell_bounds[1], 0.0) && IS_EQUAL(xcell_bounds[2], 0.0))
      || (IS_EQUAL(ycell_bounds[1], 0.0) && IS_EQUAL(ycell_bounds[2], 0.0)))
    {
      Varray2D<double> halflons(xlength + 1, Varray<double>(ylength));
      Varray2D<double> halflats(xlength, Varray<double>(ylength + 1));
      Varray2D<double> halflonsOnhalflats(xlength + 1, Varray<double>(ylength + 1));
      Varray2D<double> halflatsOnhalflons(xlength + 1, Varray<double>(ylength + 1));

      /**/
      /*************Half-lons with 360-0 transmission check**************/
      /**/
      for (int j = 0; j < ylength; ++j)
        {
          for (int i = 1; i < xlength; ++i)
            halflons[i][j] = lonbnds_mids_trans_check(xcoord_vals[i - 1 + j * xlength], xcoord_vals[i + j * xlength]);
          /*left and right boundary: */
          halflons[0][j] = lonbnds_bnds_trans_check(xcoord_vals[j * xlength], halflons[1][j]);
          halflons[xlength][j] = lonbnds_bnds_trans_check(xcoord_vals[j * xlength - 1], halflons[xlength - 1][j]);
        }
      /**/
      /*************Half-lats **************/
      /**/
      for (int i = 0; i < xlength; ++i)
        {
          for (int j = 1; j < ylength; ++j)
            halflats[i][j] = (ycoord_vals[i + (j - 1) * xlength] + ycoord_vals[i + j * xlength]) * 0.5;
          /*upper and lower boundary: */
          halflats[i][0] = 2 * ycoord_vals[i] - halflats[i][1];
          halflats[i][ylength] = 2 * ycoord_vals[i + (ylength - 1) * xlength] - halflats[i][ylength - 1];
        }
      /**/
      /****************Half-lons-on-half-lats with 0-360 transmission check**********/
      /****************Half-lats-on-half-lons                              **********/
      /**/

      for (int i = 1; i < xlength; ++i)
        {
          for (int j = 1; j < ylength; ++j)
            {
              halflonsOnhalflats[i][j] = lonbnds_mids_trans_check(halflons[i][j - 1], halflons[i][j]);
              halflatsOnhalflons[i][j] = (halflats[i - 1][j] + halflats[i][j]) * 0.5;
            }
          /*upper and lower boundary: */
          halflonsOnhalflats[i][0] = lonbnds_bnds_trans_check(halflons[i][0], halflonsOnhalflats[i][1]);
          halflonsOnhalflats[i][ylength]
              = lonbnds_bnds_trans_check(halflons[i][ylength - 1], halflonsOnhalflats[i][ylength - 1]);
          halflatsOnhalflons[i][0] = (halflats[i - 1][0] + halflats[i][0]) * 0.5;
          halflatsOnhalflons[i][ylength] = (halflats[i - 1][ylength] + halflats[i][ylength]) * 0.5;
        }

      /*left and right boundary: */
      for (int j = 1; j < ylength; ++j)
        {
          halflonsOnhalflats[0][j] = lonbnds_mids_trans_check(halflons[0][j - 1], halflons[0][j]);
          halflonsOnhalflats[xlength][j] = lonbnds_mids_trans_check(halflons[xlength][j - 1], halflons[xlength][j]);

          halflatsOnhalflons[0][j] = 2 * halflats[0][j] - halflatsOnhalflons[1][j];
          halflatsOnhalflons[xlength][j] = 2 * halflats[xlength - 1][j] - halflatsOnhalflons[xlength - 1][j];
        }
      halflatsOnhalflons[0][0] = 2 * halflats[0][0] - halflatsOnhalflons[1][0];
      halflatsOnhalflons[0][ylength] = 2 * halflats[0][ylength] - halflatsOnhalflons[1][ylength];
      halflatsOnhalflons[xlength][0] = 2 * halflats[xlength - 1][0] - halflatsOnhalflons[xlength - 1][0];
      halflatsOnhalflons[xlength][ylength]
          = 2 * halflats[xlength - 1][ylength] - halflatsOnhalflons[xlength - 1][ylength];

      halflonsOnhalflats[0][0] = lonbnds_bnds_trans_check(halflons[0][0], halflonsOnhalflats[0][1]);
      halflonsOnhalflats[0][ylength]
          = lonbnds_bnds_trans_check(halflons[0][ylength - 1], halflonsOnhalflats[0][ylength - 1]);
      halflonsOnhalflats[xlength][0] = lonbnds_bnds_trans_check(halflons[xlength][0], halflonsOnhalflats[xlength][1]);
      halflonsOnhalflats[xlength][ylength]
          = lonbnds_bnds_trans_check(halflons[xlength][ylength - 1], halflonsOnhalflats[xlength - 1][ylength]);

      for (int i = 0; i < xlength; ++i)
        for (int j = 0; j < ylength; ++j)
          {
            xcell_bounds[4 * (j * xlength + i)] = halflonsOnhalflats[i][j + 1];
            xcell_bounds[4 * (j * xlength + i) + 1] = halflonsOnhalflats[i][j];
            xcell_bounds[4 * (j * xlength + i) + 2] = halflonsOnhalflats[i + 1][j];
            xcell_bounds[4 * (j * xlength + i) + 3] = halflonsOnhalflats[i + 1][j + 1];
            ycell_bounds[4 * (j * xlength + i)] = halflatsOnhalflons[i][j + 1];
            ycell_bounds[4 * (j * xlength + i) + 1] = halflatsOnhalflons[i][j];
            ycell_bounds[4 * (j * xlength + i) + 2] = halflatsOnhalflons[i + 1][j];
            ycell_bounds[4 * (j * xlength + i) + 3] = halflatsOnhalflons[i + 1][j + 1];
          }
      gridDefNvertex(gridID, 4);
      gridDefXbounds(gridID, xcell_bounds);
      gridDefYbounds(gridID, ycell_bounds);
    }
}
static void
register_lon_axis(int gridID, int xlength, int *axis_ids)
{
  Varray<double> xcoord_vals(xlength);
  if (gridInqXvals(gridID, xcoord_vals.data()))
    {
      Varray<double> xcell_bounds(2 * xlength);
      int xnbounds = gridInqXbounds(gridID, xcell_bounds.data());
      check_and_gen_bounds(gridID, xnbounds, xlength, xcoord_vals.data(), xcell_bounds.data(), 1);
      int cmf = cmor_axis(new_axis_id(axis_ids), (char *) "longitude", (char *) "degrees_east", xlength,
                          (void *) xcoord_vals.data(), 'd', (void *) xcell_bounds.data(), 2, nullptr);
      if (cmf != 0) cdo_abort("ERROR (infile: '%s')! Function cmor_axis failed!", cdo_get_stream_name(0));
    }
}

static void
register_lat_axis(int gridID, int ylength, int *axis_ids)
{
  Varray<double> ycoord_vals(ylength);
  if (gridInqYvals(gridID, ycoord_vals.data()))
    {
      Varray<double> ycell_bounds(2 * ylength);
      int ynbounds = gridInqYbounds(gridID, ycell_bounds.data());
      check_and_gen_bounds(gridID, ynbounds, ylength, ycoord_vals.data(), ycell_bounds.data(), 0);
      int cmf = cmor_axis(new_axis_id(axis_ids), (char *) "latitude", (char *) "degrees_north", ylength,
                          (void *) ycoord_vals.data(), 'd', (void *) ycell_bounds.data(), 2, nullptr);
      if (cmf != 0) cdo_abort("ERROR (infile: '%s')! Function cmor_axis failed!", cdo_get_stream_name(0));
    }
}

static void
register_projection(int *grid_ids, int projID, double *ycoord_vals, double *xcoord_vals, double *ycell_bounds,
                    double *xcell_bounds, int xlength, int ylength, int projtype)
{
  int cmf = 0;
  int pxnbounds;
  int pynbounds;
  int pylength = gridInqYsize(projID);
  int pxlength = gridInqXsize(projID);
  double *pxcoord_vals = (double *) Malloc(pxlength * sizeof(double));
  double *pycoord_vals = (double *) Malloc(pylength * sizeof(double));
  double *pxcell_bounds = (double *) Malloc(2 * pxlength * sizeof(double));
  double *pycell_bounds = (double *) Malloc(2 * pylength * sizeof(double));
  inquire_vals_and_bounds(projID, &pxnbounds, &pynbounds, pxcoord_vals, pycoord_vals, pxcell_bounds, pycell_bounds);
  check_and_gen_bounds(projID, pxnbounds, pxlength, pxcoord_vals, pxcell_bounds, 1);
  check_and_gen_bounds(projID, pynbounds, pylength, pycoord_vals, pycell_bounds, 0);

  char p_rll_cmor[CMOR_MAX_STRING];
  int l_p_rll = strlen("grid_north_pole_longitude") + 1;
  memcpy(p_rll_cmor, "grid_north_pole_latitude\0 "
                     "grid_north_pole_longitude\0north_pole_grid_longitude\0",
         3 * l_p_rll);

  char u_rll_cmor[CMOR_MAX_STRING];
  int l_u_rll = strlen("degrees_north") + 1;
  memcpy(u_rll_cmor, "degrees_north\0degrees_east\0 degrees_east\0 ", 3 * l_u_rll);

  char p_lcc_cmor[CMOR_MAX_STRING];
  int l_p_lcc = strlen("longitude_of_central_meridian") + 1;
  memcpy(p_lcc_cmor, "standard_parallel1\0           "
                     "longitude_of_central_meridian\0latitude_of_projection_"
                     "origin\0standard_parallel2\0           ",
         4 * l_p_lcc);

  char u_cmor[CMOR_MAX_STRING];
  int l_u_cmor = 6;
    
  const char *p_rll[] = { "grid_north_pole_latitude", "grid_north_pole_longitude", "north_pole_grid_longitude", nullptr };

  const char *p_lcc[] = { "standard_parallel1", "longitude_of_central_meridian", "latitude_of_projection_origin",
                          "standard_parallel2", nullptr };
  const char *p_ps_falselsp[] = {  "straight_vertical_longitude_from_pole",
                                "latitude_of_projection_origin",                               
                                "scale_factor_at_projection_origin",
                                "false_easting",
                                "false_northing",nullptr };
  const char *p_ps_falsestandard[] = {  "straight_vertical_longitude_from_pole",
                                "latitude_of_projection_origin",                               
                                "standard_parallel",
                                "false_easting",
                                "false_northing",nullptr };
  const char *p_ps_truestandard[] = {  "straight_vertical_longitude_from_pole",
                    "latitude_of_projection_origin",      
                    "standard_parallel", nullptr };  
  const char **p_ps = nullptr;
  int sizeof_p_ps=0;
    
  int atttype, attlen;
  char attname[CDI_MAX_NAME];

  int natts;
  cdiInqNatts(projID, CDI_GLOBAL, &natts);

    
  char p_ps_cmor[CMOR_MAX_STRING];
  int l_p_ps = strlen("straight_vertical_longitude_from_pole") + 1;

  if (projtype == CDI_PROJ_STERE)
    {
      bool lsp = false;
      bool lfalse = false ;
      for (int iatt = 0; iatt < natts; ++iatt)
        {
          cdiInqAtt(projID, CDI_GLOBAL, iatt, attname, &atttype, &attlen);
          if (strcmp(attname, "standard_parallel") == 0)
            lsp = true ;
          else if (strcmp(attname, "false_easting") == 0)
            lfalse = true;
        }
      
      if ( lfalse == true )
        {
          memcpy(u_cmor, "degrees_east\0 "
                         "degrees_north\0"
                         "degrees_north\0"
                         "degrees_east\0 "
                         "degrees_north\0", 3 * l_u_rll);
          if ( lsp == false )
            {
              sizeof_p_ps=sizeof(p_ps_falselsp);
              p_ps=p_ps_falselsp;
              memcpy(p_ps_cmor, "straight_vertical_longitude_from_pole\0"
                                "latitude_of_projection_origin\0        "
                                "scale_factor_at_projection_origin\0    "
                                "false_easting\0                        "
                                "false_northing\0                       ",

                     5 * l_p_ps);

            }
          else 
            {
              sizeof_p_ps=sizeof(p_ps_falsestandard);
              p_ps=p_ps_falsestandard;
              memcpy(p_ps_cmor, "straight_vertical_longitude_from_pole\0"
                                "latitude_of_projection_origin\0        "
                                "standard_parallel\0                    "
                                "false_easting\0                        "
                                "false_northing\0                       ",

                     5 * l_p_ps);
            }
        }
      else
        {
          sizeof_p_ps=sizeof(p_ps_truestandard);
          p_ps=p_ps_truestandard;
          memcpy(p_ps_cmor, "straight_vertical_longitude_from_pole\0"
                            "latitude_of_projection_origin\0        "      
                            "standard_parallel\0                    ",
                     3 * l_p_ps);
          memcpy(u_cmor, "degrees_east\0 "
                         "degrees_north\0"
                         "degrees_north\0", 3 * l_u_rll);
        }
    }
  double *parameter_values = nullptr;

  char mapping[CDI_MAX_NAME];
  int length = CDI_MAX_NAME;
  cdiInqKeyString(projID, CDI_GLOBAL, CDI_KEY_GRIDMAP_VARNAME, mapping, &length);

  int p_len = 0;
  switch (projtype)
    {
    case CDI_PROJ_RLL: p_len = (sizeof(p_rll) / sizeof(p_rll[0])) - 1; break;
    case CDI_PROJ_LAEA:
      cdo_abort("ERROR (infile: '%s')! In grid registration:\n          This grid projection is not yet enabled.", cdo_get_stream_name(0));
      break;
    case CDI_PROJ_LCC: p_len = (sizeof(p_lcc) / sizeof(p_lcc[0])) - 1; break;
    case CDI_PROJ_SINU:
      cdo_abort("ERROR (infile: '%s')! In grid registration:\n          This grid projection is not yet enabled.",cdo_get_stream_name(0));
      break;
    case CDI_PROJ_STERE: p_len = (sizeof_p_ps / sizeof(p_ps[0])) - 1; break;
    }
  if (natts < p_len)
    cdo_warning("In grid registration:\n          Number of required grid mapping attributes '%d' is larger than the "
               "number of given grid mapping attributes '%d'.\n          Note that all required mapping attributes are "
               "set to 0.0 by default in case they are not given.",
               p_len, natts);
    
  parameter_values = (double *) Malloc(p_len * sizeof(double));
  for (int i = 0; i < p_len; ++i) parameter_values[i] = 0.0;

  for (int iatt = 0; iatt < natts; ++iatt)
    {
      cdiInqAtt(projID, CDI_GLOBAL, iatt, attname, &atttype, &attlen);
      if (atttype == CDI_DATATYPE_FLT32 || atttype == CDI_DATATYPE_FLT64)
        {
          if (attlen > 1)
            cdo_abort("ERROR (infile: '%s')! In grid registration:\n          Dont know what to do with grid mapping attribute '%s'.",
                     cdo_get_stream_name(0), attname);
          Varray<double> attflt(attlen);
          cdiInqAttFlt(projID, CDI_GLOBAL, attname, attlen, attflt.data());
          int i = 0;
          for (i = 0; i < p_len; ++i)
            {
              if (projtype == CDI_PROJ_RLL)
                {
                  if (strcmp(attname, p_rll[i]) == 0)
                    {
                      parameter_values[i] = attflt[0];
                      break;
                    }
                }
              else if (projtype == CDI_PROJ_LCC)
                {
                  if (strcmp(attname, p_lcc[i]) == 0)
                    {
                      parameter_values[i] = attflt[0];
                      break;
                    }
                }
              else if (projtype == CDI_PROJ_STERE)
                {
                  if (strcmp(attname, p_ps[i]) == 0)
                    {
                      parameter_values[i] = attflt[0];
                      break;
                    }
                }
            }
          if (i == p_len)
            cdo_warning("In grid registration:\n          grid mapping attribute '%s' is neglected.", attname);
        }
      else if (atttype == CDI_DATATYPE_TXT)
        {
          std::vector<char> atttxt(attlen + 1);
          cdiInqAttTxt(projID, CDI_GLOBAL, attname, attlen, atttxt.data());
          atttxt[attlen] = 0;
        }
    }

  int grid_axis[2];
  if (projtype == CDI_PROJ_RLL)
    {
      cmf = cmor_axis(&grid_axis[0], (char *) "grid_latitude", (char *) "degrees_north", pylength,
                      (void *) pycoord_vals, 'd',
                      (void *) pycell_bounds, 2, nullptr);
      cmf = cmor_axis(&grid_axis[1], (char *) "grid_longitude", (char *) "degrees_east", pxlength,
                      (void *) pxcoord_vals, 'd',
                      (void *) pxcell_bounds, 2, nullptr);
      cmf = cmor_grid(&grid_ids[0], 2, grid_axis, 'd', (void *) ycoord_vals, (void *) xcoord_vals, 4,
                      (void *) ycell_bounds, (void *) xcell_bounds);
#if (CMOR_VERSION_MAJOR == 2)
      cmf = cmor_set_grid_mapping(grid_ids[0], "rotated_latitude_longitude", p_len, (char **) p_rll_cmor, l_p_rll,
                                  parameter_values, (char **) u_rll_cmor, l_u_rll);
#elif (CMOR_VERSION_MAJOR == 3)
      cmf = cmor_set_grid_mapping(grid_ids[0], (char *) "rotated_latitude_longitude", p_len, p_rll_cmor, l_p_rll,
                                  parameter_values, u_rll_cmor, l_u_rll);
#endif
    }
  else if (projtype == CDI_PROJ_LCC)
    {
      memcpy(u_cmor, "      \0      \0      \0      \0", 4 * l_u_cmor);
      double *xii = (double *) Malloc(xlength * sizeof(double));
      double *yii = (double *) Malloc(ylength * sizeof(double));
      for (int i = 0; i < xlength; ++i) xii[i] = (double) i;
      for (int i = 0; i < ylength; ++i) yii[i] = (double) i;
      cmf = cmor_axis(&grid_axis[0], (char *) "x", (char *) "m", ylength, (void *) yii, 'd', 0, 0, nullptr);
      cmf = cmor_axis(&grid_axis[1], (char *) "y", (char *) "m", xlength, (void *) xii, 'd', 0, 0, nullptr);
      cmf = cmor_grid(&grid_ids[0], 2, grid_axis, 'd', (void *) ycoord_vals, (void *) xcoord_vals, 4,
                      (void *) ycell_bounds, (void *) xcell_bounds);
#if (CMOR_VERSION_MAJOR == 2)
      cmf = cmor_set_grid_mapping(grid_ids[0], mapping, p_len, (char **) p_lcc_cmor, l_p_lcc, parameter_values,
                                  (char **) u_cmor, l_u_cmor);
#elif (CMOR_VERSION_MAJOR == 3)
      cmf = cmor_set_grid_mapping(grid_ids[0], mapping, p_len, p_lcc_cmor, l_p_lcc, parameter_values, u_cmor,
                                  l_u_cmor);
#endif
      Free(xii);
      Free(yii);
    }
  else if (projtype == CDI_PROJ_STERE)
    { 
      cmf = cmor_axis(&grid_axis[0], (char *) "y", (char *) "m", ylength, (void *) pycoord_vals, 'd',
                      (void *) pycell_bounds, 2, nullptr);
      cmf = cmor_axis(&grid_axis[1], (char *) "x", (char *) "m", xlength, (void *) pxcoord_vals, 'd', 
                      (void *) pxcell_bounds, 2, nullptr);
      cmf = cmor_grid(&grid_ids[0], 2, grid_axis, 'd', (void *) ycoord_vals, (void *) xcoord_vals, 4,
                      (void *) ycell_bounds, (void *) xcell_bounds);
#if (CMOR_VERSION_MAJOR == 2)
      cmf = cmor_set_grid_mapping(grid_ids[0], "polar_stereographic", p_len, (char **) p_ps_cmor, l_p_ps, parameter_values,
                                  (char **) u_ps_cmor, l_u_ps);
#elif (CMOR_VERSION_MAJOR == 3)
      cmf = cmor_set_grid_mapping(grid_ids[0], (char *)"polar_stereographic", p_len, p_ps_cmor, l_p_ps, parameter_values, u_cmor,
                                  l_u_rll);
#endif
    }
    
  Free(parameter_values);
  Free(pxcell_bounds);
  Free(pycell_bounds);
  Free(pxcoord_vals);
  Free(pycoord_vals);
  if (cmf != 0) cdo_abort("ERROR (infile: '%s')! Function cmor_axis or cmor_set_grid_mapping failed!", cdo_get_stream_name(0));
}

static void
register_grid(KVList *kvl, int vlistID, int varID, int *axis_ids, int *grid_ids, char *project_id, char *cmor_name)
{
  int cmf = 0;
  int gridID = vlistInqVarGrid(vlistID, varID);

  char *chardimatt = kv_get_a_val(kvl, "ca", nullptr);
  char *movelons = kv_get_a_val(kvl, "ml", "y");
  char *chardim = get_txtatt(vlistID, varID, "character_axis");
  check_compare_set(&chardim, chardimatt, "character_axis", "notSet");

  char xname[CDI_MAX_NAME], yname[CDI_MAX_NAME];
  int length = CDI_MAX_NAME;
  cdiInqKeyString(gridID, CDI_XAXIS, CDI_KEY_NAME, xname, &length);
  length = CDI_MAX_NAME;
  cdiInqKeyString(gridID, CDI_YAXIS, CDI_KEY_NAME, yname, &length);

  char xdimname[CDI_MAX_NAME], ydimname[CDI_MAX_NAME];
  length = CDI_MAX_NAME;
  cdiInqKeyString(gridID, CDI_XAXIS, CDI_KEY_DIMNAME, xdimname, &length);
  length = CDI_MAX_NAME;
  cdiInqKeyString(gridID, CDI_YAXIS, CDI_KEY_DIMNAME, ydimname, &length);

  auto totalsize = gridInqSize(gridID);

  if (totalsize > 1)
    {
      int projID = gridInqProj(gridID);
      int projtype = CDI_UNDEFID;
      if ( projID != CDI_UNDEFID )
        projtype=gridInqProjType(projID);
      int type = gridInqType(gridID);
      int ylength = gridInqYsize(gridID);
      int xlength = gridInqXsize(gridID);
      double *xcoord_vals = nullptr;
      double *ycoord_vals = nullptr;
      double *xcell_bounds = nullptr;
      double *ycell_bounds = nullptr;
      int xnbounds;
      int ynbounds;

      if (type == GRID_GAUSSIAN || type == GRID_LONLAT)
        {
          grid_ids[0] = 0;
          xcoord_vals = (double *) Malloc(xlength * sizeof(double));
          ycoord_vals = (double *) Malloc(ylength * sizeof(double));
          xcell_bounds = (double *) Malloc(2 * xlength * sizeof(double));
          ycell_bounds = (double *) Malloc(2 * ylength * sizeof(double));
          inquire_vals_and_bounds(gridID, &xnbounds, &ynbounds, xcoord_vals, ycoord_vals, xcell_bounds, ycell_bounds);

          check_and_gen_bounds(gridID, xnbounds, xlength, xcoord_vals, xcell_bounds, 1);
          check_and_gen_bounds(gridID, ynbounds, ylength, ycoord_vals, ycell_bounds, 0);

          if ( ylength > 1 )
            cmf = cmor_axis(new_axis_id(axis_ids), (char *) "latitude", (char *) "degrees_north", ylength,
                          (void *) ycoord_vals, 'd', (void *) ycell_bounds, 2, nullptr);
          if ( xlength > 1 )
            cmf = cmor_axis(new_axis_id(axis_ids), (char *) "longitude", (char *) "degrees_east", xlength,
                          (void *) xcoord_vals, 'd', (void *) xcell_bounds, 2, nullptr);

          Free(xcell_bounds);
          Free(ycell_bounds);
          Free(xcoord_vals);
          Free(ycoord_vals);
        }
      else if ( type == GRID_UNSTRUCTURED )
        {
          int nvertex = gridInqNvertex(gridID);
          xcoord_vals = (double *) Malloc(totalsize * sizeof(double));
          ycoord_vals = (double *) Malloc(totalsize * sizeof(double));
          /* maximal 4 gridbounds per gridcell permitted */
          if ( nvertex )
            {
              xcell_bounds = (double *) Malloc(nvertex * totalsize * sizeof(double));
              ycell_bounds = (double *) Malloc(nvertex * totalsize * sizeof(double));
            }
          inquire_vals_and_bounds(gridID, &xnbounds, &ynbounds, xcoord_vals, ycoord_vals, xcell_bounds, ycell_bounds);
          /* In a projection, this is done by setting mapping parameter */
          if ( strcmp(movelons, "y") == 0 )
            move_lons(xcoord_vals, xcell_bounds, totalsize, nvertex * totalsize, xnbounds);
          int grid_axis[2];
          double *coord_vals;
          coord_vals = (double *) Malloc(xlength * sizeof(double));
          for (int j = 0; j < xlength; ++j) coord_vals[j] = (double) j;
          if ( strcmp(chardim, "site") == 0 )
            {
              cmf = cmor_axis(&grid_axis[0], (char *) "site", (char *) "1", xlength, (void *) coord_vals, 'd', 0,
                          0, nullptr);
              get_cmor_table(kvl, project_id);
              cmf = cmor_grid(&grid_ids[0], 1, grid_axis, 'd', (void *) ycoord_vals, (void *) xcoord_vals, 0, nullptr, nullptr);
              kv_insert_vals(kvl, "capr", (char *)"y", true, false);
            }
          else
            {
              get_cmor_table(kvl, project_id);
              cmf = cmor_axis(&grid_axis[0], (char *) "i_index", (char *) "1", xlength, (void *) coord_vals, 'd', 0,
                          0, nullptr);
              cmf = cmor_grid(&grid_ids[0], 1, grid_axis, 'd', (void *) ycoord_vals, (void *) xcoord_vals, nvertex,
                          (void *) ycell_bounds, (void *) xcell_bounds);
            }
          Free(coord_vals);
          Free(xcoord_vals);
          Free(ycoord_vals);
          if (xcell_bounds) Free(xcell_bounds);
          if (ycell_bounds) Free(ycell_bounds);
        }
      else if (type == GRID_CURVILINEAR )
        {
          xcoord_vals = (double *) Malloc(totalsize * sizeof(double));
          ycoord_vals = (double *) Malloc(totalsize * sizeof(double));
          /* maximal 4 gridbounds per gridcell permitted */
          xcell_bounds = (double *) Malloc(4 * totalsize * sizeof(double));
          ycell_bounds = (double *) Malloc(4 * totalsize * sizeof(double));
          inquire_vals_and_bounds(gridID, &xnbounds, &ynbounds, xcoord_vals, ycoord_vals, xcell_bounds, ycell_bounds);
          /* In a projection, this is done by setting mapping parameter */
          if ( strcmp(movelons, "y") == 0 )
            move_lons(xcoord_vals, xcell_bounds, totalsize, 4 * totalsize, xnbounds);
          get_cmor_table(kvl, project_id);
          int grid_axis[2];
          check_and_gen_bounds_curv(gridID, totalsize, xnbounds, xlength, xcoord_vals, xcell_bounds, ynbounds, ylength,
                                    ycoord_vals, ycell_bounds);
          if ( projID == CDI_UNDEFID || projtype == CDI_UNDEFID )
            {
              double *xncoord_vals;
              double *yncoord_vals;
              xncoord_vals = (double *) Malloc(xlength * sizeof(double));
              yncoord_vals = (double *) Malloc(ylength * sizeof(double));
              for (int j = 0; j < ylength; ++j) yncoord_vals[j] = (double) j;
              for (int j = 0; j < xlength; ++j) xncoord_vals[j] = (double) j;
              cmf = cmor_axis(&grid_axis[0], (char *) "j_index", (char *) "1", ylength, (void *) yncoord_vals, 'd', 0,
                              0, nullptr);
              cmf = cmor_axis(&grid_axis[1], (char *) "i_index", (char *) "1", xlength, (void *) xncoord_vals, 'd', 0,
                              0, nullptr);
              cmf = cmor_grid(&grid_ids[0], 2, grid_axis, 'd', (void *) ycoord_vals, (void *) xcoord_vals, 4,
                              (void *) ycell_bounds, (void *) xcell_bounds);
              Free(xncoord_vals);
              Free(yncoord_vals);
              Free(xcoord_vals);
              Free(ycoord_vals);
              Free(xcell_bounds);
              Free(ycell_bounds);
            }
          /*else
            {
              cmf = cmor_axis(&grid_axis[0],    "grid_longitude",   "degrees",    xlength,    (void *)xcoord_vals,
            'd', 0, 0, nullptr);
              cmf = cmor_axis(&grid_axis[1],    "grid_latitude",    "degrees",    ylength,    (void *)ycoord_vals,
            'd', 0, 0, nullptr);
              cmf = cmor_grid(&grid_ids[0],    2,    grid_axis,    'd',    (void *)ycoord_vals,    (void *)xcoord_vals,
            2,     (void *)ycell_bounds,    (void *)xcell_bounds);
            }*/
        }
      else if (type == GRID_GENERIC && chardim && ( (!strstr(xname, "lon") && !strstr(xdimname, "lon")) ||
                                                  (!strstr(yname, "lat") && !strstr(ydimname, "lat")) )


 )/* && (strcmp(chardim, "oline") == 0 || strcmp(chardim, "basin") == 0 || strcmp(chardim, "siline") == 0)) */
        {
          if (Options::cdoVerbose) cdo_print("Unknown grid type.");
          if (Options::cdoVerbose) cdo_print("Start to define a character axis '%s' instead of a grid axis'.", chardim);
          grid_ids[0] = 0;
          int numchar = 0;
          char *charvalstring = (char *) Malloc(CMOR_MAX_STRING * sizeof(char));
          sprintf(charvalstring, "char_axis_%s_%s", chardim, cmor_name);
          std::vector<std::string> charvals = kv_get_vals(kvl, charvalstring, &numchar);
          if ( numchar == 0 )
            {
              sprintf(charvalstring, "char_axis_%s", chardim);
              charvals = kv_get_vals(kvl, charvalstring, &numchar);
            }
          Free(charvalstring);
          if ((xlength > 0 && xlength != numchar) && (ylength > 0 && ylength != numchar))
            cdo_abort("ERROR (infile: '%s')! In registration of a character coordinate as substitution for a horizontal axis:\n        "
                     "  You configured a character coordinate '%s' with '%d' string values but you also registered a "
                     "grid with '%d' numerical values on X axis and '%d' numerical values on Y axis. One axis must "
                     "match the number of string values.",
                     cdo_get_stream_name(0), chardim, numchar, xlength, ylength);
          if (!charvals.size())
            cdo_abort("ERROR (infile: '%s')! In registration of a character coordinate as substitution for a horizontal axis:\n        "
                     "  You configured a character coordinate '%s' but no values are found! Configure values via "
                     "attribute 'char_axis_vals'!",
                     cdo_get_stream_name(0), chardim);
          if (charvals.size() && (xlength == numchar || xlength == 0))
            {
              register_char_axis(numchar, charvals, axis_ids, chardim);
              if (ylength > 1) register_lat_axis(gridID, ylength, axis_ids);
            }
          else
            {
              if (xlength > 1 ) register_lon_axis(gridID, xlength, axis_ids);
              register_char_axis(numchar, charvals, axis_ids, chardim);
            }
          if (Options::cdoVerbose) cdo_print("Successfully defined a character axis '%s' instead of a grid axis.", chardim);
          kv_insert_vals(kvl, "capr", (char *)"y", true, false);
        }
      else if (type == GRID_CHARXY)
        {
          grid_ids[0] = 0;
          if (strcmp(xdimname, "line") == 0) strcpy(xdimname, "oline");
          int dimstrlen;
          if ((dimstrlen = gridInqXIsc(gridID)))
            {
              std::vector<std::string> xcharspp;
              char **xchars = (char **) Malloc((xlength + 1) * sizeof(char *));
              for (int i = 0; i < xlength; ++i) xchars[i] = (char *) Malloc((dimstrlen + 1) * sizeof(char));
              gridInqXCvals(gridID, xchars);
              for (int j = 0; j < xlength; ++j) xchars[j][dimstrlen] = 0;
              xchars[xlength] = nullptr;
              for (int j = 0; j < xlength; ++j) xcharspp.push_back(xchars[j]);
              register_char_axis(xlength, xcharspp, axis_ids, xdimname);
              free_array(xchars);
            }
          else if (xlength)
            register_lon_axis(gridID, xlength, axis_ids);

          if ((dimstrlen = gridInqYIsc(gridID)))
            {
              std::vector<std::string> ycharspp;
              char **ychars = (char **) Malloc((ylength + 1) * sizeof(char));
              for (int i = 0; i < ylength; ++i) ychars[i] = (char *) Malloc((dimstrlen + 1) * sizeof(char));
              gridInqYCvals(gridID, ychars);
              for (int j = 0; j < ylength; ++j) ychars[j][dimstrlen] = 0;
              ychars[ylength] = nullptr;
              for (int j = 0; j < ylength; ++j) ycharspp.push_back(ychars[j]);
              register_char_axis(ylength, ycharspp, axis_ids, ydimname);
              free_array(ychars);
            }
          else if (ylength)
            register_lat_axis(gridID, ylength, axis_ids);
        }
      else if (type == GRID_PROJECTION)
        {
          cdo_abort("ERROR (infile: '%s')! In grid registration:\n          For a 'rotated_lat_lon' projection, both grids, the "
                   "unprojected lat/lon and the projected rlat/rlon are required.", cdo_get_stream_name(0));
        }
      else
        {
          grid_ids[0] = 0;
          cdo_warning(
              "Registration of a grid is skipped. Either the grid type is unknown or a registration is not necessary.");
        }

      if (projtype != CDI_UNDEFID)
        {
          register_projection(grid_ids, projID, ycoord_vals, xcoord_vals, 
                              ycell_bounds, xcell_bounds, xlength, ylength, projtype);
          Free(xcoord_vals);
          Free(ycoord_vals);
          Free(xcell_bounds);
          Free(ycell_bounds);
        }
    }
  else
    grid_ids[0] = 0;
  Free(chardim);
  if (cmf != 0) cdo_abort("ERROR (infile: '%s')! Function cmor_axis failed!", cdo_get_stream_name(0));
}

static void
register_variable(KVList *kvl, int vlistID, int varID, int *axis_ids, struct mapping *var, int *grid_ids, char *name)
{
  int cmf = 0;
  if (Options::cdoVerbose) cdo_print("8.5.1. Start to retrieve 'positive' and 'units'.");
  char *positive = get_txtatt(vlistID, varID, "positive");
  char *origname = get_txtatt(vlistID, varID, "original_name");
  char *history = get_txtatt(vlistID, varID, "history");
  char *varcom = get_txtatt(vlistID, varID, "variable_comment");
  char *units = (char *) Malloc(CDI_MAX_NAME * sizeof(char));
  vlistInqVarUnits(vlistID, varID, units);
  char *attunits = kv_get_a_val(kvl, "u", nullptr);
  char *attp = kv_get_a_val(kvl, "p", nullptr);
  char *attorigname = kv_get_a_val(kvl, "original_name", nullptr);
  char *attvarcom = kv_get_a_val(kvl, "vc", nullptr);
  check_compare_set(&positive, attp, "positive", "");
  if (strcmp(positive, " ") == 0) strcpy(positive, "");
  check_compare_set(&units, attunits, "units", nullptr);
  check_compare_set(&origname, attorigname, "original_name", "");
  if (strcmp(origname, "") == 0 || strstr(origname, "var"))
    {
      Free(origname);
      origname = nullptr;
    }
  check_compare_set(&varcom, attvarcom, "variable_comment", "");
  if (strcmp(varcom, "") == 0)
    {
      Free(varcom);
      varcom = nullptr;
    }
  if (Options::cdoVerbose) cdo_print("8.5.1. Successfully retrieved 'positive': '%s' and 'units' : '%s' and 'variable_comment': '%s'", positive, units, varcom);
  if ( strcmp(units, "") == 0 )
    cdo_abort("ERROR (infile: '%s')! No units found for CMOR variable '%s'.", cdo_get_stream_name(0), name);
  char missing_value[sizeof(double)];
  double tolerance = 1e-4;
  size_t gridsize = vlistGridsizeMax(vlistID);
  int zsize = zaxisInqSize(vlistInqVarZaxis(vlistID, varID));
  var->help_var = 0;
  if ( vlistInqVarDatatype(vlistID, varID) == CDI_DATATYPE_FLT64 )
    {
      if ( !var->data )
        {
          var->charvars = 0;
          var->datatype = 'd';
          var->data = Malloc(gridsize * zsize * sizeof(double));
        }
      *(double *) missing_value = vlistInqVarMissval(vlistID, varID);
    }
  else
    {
      if ( !var->data )
        {
          var->charvars = 0;
          var->datatype = 'f';
          var->data = Malloc(gridsize * zsize * sizeof(float));
        }
      *(float *) missing_value = vlistInqVarMissval(vlistID, varID);
    }
  if (Options::cdoVerbose) cdo_print("8.5.2. Start to call cmor_variable.");
  if ( (zaxisInqType(vlistInqVarZaxis(vlistID, varID)) != ZAXIS_HYBRID && 
        zaxisInqType(vlistInqVarZaxis(vlistID, varID)) != ZAXIS_HYBRID_HALF) && grid_ids[0] != 0)
    {
      int *tmp_id = new_axis_id(axis_ids);
      *tmp_id = grid_ids[0];
      cmf = cmor_variable(&var->cmor_varID, name, units, (count_axis_ids(axis_ids)), axis_ids, var->datatype,
                          (void *) missing_value, &tolerance, positive, origname, history,
                          varcom);
    }
  else
    {
      cmf = cmor_variable(&var->cmor_varID, name, units, count_axis_ids(axis_ids), axis_ids, var->datatype,
                          (void *) missing_value, &tolerance, positive, origname, history,
                          varcom);
    }
  if (cmf != 0) cdo_abort("ERROR (infile: '%s')! Function cmor_variable failed!", cdo_get_stream_name(0));
  if (Options::cdoVerbose) cdo_print("8.5.2. Successfully called cmor_variable.");
#if (CMOR_VERSION_MAJOR == 3 && CMOR_VERSION_MINOR >= 3)
  char *deflates = NULL;
  long int deflate = 0;
  if ( ( deflates = kv_get_a_val(kvl, "dl", nullptr) ) )
    {
      if (Options::cdoVerbose) cdo_print("8.5.3. Start to set deflate for variable '%s'.", name);
      if ( ( deflate = atol((const char *)deflates) ) )
        {
          if ( deflate == -1 )
            cmf = cmor_set_deflate(var->cmor_varID, 0, 0, 0);
          else
            cmf = cmor_set_deflate(var->cmor_varID, 1, 1, deflate);
        }
      else
        cmf = cmor_set_deflate(var->cmor_varID, 1, 1, 1);
      if (cmf != 0) cdo_abort("ERROR (infile: '%s')! Function cmor_variable failed!", cdo_get_stream_name(0));
      if (Options::cdoVerbose) cdo_print("8.5.3. Successfully set deflate for variable '%s'.", name);
    }
#endif
  if (positive) Free(positive);
  if (origname) Free(origname);
  if (history) Free(history);
  if (units) Free(units);
  if (varcom) Free(varcom);
}

static void
switch_grid_info(KVList *kvl, CdoStreamID streamID, char *grid_file, int varID)
{
  if (Options::cdoVerbose)
    cdo_print("You configured a grid_info file: '%s'. It is tested for a valid use as substitution.\n", grid_file);
  int vlistID = cdo_stream_inq_vlist(streamID);
  int nvars = vlistNvars(vlistID);
  if (nvars > 1) cdo_print("Note that the grids of the variables found first in both files are switched.");

  int byteorder = 0;
  int filetype = cdiGetFiletype(grid_file, &byteorder);
  if ((filetype == CDI_FILETYPE_NC) || (filetype == CDI_FILETYPE_NC2) || (filetype == CDI_FILETYPE_NC4)
      || (filetype == CDI_FILETYPE_NC4C))
    {
      int streamID2 = stream_open_read_locked(grid_file);
      int vlistID2 = streamInqVlist(streamID2);
      int gridID2 = vlistInqVarGrid(vlistID2, 0);
      int zaxisID2 = vlistInqVarZaxis(vlistID2, 0);

      if (strcmp(kv_get_a_val(kvl, "switch_xy", "y"), "y") == 0)
        {
          int gridID = vlistInqVarGrid(vlistID, varID);
          change_grid(vlistID, gridID, gridID2, grid_file);
        }

      if (strcmp(kv_get_a_val(kvl, "switch_z", "y"), "y") == 0)
        {
          int zaxisID = vlistInqVarZaxis(vlistID, varID);
          change_zaxis(kvl, vlistID, zaxisID, zaxisID2, grid_file);
        }

      streamClose(streamID2);
    }
  else
    {
      if (parse_kv_file(kvl, grid_file) == 0) cdo_abort("ERROR (infile: '%s')! File '%s' does not exist.", cdo_get_stream_name(0), grid_file);
    }
}

static void
register_all_dimensions(KVList *kvl, CdoStreamID streamID, struct mapping vars[], int table_id, char *project_id,
                        int miptab_freq, int *time_axis, int *mergeIDs)
{
  int cmf = 0;
  int vlistID = cdo_stream_inq_vlist(streamID);

  char *time_units = kv_get_a_val(kvl, "rtu", nullptr);

  if (Options::cdoVerbose) cdo_print("7. Start to retrieve requested variables.");

  int numvals = 0;
  std::vector<std::string> cmor_names = kv_get_vals(kvl, "cn", &numvals);

  /* Cmdlinemapping: */
  char *mapname, *mapcode;
  if (!kv_get_a_val(kvl, "mt", nullptr) && numvals)
    {
      if (Options::cdoVerbose) cdo_print("7.1. Start to search for a command line mapping");
      if ((mapname = kv_get_a_val(kvl, "n", nullptr)))
        {
          if ( change_name_via_name(vlistID, mapname, cmor_names[0].c_str()) )
            {if (Options::cdoVerbose) cdo_print("7.1. Successfully mapped '%s' on '%s' via command line attribute name.", mapname, cmor_names[0]);}
        }
      else if ((mapcode = kv_get_a_val(kvl, "c", nullptr)))
        {
          if ( change_name_via_code(vlistID, mapcode, cmor_names[0].c_str()) )
            {if (Options::cdoVerbose) cdo_print("7.1. Successfully mapped via command line attribute name.");}
        }
      else
        if (Options::cdoVerbose) cdo_print("7.1. No command line mapping given.");
    }

  if (!cmor_names.size() && vlistNvars(vlistID) > 1)
    cdo_print("Function 'all axes registration':\n          You have not requested a particular variable via "
             "'cmor_name'.\n          There are several in infile and all will be processed.\n          Notice that "
             "attributes specified in the cmdline will be used for all infile variables.");
  if (Options::cdoVerbose) cdo_print("7. Successfully retrieved requested variables");
  int foundName = 0, psRequired = 0;
  std::map<int,int> timeAxes;

  for (int varID = 0; varID < vlistNvars(vlistID); ++varID)
    {
      char name[CDI_MAX_NAME];
      vlistInqVarName(vlistID, varID, name);
      if (!cmor_names.size() || in_list(cmor_names, name, numvals))
        {
          struct mapping *var = map_var(varID, vars);
          if (!var) var = new_var_mapping(vars);
          else if (Options::cdoVerbose) cdo_print("Already mapped '%d'", varID);
          var->cdi_varID = varID;
          char *grid_file = kv_get_a_val(kvl, "gi", nullptr);
          if (grid_file) switch_grid_info(kvl, streamID, grid_file, varID);
          int axis_ids[CMOR_MAX_AXES];
          axis_ids[0] = CMOR_UNDEFID;
          int zaxisID = vlistInqVarZaxis(vlistID, varID);
          int zsize = zaxisInqSize(zaxisID);
          if (Options::cdoVerbose) cdo_print("8. Start to define variable with ID: '%d' and name: '%s'", varID, name);
          if ( (zaxisInqType(zaxisID) == ZAXIS_HYBRID || zaxisInqType(zaxisID) == ZAXIS_HYBRID_HALF ) && zsize > 1 )
            {
              if (Options::cdoVerbose)
                cdo_print("Since the zaxis of variable '%s' is of type HYBRID, surface pressure needs to be available in the input file to fully describe the axis.",
                       name);
              psRequired++;
            }
          /* Time-Axis */
          if (Options::cdoVerbose) cdo_print("8.1. Start to register time axis of '%s'", name);
          char cmor_time_name[CMOR_MAX_STRING];
          cmor_time_name[0] = '\0';
          get_time_method(kvl, vlistID, varID, cmor_time_name, project_id, miptab_freq, time_axis);
          if ( strcmp(cmor_time_name, "none") != 0 )
            {
              auto search = timeAxes.find(*time_axis);
              int *cmorTimeAxisID = new_axis_id(axis_ids);
              if ( foundName )
                {
                  if (search != timeAxes.end())
                    {
                      *cmorTimeAxisID = search->second;
                      if (Options::cdoVerbose)
                        cdo_print("8.1. Use already defined time axis '%d'", *cmorTimeAxisID);
                    }
                }
              if ( *cmorTimeAxisID == CMOR_UNDEFID )
                {
                  cmf = cmor_axis(cmorTimeAxisID, cmor_time_name, time_units, 0, nullptr, 0, nullptr, 0, nullptr);
                  if ( (strcmp(project_id, "CMIP5") == 0 || strcmp(project_id, "CORDEX") == 0 ) &&
                   (strcmp(kv_get_a_val(kvl, "realization", nullptr), "0") == 0 ||
                    strcmp(kv_get_a_val(kvl, "initialization_method", nullptr), "0") == 0 ||
                    strcmp(kv_get_a_val(kvl, "physics_version", nullptr), "0") == 0 ) )
                    cdo_warning("At least one ensemble index is set to '0' while cell_methods is not 'none'!\n"

                           "          '0' is usually reserved for fixed fields!");
                }
              if (search == timeAxes.end())
                timeAxes.insert(std::pair<int,int>(*time_axis, axis_ids[count_axis_ids(axis_ids)-1]));
            }
          if (Options::cdoVerbose && cmf == 0)
            cdo_print("8.1. Successfully handled time axis registration.");
          else if (cmf != 0)
            cdo_abort("ERROR (infile: '%s')! Function cmor_axis failed!", cdo_get_stream_name(0));
          /* Grid: */
          if (Options::cdoVerbose) cdo_print("8.2. Start to register grid of '%s'", name);
          int grid_ids[CMOR_MAX_GRIDS];
          register_grid(kvl, vlistID, varID, axis_ids, grid_ids, project_id, name);
          if ( (zaxisInqType(zaxisID) == ZAXIS_HYBRID || zaxisInqType(zaxisID) == ZAXIS_HYBRID_HALF) && grid_ids[0] != 0 )
            {
              int *tmp_id = new_axis_id(axis_ids);
              *tmp_id = grid_ids[0];
            }
          cmf = cmor_set_table(table_id);
          if (cmf != 0) cdo_abort("ERROR (infile: '%s')! Function cmor_set_table failed!", cdo_get_stream_name(0));
          if (Options::cdoVerbose) cdo_print("8.2. Successfully handled grid registration.");
          if ( miptab_freq == 19 )
            {
          /* Possible 4th spatial dimension */
          /* Z-Axis */
              if (Options::cdoVerbose) cdo_print("8.4. Start to register 4th spatial dimension of '%s'", name);
              register_fourth_axis(kvl, vlistID, varID, name, axis_ids, project_id, miptab_freq, mergeIDs);
              if (Options::cdoVerbose) cdo_print("8.4. Successfully handled 4th spatial dimension registration.");
          /* Z-Axis */
              if (Options::cdoVerbose) cdo_print("8.3. Start to register zaxis of '%s'", name);
              register_z_axis(kvl, vlistID, varID, zaxisID, name, axis_ids, &var->zfactor_id, project_id, miptab_freq, mergeIDs, vars, streamID);
              if (Options::cdoVerbose) cdo_print("8.3. Successfully handled zaxis registration.");
            }
          else
            {
          /* Z-Axis */
              if (Options::cdoVerbose) cdo_print("8.3. Start to register zaxis of '%s'", name);
              register_z_axis(kvl, vlistID, varID, zaxisID, name, axis_ids, &var->zfactor_id, project_id, miptab_freq, mergeIDs, vars, streamID);
              if (Options::cdoVerbose) cdo_print("8.3. Successfully handled zaxis registration.");
          /* Possible 4th spatial dimension */
          /* Z-Axis */
              if (Options::cdoVerbose) cdo_print("8.4. Start to register 4th spatial dimension of '%s'", name);
              register_fourth_axis(kvl, vlistID, varID, name, axis_ids, project_id, miptab_freq, mergeIDs);
              if (Options::cdoVerbose) cdo_print("8.4. Successfully handled 4th spatial dimension registration."); 
            }
          /* Variable */
          if (Options::cdoVerbose) cdo_print("8.5. Start to register variable '%s'", name);
          register_variable(kvl, vlistID, varID, axis_ids, var, grid_ids, name);
          if (Options::cdoVerbose) cdo_print("8.5. Successfully handled variable registration.");
          if (Options::cdoVerbose) cdo_print("8. Successfully defined variable with ID: '%d' and name: '%s'.", varID, name);
          kv_insert_vals(kvl, "capr", (char *)"n", true, false);
          foundName++;
        }
    }
  if ( mergeIDs[0] != CMOR_UNDEFID )
    {
      int refdatatype = vlistInqVarDatatype(vlistID, mergeIDs[0]);
      int mergeID = 0;
      while (mergeIDs[mergeID] != CMOR_UNDEFID)
        {
          struct mapping *var = map_var(mergeIDs[mergeID], vars);
          if (!var) 
            {
              var = new_var_mapping(vars);  
              var->charvars = 0;
              var->cdi_varID = mergeIDs[mergeID];
              if ( vlistInqVarDatatype(vlistID, mergeIDs[mergeID]) != refdatatype )
                cdo_abort("ERROR (infile: '%s')! Variable with ID '%d' has datatype '%d' but"
                         " variable with id '%d' has datatype '%d'."
                         " All variables that should be merged needs to be"
                         " of the same datatype.", cdo_get_stream_name(0), mergeIDs[0], refdatatype, mergeIDs[mergeID], vlistInqVarDatatype(vlistID, mergeIDs[mergeID])); 
              if ( vlistInqVarDatatype(vlistID, mergeIDs[mergeID]) == CDI_DATATYPE_FLT64 )
                {
                  var->datatype = 'd';
                  var->data = Malloc(gridInqSize(vlistInqVarGrid(vlistID, mergeIDs[mergeID])) * 
                                 zaxisInqSize(vlistInqVarZaxis(vlistID, mergeIDs[mergeID])) * 
                                 sizeof(double));
                }
              else
                {
                  var->datatype = 'f';
                  var->data = Malloc(gridInqSize(vlistInqVarGrid(vlistID, mergeIDs[mergeID])) * 
                                 zaxisInqSize(vlistInqVarZaxis(vlistID, mergeIDs[mergeID])) * 
                                 sizeof(float));
                }
            }
          mergeID++;
        }
    }
  if ( psRequired && cmor_names.size() && !in_list(cmor_names, "ps", numvals) )
    {
      int psindex = getVarIDToMap(vlistID, vlistNvars(vlistID), "name", "ps");
      if (psindex == CDI_UNDEFID)
        psindex = getVarIDToMap(vlistID, vlistNvars(vlistID), "code", "134");
      if (psindex != CDI_UNDEFID)
        {    
          int psID = getRegisteredPsid(vars, psindex);
          if ( psID == CDI_UNDEFID )
            cdo_abort("ERROR (infile: '%s')! Could not find ps for hybrid axis."); 
          else if (Options::cdoVerbose)
            cdo_print("9. Set ps as a auxilliary variable.");
          vars[psID].help_var = 1;
        }
    }
  if (!foundName && cmor_names.size())
    cdo_abort("ERROR (infile: '%s')! After registration of all dimensions for all variables:\n          None of the given variables to "
             "process by attribute 'cmor_name' is found in infile.", cdo_get_stream_name(0));
  if (Options::cdoVerbose) cdo_print("Successfully registered all dimensions for %d variables successfully.", foundName);
}

static char *
get_frequency(/*KVList *kvl,*/ int vlistID, int miptab_freq)
{
  char *frequency = (char *) Malloc(CMOR_MAX_STRING * sizeof(char));
  strcpy(frequency, "no");
  int ntsteps = vlistNtsteps(vlistID);
  int reccounter = 0;
  int recdummy = 0;

  switch (miptab_freq)
    {
    case 11: strcpy(frequency, "yr"); break;
    case 2: strcpy(frequency, "yr"); break;
    case 12: strcpy(frequency, "mon"); break;
    case 3: strcpy(frequency, "mon"); break;
    case 13: strcpy(frequency, "day"); break;
    case 4: strcpy(frequency, "day"); break;
    case 14: strcpy(frequency, "6hr"); break;
    case 5: strcpy(frequency, "6hr"); break;
    case 6: strcpy(frequency, "6hr"); break;
    case 15: strcpy(frequency, "3hr"); break;
    case 7: strcpy(frequency, "1hr"); break;
    case 8: strcpy(frequency, "3hr"); break;
    case 16: strcpy(frequency, "1hr"); break;
    case 17: strcpy(frequency, "sem"); break;
    case 18: strcpy(frequency, "dec"); break;
    case 19: strcpy(frequency, "subhr"); break;
    default:
      {
        if (cdo_assert_files_only() == false)
          {
            cdo_abort("ERROR (infile: '%s')! No frequency could be determined from MIP-table and, additionally,"
                     " cdo cmor cannot check frequency of "
                     "Ifile recs since you piped several cdo operators.", cdo_get_stream_name(0));
            /*          char *dummy;
                      cdo_warning("Cdo cmor cannot check frequency of Ifile recs since you piped several cdo
               operators.\nIt is tried to use a configuration attribute frequency.");
                      if ( !(dummy = kv_get_a_val(kvl, "frequency", nullptr)) )
                        cdo_abort("ERROR (infile: '%s')! No attribute frequency is found.");
                      else
                        {
                          strcpy(frequency, dummy);
                          return frequency;
                        }
            */
          }

        CdiStreamID streamID2 = streamOpenRead(cdo_get_stream_name(0));
        int vlistID2 = streamInqVlist(streamID2);
        int taxisID2 = vlistInqTaxis(vlistID2);
        if (ntsteps < 0)
          {
            while ((recdummy = streamInqTimestep(streamID2, reccounter++)))
              ;
            ntsteps = reccounter;
          }
        ntsteps -= 1;
        int fyear, lyear, fmonth, lmonth, dummytwo;

        if (ntsteps > 2)
          {
            streamInqTimestep(streamID2, 0);
            cdiDate_decode(taxisInqVdatetime(taxisID2).date, &fyear, &fmonth, &dummytwo);
            streamInqTimestep(streamID2, ntsteps);
            cdiDate_decode(taxisInqVdatetime(taxisID2).date, &lyear, &lmonth, &dummytwo);

            double covered_years = lyear - fyear + 1.0;
            double ntperyr = (double) ((ntsteps + 1) / covered_years);
            if (DBL_IS_EQUAL(ntperyr, (double) 1))
              strcpy(frequency, "yr");
            else if (DBL_IS_EQUAL(ntperyr, (double) 12))
              strcpy(frequency, "mon");
            else if (DBL_IS_EQUAL(ntperyr, (double) 365) || DBL_IS_EQUAL(ntperyr, (double) 365.25)
                     || DBL_IS_EQUAL(ntperyr, (double) 366))
              strcpy(frequency, "day");
            else if (DBL_IS_EQUAL(ntperyr, (double) 365 * 4) || DBL_IS_EQUAL(ntperyr, (double) 365.25 * 4)
                     || DBL_IS_EQUAL(ntperyr, (double) 366 * 4))
              strcpy(frequency, "6hr");
            else if (DBL_IS_EQUAL(ntperyr, (double) 365 * 8) || DBL_IS_EQUAL(ntperyr, (double) 365.25 * 8)
                     || DBL_IS_EQUAL(ntperyr, (double) 366 * 8))
              strcpy(frequency, "3hr");
            else
              {
                int step_per_year = 0;
                reccounter = 0;
                if (Options::cdoVerbose)
                  cdo_print("Frequency could not be determined by comparing all time steps (%d) divided by covered "
                           "years (%f).\n          It is now calculated by counting all timesteps in year %d\n         "
                           " in order to calculate time bounds in case they are not given.",
                           ntsteps, covered_years, fyear);
                while ((recdummy = streamInqTimestep(streamID2, reccounter++)))
                  {
                    int reqyear;
                    cdiDate_decode(taxisInqVdatetime(taxisID2).date, &reqyear, &lmonth, &dummytwo);
                    if (reqyear == (fyear + 1)) break;
                    step_per_year++;
                  }
                int covered_months = lmonth - fmonth + 1;
                if (step_per_year > 366 * 8)
                  cdo_abort("ERROR (infile: '%s')! In estimating frequency:\n          Frequency is sub-3hourly! Not yet enabled.", cdo_get_stream_name(0));
                else
                  {
                    if ((double) step_per_year / (double) covered_months > 31 * 8)
                      cdo_abort("ERROR (infile: '%s')! Frequency is sub-3hourly! Not yet enabled.", cdo_get_stream_name(0));
                    else if ((double) step_per_year / (double) covered_months > 31 * 4)
                      strcpy(frequency, "3hr");
                    else if ((double) step_per_year / (double) covered_months > 31)
                      strcpy(frequency, "6hr");
                    else if ((double) step_per_year / (double) covered_months > 1)
                      strcpy(frequency, "day");
                    else
                      strcpy(frequency, "mon");
                  }
                if (Options::cdoVerbose)
                  cdo_print("Found %d time steps in year %d.\n          Therefore, the frequency is %s.", step_per_year,
                           fyear, frequency);
              }
          }
        else
          {
            if (!taxisHasBounds(taxisID2) && ntsteps > 0)
              cdo_abort("ERROR (infile: '%s')! In estimating frequency:\n          No time bounds are found in infile and for %d found "
                       "timesteps no frequency can be computed - at least 3 timesteps are required.\n          Define "
                       "time bounds before cdo cmor.",
                       cdo_get_stream_name(0), ntsteps);
            else
              cdo_warning("In frequency estimation:\n          For %d found timesteps no frequency can be computed - at "
                         "least 3 timesteps are required.\n          Time bounds of the rec are used.",
                         ntsteps);
          }
        streamClose(streamID2);
      }
    }
  return frequency;
}

static int
get_tunitsec(int tunit)
{
  switch (tunit)
    {
    case TUNIT_MINUTE: return 60;
    case TUNIT_HOUR: return 3600;
    case TUNIT_DAY: return 86400;
    default: return 3600;
    }
}

static JulianDate
get_cmor_time_val(KVList *kvl, int taxisID, JulianDate ref_date, int /*tunitsec*/, int calendar, char *frequency, int ts_id, int time_axis)
{
  auto vDateTime = taxisInqVdatetime(taxisID);
  int year, month, day, hour, min, sec, ms;
  cdiDate_decode(vDateTime.date, &year, &month, &day);
  auto juldate = julianDate_encode(calendar, vDateTime);

  if (month == 0 || day == 0)
    {
      int timeoffset;
      if ((timeoffset = atol(kv_get_a_val(kvl, "firsttimeval", "-99"))) < 0)
        cdo_abort("ERROR (infile: '%s')! Time axis is broken (month or day = 0).\n          Provide 'timeoffset' and the operator "
                 "tries to calculate time values with frequency.", cdo_get_stream_name(0));
      else
        {
          int ryear, rmonth, rday, addseconds = 0;
          auto rDateTime = julianDate_decode(calendar, ref_date);
          cdiDate_decode(rDateTime.date, &ryear, &rmonth, &rday);
          /* Only print this for the first time step */
          if (ts_id < 2)
            cdo_warning("In writing the data:\n          Time axis is broken (month or day = 0). It is tried to "
                       "calculate time values with frequency and timeoffset ignoring the time stamp year.\n          "
                       "Note: These can only be valid if\n           - cm=m \n           - a equally spaced "
                       "monotonical time axis exist according to the frequency \n           - a correct calendar "
                       "exist!");
          /**/
          /* First record is valid for correfdate = refdate + timeoffset. */
          /**/
          while (timeoffset != 0)
            {
              if (timeoffset > 11)
                {
                  ryear += 1;
                  timeoffset -= 12;
                }
              else if (timeoffset != 0)
                {
                  rmonth += timeoffset;
                  if (rmonth > 12)
                    {
                      ryear += 1;
                      rmonth -= 12;
                    }
                  timeoffset = 0;
                }
            }
          rDateTime.date = cdiDate_encode(ryear, rmonth, rday);
          ref_date = julianDate_encode(calendar, rDateTime);
          /**/
          /* Add time index * frequency on correfdate */
          /**/
          if (strcmp(frequency, "yr") == 0)
            {
              year = ryear + ts_id;
              month = 6; /* Is set to mid point by CMOR */
              day = 14;  /* Is set to mid point by CMOR */
            }
          else if (strcmp(frequency, "mon") == 0)
            {
              year = ryear + std::floor(((double) (ts_id - 1)) / 12);
              month = (ts_id % 12);
              if (month == 0) month = 12;
              day = 14; /* Is set to mid point by CMOR */
            }
          else if (strcmp(frequency, "day") == 0)
            {
              addseconds = ts_id * 24 * 60 * 60 + 60 * 60 * 12;
              juldate = julianDate_add_seconds(ref_date, addseconds);
            }
          else if (strcmp(frequency, "6hr") == 0)
            {
              addseconds = ts_id * 6 * 60 * 60;
              juldate = julianDate_add_seconds(ref_date, addseconds);
            }
          else if (strcmp(frequency, "3hr") == 0)
            {
              addseconds = ts_id * 3 * 60 * 60;
              juldate = julianDate_add_seconds(ref_date, addseconds);
            }
          else if (strcmp(frequency, "1hr") == 0)
            {
              addseconds = ts_id * 60 * 60;
              juldate = julianDate_add_seconds(ref_date, addseconds);
            }
          if (addseconds == 0)
            {
              CdiDateTime dt{};
              dt.date = cdiDate_encode(year, month, 1);
              juldate = julianDate_encode(calendar, dt);
            }
        }
    }
  if ( time_axis == 1 && strcmp(kv_get_a_val(kvl, "ta", "n"), "cmip") == 0 )
    {
      auto julDateTime = julianDate_decode(calendar, juldate);
      cdiTime_decode(julDateTime.time, &hour, &min, &sec, &ms);
      if (strcmp(frequency, "6hr") == 0)
        {
          int iv[]= {0,6,12,18,24};
          size_t ivsize = sizeof(iv)/sizeof(iv[0]);
          int minid = 0;
          julDateTime.time = cdiTime_encode(0, 0, 0, 0);
          double diff = julianDate_to_seconds(julianDate_sub(juldate, julianDate_encode(calendar, julDateTime)));

          for ( size_t loopid = 1; loopid < ivsize; loopid++ )
            {
              julDateTime.time = cdiTime_encode(iv[loopid], 0, 0, 0);
              if ( std::fabs(julianDate_to_seconds(julianDate_sub(juldate, julianDate_encode(calendar, julDateTime)))) < diff )
                {
                  minid = loopid;
                  julDateTime.time = cdiTime_encode(iv[minid], 0, 0, 0);
                  diff = std::fabs(julianDate_to_seconds(julianDate_sub(juldate, julianDate_encode(calendar, julDateTime))));
                }
            }

          vDateTime.time = cdiTime_encode(iv[minid], 0, 0, 0);
          juldate = julianDate_encode(calendar, vDateTime);
        }
      else if (strcmp(frequency, "3hr") == 0)
        {
          int iv[]= {0,3,6,9,12,15,18,21,24};
          size_t ivsize = sizeof(iv)/sizeof(iv[0]);
          int minid = 0;
          julDateTime.time = cdiTime_encode(0, 0, 0, 0);
          double diff = julianDate_to_seconds(julianDate_sub(juldate, julianDate_encode(calendar, julDateTime)));

          for ( size_t loopid = 1; loopid < ivsize; loopid++ )
            {
              julDateTime.time = cdiTime_encode(iv[loopid], 0, 0, 0);
              if ( std::fabs(julianDate_to_seconds(julianDate_sub(juldate, julianDate_encode(calendar, julDateTime)))) < diff )
                {
                  minid = loopid;
                  julDateTime.time = cdiTime_encode(iv[minid], 0, 0, 0);
                  diff = std::fabs(julianDate_to_seconds(julianDate_sub(juldate, julianDate_encode(calendar, julDateTime))));
                }
            }

          vDateTime.time = cdiTime_encode(iv[minid], 0, 0, 0);
          juldate = julianDate_encode(calendar, vDateTime);
        }
      else if (strcmp(frequency, "1hr") == 0)
        {
          size_t ivsize = 25;
          int minid = 0;
          julDateTime.time = cdiTime_encode(0, 0, 0, 0);
          double diff = julianDate_to_seconds(julianDate_sub(juldate, julianDate_encode(calendar, julDateTime)));

          for ( size_t loopid = 1; loopid < ivsize; loopid++ )
            {
              julDateTime.time = cdiTime_encode(loopid, 0, 0, 0);
              if ( std::fabs(julianDate_to_seconds(julianDate_sub(juldate, julianDate_encode(calendar, julDateTime)))) < diff )
                {
                  minid = loopid;
                  julDateTime.time = cdiTime_encode(minid, 0, 0, 0);
                  diff = std::fabs(julianDate_to_seconds(julianDate_sub(juldate, julianDate_encode(calendar, julDateTime))));
                }
            }

          julDateTime.time = cdiTime_encode(minid, 0, 0, 0);
          juldate = julianDate_encode(calendar, julDateTime);
        }
    }
  return juldate;
}

static double *
get_time_bounds(KVList *kvl, int taxisID, int ifreq, JulianDate ref_date, JulianDate jtime_val, int calendar,
                int tunitsec, double *time_bnds, int time_axis, int /*vlistID*/)
{
  double time_val = julianDate_to_seconds(julianDate_sub(jtime_val, ref_date)) / tunitsec;
  const auto vDateTime = taxisInqVdatetime(taxisID);
  auto vDateTimeCorr = vDateTime;
  CdiDateTime vDateTime0b{};
  CdiDateTime vDateTime1b{};
  int year, month, day;
  int hour, min, sec, ms;
  cdiDate_decode(vDateTime.date, &year, &month, &day);
  if (month == 0 || day == 0)
    {
      vDateTimeCorr = julianDate_decode(calendar, jtime_val);
      cdiDate_decode(vDateTimeCorr.date, &year, &month, &day);
    }
  /***/
  /* If file time axis has bounds use them, otherwise use cmor time axis deduced from miptable frequency and
   * cell_methods or frequency itself*/
  /***/

  if (!taxisHasBounds(taxisID) || strcmp(kv_get_a_val(kvl, "ta", "n"), "cmip") == 0
      || time_axis == 2 || ifreq == 8 )
    {
      /***/
      /* Climatologies */
      /***/

      if (time_axis == 2)
        {
          if (Options::cdoVerbose) cdo_print("10.4.1. Get climatology interval.");
          int numdates = 0;
          std::vector<std::string> climyears = kv_get_vals(kvl, "ci", &numdates);
          if (!climyears.size())
            cdo_abort("ERROR (infile: '%s')! In writing model output:\n          Could not calculate time bounds for climatology time "
                     "axis because attribute 'climatology_interval' is not available.", cdo_get_stream_name(0));
          if (numdates != 2)
            cdo_abort("ERROR (infile: '%s')! In writing model output:\n          Could not calculate time bounds for climatology time "
                     "axis because attribute 'climatology_interval' has not two values.", cdo_get_stream_name(0));
          int expstartyear = atol(climyears[0].c_str());
          int expendyear = atol(climyears[1].c_str());

          vDateTime0b.date = cdiDate_encode(expstartyear, month, 1);
          month++;
          if (month != 12)
            vDateTime1b.date = cdiDate_encode(expendyear, month, 1);
          else
            vDateTime1b.date = cdiDate_encode(expendyear + 1, 1, 1);
        }
      /***/
      /* Diurnal cycle */
      /***/
      else if (time_axis == 3)
        {
          vDateTime0b.date = cdiDate_encode(year, month, 1);
          cdiTime_decode(vDateTime.time, &hour, &min, &sec, &ms);
          vDateTime0b.time = cdiTime_encode(hour, 0, 0, 0);

          hour++;
          month++;
          vDateTime1b.time = cdiTime_encode(hour, 0, 0, 0);
          vDateTime1b.date = (hour > 23) ? cdiDate_encode(year, month, 2) : cdiDate_encode(year, month, 1);
        }
      else
        {
          /***/
          /* Frequency dependent: */
          /***/
          if ( ifreq == 1)
            {
              vDateTime0b.date = cdiDate_encode(year, 1, 1);
              vDateTime1b.date = cdiDate_encode(year + 1, 1, 1);
            }
          else if ( ifreq == 7)
            {
              if ( month == 12 ||  month == 1 ||  month == 2 )
                {
                  if ( month == 12 )
                    {
                      vDateTime0b.date = cdiDate_encode(year, 12, 1);
                      vDateTime1b.date = cdiDate_encode(year+1, 3, 1);
                    }
                  else
                    {
                      vDateTime0b.date = cdiDate_encode(year-1, 12, 1);
                      vDateTime1b.date = cdiDate_encode(year, 3, 1);
                    }
                }
              else if ( month == 3 ||  month == 4 ||  month == 5 )
                {
                  vDateTime0b.date = cdiDate_encode(year, 3, 1);
                  vDateTime1b.date = cdiDate_encode(year, 6, 1);
                }
              else if ( month == 6 ||  month == 7 ||  month == 8 )
                {
                  vDateTime0b.date = cdiDate_encode(year, 6, 1);
                  vDateTime1b.date = cdiDate_encode(year, 9, 1);
                }
              else
                {
                  vDateTime0b.date = cdiDate_encode(year, 9, 1);
                  vDateTime1b.date = cdiDate_encode(year, 12, 1);
                }
            }
          else if ( ifreq == 2)
            {
              vDateTime0b.date = cdiDate_encode(year, month, 1);
              month++;
              if (month > 12)
                {
                  month = 1;
                  year++;
                }
              vDateTime1b.date = cdiDate_encode(year, month, 1);
            }
          else if ( ifreq == 3)
            {
/*Assuming that the requested axis is always "days since.. 00:00:00" */
              time_bnds[0] = std::floor(time_val);
              time_bnds[1] = std::ceil(time_val);
/*Assuming that the averaged value is written for and at the end of the interval */
              if ( DBL_IS_EQUAL(time_bnds[0], time_bnds[1]) )
                time_bnds[0] -= 1.;
              return time_bnds;
            }
          /***/
          /* Note that time_val must be correct in Infile for subdaily frequencies */
          /***/
          else if ( ifreq == 4)
            {
              cdiTime_decode(vDateTime.time, &hour, &min, &sec, &ms);
              int iv[]= {0,6,12,18,24};
              size_t ivsize = sizeof(iv)/sizeof(iv[0]);
              int minid = 0, newhour[2];
              vDateTimeCorr.time = cdiTime_encode(0, 0, 0, 0);
              double diff = julianDate_to_seconds(julianDate_sub(jtime_val, julianDate_encode(calendar, vDateTimeCorr)));
              for ( size_t loopid = 1; loopid < ivsize; loopid++ )
                {
                  vDateTimeCorr.time = cdiTime_encode(iv[loopid], 0, 0, 0);
                  if ( std::fabs(julianDate_to_seconds(julianDate_sub(jtime_val, julianDate_encode(calendar, vDateTimeCorr)))) < diff )
                    {
                      minid = loopid;
                      vDateTimeCorr.time = cdiTime_encode(iv[minid], 0, 0, 0);
                      diff = std::fabs(julianDate_to_seconds(julianDate_sub(jtime_val, julianDate_encode(calendar, vDateTimeCorr))));
                    }
                }
              newhour[0] = ( hour-iv[minid] < 0 ) ? iv[minid-1] : iv[minid];
              newhour[1] = ( hour-iv[minid] < 0 ) ? iv[minid] : iv[minid+1];

              vDateTime0b.time = cdiTime_encode(newhour[0], 0 , 0, 0);
              vDateTime1b.time = cdiTime_encode(newhour[1], 0 , 0, 0);

              vDateTime0b.date = cdiDate_encode(year, month, day);
              vDateTime1b.date = vDateTime0b.date;
            }
          else if ( ifreq == 5)
            {
              cdiTime_decode(vDateTime.time, &hour, &min, &sec, &ms);
              int iv[]= {0,3,6,9,12,15,18,21,24};
              size_t ivsize = sizeof(iv)/sizeof(iv[0]);
              int minid = 0, newhour[2];
              vDateTimeCorr.time = cdiTime_encode(0, 0, 0, 0);
              double diff = julianDate_to_seconds(julianDate_sub(jtime_val, julianDate_encode(calendar, vDateTimeCorr)));
              for ( size_t loopid = 1; loopid < ivsize; loopid++ )
                {
                  vDateTimeCorr.time = cdiTime_encode(iv[loopid], 0, 0, 0);
                  if ( std::fabs(julianDate_to_seconds(julianDate_sub(jtime_val, julianDate_encode(calendar, vDateTimeCorr)))) < diff )
                    {
                      minid = loopid;
                      vDateTimeCorr.time = cdiTime_encode(iv[minid], 0, 0, 0);
                      diff = std::fabs(julianDate_to_seconds(julianDate_sub(jtime_val, julianDate_encode(calendar, vDateTimeCorr))));
                    }
                }
              newhour[0] = ( hour-iv[minid] < 0 ) ? iv[minid-1] : iv[minid];
              newhour[1] = ( hour-iv[minid] < 0 ) ? iv[minid] : iv[minid+1];

              vDateTime0b.time = cdiTime_encode(newhour[0], 0 , 0, 0);
              vDateTime1b.time = cdiTime_encode(newhour[1], 0 , 0, 0);

              vDateTime0b.date = cdiDate_encode(year, month, day);
              vDateTime1b.date = vDateTime0b.date;
            }
          else if ( ifreq == 6)
            {
              cdiTime_decode(vDateTime.time, &hour, &min, &sec, &ms);
              vDateTime0b.time = cdiTime_encode(hour, 0 , 0, 0);
              vDateTime1b.time = cdiTime_encode(hour+1, 0 , 0, 0);

              vDateTime0b.date = cdiDate_encode(year, month, day);
              vDateTime1b.date = vDateTime0b.date;
            }
          else if ( ifreq == 8)
            {
              if (Options::cdoVerbose) cdo_print("10.4.1. Get decadal interval.");
              int numdates = 0;
              std::vector<std::string> decyears = kv_get_vals(kvl, "di", &numdates);
              if (!decyears.size())
                cdo_abort("ERROR (infile: '%s')! In writing model output:\n          Could not calculate time bounds for decadal time "
                     "axis because attribute 'decadal_interval' is not available.", cdo_get_stream_name(0));
              if (numdates != 2)
                cdo_abort("ERROR (infile: '%s')! In writing model output:\n          Could not calculate time bounds for decadal time "
                     "axis because attribute 'decadal_interval' has not two values.", cdo_get_stream_name(0));
              int expstartyear = atol(decyears[0].c_str());
              int expendyear = atol(decyears[1].c_str());
              
              vDateTime0b.date = cdiDate_encode(expstartyear, 1, 1);
              vDateTime1b.date = cdiDate_encode(expendyear, 1, 1);
              
            }
        }
    }
  else
    {
      taxisInqVdatetimeBounds(taxisID, &vDateTime0b, &vDateTime1b);
    }

  auto juldate = julianDate_encode(calendar, vDateTime0b);
  time_bnds[0] = julianDate_to_seconds(julianDate_sub(juldate, ref_date)) / tunitsec;

  juldate = julianDate_encode(calendar, vDateTime1b);
  time_bnds[1] = julianDate_to_seconds(julianDate_sub(juldate, ref_date)) / tunitsec;
  if (time_axis == 3) time_bnds[1] -= 1;

  return time_bnds;
}

static void
read_record(CdoStreamID streamID, struct mapping vars[], int vlistID)
{
  int varID, levelID;
  cdo_inq_record(streamID, &varID, &levelID);

  int gridID = vlistInqVarGrid(vlistID, varID);
  int type = gridInqType(gridID);
  auto gridsize = gridInqSize(gridID);
  double *buffer = (double *) Malloc(gridsize * sizeof(double));

  struct mapping *var = map_var(varID, vars);
  if (var && var->charvars != 1)
    {
      int zaxisID = vlistInqVarZaxis(vlistID, varID);
      int ztype = zaxisInqType(zaxisID) ;
/*      int latdim = gridInqYsize(gridID); */
      int levdim = zaxisInqSize(zaxisID);
      size_t nmiss;
      cdo_read_record(streamID, buffer, &nmiss);
      for (size_t i = 0; i < gridsize; ++i)
        {
          // Wrong:  (lat x basin, lev ) gridsize * levelID + i
          // Wrong:  (basin x lat, lev) gridsize * levelID + i * chardim - ( int ) std::floor(i / latdim) * gridsize + ( int
          // ) std::floor(i/latdim)
          // Wrong:  (basin x lev, lat ) gridsize/latdim * levdim * ( i - ( int ) std::floor(i/latdim) * latdim ) + ( int )
          // std::floor(i/latdim) + gridsize/latdim * levelID;
          // Wrong:  (lat x lev, basin ) latdim * levdim * ( int ) std::floor(i/latdim) + ( i - ( int ) std::floor(i/latdim) *
          // latdim ) + levelID * latdim
          // (lev x lat, basin )
          int newIndex;
          if (levdim > 1 )
            {
              if ( (type == GRID_UNSTRUCTURED || type == GRID_CURVILINEAR ) && ztype != ZAXIS_HYBRID)
                newIndex = i + gridsize * levelID; 
//              else if ( type == GRID_CURVILINEAR )
//                newIndex = i + gridsize * levelID;
              else
                newIndex = i * levdim + levelID;
            }
          else
            newIndex = i;
          if (var->datatype == 'f')
            {
              ((float *) var->data)[newIndex] = (float) buffer[i];
            }
          else
            {
              ((double *) var->data)[newIndex] = (double) buffer[i];
            }
        }
    }
  Free(buffer);
}

static int
check_append_and_size(KVList *kvl, int /*vlistID*/, char *testIn, int ifreq, int calendar)
{
  char *test = testIn;
  size_t filesize = FileUtils::size((const char *) testIn);
  char old_start_date[CMOR_MAX_STRING];
  char old_end_date[CMOR_MAX_STRING];
  int i = 0, j = 0;
  /* Get dates from chunk string */
  if (Options::cdoVerbose) cdo_print("Start to retrieve dates from chunk string.");
  while (*(test + i) != 0)
    {
      if (*(test + i) == '_')
        {
          test += (i + 1);
          i = 0;
        }
      if (*(test + i) == '-') j = i;
      i++;
    }
  if (!i || !j || *(test + j + 1) == 0 || *(test + 2 * j) == 0)
    {
      if (Options::cdoVerbose) cdo_print("In checking the last chunk:\n          Date from filename of the chunk cannot be read.\n          "
                 "Switched to replace mode for this variable.");
      return 0;
    }

  strncpy(old_start_date, test, j);
  old_start_date[j] = 0;
  test += (j + 1);
  strncpy(old_end_date, test, j);
  old_end_date[j] = 0;

  if (Options::cdoVerbose)
    cdo_print("Successfully retrieved start date: '%s' and end date: '%s' chunk string.", old_start_date,
             old_end_date);
  /* Check frequency of chunk with frequency of file */

  if ((j == 12 && ifreq < 4) ||
      (j == 8 && ifreq != 3) ||
      (j == 6 && ifreq != 2) ||
      (j == 4 && ifreq != 1 && ifreq != 8) )
    {
      if (Options::cdoVerbose) cdo_print("In checking last chunk:\n          Frequency of chunk file does not agree with frequency of the "
                 "working file.\n         Switched to replace mode for this variable.");
      return 0;
    }

  /* Encode in julseconds depending on frequency */
  if (Options::cdoVerbose) cdo_print("Start to encode dates with frequencies to julseconds.");

  int old_start_year = 0, old_start_month = 1, old_start_day = 1, old_start_hr = 1, old_start_min = 1;
  int old_end_year = 0, old_end_month = 1, old_end_day = 1, old_end_hr = 1, old_end_min = 1;

  switch (j)
    {
    case (12):
      std::sscanf(old_start_date, "%04d%02d%02d%02d%02d", &old_start_year, &old_start_month, &old_start_day, &old_start_hr, &old_start_min);
      std::sscanf(old_end_date, "%04d%02d%02d%02d%02d", &old_end_year, &old_end_month, &old_end_day, &old_end_hr, &old_end_min);
      break;
    case (8):
      std::sscanf(old_start_date, "%04d%02d%02d", &old_start_year, &old_start_month, &old_start_day);
      std::sscanf(old_end_date, "%04d%02d%02d", &old_end_year, &old_end_month, &old_end_day);
      break;
    case (6):
      std::sscanf(old_start_date, "%04d%02d", &old_start_year, &old_start_month);
      std::sscanf(old_end_date, "%04d%02d", &old_end_year, &old_end_month);
      break;
    case (4):
      old_start_year = atol(old_start_date);
      old_end_year = atol(old_end_date);
      break;
    default:
      {
        if ( strcmp(kv_get_a_val(kvl, "om", "r"),"a") == 0 )
          {
            cdo_warning("Last chunk's frequency cannot yet be tested for being suitable.");
            return 1;
          }
        else if (Options::cdoVerbose)
          {
            cdo_print("In checking last chunk:\n          Last chunk has unknown frequency "
                   "which is yet not enabled to be appended by "
                   "cdo cmor.\n          Switched to replace mode for this variable.");
            return 0;
          }
      }
    }

  CdiDateTime startDateTime{};
  CdiDateTime endDateTime{};
  startDateTime.date = cdiDate_encode(old_start_year, old_start_month, old_start_day);
  endDateTime.date = cdiDate_encode(old_end_year, old_end_month, old_end_day);

  startDateTime.time = cdiTime_encode(old_start_hr, old_start_min, 0, 0);
  endDateTime.time = cdiTime_encode(old_end_hr, old_end_min, 0, 0);
  auto julostart = julianDate_encode(calendar, startDateTime);
  auto juloend = julianDate_encode(calendar, endDateTime);

  if (Options::cdoVerbose) cdo_print("Successfully calculated juldates.");
  /* Read in first vdate in case not piped */
  if (Options::cdoVerbose) cdo_print("Start to calculate temporal gap between chunk and working file.");
  if (cdo_assert_files_only() == false)
    {
      if (Options::cdoVerbose) cdo_print("Cdo cmor cannot enable append mode since you piped several cdo operators.\n          Switched to "
                 "replace mode for this variable.");
      if ( strcmp(kv_get_a_val(kvl, "om", "r"),"a") == 0 )
        {
          cdo_warning("Could not check whether chunk is suitable to be appended since you piped operators.\n"
                     "          Note that the operator did not check 1. for time gaps and 2. for the max size.");
          cdo_print("Output mode: (A)pend.");
          return 1;
        }
      return 0;
    }

  CdiStreamID streamID2 = streamOpenRead(cdo_get_stream_name(0));
  int vlistID2 = streamInqVlist(streamID2);
  int taxisID2 = vlistInqTaxis(vlistID2);
  const auto vDateTime2 = taxisInqVdatetime(taxisID2);
  auto firstdate = julianDate_encode(calendar, vDateTime2);

  int fyear, fmonth, dummy;
  cdiDate_decode(vDateTime2.date, &fyear, &fmonth, &dummy);
  if ( ifreq == 1 && fyear == old_end_year )
    {
      if (Options::cdoVerbose) cdo_print("In checking the last chunk:\n          The years of the end date of the chunk file "
                 "and the first date of the working file are the same: '%d'."
                 "   Switched to replace mode for this variable.", fyear);
      return 0;
    }

  /* Check temporal distance between last chunk date and first file date */
  double append_distance = julianDate_to_seconds(julianDate_sub(firstdate, juloend)) / 3600.0;
  if (   (ifreq == 6 && (append_distance > 2.0 || append_distance < 0))
      || (ifreq == 5 && (append_distance > 6.0 || append_distance < 0))
      || (ifreq == 4 && (append_distance > 12.0 || append_distance < 0))
      || (ifreq == 3 && (append_distance > 48.0 || append_distance < 0))
      || (ifreq == 2 && (append_distance / 24.0 > 62.0 || append_distance < 0))
      || (ifreq == 1 && (append_distance / 24.0 / 30.5 > 24.0 || append_distance < 0))
      || (ifreq == 1 && (append_distance / 24.0 / 30.5 < 1.0 ))
      || (ifreq == 8 && (append_distance / 24.0 / 30.5 / 12.0 < 1.0 || append_distance < 0))
      || (ifreq == 8 && (append_distance / 24.0 / 30.5 / 12.0 > 20.0 ))
)
    {
      if (Options::cdoVerbose) cdo_print("In checking the last chunk:\n          A temporal gap is diagnosed between end date of chunk file "
                 "and first date of working file of: '%f' hours ( '%f' days, '%f' months, '%f' years)"
                 ". Maximal valid gaps are:\n"
                 "          2 hours for 1-hourly frequency\n"
                 "          6 hours for 3-hourly frequency\n"
                 "          12 hours for 6-hourly frequency\n"
                 "          48 hours for daily frequency\n"
                 "          62 days for monthly frequency\n"
                 "          24 months for yearly frequency\n"
                 "          20 years for decadal frequency\n"
                 "Minimal valid gaps are:\n"
                 "          1 month for yearly frequency\n"
                 "          1 year for decadal frequency\n"
                 "          Switched to replace mode for this variable.",
                 append_distance, append_distance/24.0, append_distance/24.0/30.5, append_distance/24.0/30.5/12.0);
      streamClose(streamID2);
      return 0;
    }
  else if ( Options::cdoVerbose )
    cdo_print("The temporal gap between end date of chunk file and first date of working file is '%f' and therefore valid.", append_distance);

  if (Options::cdoVerbose) cdo_print("Successfully checked temporal gap.");
  /* Check file size */
  if (Options::cdoVerbose) cdo_print("Start to check file size of chunk + working file.");
  double old_interval_sec = julianDate_to_seconds(julianDate_sub(juloend, julostart));
  double size_per_sec = (double) filesize / old_interval_sec;

  int maxsizegb = atol(kv_get_a_val(kvl, "ms", "2"));
  int maxsizeb = maxsizegb * 1024 * 1024 * 1024;

  int ntsteps = vlistNtsteps(vlistID2);
  if (ntsteps < 0)
    {
      ntsteps = 0;
      while (streamInqTimestep(streamID2, ntsteps++))
        ;
      if (ntsteps == 0)
        {
          if (Options::cdoVerbose) cdo_print("In checking whether append mode is possible:\n          No time steps found in infile.\n         "
                     " Switched to replace mode for this variable.");
          streamClose(streamID2);
          return 0;
        }
    }

  double estimated_size;
  switch (ifreq)
    {
    case (6): estimated_size = ntsteps * 60 * 60 * 1 * size_per_sec + (double) filesize; break;
    case (5): estimated_size = ntsteps * 60 * 60 * 3 * size_per_sec + (double) filesize; break;
    case (4): estimated_size = ntsteps * 60 * 60 * 6 * size_per_sec + (double) filesize; break;
    case (3): estimated_size = ntsteps * 60 * 60 * 24 * size_per_sec + (double) filesize; break;
    case (2): estimated_size = ntsteps * 60 * 60 * 24 * 30.5 * size_per_sec + (double) filesize; break;
    case (1): estimated_size = ntsteps * 60 * 60 * 24 * 365.25 * size_per_sec + (double) filesize; break;
    case (8): estimated_size = ntsteps * 10 * 60 * 60 * 24 * 365.25 * size_per_sec + (double) filesize; break;
    default:
      {
        if (Options::cdoVerbose) cdo_print("In checking whether append mode is valid:\n          Selected chunk to append data has subdaily frequency which is yet not enabled by cdo cmor.\n          Switched to replace mode for this variable.");
        streamClose(streamID2);
        return 0;
      }
    }

  if ( maxsizeb != 0 && (unsigned int) estimated_size > (unsigned int) maxsizeb)
    {
      if (Options::cdoVerbose) cdo_print("In checking whether append mode is valid:\n          Estimated file size of appended file is : '%f'gb and exceeds maximal allowed file size: '%d'gb.\n          Switched to replace mode for this variable.",
                 estimated_size / 1024.0 / 1024.0 / 1024.0, maxsizegb);
      streamClose(streamID2);
      return 0;
    }
  streamClose(streamID2);
  if (Options::cdoVerbose) cdo_print("Successfully checked file size of chunk + working file.");
  return 1;
}

static char *
use_chunk_des_files(KVList *kvl, int vlistID, int /*var_id*/, char *chunk_des_file, int ifreq, int calendar)
{
  char *chunk_file = (char *) Malloc(4096 * sizeof(char));
  if (file_exist(chunk_des_file, false, "chunk_description", false))
    {
      auto *fp = std::fopen(chunk_des_file, "r");
      if (fp == nullptr)
        cdo_abort("Could not open chunk description file '%s'", chunk_des_file);
      ListBuffer listBuffer;
      auto status = listBuffer.read(fp, chunk_des_file);
      if (status) cdo_abort("Read error on chunk_description %s!", chunk_des_file);

      char *cbuffer = (char *)(listBuffer.buffer.data());
      size_t pos = 0;
      for ( pos = 0; pos < listBuffer.buffer.size(); pos++ )
        if ( cbuffer[pos] == '\n' )
          {
            cbuffer[pos] = 0;
            break;
          }
      if ( pos == listBuffer.buffer.size() )
        {
          cbuffer[pos-1] = 0;
          pos--;
        }
      std::snprintf(chunk_file, pos+1, "%s", cbuffer);

      if (file_exist(chunk_file, false, "chunk_description", false) )
        {
          if ( check_append_and_size(kvl, vlistID, chunk_file, ifreq, calendar) )
            return chunk_file;
          else
            {
              if (Options::cdoVerbose) cdo_print("In checking the last chunk:\n          Chunk '%s' configured via chunk description file is not suitable to be appended.",
                   chunk_file);
              cdo_print("Output mode: (R)eplace.");
            }
        }
      else
        {
          if ( chunk_file[0] )
            {
              if (Options::cdoVerbose) cdo_print("In checking the last chunk:\n          Chunk '%s' configured via chunk description file does not exist.", chunk_file);
              cdo_print("Output mode: (R)eplace.");
            }
          else
            {
              if (Options::cdoVerbose) cdo_print("In checking the last chunk:\n          No name found in chunk description file.");
              cdo_print("Output mode: (R)eplace.");
            }
        }
    }
  else
    {
      if (Options::cdoVerbose) cdo_print("Chunk description file '%s' could not be opened.", chunk_des_file);
    }
  strcpy(chunk_file, " \0");
  return chunk_file;
}

static char **
empty_array(struct mapping vars[], char ***chunk_files)
{
  for (int i = 0; vars[i].cmor_varID != CMOR_UNDEFID; ++i) (*chunk_files)[i] = nullptr;
  return *chunk_files;
}

static char **
get_chunk_des_files(KVList *kvl, struct mapping vars[], char *miptab_freqptr, int nreq, int vlistID, char *charname, char *project_id)
{
  char **chunk_des_files = (char **) Malloc((nreq + 1) * sizeof(char *));
  chunk_des_files[nreq] = nullptr;

  char trunk[CMOR_MAX_STRING];
  if (strcmp(project_id, "CMIP6") == 0)
    sprintf(trunk,"%s_", kv_get_a_val(kvl, "source_id", ""));
  else
    sprintf(trunk,"%s_", kv_get_a_val(kvl, "model_id", ""));
  const char *description_atts[] = { "experiment_id", "member", "sub_experiment_id", nullptr };
  strcpy(trunk, miptab_freqptr);
  for (int i = 0; description_atts[i]; ++i)
    {
      strcat(trunk, "_");
      strcat(trunk, kv_get_a_val(kvl, description_atts[i], ""));
    }

  for (int j = 0; vars[j].cmor_varID != CMOR_UNDEFID; ++j)
    {
      char *name = (char *) Malloc(CDI_MAX_NAME * sizeof(char));
      if (charname)
        strcpy(name, charname);
      else
        vlistInqVarName(vlistID, vars[j].cdi_varID, name);
      chunk_des_files[j] = (char *) Malloc(CMOR_MAX_STRING * sizeof(char));
      sprintf(chunk_des_files[j], ".CHUNK_FILE_%s_%s.txt", name, trunk);
      Free(name);
    }
  return chunk_des_files;
}

static char **
get_chunk_files(KVList *kvl, struct mapping vars[], int vlistID, int ifreq, int time_axis, int calendar,
                char *miptab_freqptr, char *project_id, int *mergeIDs, int psID)
{
  int i = 0;
  for (i = 0; vars[i].cmor_varID != CMOR_UNDEFID; ++i) ;
  if ( mergeIDs[0] != CMOR_UNDEFID )
    i = 1;
  char **chunk_files = (char **) Malloc((i + 1) * sizeof(char *));
  chunk_files[i] = nullptr;

  if (Options::cdoVerbose) cdo_print("10.2.1. Start to validate append mode.");
  char *dummy = kv_get_a_val(kvl, "om", "a");
  if (strcmp(dummy, "a") != 0)
    return empty_array(vars, &chunk_files);
  else if (time_axis == 4)
    {
      if (Options::cdoVerbose) cdo_print("In validating append mode:\n          CMOR APPEND mode not possible for time independent "
                 "variables.");
      cdo_print("Output mode: (R)eplace.");
      return empty_array(vars, &chunk_files);
    }
  if (Options::cdoVerbose) cdo_print("10.2.1. Successfully validated append mode.");


  if (Options::cdoVerbose) cdo_print("10.2.2. Start to get chunk names.");
  int num_aaf = 0;
  std::vector<std::string> chunk_att_files = kv_get_vals(kvl, "lc", &num_aaf);
  char **chunk_des_files = nullptr;
  if (num_aaf != i && num_aaf > 0)
    {
      if (Options::cdoVerbose) cdo_print(
          "Number of chunk files '%d' disagree with number of requested variables '%d'.\n Switched to replace mode.\n",
          num_aaf, i);
      cdo_print("Output mode: (R)eplace.");
      return empty_array(vars, &chunk_files);
    }
  else if (num_aaf == 0)
    {
      char *nd = kv_get_a_val(kvl, "d", "y");
      /* For chunk description file : */
      if (nd[0] == 'y')
        chunk_des_files = get_chunk_des_files(kvl, vars, miptab_freqptr, i, vlistID, nullptr, project_id);
      else
        {
          if (Options::cdoVerbose) cdo_print("In getting chunk names:\n          Automatic chunk configuration via file not possible if DRS is "
                     "not created.");
          cdo_print("Output mode: (R)eplace.");
          return empty_array(vars, &chunk_files);
        }
    }
  if (Options::cdoVerbose) cdo_print("10.2.2. Successfully retrieved chunk names.");


  for (int j = 0; vars[j].cmor_varID != CMOR_UNDEFID; ++j)
    {
      if (vars[j].cmor_varID == psID and vars[j].help_var)
        {
          if (Options::cdoVerbose)
            cdo_print("Chunkfile for ps which was registered for hybrid axis is skipped.");
          continue ;
        }
      if (num_aaf != 0)
        {
          if (file_exist(chunk_att_files[j].c_str(), false, "chunk file", false)
              && check_append_and_size(kvl, vlistID, (char *)chunk_att_files[j].c_str(), ifreq, calendar))
            chunk_files[j] = strdup(chunk_att_files[j].c_str());
          else
            {
              if (Options::cdoVerbose) cdo_print("Chunk '%s' could not be used.", chunk_att_files[j]);
              cdo_print("Output mode: (R)eplace.");
              chunk_files[j] = strdup(" ");
            }
        }
      else
        {
          if (Options::cdoVerbose)
            cdo_print("It is tried to open a chunk description file for varID: '%d': '%s'.", vars[j].cdi_varID,
                     chunk_des_files[j]);
          chunk_files[j] = use_chunk_des_files(kvl, vlistID, vars[j].cdi_varID, chunk_des_files[j], ifreq, calendar);
          printf("%d %s\n", j, chunk_files[j]);
        }
      if (strcmp(chunk_files[j], " ") != 0)
        cdo_print("Output mode: (A)ppend.\n          (Chunk file for var ID %d is: '%s')", vars[j].cdi_varID, chunk_files[j]);
    }
  if (chunk_des_files) free_array(chunk_des_files);
  if (Options::cdoVerbose) cdo_print("Successfully processed chunk file retrieval.");
  return chunk_files;
}

static void
write_variables(KVList *kvl, CdoStreamID streamID, struct mapping vars[], int miptab_freq, int time_axis, int calendar,
                char *miptab_freqptr, char *project_id, int *mergeIDs)
{
  int cmf = 0;
  int vlistID = cdo_stream_inq_vlist(streamID);
  int taxisID = vlistInqTaxis(vlistID);
  int tsID = 0;
  int nrecs;
  size_t gridsize = vlistGridsizeMax(vlistID);

  if (Options::cdoVerbose) cdo_print("10. Start to write variables via cmor_write.");
  if (Options::cdoVerbose) cdo_print("10.1. Start to get frequency.");
  int time_unit;
  const auto sDateTime = get_taxis(kv_get_a_val(kvl, "rtu", nullptr), &time_unit);
  int tunitsec = get_tunitsec(time_unit);
  const auto ref_date = julianDate_encode(calendar, sDateTime);
  char *frequency = nullptr;
  if (time_axis != 4)
    {
      frequency = get_frequency(/*kvl,*/ vlistID, miptab_freq);
      if (Options::cdoVerbose && strcmp(frequency, "no") != 0 ) cdo_print("10.1. Successfully retrieved frequency %s.", frequency);
    }
  else if (Options::cdoVerbose ) cdo_print("10.1. Successfully retrieved fixed time axis.");

  int ifreq = 0;
  if (frequency)
    {
      if (strcmp(frequency, "yr") == 0) ifreq = 1;
      if (strcmp(frequency, "mon") == 0) ifreq = 2;
      if (strcmp(frequency, "day") == 0) ifreq = 3;
      if (strstr(frequency, "6hr") ) ifreq = 4;
      if (strstr(frequency, "3hr") ) ifreq = 5;
      if (strstr(frequency, "1hr") ) ifreq = 6;
      if (strstr(frequency, "sem") ) ifreq = 7;
      if (strstr(frequency, "dec") ) ifreq = 8;
      if (strstr(frequency, "clim") ) ifreq = 9;
      if (strstr(frequency, "subhr") ) ifreq = 10;
    }

  int ps_index = getVarIDToMap(vlistID, vlistNvars(vlistID), "name", "ps");
  int psID = getRegisteredPsid(vars, ps_index);

  if (Options::cdoVerbose) cdo_print("10.2. Start to get chunk files.");
  char **chunk_files = nullptr;
  if ( ifreq > 0 && ifreq != 7 )
    chunk_files = get_chunk_files(kvl, vars, vlistID, ifreq, time_axis, calendar, miptab_freqptr, project_id, mergeIDs, psID);
  else
    {
      int number = 0;
      for (number = 0; vars[number].cmor_varID != CMOR_UNDEFID; number++)
      chunk_files = (char **) Malloc((number + 1) * sizeof(char *));
      empty_array(vars, &chunk_files);
    }
  if ( ifreq == 7 ) cdo_print("10.2. Append mode not possible for frequency '%s'. Switch to replace mode.", frequency);
  if (Options::cdoVerbose) cdo_print("10.2. Successfully retrieved chunk files.");
  int i = 0;
  if ( chunk_files[0] && chunk_files[0][0] != ' ' && strcmp(kv_get_a_val(kvl, "sc", "n"), "y") == 0 )
    {
      while ( chunk_files[i] )
        {
          char command[CDI_MAX_NAME];
          sprintf(command, "cp %s %s.save", chunk_files[i], chunk_files[i]);
          int dir_err = system(command);
          if (dir_err != 0)
            cdo_warning("Could not create a .save file out of the previous chunk '%s'.", chunk_files[i]);
          i++;
        }
    }
  i = 0;

  int zaxisID, zsize = 0, pscheck = 1;
  char *charname = nullptr;
  CdoStreamID newstreamID = nullptr;
  for (i = 0; vars[i].cdi_varID != CDI_UNDEFID; ++i)
    if (vars[i].charvars)
      {
        if (Options::cdoVerbose) cdo_print("10.3. Start to get auxiliary variables.");
        zaxisID = vlistInqVarZaxis(vlistID, vars[i].cdi_varID);
        zsize = zaxisInqSize(zaxisID);
        charname = (char *) Malloc(CDI_MAX_NAME * sizeof(char));
        vlistInqVarName(vlistID, vars[i].cdi_varID, charname);

        cdo_stream_close(streamID);
        newstreamID = cdo_open_read(0);

        vlistID = cdo_stream_inq_vlist(newstreamID);
        taxisID = vlistInqTaxis(vlistID);

        pscheck = 0;
        if (Options::cdoVerbose) cdo_print("10.3. Successfully retrieved auxiliary variables.");
        break;
      }
  if (pscheck == 0)
    if (Options::cdoVerbose)
      cdo_print("Since you defined a variable with character coordinate axis you cannot write another variable with zaxis "
             "of type ZAXIS_HYBRID.");
  if ( !newstreamID )
    newstreamID = streamID;

  int fsize = 0;
  int *mergeIdx = nullptr;
  if ( mergeIDs[0] != CMOR_UNDEFID )
    {
      while ( mergeIDs[fsize] != CMOR_UNDEFID ) fsize++;;
      if (Options::cdoVerbose)
        cdo_print("10.3. '%d' Variables will be merged.", fsize);
      zaxisID = vlistInqVarZaxis(vlistID, mergeIDs[0]);
      zsize = zaxisInqSize(zaxisID);

      mergeIdx = (int *) Malloc(fsize*sizeof(int));
      mergeIdx[0] = -1;
      for ( int j = 0; j < fsize; j++ )
        {
          for (i = 0; vars[i].cdi_varID != CDI_UNDEFID; ++i)
            {
              if ( vars[i].cdi_varID == mergeIDs[j] )
                mergeIdx[j] = i;
            }
        }
      if ( mergeIdx[0] == -1 )
        cdo_abort("Could not find registered CMOR varID.");
    }

  if (Options::cdoVerbose) cdo_print("10.4. Start to loop over time steps.");
  while ((nrecs = cdo_stream_inq_timestep(newstreamID, tsID++)))
    {
      double time_bnds[2];
      double *time_bndsp = nullptr;
      JulianDate jtime_val;
      double time_val;
      if (time_axis != 4)
        {
          jtime_val = get_cmor_time_val(kvl, taxisID, ref_date, tunitsec, calendar, frequency, tsID, time_axis);
          time_val = julianDate_to_seconds(julianDate_sub(jtime_val, ref_date)) / tunitsec;
          time_bndsp = (time_axis != 1) ? get_time_bounds(kvl, taxisID, ifreq, ref_date, jtime_val, calendar,
                                                          tunitsec, time_bnds, time_axis, vlistID)
                                        : 0;
        }
      while (nrecs--) read_record(newstreamID, vars, vlistID);
      if ( mergeIDs[0] != CMOR_UNDEFID )
        {
          void *dataslice;
          if ( vars[mergeIdx[0]].datatype == 'd' )
            dataslice = (void *) Malloc(gridsize * zsize * fsize * sizeof(double));
          else
            dataslice = (void *) Malloc(gridsize * zsize * fsize * sizeof(float));

          for (i = 0; i<fsize; i++ )
            {
              for (int j = 0; j < (int) gridsize * zsize; ++j)
                 {
                   if ( miptab_freq == 8 )
                     {
                       if ( vars[mergeIdx[0]].datatype == 'd' )
                         ((double *) dataslice)[j*fsize+i] = ((double *) vars[mergeIdx[i]].data)[j];
                       else
                         ((float *) dataslice)[j*fsize+i] = ((float *) vars[mergeIdx[i]].data)[j];
                     }
                   else
                     {
                       if ( vars[mergeIdx[0]].datatype == 'd' )
                         ((double *) dataslice)[j+i*gridsize*zsize] = ((double *) vars[mergeIdx[i]].data)[j];
                       else
                         ((float *) dataslice)[j+i*gridsize*zsize] = ((float *) vars[mergeIdx[i]].data)[j];
                     }
                 }
            }
#if (CMOR_VERSION_MAJOR == 3 && CMOR_VERSION_MINOR <= 2 && CMOR_VERSION_PATCH <= 7)
          cmf = cmor_write(vars[mergeIdx[0]].cmor_varID, dataslice, vars[mergeIdx[0]].datatype, 1, &time_val, time_bndsp, NULL);
#else
          cmf = cmor_write(vars[mergeIdx[0]].cmor_varID, dataslice, vars[mergeIdx[0]].datatype, chunk_files[0], 1, &time_val,
                                       time_bndsp, NULL);
#endif
          Free(dataslice);
        }
      for (i = 0; vars[i].cmor_varID != CMOR_UNDEFID; ++i)
        {
          /*          char name[CDI_MAX_NAME];
                    vlistInqVarName(vlistID, vars[i].cdi_varID, name); */
              if (!vars[i].help_var)
                {
                  if (time_axis != 4)
                    {
                      if (vars[i].charvars)
                        {
                          void *dataslice;
                          if ( vars[i].datatype == 'd' )
                            {
                              dataslice = (void *) Malloc(gridsize * zsize * sizeof(double));
                              for (int j = 0; j < (int) gridsize * zsize; ++j)
                                ((double *) dataslice)[j] = ((double *) vars[i].data)[(tsID - 1) * gridsize * zsize + j];
                            }
                          else
                            {
                              dataslice = (void *) Malloc(gridsize * zsize * sizeof(float));
                              for (int j = 0; j < (int) gridsize * zsize; ++j)
                                ((float *) dataslice)[j] = ((float *) vars[i].data)[(tsID - 1) * gridsize * zsize + j];
                            }
#if (CMOR_VERSION_MAJOR == 3 && CMOR_VERSION_MINOR <= 2 && CMOR_VERSION_PATCH <= 7)
                          cmf = cmor_write(vars[i].cmor_varID, dataslice, vars[i].datatype, 1, &time_val, time_bndsp, nullptr);
#else
                          cmf = cmor_write(vars[i].cmor_varID, dataslice, vars[i].datatype, chunk_files[i], 1, &time_val,
                                       time_bndsp, nullptr);
#endif
                          Free(dataslice);
                        }
                      else if ( vars[i].cdi_varID != mergeIDs[0] )
                        {
#if (CMOR_VERSION_MAJOR == 3 && CMOR_VERSION_MINOR <= 2 && CMOR_VERSION_PATCH <= 7)

                          cmf = cmor_write(vars[i].cmor_varID, vars[i].data, vars[i].datatype, 1, &time_val, time_bndsp,
                                       nullptr);
#else
                          cmf = cmor_write(vars[i].cmor_varID, vars[i].data, vars[i].datatype, chunk_files[i], 1, &time_val,
                                       time_bndsp, nullptr);
#endif
                        }
                      if (vars[i].zfactor_id > 0)
                        {
#if (CMOR_VERSION_MAJOR == 3 && CMOR_VERSION_MINOR <= 2 && CMOR_VERSION_PATCH <= 7)
                          cmf = cmor_write(vars[i].zfactor_id, vars[psID].data, vars[psID].datatype, 1, &time_val,
                                       time_bndsp, &vars[i].cmor_varID);
#else
                          cmf = cmor_write(vars[i].zfactor_id, vars[psID].data, vars[psID].datatype, chunk_files[i],
                                       1, &time_val, time_bndsp, &vars[i].cmor_varID);
#endif
                        }
                    }
                  else
                    {
#if (CMOR_VERSION_MAJOR == 3 && CMOR_VERSION_MINOR <= 2 && CMOR_VERSION_PATCH <= 7)
                      cmf = cmor_write(vars[i].cmor_varID, vars[i].data, vars[i].datatype, 0, 0, 0, nullptr);
#else
                      cmf = cmor_write(vars[i].cmor_varID, vars[i].data, vars[i].datatype, chunk_files[i], 0, 0, 0, nullptr);
#endif
                    }
                }
              if (cmf != 0) cdo_abort("ERROR (infile: '%s')! Function cmor_write failed!", cdo_get_stream_name(0));
            }
    }

  if (Options::cdoVerbose) cdo_print("10.4. Successfully looped over time steps.");
  if (Options::cdoVerbose) cdo_print("10. Successfully written variables via cmor_write.");
  if (Options::cdoVerbose) cdo_print("11. Start to close files, free allocated memory and, if necessary, write chunk files.");
  char **chunkdf = nullptr;
  if (mergeIDs[0]!=CMOR_UNDEFID)
    i = 1;
  if (strcmp(kv_get_a_val(kvl, "d", "y"), "y") == 0)
    chunkdf = get_chunk_des_files(kvl, vars, miptab_freqptr, i, vlistID, charname, project_id);

  char file_name[CMOR_MAX_STRING];
  for (i = 0; vars[i].cmor_varID != CMOR_UNDEFID; ++i)
    {
      if (!vars[i].help_var)
        {
          cmf = cmor_close_variable(vars[i].cmor_varID, file_name, nullptr);
          if (strcmp(file_name, "") == 0) cdo_abort("ERROR (infile: '%s')! Function cmor_write failed!", cdo_get_stream_name(0));
#if (CMOR_VERSION_MAJOR == 2)
          if (  strcmp(project_id, "CORDEX") == 0 )
            {
              if ( strcmp(kv_get_a_val(kvl, "tp", "y"), "y") == 0 )
                {
                  if (Options::cdoVerbose) cdo_print("11.1. Start to set a prefix for tracking_id.");
                  int ncid, status;
                  status = nc_open((const char *)file_name, NC_WRITE, &ncid);
                  status = nc_redef(ncid);
                  char prelim[CMOR_MAX_STRING];
                  status = nc_get_att_text(ncid,NC_GLOBAL, "tracking_id", prelim);
                  char *prefixCordex = strdup("hdl:21.14103/");
                  int lengthCombi = strlen(prelim) + strlen(prefixCordex);
                  char *track = (char *) Malloc( lengthCombi *sizeof(char));
                  sprintf(track, "%s%s", prefixCordex, prelim);
                  status = nc_put_att_text (ncid, NC_GLOBAL, "tracking_id", (size_t)lengthCombi, (const char *) track);
                  status = nc_enddef(ncid);
                  status = nc_close(ncid);
                  if ( status != NC_NOERR )
                    cdo_abort("ERROR (infile: '%s')! Could not set a prefix for tracking_id", cdo_get_stream_name(0));
                  if (Options::cdoVerbose) cdo_print("11.1. Successfully set a prefix for tracking_id.");
                  Free(track);
                  Free(prefixCordex);
                }
            }
#endif
/*
              if ( strcmp(kv_get_a_val(kvl, "tracking_prefix", "n"), "y") == 0 )
                {
                  if (Options::cdoVerbose) cdo_print("11.1. Start to set a prefix for tracking_id.");
                  CdiStreamID streamIDF = streamOpenRead(file_name);
                  int vlistIDF = streamInqVlist(streamIDF);
                  char *prelim = get_txtatt(vlistIDF, CDI_GLOBAL, "tracking_id");
                  char *prefixCordex = strdup("21.14103/");
                  size_t lengthCombi = (size_t) (strlen(prelim) + strlen(prefixCordex));
                  char *track = (char *) Malloc( lengthCombi *sizeof(char));
                  sprintf(track, "%s%s", prefixCordex, prelim);
                  cdiDefAttTxt(vlistIDF, CDI_GLOBAL, "tracking_id", lengthCombi, (const char *)track);
                  streamClose(streamIDF);
                  if (Options::cdoVerbose) cdo_print("11.1. Successfully set a prefix for tracking_id.");
                }
*/
          bool isCordexName = false;
          char cordex_file_name[CMOR_MAX_STRING];
          if (strcmp(project_id, "CORDEX") == 0 && kv_get_a_val(kvl, "cordexDir", nullptr)
              && kv_get_a_val(kvl, "cordexFileTem", nullptr))
            {
              char varname[CMOR_MAX_STRING], timename[CMOR_MAX_STRING];
              char *dummy = file_name;
              int count = 0, firsts = 0, lasts = 0;
              while (file_name[count])
                {
                  if (file_name[count] == '_')
                    {
                      if (firsts == 0) firsts = count;

                      lasts = count;
                    }
                  count++;
                }
              strncpy(varname, file_name, firsts);
              varname[firsts] = '\0';
              dummy += lasts-1;
/* Check for CMOR-bug */
              while ( *dummy != '_' )
                {
                  if ( *dummy == '-' )
                    break;
                  dummy--;
                }
              if ( *dummy != '_' )
                {
                  while ( *dummy != '_' )
                    dummy--;
                  strcpy(timename, dummy);
                  lasts = 1;
                  dummy++;
                  while ( *dummy != '_' )
                    {
                      lasts++;
                      dummy++;
                    }
                  timename[lasts] = '.';
                  timename[lasts+1] = 'n';
                  timename[lasts+2] = 'c';
                  timename[lasts+3] = '\0';
                }
/* end check for CMOR-bug */
              else
                {
                  dummy = file_name;
                  dummy += lasts;
                  strcpy(timename, dummy);
                }
              if ( ifreq == 7 )
                {
                  char smon1[12], smon2[12];
                  strncpy(smon1, &timename[5], 2);
                  strncpy(smon2, &timename[12], 2);
                  smon1[2] = '\0';
                  smon2[2] = '\0';
                  if ( atol(smon1) != 1 )
                    sprintf(smon1, "%02ld", atol(smon1)-1);
                  else
                    {
                      char syr[5];
                      strncpy(syr, &timename[1], 4);
                      syr[4] = '\0';
                      sprintf(syr, "%04ld", atol(syr)-1);
                      strncpy(&timename[1], syr, 4);

                      sprintf(smon1, "12");
                    }
                  if ( atol(smon2) != 1 )
                    sprintf(smon2, "%02ld", atol(smon2)+1);
                  else
                    {
                      char syr[12];
                      strncpy(syr, &timename[8], 4);
                      syr[4] = '\0';
                      sprintf(syr, "%04ld", atol(syr)-1);
                      strncpy(&timename[8], syr, 4);

                      sprintf(smon2, "12");
                    }
                  strncpy(&timename[5], smon1, 2);
                  strncpy(&timename[12], smon2, 2);
                }

              int cmdlen = 11+strlen(kv_get_a_val(kvl, "cordexDir", nullptr))+strlen(varname);
              char command1[cmdlen];
              std::snprintf(command1, cmdlen, "mkdir -p %s/%s", kv_get_a_val(kvl, "cordexDir", nullptr), varname);

              int dir_err = system(command1);
              if (dir_err != 0)
                {
                  cdo_warning("Could not create CORDEX compliant path for output files of cdo cmor. Files are created "
                             "in current working directory.");
                }

              sprintf(cordex_file_name, "%s/%s/%s_%s%s", kv_get_a_val(kvl, "cordexDir", nullptr), varname, varname,
                      kv_get_a_val(kvl, "cordexFileTem", nullptr), timename);

              cmdlen=5+strlen(file_name)+strlen(cordex_file_name);
              char command2[cmdlen];

              sprintf(command2, "mv %s %s", file_name, cordex_file_name);
              dir_err = system(command2);
              if (dir_err != 0)
                {
                  cdo_warning("Could not move cdo cmor output file to CORDEX compliant path.");
                  cdo_print("     File stored in:  '%s' with cmor!", file_name);
                  if (Options::silentMode)
                    cdo_warning("     File stored in:  '%s' with cmor!", file_name);

                }
              else
                {
                  isCordexName = true;
                  cdo_print("     File stored in:  '%s' with cmor!", cordex_file_name);
                  if (Options::silentMode)
                    cdo_warning("     File stored in:  '%s' with cmor!", cordex_file_name);
                }
            }
          else
            {
              char *realization = kv_get_a_val(kvl, "realization", nullptr);
              if ( !realization )
                realization = kv_get_a_val(kvl, "realization_index", nullptr);

              if ( realization[0] == '0' && realization[1] )
                {
                  char newname[CDI_MAX_NAME], oldmember[CDI_MAX_NAME],
                       newmember[CDI_MAX_NAME], chunkpath[CDI_MAX_NAME], oldchunkpath[CDI_MAX_NAME];
                  sprintf(oldmember, "r%ldi", atol(realization));
                  sprintf(newmember, "r%si", realization);

                  char *startcmp = file_name;
                  int startpattern = 0, lastSlash = 0;
                  strcpy(chunkpath,file_name);
                  int patternlength = strlen(oldmember);
                  bool oldchunkcopied = false;
/* member is once in the path, once in the file name */
                  while ( file_name[startpattern] )
                    {
                      if ( strncmp(startcmp, oldmember, patternlength) == 0 )
                            {
                              if ( !oldchunkcopied )
                                chunkpath[startpattern] = '\0';
                              else
                                chunkpath[startpattern+strlen(newmember)-patternlength] = '\0';
                              startcmp += patternlength;
                              startpattern += patternlength;
                              sprintf(newname, "%s%s%s", chunkpath, newmember, startcmp);
                              if ( !oldchunkcopied )
                                {
                                  sprintf(oldchunkpath, "%s%s", chunkpath, oldmember);
                                  oldchunkcopied = true;
                                }
                              strcpy(chunkpath, newname);
                            }
                      startcmp++; startpattern++;
                    }
                  startpattern = 0;
                  while ( newname[startpattern] )
                    {
                      if ( newname[startpattern] == '/' )
                        lastSlash = startpattern;
                      startpattern++;
                    }

                  strcpy(chunkpath, newname);
                  chunkpath[lastSlash] = '\0';
                  char command[CDI_MAX_NAME];
                  sprintf(command, "mkdir -p %s; mv %s %s;",
                          chunkpath, file_name, newname);

                  int dir_err = system(command);
                  if (dir_err != 0)
                    {
                      cdo_warning("Could not move cdo cmor output file to a path with realization=0*.");
                      cdo_print("     File stored in:  '%s' with cmor!", file_name);
                      if (Options::silentMode)
                        cdo_warning("     File stored in:  '%s' with cmor!", file_name);
                    }
                  else
                    {
                      cdo_print("     File stored in:  '%s' with cmor!", newname);
                      if (Options::silentMode)
                        cdo_warning("     File stored in:  '%s' with cmor!", newname);
                  }

                  sprintf(command, "rmdir %s*;", oldchunkpath);
                  if (Options::cdoVerbose) cdo_print("Start to remove wrong data path (r1 instead of r0*).");
                  dir_err = system(command);
                  if (dir_err != 0)
                    if (Options::cdoVerbose) cdo_print("Failed to remove '%s'", oldchunkpath);
                }
              else
                {
                  cdo_print("     File stored in:  '%s' with cmor!", file_name);
                  if (Options::silentMode)
                    cdo_warning("     File stored in:  '%s' with cmor!", file_name);
                }
            }
          if (chunkdf)
            {
              if (Options::cdoVerbose) cdo_print("11.2. Start to write a chunk description file.");
              FILE *fp = std::fopen(chunkdf[i], "w+");
              if (fp)
                {
                  if ( isCordexName )
                    fprintf(fp, "%s\n", cordex_file_name);
                  else
                    fprintf(fp, "%s\n", file_name);
                }
              else
                {
                  cdo_print("Could not open a chunk description file '%s'.", chunkdf[i]);
                  continue;
                }
              std::fclose(fp);
              if (Options::cdoVerbose) cdo_print("11.2. Successfully written a chunk description file '%s'.", chunkdf[i]);
            }
        }
    }
  if (cmf != 0) cdo_abort("ERROR (infile: '%s')! Function cmor_close_variable failed!", cdo_get_stream_name(0));
  if ( mergeIDs[0] != CMOR_UNDEFID ) Free(mergeIdx);
  if (frequency) Free(frequency);
  if (chunk_files) free_array(chunk_files);
  if (chunkdf) free_array(chunkdf);
  if (charname) Free(charname);
  cdo_stream_close(newstreamID);
  if (Options::cdoVerbose) cdo_print("11. Successfully closed files and freed allocated memory.");
}
/*
static int find_cmorvar(list_t *tester, char *cn, char * miptabfreq)
{
  KeyValues *kvcharcn = nullptr, *kvmip = nullptr;
  kvcharcn = tester->search("cmor_name");

  if ( kvcharcn )
    {
      if ( strcmp(kvcharcn->values[0], cn) == 0 )
        {
          kvmip = tester->search("project_mip_table");
          if ( kvmip )
            {
              if ( strcmp(kvmip->values[0], miptabfreq) == 0 )
                return 1;
              else
                return 0;
            }
          return 2;
        }
    }
  return 0;
}
*/
static void
read_maptab(KVList *kvl, CdoStreamID streamID, char *miptabfreq, struct mapping vars[])
{
  /***/
  /* Build mapping table from a combination of two attributes if mt does not begin with / and a directory path is given
   */
  /***/
  if (Options::cdoVerbose) cdo_print("5. Start to find, read and apply mapping table.");
  char *maptab = kv_get_a_val(kvl, "mt", nullptr);
  char *maptabdir = kv_get_a_val(kvl, "mapping_table_dir", nullptr);
  char *maptabbuild = nullptr;
  const KeyValues *kvn = kvl->search("n");
  const KeyValues *kvc = kvl->search("c");
  const KeyValues *kvcn = kvl->search("cn");
  int filetype = cdo_inq_filetype(streamID);

  if (maptab && maptabdir)
    if (maptab[0] != '/')
      {
        maptabbuild = (char *) Malloc((strlen(maptab) + strlen(maptabdir) + 2) * sizeof(char));
        sprintf(maptabbuild, "%s/%s", maptabdir, maptab);
      }
  if (maptab)
    {
      if (maptabbuild) maptab = maptabbuild;
      int vlistID = cdo_stream_inq_vlist(streamID);

      /***/
      /* Parse the table as a fortran namelist wich contains lists (=lines) of keyvalues */
      /***/
      if (Options::cdoVerbose) cdo_print("5.1 Try to read mapping table: '%s'", maptab);
      kv_insert_vals(kvl, "workfile4err", maptab, true, false);
      PMList pml = cdo_parse_cmor_file(maptab, true);
      if (!pml.size())
        {
          cdo_warning("5.1. In parsing the mapping table '%s':\n          Mapping table could not be parsed. Operator "
                     "continues.",
                     maptab);
          return;
        }
      /***/
      /* If a variable selector name or code is given in cmdline, the corresponding variable is picked from Infile and
       * mapped. */
      /* Only the first value of name/code given in the cmdline is processed */
      /* If no variable selector is given, process all variables and map via name and code */
      /***/
      /* However, if the mapping table contains a keyvalue pair for name or code with more than one value, */
      /* the corresponding variable has a character coordinate and requires special treatment */
      /* This is tested once before mapping. If the special variable equals the variable which is to map, */
      /* the special treatment begins with fct addcharvar */
      /***/
      /* Different CMOR variables are built with one model variable. */
      /* Consequently, for one model variable more than one mapping table entry can exist */
      /* As a second identification argument, the mapping table name (miptabfreq) is used */
      /***/
      /* If no variable selector is given in the mapping table, it is assumed that the infile variable is already named
       * like cmor_name */
      /***/

      if (kvn)
        {
          if (Options::cdoVerbose) cdo_print("5. No character axis is built from several variables"
                                   "(only possible if values for 'name' are provided in the mapping table).");
          if (kvn->nvalues > 1)
            cdo_warning("5.1. In applying the mapping table '%s':\n          Only the first value of commandline "
                        "variable selection key 'name' is processed.",
                        maptab);
          maptab_via_cmd(maptab, pml, kvn->values[0].c_str(), vlistID, "name", (char *) kvcn->values[0].c_str(),
                         miptabfreq, filetype, maptab);
        }
      else if (kvc)
        {
          if (Options::cdoVerbose) cdo_print("5. No character axis is built from several variables"
                                   "(only possible if values for 'code' are provided in the mapping table).");
          if (kvc->nvalues > 1)
            cdo_warning("5.1. In applying the mapping table '%s':\n          Only the first value of commandline "
                       "variable selection key 'code' is processed.",
                       maptab);
          maptab_via_cmd(maptab, pml, kvc->values[0].c_str(), vlistID, "code", (char *) kvcn->values[0].c_str(),
                         miptabfreq, filetype, maptab);
        }
      else if (kvcn)
        {
          maptab_via_cn(maptab, pml, kvcn->values, vlistID, kvcn->nvalues, miptabfreq, filetype, vars, true);
        }
      else
        {
          if (Options::cdoVerbose) cdo_print("5. No character axis is built from several variables"
                                   "(only possible if cmor_name is provided in command line).");
          for (int varID = 0; varID < vlistNvars(vlistID); ++varID)
            {
              /***/
              /* Begin with Code in case infile is of type GRB */
              /***/
              if (filetype == CDI_FILETYPE_GRB || filetype == CDI_FILETYPE_GRB2)
                if (maptab_via_key(maptab, pml, vlistID, varID, "code", miptabfreq))
                  {
                    if (Options::cdoVerbose) cdo_print("5.1. Successfully mapped varID '%d' via code.", varID);
                    continue;
                  }
              if (maptab_via_key(maptab, pml, vlistID, varID, "name", miptabfreq))
                {
                  if (Options::cdoVerbose) cdo_print("5.1. Successfully mapped varID '%d' via name.", varID);
                  continue;
                }
              if (maptab_via_key(maptab, pml, vlistID, varID, "code", miptabfreq))
                {
                  if (Options::cdoVerbose) cdo_print("5.1. Successfully mapped varID '%d' via code.", varID);
                  continue;
                }
              /***/
              /* In case corresponding mapping table entry does not contain a variable selector attribute */
              /***/
              if (maptab_via_key(maptab, pml, vlistID, varID, "cmor_name", miptabfreq))
                {
                  if (Options::cdoVerbose) cdo_print("5.1. Successfully mapped varID '%d' via cmor_name.", varID);
                  continue;
                }
              cdo_warning("5.1. In applying the mapping table '%s':\n          Could not map variable with id '%d'.",
                         maptab, varID);
            }
        }
      /***/
      /* In case a requested variable needs an auxilliary variable, the latter may be mapped later. */
      /* If a mapping table exists is saved here */
      /***/
      kv_insert_vals(kvl, "mtproof", maptab, true, false);
      cdo_print("Mapping Table = '%s'.", maptab);
      if (maptabbuild) Free(maptabbuild);
    }
  else if (Options::cdoVerbose)
    cdo_print("5. No mapping table found.");
}

static void
replace_key(KVList *kvl, KeyValues kv, const char *newkey)
{
  char **values = (char **) Malloc((kv.nvalues+1) * sizeof(char *));
  int k = 0;
  for ( k = 0; k<kv.nvalues; k++ )
    values[k] = strdup(kv.values[k].c_str());
  values[kv.nvalues] = nullptr;
  kvl->remove(kv.key);
  kvl->append(newkey, values, k);
  free_array(values);
}

static void
parse_cmdline(KVList *kvl, std::vector<std::string> &params, int nparams, const char *ventry)
{
/* Already set params++ in main function */
  if (kvl->parse_arguments(nparams-1, params) != 0 )
    cdo_abort("ERROR (infile: '%s')! Could not parse command line.",cdo_get_stream_name(0));

  std::vector<KeyValues> keystorm, keystosubs;
  for (const auto &kv : *kvl)
    {
              const char *short_key = check_short_key((char *)kv.key.c_str());
              if (short_key)
                {
                  if (kv.key != short_key) keystosubs.push_back(kv);
                }
              else
                {
                  cdo_warning("Unknown commandline keyword: '%s'\n", kv.key);
                  keystorm.push_back(kv);
                }
    }
  for ( size_t i = 0; i < keystosubs.size(); ++i)
    {
      replace_key(kvl, keystosubs[i], check_short_key((char *)keystosubs[i].key.c_str()));
    }
  for ( size_t i = 0; i < keystorm.size(); ++i)
    {
      kvl->remove(keystorm[i].key);
    }
}

static char *
get_mip_table(char *params, KVList *kvl, char *project_id, bool print)
{
  char *miptab;
  if (print && Options::cdoVerbose) cdo_print("2.2. Start to find a MIP table file.");
  if (!params) cdo_abort("ERROR (infile: '%s')! First parameter not passed. A MIP table file is required.", cdo_get_stream_name(0));
  if (file_exist(params, false, "MIP table", print))
    {
      miptab = strdup(params);
      int i = 0, j = 0;
      while (params[i])
        {
          if (params[i] == '/') j = i;
          i++;
        }
      char miptabdir[1024];
      char cwd[1024];
      getcwd(cwd, sizeof(cwd));
      cwd[strlen(cwd)] = '\0';
      if (params[0] == '/')
        {
          strncpy(miptabdir, params, j + 1);
          miptabdir[j+1] = '\0';
        }
      else if (j == 0 )
        strcpy(miptabdir, cwd);
      else
        {
          strcpy(miptabdir, cwd);
          strcat(miptabdir, "/");
          strncat(miptabdir, params, j + 1);
          miptabdir[strlen(cwd)+j + 1] = '\0';
        }
      kv_insert_vals(kvl, "mip_table_dir", miptabdir, true, false);
      if (print) cdo_print("MIP table file = '%s'.", miptab);
      return miptab;
    }
  else
    {
      if (print && Options::cdoVerbose)
        cdo_print("Try to build a path with additional configuration attributes:\n          'mip_table_dir' and "
                 "'project_id'\n          in order to use '%s' as MIP-table.",
                 params);
      char *miptabdir = kv_get_a_val(kvl, "mip_table_dir", nullptr);
      if (miptabdir && project_id)
        {
#if (CMOR_VERSION_MAJOR == 2)
          {
            miptab = (char *) Malloc((strlen(miptabdir) + strlen(project_id) + strlen(params) + 3) * sizeof(char));
            sprintf(miptab, "%s/%s_%s", miptabdir, project_id, params);
          }
#elif (CMOR_VERSION_MAJOR == 3)
          {
            miptab = (char *) Malloc((strlen(miptabdir) + strlen(project_id) + strlen(params) + 8) * sizeof(char));
            sprintf(miptab, "%s/%s_%s.json", miptabdir, project_id, params);
          }
#endif
          file_exist(miptab, true, "MIP table", print);
          if (print) cdo_print("MIP table file = '%s'", miptab);
          return miptab;
        }
      else
        cdo_abort("ERROR (infile: '%s')! In finding the MIP table:\n          Could not find attribute 'mip_table_dir'.", cdo_get_stream_name(0));
    }

  return nullptr;
}

static char *
freq_from_path(char *mip_table)
{
  char *freq = mip_table;
  int fpos = 0, k = 0, j = 0;
  while (*(mip_table + j))
    {
      j++;
      if (*(mip_table + j) == '/') k = j + 1;
      if (*(mip_table + j) == '_' && *(mip_table + j + 1)) fpos = j + 1;
    }
  freq += k;
  if (fpos > k) freq += fpos - k;
  return freq;
}

static int
get_miptab_freq(char *mip_table, char *project_id)
{
  int miptab_freq = 0;
  char *freq = freq_from_path(mip_table);
  if (freq != nullptr)
    {
      if (strstr(freq, "yr") || strstr(freq, "Yr"))
        miptab_freq = 11;
      else if (strstr(freq, "mon") || strstr(freq, "Mon"))
        miptab_freq = 12;
      else if (strstr(freq, "day") || strstr(freq, "Day"))
        miptab_freq = 13;
      else if (strstr(freq, "6h"))
        miptab_freq = 14;
      else if (strstr(freq, "3h"))
        miptab_freq = 15;
      else if (strstr(freq, "1h") || strstr(freq, "AERhr"))
        miptab_freq = 16;
      else if (strstr(freq, "sem"))
        miptab_freq = 17;
      else if (strstr(freq, "dec"))
        miptab_freq = 18;
      else if (strstr(freq, "subhr"))
        miptab_freq = 19;

      if (strcmp(freq, "Oclim") == 0)
        miptab_freq = 1;
      else if (strcmp(freq, "Oyr") == 0)
        miptab_freq = 2;
      else if (strcmp(freq, "cfMon") == 0)
        miptab_freq = 3;
      else if (strcmp(freq, "day") == 0)
        miptab_freq = 4;
      else if (strcmp(freq, "6hrPlev") == 0 && strcmp(project_id, "CMIP5") == 0)
        miptab_freq = 5;
      else if (strcmp(freq, "6hrPlevPt") == 0)
        miptab_freq = 5;
      else if (strcmp(freq, "6hrLev") == 0)
        miptab_freq = 6;
      else if (strcmp(freq, "E1hrClimMon") == 0)
        miptab_freq = 7;
      else if (strcmp(freq, "E3hrPt") == 0)
        miptab_freq = 8;
    }
  return miptab_freq;
}

static void
check_cmdline_mapping(KVList *kvl)
{
  char *name = kv_get_a_val(kvl, "n", nullptr);
  char *code = kv_get_a_val(kvl, "c", nullptr);
  char *cn = kv_get_a_val(kvl, "cn", nullptr);
  if ((name && code))
    cdo_abort("ERROR (infile: '%s')! Mapping via command line failed. Only one variable selector of 'name' and 'code' is allowed.", cdo_get_stream_name(0));
  if ((name && !cn) || (code && !cn))
    cdo_abort("ERROR (infile: '%s')! Mapping via command line failed. A corresponding 'cmor_name' is needed.", cdo_get_stream_name(0));
}

static char *
get_project_id(KVList *kvl, char *params)
{
  if (Options::cdoVerbose) cdo_print("2.1. Start to check whether 'project_id' or 'mip_era' is denoted.");
  char *project_id = nullptr, *dummy, *dummy2;
  dummy = kv_get_a_val(kvl, "project_id", nullptr);
  dummy2 = kv_get_a_val(kvl, "mip_era", nullptr);

  char tester[CDI_MAX_NAME];
  strcpy(tester, params);
  char *testerpointer = tester;
  int underscore = 0, slash = 0, testint = 0;
  while ( tester[testint] )
    {
      if ( tester[testint] == '_' )
        underscore = testint;
      if ( tester[testint] == '/' )
        slash = testint;
      testint++;
    }
  if ( underscore > slash )
    {
       if ( slash )
         testerpointer += slash+1;
       testerpointer[underscore-slash-1] = '\0';
    }
#if defined(CMOR_VERSION_MAJOR)
#if (CMOR_VERSION_MAJOR == 2)
  {
    if (!dummy && !dummy2)
      {
        if ( strcmp(testerpointer, tester) != 0 )
          {
            if (Options::cdoVerbose) 
              cdo_print("Could not find attribute 'project_id'.\n          "
                   "Try to use substring from MIP-table input '%s' as project_id.", testerpointer);
            project_id = strdup(testerpointer);
          }
        else
          cdo_abort("ERROR (infile: '%s')! Attribute 'project_id' is required.", cdo_get_stream_name(0));
      }
    else if (!dummy)
      cdo_abort("ERROR (infile: '%s')! Cannot produce CMIP6 standard with CMOR2.\n          "
               "Value for attribute 'project_id' is required.", cdo_get_stream_name(0));
    else
      project_id = strdup(dummy);
  }
#elif (CMOR_VERSION_MAJOR == 3)
  {
    if (!dummy && !dummy2)
      {
        if ( strcmp(testerpointer, tester) != 0 )
          {
            if (Options::cdoVerbose) 
              cdo_print("Could not find attribute 'project_id'.\n          "
                   "Try to use substring from MIP-table input '%s' as project_id.", testerpointer);
            project_id = strdup(testerpointer);
          }
        else
          cdo_abort("ERROR (infile: '%s')! Attribute 'mip_era' or 'project_id' is required.", cdo_get_stream_name(0));
      }
    else if (!dummy2)
      {
        if (Options::cdoVerbose) cdo_print("You have not provided 'mip_era' but only 'project_id'."
          " If you try to produce CMIP5 standard,\n          It is recommended to use CMOR2 for this job instead.");
        project_id = strdup(dummy);
      }
    else
      project_id = strdup(dummy2);
  }
#endif
#else
  cdo_abort("ERROR (infile: '%s')! Cannot check CMOR version: Missing makro CMOR_VERSION_MAJOR", cdo_get_stream_name(0));
#endif

  if (Options::cdoVerbose) cdo_print("2.1. Successfully found project_id / mip_era: '%s'.", project_id);
  return project_id;
}

static int
cmor_load_and_set_table(KVList *kvl, char *param0, char *project_id, char **mip_table)
{
  int table_id = 0, cmf = 0;
#if (CMOR_VERSION_MAJOR == 3)
  Free(*mip_table);
  *mip_table = get_mip_table(param0, kvl, project_id, false);
#endif
  cmf = cmor_load_table(*mip_table, &table_id);
  if (cmf != 0) cdo_abort("ERROR (infile: '%s')! Function cmor_load_table failed!", cdo_get_stream_name(0));
  cmf = cmor_set_table(table_id);
  if (cmf != 0) cdo_abort("ERROR (infile: '%s')! Function cmor_set_table failed!", cdo_get_stream_name(0));
  return table_id;
}

#endif

void *
CMOR(void *process)
{
  cdo_initialize(process);

#ifdef HAVE_LIBCMOR
  int nparams = cdo_operator_argc();

  if (nparams < 1) cdo_abort("ERROR (infile: '%s')! No parameter found. Need at least a MIP-table.", cdo_get_stream_name(0));

  auto params = cdo_get_oper_argv();
  char *miptableInput = strdup(params[0].c_str());
  //copy everything from params except the first
  params = std::vector<std::string>(params.begin() + 1, params.end());

  /* Define cmdline list and read cmdline */
  const char *pmlistHelper[] = { "cmdline" };
  KVList kvl;
  parse_cmdline(&kvl, params, nparams, pmlistHelper[0]);

  /* Check whether a command line mapping is active */
  check_cmdline_mapping(&kvl);

  /* Config files are read with descending priority. */
  read_config_files(&kvl);

  /* Get project_id, mip_table and mip_table frequency*/
  if (Options::cdoVerbose) cdo_print("2. Start to find a MIP table and to deduce a frequency from MIP table file.");
  char *project_id = get_project_id(&kvl, miptableInput);
  char *mip_table = get_mip_table(miptableInput, &kvl, project_id, true);
#if (CMOR_VERSION_MAJOR == 3)
  int lenmt = strlen(mip_table) - 4 ;
  std::vector<char> miptemp(lenmt);
  strncpy(miptemp.data(), mip_table, lenmt);
  miptemp[strlen(mip_table) - 5] = '\0';
  Free(mip_table);
  mip_table = strdup(miptemp.data());
#endif
  int miptab_freq = get_miptab_freq(mip_table, project_id);

  char *miptab_freqptr = strdup(freq_from_path(mip_table));
  kv_insert_vals(&kvl, "miptab_freq", miptab_freqptr, true, false);

  if (Options::cdoVerbose)
    cdo_print("2. Successfully found a MIP table '%s' and deduced a MIP table frequency '%s'.", mip_table,
             miptab_freqptr);

  if (Options::cdoVerbose) cdo_print("3. Start to open infile '%s'.", cdo_get_stream_name(0));
  const auto streamID = cdo_open_read(0);
  const auto vlistID = cdo_stream_inq_vlist(streamID);

  if (Options::cdoVerbose) cdo_print("3. Successfully opened infile '%s'.", cdo_get_stream_name(0));

  if (Options::cdoVerbose) cdo_print("4. Start to check attributes.");
  /* Short keys from rtu, mt, gi must be included similar to global atts */
  add_globalhybrids(&kvl, vlistID);
  /* Allow time units from infile */
  check_required_time_units(&kvl, vlistInqTaxis(vlistID));

  /* Check for attributes and member name */
  check_attr(&kvl, project_id, vlistID);
  check_mem(&kvl, project_id);
  if (Options::cdoVerbose) cdo_print("4. Successfully checked global attributes.");

  /* dump_global_attributes(pml, streamID); */

  struct mapping *vars = construct_var_mapping(vlistID);

  /* read mapping table */
  read_maptab(&kvl, streamID, miptab_freqptr, vars);

  int time_axis = 0, calendar = 0;

  setup_dataset(&kvl, streamID, &calendar, project_id);

  int table_id = cmor_load_and_set_table(&kvl, miptableInput, project_id, &mip_table);

  int mergeIDs[150];
  mergeIDs[0] = CMOR_UNDEFID;
  register_all_dimensions(&kvl, streamID, vars, table_id, project_id, miptab_freq, &time_axis, mergeIDs);
  write_variables(&kvl, streamID, vars, miptab_freq, time_axis, calendar, miptab_freqptr, project_id, mergeIDs);

  destruct_var_mapping(vars);
  Free(mip_table);
  Free(project_id);
  Free(miptab_freqptr);
  /* Free(miptableInput); */

#else
  cdo_abort("CMOR support not compiled in!");
#endif
  cdo_finish();
  return 0;
}
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
