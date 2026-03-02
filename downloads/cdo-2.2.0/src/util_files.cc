#include <cdi.h>

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <sys/stat.h>

#include "util_files.h"
#include "cdo_options.h"
#include "cdo_vlist.h"
#include "readline.h"

#include "cdo_default_values.h"
bool
FileUtils::file_exists(const char *filename)
{
  struct stat buf;
  const auto status = stat(filename, &buf);
  return (status == 0) && (S_ISREG(buf.st_mode) && buf.st_size > 0);
}

bool
FileUtils::file_exists(const std::string &filename)
{
  return FileUtils::file_exists(filename.c_str());
}

bool
FileUtils::user_file_overwrite(const char *filename)
{
  auto status = false;

  if (!Options::silentMode && cdo::stdinIsTerminal && cdo::stderrIsTerminal)
    {
      fprintf(stderr, "File %s already exists, overwrite? (yes/no): ", filename);
      char line[1024];
      cdo::readline(stdin, line, 1024);
      char *pline = line;
      while (isspace((int) *pline)) pline++;
      const auto len = strlen(pline);
      if (len == 3)
        {
          if (strncmp(pline, "yes", 3) == 0 || strncmp(pline, "YES", 3) == 0) status = true;
        }
      else if (len == 1)
        {
          if (pline[0] == 'y' || pline[0] == 'Y') status = true;
        }
    }

  return status;
}

off_t
FileUtils::size(const char *filename)
{
  off_t filesize = 0;

  if (filename[0] != '(') /* && filename[1] != 'p') */
    {
      struct stat buf;
      if (stat(filename, &buf) == 0) filesize = buf.st_size;
    }

  return filesize;
}

static void
gen_file_suffix(int filetype, const char *refname, size_t maxlen, char *filesuffix, int vlistID)
{
  auto lready = false;
  auto lcompsz = false;

  if (filetype == CdoDefault::FileType && CdoDefault::DataType == -1 && CdoDefault::Byteorder == -1)
    {
      size_t len = 0;
      if (refname != nullptr && *refname != 0 && *refname != '-' && *refname != '.') len = strlen(refname);

      if (len > 2)
        {
          const char *result = strrchr(refname, '.');
          if (result != nullptr && result[1] != 0)
            {
              const int firstchar = tolower(result[1]);
              switch (firstchar)
                {
                case 'g':
                  if (CdoDefault::FileType == CDI_FILETYPE_GRB || CdoDefault::FileType == CDI_FILETYPE_GRB2) lready = true;
                  break;
                case 'n':
                  if (CdoDefault::FileType == CDI_FILETYPE_NC || CdoDefault::FileType == CDI_FILETYPE_NC2
                      || CdoDefault::FileType == CDI_FILETYPE_NC4 || CdoDefault::FileType == CDI_FILETYPE_NC4C
                      || CdoDefault::FileType == CDI_FILETYPE_NC5)
                    lready = true;
                  break;
                case 's':
                  if (CdoDefault::FileType == CDI_FILETYPE_SRV) lready = true;
                  break;
                case 'e':
                  if (CdoDefault::FileType == CDI_FILETYPE_EXT) lready = true;
                  break;
                case 'i':
                  if (CdoDefault::FileType == CDI_FILETYPE_IEG) lready = true;
                  break;
                }
            }

          // if ( lready )  strncat(filesuffix, result, maxlen-1);
          if (lready && ((len = strlen(result)) < (maxlen - 1)))
            {
              while (len--)
                {
                  if (*result == '.' || isalnum(*result)) strncat(filesuffix, result, 1);
                  result++;
                }
            }
        }
    }

  if (!lready)
    {
      strncat(filesuffix, streamFilesuffix(CdoDefault::FileType), maxlen - 1);
      if (CdoDefault::FileType == CDI_FILETYPE_GRB && vlist_is_szipped(vlistID)) lcompsz = true;
    }

  if (CdoDefault::FileType == CDI_FILETYPE_GRB && Options::cdoCompType == CDI_COMPRESS_SZIP) lcompsz = true;
  if (lcompsz) strncat(filesuffix, ".sz", maxlen - 1);
}

void
FileUtils::gen_suffix(char *filesuffix, size_t maxlen, int filetype, int vlistID, const char *refname)
{
  if (strncmp(cdo::File_Suffix, "NULL", 4) != 0)
    {
      if (cdo::File_Suffix[0] != 0) { strncat(filesuffix, cdo::File_Suffix, maxlen - 1); }
      else { gen_file_suffix(filetype, refname, maxlen, filesuffix, vlistID); }
    }
}
