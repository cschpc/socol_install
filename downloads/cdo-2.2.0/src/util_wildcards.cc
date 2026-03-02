/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <cstdlib>
#include <cstddef>
#include <cstring>
#include <vector>
#include <string>
#ifdef HAVE_WORDEXP_H
#include <wordexp.h>
#endif
#ifdef HAVE_GLOB_H
#include <glob.h>
#endif
#ifdef HAVE_FNMATCH_H
#include <fnmatch.h>
#endif

#include "util_wildcards.h"

#ifdef HAVE_GLOB_H
static int
get_glob_flags(void)
{
  int glob_flags = 0;

#ifdef GLOB_NOCHECK
  glob_flags |= GLOB_NOCHECK;
#endif
#ifdef GLOB_TILDE
  glob_flags |= GLOB_TILDE;
#endif

  return glob_flags;
}
#endif

static int
find_wildcard(const char *string, size_t len)
{
  int status = 0;

  if (len > 0)
    {
      if (string[0] == '~') status = 1;

      if (status == 0)
        {
          for (size_t i = 0; i < len; ++i)
            if (string[i] == '?' || string[i] == '*' || string[i] == '[')
              {
                status = 1;
                break;
              }
        }
    }

  return status;
}

// used in griddes.cc
char *
expand_filename(const char *string)
{
  char *filename = nullptr;

  if (find_wildcard(string, strlen(string)))
    {
#ifdef HAVE_GLOB_H
      const auto glob_flags = get_glob_flags();
      glob_t glob_results;
      glob(string, glob_flags, 0, &glob_results);
      if (glob_results.gl_pathc == 1) filename = strdup(glob_results.gl_pathv[0]);
      globfree(&glob_results);
#endif
    }

  return filename;
}

#ifdef HAVE_WORDEXP_H
// Expands all input file wildcards and removes the wildcard while inserting all expanded files into argv
std::vector<std::string>
expand_wild_cards(std::vector<std::string> argv)
{
  int flags = WRDE_UNDEF;
  wordexp_t glob_results;

  bool applyActive = false;
  int bracketsOpen = 0;
  for (size_t idx = 1; idx < argv.size(); idx++)
    {
      // if argv[idx] contains wildcard (* or [?]+), multiple ** are ignored
      if (argv[idx].compare(0, 6, "-apply") == 0)
        {
          applyActive = true;
          continue;
        }
      if (argv[idx].size() == 1 && argv[idx][0] == '[')
        {
          bracketsOpen++;
          continue;
        }
      if (argv[idx].size() == 1 && argv[idx][0] == ']')
        {
          bracketsOpen--;
          if (bracketsOpen == 0) { applyActive = false; }
          continue;
        }
      if (argv[idx][0] != '-' && argv[idx].find_first_of("*?[ ") != std::string::npos)
        {
          const auto status = wordexp(argv[idx].c_str(), &glob_results, flags);
          if (status != 0)
            {
              fprintf(stderr, "%s: ", __func__);
              if (status == WRDE_BADCHAR)
                fprintf(stderr,
                        "Argument '%s' contains one of the following unsupported unquoted characters: <newline>, `|', "
                        "`&', `;', `<', `>', `(', `)', `{', `}'.\n",
                        argv[idx].c_str());
              else if (status == WRDE_NOSPACE)
                fprintf(stderr, "Not enough memory to store the result.\n");
              else if (status == WRDE_SYNTAX)
                fprintf(stderr, "Shell syntax error in '%s'\n", argv[idx].c_str());
              else if (status == WRDE_BADVAL)
                fprintf(stderr, "Undefined shell variable in '%s'\n", argv[idx].c_str());
              else
                fprintf(stderr, "wordexp() returns an error.\n");
              exit(EXIT_FAILURE);
            }
          // range based insert (glob_results.we_wordv is inserted before wildcard
          if (std::string(glob_results.we_wordv[0]).find_first_of("*?[ ") == std::string::npos)
            {
              auto insertAt = idx + 1;
              if (applyActive == false)
                {
                  argv.insert(argv.begin() + insertAt, "]");
                  argv.insert(argv.begin() + insertAt, "[");
                  insertAt += 1;
                }
              argv.insert(argv.begin() + insertAt, glob_results.we_wordv, glob_results.we_wordv + glob_results.we_wordc);
              argv.erase(argv.begin() + idx);
            }
          // delete wildcard
          wordfree(&glob_results);
        }
    }

  return argv;
}
#else
std::vector<std::string>
expand_wild_cards(std::vector<std::string> argv)
{
  return argv;
}
#endif

#ifdef HAVE_FNMATCH_H
int
wildcardmatch(const char *pattern, const char *string)
{
  return fnmatch(pattern, string, 0);
}
#else
// The wildcardmatch function checks if two given strings match.
// The first string may contain wildcard characters
// * --> Matches with 0 or more instances of any character or set of characters.
// ? --> Matches with any one character.
// source code from http://www.geeksforgeeks.org/wildcard-character-matching/
int
wildcardmatch(const char *w, const char *s)
{
  // If we reach at the end of both strings, we are done
  if (*w == '\0' && *s == '\0') return 0;

  // Make sure that the characters after '*' are present in second string.
  // This function assumes that the first string will not contain two consecutive '*'
  if (*w == '*' && *(w + 1) != '\0' && *s == '\0') return 1;

  // If the first string contains '?', or current characters of both strings match
  if ((*w == '?' && *s != '\0') || *w == *s) return wildcardmatch(w + 1, s + 1);

  // If there is *, then there are two possibilities
  // a) We consider current character of second string
  // b) We ignore current character of second string.
  if (*w == '*') return wildcardmatch(w + 1, s) || wildcardmatch(w, s + 1);

  return 1;
}
#endif

int
wildcardmatch(const std::string &pattern, const std::string &string)
{
  return wildcardmatch(pattern.c_str(), string.c_str());
}
