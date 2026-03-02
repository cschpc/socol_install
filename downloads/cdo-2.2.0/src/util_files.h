#ifndef UTIL_FILES_H
#define UTIL_FILES_H

#include <sys/types.h>
#include <string>
namespace FileUtils
{
bool file_exists(const std::string &filename);
bool file_exists(const char *filename);
bool user_file_overwrite(const char *filename);
off_t size(const char *filename);
void gen_suffix(char *filesuffix, size_t maxlen, int filetype, int vlistID, const char *refname);
}  // namespace FileUtils

#endif
