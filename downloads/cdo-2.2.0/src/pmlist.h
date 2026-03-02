/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/
#ifndef PMLIST_H
#define PMLIST_H

#include <list>
#include <vector>
#include <string>

#include "listbuffer.h"
#include "namelist.h"

struct KeyValues
{
  int nvalues = 0;
  std::string key;
  std::vector<std::string> values;
};

// clang-format off
class  // KVList
#ifdef WARN_UNUSED
[[gnu::warn_unused]]
#endif
KVList : public std::list<KeyValues>
// clang-format on
{
public:
  std::string name;

  void print(FILE *fp = stderr) const;
  int parse_arguments(int argc, const std::vector<std::string> &argv);
  const KeyValues *search(const std::string &key) const;
  void remove(const std::string &inkey);
  void append(const char *, const char *const *, int);
  char *get_first_value(const char *key, const char *replacer);
};

// clang-format off
class  // PMList
#ifdef WARN_UNUSED
[[gnu::warn_unused]]
#endif
PMList : public std::list<KVList>
// clang-format on
{
public:
  const KVList *searchKVListVentry(const std::string &key, const std::string &value, const std::vector<std::string> &entry);
  const KVList *getKVListVentry(const std::vector<std::string> &entry);
  void print(FILE *fp = stderr);
  void read_namelist(FILE *fp, const char *name);
  void read_cmor_table(FILE *fp, const char *name);
};

int parse_namelist(PMList &pmlist, NamelistParser &parser, char *buf, bool cdocmor);
int parse_list_buffer(NamelistParser &p, ListBuffer &listBuffer);

#endif
