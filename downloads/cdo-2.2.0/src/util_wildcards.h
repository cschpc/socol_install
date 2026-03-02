/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/
#ifndef UTIL_WILDCARDS_H
#define UTIL_WILDCARDS_H

#include <vector>
#include <string>

char *expand_filename(const char *args);

std::vector<std::string> expand_wild_cards(std::vector<std::string> argv);

int wildcardmatch(const char *w, const char *s);
int wildcardmatch(const std::string &w, const std::string &s);

#endif
