#ifndef CDO_QUERY_H
#define CDO_QUERY_H

#include <string>

#include <cdi.h>

#include "pmlist.h"

std::string set_query_parameter(const KVList &kvlist, CdiQuery *query);
std::string set_query_parameter(const std::string &params, CdiQuery *query);

#endif
