/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <stdlib.h>

#include <cdi.h>
#include "table.h"
#include "cdo_output.h"
#include "util_files.h"

std::string getenv_string(const std::string &envVar);

namespace cdo
{

int
define_table(const std::string &tablearg)
{
  auto tablename = tablearg.c_str();

  auto tableID = FileUtils::file_exists(tablename) ? tableRead(tablename) : CDI_UNDEFID;

  if (tableID == CDI_UNDEFID)
    {
      const auto &tablepath = getenv_string("CD_TABLEPATH");
      if (tablepath.size())
        {
          auto tablefile = tablepath + "/" + tablename;
          if (FileUtils::file_exists(tablefile)) tableID = tableRead(tablefile.c_str());
        }
    }

  if (tableID == CDI_UNDEFID) tableID = tableInq(-1, 0, tablename);

  if (tableID == CDI_UNDEFID) cdo_abort("table <%s> not found", tablename);

  return tableID;
}

}
