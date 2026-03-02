#include "module_static_maps.h"
#include "modules.h"

std::map<std::string, module_t> &
get_modules()
{
  static std::map<std::string, module_t> modules;
  return modules;
};

std::map<std::string, std::string> &
get_aliases()
{
  static std::map<std::string, std::string> aliases;
  return aliases;
}

std::map<std::string, std::string> &
get_module_map()
{
  static std::map<std::string, std::string> modules_map;
  return modules_map;
}
