#ifndef MODULE_INFO_H
#define MODULE_INFO_H

#include <string>
#include <map>

const std::string s_obase = "obase";
const std::string s_arbIn = "arbitrary";
const std::string s_filesOnly = "filesOnly";
const std::string s_onlyFirst = "onlyFirst";
const std::string s_noOutput = "noOutput";

struct ModListOptions
{
  bool printAll = false;
  bool operInfoRequested = false;
  std::map<const std::string, int> opt
      = { { s_obase, false }, { s_arbIn, false }, { s_filesOnly, false }, { s_onlyFirst, false }, { s_noOutput, false } };

  bool requested(const std::string &name);
  bool mod_info_requested();
  bool parse_request(const std::string &requestString);
};

void operator_print_list(ModListOptions &p_modListOpt);
#endif
