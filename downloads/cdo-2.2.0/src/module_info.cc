#include "module_info.h"
#include "modules.h"
#include "mpmo_color.h"
#include "modules.h"
#include "util_string.h"
#include <algorithm>
#include <iostream>
#include <functional>
#include <string>

typedef std::function<bool(module_t &mod)> ModuleQuery;

bool
ModListOptions::requested(const std::string &name)
{
  return opt[name];
}

bool
ModListOptions::mod_info_requested()
{
  return (operInfoRequested || printAll || requested(s_obase) || requested(s_arbIn) || requested(s_filesOnly)
          || requested(s_onlyFirst) || requested(s_arbIn) || requested(s_noOutput));
}

bool
ModListOptions::parse_request(const std::string &requestString)
{
  auto all = true;
  const auto splitString = split_with_seperator(requestString, ',');

  if (requestString.size() > 0)
    {
      all = false;
      for (size_t i = 0; i < splitString.size(); ++i)
        {
          auto it = find_module(splitString[i]);
          if (it != get_modules().end())
            {
              operInfoRequested = true;
              std::cerr << splitString[i] << ": " << it->second.toString() << std::endl;
            }
          else
            {
              if (opt.find(splitString[i]) != opt.end()) { opt[splitString[i]] = 1; }
              else
                {
                  std::cerr << "option " << splitString[i] << " not found" << std::endl;
                  return false;
                }
            }
        }
    }
  printAll = all;

  return true;
}

static std::string
get_operator_description(const std::string &p_current_op_name, const char **help)
{
  std::string description = "";
  unsigned long cur_help_idx = 0;
  std::string line;
  unsigned long operator_section = 0;

  // search for operator section
  size_t help_size = 0;
  while (help[help_size]) help_size++;
  if (!help_size) return description;
  while (operator_section == 0 && cur_help_idx < help_size - 1)
    {
      line = help[++cur_help_idx];
      if (line.find("OPERATORS") != std::string::npos) { operator_section = cur_help_idx; }
    }
  // if no operator section is found
  if (operator_section == 0)
    {
      cur_help_idx = 0;
      line = help[0];
      std::string name_section = help[0];
      bool help_contains_name = false;
      // search for the operator name in the description
      while (!line.empty())
        {
          line = help[++cur_help_idx];
          if (line.find(p_current_op_name) != std::string::npos) { help_contains_name = true; }
          name_section += line;
        }
      // if the name was found save description for later use
      if (help_contains_name) { description = name_section.substr(name_section.find_first_of('-') + 2, name_section.size()); }
    }
  else
    {
      line = help[++operator_section];
      // search the operator section for current operator line
      while (line.find(p_current_op_name + " ") == std::string::npos && !line.empty() && operator_section < help_size - 1)
        {
          line = help[++operator_section];
        }
      // if operator line found save description for later use
      if (!line.empty() && line.find("    " + p_current_op_name + " ") != std::string::npos)
        {
          auto op_name_start = line.find_first_not_of(" \t");

          description = line.substr(line.find_first_not_of(" \t", op_name_start + p_current_op_name.size()), line.size());
        }
    }

  return description;
}

// helper function for setting the spacing in operator_print_list
static std::string
get_spacing_for(int p_space, const std::string &str)
{
  std::string spacing = "";
  for (int i = str.size(); i <= p_space; ++i) spacing += " ";
  return spacing;
}

static std::string
operatorGetShortInfoString(std::string &current_op_name, const module_t &p_module)
{
  std::string shortInfo = current_op_name;
  if (get_aliases().find(current_op_name) != get_aliases().end())
    {
      shortInfo += std::string(get_spacing_for(16, current_op_name) + "--> " + get_aliases()[current_op_name]);
    }
  else if (p_module.help)
    {
      // add spaceing and saving output line to the output list
      auto description = get_operator_description(current_op_name, p_module.help);
      shortInfo += get_spacing_for(16, current_op_name) + description;
    }
  std::string in_out_info
      = "(" + std::to_string(p_module.get_stream_in_cnt()) + "|" + std::to_string(p_module.get_stream_out_cnt()) + ")";
  shortInfo += get_spacing_for(90, shortInfo) + in_out_info;
  return shortInfo;
}

void
operator_print_list(std::function<bool(module_t &)> selectionCriteria)
{
  std::vector<std::string> output_list;

  // for (size_t out_list_idx = 0; out_list_idx < list_length; out_list_idx++)
  for (auto &current_op_name : get_sorted_operator_name_list())
    {
      module_t &current_module = get_modules()[get_module_name_to(current_op_name)];
      if (selectionCriteria(current_module)) { output_list.push_back(operatorGetShortInfoString(current_op_name, current_module)); }
    }
  // print generated output list
  for (const std::string &str : output_list) { std::cout << str << std::endl; }
}

void
operator_print_list(ModListOptions &p_opt)
{
  set_text_color(stderr, GREEN);

  if (p_opt.printAll == true)
    {
      operator_print_list([](module_t &) { return true; });
    }
  else
    {

      ModuleQuery defaultModuleQuery = [](module_t &) -> bool { return false; };
      ModuleQuery runquestDefaultModuleQuery = [](module_t &) -> bool { return true; };

      // clang-format off
      ModuleQuery hasObase  = p_opt.requested(s_obase)     ? [](module_t &mod) -> bool { return mod.get_stream_out_cnt() == -1;        } : defaultModuleQuery;
      ModuleQuery hasNoOut  = p_opt.requested(s_noOutput)  ? [](module_t &mod) -> bool { return mod.get_stream_out_cnt() ==  0;        } : defaultModuleQuery;
      ModuleQuery hasArb    = p_opt.requested(s_arbIn)     ? [](module_t &mod) -> bool { return mod.get_stream_in_cnt()  == -1;        } : defaultModuleQuery;
      ModuleQuery filesOnly = p_opt.requested(s_filesOnly) ? [](module_t &mod) -> bool { return mod.get_pos_restriction() == FilesOnly; } : defaultModuleQuery;
      ModuleQuery onlyFirst = p_opt.requested(s_onlyFirst) ? [](module_t &mod) -> bool { return mod.get_pos_restriction() == OnlyFirst; } : defaultModuleQuery;
      // clang-format on

      operator_print_list(
          [&](module_t &mod) { return (hasObase(mod) || hasArb(mod) || hasNoOut(mod) || filesOnly(mod) || onlyFirst(mod)); });
    }

  reset_text_color(stderr);

  return;
}
