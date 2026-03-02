/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida
          Oliver Heidmann

*/

#include "modules.h"
#include "util_string.h"
#include "cdo_output.h"

#include <cstring>
#include <algorithm>  // for std::sort()

int add_alias(const std::string &alias, const std::string &original);
void add_module(const std::string &module_name, const module_t &new_module);

/* removes '-' from operator string and returns copy of parameter operatorCommand */
/***
  vector for library handles for loaded custom modules
  */
std::vector<void *> custom_modules_lib_handles;

/***
  Contains added modules as values and their names as key
  */

Alias::Alias(const std::string &_alias, const std::string &_original) : alias(_alias), original(_original) {}

// clang-format off
module_t::module_t(const std::string &p_mod_name,
                   void *(*p_func)(void *),
                   const char **p_help,
                   const std::vector<std::string> &p_opers,
                   short p_m, short p_n, short p_siC, short p_soC,
                   PositionRestrictions p_pos_restriction,
                   const std::vector<Alias> &p_aliases)

                 : mod_name(p_mod_name),
                   func(p_func),
                   help(p_help),
                   mode(p_m),
                   number(p_n),
                   operators(p_opers),
                   constraints({p_siC, p_soC, p_pos_restriction}),
                   aliases(p_aliases)
{
  add_module(p_mod_name, *this);
  for (const auto &al : aliases) add_alias(al.alias, al.original);
}
// clang-format on

std::string
module_t::toString()
{
  std::string inp = (get_stream_in_cnt() >= 0) ? std::to_string(get_stream_in_cnt()) : "Arbitrary";
  std::string out = (get_stream_out_cnt() >= 0) ? std::to_string(get_stream_out_cnt()) : "Output base";
  std::string restriction = "none";

  if (get_pos_restriction() == OnlyFirst) restriction = "Can only be the first operator";
  if (get_pos_restriction() == FilesOnly)
    (restriction != Green("none")) ? restriction = Yellow("Can only use files as input.")
                                   : restriction += Yellow(", Can only use files as input");

  std::string desc = "Input: " + inp + ", Ouput: " + out + ", Restricton: " + restriction;

  return desc;
}

void
extract_name_and_argument(const std::string &command, std::string &operatorName, std::string &operatorArgument)
{
  const std::string delimiter = ",";

  const size_t start = (command[0] == '-') ? 1 : 0;

  size_t len = command.find(delimiter);
  if (len == std::string::npos)
    {
      len = command.size();
      operatorArgument = "";
    }
  else { operatorArgument = command.substr(len + 1, std::string::npos); }

  operatorName = command.substr(start, len - start);
}

std::string
extract_operator_name(const std::string &operatorCommand)
{
  const std::string delimiter = ",";

  const size_t start = (operatorCommand[0] == '-') ? 1 : 0;

  size_t len = operatorCommand.find(delimiter);
  if (len == std::string::npos) len = operatorCommand.size();

  std::string oper_name = operatorCommand.substr(start, len - start);
  return oper_name;
}

bool
is_alias(const std::string &operatorName)
{
  return (get_aliases().find(operatorName) != get_aliases().end());
}

std::string
get_original(const std::string &operatorName)
{
  std::string name = extract_operator_name(operatorName);
  return (is_alias(name) ? get_aliases()[name] : name);
}
/**
 * @param a pointer to a string/substring
 * @param b pointer to a string/substring
 * @param alen length of string a
 * @param blen length of string b
 * @retval true if a is similar to b
 * @retval false if a is not similar to b
 *
 * Recursive function for finding substrings of a operator name that match other operators.
 */
static bool
similar(const char *a, const char *b, unsigned long alen, unsigned long blen)
{
  if (alen > 2 && blen > 2 && strstr(b, a)) return true;

  while (*a && *b && *a == *b)
    {
      a++;
      b++;
    }
  if (!*a && !*b) return true;

  //  printf("%d %d %s %s\n", alen, blen, a, b);

  if (alen >= 2 && blen >= 1 && *a && similar(a + 1, b, alen - 2, blen - 1)) return true;

  if (alen >= 1 && blen >= 2 && *b && similar(a, b + 1, alen - 1, blen - 2)) return true;

  return false;
}

/**
 * @param original string tested for similarity to \p other
 * @param other string that \p original will be compared to
 * @retval true if original and other are similar
 * @retval false if not
 *
 * Wrapper function for #similar() to parse c++ strings to c strings
 */
static bool
similar(const std::string &original, const std::string &other)
{
  return (similar(original.c_str(), other.c_str(), original.size(), other.size()));
}

/**
 * @param operatorName operator name
 * @retval true if #modules_map contains \p operatorName
 * @retval false if not
 */
bool
operator_name_exists(const std::string &operatorName)
{
  if (get_module_map().find(operatorName) != get_module_map().end()) return true;

  if (get_aliases().find(operatorName) != get_aliases().end()) return true;

  return false;
}

/**
 * @param moduleName module name
 * @retval true if #modules contains \a moduleName
 * @retval false if not
 */
/*
static bool
module_map_contains(std::string moduleName)
{
  if (modules.find(moduleName) != modules.end())
    {
      return true;
    }
  else
    {
      cdo_abort("Module %s not found", moduleName);
    }
  return false;
}
*/
/***
 * function for finding similar operator names for the given string
 * @param operatorName operator name to find similar operators for
 * @returns A string with all found names. The string is seqmented into lines
 * with a max length of 75 characters
 */
std::string
find_similar_operators(const std::string &operatorName)
{
  std::string found_similar_operators = "";
  size_t lines = 1;
  constexpr size_t line_length = 105;

  if (operatorName != "")
    {
      // Searching for similar operator names in operator to module map
      for (const auto &str : get_module_map())
        {
          if (similar(string_to_lower(operatorName), str.first))
            {
              if (found_similar_operators.size() + str.first.size() > lines * line_length)
                {
                  found_similar_operators += "\n";
                  lines++;
                }
              found_similar_operators += str.first;
              found_similar_operators += " ";
            }
        }
      // Searching for similar operator names in aliases to original map
      for (const auto &str : get_aliases())
        {
          if (similar(string_to_lower(operatorName), str.first))
            {
              if (found_similar_operators.size() + str.first.size() > lines * line_length)
                {
                  found_similar_operators += "\n";
                  lines++;
                }
              found_similar_operators += str.first;
              found_similar_operators += " ";
            }
        }
    }

  return found_similar_operators;
}

/**
 * @param operatorName operator name.
 * @retval true if \p operatorName exists.
 * @retval false if \p operatorName is not in #modules_map
 *
 * Checks if given \p operatorName is in #modules_map. Else returns false.

 * Checks if \p operatorName is not a file.

 * If no matching operator is found checks for similar operators using find_similar_operators().
 *
 *  \note If \p operatorName is a file name the program will exit.
 */
static bool
check_operator(const std::string &operatorName)
{
  if (operator_name_exists(operatorName)) { return true; }
  else if (operatorName == "") { cdo_abort("Operator name missing!"); }
  else
    {
      // Checking if the operatorname is an existing file name
      auto fp = std::fopen(operatorName.c_str(), "r");
      if (fp)
        {
          std::fclose(fp);
          fprintf(stderr, "Use commandline option -h for help.\n");
          cdo_abort("Operator missing, %s is a file on disk!", operatorName);
        }
      // Operator is no filename
      // Checking for similar operators
      fprintf(stderr, "Operator >%s< not found!\n", operatorName.c_str());
      fprintf(stderr, "Similar operators are:\n");
      const auto found_similar_operators = find_similar_operators(operatorName);

      if (found_similar_operators.size() > 0) { std::cerr << found_similar_operators << std::endl; }
      else { fprintf(stderr, "(not found)\n"); }

      exit(EXIT_FAILURE);
    }

  return false;
}

void
register_operators(const std::string &mod_name, const module_t &mod)
{
  for (const auto &y : mod.operators)
    {
      if (!operator_name_exists(y)) { get_module_map()[y] = mod_name; }
      else { std::cout << "Tried to add operator but the operator name already exists" << std::endl; }
    }
}

void
register_operators_from_modules()
{
  for (const auto &x : get_modules()) { register_operators(x.first, x.second); }
}

/***
 * Adds a module and its operators to cdo.
 * Adds the module to modules
 * Adds the operators of modules to modules_map
 * @param new_module newly constructed module
 * @note: if an error happens while adding the new module cdo will exit.
 */
void
add_module(const std::string &module_name, const module_t &new_module)
{
  if (get_modules().find(module_name) == get_modules().end())
    {
      get_modules()[module_name] = new_module;
      for (const auto &operatorName : new_module.operators)
        {
          // if the operator name is not already in the map or in the aliases
          if (!operator_name_exists(operatorName)) { get_module_map()[operatorName] = module_name; }
          else { cdo_abort("Tried to add operator but the operator name already exists"); }
        }
    }
  else { cdo_abort("Module %s name already exists", module_name); }
}

/**
 * adds an key value pair to #modules_map with alias as key and originals name as value
 * @param alias new alias to be added
 * @param original original operator name
 */
int
add_alias(const std::string &alias, const std::string &original)
{
  auto iter_original = get_module_map().find(original);
  auto iter_alias = get_aliases().find(alias);

  if (iter_alias != get_aliases().end())
    {
      printf("WARNING: alias %s could not be added: it already exists\n", alias.c_str());
      return -1;
    }

  if (iter_original == get_module_map().end())
    {
      printf("ERROR: alias %s could not be added: operator %s does not exist\n", alias.c_str(), original.c_str());
      return -2;
    }

  if (get_module_map().find(alias) != get_module_map().end())
    {
      printf("ERROR alias %s could not be added: alias name already exists as an operator\n", alias.c_str());
      return -3;
    }

  get_aliases()[alias] = original;

  return 0;
}

/***
 * returns the module of given operator
 *
 * parameters:
 *      std::string operatorName -> name of the operator
 * return value:
 *      std::string -> name of the operators module
 */
std::string
get_module_name_to(const std::string &operatorName)
{
  // if not true the programm will exit in function check_operator
  if (check_operator(operatorName))
    {
      if (get_module_map().count(operatorName) > 0) { return get_module_map()[operatorName]; }
      else if (get_aliases().count(operatorName) > 0) { return get_module_map()[get_aliases()[operatorName]]; }
    }
  // Only to quell the warning. Function should never reach this.
  return "";
}

/**
 * @fn void *(*operatorModule(const char *operatorName))(void *)

 * returns the function of the module of \a operatorName
 * @param operatorName name of the operator
 * @returns the function of the #module_t of \a operatorName
 */
/* not used
void *(*operatorModule(const std::string &operatorName))(void *)
{
  std::string operator_name = std::string(operatorName);
  return get_modules()[get_module_name_to(operatorName)].func;
}
*/

/***
 * returns help for given operator name
 * if there is no help nullptr will be returned
 * @param operatorName operator name
 * @return vector of strings containing the help
 */
const char **
operator_help(const std::string &operatorName)
{
  return get_modules()[get_module_name_to(get_original(operatorName))].help;
}

/***
 * Returns the number of input streams a operator takes.
 * returns -1 for a unlimited number of input streams.
 * @param operatorName operator name
 * @retval short
 */
/* not used
int
operatorStreamInCnt(const char *operatorName)
{
  std::string operator_name = std::string(operatorName);
  return get_modules()[get_module_name_to(operator_name)].streamInCnt;
}
*/
/***
 * Returns the number of output streams a given operator has.
 * returns -1 if obase is used
 * @param operatorName operator name
 * @return 1 for CDI_REAL, 2 for CDI_COMP (complex), 3 for CDI_REAL and CDI_COMP
 * @retval short
 */
/* not used
int
operatorStreamOutCnt(const char *operatorName)
{
  std::string operator_name = std::string(operatorName);
  return get_modules()[get_module_name_to(operator_name)].streamOutCnt;
}
*/
/***
 * Returns the number type this operator works with
 * @param operatorName operator name
 * @reval short
 */
int
operator_stream_number(std::string &p_operatorName)
{
  return get_modules()[get_module_name_to(p_operatorName)].get_number();
}

/***
 * Creates a sorted vector with all operator names and alisases excluding all modules that are marked as internal
 * @return a sorted std::vector containing all operator names and aliases
 * excluding all operators which modules are marked as internal
 */
std::vector<std::string>
get_sorted_operator_name_list()
{
  std::vector<std::string> names;
  for (const auto &operator_module_names_pair : get_module_map())
    {
      if (get_modules()[operator_module_names_pair.second].get_mode() == 1) { names.push_back(operator_module_names_pair.first); }
    }
  // adding operators names from alias_map
  for (const auto &alias : get_aliases()) { names.push_back(alias.first); }

  std::sort(names.begin(), names.end());

  return names;
}

std::vector<std::string>
get_no_output_operator_list()
{
  std::vector<std::string> names;
  for (const auto &operator_module_names_pair : get_module_map())
    {
      if (get_modules()[operator_module_names_pair.second].get_mode() == 1
          && get_modules()[operator_module_names_pair.second].get_stream_out_cnt() == 0)
        {
          names.push_back(operator_module_names_pair.first);
        }
    }
  // adding operators names from alias_map
  for (const auto &alias : get_aliases())
    {
      const auto &original = alias.second;
      if (get_modules()[get_module_map()[original]].get_mode() == 1
          && get_modules()[get_module_map()[original]].get_stream_out_cnt() == 0)
        {
          names.push_back(alias.first);
        }
    }

  std::sort(names.begin(), names.end());

  return names;
}

void
operatorPrintAll(void)
{
  int number_of_chars = 0;
  std::string tab = "   ";
  int tab_width = tab.size();
  // using a set because it sorts the operators alphabetically on its own
  std::vector<std::string> sorted_operator_names = get_sorted_operator_name_list();

  std::cout << tab;
  for (const auto &operatorName : sorted_operator_names)
    {
      if (number_of_chars > 85)
        {
          number_of_chars = tab_width;
          std::cerr << std::endl << tab;
        }

      std::cerr << " " << operatorName;
      number_of_chars += 1 + operatorName.size();
    }

  std::cerr << std::endl;
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

/***
 * Prints all operator names and their short descriptions
 * Aliases are listed and point to their original operator name.
 * If the module is not documented the description is empty
 * If a module has only one operator the short module description is listed
 * If the operator is not documented the description is empty
 */
void
operator_print_list(bool print_no_output)
{
  std::vector<std::string> output_list = print_no_output ? get_no_output_operator_list() : get_sorted_operator_name_list();

  auto list_length = output_list.size();

  // help variables

  for (size_t out_list_idx = 0; out_list_idx < list_length; out_list_idx++)
    {
      auto &current_op_name = output_list[out_list_idx];
      const auto *current_module = &get_modules()[get_module_name_to(current_op_name)];
      if (get_aliases().find(current_op_name) != get_aliases().end())
        {
          output_list[out_list_idx] += get_spacing_for(16, current_op_name) + "--> " + get_aliases()[current_op_name];
        }
      else if (current_module->help)
        {
          // add spaceing and saving output line to the output list
          auto description = get_operator_description(current_op_name, current_module->help);
          output_list[out_list_idx] += get_spacing_for(16, current_op_name) + description;
        }
      std::string in_out_info = " (" + std::to_string(current_module->get_stream_in_cnt()) + "|"
                                + std::to_string(current_module->get_stream_out_cnt()) + ")";
      output_list[out_list_idx] += get_spacing_for(90, output_list[out_list_idx]) + in_out_info;
    }
  // print generated output list
  for (const std::string &str : output_list) { std::cout << str << std::endl; }
}

std::map<std::string, module_t>::iterator
find_module(const std::string &operator_name)
{
  std::string operName = operator_name;

  auto it = get_aliases().find(operName);
  if (it != get_aliases().end()) { operName = it->second; }

  auto modMapIter = get_module_map().find(operName);
  auto modIter = get_modules().end();
  if (!(modMapIter == get_module_map().end())) { modIter = get_modules().find(modMapIter->second); }
  return modIter;
}

std::vector<std::string>
get_module_operator_names(std::string module_name)
{
  auto &modules = get_modules();
  if (std::islower(static_cast<unsigned char>(module_name[0]))) {module_name[0] = std::toupper(module_name[0]); }
  if (modules.find(module_name) == modules.end()) { return std::vector<std::string>(); }
  return modules[module_name].operators;
}

module_t &
get_module(const std::string &operator_name)
{
  return get_modules()[get_module_name_to(get_original(operator_name))];
}
