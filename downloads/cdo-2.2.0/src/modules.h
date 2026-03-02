/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/
#ifndef MODULES_H
#define MODULES_H

#include <iostream>
#include <map>
#include <array>
#include <vector>

#include "module_static_maps.h"

// Obase uses input name as base name for files e.g 'test' gets used as test_001 test_002 which are created inside the operator
#define OBASE -1
#define INTERNAL 0
#define EXPOSED 1

/***
  type definition for module functions loaded from a custom module
  */
using dyn_oper_t = void (*)(void *arg);

enum PositionRestrictions
{
  NoRestriction = 0,
  FilesOnly = 1,
  OnlyFirst = 2
};

struct Alias
{
  Alias(const std::string &_alias, const std::string &_original);
  std::string alias;
  std::string original;
};

struct module_constraints
{
  short streamInCnt;   // Number of input streams
  short streamOutCnt;  // Number of output streams
  PositionRestrictions pos_restriction = NoRestriction;
};

struct module_t
{
  std::string mod_name;
  void *(*func)(void *);               // Module
  const char **help;                   // Help
  short mode;          // Module mode: 0:intern 1:extern
  short number;        // Allowed number type
  std::vector<std::string> operators;  // Operator names
  module_constraints constraints;
  std::vector<Alias> aliases;

  module_t(const std::string &mod_name, void *(*p_func)(void *), const char **p_help, const std::vector<std::string> &p_opers,
           short p_m, short p_n, short p_siC, short p_soC, PositionRestrictions p_pos_restriction,
           const std::vector<Alias> &p_aliases = {});
  module_t(){};
  std::string toString();
  int
  get_stream_in_cnt() const
  {
    return constraints.streamInCnt;
  }
  int
  get_stream_out_cnt() const
  {
    return constraints.streamOutCnt;
  }
  int
  get_number() const
  {
    return number;
  }
  int
  get_mode() const
  {
    return mode;
  }
  int
  get_pos_restriction()
  {
    return constraints.pos_restriction;
  }
};

/***
  vector for library handles for loaded custom modules
  */
extern std::vector<void *> custom_modules_lib_handles;

std::string extract_operator_name(const std::string &operatorCommand);
void extract_name_and_argument(const std::string &command, std::string &operatorName, std::string &operatorArgument);

std::string find_similar_operators(const std::string &operatorName);

/***
  Key: operator alias / Value: operator original name
 */

void register_operators(const std::string &p_mod_name, const module_t &p_mod);

// void *(*operatorModule(const char *operatorName))(void *);

std::map<std::string, module_t>::iterator find_module(const std::string &operatorName);
std::vector<std::string> get_module_operator_names(std::string module_name);
module_t &get_module(const std::string &operatorName);

std::string get_original(const std::string &operatorName);

void init_modules();

const char **operator_help(const std::string &operatorName);
// int operatorStreamInCnt(const char *operatorName);
// int operatorStreamOutCnt(const char *operatorName);
int operator_stream_number(std::string &operatorName);

std::string get_module_name_to(const std::string &operatorName);
std::vector<std::string> get_sorted_operator_name_list();

#endif /* MODULES_H */
