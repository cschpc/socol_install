/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <string.h>
#include <assert.h>
#include "cdo_getopt.h"
#include "util_string.h"
#include <cstdio>
#include <cstdlib>

#include <vector>
#include <set>
#include <map>

#define CLIOP_DBG false

#include <sys/ioctl.h>
#include <stdio.h>
#include <unistd.h>

std::string
CLIOptions::pad_size_terminal(char padWith, const std::string &sourround)
{
  size_t terminal_width = 120;
  if (stderr_is_tty)
    {
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
      CONSOLE_SCREEN_BUFFER_INFO csbi;
      int columns, rows;

      GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &csbi);
      columns = csbi.srWindow.Right - csbi.srWindow.Left + 1;
#elif __APPLE__
#include <TargetConditionals.h>
#if TARGET_IPHONE_SIMULATOR
      // iOS, tvOS, or watchOS Simulator
#elif TARGET_OS_MACCATALYST
      // Mac's Catalyst (ports iOS API into Mac, like UIKit).
#elif TARGET_OS_IPHONE
      // iOS, tvOS, or watchOS device
#elif TARGET_OS_MAC
      // Other kinds of Apple platforms
#else
#error "Unknown Apple platform"
#endif
#elif __linux__ || __unix__ || defined(_POSIX_VERSION)
      struct winsize w;
      ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
      terminal_width = w.ws_col;
#endif
    }

  std::string line;
  if (sourround.empty()) { line = std::string(terminal_width, padWith); }
  else
    {
      size_t front_pad = 3;
      size_t num_spaces = 2;
      size_t space_left = terminal_width - front_pad - num_spaces - sourround.size();

      line = std::string(front_pad, padWith);
      line += " " + sourround + " ";
      line += std::string(space_left, padWith);
    }
  return line + "\n";
}

void
CLIOptions::usage()
{
  std::cerr << pad_size_terminal('-');
  fprintf(stderr, "  Usage : cdo  [Options]  Operator1  [-Operator2  [-OperatorN]]\n");
  std::cerr << pad_size_terminal('-');

  fprintf(stderr, "\n");
  std::cerr << pad_size_terminal('=', "Options");

  set_text_color(stderr, BLUE);
  print_options_help();
  reset_text_color(stderr);
  std::cerr << pad_size_terminal('-');

  std::cerr << pad_size_terminal('=', "Environment Variables");
  set_text_color(stderr, BLUE);
  print_envvar_help();
  reset_text_color(stderr);
  std::cerr << pad_size_terminal('-');

  fprintf(stderr, "\n");

  fprintf(stderr, "  Operators:\n");
  fprintf(stderr, "    Use option --operators for a list of all operators.\n");

  fprintf(stderr, "\n");
  fprintf(stderr, "  CDO version %s, Copyright (C) 2002-2023 MPI für Meteorologie\n", VERSION);
  fprintf(stderr, "  This is free software and comes with ABSOLUTELY NO WARRANTY\n");
  fprintf(stderr, "  Report bugs to <https://mpimet.mpg.de/cdo>\n");
}

std::vector<std::shared_ptr<cdo_option_2>> CLIOptions::optionRegistry;
std::map<std::string, std::shared_ptr<cdo_option_2>> CLIOptions::optionMap;

std::vector<std::shared_ptr<cdo_option_2>> CLIOptions::envvarRegistry;
std::map<std::string, std::shared_ptr<cdo_option_2>> CLIOptions::envvarMap;

const int CLIOptions::EXIT_REQUESTED = -2;
const int CLIOptions::ABORT_REQUESTED = -1;
bool CLIOptions::print_settings = false;
bool CLIOptions::print_envvars = false;
const int CLIOptions::padding = 46;

bool CLIOptions::stderr_is_tty = true;

void
CLIOptions::set_tty_status(bool tty_status)
{
  stderr_is_tty = tty_status;
}

int
CLIOptions::parse(std::vector<std::string> p_argv)
{
  int retval = p_argv.size();
  for (size_t i = 1; i < p_argv.size(); ++i)
    {
      Debug(CLIOP_DBG, "Checking: %s", p_argv[i]);
      const std::map<std::string, std::shared_ptr<cdo_option_2>>::iterator &it = optionMap.find(p_argv[i]);
      if (it == optionMap.end())
        {
          retval = i;
          break;
        }
      if (it->second->hasArgument)
        {
          i++;
          if (i >= p_argv.size() || ((p_argv[i][0] == '-') && !std::get<0>(tokenize_comma_seperated_int_list(p_argv[i]))))
            {
              if (it->second->default_value.size() > 0) { it->second->argument = it->second->default_value; }
              else if (!it->second->argumentIsOptional)
                cdo_abort("missing argument for %s", p_argv[i - 1]);
            }
          else { it->second->argument = p_argv[i]; }
        }
      Debug(CLIOP_DBG, "executing option %s", (*it).first);
      it->second->execute();
      if (it->second->abortOnExecute)
        {
          retval = EXIT_REQUESTED;
          break;
        }
      if (i >= p_argv.size() - 1)
        {
          // TODO: missing err msg
          retval = ABORT_REQUESTED;
          break;
        }
    }
  if (print_settings) print_registry(optionRegistry);
  if (print_envvars) print_registry(envvarRegistry);

  return retval;
}

void
CLIOptions::get_env_vars()
{
  for (auto &setting_ptr : envvarRegistry)
    {
      const char *envVarValue = getenv(setting_ptr->name.c_str());
      if (envVarValue)
        {
          if (!setting_ptr->hasArgument || (setting_ptr->hasArgument && *envVarValue))
            {
              Debug(CLIOP_DBG, "Executing envvar %s", setting_ptr->name);
              setting_ptr->argument = envVarValue ? std::string(envVarValue) : std::string();
              setting_ptr->execute();
            }
        }
    }
}

void
CLIOptions::print_registry(const std::vector<std::shared_ptr<cdo_option_2>> &p_registry)
{
  for (auto &it : p_registry)
    {
      if (it->argument.size() > 0) fprintf(stderr, "%s = %s\n", it->name.c_str(), it->argument.c_str());
    }
}

std::shared_ptr<cdo_option_2> &
CLIOptions::envvar(const std::string &p_name)
{
  if (envvarMap.find(p_name) == envvarMap.end())
    {
      envvarRegistry.push_back(std::make_shared<cdo_option_2>());
      auto &newEnvVar = envvarRegistry.back();
      envvarMap[p_name] = newEnvVar;
      newEnvVar->name = p_name;
    }
  else { cdo_abort("Environment Variable already registered!"); }
  return envvarRegistry.back();
}

std::shared_ptr<cdo_option_2> &
CLIOptions::option(const std::string &p_name, const std::string &p_shortForm)
{
  std::string name = "--" + p_name;
  std::string shortForm = p_shortForm;
  Debug(CLIOP_DBG, "registering key: %s", name);
  if (optionMap.find(name) != optionMap.end()) { cdo_abort("option name already exists: %s", name); }

  optionRegistry.push_back(std::make_shared<cdo_option_2>());
  optionMap[name] = optionRegistry.back();
  optionMap[name]->name = name;

  if (!shortForm.empty())
    {
      if (optionMap.find(shortForm) == optionMap.end())
        {
          shortForm = "-" + p_shortForm;
          optionMap[shortForm] = optionRegistry.back();
          optionMap[name]->shortName = shortForm;
        }
      else { cdo_abort("Option short form already exists: %s", shortForm); }
    }
  return optionRegistry.back();
}
void
CLIOptions::print_envvar_help()
{
  for (auto &it : envvarMap)
    {
      if (!it.second->isInternal)
        {
          std::string helpString = std::string(4, ' ');
          helpString += it.second->name;
          if (it.second->hasArgument) helpString += " <" + it.second->argumentDescription + "> ";
          int spaceLeft = padding - helpString.size();
          for (auto &line : it.second->description) { helpString += std::string(spaceLeft, ' ') + line + "\n"; }
          std::cerr << helpString;
        }
    }
}
void
CLIOptions::print_options_help()
{
  for (auto &iter : optionMap)
    {
      if (iter.first.size() != 2 && !iter.second->isInternal)
        {
          auto option = iter.second;
          std::string helpString = "    ";
          if (!option->shortName.empty()) { helpString += option->shortName + ", "; }
          else { helpString += "    "; }
          helpString += option->name + " ";
          if (option->hasArgument) helpString += " <" + option->argumentDescription + "> ";
          if (option->hasArgument && option->argumentDescription.empty())
            {
              std::cerr << "error: help argument of " << option->name << " has no description!" << std::endl;
              exit(0);
            }
          int spaceLeft = padding - helpString.size();
          if (spaceLeft <= 0)
            {
              helpString += "\n";
              spaceLeft = padding;
            }
          for (auto &line : option->description)
            {
              helpString += std::string(spaceLeft, ' ') + line + "\n";
              spaceLeft = padding;
            }
          if (option->description.empty()) helpString += "\n";
          std::cerr << helpString;
        }
    }
}
#include <algorithm>

void
CLIOptions::option_from_envvar(const std::string &p_envvarName)
{
  if (envvarMap.find(p_envvarName) == envvarMap.end()) { cdo_abort("Error envvar %s does not exist", p_envvarName); }
  std::string ENVVAR_SUFFIX = "CDO_";
  std::string optionName = p_envvarName.substr(ENVVAR_SUFFIX.size(), p_envvarName.size());

  std::transform(optionName.begin(), optionName.end(), optionName.begin(), [](unsigned char c) { return std::tolower(c); });

  optionName = "--" + optionName;

  if (optionMap.find(optionName) != optionMap.end())
    {
      cdo_abort("Error autogenerated name %s for envvar option %s does already exist!", optionName, p_envvarName);
    }
  optionRegistry.push_back(std::make_shared<cdo_option_2>());
  auto &newOption = optionRegistry.back();
  auto &envVar = envvarMap[p_envvarName];
  newOption->name = optionName;
  newOption->hasArgument = envVar->hasArgument;
  newOption->description = envVar->description;
  newOption->argumentDescription = envVar->argumentDescription;
  newOption->description = { "This option is generated from " + p_envvarName + " and will overwrite it.",
                             "See help of corresponding environment variable." };
  newOption->effect = envVar->effect;
  newOption->isInternal = envVar->isInternal;
  newOption->default_value = envVar->default_value;
  optionMap[optionName] = newOption;
}

void
CLIOptions::print_available_options()
{
  for (auto &iter : optionMap) { std::cout << iter.first << std::endl; }
  std::cout << "_---------------------------------_" << std::endl;
  for (auto &iter : optionMap)
    {
      if (iter.first.size() == 2) std::cout << iter.second->shortName[1] << (iter.second->hasArgument ? ":" : "");
    }
}
