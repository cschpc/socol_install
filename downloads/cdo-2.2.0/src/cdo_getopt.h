#ifndef _CDO_GETOPT_H
#define _CDO_GETOPT_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <map>
#include <set>
#include <vector>
#include <string>
#include <functional>
#include <cstdio>
#include <memory>

#include "cdo_output.h"

struct cdo_option_2
{

  std::function<void(std::string argument)> effect = [](std::string) {
    printf("ERROR");
    exit(0);
  };
  bool hasArgument = false;
  bool abortOnExecute = false;
  std::string argument = std::string();
  std::vector<std::string> description = std::vector<std::string>();
  std::string argumentDescription = std::string();
  std::string name = std::string();
  std::string shortName = std::string();
  bool isInternal = false;
  bool argumentIsOptional = false;
  std::string linkedEnvironmentVariable = std::string();
  std::string default_value = "";

  cdo_option_2 *
  add_effect(const std::function<void(std::string argument)> p_effect)
  {
    hasArgument = true;
    effect = p_effect;
    return this;
  }
  cdo_option_2 *
  set_internal(bool p_is_internal)
  {
    isInternal = p_is_internal;
    return this;
  }

  cdo_option_2 *
  add_default(const std::string &p_default_value)
  {
    default_value = p_default_value;
    return this;
  }

  cdo_option_2 *
  add_effect(const std::function<void(void)> p_effect)
  {
    effect = [p_effect](std::string) { p_effect(); };
    return this;
  }

  cdo_option_2 *
  describe_argument(const std::string &desc)
  {
    hasArgument = true;
    argumentDescription = desc;
    return this;
  }
  cdo_option_2 *
  set_argument_optional(bool p_isOptional)
  {
    argumentIsOptional = p_isOptional;
    return this;
  }

  cdo_option_2 *
  add_help(std::string desc...)
  {

    description.push_back(desc);
    return this;
  }

  template <typename... T>
  cdo_option_2 *
  add_help(const std::string &arg, T... args)
  {
    add_help(arg);
    add_help(args...);
    return this;
  }

  std::string
  get_arg()
  {
    return argument;
  }
  void
  execute()
  {
    if (effect == nullptr) { cdo_abort("effect not set"); }
    effect(get_arg());
  }

  cdo_option_2 *
  aborts_program(bool aborts)
  {
    abortOnExecute = aborts;
    return this;
  }

  cdo_option_2 *
  add_env_var(const std::string &p_envVarName)
  {
    linkedEnvironmentVariable = p_envVarName;
    return this;
  }

  void
  evaluate_env_var()
  {
    const char *env_var = getenv(linkedEnvironmentVariable.c_str());
    if (env_var != nullptr)
      {
        if (hasArgument)
          {
            if (*env_var == '\0') { cdo_abort("Error: %s defined but has no value", linkedEnvironmentVariable); }
            effect(std::string(env_var));
          }
        else { effect(std::string()); }
      }
  }
};

class CLIOptions
{
private:
  static std::vector<std::shared_ptr<cdo_option_2>> optionRegistry;
  static std::map<std::string, std::shared_ptr<cdo_option_2>> optionMap;

  static std::vector<std::shared_ptr<cdo_option_2>> envvarRegistry;
  static std::map<std::string, std::shared_ptr<cdo_option_2>> envvarMap;
  static bool stderr_is_tty;

public:
  static std::string pad_size_terminal(char padWith, const std::string &sourround = std::string());
  const static int EXIT_REQUESTED;
  const static int ABORT_REQUESTED;
  const static int padding;
  static bool print_settings;
  static bool print_envvars;

  static void set_tty_status(bool tty_status);
  static int parse(std::vector<std::string> p_argv);
  static void get_env_vars();
  static std::shared_ptr<cdo_option_2> &option(const std::string &p_name, const std::string &p_short_form = std::string());
  static std::shared_ptr<cdo_option_2> &envvar(const std::string &p_name);
  static void usage();
  static void print_available_options();
  static void print_options_help();
  static void print_envvar_help();
  static void print_registry(const std::vector<std::shared_ptr<cdo_option_2>> &p_registry);
  static void option_from_envvar(const std::string &p_envvarName);
};
#endif /* _CDO_GETOPT_H */
