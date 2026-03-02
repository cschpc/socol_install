/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida
          Oliver Heidmann

*/

#include <vector>
#include <iostream>
#include "modules.h"
#include "cdo_output.h"
#include "cdo_apply.h"

static ApplyStatus errState = ApplyStatus::OK;

static std::string
parse_arg(const std::string &oper)
{
  if (oper == "")
    {
      fprintf(stderr, "cdo apply: no operator given for apply.\n");
      errState = ApplyStatus::MISSING_ARG;
    }

  const auto mod = get_module(oper);
  if (mod.get_stream_in_cnt() != 1)
    {
      fprintf(stderr, "cdo apply: operator %s can not be used with apply.\n", oper.c_str());
      if (mod.get_stream_in_cnt() == -1)
        {
          fprintf(stderr, "           %s has variable input.\n", oper.c_str());
          errState = ApplyStatus::ARG_VARIABLE_INPUT;
        }

      if (mod.get_stream_in_cnt() == 0)
        {
          fprintf(stderr, "           %s has no input.\n", oper.c_str());
          errState = ApplyStatus::ARG_NO_INPUT;
        }
    }
  if (mod.get_stream_in_cnt() > 1)
    {
      fprintf(stderr, "           %s has more than one output.\n", oper.c_str());
      errState = ApplyStatus::ARG_TOO_MANY_OUT;
    }

  return oper;
}

static std::vector<std::string>
generate_tokens(const std::string &p_oper)
{
  std::vector<std::string> result;

  auto end = p_oper.find(' ');
  auto start = 0;
  while (end != std::string::npos && errState == ApplyStatus::OK)
    {
      auto oper = p_oper.substr(start, end - start);
      if (oper[0] == '-') oper = parse_arg(oper);
      result.push_back(oper);
      start = end + 1;
      end = p_oper.find(' ', start);
    }
  auto oper = p_oper.substr(start, end - start);
  if (oper[0] == '-') oper = parse_arg(oper);
  result.push_back(oper);

  return result;
}

static std::vector<std::string>::iterator
expand(const std::vector<std::string> &p_oper, std::vector<std::string> &p_result, std::vector<std::string> &p_argv,
       std::vector<std::string>::iterator p_start)
{
  auto argvIter = p_start;
  while (argvIter != p_argv.end() && (*(argvIter))[0] != ']')
    {
      auto currentArg = (*(argvIter));
      if (currentArg[0] == '[')
        {
          fprintf(stderr, "cdo apply: brackets not allowed for command apply.\n");
          errState = ApplyStatus::BRACKET_USED;
          return argvIter;
        }
      if (currentArg[0] == '-')
        {
          const auto mod = get_module(currentArg.c_str());
          if (mod.get_stream_in_cnt() > 0)
            {
              fprintf(stderr, "cdo apply: operators with inputs are not allowed for command apply.\n");
              errState = ApplyStatus::OPER_WITH_INPUT_USED;
              return argvIter;
            }
        }
      for (auto &oper : p_oper) { p_result.push_back(oper); }
      p_result.push_back(*argvIter);
      ++argvIter;
    }

  if (argvIter == p_argv.end() || (*(argvIter))[0] != ']')
    {
      fprintf(stderr, "cdo apply: missing closing bracket.\n");
      errState = ApplyStatus::MISSING_CLOSING_BRACKET;
    }

  return argvIter;
}

static std::vector<std::string>
scan(std::vector<std::string> p_argv)
{
  std::vector<std::string> newArgv = {};
  if (p_argv[0].compare(0, strlen("-apply,"), "-apply,") == 0)
    {
      fprintf(stderr, "cdo apply: apply can not be first.\n");
      errState = ApplyStatus::APPLY_USED_FIRST;
      return p_argv;
    }
  for (auto argvIter = p_argv.begin(); argvIter < p_argv.end() && errState == ApplyStatus::OK; argvIter++)
    {
      std::string currentArgv = *argvIter;
      if (currentArgv.compare(0, 6, "-apply") == 0)
        {
          ++argvIter;
          if ((*(argvIter))[0] == '[')
            {
              const auto pos = currentArgv.find(',');
              if (pos != std::string::npos)
                {
                  auto parameter = currentArgv.substr(pos + 1);
                  if (parameter.empty())
                    {
                      fprintf(stderr, "cdo apply: missing argument for apply.\n");
                      errState = ApplyStatus::MISSING_ARG;
                      break;
                    }
                  else if (parameter[0] != '-')
                    {
                      fprintf(stderr, "cdo apply: missing pipe symbol in apply argument: %s\n", parameter.c_str());
                      errState = ApplyStatus::MISSING_PIPE_SYM;
                      break;
                    }
                  auto tokens = generate_tokens(parameter);
                  if (errState != ApplyStatus::OK) break;
                  argvIter = expand(tokens, newArgv, p_argv, ++argvIter);
                  if (errState != ApplyStatus::OK) break;
                }
              else
                {
                  fprintf(stderr, "cdo apply: missing argument for apply.\n");
                  errState = ApplyStatus::MISSING_ARG;
                  break;
                }
            }
          else
            {
              fprintf(stderr, "cdo apply: Missing bracket after apply, apply can only be used with [].\n");
              errState = ApplyStatus::MISSING_BRACKET;
            }
        }
      else { newArgv.push_back(*argvIter); }
    }

  return newArgv;
}

std::vector<std::string>
expand_apply(const std::vector<std::string> &p_argv, ApplyStatus &expandSuccess)
{
  errState = ApplyStatus::OK;
  auto result = scan(p_argv);
  expandSuccess = errState;
  return result;
}
