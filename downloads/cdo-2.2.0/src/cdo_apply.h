#ifndef CDO_APPLY_H
#define CDO_APPLY_H

#include <vector>

/*
 * this feature of cdo allows the use of the -apply keyword.
 * When used all files that are inside the ][ brackets of apply are prepended with given operator.
 */

/*
 * checks for:
 *  missing bracket at pos after -apply -> error
 *  checks for new bracket -> error
 *  checks for no closing bracket -> error
 *  checks for apply not beeing at the first position in the argument string -> error
 *  checks for missing argument for apply
 *  checks for operators that are used within the apply brackets -> error
 */

/*
 * supports:
 * operators with number of inputs = 1 (no others)
 * operators with arguments: written as -apply,-OPERNAME,arguments
 */
enum class ApplyStatus : int
{
  OK = 0,
  ARG_TOO_MANY_OUT = 1,
  ARG_NO_INPUT = 2,
  ARG_VARIABLE_INPUT = 3,
  BRACKET_USED = 4,
  OPER_WITH_INPUT_USED = 5,
  MISSING_CLOSING_BRACKET = 6,
  APPLY_USED_FIRST = 7,
  MISSING_ARG = 8,
  MISSING_BRACKET = 9,
  MISSING_PIPE_SYM = 10
};

std::vector<std::string> expand_apply(const std::vector<std::string> &p_argv, ApplyStatus &p_status);

#endif
