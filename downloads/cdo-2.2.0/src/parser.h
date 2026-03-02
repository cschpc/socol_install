#ifndef PARSER_H
#define PARSER_H

#include <vector>
#include <stack>
#include "node.h"
#include "cdo_syntax_error.h"

namespace Parser
{

//'Regular' parser messages
static std::string errmsg_multiple_variable = "Using two operators with variable input without using sub groups [] is not allowed";
static std::string errmsg_missing_outputs = "Missing outputs";
static std::string errmsg_missing_inputs = "Missing inputs";
static std::string errmsg_unprocessed_inputs = "Unprocessed Inputs";
static std::string errmsg_keyword_output = "Keywords cannot be used as file names";

//Subgroup errors
static std::string errmsg_mixed_input = "Mixing of normal inputs and subgroups is not allowed";
static std::string errmsg_missing_sub_group = "Closing bracket without open subgroup";
static std::string errmsg_empty_subgroup = "Empty Subgroup";
static std::string errmsg_bracket_not_closed = "Bracket not closed";
static std::string errmsg_malformed_subgroup = "Malformed Subgroup";

//Apply error messages
static std::string errmsg_only_1_to_1_operators = "Only operators with a single in and output allowed";
static std::string errmsg_apply_missing_argument = "Missing arguments";
static std::string errmsg_apply_multiple_roots = "Apply can only process chains with a single in and out put";
static std::string errmsg_apply_requires_bracket = "Apply requires brackets";
static std::string errmsg_apply_no_inputs = "Apply content has no available free inputs";
static std::string errmsg_apply_in_first_pos = "Apply can not be in first position";

std::string err_msg_oper_not_found(std::string name);

std::vector<std::shared_ptr<Node>> run(std::vector<std::string> &p_argv);
std::vector<std::shared_ptr<Node>> parse(std::vector<std::string> p_argv, const char *(*context)(void) );
std::vector<std::shared_ptr<Node>> _parse(std::vector<std::string> p_argv, const char *(*context)(void) );

namespace Util
{

std::string result_to_string(std::vector<std::shared_ptr<Node>> p_roots, std::string p_text = "returning: ");

std::string build_err_msg(std::vector<std::string> &p_argv, const std::vector<std::string>::const_iterator &iter,
                          const std::string &prompt, int cdo_abort_prompt_spacing = 10);
}  // namespace Util

struct MissingOutFileException : public std::invalid_argument
{
  explicit MissingOutFileException(const std::string &p_msg) : std::invalid_argument(p_msg) {}
};

}  // namespace Parser

#endif
