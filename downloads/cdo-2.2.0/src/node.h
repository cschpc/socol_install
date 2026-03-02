#ifndef CDO_NODE_HPP
#define CDO_NODE_HPP

#include <vector>
#include <string>
#include <memory>
#include <iostream>
#include "modules.h"

static std::string errmsg_node_to_many_inputs = "To many inputs";
static std::string errmsg_node_no_output= "Operator has no output, cannot be used with pipes unless used first";
static std::string errmsg_node_unassigned  ="Could not be assigned, leftover input";
static std::string errmsg_node_file_to_file = "Attempted to attach file to file";
static std::string errmsg_node_only_accepts_files = "Operator cannot be piped into an operator that takes only files";

class Node
{
public:
  std::vector<std::string>::const_iterator iter;
  const std::string oper;
  const std::string arguments;
  module_constraints constraints;
  std::vector<std::shared_ptr<Node>> children = {};

  Node(std::vector<std::string>::const_iterator p_iter, const std::string &p_operName, const std::string &p_args,
       module_constraints p_constraints);
  Node(std::vector<std::string>::const_iterator p_iter, bool p_isOutFile);
  explicit Node(Node *p_nodePtr);
  std::shared_ptr<Node> copy();

  // Ready to be returned and process
  bool has_missing_input();
  // Done in terms of beein on the stack
  bool is_done();
  bool is_temporary_leaf();
  bool is_leaf();

  void add_leaf(std::shared_ptr<Node> &p_newNode);
  void append(std::vector<std::shared_ptr<Node>> &p_node);
  void append(std::shared_ptr<Node> &p_node);

  int
  numMaxChildren()
  {
    return constraints.streamInCnt;
  }
  int
  numOut()
  {
    return constraints.streamOutCnt;
  }
  PositionRestrictions
  get_restriction()
  {
    return constraints.pos_restriction;
  }
  const bool isFile = false;
  const bool isOutFile = false;

  std::string to_string();
};


#endif
