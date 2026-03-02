#include "node.h"
#include "cdo_output.h"
#include "cdo_node_attach_exception.h"


Node::Node(std::vector<std::string>::const_iterator p_iter, const std::string &p_oper, const std::string &args,
           module_constraints p_constraints)
    : iter(p_iter), oper(p_oper), arguments(args), constraints(p_constraints), isFile(false), isOutFile(false)
{
}

Node::Node(std::vector<std::string>::const_iterator p_iter, bool p_isOutFile = false)
    : iter(p_iter), oper(*p_iter), arguments(""),
      constraints(
          { p_isOutFile ? (short) 1 : (short) 0, p_isOutFile ? (short) 0 : (short) 1, PositionRestrictions::NoRestriction }),
      isFile(true), isOutFile(p_isOutFile)
{
}

Node::Node(Node *node_ptr)
    : iter(node_ptr->iter), oper(node_ptr->oper), arguments(node_ptr->arguments), constraints(node_ptr->constraints),
      isFile(node_ptr->isFile), isOutFile(node_ptr->isOutFile)
{
}

bool
Node::has_missing_input()
{
  if (children.size() == 0 && constraints.streamInCnt != 0) return true;
  if (isFile || children.size() == (size_t) constraints.streamInCnt || constraints.streamInCnt == -1) return false;
  return true;
}

bool
Node::is_done()
{
  bool done = false;
  if (isFile && !isOutFile) { done = true; }
  else if ((int) children.size() == constraints.streamInCnt)
    {
      done = true;
    }
  Debug(CDO_NODE, "%s is done: %s", oper, done ? "true" : "false");
  return done;
}

void
Node::append(std::shared_ptr<Node> &node)
{

  Debug(CDO_NODE, "appending  %s to %s", node->oper, oper);
  if (isFile && !isOutFile && node->isFile) { throw NodeAttachException(node, errmsg_node_file_to_file); }
  if (isOutFile && is_done()) { throw NodeAttachException(node, errmsg_node_unassigned); }
  if (constraints.streamInCnt >= 0 && (int) children.size() == constraints.streamInCnt)
    {
      throw NodeAttachException(iter, errmsg_node_to_many_inputs);
    }
  if (node->numOut() == 0)
    {
      throw NodeAttachException(node, errmsg_node_no_output);
    }
  if (get_restriction() == PositionRestrictions::FilesOnly && !node->isFile)
    {
      size_t start = ((*iter)[0] == '-') ? 1 : 0;
      std::string name = (*iter).substr(start, (*iter).size() - 0);
      throw NodeAttachException(node, errmsg_node_only_accepts_files);
    }
  children.push_back(node);
}

void
Node::append(std::vector<std::shared_ptr<Node>> &n)
{
  for (auto &x : n) append(x);
}

bool
Node::is_temporary_leaf()
{
  return (children.size() != (size_t) constraints.streamInCnt) && !isFile && (constraints.streamInCnt != 0);
}

void
Node::add_leaf(std::shared_ptr<Node> &new_node)
{
  Debug(CDO_NODE, "add_leaf of node: %s", new_node->oper);
  if (is_temporary_leaf())
    {
      Debug(CDO_NODE, "adding leaf to %s", oper);
      append(new_node);
    }
  else
    {
      for (auto &c : children) { c->add_leaf(new_node); }
    }
}
std::shared_ptr<Node>
Node::copy()
{
  auto copiedNode = std::make_shared<Node>(this);
  copiedNode->children.clear();
  for (auto &child : children) { copiedNode->children.push_back(child->copy()); }
  return copiedNode;
}

std::string
Node::to_string()
{
  std::string r = oper;
  if (!arguments.empty()) { r += "," + arguments; }
  if (children.size() > 0) r += " [";
  for (auto &c : children) { r += " " + c->to_string(); }
  if (children.size() > 0) r += " ]";
  return r;
}
