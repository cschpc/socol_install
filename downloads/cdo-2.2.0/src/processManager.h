#ifndef PROCESS_MANAGER_H
#define PROCESS_MANAGER_H

// Stdlib includes
#include <map>
#include <set>
#include <string>
#include <vector>
#include <string>
#include <memory>

#include <pthread.h>
#include "node.h"
#include "modules.h"

// cdo includes

// Froward declarations
class Process;

// Error codes
enum class ParseStatus
{
  Ok = 0,
  OpenBracketMissing = -1,
  ClosingBracketMissing = -2,
  UnprocessedInput = -3,
  MissingOutFile = -4,
  MissingObase = -5,
  OperatorNotFirst = -6,
  FileIsInAndOutput = -7

};
class ProcessManager
{

private:
  std::map<int, std::shared_ptr<Process>> m_processes;
  std::vector<pthread_t> m_threadIDs;

  int m_numProcesses = 0;
  int m_numProcessesActive = 0;

  std::vector<std::string> get_operator_argv(std::string operatorArguments);
  std::vector<std::string> split_args(std::string operatorArguments);
  const std::shared_ptr<Process> create_process(const std::string &operatorName, const std::vector<std::string> &arguments);

public:
  void run_processes();
  void kill_processes();
  void clear_processes();
  int get_num_processes();
  int get_num_active_processes();
  const std::shared_ptr<Process> &get_process_from_id(int p_processID);

  void buildProcessTree(std::vector<std::shared_ptr<Node>> root);

  std::shared_ptr<Process> build_node(std::shared_ptr<Node> ptr);
};

#endif
