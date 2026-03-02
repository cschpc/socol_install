/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Oliver Heidmann

*/

#include "processManager.h"
#include "process.h"
#include "cdo_output.h"
#include "cdo_options.h"
#include "fileStream.h"

#include <stack>
#include <mutex>

static std::mutex processMutex;

static std::string parse_err_msg = "";

static const int IS_OBASE = -1;

void
ProcessManager::buildProcessTree(std::vector<std::shared_ptr<Node>> roots)
{
  Debug(PROCESS,"Building process Tree");
  std::shared_ptr<Node> node = roots[0]->isFile ? roots[0]->children[0] : roots[0];

  auto first_process = create_process(node->oper, split_args(node->arguments));
  if (node->numOut() == IS_OBASE)
    {
      Debug(PROCESS,"Setting obase for %s", node->oper);
      first_process->set_obase(roots[0]->oper);
    }
  else if (node->numOut() > 0)
    {
      for (const auto &n : roots)
        {
          Debug(PROCESS,"adding out files to %s", node->oper);
          first_process->add_file_out_stream(n->oper);
        }
    }

  for (const auto &c : node->children)
    {
      if (c->isFile)
        {
          Debug(PROCESS,"Adding file in stream: %s", c->oper);
          first_process->add_file_in_stream(c->oper);
        }
      else
        {
          auto c_ptr = build_node(c);
          first_process->add_child(c_ptr);
          c_ptr->add_parent(first_process);
        }
    }

  set_process_num(m_processes.size());
  FileStream::enableTimers(m_processes.size() == 1 && Threading::ompNumThreads == 1);
}

std::shared_ptr<Process>
ProcessManager::build_node(std::shared_ptr<Node> parent_node)
{
  Debug(PROCESS,"Building process for %s", parent_node->oper);
  auto parent_process = create_process(parent_node->oper, split_args(parent_node->arguments));
  for (auto child_node : parent_node->children)
    {
      if (child_node->isFile)
        {
          Debug(PROCESS,"Adding file in stream: %s", child_node->oper);
          parent_process->add_file_in_stream(child_node->oper);
        }
      else
        {
          Debug(PROCESS,"Building Process for %s", child_node->oper);
          auto child_process = build_node(child_node);
          parent_process->add_child(child_process);
          child_process->add_parent(parent_process);
        }
    }
  return parent_process;
}

void
ProcessManager::run_processes()
{
  for (auto &idProcessPair : m_processes)
    {
      if (idProcessPair.first)
        {
          /*TEMP*/
          if (!Options::silentMode && (cdo::stdoutIsTerminal || Options::cdoVerbose))
            {
              // MpMO::Print(Green("%s: ") + "Process started", idProcessPair.second->prompt);
              set_text_color(stdout, GREEN);
              fprintf(stdout, "%s: ", idProcessPair.second->prompt);
              reset_text_color(stdout);
              fprintf(stdout, "Process started\n");
            }
          m_threadIDs.push_back(idProcessPair.second->run());
        }
    }
  m_threadIDs.push_back(pthread_self());
  // MpMO::PrintCerr(Green("%s: ") + "xProcess started", get_process_from_id(0).inq_prompt());
  get_process_from_id(0)->m_module.func(get_process_from_id(0).get());
}

void
ProcessManager::kill_processes()
{
  for (auto threadID : m_threadIDs)
    {
      if (threadID != pthread_self())
        {
          pthread_cancel(threadID);
          Debug(PROCESS_MANAGER, "process killed: %ld", threadID);
        }
    }
}

void
ProcessManager::clear_processes()
{
  Debug(PROCESS_MANAGER, "Deleting Processes");
  m_processes.clear();
  m_numProcesses = 0;
  m_numProcessesActive = 0;
}
const std::shared_ptr<Process>
ProcessManager::create_process(const std::string &operatorName, const std::vector<std::string> &arguments)
{
  auto processID = m_numProcesses++;
  if (processID >= MAX_PROCESS) cdo_abort("Limit of %d processes reached!", MAX_PROCESS);
  auto success = m_processes.insert(std::make_pair(processID, std::make_shared<Process>(processID, operatorName, arguments)));
  if (!success.second) cdo_abort("Process %d could not be created", processID);
  m_numProcessesActive++;
  return success.first->second;
}

int
ProcessManager::get_num_processes(void)
{
  std::scoped_lock lock(processMutex);
  int pnums = m_processes.size();
  return pnums;
}

int
ProcessManager::get_num_active_processes(void)
{
  std::scoped_lock lock(processMutex);
  int pnums = m_numProcessesActive;
  return pnums;
}
std::vector<std::string>
ProcessManager::get_operator_argv(std::string operatorArguments)
{
  std::vector<std::string> argument_vector;
  Debug(PROCESS && strchr(operatorArguments.c_str(), ',') != nullptr, "Setting operator arguments: %s", operatorArguments);

  constexpr char delimiter = ',';

  auto pos = operatorArguments.find(delimiter);
  if (pos != std::string::npos)
    {
      // remove operator name
      operatorArguments.erase(0, pos + 1);

      while ((pos = operatorArguments.find(delimiter)) != std::string::npos)
        {
          argument_vector.push_back(operatorArguments.substr(0, pos));
          Debug(PROCESS,"added argument %s", argument_vector.back());
          operatorArguments.erase(0, pos + 1);
        }
      argument_vector.push_back(operatorArguments);
    }
  return argument_vector;
}

std::vector<std::string>
ProcessManager::split_args(std::string operatorArguments)
{
  if (operatorArguments.empty()) return {};
  Debug(PROCESS, "Setting operator arguments: '%s'", operatorArguments);
  std::vector<std::string> argument_vector = {};
  constexpr char delimiter = ',';
  size_t pos;
  while ((pos = operatorArguments.find(delimiter)) != std::string::npos)
    {
      argument_vector.push_back(operatorArguments.substr(0, pos));
      Debug(PROCESS,"added argument %s", argument_vector.back());
      operatorArguments.erase(0, pos + 1);
    }
  argument_vector.push_back(operatorArguments);
  return argument_vector;
}

const std::shared_ptr<Process> &
ProcessManager::get_process_from_id(int p_processID)
{
  std::scoped_lock lock(processMutex);

  auto process = m_processes.find(p_processID);
  if (process == m_processes.end()) cdo_abort("Process with ID: %d not found", p_processID);

  return process->second;
}
