/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida
          Oliver Heidmann

*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_PTHREAD_H
#include <pthread.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#include <stdio.h>
#include <iostream>
#include <string>
#include <sstream>

#include "process.h"
#include "cdo_wtime.h"
#include "util_string.h"
#include "cdo_options.h"
#include "fileStream.h"
#include "pipeStream.h"

static int processNum = 0;

int
get_process_num()
{
  return processNum;
}

void
set_process_num(int p_num)
{
  processNum = p_num;
}

Process::Process(int p_ID, const std::string &p_operatorName, const std::vector<std::string> &p_arguments)
    : m_ID(p_ID), operatorName(p_operatorName)
{
  m_isActive = true;
  init_process(p_operatorName, p_arguments);
}

bool
Process::input_is_variable()
{
  return (m_module.get_stream_in_cnt() == -1);
}

void
Process::init_process(const std::string &p_operatorName, const std::vector<std::string> &p_arguments)
{
  startTime = cdo_get_wtime();
#ifdef HAVE_LIBPTHREAD
  threadID = pthread_self();
#endif
  m_oargv = p_arguments;
  operatorName = get_original(p_operatorName);

  m_module = get_module(p_operatorName);

  def_prompt();  // has to be called after get operatorName
}

int
Process::get_stream_cnt_in()
{
  return inputStreams.size();
}

int
Process::get_stream_cnt_out()
{
  return outputStreams.size();
}

void
Process::def_prompt()
{
  if (m_ID == 0)
    std::snprintf(prompt, sizeof(prompt), "%s    %s", cdo::progname, operatorName.c_str());
  else
    std::snprintf(prompt, sizeof(prompt), "%s(%d) %s", cdo::progname, m_ID, operatorName.c_str());
}

const char *
Process::inq_prompt() const
{
  return prompt;
}

void
Process::handle_process_err()
{
  switch (m_status)
    {
    case ProcessStatus::UnlimitedIOCounts:
      {
        cdo_abort("I/O stream counts unlimited no allowed!");
        break;
      }
    case ProcessStatus::MissInput:
      {
        cdo_abort("Input streams missing!");
        break;
      }
    case ProcessStatus::MissOutput:
      {
        cdo_abort("Output streams missing!");
        break;
      }
    case ProcessStatus::TooManyStreams:
    case ProcessStatus::TooFewStreams:
      {
        const auto inCnt = m_module.get_stream_in_cnt();
        auto outCnt = m_module.get_stream_out_cnt();
        const bool lobase = (outCnt == -1);
        if (lobase) outCnt = 1;

        std::string caseCount = (m_status == ProcessStatus::TooManyStreams) ? "many" : "few";
        std::string pluralIn = (inCnt > 1) ? "s" : "";
        std::string pluralOut = (outCnt > 1) ? "s" : "";

        std::stringstream errMsg;
        errMsg << "Too " << caseCount << " streams specified! Operator " << operatorName << " needs " << inCnt << " input stream"
               << pluralIn << " and " << outCnt << " output " << (lobase ? "basename" : "stream") << pluralOut << "!";
        cdo_abort(errMsg.str());
        break;
      }
    case ProcessStatus::Ok: break;
    }
}

void
Process::validate()
{
  check_stream_cnt();
  if (m_status != ProcessStatus::Ok) handle_process_err();
}

void
Process::check_stream_cnt()
{
  int streamCnt = 0;

  int wantedStreamInCnt = m_module.get_stream_in_cnt();
  int wantedStreamOutCnt = m_module.get_stream_out_cnt();

  int streamInCnt0 = wantedStreamInCnt;

  bool obase = false;
  if (wantedStreamOutCnt == -1)
    {
      wantedStreamOutCnt = 1;
      obase = true;
    }

  if (wantedStreamInCnt == -1 && wantedStreamOutCnt == -1) { m_status = ProcessStatus::UnlimitedIOCounts; }

  // printf(" wantedStreamInCnt,wantedStreamOutCnt %d %d\n",
  // wantedStreamInCnt,wantedStreamOutCnt);
  else if (wantedStreamInCnt == -1)
    {
      wantedStreamInCnt = m_streamCnt - wantedStreamOutCnt;
      if (wantedStreamInCnt < 1) m_status = ProcessStatus::MissInput;
    }

  else if (wantedStreamOutCnt == -1)
    {
      wantedStreamOutCnt = m_streamCnt - wantedStreamInCnt;
      if (wantedStreamOutCnt < 1) m_status = ProcessStatus::MissOutput;
    }
  else
    {
      // printf(" wantedStreamInCnt,wantedStreamOutCnt %d %d\n",
      // wantedStreamInCnt,wantedStreamOutCnt);

      streamCnt = wantedStreamInCnt + wantedStreamOutCnt;
      // printf(" streamCnt %d %d\n", m_streamCnt, streamCnt);

      if (m_streamCnt > streamCnt)
        m_status = ProcessStatus::TooManyStreams;
      else if (m_streamCnt < streamCnt && !obase)
        m_status = ProcessStatus::TooFewStreams;
      else if (wantedStreamInCnt > (int) inputStreams.size())
        m_status = ProcessStatus::TooFewStreams;

      else if (wantedStreamInCnt == 1 && streamInCnt0 == -1)
        m_status = ProcessStatus::Ok;
    }
}

bool
Process::has_hall_inputs()
{
  if (m_module.get_stream_in_cnt() == -1) return false;

  return (m_module.get_stream_in_cnt() == static_cast<short>(inputStreams.size()));
}

void
Process::set_inactive()
{
  m_isActive = false;
}

int
Process::operator_add(const char *name, int f1, int f2, const char *enter)
{
  const int operID = m_noper;

  if (operID >= MAX_OPERATOR) cdo_abort("Maximum number of %d operators reached!", MAX_OPERATOR);

  oper[m_noper] = { f1, f2, name, enter };

  m_noper++;

  return operID;
}

int
Process::get_operator_id()
{
  if (m_noper > 0)
    {
      for (int operID = 0; operID < m_noper; operID++)
        {
          if (operatorName == oper[operID].name) return operID;
        }

      cdo_abort("Operator not callable by this name! Name is: %s", operatorName);
    }

  cdo_abort("Operator not initialized!");

  return -1;
}

void
Process::add_file_in_stream(const std::string &file)
{
  inputStreams.push_back(std::make_shared<FileStream>(file));
  m_streamCnt++;
}

void
Process::add_file_out_stream(const std::string &file)
{
  if (file[0] == '-') { cdo_abort("Missing output file. Found an operator instead of filename: %s", file); }
  outputStreams.push_back(std::make_shared<FileStream>(file));
  m_streamCnt++;
}

void
Process::add_child(const std::shared_ptr<Process> &childProcess)
{
  childProcesses.push_back(childProcess);
  nchild = childProcesses.size();
  add_pipe_in_stream();
}

void
Process::add_pipe_in_stream()
{
#ifdef HAVE_LIBPTHREAD
  inputStreams.push_back(std::make_shared<PipeStream>(m_ID));
  m_streamCnt++;
#else
  cdo_abort("Cannot use pipes, pthread support not compiled in!");
#endif
}

void
Process::add_parent(const std::shared_ptr<Process> &parentProcess)
{
  parentProcesses.push_back(parentProcess);
  m_posInParent = parentProcess->inputStreams.size() - 1;
  add_pipe_out_stream();
}

void
Process::add_pipe_out_stream()
{
  outputStreams.push_back(parentProcesses[0]->inputStreams[m_posInParent]);
  m_streamCnt++;
}

pthread_t
Process::run()
{
  Debug(PROCESS, "starting new thread for process %d", m_ID);
  pthread_attr_t attr;
  auto status = pthread_attr_init(&attr);
  if (status) cdo_sys_error("pthread_attr_init failed for '%s'", operatorName.c_str());
  status = pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
  if (status) cdo_sys_error("pthread_attr_setdetachstate failed for '%s'", operatorName.c_str());
  /*
    param.sched_priority = 0;
    status = pthread_attr_setschedparam(&attr, &param);
    if ( status ) cdo_sys_error("pthread_attr_setschedparam failed for '%s'", newarg+1);
  */
  /* status = pthread_attr_setinheritsched(&attr, PTHREAD_EXPLICIT_SCHED); */
  /* if ( status ) cdo_sys_error("pthread_attr_setinheritsched failed for '%s'", newarg+1); */

  int pthreadScope;
  pthread_attr_getscope(&attr, &pthreadScope);

  /* status = pthread_attr_setscope(&attr, PTHREAD_SCOPE_PROCESS); */
  /* if ( status ) cdo_sys_error("pthread_attr_setscope failed for '%s'", newarg+1); */
  /* If system scheduling scope is specified, then the thread is scheduled against all threads in the system */
  /* pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM); */

  size_t stacksize = 0;
  status = pthread_attr_getstacksize(&attr, &stacksize);
  if (stacksize < 2097152)
    {
      stacksize = 2097152;
      pthread_attr_setstacksize(&attr, stacksize);
    }

  pthread_t thrID;
  const auto rval = pthread_create(&thrID, &attr, m_module.func, this);
  if (rval != 0)
    {
      errno = rval;
      cdo_sys_error("pthread_create failed for '%s'", operatorName.c_str());
    }

  m_isActive = true;

  return thrID;
}

// local helper function
extern "C" size_t getPeakRSS();
static void
get_max_memstring(char *p_memstring, size_t memstringLen)
{
  auto memmax = getPeakRSS();
  if (memmax)
    {
      size_t muindex = 0;
      const char *mu[] = { "B", "KB", "MB", "GB", "TB", "PB" };
      const size_t nmu = sizeof(mu) / sizeof(char *);
      while (memmax > 9999 && muindex < nmu - 1)
        {
          memmax /= 1024;
          muindex++;
        }
      std::snprintf(p_memstring, memstringLen, " %zu%s", memmax, mu[muindex]);
    }
}

void
Process::print_benchmarks(double p_walltime, const char *p_memstring)
{
  const auto numberOfUsedThreads = get_process_num();
  if (Options::test)
    fprintf(stdout, " [%.2fs%s %dthread%s]", p_walltime, p_memstring, numberOfUsedThreads, ADD_PLURAL(numberOfUsedThreads));
  else
    fprintf(stdout, " [%.2fs%s]", p_walltime, p_memstring);
}

void
Process::print_processed_values()
{
  set_text_color(stdout, GREEN);
  fprintf(stdout, "%s: ", prompt);
  reset_text_color(stdout);

  const auto nvals = inq_nvals();

  if (nvals > 0)
    {
      fprintf(stdout, "Processed %zu value%s from %d variable%s", nvals, ADD_PLURAL(nvals), nvars, ADD_PLURAL(nvars));
    }
  else if (nvars > 0) { fprintf(stdout, "Processed %d variable%s", nvars, ADD_PLURAL(nvars)); }

  if ((nvals || nvars) && ntimesteps > 0) fprintf(stdout, " over %d timestep%s", ntimesteps, ADD_PLURAL(ntimesteps));

  if (m_ID == 0)
    {
      char memstring[32] = { "" };
      get_max_memstring(memstring, sizeof(memstring));
      print_benchmarks(cdo_get_wtime() - startTime, memstring);
    }

  // if (nvars > 0 || nvals > 0 || ntimesteps > 0 || m_ID == 0) fprintf(stdout, ".");
  fprintf(stdout, "\n");
}

bool
Process::has_out_stream(const CdoStreamID p_streamID)
{
  for (const CdoStreamID &streamID : outputStreams)
    {
      if (streamID == p_streamID) return true;
    }
  return false;
}

bool
Process::has_in_stream(const CdoStreamID p_streamID)
{
  for (const CdoStreamID &streamID : inputStreams)
    {
      if (streamID == p_streamID) return true;
    }
  return false;
}

size_t
Process::inq_nvals()
{
  size_t nvals = 0;
  for (size_t i = 0; i < inputStreams.size(); ++i)
    {
      Debug(PROCESS, "Inquiring nvals from instream %s", inputStreams[i]->m_name);
      nvals += inputStreams[i]->getNvals();
    }
  return nvals;
}

bool
Process::has_no_pipes()
{
  return (childProcesses.size() == 0);
}

const char *
Process::get_out_stream_name()
{
  return outputStreams[0]->m_name.c_str();
}

size_t
Process::get_oper_argc()
{
  return m_oargv.size();
}

std::string
Process::get_argv(int p_idx)
{
  if (!(p_idx > (int) get_oper_argc() && p_idx > 0))
    cdo_abort("Process Argv not found. Idx: %d, Process argc: %d", p_idx, m_oargv.size());

  return m_oargv[p_idx];
}

const std::string
Process::get_obase()
{
  return m_obase;
}
void
Process::close_streams()
{
  for (auto s : inputStreams) { s->close(); }
  for (auto s : outputStreams) { s->close(); }
}

int
Process::get_id()
{
  return m_ID;
}
