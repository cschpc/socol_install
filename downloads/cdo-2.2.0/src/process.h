/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida
          Oliver Heidmann

*/
#ifndef PROCESS_H
#define PROCESS_H

#include "cdoStream.h"
#include "modules.h"

#include <vector>
#include <iostream>
#include <string>
#include <set>

constexpr int MAX_PROCESS = 65536;
constexpr int MAX_OPERATOR = 128;

enum class ProcessStatus
{
  Ok = 0,
  UnlimitedIOCounts = -1,
  MissInput = -2,
  MissOutput = -3,
  TooManyStreams = -4,
  TooFewStreams = -5,
};

class oper_t
{
public:
  int f1 = 0;
  int f2 = 0;
  std::string name;
  const char *enter = nullptr;

  oper_t() {}
  oper_t(int _f1, int _f2, const char *_name, const char *_enter) : f1(_f1), f2(_f2), name(_name), enter(_enter) {}
};

class Process
{
public:
  int m_ID;
  int m_posInParent;
#ifdef HAVE_LIBPTHREAD
  pthread_t threadID;
#endif
  short nchild = 0;
  std::vector<std::shared_ptr<Process>> childProcesses;
  std::vector<std::shared_ptr<Process>> parentProcesses;
  std::vector<CdoStreamID> inputStreams;
  std::vector<CdoStreamID> outputStreams;
  int nChildActive = 0;

  double startTime;

  int nvars = 0;

  int ntimesteps = 0;
  int m_streamCnt = 0;
  std::string m_operatorCommand = "UNINITALIZED";
  std::string operatorName;
  std::string m_obase;
  char prompt[64];

  short m_noper = 0;
  bool m_isActive;

  module_t m_module;
  std::vector<std::string> m_oargv;
  oper_t oper[MAX_OPERATOR];

  Process(int p_ID, const std::string &p_operatorNamme, const std::vector<std::string> &operatorArguments);

  pthread_t run();
  ProcessStatus m_status = ProcessStatus::Ok;

  /**
   * returns the number of in streams this process currently has.
   **/
  int get_stream_cnt_in();
  /**
   * returns the number of out streams this process currently has.
   */
  int get_stream_cnt_out();

  /**
   * Adds a Process as child and creates and adds a new pipe stream.
   */
  void add_child(const std::shared_ptr<Process> &child_process);
  /**
   * Adds a Process as parent and adds the parents input stream to out streams.
   */
  void add_parent(const std::shared_ptr<Process> &parent_process);
  /**
   * Compares the wanted and current stream counts.
   * @return if the wanted count is -1 this function always returns false.
   * Are the current and wanted stream counts equal 1 and if they differ false.
   */
  bool has_hall_inputs();
  /**
   * Adds and creates a new file pstream to the in streams
   */
  void add_file_in_stream(const std::string &file);
  /**
   * Adds and creates a new file pstream to the out streams
   */
  void add_file_out_stream(const std::string &file);
  /**
   * Adds and creates a new pipe pstream to the in streams
   */
  void add_pipe_in_stream();
  /**
   * Adds and creates a new file pstream to the out streams
   */
  void add_pipe_out_stream();
  /**
   * Adds an operator to the process
   */
  int operator_add(const char *name, int f1, int f2, const char *enter);
  /**
   * returns the operatorID of the currently in use operator
   */
  int get_operator_id();
  void set_inactive();
  const char *inq_prompt() const;
  void print_processed_values();
  void print_benchmarks(double p_walltime, const char *p_memstring);
  void check_stream_cnt();
  void validate();
  void handle_process_err();

  bool has_out_stream(const CdoStreamID p_streamPtr);
  bool has_in_stream(const CdoStreamID p_streamPtr);
  bool has_no_pipes();
  size_t inq_nvals();
  const char *get_out_stream_name();
  size_t get_oper_argc();
  std::string get_argv(int idx);
  const std::string get_obase();
  void
  set_obase(const std::string &obase)
  {
    m_obase = obase;
  }  // TODO into cc
  bool input_is_variable();
  int get_id();

  void init_process(const std::string &p_operatorName, const std::vector<std::string> &p_arguments);
  void close_streams();

private:
  Process();
  void def_prompt();
};

int get_process_num();
void set_process_num(int p_num);

#endif /* PROCESS_H */
