/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/
#ifndef PIPE_H
#define PIPE_H

#include <mutex>
#include <condition_variable>

#include "field.h"

struct pipe_t
{
public:
  pipe_t();
  void pipe_init();

  int pipe_inq_vlist(int &vlistID);
  void pipe_def_vlist(int &target_vlistID, int new_vlistID);

  int pipe_inq_timestep(int p_tsID);
  void pipe_def_timestep(int p_vlistID, int tsID);

  int pipe_inq_record(int *varID, int *levelID);
  void pipe_def_record(int p_varId, int p_levelID);

  void pipe_write_record(double *p_data, size_t p_nmiss);
  void pipe_write_record(float *p_data, size_t p_nmiss);
  void pipe_write_record(Field *p_data, size_t p_nmiss);

  size_t pipe_read_record(int p_vlistID, double *data, size_t *nmiss);
  size_t pipe_read_record(int p_vlistID, float *data, size_t *nmiss);
  size_t pipe_read_record(int p_vlistID, Field *data, size_t *nmiss);

  size_t pipe_read_pipe_record(double *data, int vlistID, size_t *p_nmiss);
  size_t pipe_read_pipe_record(float *data, int vlistID, size_t *p_nmiss);

  void pipe_set_name(int processID, int inputIDX);
  void close();

  bool EOP;
  bool usedata;
  bool hasdata;

  int varID, levelID;
  int recIDr, recIDw, tsIDr, tsIDw;

  size_t nmiss;
  int nrecs;

  bool data_is_float;
  double *data_d;
  float *data_f;

  std::mutex m_mutex;
  std::condition_variable tsDef_cond, tsInq_cond, vlistDef_cond, isClosed_cond;
  std::condition_variable recDef_cond, recInq_cond;
  std::condition_variable write_cond, read_cond;

  std::string name;

private:
  void wait_for_read();
};

#endif /* PIPE_H */
