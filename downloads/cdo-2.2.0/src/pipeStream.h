/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida
          Oliver Heidmann

*/
#ifndef PIPESTREAM_H
#define PIPESTREAM_H

#include "cdoStream.h"
#include "pipe.h"

#ifdef HAVE_LIBPTHREAD

class FileStream;  // Predeclaration only for copyRecord(...)

class PipeStream : public CdoStream
{
public:
  // Constructors
  explicit PipeStream(int p_processID);
  // ---

  // CdoStream Interface functions
  int open_read();
  int open_write(int p_filetype);
  int open_append();

  int inq_vlist();
  void def_vlist(int p_vlistID);

  void inq_record(int *varID, int *levelID);
  void defRecord(int varID, int levelID);

  void read_record(float *p_data, size_t *nmiss);
  void read_record(double *p_data, size_t *nmiss);
  void read_record(Field *p_field, size_t *nmiss);

  void write_record(float *p_data, size_t nmiss);
  void write_record(double *p_data, size_t nmiss);
  void write_record(Field *p_field, size_t nmiss);

  void copyRecord(CdoStreamID p_fileStream);

  int inq_timestep(int tsID);
  void def_timestep(int tsID);

  int inqFileType();
  int inqByteorder();

  void close();

  size_t getNvals();
  // ---

  // FileStreamOnly
  // ---

private:
  PipeStream() = delete;
  std::shared_ptr<pipe_t> m_pipe = std::make_shared<pipe_t>();
  pthread_t rthreadID;  // read  thread ID
  pthread_t wthreadID;  // write thread ID
  void waitForPipe();
};
#endif

#endif
