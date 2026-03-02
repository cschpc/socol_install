#ifndef FILESTREAM_H
#define FILESTREAM_H

#include <string>
#include "cdoStream.h"

class FileStream : public CdoStream
{
public:
  // Constructors
  explicit FileStream(const std::string &p_fileStream);
  // ---

  // CdoStream Interface functions
  int open_read() override;
  int open_write(int p_filetype) override;
  int open_append() override;

  int inq_vlist() override;
  void def_vlist(int p_vlistID) override;

  void inq_record(int *varID, int *levelID) override;
  void defRecord(int varID, int levelID) override;

  void read_record(float *p_data, size_t *nmiss) override;
  void read_record(double *p_data, size_t *nmiss) override;
  void read_record(Field *p_field, size_t *nmiss) override;

  void write_record(float *p_data, size_t nmiss) override;
  void write_record(double *p_data, size_t nmiss) override;
  void write_record(Field *p_field, size_t nmiss) override;

  void copyRecord(CdoStreamID p_fileStream) override;

  int inq_timestep(int tsID) override;
  void def_timestep(int tsID) override;

  int inqFileType() override;
  int inqByteorder() override;

  void close() override;

  size_t getNvals() override;
  // ---

  // FileStreamOnly
  int getFileID();
  static void
  enableTimers(const bool p_enable)
  {
    FileStream::TimerEnabled = p_enable;
  }
  static bool
  timersEnabled()
  {
    return FileStream::TimerEnabled;
  }
  // ---

protected:
  static bool TimerEnabled;

private:
  FileStream() = delete;
  std::string m_filename;
  void checkDatarange(int varID, double *array, size_t nmiss);
};

#endif
