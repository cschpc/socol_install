#ifndef CDOSTREAM_H
#define CDOSTREAM_H

#ifdef HAVE_CONFIG_H
#include "config.h" /* _FILE_OFFSET_BITS influence off_t */
#endif

#include <sys/types.h> /* off_t */
#include <vector>
#include <memory>

#include "datarangelist.h"
#include "field.h"

class FileStream;

class CdoStream;

using CdiStreamID = int;
using CdoStreamID = std::shared_ptr<CdoStream>;
#define CDO_STREAM_UNDEF nullptr

class CdoStream
{
public:
  // Constructors
  virtual int open_read() = 0;
  virtual int open_write(int p_filetype) = 0;
  virtual int open_append() = 0;

  virtual int inq_vlist() = 0;
  virtual void def_vlist(int p_vlistID) = 0;

  virtual void inq_record(int *varID, int *levelID) = 0;
  virtual void defRecord(int varID, int levelID) = 0;

  virtual void read_record(float *p_data, size_t *nmiss) = 0;
  virtual void read_record(double *p_data, size_t *nmiss) = 0;
  virtual void read_record(Field *p_field, size_t *nmiss) = 0;

  virtual void write_record(float *p_data, size_t nmiss) = 0;
  virtual void write_record(double *p_data, size_t nmiss) = 0;
  virtual void write_record(Field *p_field, size_t nmiss) = 0;

  virtual void copyRecord(CdoStreamID dest) = 0;

  virtual int inq_timestep(int tsID) = 0;
  virtual void def_timestep(int tsID) = 0;

  virtual int inqFileType() = 0;
  virtual int inqByteorder() = 0;

  virtual void close() = 0;

  virtual size_t getNvals() = 0;

  int get_id();
  int getTsID();
  int getVarID();
  int getVlistID();

  int m_cdoStreamID = -1;  // aka the id of the pstream
  int m_filetype = CDI_UNDEFID;
  size_t m_nvals = 0;
  bool isopen = false;

  std::string m_name;
  std::vector<Datarange> m_datarangelist;

  int m_vlistID = -1;
  int m_tsID = -1;
  int m_varID = -1;  // next varID defined with streamDefVar

  // to be removed or to be moved to FileStream! // some operators need some refactoring for these to be able to be moved
  int m_fileID = 0;
  bool ispipe = false;

protected:
  CdoStream();
  virtual ~CdoStream() = 0;
  void defDatarangeList(int p_vlistID);

private:
};

#endif /* CDOSTREAM_H */
