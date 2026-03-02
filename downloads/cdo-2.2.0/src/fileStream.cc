/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#ifdef HAVE_CONFIG_H
#include "config.h" /* HAVE_NC4HDF5_THREADSAFE */
#endif

#include <sys/stat.h> /* stat */
#include <cdi.h>

#include "fileStream.h"

#include "cdi_lockedIO.h"
#include "cdo_options.h"
#include "cdo_output.h"
#include "cdo_default_values.h"
#include "cdo_history.h"
#include "cdo_query.h"
#include "process.h"

#include "timer.h"
#include "commandline.h"

bool FileStream::TimerEnabled = false;
#ifndef HAVE_NC4HDF5_THREADSAFE
static auto inputFileTypeIsNetCDF4 = false;
#endif

FileStream::FileStream(const std::string &p_fileName) : m_filename(p_fileName)
{
  m_name = p_fileName;
  m_fileID = CDI_UNDEFID;
}

int
FileStream::open_read()
{
  int fileID = -1;

  open_lock();

  if (m_filename.size() > 6 && m_filename.rfind("query:", 0) == 0)
    {
      CdiQuery *query = cdiQueryCreate();
      auto path = set_query_parameter(m_filename.substr(6), query);
      if (Options::cdoVerbose) cdiQueryPrint(query);

      fileID = streamOpenReadQuery(path.c_str(), query);
    }
  else
    {
      fileID = streamOpenRead(m_filename.c_str());
    }

  if (fileID < 0) cdi_open_error(fileID, "Open failed on >%s<", m_filename.c_str());
  isopen = true;

  m_filetype = streamInqFiletype(fileID);
  if (CdoDefault::FileType == CDI_UNDEFID) CdoDefault::FileType = m_filetype;
  m_fileID = fileID;

#ifndef HAVE_NC4HDF5_THREADSAFE
  if (m_filetype == CDI_FILETYPE_NC4 || m_filetype == CDI_FILETYPE_NC4C) inputFileTypeIsNetCDF4 = true;
#endif

  Debug(FILE_STREAM, "Set number of worker to %d", Options::numStreamWorker);
  if (Options::numStreamWorker > 0) streamDefNumWorker(fileID, Options::numStreamWorker);

  open_unlock();

  return fileID;
}

int
FileStream::open_write(int p_filetype)
{
  if (Options::cdoInteractive)
    {
      struct stat stbuf;
      auto rstatus = stat(m_name.c_str(), &stbuf);
      // If permanent file already exists, query user whether to overwrite or exit
      if (rstatus != -1) query_user_exit(m_name.c_str());
    }
  if (p_filetype == CDI_UNDEFID) p_filetype = CDI_FILETYPE_GRB;

#ifndef HAVE_NC4HDF5_THREADSAFE
  auto outputFileTypeIsNetCDF4 = (p_filetype == CDI_FILETYPE_NC4 || p_filetype == CDI_FILETYPE_NC4C);
  if (inputFileTypeIsNetCDF4 && outputFileTypeIsNetCDF4 && get_process_num() > 1 && Threading::cdoLockIO == false)
    {
      cdo_warning("Using a non-thread-safe NetCDF4/HDF5 library in a multi-threaded environment may lead to erroneous results!");
      cdo_warning("Use a thread-safe NetCDF4/HDF5 library or the CDO option -L to avoid such errors.");
    }
#endif

  // TODO FIX THIS: if (FileStream::timersEnabled()) timer_start(timer_write);

  open_lock();
  auto fileID = streamOpenWrite(m_filename.c_str(), p_filetype);
  open_unlock();

  // TODO FIX THIS: if(FileStream::timersEnabled()) timer_stop(timer_write);
  if (fileID < 0) cdi_open_error(fileID, "Open failed on >%s<", m_name.c_str());
  isopen = true;

  if (CdoDefault::Byteorder != CDI_UNDEFID) streamDefByteorder(fileID, CdoDefault::Byteorder);

  set_compression(fileID, p_filetype);

  m_fileID = fileID;
  m_filetype = p_filetype;

  return m_cdoStreamID;
}

int
FileStream::open_append()
{
  if (FileStream::timersEnabled()) timer_start(timer_write);

  open_lock();
  auto fileID = streamOpenAppend(m_filename.c_str());
  open_unlock();

  if (FileStream::timersEnabled()) timer_stop(timer_write);

  if (fileID < 0) cdi_open_error(fileID, "Open failed on >%s<", m_filename.c_str());

  isopen = true;

  m_filetype = streamInqFiletype(fileID);
  set_compression(fileID, m_filetype);

  m_fileID = fileID;

  return m_fileID;
}

void
FileStream::def_vlist(int p_vlistID)
{
  if (m_filetype == CDI_FILETYPE_NCZARR)
    {
      auto maxSteps = vlistNtsteps(p_vlistID);
      if (maxSteps >= 0)
        streamDefMaxSteps(m_fileID, maxSteps);
      else
        cdo_warning("Unknown number of timesteps, use operator setmaxsteps to set max. number of timesteps!");
    }

  cdo_append_history(p_vlistID, command_line());

  if (CdoDefault::DataType != CDI_UNDEFID)
    {
      auto nvars = vlistNvars(p_vlistID);
      for (int varID = 0; varID < nvars; ++varID) vlistDefVarDatatype(p_vlistID, varID, CdoDefault::DataType);

      if (CdoDefault::DataType == CDI_DATATYPE_FLT64 || CdoDefault::DataType == CDI_DATATYPE_FLT32)
        {
          for (int varID = 0; varID < nvars; ++varID) cdiDeleteKey(p_vlistID, varID, CDI_KEY_ADDOFFSET);
          for (int varID = 0; varID < nvars; ++varID) cdiDeleteKey(p_vlistID, varID, CDI_KEY_SCALEFACTOR);
        }
    }

  if (Options::nsb > 0)
    {
      auto nvars = vlistNvars(p_vlistID);
      for (int varID = 0; varID < nvars; ++varID) vlistDefVarNSB(p_vlistID, varID, Options::nsb);
    }

  if (Options::cdoChunkType != CDI_UNDEFID)
    {
      auto nvars = vlistNvars(p_vlistID);
      for (int varID = 0; varID < nvars; ++varID) cdiDefKeyInt(p_vlistID, varID, CDI_KEY_CHUNKTYPE, Options::cdoChunkType);
      auto ngrids = vlistNgrids(p_vlistID);
      for (int index = 0; index < ngrids; ++index)
        cdiDefKeyInt(vlistGrid(p_vlistID, index), CDI_GLOBAL, CDI_KEY_CHUNKTYPE, Options::cdoChunkType);
    }

  if (Options::cdoChunkSize != CDI_UNDEFID)
    {
      auto nvars = vlistNvars(p_vlistID);
      for (int varID = 0; varID < nvars; ++varID) cdiDefKeyInt(p_vlistID, varID, CDI_KEY_CHUNKSIZE, Options::cdoChunkSize);
      auto ngrids = vlistNgrids(p_vlistID);
      for (int index = 0; index < ngrids; ++index)
        cdiDefKeyInt(vlistGrid(p_vlistID, index), CDI_GLOBAL, CDI_KEY_CHUNKSIZE, Options::cdoChunkSize);
    }

  if (Options::CMOR_Mode)
    {
      cdo_def_tracking_id(p_vlistID, "tracking_id");
      cdo_def_creation_date(p_vlistID);
    }

  if (Options::VersionInfo) cdiDefAttTxt(p_vlistID, CDI_GLOBAL, "CDO", (int) strlen(cdo_comment()), cdo_comment());

#ifdef _OPENMP
  if (Threading::ompNumThreads > 1)
    cdiDefAttInt(p_vlistID, CDI_GLOBAL, "cdo_openmp_thread_number", CDI_DATATYPE_INT32, 1, &Threading::ompNumThreads);
#endif
  defDatarangeList(p_vlistID);

  if (FileStream::timersEnabled()) timer_start(timer_write);
  stream_def_vlist_locked(m_fileID, p_vlistID);
  if (FileStream::timersEnabled()) timer_stop(timer_write);
}

int
FileStream::inq_vlist()
{
  if (FileStream::timersEnabled()) timer_start(timer_read);
  auto vlistID = stream_inq_vlist_locked(m_fileID);
  if (FileStream::timersEnabled()) timer_stop(timer_read);
  if (vlistID == -1) cdo_abort("Couldn't read data from input fileID %d!", m_fileID);

  auto nsubtypes = vlistNsubtypes(vlistID);
  if (nsubtypes > 1) cdo_warning("Subtypes are unsupported, the processing results are possibly wrong!");

  if (CdoDefault::TaxisType != CDI_UNDEFID) taxisDefType(vlistInqTaxis(vlistID), CdoDefault::TaxisType);

  m_vlistID = vlistID;
  return m_vlistID;
}

void
FileStream::inq_record(int *const varID, int *const levelID)
{
  if (FileStream::timersEnabled()) timer_start(timer_read);
  stream_inq_rec_locked(m_fileID, varID, levelID);
  if (FileStream::timersEnabled()) timer_stop(timer_read);
  m_varID = *varID;
}

void
FileStream::defRecord(int varID, int levelID)
{
  if (FileStream::timersEnabled()) timer_start(timer_write);
  stream_def_rec_locked(m_fileID, varID, levelID);
  if (FileStream::timersEnabled()) timer_stop(timer_write);
  m_varID = varID;
}

void
FileStream::read_record(float *const p_data, size_t *const nmiss)
{
  if (FileStream::timersEnabled()) timer_start(timer_read);
  stream_readrecord_float_locked(m_fileID, p_data, nmiss);
  if (FileStream::timersEnabled()) timer_stop(timer_read);
}

void
FileStream::read_record(double *const p_data, size_t *const nmiss)
{
  if (FileStream::timersEnabled()) timer_start(timer_read);
  stream_readrecord_double_locked(m_fileID, p_data, nmiss);
  if (FileStream::timersEnabled()) timer_stop(timer_read);
}

void
FileStream::read_record(Field *const p_field, size_t *const nmiss)
{
  read_record(p_field->vec_d.data(), nmiss);
}

void
FileStream::write_record(float *const p_data, size_t p_nmiss)
{
  if (FileStream::timersEnabled()) timer_start(timer_write);

  auto varID = m_varID;
  if (varID < (int) m_datarangelist.size())
    if (m_datarangelist[varID].checkDatarange) m_datarangelist[varID].check_datarange(p_data, p_nmiss);

  stream_write_record_float_locked(m_fileID, p_data, p_nmiss);

  if (FileStream::timersEnabled()) timer_stop(timer_write);
}

void
FileStream::write_record(double *const p_data, size_t p_nmiss)
{
  if (FileStream::timersEnabled()) timer_start(timer_write);

  auto varID = m_varID;
  if (varID < (int) m_datarangelist.size())
    if (m_datarangelist[varID].checkDatarange) m_datarangelist[varID].check_datarange(p_data, p_nmiss);

  stream_write_record_double_locked(m_fileID, p_data, p_nmiss);

  if (FileStream::timersEnabled()) timer_stop(timer_write);
}

void
FileStream::write_record(Field *const p_field, size_t p_nmiss)
{
  write_record(p_field->vec_d.data(), p_nmiss);
}

void
FileStream::copyRecord(CdoStreamID p_destination)
{
  FileStream *fStream = dynamic_cast<FileStream *>(p_destination.get());
  stream_copy_record_locked(m_fileID, fStream->getFileID());
}
/*
 * FileStream::inq_timestep(int p_tsID)
 * stets internal state of the cdi datastructure to work on the given timestep (p_tsID) and returns the number of records that the
 * timestep contains.
 * Inquires and defines the time axis type if the timestep ID is 0 AND the taxis type is yet to be defined.
 * When the timestep inquiry was successfull m_tsID is set to the wanted p_tsID IF p_tsID != m_tsID
 * When only one process is running the timers are enabled.
 * -- last Documentation update(2019-06-14) --
 */
int
FileStream::inq_timestep(int p_tsID)
{
  if (FileStream::timersEnabled()) timer_start(timer_read);
  auto nrecs = stream_inq_time_step_locked(m_fileID, p_tsID);
  if (FileStream::timersEnabled()) timer_stop(timer_read);

  if (p_tsID == 0 && CdoDefault::TaxisType != CDI_UNDEFID) taxisDefType(vlistInqTaxis(m_vlistID), CdoDefault::TaxisType);

  if (nrecs && p_tsID != m_tsID) m_tsID = p_tsID;
  Debug(FILE_STREAM, "Current TsID: %d,  nrecs: %d", m_tsID, nrecs);

  return nrecs;
}

void
FileStream::def_timestep(int p_tsID)
{
  if (FileStream::timersEnabled()) timer_start(timer_write);
  // don't use sync -> very slow on GPFS
  //  if ( p_tsID > 0 ) streamSync(fileID);

  stream_def_time_step_locked(m_fileID, p_tsID);

  if (FileStream::timersEnabled()) timer_stop(timer_write);
}

int
FileStream::inqFileType()
{
  return streamInqFiletype(m_fileID);
}

int
FileStream::inqByteorder()
{
  return streamInqByteorder(m_fileID);
}

size_t
FileStream::getNvals()
{
  // set when the stream is closed
  // see: FileStream::close()
  return m_nvals;
}

int
FileStream::getFileID()
{
  return m_fileID;
}

void
FileStream::close()
{
  Debug(FILE_STREAM, "%s fileID %d", m_name, m_fileID);

  m_nvals = streamNvals(m_fileID);
  stream_close_locked(m_fileID);

  isopen = false;
  m_vlistID = -1;

  if (m_datarangelist.size())
    {
      m_datarangelist.clear();
      m_datarangelist.shrink_to_fit();
    }
}
