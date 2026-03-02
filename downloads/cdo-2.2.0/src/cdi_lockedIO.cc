/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include "cdi_lockedIO.h"
#include "cdo_options.h"
#include "cdo_output.h"
#include "cthread_debug.h"

#include <mutex>

static std::mutex streamOpenReadMutex;
static std::mutex streamMutex;

int
stream_open_read_locked(const char *const p_filename)
{
  open_lock();
  const auto streamID = streamOpenRead(p_filename);
  open_unlock();
  if (streamID < 0) cdi_open_error(streamID, "Open failed on >%s<", p_filename);

  return streamID;
}

void
stream_close_locked(const int p_fileID)
{
  if (Threading::cdoLockIO) cthread_mutex_lock(streamMutex);
  streamClose(p_fileID);
  if (Threading::cdoLockIO) cthread_mutex_unlock(streamMutex);
}

void
stream_inq_rec_locked(const int p_fileID, int *const p_varID, int *const p_levelID)
{
  if (Threading::cdoLockIO) cthread_mutex_lock(streamMutex);
  streamInqRecord(p_fileID, p_varID, p_levelID);
  if (Threading::cdoLockIO) cthread_mutex_unlock(streamMutex);
}

void
stream_def_rec_locked(const int p_fileID, const int p_varID, const int p_levelID)
{
  if (Threading::cdoLockIO) cthread_mutex_lock(streamMutex);
  streamDefRecord(p_fileID, p_varID, p_levelID);
  if (Threading::cdoLockIO) cthread_mutex_unlock(streamMutex);
}

void
stream_readrecord_float_locked(const int p_fileID, float *const p_data, size_t *const p_nmiss)
{
  if (Threading::cdoLockIO) cthread_mutex_lock(streamMutex);
  streamReadRecordF(p_fileID, p_data, p_nmiss);
  if (Threading::cdoLockIO) cthread_mutex_unlock(streamMutex);
}

void
stream_readrecord_double_locked(const int p_fileID, double *const p_data, size_t *const p_nmiss)
{
  if (Threading::cdoLockIO) cthread_mutex_lock(streamMutex);
  streamReadRecord(p_fileID, p_data, p_nmiss);
  if (Threading::cdoLockIO) cthread_mutex_unlock(streamMutex);
}

void
stream_def_vlist_locked(const int p_fileID, const int p_vlistID)
{
  if (Threading::cdoLockIO) cthread_mutex_lock(streamMutex);
  streamDefVlist(p_fileID, p_vlistID);
  if (Threading::cdoLockIO) cthread_mutex_unlock(streamMutex);
}

int
stream_inq_vlist_locked(const int p_fileID)
{
  if (Threading::cdoLockIO) cthread_mutex_lock(streamMutex);
  const auto vlistID = streamInqVlist(p_fileID);
  if (Threading::cdoLockIO) cthread_mutex_unlock(streamMutex);

  return vlistID;
}

void
stream_write_record_double_locked(const int p_fileID, const double *const p_data, const size_t p_nmiss)
{
  if (Threading::cdoLockIO) cthread_mutex_lock(streamMutex);
  streamWriteRecord(p_fileID, p_data, p_nmiss);
  if (Threading::cdoLockIO) cthread_mutex_unlock(streamMutex);
}

void
stream_write_record_float_locked(const int p_fileID, const float *const p_data, const size_t p_nmiss)
{
  if (Threading::cdoLockIO) cthread_mutex_lock(streamMutex);
  streamWriteRecordF(p_fileID, p_data, p_nmiss);
  if (Threading::cdoLockIO) cthread_mutex_unlock(streamMutex);
}

int
stream_inq_time_step_locked(const int p_fileID, const int p_tsID)
{
  if (Threading::cdoLockIO) cthread_mutex_lock(streamMutex);
  const auto nrecs = streamInqTimestep(p_fileID, p_tsID);
  if (Threading::cdoLockIO) cthread_mutex_unlock(streamMutex);

  return nrecs;
}

int
stream_def_time_step_locked(const int p_fileID, const int p_tsID)
{
  if (Threading::cdoLockIO) cthread_mutex_lock(streamMutex);
  const auto success = streamDefTimestep(p_fileID, p_tsID);
  if (Threading::cdoLockIO) cthread_mutex_unlock(streamMutex);
  return success;
}

int
stream_copy_record_locked(const int p_fileID, const int p_targetFileID)
{
  if (Threading::cdoLockIO) cthread_mutex_lock(streamMutex);
  streamCopyRecord(p_fileID, p_targetFileID);
  if (Threading::cdoLockIO) cthread_mutex_unlock(streamMutex);
  return p_targetFileID;
}

void
vlist_copy_flag_locked(const int p_vlistID2, const int p_vlistID1)
{
  cthread_mutex_lock(streamMutex);
  vlistCopyFlag(p_vlistID2, p_vlistID1);
  cthread_mutex_unlock(streamMutex);
}

void
open_lock(void)
{
  cthread_mutex_lock(Threading::cdoLockIO ? streamMutex : streamOpenReadMutex);
}

void
open_unlock(void)
{
  cthread_mutex_unlock(Threading::cdoLockIO ? streamMutex : streamOpenReadMutex);
}

void
cdo_vlist_copy_flag(const int vlistID2, const int vlistID1)
{
  vlist_copy_flag_locked(vlistID2, vlistID1);
}
