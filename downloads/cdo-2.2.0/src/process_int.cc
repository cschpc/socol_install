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

#include <cassert>

// Debug and message includes
#include <cdi.h>

#include "cdo_options.h"
#include "cdo_default_values.h"
#include "process.h"

thread_local Process *localProcess;
static int NumProcessActive = 0;

int
cdo_filetype(void)
{
  if (CdoDefault::FileType == CDI_UNDEFID)
    {
      CdoDefault::FileType = CDI_FILETYPE_GRB;
      if (Options::cdoVerbose) cdo_print("Set default filetype to GRIB1");
    }

  return CdoDefault::FileType;
}

Process &
process_self(void)
{
  return *localProcess;
}

void
process_def_var_num(int nvars)
{
  localProcess->nvars += nvars;
}

/* cdoStreamInqTimeStep(...)
 * sets the given p_pstreamptr to given tsID and returns the number of records this timestep contains
 */
int
cdo_stream_inq_timestep(CdoStreamID p_pstreamptr, int tsID)
{
  Debug(PROCESS_INT, "%s pstreamID %d", p_pstreamptr->m_name, p_pstreamptr->get_id());
  const auto nrecs = p_pstreamptr->inq_timestep(tsID);
  if (nrecs && tsID == p_pstreamptr->getTsID())
    {
      localProcess->ntimesteps++;
      Debug(PROCESS_INT, "Timestep cnt for localProcess: %s: %d", localProcess->prompt, localProcess->ntimesteps);
    }
  return nrecs;
}

void
cdo_add_steps(int numSteps)
{
  localProcess->ntimesteps += numSteps;
}

int
cdo_operator_argc(void)
{
  return localProcess->get_oper_argc();
}

const std::string &
cdo_operator_argv(size_t idx)
{
  const auto argc = localProcess->get_oper_argc();
  assert(((void) "Internal error: argument index out of bounds!", (idx < argc)));

  return localProcess->m_oargv[idx];
}

const std::vector<std::string> &
cdo_get_oper_argv()
{
  return localProcess->m_oargv;
}

void
operator_check_argc(int numargs)
{
  const int argc = localProcess->get_oper_argc();

  if (argc < numargs)
    cdo_abort("Too few arguments! Need %d found %d.", numargs, argc);
  else if (argc > numargs)
    cdo_abort("Too many arguments! Need %d found %d.", numargs, argc);
}

static void
print_enter(const char *prompt, const char *enter)
{
  set_text_color(stderr, BRIGHT, MAGENTA);
  fprintf(stderr, "%-16s : ", prompt);
  reset_text_color(stderr);
  fprintf(stderr, "Enter %s > ", enter);
}

static std::string
get_input()
{
  std::string line;
  std::string fline = "";

  do {
      std::getline(std::cin, line);
      fline += line;
    }
  while (std::string(line).find("\\") != std::string::npos);

  return fline;
}

void
operator_input_arg(const char *enter)
{
  while (true)
    {
      if (localProcess->get_oper_argc() != 0) return;

      if (enter) print_enter(localProcess->prompt, enter);
      std::stringstream stringStream(get_input());

      std::string line;
      while (std::getline(stringStream, line))
        {
          std::size_t prev = 0, pos;
          while ((pos = line.find_first_of(", \\", prev)) != std::string::npos)
            {
              if (pos > prev) localProcess->m_oargv.push_back(line.substr(prev, pos - prev));
              prev = pos + 1;
            }
          if (prev < line.length()) localProcess->m_oargv.push_back(line.substr(prev, std::string::npos));
        }
    }
}

int
cdo_operator_add(const char *name, int f1, int f2, const char *enter)
{
  return localProcess->operator_add(name, f1, f2, enter);
}

int
cdo_operator_id(void)
{
  return localProcess->get_operator_id();
}

int
cdo_operator_f1(int operID)
{
  return localProcess->oper[operID].f1;
}

int
cdo_operator_f2(int operID)
{
  return localProcess->oper[operID].f2;
}

const char *
cdo_operator_name(int operID)
{
  return localProcess->oper[operID].name.c_str();
}

std::string
cdo_module_name()
{
  return get_module_name_to(cdo_operator_name(cdo_operator_id()));
}

const char *
cdo_operator_enter(int operID)
{
  return localProcess->oper[operID].enter;
}

int
cdo_stream_number()
{
  return operator_stream_number(localProcess->operatorName);
}

int
cdo_stream_cnt(void)
{
  return localProcess->m_streamCnt;
}

CdoStreamID
cdo_open_read(int inStreamIDX)
{
  Debug(PROCESS_INT, "Getting in stream %d of process %d", inStreamIDX, localProcess->m_ID);

  if (localProcess->get_stream_cnt_in() < inStreamIDX || inStreamIDX < 0)
    cdo_abort("instream %d of process %d not found", inStreamIDX, localProcess->m_ID);

  const auto inStream = localProcess->inputStreams[inStreamIDX];
  inStream->open_read();

  return inStream;
}

/*parameters:
 *  p_outStreamIDX: In operator defined out stream ID
 *  filetype      : Currently not used in operators! Default value is CDI_UNDEFID.
 */
CdoStreamID
cdo_open_write(int p_outStreamIDX, int filetype)
{
  if (filetype == CDI_UNDEFID) filetype = cdo_filetype();
  Debug(PROCESS_INT, "Getting out stream %d of process %d", p_outStreamIDX, localProcess->m_ID);

  const int outStreamIDX = p_outStreamIDX - localProcess->inputStreams.size();
  if (outStreamIDX > localProcess->get_stream_cnt_out() || outStreamIDX < 0)
    {
      cdo_abort("outstream %d of %d not found. Called with streamIdx = %d", outStreamIDX, localProcess->m_ID, p_outStreamIDX);
    }

  const auto outStream = localProcess->outputStreams[outStreamIDX];
  outStream->open_write(filetype);

  return outStream;
}

CdoStreamID
cdo_open_write(const std::string &p_filename, int filetype)
{
  if (filetype == CDI_UNDEFID) filetype = cdo_filetype();

  localProcess->add_file_out_stream(p_filename);

  const auto pstreamID = localProcess->outputStreams.back()->open_write(filetype);
  if (pstreamID == -1) cdo_abort("Could not create pstream for file: %s", p_filename);

  return localProcess->outputStreams.back();
}

CdoStreamID
cdo_open_append(int p_outFileIndex)
{
  const auto streamIndex = p_outFileIndex - localProcess->inputStreams.size();
  const auto outStream = localProcess->outputStreams[streamIndex];

  const auto fileID = outStream->open_append();
  if (fileID < 0) cdi_open_error(fileID, "Open failed on >%s<", outStream->m_name.c_str());

  return outStream;
}

/**
 * Returns the output stream name as std::string for \p outStreamID.
 */
static const char *
cdoGetOutStreamName()
{
  return localProcess->get_out_stream_name();
}

/**
 * Returns the input stream name as std::string for inStreamID.
 */

const std::shared_ptr<Process>
cdoGetInputChild(int p_inID)
{
  for (const std::shared_ptr<Process> &processPtr : localProcess->childProcesses)
    {
      if (processPtr->has_out_stream(localProcess->inputStreams[p_inID])) return processPtr;
    }
  return nullptr;
}

static const char *
cdoGetInStreamName(int p_inStream)
{
  return localProcess->inputStreams[p_inStream]->m_name.c_str();
}

const char *
cdo_get_stream_name(int p_streamIndex)
{
  const char *streamName = nullptr;
  Debug(PROCESS_INT, "stridx %d", p_streamIndex);
  if (p_streamIndex >= static_cast<int>(localProcess->inputStreams.size()))
    {
      Debug(PROCESS_INT, "Getting output stream name %d", p_streamIndex);
      streamName = cdoGetOutStreamName();
    }
  else
    {
      Debug(PROCESS_INT, "Getting input stream name %d", p_streamIndex);
      streamName = cdoGetInStreamName(p_streamIndex);
    }

  Debug(PROCESS_INT, "StreamName is: %s", streamName);

  return streamName;
}

std::string
cdo_get_command_from_in_stream(int p_streamIndex)
{
  const std::shared_ptr<Process> process = cdoGetInputChild(p_streamIndex);
  if (!process)  // then this should be a file and cdo_get_stream_name can handle this
    {
      return cdo_get_stream_name(p_streamIndex);
    }
  return process->operatorName;
}

bool
cdo_assert_files_only()
{
  return localProcess->has_no_pipes();
}

std::string
cdo_get_obase()
{
  return localProcess->get_obase();
}

void
cdo_initialize(void *p_process)
{
#if defined(_OPENMP)
  omp_set_num_threads(Threading::ompNumThreads);  // Has to be called for every module (pthread)!
#endif
  localProcess = (Process *) p_process;
#ifdef HAVE_LIBPTHREAD
  localProcess->threadID = pthread_self();
#endif
  Debug(PROCESS_INT, "Initializing process: %s (id: %d)", localProcess->operatorName, localProcess->get_id());

#ifdef HAVE_LIBPTHREAD
  Debug(PROCESS_INT, "process %d thread %ld", localProcess->m_ID, pthread_self());
#endif
}

const char *
process_inq_prompt(void)
{
  return localProcess ? localProcess->inq_prompt() : cdo::progname;
}

extern "C" size_t getPeakRSS();

void
cdo_finish(void)
{
  Debug(PROCESS_INT, "Finishing process: %d", localProcess->m_ID);

#ifdef HAVE_LIBPTHREAD
  Debug(PROCESS_INT, "process %d thread %ld", localProcess->m_ID, pthread_self());
#endif

  if (!Options::silentMode && (cdo::stdoutIsTerminal || Options::cdoVerbose)) localProcess->print_processed_values();

  // localProcess->close_streams();
  localProcess->set_inactive();
  NumProcessActive--;
}

/* TODO move cdoClose internals into Pstream::close/ CdoStream::close */
void
cdo_stream_close(CdoStreamID pstreamPtr)
{
  Debug(PROCESS_INT, "Adding %d to pstream %d %s", localProcess->inq_nvals(), pstreamPtr->get_id(), pstreamPtr->m_name);
  pstreamPtr->close();
}

void
cdo_close_streams()
{
  localProcess->close_streams();
}

/******************************************************************/
/* Functions that only work on streamID's and do not need process */
/******************************************************************/

int
cdo_stream_inq_vlist(CdoStreamID p_pstreamPtr)
{
  if (p_pstreamPtr == nullptr) return -1;

  Debug(PROCESS_INT, "Inquiring Vlist from pstream %d", p_pstreamPtr->get_id());
  const auto vlistID = p_pstreamPtr->inq_vlist();
  if (vlistNumber(vlistID) == CDI_COMP && cdo_stream_number() == CDI_REAL)
    cdo_abort("Fields with complex numbers are not supported by this operator!");

  if (vlistNumber(vlistID) == CDI_REAL && cdo_stream_number() == CDI_COMP)
    cdo_abort("This operator needs fields with complex numbers!");

  process_def_var_num(vlistNvars(vlistID));

  return vlistID;
}

// - - - - - - -
void
cdo_write_record_f(CdoStreamID p_pstreamPtr, float *data, size_t nmiss)
{
  p_pstreamPtr->write_record(data, nmiss);
}

void
cdo_write_record(CdoStreamID p_pstreamPtr, double *data, size_t nmiss)
{
  p_pstreamPtr->write_record(data, nmiss);
}

void
cdo_write_record(CdoStreamID p_pstreamPtr, Field &field)
{
  if (field.memType == MemType::Float)
    cdo_write_record_f(p_pstreamPtr, field.vec_f.data(), field.nmiss);
  else
    cdo_write_record(p_pstreamPtr, field.vec_d.data(), field.nmiss);
}

void
cdo_write_record(CdoStreamID p_pstreamPtr, Field3D &field, int levelID, size_t nmiss)
{
  const auto offset = levelID * field.gridsize * field.nwpv;
  if (field.memType == MemType::Float)
    cdo_write_record_f(p_pstreamPtr, field.vec_f.data() + offset, nmiss);
  else
    cdo_write_record(p_pstreamPtr, field.vec_d.data() + offset, nmiss);
}
// - - - - - - -
void
cdo_def_vlist(CdoStreamID p_pstreamPtr, int vlistID)
{
  p_pstreamPtr->def_vlist(vlistID);
}
// - - - - - - -
void
cdo_def_timestep(CdoStreamID p_pstreamPtr, int tsID)
{
  p_pstreamPtr->def_timestep(tsID);
}
// - - - - - - -
//
int
cdo_inq_filetype(CdoStreamID p_pstreamPtr)
{
  return p_pstreamPtr->inqFileType();
}
// - - - - - - -
void
cdo_inq_grib_info(CdoStreamID p_streamPtr, int *intnum, float *fltnum, off_t *bignum)
{
  streamInqGRIBinfo(p_streamPtr->m_fileID, intnum, fltnum, bignum);
}

// - - - - - - -
int
cdo_inq_byteorder(CdoStreamID p_pstreamPtr)
{
  return p_pstreamPtr->inqByteorder();
}

// - - - - - - -
void
cdo_read_record_f(CdoStreamID p_pstreamPtr, float *data, size_t *nmiss)
{
  if (data == nullptr) cdo_abort("Data pointer not allocated (cdo_read_record)!");
  p_pstreamPtr->read_record(data, nmiss);
}

void
cdo_read_record(CdoStreamID p_pstreamPtr, double *data, size_t *nmiss)
{
  if (data == nullptr) cdo_abort("Data pointer not allocated (cdo_read_record)!");
  p_pstreamPtr->read_record(data, nmiss);
}

void
cdo_read_record(CdoStreamID streamID, Field &field)
{
  if (field.memType == MemType::Float)
    cdo_read_record_f(streamID, field.vec_f.data(), &field.nmiss);
  else
    cdo_read_record(streamID, field.vec_d.data(), &field.nmiss);
}

void
cdo_read_record(CdoStreamID streamID, Field3D &field, int levelID, size_t *nmiss)
{
  const auto offset = levelID * field.gridsize * field.nwpv;
  if (field.memType == MemType::Float)
    cdo_read_record_f(streamID, field.vec_f.data() + offset, nmiss);
  else
    cdo_read_record(streamID, field.vec_d.data() + offset, nmiss);
}
// - - - - - - -

void
cdo_copy_record(CdoStreamID pstreamPtrDest, CdoStreamID pstreamPtrSrc)
{
  Debug(PROCESS_INT, "pstreamIDdest = %d pstreamIDsrc = %d", pstreamPtrDest->get_id(), pstreamPtrSrc->get_id());
  pstreamPtrDest->copyRecord(pstreamPtrSrc);
}

// - - - - - - -

void
cdo_inq_record(CdoStreamID pstreamptr, int *varID, int *levelID)
{
  pstreamptr->inq_record(varID, levelID);
}

void
cdo_def_record(CdoStreamID pstreamptr, int varID, int levelID)
{
  pstreamptr->defRecord(varID, levelID);
}
// - - - - - - -

void
cdo_def_comp_type(CdoStreamID p_streamID, int p_cdi_compression_type)
{
  streamDefCompType(p_streamID->m_fileID, p_cdi_compression_type);
}

bool
data_is_unchanged()
{
  const auto unchanged
      = (localProcess->m_ID == 0 && localProcess->has_no_pipes() && Options::cdoRegulargrid == false && CdoDefault::FileType == -1
         && CdoDefault::DataType == -1 && CdoDefault::Byteorder == -1 && Options::CDO_Memtype == MemType::Native);
  return unchanged;
}

void
cdo_set_nan(double missval, size_t gridsize, double *array)
{
  if (std::isnan(missval))
    {
      constexpr double newmissval = -9.e33;
      for (size_t i = 0; i < gridsize; ++i)
        if (DBL_IS_EQUAL(array[i], missval)) array[i] = newmissval;
    }
}
