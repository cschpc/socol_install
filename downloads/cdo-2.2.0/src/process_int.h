/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida
          Oliver Heidmann

*/
#ifndef PROCESS_INT_H
#define PROCESS_INT_H

#include "processManager.h"
#include "cdo_cdi_wrapper.h"
#include "process.h"

Process &process_self(void);
bool cdo_assert_files_only();
bool cdo_stream_is_pipe(CdoStreamID streamID);
bool data_is_unchanged();
const std::string &cdo_operator_argv(size_t p_idx);
const std::vector<std::string> &cdo_get_oper_argv();
std::string cdo_get_obase();
const char *cdo_get_stream_name(int p_streamIndex);
const char *cdo_operator_enter(int operID);
const char *cdo_operator_name(int operID);
std::string cdo_module_name();
const char *process_inq_prompt(void);
std::string cdo_get_command_from_in_stream(int p_streamIndex);
int cdo_filetype(void);
int cdo_inq_byteorder(CdoStreamID streamID);
int cdo_inq_filetype(CdoStreamID streamID);
int cdo_stream_cnt(void);
int cdo_stream_number();
int cdo_operator_argc(void);
void operator_check_argc(int numargs);
void operator_input_arg(const char *enter);
void process_def_var_num(int nvars);
void cdo_inq_grib_info(CdoStreamID streamID, int *intnum, float *fltnum, off_t *bignum);

void cdo_set_nan(double missval, size_t gridsize, double *array);
// ***********************************************************

void cdo_inq_record(CdoStreamID streamID, int *varID, int *levelID);
void cdo_def_record(CdoStreamID streamID, int varID, int levelID);

void cdo_read_record_f(CdoStreamID streamID, float *data, size_t *nmiss);
void cdo_read_record(CdoStreamID streamID, double *data, size_t *nmiss);
void cdo_read_record(CdoStreamID streamID, Field &field);
void cdo_read_record(CdoStreamID streamID, Field3D &field, int levelID, size_t *nmiss);

void cdo_write_record_f(CdoStreamID streamID, float *data, size_t nmiss);
void cdo_write_record(CdoStreamID streamID, double *data, size_t nmiss);
void cdo_write_record(CdoStreamID streamID, Field &data);
void cdo_write_record(CdoStreamID streamID, Field3D &data, int levelID, size_t nmiss);

void cdo_copy_record(CdoStreamID streamIDdest, CdoStreamID streamIDsrc);

void cdo_add_steps(int numSteps);
int cdo_stream_inq_timestep(CdoStreamID streamID, int tsID);
void cdo_def_timestep(CdoStreamID streamID, int tsID);
int cdo_stream_inq_vlist(CdoStreamID streamID);
void cdo_def_vlist(CdoStreamID streamID, int vlistID);

void cdo_stream_close(CdoStreamID streamID);

// ***********************************************************
void cdo_def_comp_type(CdoStreamID p_streamID, int p_cdi_compression_type);
// ***********************************************************
CdoStreamID cdo_open_append(int outStreamIDX);
CdoStreamID cdo_open_read(int inStreamIDX);
CdoStreamID cdo_open_write(int outStreamIDX, int filetype = CDI_UNDEFID);
CdoStreamID cdo_open_write(const std::string &p_filename, int filetype = CDI_UNDEFID);
int cdo_operator_add(const char *name, int f1, int f2, const char *enter);
int cdo_operator_f1(int operID);
int cdo_operator_f2(int operID);
int cdo_operator_id(void);
void cdo_finish();
void cdo_initialize(void *process);

#endif
