/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_LIBPTHREAD
#include <pthread.h>
#endif

#include "process_int.h"
#include "par_io.h"

void *
read_record(void *arg)
{
  auto *read_arg = (read_arg_t *) arg;
  auto streamID = read_arg->streamID;
  auto p_varID = read_arg->varID;
  auto p_levelID = read_arg->levelID;
  auto p_nmiss = read_arg->nmiss;
  auto p_array = read_arg->array;

  // fprintf(stderr, "read_record: streamID = %d\n", streamID);
  cdo_inq_record(streamID, p_varID, p_levelID);
  cdo_read_record(streamID, p_array, p_nmiss);
  // fprintf(stderr, "read_record: varID %d levelID %d\n", *p_varID, *p_levelID);

  return nullptr;
}

void
par_read_record(CdoStreamID streamID, int *varID, int *levelID, double *array, size_t *nmiss, par_io_t *parIO)
{
  bool lpario = false;
  int recID = 0, nrecs = 0;
#ifdef HAVE_LIBPTHREAD
  pthread_t thrID = 0;
  // pthread_attr_t attr;
  int rval;

  if (parIO)
    {
      lpario = true;
      recID = parIO->recID;
      nrecs = parIO->nrecs;
      thrID = parIO->thrID;
    }
#endif

  if (recID == 0 || !lpario)
    {
      read_arg_t read_arg;
      read_arg.streamID = streamID;
      read_arg.varID = varID;
      read_arg.levelID = levelID;
      read_arg.nmiss = nmiss;
      read_arg.array = array;

      read_record(&read_arg);
    }
#ifdef HAVE_LIBPTHREAD
  else
    {
      // fprintf(stderr, "parIO1: %ld streamID %d %d %d\n", (long)thrID, streamID, recID, nrecs);
      rval = pthread_join(thrID, nullptr);
      if (rval != 0) cdo_abort("pthread_join failed!");

      *varID = parIO->varID;
      *levelID = parIO->levelID;
      *nmiss = parIO->nmiss;
      // fprintf(stderr, "parIO2: %ld streamID %d %d %d\n", (long)thrID, streamID, *varID, *levelID);
      array_copy(parIO->array_size, parIO->array, array);
    }

  if (lpario && nrecs > 1)
    {
      read_arg_t *read_arg = &(parIO->read_arg);
      if ((recID + 1) < nrecs)
        {
          if (recID == 0)
            {
              pthread_attr_init(&parIO->attr);
              pthread_attr_setdetachstate(&parIO->attr, PTHREAD_CREATE_JOINABLE);
            }

          read_arg->streamID = streamID;
          read_arg->varID = &parIO->varID;
          read_arg->levelID = &parIO->levelID;
          read_arg->nmiss = &parIO->nmiss;
          read_arg->array = parIO->array;

          // fprintf(stderr, "pthread_create: streamID %d %d\n", read_arg->streamID,streamID);
          rval = pthread_create(&thrID, &parIO->attr, read_record, read_arg);
          if (rval != 0) cdo_abort("pthread_create failed!");

          // fprintf(stderr, "thrID = %ld\n", (long) thrID);
          parIO->thrID = thrID;
        }
      else
        pthread_attr_destroy(&parIO->attr);
    }
#endif
}
