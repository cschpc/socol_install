/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#ifndef CDO_TASK_H
#define CDO_TASK_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_LIBPTHREAD
#include <pthread.h>
#endif

namespace cdo
{
enum class TaskState
{
  SETUP,
  IDLE,
  JOB,
  DIE
};

class Task
{
public:
  void *(*routine)(void *) = nullptr;
  void *arg = nullptr;
  void *result = nullptr;
  enum TaskState state = cdo::TaskState::SETUP;
#ifdef HAVE_LIBPTHREAD
  pthread_t thread;
  pthread_cond_t work_cond;
  pthread_mutex_t work_mtx;
  pthread_cond_t boss_cond;
  pthread_mutex_t boss_mtx;
#endif

  Task();
  ~Task();
  void start(void *(*task_routine)(void *), void *task_arg);
  void *wait();
};
}  // namespace cdo

#endif /* CDO_TASK_H */
