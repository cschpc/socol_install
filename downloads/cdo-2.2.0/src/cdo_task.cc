/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "cdo_task.h"

#include <cstdlib>
#include <cstdio>

#ifdef HAVE_LIBPTHREAD
static void *
cdo_task(void *task)
{
  cdo::Task *task_info = (cdo::Task *) task;

  {
    // cond_wait mutex must be locked before we can wait
    pthread_mutex_lock(&(task_info->work_mtx));

    // printf("<worker> start\n");

    // ensure boss is waiting
    pthread_mutex_lock(&(task_info->boss_mtx));

    // signal to boss that setup is complete
    task_info->state = cdo::TaskState::IDLE;

    // wake-up signal
    pthread_cond_signal(&(task_info->boss_cond));
    pthread_mutex_unlock(&(task_info->boss_mtx));

    while (1)
      {
        pthread_cond_wait(&(task_info->work_cond), &(task_info->work_mtx));

        if (cdo::TaskState::DIE == task_info->state) break;  // kill thread

        if (cdo::TaskState::IDLE == task_info->state) continue;  // accidental wake-up

        // do blocking task
        // printf("<worker> JOB start\n");
        task_info->result = task_info->routine(task_info->arg);
        // printf("<worker> JOB end\n");

        // ensure boss is waiting
        pthread_mutex_lock(&(task_info->boss_mtx));

        // indicate that job is done
        task_info->state = cdo::TaskState::IDLE;

        // wake-up signal
        pthread_cond_signal(&(task_info->boss_cond));
        pthread_mutex_unlock(&(task_info->boss_mtx));
      }

    pthread_mutex_unlock(&(task_info->work_mtx));
  }

  pthread_exit(nullptr);

  return nullptr;
}
#endif

namespace cdo
{
void
Task::start(void *(*task_routine)(void *), void *task_arg)
{
  // ensure worker is waiting
#ifdef HAVE_LIBPTHREAD
  pthread_mutex_lock(&(this->work_mtx));
#endif

  // set job information & state
  this->routine = task_routine;
  this->arg = task_arg;
  this->state = cdo::TaskState::JOB;

#ifndef HAVE_LIBPTHREAD
  this->result = this->routine(this->arg);
#endif

  // wake-up signal
#ifdef HAVE_LIBPTHREAD
  pthread_cond_signal(&(this->work_cond));
  pthread_mutex_unlock(&(this->work_mtx));
#endif
}

void *
Task::wait()
{
#ifdef HAVE_LIBPTHREAD
  while (1)
    {
      if (cdo::TaskState::IDLE == this->state) break;

      pthread_cond_wait(&(this->boss_cond), &(this->boss_mtx));

      // if (cdo::TaskState::IDLE == task_info->state) break;
    }
#endif

  return this->result;
}

Task::Task()
{
#ifdef HAVE_LIBPTHREAD
  pthread_attr_t attr;
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
  size_t stacksize;
  pthread_attr_getstacksize(&attr, &stacksize);
  if (stacksize < 2097152)
    {
      stacksize = 2097152;
      pthread_attr_setstacksize(&attr, stacksize);
    }

  pthread_cond_init(&(this->work_cond), nullptr);
  pthread_mutex_init(&(this->work_mtx), nullptr);
  pthread_cond_init(&(this->boss_cond), nullptr);
  pthread_mutex_init(&(this->boss_mtx), nullptr);

  pthread_mutex_lock(&(this->boss_mtx));

  pthread_create(&(this->thread), &attr, cdo_task, (void *) this);

  this->wait();
#endif
}

Task::~Task()
{
#ifdef HAVE_LIBPTHREAD
  // ensure the worker is waiting
  pthread_mutex_lock(&(this->work_mtx));

  // printf("Task::delete: send DIE to <worker>\n");
  this->state = cdo::TaskState::DIE;

  // wake-up signal
  pthread_cond_signal(&(this->work_cond));
  pthread_mutex_unlock(&(this->work_mtx));

  // wait for thread to exit
  pthread_join(this->thread, nullptr);

  pthread_mutex_destroy(&(this->work_mtx));
  pthread_cond_destroy(&(this->work_cond));

  pthread_mutex_unlock(&(this->boss_mtx));
  pthread_mutex_destroy(&(this->boss_mtx));
  pthread_cond_destroy(&(this->boss_cond));
#endif
}
}  // namespace cdo

#ifdef TEST_CDO_TASK
// g++ -g -DTEST_CDO_TASK -DHAVE_LIBPTHREAD cdo_task.cc

void *
mytask(void *arg)
{
  printf("run mytask\n");
  return nullptr;
}

int
main(int argc, char **argv)
{
  cdo::Task task;

  void *myarg = nullptr;
  void *myresult;

  task.start(mytask, myarg);

  myresult = task.wait();

  task.start(mytask, myarg);

  myresult = task.wait();

  return 0;
}
#endif
