#ifndef CDO_WTIME_H
#define CDO_WTIME_H

#ifdef _OPENMP
//#define CDO_USE_OMP_WTIME
#endif

#ifdef CDO_USE_OMP_WTIME

#ifdef _OPENMP
#include <omp.h>  // omp_get_wtime
#endif

inline double
cdo_get_wtime()
{
#ifdef _OPENMP
  return omp_get_wtime();
#else
  return 0;
#endif
}

#else

#include <chrono>

inline double
cdo_get_wtime()
{
  using namespace std::chrono;
  return duration_cast<duration<double>>(steady_clock::now().time_since_epoch()).count();
}

#endif

#endif
