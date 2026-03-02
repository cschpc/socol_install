#ifndef YAC_LD_INTERFACE_H
#define YAC_LD_INTERFACE_H

#ifdef YAC_FOR_CDO
#define YAC_LONG_DOUBLE_INTERFACE_ID 1
#endif

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifndef YAC_LONG_DOUBLE_INTERFACE_ID
#error None of the supported long double arithmetic interfaces is available
#endif

#if YAC_LONG_DOUBLE_INTERFACE_ID == 1 // long double C type

#elif YAC_LONG_DOUBLE_INTERFACE_ID == 2 // GNU MPFR

#include <gmp.h>
#include <mpfr.h>
#include <mpf2mpfr.h>

#undef USE_MPF
#define USE_MPF

#elif YAC_LONG_DOUBLE_INTERFACE_ID == 3 // GNU MP

#include <gmp.h>

#undef USE_MPF
#define USE_MPF

#else

#error Unexpected value for YAC_LONG_DOUBLE_INTERFACE_ID

#endif

#ifdef USE_MPF

#ifndef YAC_MPF_PRECISION
// This is what we usually get for the long double C type on systems that
// support it:
#define YAC_MPF_PRECISION 64
#endif

static void yac_mpf_init_vector(mpf_t x[3]) {
  mpf_init2(x[0], YAC_MPF_PRECISION);
  mpf_init2(x[1], YAC_MPF_PRECISION);
  mpf_init2(x[2], YAC_MPF_PRECISION);
}

static void yac_mpf_set_vector_d(mpf_t rop[3], const double op[3]) {
  yac_mpf_init_vector(rop);
  mpf_set_d(rop[0], op[0]);
  mpf_set_d(rop[1], op[1]);
  mpf_set_d(rop[2], op[2]);
}

static void yac_mpf_clear_vector(mpf_t x[3]) {
  mpf_clear(x[0]);
  mpf_clear(x[1]);
  mpf_clear(x[2]);
}

#endif // USE_MPF

#endif // YAC_LD_INTERFACE_H

