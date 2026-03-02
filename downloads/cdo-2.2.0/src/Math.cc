/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Math       abs             Absolute value
      Math       sqr             Square
      Math       sqrt            Square root
      Math       exp             Exponential
      Math       ln              Natural logarithm
      Math       log10           Base 10 logarithm
      Math       sin             Sine
      Math       cos             Cosine
      Math       tan             Tangent
      Math       asin            Arc sine
      Math       acos            Arc cosine
      Math       atan            Arc tangent
      Math       pow             Power
      Math       reci            Reciprocal
*/

#include <cstdlib>
#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "param_conversion.h"

static void
check_out_of_range(size_t &nmiss, const size_t len, double missval, Varray<double> &v, double rmin, double rmax)
{
  if (nmiss)
    {
      for (size_t i = 0; i < len; ++i)
        if (!DBL_IS_EQUAL(v[i], missval) && (v[i] < rmin || v[i] > rmax))
          {
            v[i] = missval;
            nmiss++;
          }
    }
  else
    {
      for (size_t i = 0; i < len; ++i)
        if (v[i] < rmin || v[i] > rmax)
          {
            v[i] = missval;
            nmiss++;
          }
    }
}

static void
check_lower_range(size_t &nmiss, const size_t len, double missval, Varray<double> &v, double rmin)
{
  if (nmiss)
    {
      for (size_t i = 0; i < len; ++i)
        if (!DBL_IS_EQUAL(v[i], missval) && v[i] < rmin)
          {
            v[i] = missval;
            nmiss++;
          }
    }
  else
    {
      for (size_t i = 0; i < len; ++i)
        if (v[i] < rmin)
          {
            v[i] = missval;
            nmiss++;
          }
    }
}

template <typename T, class UnaryOperation>
void
math_varray_transform(const size_t nmiss, const size_t len, double missval1, const Varray<T> &array1, Varray<T> &array2,
                      UnaryOperation unary_op)
{
  if (nmiss)
    for (size_t i = 0; i < len; ++i) array2[i] = DBL_IS_EQUAL(array1[i], missval1) ? missval1 : unary_op(array1[i]);
  else
    for (size_t i = 0; i < len; ++i) array2[i] = unary_op(array1[i]);
}

static void
math_varray_sqr_cplx(const size_t len, const Varray<double> &array1, Varray<double> &array2)
{
  for (size_t i = 0; i < len; ++i)
    {
      array2[i * 2] = array1[i * 2] * array1[i * 2] + array1[i * 2 + 1] * array1[i * 2 + 1];
      array2[i * 2 + 1] = 0.0;
    }
}

static void
math_varray_sqrt_cplx(const size_t len, double missval1, const Varray<double> &array1, Varray<double> &array2)
{
  auto missval2 = missval1;
  auto rsqrt2 = 1.0 / std::sqrt(2.0);
  for (size_t i = 0; i < len; ++i)
    {
      double abs = SQRTMN(ADDMN(MULMN(array1[2 * i], array1[2 * i]), MULMN(array1[2 * i + 1], array1[2 * i + 1])));
      array2[i * 2] = MULMN(rsqrt2, SQRTMN(ADDMN(array1[i * 2], abs)));
      array2[i * 2 + 1] = MULMN(rsqrt2, DIVMN(array1[2 * i + 1], SQRTMN(ADDMN(array1[2 * i], abs))));
    }
}

static void
math_varray_conj_cplx(const size_t len, const Varray<double> &array1, Varray<double> &array2)
{
  for (size_t i = 0; i < len; ++i)
    {
      array2[i * 2] = array1[i * 2];
      array2[i * 2 + 1] = -array1[i * 2 + 1];
    }
}

static void
math_varray_abs_cplx(const size_t len, double missval1, const Varray<double> &array1, Varray<double> &array2)
{
  auto missval2 = missval1;
  for (size_t i = 0; i < len; ++i)
    {
      array2[i] = SQRTMN(ADDMN(MULMN(array1[2 * i], array1[2 * i]), MULMN(array1[2 * i + 1], array1[2 * i + 1])));
    }
}

static void
math_varray_arg_cplx(const size_t len, double missval1, const Varray<double> &array1, Varray<double> &array2)
{
  for (size_t i = 0; i < len; ++i)
    {
      array2[i] = (DBL_IS_EQUAL(array1[2 * i], missval1) || DBL_IS_EQUAL(array1[2 * i + 1], missval1))
                      ? missval1
                      : std::atan2(array1[2 * i + 1], array1[2 * i]);
    }
}

enum struct Oper
{
  Abs,
  Int,
  Nint,
  Sqr,
  Sqrt,
  Exp,
  Ln,
  Log10,
  Sin,
  Cos,
  Tan,
  Asin,
  Acos,
  Atan,
  Pow,
  Rand,
  Reci,
  Not,
  Conj,
  Re,
  Im,
  Arg
};

static void
addOperators(void)
{
  // clang-format off
  cdo_operator_add("abs",   (int)Oper::Abs,   0, nullptr);
  cdo_operator_add("int",   (int)Oper::Int,   0, nullptr);
  cdo_operator_add("nint",  (int)Oper::Nint,  0, nullptr);
  cdo_operator_add("sqr",   (int)Oper::Sqr,   0, nullptr);
  cdo_operator_add("sqrt",  (int)Oper::Sqrt,  0, nullptr);
  cdo_operator_add("exp",   (int)Oper::Exp,   0, nullptr);
  cdo_operator_add("ln",    (int)Oper::Ln,    0, nullptr);
  cdo_operator_add("log10", (int)Oper::Log10, 0, nullptr);
  cdo_operator_add("sin",   (int)Oper::Sin,   0, nullptr);
  cdo_operator_add("cos",   (int)Oper::Cos,   0, nullptr);
  cdo_operator_add("tan",   (int)Oper::Tan,   0, nullptr);
  cdo_operator_add("asin",  (int)Oper::Asin,  0, nullptr);
  cdo_operator_add("acos",  (int)Oper::Acos,  0, nullptr);
  cdo_operator_add("atan",  (int)Oper::Atan,  0, nullptr);
  cdo_operator_add("pow",   (int)Oper::Pow,   0, nullptr);
  cdo_operator_add("rand",  (int)Oper::Rand,  0, nullptr);
  cdo_operator_add("reci",  (int)Oper::Reci,  0, nullptr);
  cdo_operator_add("not",   (int)Oper::Not,   0, nullptr);
  cdo_operator_add("conj",  (int)Oper::Conj,  0, nullptr);
  cdo_operator_add("re",    (int)Oper::Re,    0, nullptr);
  cdo_operator_add("im",    (int)Oper::Im,    0, nullptr);
  cdo_operator_add("arg",   (int)Oper::Arg,   0, nullptr);
  // clang-format on
}
class ModuleMath
{
private:
  Oper operfunc;

  CdoStreamID streamID1;
  CdoStreamID streamID2;

  int taxisID1;
  int taxisID2;
  int vlistID1;
  int vlistID2;

  VarList varList1;
  Varray<double> array1;
  Varray<double> array2;

  double rc = 0.0;

public:
  void
  init(void *process)
  {

    cdo_initialize(process);

    addOperators();

    auto operatorID = cdo_operator_id();
    operfunc = (Oper) cdo_operator_f1(operatorID);

    if (operfunc == Oper::Pow)
      {
        operator_input_arg("value");
        rc = parameter_to_double(cdo_operator_argv(0));
      }
    else { operator_check_argc(0); }

    if (operfunc == Oper::Rand) std::srand(Options::Random_Seed);

    streamID1 = cdo_open_read(0);

    vlistID1 = cdo_stream_inq_vlist(streamID1);
    vlistID2 = vlistDuplicate(vlistID1);

    if (operfunc == Oper::Re || operfunc == Oper::Im || operfunc == Oper::Abs || operfunc == Oper::Arg)
      {
        auto nvars = vlistNvars(vlistID2);
        for (int varID = 0; varID < nvars; ++varID)
          {
            if (vlistInqVarDatatype(vlistID2, varID) == CDI_DATATYPE_CPX32)
              vlistDefVarDatatype(vlistID2, varID, CDI_DATATYPE_FLT32);
            if (vlistInqVarDatatype(vlistID2, varID) == CDI_DATATYPE_CPX64)
              vlistDefVarDatatype(vlistID2, varID, CDI_DATATYPE_FLT64);
          }
      }

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    auto gridsizemax = vlistGridsizeMax(vlistID1);
    if (vlistNumber(vlistID1) != CDI_REAL) gridsizemax *= 2;

    array1 = Varray<double>(gridsizemax);
    array2 = Varray<double>(gridsizemax);

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);

    varListInit(varList1, vlistID1);
  }
  void
  run()
  {
    int tsID = 0;
    while (true)
      {
        auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
        if (nrecs == 0) break;

        cdo_taxis_copy_timestep(taxisID2, taxisID1);
        cdo_def_timestep(streamID2, tsID);

        for (int recID = 0; recID < nrecs; ++recID)
          {
            int varID, levelID;
            cdo_inq_record(streamID1, &varID, &levelID);
            size_t nmiss;
            cdo_read_record(streamID1, &array1[0], &nmiss);

            auto missval1 = varList1[varID].missval;
            auto n = varList1[varID].gridsize;
            auto number = varList1[varID].nwpv;

            if (number == CDI_REAL)
              {
                // clang-format off
              switch (operfunc)
                {
                case Oper::Abs:   math_varray_transform(nmiss, n, missval1, array1, array2, unary_op_abs); break;
                case Oper::Int:   math_varray_transform(nmiss, n, missval1, array1, array2, unary_op_int); break;
                case Oper::Nint:  math_varray_transform(nmiss, n, missval1, array1, array2, unary_op_nint); break;
                case Oper::Sqr:   math_varray_transform(nmiss, n, missval1, array1, array2, unary_op_sqr); break;
                case Oper::Sqrt:
                  for (size_t i = 0; i < n; ++i) array2[i] = DBL_IS_EQUAL(array1[i], missval1) ? missval1 : SQRTMN(array1[i]);
                  break;
                case Oper::Exp:   math_varray_transform(nmiss, n, missval1, array1, array2, unary_op_exp); break;
                case Oper::Ln:    check_lower_range(nmiss, n, missval1, array1, -1);
                                  math_varray_transform(nmiss, n, missval1, array1, array2, unary_op_log); break;
                case Oper::Log10: check_lower_range(nmiss, n, missval1, array1, -1);
                                  math_varray_transform(nmiss, n, missval1, array1, array2, unary_op_log10); break;
                case Oper::Sin:   math_varray_transform(nmiss, n, missval1, array1, array2, unary_op_sin); break;
                case Oper::Cos:   math_varray_transform(nmiss, n, missval1, array1, array2, unary_op_cos); break;
                case Oper::Tan:   math_varray_transform(nmiss, n, missval1, array1, array2, unary_op_tan); break;
                case Oper::Asin:  check_out_of_range(nmiss, n, missval1, array1, -1, 1);
                                  math_varray_transform(nmiss, n, missval1, array1, array2, unary_op_asin); break;
                case Oper::Acos:  check_out_of_range(nmiss, n, missval1, array1, -1, 1);
                                  math_varray_transform(nmiss, n, missval1, array1, array2, unary_op_acos); break;
                case Oper::Atan:  math_varray_transform(nmiss, n, missval1, array1, array2, unary_op_atan); break;
                case Oper::Pow:
                  for (size_t i = 0; i < n; ++i) array2[i] = DBL_IS_EQUAL(array1[i], missval1) ? missval1 : std::pow(array1[i], rc);
                  break;
                case Oper::Rand:
                  for (size_t i = 0; i < n; ++i) array2[i] = DBL_IS_EQUAL(array1[i], missval1) ? missval1 : ((double) std::rand()) / ((double) RAND_MAX);
                  break;
                case Oper::Reci:  math_varray_transform(nmiss, n, missval1, array1, array2, unary_op_reci); break;
                case Oper::Not:   math_varray_transform(nmiss, n, missval1, array1, array2, unary_op_not); break;
                case Oper::Re:
                case Oper::Arg:   math_varray_transform(nmiss, n, missval1, array1, array2, unary_op_nop); break;
                default: cdo_abort("Operator not implemented for real data!"); break;
                }
                // clang-format on

                nmiss = varray_num_mv(n, array2, missval1);
              }
            else
              {
                // clang-format off
              switch (operfunc)
                {
                case Oper::Sqr:   math_varray_sqr_cplx(n, array1, array2); break;
                case Oper::Sqrt:  math_varray_sqrt_cplx(n, missval1, array1, array2); break;
                case Oper::Conj:  math_varray_conj_cplx(n, array1, array2); break;
                case Oper::Re:    for (size_t i = 0; i < n; ++i) array2[i] = array1[i * 2]; break;
                case Oper::Im:    for (size_t i = 0; i < n; ++i) array2[i] = array1[i * 2 + 1]; break;
                case Oper::Abs:   math_varray_abs_cplx(n, missval1, array1, array2); break;
                case Oper::Arg:   math_varray_arg_cplx(n, missval1, array1, array2); break;
                default: cdo_abort("Fields with complex numbers are not supported by this operator!"); break;
                }
                // clang-format on

                nmiss = 0;
              }

            cdo_def_record(streamID2, varID, levelID);
            cdo_write_record(streamID2, array2.data(), nmiss);
          }

        tsID++;
      }
  }
  void
  close()
  {
    cdo_stream_close(streamID2);
    cdo_stream_close(streamID1);

    vlistDestroy(vlistID2);

    cdo_finish();
  }
};

void *
Math(void *process)
{
  ModuleMath math;
  math.init(process);
  math.run();
  math.close();
  return nullptr;
}
