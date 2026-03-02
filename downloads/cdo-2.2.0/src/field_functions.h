/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/
#ifndef FIELD_FUNCTIONS_H
#define FIELD_FUNCTIONS_H

#include "field.h"

double var_to_std(double rvar, double missval);

enum FieldFunc
{
  FieldFunc_Min = 100,
  FieldFunc_Max,
  FieldFunc_Range,
  FieldFunc_Sum,
  FieldFunc_Avg,
  FieldFunc_Mean,
  FieldFunc_Var,
  FieldFunc_Var1,
  FieldFunc_Std,
  FieldFunc_Std1,
  FieldFunc_Skew,
  FieldFunc_Kurt,
  FieldFunc_Median,
  FieldFunc_Count,
  FieldFunc_Pctl,

  FieldFunc_Cor,
  FieldFunc_Covar,
  FieldFunc_Avgw,
  FieldFunc_Meanw,
  FieldFunc_Stdw,
  FieldFunc_Std1w,
  FieldFunc_Varw,
  FieldFunc_Var1w,
  FieldFunc_Minidx,
  FieldFunc_Maxidx,
  FieldFunc_Rmsd,

  FieldFunc_Add,
  FieldFunc_Sub,
  FieldFunc_Mul,
  FieldFunc_Div,
  FieldFunc_Mod,

  FieldFunc_EQ,
  FieldFunc_NE,
  FieldFunc_LE,
  FieldFunc_LT,
  FieldFunc_GE,
  FieldFunc_GT,

  FieldFunc_Atan2,
  FieldFunc_Setmiss,
};

// fieldmem.cc
void fields_from_vlist(int vlistID, FieldVector2D &field2D);
void fields_from_vlist(int vlistID, FieldVector2D &field2D, int ptype);
void fields_from_vlist(int vlistID, FieldVector2D &field2D, int ptype, double fillValue);

// field.cc
double field_function(const Field &field, int function);

double field_min(const Field &field);
double field_range(const Field &field);
double field_max(const Field &field);
double field_sum(const Field &field);
double field_mean(const Field &field);
double field_meanw(const Field &field);
double field_avg(const Field &field);
double field_avgw(const Field &field);
double field_std(const Field &field);
double field_std1(const Field &field);
double field_var(const Field &field);
double field_var1(const Field &field);
double field_stdw(const Field &field);
double field_std1w(const Field &field);
double field_varw(const Field &field);
double field_var1w(const Field &field);
double field_skew(const Field &field);
double field_kurt(const Field &field);
double field_median(const Field &field);
double field_count(const Field &field);

// ENS validation
double field_rank(Field &field);

double field_pctl(Field &field, double pn);

// field_zonal.cc
void zonal_function(const Field &field1, Field &field2, int function);
void zonal_min(const Field &field1, Field &field2);
void zonal_max(const Field &field1, Field &field2);
void zonal_range(const Field &field1, Field &field2);
void zonal_sum(const Field &field1, Field &field2);
void zonal_avg(const Field &field1, Field &field2);
void zonal_mean(const Field &field1, Field &field2);
void zonal_std(const Field &field1, Field &field2);
void zonal_std1(const Field &field1, Field &field2);
void zonal_var(const Field &field1, Field &field2);
void zonal_var1(const Field &field1, Field &field2);
void zonal_skew(const Field &field1, Field &field2);
void zonal_kurt(const Field &field1, Field &field2);
void zonal_median(const Field &field1, Field &field2);
void zonal_pctl(const Field &field1, Field &field2, double pn);

// field_meridional.cc
void meridional_function(const Field &field1, Field &field2, int function);
void meridional_pctl(const Field &field1, Field &field2, double pn);

void field_rms(const Field &field1, const Field &field2, Field &field3);

// fieldc.cc
void fieldc_function(Field &field, double rconst, int function);

void fieldc_mul(Field &field, double rconst);
void fieldc_div(Field &field, double rconst);
void fieldc_add(Field &field, double rconst);
void fieldc_sub(Field &field, double rconst);
void fieldc_min(Field &field, double rconst);
void fieldc_max(Field &field, double rconst);
void fieldc_mod(Field &field, double divisor);

// fieldc_complex.cc
void fieldc_function_complex(Field &field, const double rconstcplx[2], int function);

// field2.cc
void field2_function(Field &field1, const Field &field2, int function);

void field2_add(Field &field1, const Field &field2);
void field2_sum(Field &field1, const Field &field2);
void field2_sub(Field &field1, const Field &field2);
void field2_mul(Field &field1, const Field &field2);
void field2_div(Field &field1, const Field &field2);
void field2_min(Field &field1, const Field &field2);
void field2_max(Field &field1, const Field &field2);
void field2_atan2(Field &field1, const Field &field2);

void field2_sumq(Field &field1, const Field &field2);
void field2_sumw(Field &field1, const Field &field2, double w);
void field2_sumqw(Field &field1, const Field &field2, double w);
void field2_sumtr(Field &field1, const Field &field2, double refval);
void field2_vinit(Field &field1, const Field &field2);
void field2_vincr(Field &field1, const Field &field2);
void field2_vinit(Field &field1, const Field &field2, int vinit);
void field2_vincr(Field &field1, const Field &field2, int vincr);

void field2_sumsumq(Field &field1, Field &field2, const Field &field3);
void field2_maxmin(Field &field1, Field &field2, const Field &field3);
void field2_minidx(Field &field1, Field &field2, const Field &field3, int idx);
void field2_maxidx(Field &field1, Field &field2, const Field &field3, int idx);
void field2_var(Field &field1, const Field &field2, const Field &field3, int divisor);
void field2_std(Field &field1, const Field &field2, const Field &field3, int divisor);
void fieldc_var(Field &field1, const Field &field2, int numSets, int divisor);
void fieldc_std(Field &field1, const Field &field2, int numSets, int divisor);
void field2_moq(Field &field1, const Field &field2);
void field2_moqw(Field &field1, const Field &field2, double w);

void field2_count(Field &field1, const Field &field2);

// field2_complex.cc
void field2_function_complex(Field &field1, const Field &field2, int function);

#endif /* FIELD_FUNCTIONS_H */
