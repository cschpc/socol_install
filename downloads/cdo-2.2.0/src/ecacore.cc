/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Copyright (C) 2006 Brockmann Consult

  Author: Ralf Quast
          Uwe Schulzweida

*/

#include <assert.h>
#include <string.h>

#include <cdi.h>

#include "cdo_vlist.h"
#include <mpim_grid.h>
#include "process_int.h"
#include "ecacore.h"
#include "ecautil.h"
#include "util_files.h"
#include "util_date.h"
#include "datetime.h"
#include "field_functions.h"

constexpr int FIRST_VAR_ID = 0;

#define IS_NOT_SET(x) (x == nullptr)
#define IS_SET(x) (x != nullptr)

static void
init_field(Field &field, int gridID, double missval)
{
  field.grid = gridID;
  field.missval = missval;
}

static void
init_field(Field &field, int gridID, double missval, size_t gridsize)
{
  field.grid = gridID;
  field.missval = missval;
  field.resize(gridsize);
}

static void
init_field(Field &field, int gridID, double missval, size_t gridsize, double value)
{
  field.grid = gridID;
  field.missval = missval;
  field.resize(gridsize, value);
}

void
eca1(const ECA_REQUEST_1 &request)
{
  const auto operatorID = cdo_operator_id();

  CdiDateTime ovDateTime{};
  CdiDateTime inDateTime21{};
  int itsID;
  int otsID;

  const auto compareDate = cdo_operator_f2(operatorID);

  const auto istreamID = cdo_open_read(0);

  const auto ivlistID = cdo_stream_inq_vlist(istreamID);
  const auto ovlistID = vlistCreate();

  const auto gridID = vlistInqVarGrid(ivlistID, FIRST_VAR_ID);
  const auto zaxisID = vlistInqVarZaxis(ivlistID, FIRST_VAR_ID);
  const auto missval = vlistInqVarMissval(ivlistID, FIRST_VAR_ID);

  auto varID = vlistDefVar(ovlistID, gridID, zaxisID, TIME_VARYING);

  vlistDefVarMissval(ovlistID, varID, missval);

  if (IS_SET(request.var1.name)) cdiDefKeyString(ovlistID, varID, CDI_KEY_NAME, request.var1.name);
  if (IS_SET(request.var1.longname)) cdiDefKeyString(ovlistID, varID, CDI_KEY_LONGNAME, request.var1.longname);
  if (IS_SET(request.var1.units)) cdiDefKeyString(ovlistID, varID, CDI_KEY_UNITS, request.var1.units);

  if (IS_SET(request.var2.h2) || IS_SET(request.var2.h3))
    {
      varID = vlistDefVar(ovlistID, gridID, zaxisID, TIME_VARYING);

      vlistDefVarMissval(ovlistID, varID, missval);

      if (IS_SET(request.var2.name)) cdiDefKeyString(ovlistID, varID, CDI_KEY_NAME, request.var2.name);
      if (IS_SET(request.var2.longname)) cdiDefKeyString(ovlistID, varID, CDI_KEY_LONGNAME, request.var2.longname);
      if (IS_SET(request.var2.units)) cdiDefKeyString(ovlistID, varID, CDI_KEY_UNITS, request.var2.units);
    }

  if (cdo_operator_f2(operatorID) == 16) vlistDefNtsteps(ovlistID, 1);

  const auto itaxisID = vlistInqTaxis(ivlistID);
  const auto otaxisID = cdo_taxis_create(TAXIS_RELATIVE);
  taxisDefTunit(otaxisID, TUNIT_DAY);
  //  taxisDefTunit(otaxisID, TUNIT_MINUTE);
  //  taxisDefCalendar(otaxisID, CALENDAR_PROLEPTIC);
  taxisDefCalendar(otaxisID, taxisInqCalendar(itaxisID));
  taxisDefRdate(otaxisID, request.var1.refdate);
  taxisDefRtime(otaxisID, 0);
  vlistDefTaxis(ovlistID, otaxisID);

  const auto ostreamID = cdo_open_write(1);
  cdo_def_vlist(ostreamID, ovlistID);

  const auto gridsize = gridInqSize(gridID);

  Field field1, field2, field3;
  field1.resize(gridsize);
  field3.resize(gridsize);
  if (IS_SET(request.var2.h2) || IS_SET(request.var2.h3)) field2.resize(gridsize);

  const auto nlevels = zaxisInqSize(zaxisID);

  FieldVector var12(nlevels), samp1(nlevels), samp2(nlevels);
  FieldVector var13, var21, var23;
  if (IS_SET(request.var1.f3)) var13.resize(nlevels);
  if (IS_SET(request.var2.h2)) var21.resize(nlevels);
  if (IS_SET(request.var2.h3)) var23.resize(nlevels);

  for (int levelID = 0; levelID < nlevels; ++levelID)
    {
      init_field(var12[levelID], gridID, missval, gridsize);
      init_field(samp1[levelID], gridID, missval, gridsize);
      init_field(samp2[levelID], gridID, missval);

      if (IS_SET(request.var1.f3)) init_field(var13[levelID], gridID, missval, gridsize);
      if (IS_SET(request.var2.h2)) init_field(var21[levelID], gridID, missval, gridsize);
      if (IS_SET(request.var2.h3)) init_field(var23[levelID], gridID, missval, gridsize);
    }

  itsID = 0;
  otsID = 0;
  while (true)
    {
      int nrecs = 0;
      long numSets = 0;
      while (true)
        {
          nrecs = cdo_stream_inq_timestep(istreamID, itsID);
          if (nrecs == 0) break;

          const auto ivDateTime = taxisInqVdatetime(itaxisID);

          if (numSets == 0) inDateTime21 = ivDateTime;

          if (date_is_neq(ivDateTime, inDateTime21, compareDate))
            {
              cdo_add_steps(-1);
              break;
            }

          for (int recID = 0; recID < nrecs; ++recID)
            {
              int levelID;
              cdo_inq_record(istreamID, &varID, &levelID);

              if (varID != FIRST_VAR_ID) continue;

              if (numSets == 0)
                {
                  if (request.var1.f2 != &vfarnum2)
                    {
                      field_fill(var12[levelID], missval);
                      var12[levelID].nmiss = gridsize;
                    }
                  field_fill(samp1[levelID], missval);
                  if (!samp2[levelID].empty()) field_fill(samp2[levelID], 0.0);
                  if (IS_SET(request.var1.f3)) field_fill(var13[levelID], missval);
                  if (IS_SET(request.var2.h2)) field_fill(var21[levelID], missval);
                  if (IS_SET(request.var2.h3)) field_fill(var23[levelID], missval);

                  samp1[levelID].nmiss = gridsize;
                  if (IS_SET(request.var1.f3)) var13[levelID].nmiss = gridsize;
                  if (IS_SET(request.var2.h2)) var21[levelID].nmiss = gridsize;
                  if (IS_SET(request.var2.h3)) var23[levelID].nmiss = gridsize;
                }

              cdo_read_record(istreamID, field1.vec_d.data(), &field1.nmiss);
              field1.grid = var12[levelID].grid;
              field1.missval = var12[levelID].missval;

              vfarnum(samp1[levelID], field1);

              if (IS_SET(request.var2.h2))
                {
                  field2.vec_d = field1.vec_d;
                  field2.nmiss = field1.nmiss;
                  field2.grid = field1.grid;
                  field2.missval = field1.missval;
                }

              if (IS_SET(request.var1.f1)) request.var1.f1(field1, request.var1.f1arg);

              if (field1.nmiss || !samp2[levelID].empty())
                {
                  if (samp2[levelID].empty())
                    {
                      samp2[levelID].resize(gridsize);
                      field_fill(samp2[levelID], numSets);
                    }
                  for (size_t i = 0; i < gridsize; ++i)
                    {
                      if (DBL_IS_EQUAL(field1.vec_d[i], field1.missval)) continue;
                      samp2[levelID].vec_d[i]++;
                    }
                }

              if (IS_NOT_EQUAL(request.var1.mulc, 0.0)) fieldc_mul(field1, request.var1.mulc);
              if (IS_NOT_EQUAL(request.var1.addc, 0.0)) fieldc_add(field1, request.var1.addc);

              if (IS_SET(request.var1.f3) && request.var1.f2 == &vfarnum2)
                {
                  varray_copy(gridsize, var12[levelID].vec_d, field3.vec_d);
                  field3.nmiss = var12[levelID].nmiss;
                  field3.grid = var12[levelID].grid;
                  field3.missval = var12[levelID].missval;
                }

              request.var1.f2(var12[levelID], field1);

              if (IS_SET(request.var2.h2) || IS_SET(request.var2.h3))
                {
                  // if h2 is null, use the output of f2 as input for h1
                  if (IS_NOT_SET(request.var2.h2))
                    {
                      varray_copy(gridsize, var12[levelID].vec_d, field2.vec_d);
                      field2.nmiss = var12[levelID].nmiss;
                      field2.grid = var12[levelID].grid;
                      field2.missval = var12[levelID].missval;
                    }

                  if (IS_SET(request.var2.h1)) request.var2.h1(field2, request.var2.h1arg);

                  if (IS_NOT_SET(request.var2.h2))
                    request.var2.h3(var23[levelID], field2);
                  else
                    {
                      request.var2.h2(var21[levelID], field2);
                      if (IS_SET(request.var2.h3)) request.var2.h3(var23[levelID], var21[levelID]);
                    }
                }

              if (IS_SET(request.var1.f3))
                {
                  if (request.var1.f2 == &vfarnum2)
                    {
                      auto &array1 = field3.vec_d;
                      const auto &array2 = var12[levelID].vec_d;
                      const auto missval2 = field2.missval;

                      const auto len = field1.size;
                      if (len != field2.size) cdo_abort("Fields have different size (%s)", __func__);

                      if (field2.nmiss)
                        {
                          for (size_t i = 0; i < len; ++i)
                            {
                              if (!DBL_IS_EQUAL(array2[i], missval2) && !DBL_IS_EQUAL(array2[i], 0.0)) array1[i] = 0.0;
                            }
                        }
                      else
                        {
                          for (size_t i = 0; i < len; ++i)
                            {
                              if (!DBL_IS_EQUAL(array2[i], 0.0)) array1[i] = 0.0;
                            }
                        }
                      request.var1.f3(var13[levelID], field3);
                    }
                  else
                    request.var1.f3(var13[levelID], var12[levelID]);
                }
            }

          ovDateTime = ivDateTime;
          numSets++;
          itsID++;
        }

      if (request.var1.f2 == &vfarnum2)
        {
          for (int levelID = 0; levelID < nlevels; ++levelID)
            {
              if (nrecs == 0) { request.var1.f3(var13[levelID], var12[levelID]); }
              else
                {
                  auto &array1 = var13[levelID].vec_d;
                  const auto &array2 = var12[levelID].vec_d;
                  const auto missval2 = field2.missval;
                  const auto len = field1.size;

                  for (size_t i = 0; i < len; ++i)
                    {
                      if (DBL_IS_EQUAL(array2[i], numSets) || array2[i] > numSets) array1[i] = missval2;
                    }
                }
            }
        }

      if (nrecs == 0 && numSets == 0) break;

      if (request.var1.epilog == MEAN || request.var1.epilog == PERCENT_OF_TIME)
        for (int levelID = 0; levelID < nlevels; ++levelID)
          {
            auto &rvar = IS_SET(request.var1.f3) ? var13[levelID] : var12[levelID];

            if (samp2[levelID].empty())
              fieldc_div(rvar, numSets);
            else
              field2_div(rvar, samp2[levelID]);

            if (request.var1.epilog == PERCENT_OF_TIME) fieldc_mul(rvar, 100.0);
          }

      if (request.var1.refdate == 19550101)
        {
          taxisDefVdatetime(otaxisID, ovDateTime);
          cdo_def_timestep(ostreamID, otsID);
        }
      else
        {
          int year, month, day;
          cdiDate_decode(inDateTime21.date, &year, &month, &day);
          define_mid_of_time(cdo_operator_f2(operatorID), otaxisID, year, month, 12);
          cdo_def_timestep(ostreamID, otsID);
        }

      if (otsID && vlistInqVarTimetype(ivlistID, FIRST_VAR_ID) == TIME_CONSTANT) continue;

      varID = 0;
      for (int levelID = 0; levelID < nlevels; ++levelID)
        {
          auto &rvar = IS_SET(request.var1.f3) ? var13[levelID] : var12[levelID];

          vfarsel(rvar, samp1[levelID]);

          cdo_def_record(ostreamID, varID, levelID);
          cdo_write_record(ostreamID, rvar.vec_d.data(), rvar.nmiss);
        }

      if (IS_SET(request.var2.h2) || IS_SET(request.var2.h3))
        {
          varID = 1;
          for (int levelID = 0; levelID < nlevels; ++levelID)
            {
              auto &rvar = IS_SET(request.var2.h3) ? var23[levelID] : var21[levelID];

              vfarsel(rvar, samp1[levelID]);

              cdo_def_record(ostreamID, varID, levelID);
              cdo_write_record(ostreamID, rvar.vec_d.data(), rvar.nmiss);
            }
        }

      if (nrecs == 0) break;
      otsID++;
    }

  cdo_stream_close(ostreamID);
  cdo_stream_close(istreamID);
}

void
eca2(const ECA_REQUEST_2 &request)
{
  const auto operatorID = cdo_operator_id();

  CdiDateTime ovDateTime{};
  CdiDateTime inDateTime21{};
  int varID;
  int itsID;
  int otsID;

  const auto compareDate = cdo_operator_f2(operatorID);

  const auto istreamID1 = cdo_open_read(0);
  const auto istreamID2 = cdo_open_read(1);

  const auto ivlistID1 = cdo_stream_inq_vlist(istreamID1);
  const auto ivlistID2 = cdo_stream_inq_vlist(istreamID2);
  const auto ovlistID = vlistCreate();

  vlist_compare(ivlistID1, ivlistID2, CMP_ALL);

  const auto gridID = vlistInqVarGrid(ivlistID1, FIRST_VAR_ID);
  const auto zaxisID = vlistInqVarZaxis(ivlistID1, FIRST_VAR_ID);
  const auto missval1 = vlistInqVarMissval(ivlistID1, FIRST_VAR_ID);
  const auto missval2 = vlistInqVarMissval(ivlistID2, FIRST_VAR_ID);

  varID = vlistDefVar(ovlistID, gridID, zaxisID, TIME_VARYING);

  vlistDefVarMissval(ovlistID, varID, missval1);

  if (IS_SET(request.var1.name)) cdiDefKeyString(ovlistID, varID, CDI_KEY_NAME, request.var1.name);
  if (IS_SET(request.var1.longname)) cdiDefKeyString(ovlistID, varID, CDI_KEY_LONGNAME, request.var1.longname);
  if (IS_SET(request.var1.units)) cdiDefKeyString(ovlistID, varID, CDI_KEY_UNITS, request.var1.units);

  if (IS_SET(request.var2.h2))
    {
      varID = vlistDefVar(ovlistID, gridID, zaxisID, TIME_VARYING);

      vlistDefVarMissval(ovlistID, varID, missval1);

      if (IS_SET(request.var2.name)) cdiDefKeyString(ovlistID, varID, CDI_KEY_NAME, request.var2.name);
      if (IS_SET(request.var2.longname)) cdiDefKeyString(ovlistID, varID, CDI_KEY_LONGNAME, request.var2.longname);
      if (IS_SET(request.var2.units)) cdiDefKeyString(ovlistID, varID, CDI_KEY_UNITS, request.var2.units);
    }

  if (cdo_operator_f2(operatorID) == 16) vlistDefNtsteps(ovlistID, 1);

  const auto itaxisID1 = vlistInqTaxis(ivlistID1);
  const auto otaxisID = cdo_taxis_create(TAXIS_RELATIVE);
  taxisDefTunit(otaxisID, TUNIT_DAY);
  //  taxisDefTunit(otaxisID, TUNIT_MINUTE);
  //  taxisDefCalendar(otaxisID, CALENDAR_PROLEPTIC);
  taxisDefCalendar(otaxisID, taxisInqCalendar(itaxisID1));
  taxisDefRdate(otaxisID, request.var1.refdate);
  taxisDefRtime(otaxisID, 0);
  vlistDefTaxis(ovlistID, otaxisID);

  const auto ostreamID = cdo_open_write(2);
  cdo_def_vlist(ostreamID, ovlistID);

  const auto gridsize = gridInqSize(gridID);

  Field field1, field2, field3;
  field1.resize(gridsize);
  field2.resize(gridsize);
  field3.resize(gridsize);
  constexpr int MaxDays = 373;
  FieldVector2D vars2[MaxDays];

  const auto nlevels = zaxisInqSize(zaxisID);

  FieldVector var14(nlevels), samp1(nlevels), samp2(nlevels), samp3(nlevels);
  FieldVector total, var15, var22;
  if (request.var1.epilog == PERCENT_OF_TOTAL_AMOUNT) total.resize(nlevels);
  if (IS_SET(request.var1.f5)) var15.resize(nlevels);
  if (IS_SET(request.var2.h2)) var22.resize(nlevels);

  for (int levelID = 0; levelID < nlevels; ++levelID)
    {
      init_field(var14[levelID], gridID, missval1, gridsize);
      init_field(samp1[levelID], gridID, missval1, gridsize);
      init_field(samp2[levelID], gridID, missval1, gridsize);
      init_field(samp3[levelID], gridID, missval1);

      if (request.var1.epilog == PERCENT_OF_TOTAL_AMOUNT) init_field(total[levelID], gridID, missval1, gridsize);
      if (IS_SET(request.var1.f5)) init_field(var15[levelID], gridID, missval1, gridsize);
      if (IS_SET(request.var2.h2)) init_field(var22[levelID], gridID, missval1, gridsize);
    }

  itsID = 0;
  while (true)
    {
      const auto nrecs = cdo_stream_inq_timestep(istreamID2, itsID);
      if (nrecs == 0) break;

      const auto ivDateTime = taxisInqVdatetime(vlistInqTaxis(ivlistID2));

      const auto dayOfYear = decode_day_of_year(ivDateTime.date);
      if (dayOfYear < 0 || dayOfYear >= MaxDays) cdo_abort("Day %d out of range!", dayOfYear);

      if (!vars2[dayOfYear].size()) fields_from_vlist(ivlistID2, vars2[dayOfYear], FIELD_VEC);

      for (int recID = 0; recID < nrecs; ++recID)
        {
          int levelID;
          cdo_inq_record(istreamID2, &varID, &levelID);
          if (varID != FIRST_VAR_ID) continue;
          auto &rvar = vars2[dayOfYear][0][levelID];
          cdo_read_record(istreamID2, rvar.vec_d.data(), &rvar.nmiss);
        }

      itsID++;
    }

  itsID = 0;
  otsID = 0;

  while (true)
    {
      int nrecs = 0;
      long numSets = 0;
      while (true)
        {
          nrecs = cdo_stream_inq_timestep(istreamID1, itsID);
          if (nrecs == 0) break;

          const auto ivDateTime = taxisInqVdatetime(itaxisID1);

          const auto dayOfYear = decode_day_of_year(ivDateTime.date);
          if (!vars2[dayOfYear].size()) cdo_abort("Input streams have different time values!");

          if (numSets == 0) inDateTime21 = ivDateTime;

          if (date_is_neq(ivDateTime, inDateTime21, compareDate))
            {
              cdo_add_steps(-1);
              break;
            }

          for (int recID = 0; recID < nrecs; ++recID)
            {
              int levelID;
              cdo_inq_record(istreamID1, &varID, &levelID);

              if (varID != FIRST_VAR_ID) continue;

              if (numSets == 0)
                {
                  if (request.var1.f4 != &vfarnum2)
                    {
                      field_fill(var14[levelID], missval1);
                      var14[levelID].nmiss = gridsize;
                    }
                  field_fill(samp1[levelID], missval1);
                  field_fill(samp2[levelID], missval1);
                  if (!samp3[levelID].empty()) field_fill(samp3[levelID], 0.0);
                  if (request.var1.epilog == PERCENT_OF_TOTAL_AMOUNT) field_fill(total[levelID], 0.0);
                  if (IS_SET(request.var1.f5)) field_fill(var15[levelID], missval1);
                  if (IS_SET(request.var2.h2)) field_fill(var22[levelID], missval1);

                  samp1[levelID].nmiss = gridsize;
                  samp2[levelID].nmiss = gridsize;
                  if (request.var1.epilog == PERCENT_OF_TOTAL_AMOUNT) total[levelID].nmiss = gridsize;
                  if (IS_SET(request.var1.f5)) var15[levelID].nmiss = gridsize;
                  if (IS_SET(request.var2.h2)) var22[levelID].nmiss = gridsize;
                }

              cdo_read_record(istreamID1, field1.vec_d.data(), &field1.nmiss);
              field1.grid = gridID;
              field1.missval = missval1;
              field2.grid = vars2[dayOfYear][0][levelID].grid;
              field2.nmiss = vars2[dayOfYear][0][levelID].nmiss;
              field2.missval = missval2;
              field_copy(vars2[dayOfYear][0][levelID], field2);

              vfarnum(samp1[levelID], field1);
              vfarnum(samp2[levelID], field2);

              if (request.var1.epilog == PERCENT_OF_TOTAL_AMOUNT) field2_sum(total[levelID], field1);

              if (IS_SET(request.var1.f1)) request.var1.f1(field1, request.var1.f1arg);
              if (IS_SET(request.var1.f2)) request.var1.f2(field2, request.var1.f2arg);

              if (field1.nmiss || !samp3[levelID].empty())
                {
                  if (samp3[levelID].empty())
                    {
                      samp3[levelID].resize(gridsize);
                      field_fill(samp3[levelID], numSets);
                    }
                  for (size_t i = 0; i < gridsize; ++i)
                    {
                      if (DBL_IS_EQUAL(field1.vec_d[i], field1.missval)) continue;
                      samp3[levelID].vec_d[i]++;
                    }
                }

              if (IS_SET(request.var1.f5) && request.var1.f4 == &vfarnum2)
                {
                  varray_copy(gridsize, var14[levelID].vec_d, field3.vec_d);
                  field3.nmiss = var14[levelID].nmiss;
                  field3.grid = var14[levelID].grid;
                  field3.missval = var14[levelID].missval;
                }

              request.var1.f3(field1, field2);
              request.var1.f4(var14[levelID], field1);

              if (IS_SET(request.var2.h2))
                {
                  varray_copy(gridsize, var14[levelID].vec_d, field2.vec_d);
                  field2.nmiss = var14[levelID].nmiss;
                  field2.grid = var14[levelID].grid;
                  field2.missval = var14[levelID].missval;

                  if (IS_SET(request.var2.h1)) request.var2.h1(field2, request.var2.h1arg);

                  request.var2.h2(var22[levelID], field2);
                }

              if (IS_SET(request.var1.f5))
                {
                  if (request.var1.f4 == &vfarnum2)
                    {
                      auto &array1 = field3.vec_d;
                      const auto &array2 = var14[levelID].vec_d;
                      const auto missvaltemp = field1.missval;

                      const auto len = field1.size;
                      if (len != field3.size) cdo_abort("Fields have different size (%s)", __func__);

                      if (field1.nmiss)
                        {
                          for (size_t i = 0; i < len; ++i)
                            {
                              if (!DBL_IS_EQUAL(array2[i], missvaltemp) && !DBL_IS_EQUAL(array2[i], 0.0)) array1[i] = 0.0;
                            }
                        }
                      else
                        {
                          for (size_t i = 0; i < len; ++i)
                            {
                              if (!DBL_IS_EQUAL(array2[i], 0.0)) array1[i] = 0.0;
                            }
                        }
                      request.var1.f5(var15[levelID], field3, request.var1.f5arg);
                    }
                  else
                    request.var1.f5(var15[levelID], var14[levelID], request.var1.f5arg);
                }
            }

          ovDateTime = ivDateTime;
          numSets++;
          itsID++;
        }

      if (request.var1.f4 == &vfarnum2)
        {
          for (int levelID = 0; levelID < nlevels; ++levelID)
            {
              if (nrecs == 0) { request.var1.f5(var15[levelID], var14[levelID], request.var1.f5arg); }
              else
                {
                  const auto &array2 = var14[levelID].vec_d;
                  auto &array1 = var15[levelID].vec_d;
                  const auto len = field1.size;
                  const auto missvaltemp = field1.missval;
                  for (size_t i = 0; i < len; ++i)
                    {
                      if (DBL_IS_EQUAL(array2[i], numSets) || array2[i] > numSets) array1[i] = missvaltemp;
                    }
                }
            }
        }

      if (nrecs == 0 && numSets == 0) break;

      if (request.var1.epilog == MEAN || request.var1.epilog == PERCENT_OF_TIME)
        for (int levelID = 0; levelID < nlevels; ++levelID)
          {
            auto &rvar = IS_SET(request.var1.f5) ? var15[levelID] : var14[levelID];

            if (samp3[levelID].empty())
              fieldc_div(rvar, numSets);
            else
              field2_div(rvar, samp3[levelID]);

            if (request.var1.epilog == PERCENT_OF_TIME) fieldc_mul(rvar, 100.0);
          }
      else if (request.var1.epilog == PERCENT_OF_TOTAL_AMOUNT)
        for (int levelID = 0; levelID < nlevels; ++levelID)
          {
            Field &rvar = IS_SET(request.var1.f5) ? var15[levelID] : var14[levelID];

            field2_div(rvar, total[levelID]);
            fieldc_mul(rvar, 100.0);
          }

      if (request.var1.refdate == 19550101)
        {
          taxisDefVdatetime(otaxisID, ovDateTime);
          cdo_def_timestep(ostreamID, otsID);
        }
      else
        {
          int year, month, day;
          cdiDate_decode(inDateTime21.date, &year, &month, &day);
          define_mid_of_time(cdo_operator_f2(operatorID), otaxisID, year, month, 12);
          cdo_def_timestep(ostreamID, otsID);
        }

      if (otsID && vlistInqVarTimetype(ivlistID1, FIRST_VAR_ID) == TIME_CONSTANT) continue;

      varID = 0;
      for (int levelID = 0; levelID < nlevels; ++levelID)
        {
          Field &rvar = IS_SET(request.var1.f5) ? var15[levelID] : var14[levelID];

          vfarsel(rvar, samp1[levelID]);
          vfarsel(rvar, samp2[levelID]);

          cdo_def_record(ostreamID, varID, levelID);
          cdo_write_record(ostreamID, rvar.vec_d.data(), rvar.nmiss);
        }

      if (IS_SET(request.var2.h2))
        {
          varID = 1;
          for (int levelID = 0; levelID < nlevels; ++levelID)
            {
              auto &rvar = var22[levelID];

              vfarsel(rvar, samp1[levelID]);
              vfarsel(rvar, samp2[levelID]);

              cdo_def_record(ostreamID, varID, levelID);
              cdo_write_record(ostreamID, rvar.vec_d.data(), rvar.nmiss);
            }
        }

      if (nrecs == 0) break;

      otsID++;
    }

  cdo_stream_close(ostreamID);
  cdo_stream_close(istreamID2);
  cdo_stream_close(istreamID1);
}

void
eca3(const ECA_REQUEST_3 &request)
{
  const auto operatorID = cdo_operator_id();

  CdiDateTime ovDateTime{};
  CdiDateTime inDateTime21{};
  int varID;
  int itsID;
  int otsID;

  const auto compareDate = cdo_operator_f2(operatorID);

  const auto istreamID1 = cdo_open_read(0);
  const auto istreamID2 = cdo_open_read(1);

  const auto ivlistID1 = cdo_stream_inq_vlist(istreamID1);
  const auto ivlistID2 = cdo_stream_inq_vlist(istreamID2);
  const auto ovlistID = vlistCreate();

  vlist_compare(ivlistID1, ivlistID2, CMP_ALL);

  const auto gridID = vlistInqVarGrid(ivlistID1, FIRST_VAR_ID);
  const auto zaxisID = vlistInqVarZaxis(ivlistID1, FIRST_VAR_ID);
  const auto missval = vlistInqVarMissval(ivlistID1, FIRST_VAR_ID);

  varID = vlistDefVar(ovlistID, gridID, zaxisID, TIME_VARYING);

  vlistDefVarMissval(ovlistID, varID, missval);

  if (IS_SET(request.name)) cdiDefKeyString(ovlistID, varID, CDI_KEY_NAME, request.name);
  if (IS_SET(request.longname)) cdiDefKeyString(ovlistID, varID, CDI_KEY_LONGNAME, request.longname);
  if (IS_SET(request.units)) cdiDefKeyString(ovlistID, varID, CDI_KEY_UNITS, request.units);

  if (cdo_operator_f2(operatorID) == 16) vlistDefNtsteps(ovlistID, 1);

  const auto itaxisID1 = vlistInqTaxis(ivlistID1);
  const auto itaxisID2 = vlistInqTaxis(ivlistID2);
  const auto otaxisID = cdo_taxis_create(TAXIS_RELATIVE);
  taxisDefTunit(otaxisID, TUNIT_DAY);
  //  taxisDefTunit(otaxisID, TUNIT_MINUTE);
  //  taxisDefCalendar(otaxisID, CALENDAR_PROLEPTIC);
  taxisDefCalendar(otaxisID, taxisInqCalendar(itaxisID1));
  taxisDefRdate(otaxisID, request.refdate);
  taxisDefRtime(otaxisID, 0);
  vlistDefTaxis(ovlistID, otaxisID);

  const auto ostreamID = cdo_open_write(2);
  cdo_def_vlist(ostreamID, ovlistID);

  const auto gridsize = gridInqSize(gridID);

  Field field1, field2;
  field1.resize(gridsize);
  field2.resize(gridsize);

  const int nlevels = zaxisInqSize(zaxisID);

  FieldVector var1(nlevels), var2(nlevels);

  for (int levelID = 0; levelID < nlevels; ++levelID)
    {
      init_field(var1[levelID], gridID, missval, gridsize);
      init_field(var2[levelID], gridID, missval, gridsize);
    }

  itsID = 0;
  otsID = 0;
  while (true)
    {
      int nrecs = 0;
      long numSets = 0;
      while (true)
        {
          nrecs = cdo_stream_inq_timestep(istreamID1, itsID);
          if (nrecs == 0) break;

          if (!cdo_stream_inq_timestep(istreamID2, itsID)) cdo_abort("Input streams have different number of time steps!");

          const auto ivDateTime1 = taxisInqVdatetime(itaxisID1);
          const auto ivdate1 = cdiDate_get(ivDateTime1.date);
          const auto ivtime1 = cdiTime_get(ivDateTime1.time);

          const auto ivDateTime2 = taxisInqVdatetime(itaxisID2);
          const auto ivdate2 = cdiDate_get(ivDateTime2.date);
          const auto ivtime2 = cdiTime_get(ivDateTime2.time);

          if (ivdate1 != ivdate2) cdo_abort("Input streams have different verification dates at time step %d!", itsID + 1);
          if (ivtime1 != ivtime2) cdo_abort("Input streams have different verification times at time step %d!", itsID + 1);

          if (numSets == 0) inDateTime21 = ivDateTime1;

          if (date_is_neq(ivDateTime1, inDateTime21, compareDate))
            {
              cdo_add_steps(-1);
              break;
            }

          for (int recID = 0; recID < nrecs; ++recID)
            {
              int levelID;
              cdo_inq_record(istreamID1, &varID, &levelID);
              cdo_inq_record(istreamID2, &varID, &levelID);

              if (varID != FIRST_VAR_ID) continue;

              if (numSets == 0)
                {
                  for (size_t i = 0; i < gridsize; ++i)
                    {
                      var1[levelID].vec_d[i] = missval;
                      var2[levelID].vec_d[i] = missval;
                    }
                  var1[levelID].nmiss = gridsize;
                  var2[levelID].nmiss = gridsize;
                }

              cdo_read_record(istreamID1, field1.vec_d.data(), &field1.nmiss);
              field1.grid = var1[levelID].grid;
              field1.missval = var1[levelID].missval;

              cdo_read_record(istreamID2, field2.vec_d.data(), &field2.nmiss);
              field2.grid = var1[levelID].grid;
              field2.missval = var1[levelID].missval;

              request.f1(var1[levelID], field1);
              request.f2(var2[levelID], field2);
            }

          ovDateTime = ivDateTime1;
          numSets++;
          itsID++;
        }

      if (nrecs == 0 && numSets == 0) break;

      for (int levelID = 0; levelID < nlevels; ++levelID) request.f3(var1[levelID], var2[levelID]);

      if (request.refdate == 19550101)
        {
          taxisDefVdatetime(otaxisID, ovDateTime);
          cdo_def_timestep(ostreamID, otsID);
        }
      else
        {
          int year, month, day;
          cdiDate_decode(inDateTime21.date, &year, &month, &day);
          define_mid_of_time(cdo_operator_f2(operatorID), otaxisID, year, month, 12);
          cdo_def_timestep(ostreamID, otsID);
        }

      if (otsID && vlistInqVarTimetype(ivlistID1, FIRST_VAR_ID) == TIME_CONSTANT) continue;

      varID = 0;
      for (int levelID = 0; levelID < nlevels; ++levelID)
        {
          cdo_def_record(ostreamID, varID, levelID);
          cdo_write_record(ostreamID, var1[levelID].vec_d.data(), var1[levelID].nmiss);
        }

      if (nrecs == 0) break;
      otsID++;
    }

  cdo_stream_close(ostreamID);
  cdo_stream_close(istreamID2);
  cdo_stream_close(istreamID1);
}

// check for non missval values
static bool
fldhvs(const FieldVector &fieldVector, const size_t nlevels)
{
  for (size_t level = 0; level < nlevels; level++)
    {
      if (fieldVector[level].nmiss != fieldVector[level].size) return true;
    }

  return false;
}

void
eca4(const ECA_REQUEST_4 &request)
{
  const auto operatorID = cdo_operator_id();

  int yearcnt = 0;
  int varID;
  bool resetAtJan = false, resetAtJul = false;
  bool isFirstYear = true;
  CdiDateTime ovDateTime{};
  CdiDateTime inDateTime21{};
  int64_t ivdate = 0, ovdate = 0;

  const auto compareDate = cdo_operator_f2(operatorID);

  const auto istreamID1 = cdo_open_read(0);
  const auto istreamID2 = cdo_open_read(1);

  const auto ivlistID1 = cdo_stream_inq_vlist(istreamID1);
  const auto ivlistID2 = cdo_stream_inq_vlist(istreamID2);
  const auto ovlistID = vlistCreate();

  int gridID = vlistInqVarGrid(ivlistID1, FIRST_VAR_ID);
  if (gridInqSize(gridID) != gridInqSize(vlistInqVarGrid(ivlistID2, FIRST_VAR_ID)))
    cdo_abort("Grid sizes of the input fields do not match!");

  const auto zaxisID = vlistInqVarZaxis(ivlistID1, FIRST_VAR_ID);
  const auto missval = vlistInqVarMissval(ivlistID1, FIRST_VAR_ID);

  const auto ovarID1 = vlistDefVar(ovlistID, gridID, zaxisID, TIME_VARYING);

  vlistDefVarMissval(ovlistID, ovarID1, missval);

  if (IS_SET(request.name)) cdiDefKeyString(ovlistID, ovarID1, CDI_KEY_NAME, request.name);
  if (IS_SET(request.longname)) cdiDefKeyString(ovlistID, ovarID1, CDI_KEY_LONGNAME, request.longname);
  if (IS_SET(request.units)) cdiDefKeyString(ovlistID, ovarID1, CDI_KEY_UNITS, request.units);

  const auto ovarID2 = vlistDefVar(ovlistID, gridID, zaxisID, TIME_VARYING);

  vlistDefVarMissval(ovlistID, ovarID2, missval);

  if (IS_SET(request.name2)) cdiDefKeyString(ovlistID, ovarID2, CDI_KEY_NAME, request.name2);
  if (IS_SET(request.longname2)) cdiDefKeyString(ovlistID, ovarID2, CDI_KEY_LONGNAME, request.longname2);
  if (IS_SET(request.units2)) cdiDefKeyString(ovlistID, ovarID2, CDI_KEY_UNITS, request.units2);

  if (cdo_operator_f2(operatorID) == 16) vlistDefNtsteps(ovlistID, 1);

  const auto itaxisID = vlistInqTaxis(ivlistID1);
  const auto otaxisID = cdo_taxis_create(TAXIS_RELATIVE);
  taxisDefTunit(otaxisID, TUNIT_DAY);
  //  taxisDefTunit(otaxisID, TUNIT_MINUTE);
  //  taxisDefCalendar(otaxisID, CALENDAR_PROLEPTIC);
  taxisDefCalendar(otaxisID, taxisInqCalendar(itaxisID));
  taxisDefRdate(otaxisID, 19550101);
  taxisDefRtime(otaxisID, 0);
  vlistDefTaxis(ovlistID, otaxisID);

  const auto ostreamID = cdo_open_write(2);
  cdo_def_vlist(ostreamID, ovlistID);

  bool lyvals = true;
  const auto gridtype = gridInqType(gridID);
  if (gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN || gridtype == GRID_PROJECTION)
    {
      gridID = gridToCurvilinear(gridID, NeedCorners::Yes);
    }
  else if (gridtype == GRID_GME) { gridID = gridToUnstructured(gridID, NeedCorners::Yes); }
  else { lyvals = false; }

  const auto gridsize = gridInqSize(gridID);
  // for later check on northern\southern hemisphere
  std::vector<double> yvals(gridsize);
  if (lyvals) { gridInqYvals(gridID, yvals.data()); }
  else
    {
      for (size_t i = 0; i < gridsize; ++i) yvals[i] = 20;  // Northern hemisphere
    }

  // Two fields are needed because of the definition of gsl for northern and southern hemisphere
  Field fieldGt, fieldLt;
  fieldGt.resize(gridsize);
  fieldLt.resize(gridsize);

  // field for the land-water-distribution
  Field mask;
  mask.resize(gridsize);

  const auto nlevels = zaxisInqSize(zaxisID);

  FieldVector startCount(nlevels), endCount(nlevels);
  FieldVector gslDuration(nlevels), gslFirstDay(nlevels);

  FieldVector2D startDateWithHist(2), endDateWithHist(2);
  /* because of the different definitions for northern and southern hemisphere,
   * the values of the last year have to be present THE LAST YEAR HAS THE INDEX 1 */
  for (int h = 0; h < 2; h++)
    {
      startDateWithHist[h].resize(nlevels);
      endDateWithHist[h].resize(nlevels);
    }

  for (int levelID = 0; levelID < nlevels; ++levelID)
    {
      init_field(startCount[levelID], gridID, missval, gridsize, 0.0);
      init_field(endCount[levelID], gridID, missval, gridsize, 0.0);
      init_field(gslDuration[levelID], gridID, missval, gridsize);
      init_field(gslFirstDay[levelID], gridID, missval, gridsize);

      for (int h = 0; h < 2; h++) init_field(startDateWithHist[h][levelID], gridID, missval, gridsize);
      for (int h = 0; h < 2; h++) init_field(endDateWithHist[h][levelID], gridID, missval, gridsize);
    }

  int itsID = 0;
  int otsID = 0;

  if (cdo_stream_inq_timestep(istreamID2, itsID))
    {
      int levelID;
      cdo_inq_record(istreamID2, &varID, &levelID);
      cdo_read_record(istreamID2, mask.vec_d.data(), &mask.nmiss);
      mask.grid = gridID;
      mask.missval = vlistInqVarMissval(ivlistID2, 0);

      request.s3(mask, request.s3arg);
    }
  else
    cdo_abort("Could not read land-water mask!");

  while (true)
    {
      int nrecs = 0;
      long numSets = 0;
      while (true)
        {
          nrecs = cdo_stream_inq_timestep(istreamID1, itsID);
          if (nrecs == 0) break;

          const auto ivDateTime = taxisInqVdatetime(itaxisID);
          ivdate = cdiDate_get(ivDateTime.date);

          int month = (ivdate % 10000) / 100;
          if (month < 1 || month > 12) cdo_abort("month %d out of range!", month);

          if (numSets == 0) inDateTime21 = ivDateTime;

          if (date_is_neq(ivDateTime, inDateTime21, compareDate))
            {
              resetAtJan = false;
              resetAtJul = false;
              cdo_add_steps(-1);
              break;
            }

          for (int recID = 0; recID < nrecs; ++recID)
            {
              int levelID;
              cdo_inq_record(istreamID1, &varID, &levelID);

              if (varID != FIRST_VAR_ID) continue;

              if (numSets == 0)
                {
                  field_fill(gslDuration[levelID], missval);
                  field_fill(gslFirstDay[levelID], missval);
                  // reinitialize the current year
                  field_fill(startDateWithHist[0][levelID], missval);
                  field_fill(endDateWithHist[0][levelID], missval);

                  gslDuration[levelID].nmiss = 0;
                  gslFirstDay[levelID].nmiss = 0;
                  // reinitialize the current year
                  startDateWithHist[0][levelID].nmiss = gridsize;
                  endDateWithHist[0][levelID].nmiss = gridsize;
                }
              // init the history ONCE
              if (0 == itsID)
                {
                  field_fill(startDateWithHist[1][levelID], missval);
                  field_fill(endDateWithHist[1][levelID], missval);

                  startDateWithHist[1][levelID].nmiss = gridsize;
                  endDateWithHist[1][levelID].nmiss = gridsize;
                }

              cdo_read_record(istreamID1, fieldGt.vec_d.data(), &fieldGt.nmiss);
              fieldLt.vec_d = fieldGt.vec_d;
              fieldLt.nmiss = fieldGt.nmiss;
              fieldGt.grid = startCount[levelID].grid;
              fieldGt.missval = startCount[levelID].missval;
              fieldLt.grid = startCount[levelID].grid;
              fieldLt.missval = startCount[levelID].missval;

              // Reinitialization of (start|end)Count variables has to be done different for norther and southern hemisphere
              if (1 == month && !resetAtJan)
                {
                  // reset northern startCount
                  for (size_t i = 0; i < gridsize; ++i)
                    {
                      if (yvals[i] >= 0.0)
                        if (!DBL_IS_EQUAL(startCount[levelID].vec_d[i], missval))
                          {
                            startCount[levelID].vec_d[i] = missval;
                            startCount[levelID].nmiss++;
                          }
                    }
                  // reset southern endCount
                  for (size_t i = 0; i < gridsize; ++i)
                    {
                      if (yvals[i] < 0.0)
                        if (!DBL_IS_EQUAL(endCount[levelID].vec_d[i], missval))
                          {
                            endCount[levelID].vec_d[i] = missval;
                            endCount[levelID].nmiss++;
                          }
                    }

                  resetAtJan = true;
                }
              if (7 == month && !resetAtJul)
                {
#ifdef _OPENMP
#pragma omp sections
#endif
                  {
#ifdef _OPENMP
#pragma omp section
#endif
                    {
                      // reset northern endCount
                      for (size_t i = 0; i < gridsize; ++i)
                        {
                          if (yvals[i] >= 0.0)
                            {
                              if (!DBL_IS_EQUAL(endCount[levelID].vec_d[i], missval))
                                {
                                  endCount[levelID].vec_d[i] = missval;
                                  endCount[levelID].nmiss++;
                                }
                            }
                        }
                    }
#ifdef _OPENMP
#pragma omp section
#endif
                    {
                      // reset southern startCount
                      for (size_t i = 0; i < gridsize; ++i)
                        {
                          if (yvals[i] < 0.0)
                            {
                              if (!DBL_IS_EQUAL(startCount[levelID].vec_d[i], missval))
                                {
                                  startCount[levelID].vec_d[i] = missval;
                                  startCount[levelID].nmiss++;
                                }
                            }
                        }
                    }
                  }
                  resetAtJul = true;
                }

                // count the day with temperature larger/smaller than the given limit
#ifdef _OPENMP
#pragma omp sections
#endif
              {
#ifdef _OPENMP
#pragma omp section
#endif
                {
                  vfarsel(fieldGt, mask);
                  request.s1(fieldGt, request.s1arg);
                  vfarnum2(startCount[levelID], fieldGt);
                }
#ifdef _OPENMP
#pragma omp section
#endif
                {
                  vfarsel(fieldLt, mask);
                  request.s2(fieldLt, request.s1arg);
                  vfarnum2(endCount[levelID], fieldLt);
                }
              }

              if (month < 7)
                {
                  for (size_t i = 0; i < gridsize; ++i)  // dictinct between northern and southern sphere
                    // start with south
                    if (yvals[i] < 0)
                      {
                        // south: periods can also start in the first half of the year, but this date has already gone into the
                        // history
                        if (DBL_IS_EQUAL(startDateWithHist[1][levelID].vec_d[i], missval)
                            && IS_EQUAL(startCount[levelID].vec_d[i], request.consecutiveDays))
                          {
                            startDateWithHist[1][levelID].vec_d[i] = ivdate;
                            // reset the endCount, because we are only interessted in the end of the eriod, if a start was found
                            endCount[levelID].vec_d[i] = missval;
                            endDateWithHist[0][levelID].vec_d[i] = missval;
                          }
                        if (DBL_IS_EQUAL(endDateWithHist[0][levelID].vec_d[i], missval)
                            && IS_EQUAL(endCount[levelID].vec_d[i], request.consecutiveDays))
                          {
                            endDateWithHist[0][levelID].vec_d[i] = ivdate;
                          }
                      }
                    else
                      {
                        if (DBL_IS_EQUAL(startDateWithHist[0][levelID].vec_d[i], missval)
                            && IS_EQUAL(startCount[levelID].vec_d[i], request.consecutiveDays))
                          {
                            startDateWithHist[0][levelID].vec_d[i] = ivdate;
                          }
                      }
                }
              else
                {
                  for (size_t i = 0; i < gridsize; ++i)
                    {
                      if (yvals[i] < 0)
                        {
                          if (DBL_IS_EQUAL(startDateWithHist[0][levelID].vec_d[i], missval)
                              && IS_EQUAL(startCount[levelID].vec_d[i], request.consecutiveDays))
                            {
                              startDateWithHist[0][levelID].vec_d[i] = ivdate;
                            }
                        }
                      else
                        {
                          // north: periods can also start in the second half of the year
                          if (DBL_IS_EQUAL(startDateWithHist[0][levelID].vec_d[i], missval)
                              && IS_EQUAL(startCount[levelID].vec_d[i], request.consecutiveDays))
                            {
                              startDateWithHist[0][levelID].vec_d[i] = ivdate;
                              // reset the endCount, because we are only interessted in the end of the eriod, if a start was found
                              endCount[levelID].vec_d[i] = missval;
                              endDateWithHist[0][levelID].vec_d[i] = missval;
                            }
                          if (DBL_IS_EQUAL(endDateWithHist[0][levelID].vec_d[i], missval)
                              && IS_EQUAL(endCount[levelID].vec_d[i], request.consecutiveDays))
                            {
                              endDateWithHist[0][levelID].vec_d[i] = ivdate;
                            }
                        }
                    }
                }
              // update nmiss for saving data in GRIB
              startCount[levelID].nmiss = field_num_miss(startCount[levelID]);
              endCount[levelID].nmiss = field_num_miss(endCount[levelID]);
              startDateWithHist[1][levelID].nmiss = field_num_miss(startDateWithHist[1][levelID]);
              startDateWithHist[0][levelID].nmiss = field_num_miss(startDateWithHist[0][levelID]);
              endDateWithHist[1][levelID].nmiss = field_num_miss(endDateWithHist[1][levelID]);
              endDateWithHist[0][levelID].nmiss = field_num_miss(endDateWithHist[0][levelID]);
            }

          ovDateTime = ivDateTime;
          ovdate = cdiDate_get(ovDateTime.date);
          numSets++;
          itsID++;
        }

      if (nrecs == 0 && numSets == 0) break;

      adjust_end_date(nlevels, gridsize, yvals, missval, ovdate, startDateWithHist, endDateWithHist);

      /*  compute and write GSL for the previous year
       *  AND
       *  write the current start/end dates into the history
       *
       *  this is the default action if more than a year is available */
      if (yearcnt != 0)
        {
          compute_gsl(nlevels, gridsize, yvals, missval, startDateWithHist, endDateWithHist, gslDuration, gslFirstDay, false);

          // values of the privous year
          ovDateTime.date = cdiDate_encode(ovdate / 10000 - 1, 12, 31);
          write_gsl_stream(ostreamID, otaxisID, otsID, ovarID1, ovarID2, ivlistID1, FIRST_VAR_ID, gslDuration, gslFirstDay,
                           ovDateTime, nlevels);
          otsID++;
        }

      // if there is a previous year
      if (ovdate != ivdate)
        {
          /*  if the first year of data was processed, the history has to
           *  be checked befor it get's updated. This is necessary, if a
           *  growing period on the southern hemisphere was found. Otherwise,
           *  it would get overwritten. */
          if (isFirstYear)
            {
              // Check for non missing values, i.e. is there any data for the previous year?
              if (fldhvs(startDateWithHist[1], nlevels))
                {
                  compute_gsl(nlevels, gridsize, yvals, missval, startDateWithHist, endDateWithHist, gslDuration, gslFirstDay,
                              false);
                  ovDateTime.date = cdiDate_encode(ovdate / 10000 - 1, 12, 31);
                  write_gsl_stream(ostreamID, otaxisID, otsID, ovarID1, ovarID2, ivlistID1, FIRST_VAR_ID, gslDuration, gslFirstDay,
                                   ovDateTime, nlevels);
                  otsID++;
                }
              isFirstYear = false;
            }
#ifdef _OPENMP
#pragma omp sections
#endif
          {
            update_hist(startDateWithHist, nlevels, gridsize, yvals, false);
#ifdef _OPENMP
#pragma omp section
#endif
            update_hist(endDateWithHist, nlevels, gridsize, yvals, true);
          }
        }
      else  // process the current year, this only happens, if the last timestep is reached OR if data for only one year is present
        {
          compute_gsl(nlevels, gridsize, yvals, missval, startDateWithHist, endDateWithHist, gslDuration, gslFirstDay, true);

          write_gsl_stream(ostreamID, otaxisID, otsID, ovarID1, ovarID2, ivlistID1, FIRST_VAR_ID, gslDuration, gslFirstDay,
                           ovDateTime, nlevels);
          otsID++;
        }
      yearcnt++;

      if (nrecs == 0) break;
    }

  cdo_stream_close(ostreamID);
  cdo_stream_close(istreamID2);
  cdo_stream_close(istreamID1);
}
