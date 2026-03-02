/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Ymoncomp   ymoneq         Equal
      Ymoncomp   ymonne         Not equal
      Ymoncomp   ymonle         Less equal
      Ymoncomp   ymonlt         Less than
      Ymoncomp   ymonge         Greater equal
      Ymoncomp   ymongt         Greater than
      Yseascomp  yseaseq        Equal
      Yseascomp  yseasne        Not equal
      Yseascomp  yseasle        Less equal
      Yseascomp  yseaslt        Less than
      Yseascomp  yseasge        Greater equal
      Yseascomp  yseasgt        Greater than
*/

#include <cdi.h>

#include "cdo_vlist.h"
#include "cdo_season.h"
#include "process_int.h"
#include "field_functions.h"


static auto func_comp
    = [](auto hasMissvals, auto n, auto mv1, auto mv2, auto &v1, const auto &v2, auto binray_operator) {
        if (hasMissvals)
          {
            if (std::isnan(mv1) || std::isnan(mv2))
              for (size_t i = 0; i < n; ++i)
                v1[i] = (dbl_is_equal(v1[i], mv1) || dbl_is_equal(v2[i], mv2)) ? mv1 : binray_operator(v1[i], v2[i]);
            else
              for (size_t i = 0; i < n; ++i)
                v1[i] = (is_equal(v1[i], mv1) || is_equal(v2[i], mv2)) ? mv1 : binray_operator(v1[i], v2[i]);
          }
        else
          {
            for (size_t i = 0; i < n; ++i) v1[i] = binray_operator(v1[i], v2[i]);
          }
      };


template <typename T1, typename T2>
void comp_function(int operFunc, bool hasMissvals, size_t ngp, T1 missval1, T2 missval2, Varray<T1> &v1, const Varray<T2> &v2)
{
  // clang-format off
  if      (operFunc == FieldFunc_EQ) func_comp(hasMissvals, ngp, missval1, missval2, v1, v2, binary_op_EQ);
  else if (operFunc == FieldFunc_NE) func_comp(hasMissvals, ngp, missval1, missval2, v1, v2, binary_op_NE);
  else if (operFunc == FieldFunc_LE) func_comp(hasMissvals, ngp, missval1, missval2, v1, v2, binary_op_LE);
  else if (operFunc == FieldFunc_LT) func_comp(hasMissvals, ngp, missval1, missval2, v1, v2, binary_op_LT);
  else if (operFunc == FieldFunc_GE) func_comp(hasMissvals, ngp, missval1, missval2, v1, v2, binary_op_GE);
  else if (operFunc == FieldFunc_GT) func_comp(hasMissvals, ngp, missval1, missval2, v1, v2, binary_op_GT);
  else cdo_abort("Operator not implemented!");
  // clang-format on
}

void comp_function(Field &f1, const Field &f2, int operFunc)
{
  auto hasMissvals = (f1.nmiss > 0 || f2.nmiss > 0);
  if (f1.memType == MemType::Float && f2.memType == MemType::Float)
    comp_function(operFunc, hasMissvals, f1.size, (float)f1.missval, (float)f2.missval, f1.vec_f, f2.vec_f);
  else if (f1.memType == MemType::Float && f2.memType == MemType::Double)
    comp_function(operFunc, hasMissvals, f1.size, (float)f1.missval, f2.missval, f1.vec_f, f2.vec_d);
  else if (f1.memType == MemType::Double && f2.memType == MemType::Float)
    comp_function(operFunc, hasMissvals, f1.size, f1.missval, (float)f2.missval, f1.vec_d, f2.vec_f);
  else
    comp_function(operFunc, hasMissvals, f1.size, f1.missval, f2.missval, f1.vec_d, f2.vec_d);
}

constexpr int MaxMonths = 12;

static const char *monthNames[]
    = { "January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December" };

static int
get_month_index(const CdiDate vDate)
{
  int year, mon, day;
  cdiDate_decode(vDate, &year, &mon, &day);
  if (mon < 1 || mon > MaxMonths) cdo_abort("Month %d out of range!", mon);
  mon--;
  return mon;
}

enum
{
  MONTHLY,
  SEASONAL
};

static void
add_operators(void)
{
  // clang-format off
  cdo_operator_add("ymoneq",  FieldFunc_EQ, MONTHLY, nullptr);
  cdo_operator_add("ymonne",  FieldFunc_NE, MONTHLY, nullptr);
  cdo_operator_add("ymonle",  FieldFunc_LE, MONTHLY, nullptr);
  cdo_operator_add("ymonlt",  FieldFunc_LT, MONTHLY, nullptr);
  cdo_operator_add("ymonge",  FieldFunc_GE, MONTHLY, nullptr);
  cdo_operator_add("ymongt",  FieldFunc_GT, MONTHLY, nullptr);
  cdo_operator_add("yseaseq", FieldFunc_EQ, SEASONAL, nullptr);
  cdo_operator_add("yseasne", FieldFunc_NE, SEASONAL, nullptr);
  cdo_operator_add("yseasle", FieldFunc_LE, SEASONAL, nullptr);
  cdo_operator_add("yseaslt", FieldFunc_LT, SEASONAL, nullptr);
  cdo_operator_add("yseasge", FieldFunc_GE, SEASONAL, nullptr);
  cdo_operator_add("yseasgt", FieldFunc_GT, SEASONAL, nullptr);
  // clang-format on
}

class ModuleYmonComp
{
  int operfunc;
  int opertype;

  CdoStreamID streamID1;
  CdoStreamID streamID2;
  CdoStreamID streamID3;

  int taxisID1;
  int taxisID2;
  int taxisID3;

  int vlistID2;

  VarList varList1;
  Field field;
  FieldVector2D vars2[MaxMonths];

  void
  process_data(int tsID, int nrecs, int mon)
  {
    cdo_taxis_copy_timestep(taxisID3, taxisID1);
    cdo_def_timestep(streamID3, tsID);

    for (int recID = 0; recID < nrecs; ++recID)
      {
        int varID, levelID;
        cdo_inq_record(streamID1, &varID, &levelID);
        field.init(varList1[varID]);
        cdo_read_record(streamID1, field);

        comp_function(field, vars2[mon][varID][levelID], operfunc);

        cdo_def_record(streamID3, varID, levelID);
        cdo_write_record(streamID3, field);
      }
  }

  void
  load_data_from_stream(int mon, int nrecs)
  {
    fields_from_vlist(vlistID2, vars2[mon], FIELD_VEC | FIELD_NAT);

    for (int recID = 0; recID < nrecs; ++recID)
      {
        int varID, levelID;
        cdo_inq_record(streamID2, &varID, &levelID);
        cdo_read_record(streamID2, vars2[mon][varID][levelID]);
      }
  }

  void
  run_monthly()
  {
    for (int nrecs, tsID = 0; (nrecs = cdo_stream_inq_timestep(streamID2, tsID) != 0); tsID++)
      {
        auto mon = get_month_index(taxisInqVdatetime(taxisID2).date);
        if (vars2[mon].size())
          cdo_abort("%s already allocated! The second input file must contain monthly mean values for a maximum of one year.",
                    monthNames[mon]);

        load_data_from_stream(mon, nrecs);
      }

    for (int nrecs, tsID = 0; (nrecs = cdo_stream_inq_timestep(streamID1, tsID) != 0); tsID++)
      {
        auto mon = get_month_index(taxisInqVdatetime(taxisID1).date);
        if (vars2[mon].size() == 0)
          cdo_abort("%s not found! The second input file must contain monthly mean values for a maximum of one year.",
                    monthNames[mon]);

        process_data(tsID, nrecs, mon);
      }
  }

  void
  run_seasonal()
  {
    const char **seasonNames = get_season_name();
    for (int nrecs, tsID = 0; (nrecs = cdo_stream_inq_timestep(streamID2, tsID) != 0); tsID++)
      {
        auto mon = get_month_index(taxisInqVdatetime(taxisID2).date);
        auto season = month_to_season(mon + 1);
        if (vars2[season].size()) cdo_abort("Season %s already allocated!", seasonNames[season]);

        load_data_from_stream(mon, nrecs);
      }

    for (int nrecs, tsID = 0; (nrecs = cdo_stream_inq_timestep(streamID1, tsID) != 0); tsID++)
      {
        auto mon = get_month_index(taxisInqVdatetime(taxisID1).date);
        auto season = month_to_season(mon + 1);
        if (vars2[season].size() == 0) cdo_abort("Season %s not found!", seasonNames[season]);

        process_data(tsID, nrecs, season);
      }
  }

public:
  void
  init(void *process)
  {
    cdo_initialize(process);

    add_operators();

    auto operatorID = cdo_operator_id();
    opertype = cdo_operator_f2(operatorID);
    operfunc = cdo_operator_f1(operatorID);

    operator_check_argc(0);

    streamID1 = cdo_open_read(0);
    streamID2 = cdo_open_read(1);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    vlistID2 = cdo_stream_inq_vlist(streamID2);
    auto vlistID3 = vlistDuplicate(vlistID1);

    vlist_compare(vlistID1, vlistID2, CMP_ALL);

    varListInit(varList1, vlistID1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = vlistInqTaxis(vlistID2);
    taxisID3 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID3, taxisID3);

    streamID3 = cdo_open_write(2);
    cdo_def_vlist(streamID3, vlistID3);
  }

  void
  run()
  {
    if (opertype == SEASONAL)
      run_seasonal();
    else
      run_monthly();
  }

  void
  close()
  {
    cdo_stream_close(streamID2);
    cdo_stream_close(streamID3);
    cdo_stream_close(streamID1);
    cdo_finish();
  }
};

void *
Ymoncomp(void *process)
{
  ModuleYmonComp ymoncomp;
  ymoncomp.init(process);
  ymoncomp.run();
  ymoncomp.close();

  return nullptr;
}
