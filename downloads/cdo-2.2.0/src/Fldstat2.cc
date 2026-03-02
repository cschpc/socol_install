/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Fldstat2    fldcor         Correlation in grid space
      Fldstat2    fldcovar       Covariance in grid space
*/

#include <cdi.h>

#include "process_int.h"
#include "cdo_vlist.h"
#include <mpim_grid.h>
#include "field_functions.h"

// routine corr copied from PINGO
// correclation in space
auto correlation_kernel = [](const auto v1, const auto mv1, const auto v2, const auto mv2, const auto w, auto &sum0, auto &sum1,
                             auto &sum00, auto &sum01, auto &sum11, auto &wsum0, auto isEQ) {
  if (!isEQ(w, mv1) && !isEQ(v1, mv1) && !isEQ(v2, mv2))
    {
      sum0 += w * v1;
      sum1 += w * v2;
      sum00 += w * v1 * v1;
      sum01 += w * v1 * v2;
      sum11 += w * v2 * v2;
      wsum0 += w;
    }
};

template <typename T1, typename T2>
static double
correlation(const Varray<T1> &v1, const Varray<T2> &v2, const Varray<double> &weight, double missval1, double missval2,
            size_t gridsize)
{
  double sum0 = 0.0, sum1 = 0.0, sum00 = 0.0, sum01 = 0.0, sum11 = 0.0, wsum0 = 0.0;

  if (std::isnan(missval1) || std::isnan(missval2))
    {
      for (size_t i = 0; i < gridsize; ++i)
        correlation_kernel(v1[i], missval1, v2[i], missval2, weight[i], sum0, sum1, sum00, sum01, sum11, wsum0, dbl_is_equal);
    }
  else
    {
      for (size_t i = 0; i < gridsize; ++i)
        correlation_kernel(v1[i], missval1, v2[i], missval2, weight[i], sum0, sum1, sum00, sum01, sum11, wsum0, is_equal);
    }

  const auto out = is_not_equal(wsum0, 0.0)
                       ? DIVMN((sum01 * wsum0 - sum0 * sum1), SQRTMN((sum00 * wsum0 - sum0 * sum0) * (sum11 * wsum0 - sum1 * sum1)))
                       : missval1;

  return out;
}

static double
correlation(const Field &field1, const Field &field2, const Varray<double> &weight)
{
  if (field1.memType == MemType::Float && field2.memType == MemType::Float)
    return correlation(field1.vec_f, field2.vec_f, weight, field1.missval, field2.missval, field1.size);
  else if (field1.memType == MemType::Float && field2.memType == MemType::Double)
    return correlation(field1.vec_f, field2.vec_d, weight, field1.missval, field2.missval, field1.size);
  else if (field1.memType == MemType::Double && field2.memType == MemType::Float)
    return correlation(field1.vec_d, field2.vec_f, weight, field1.missval, field2.missval, field1.size);
  else
    return correlation(field1.vec_d, field2.vec_d, weight, field1.missval, field2.missval, field1.size);
}

// covariance in space
auto covariance_kernel = [](const auto v1, const auto mv1, const auto v2, const auto mv2, const auto w, auto &sum0, auto &sum1,
                            auto &sum01, auto &wsum0, auto isEQ) {
  if (!isEQ(w, mv1) && !isEQ(v1, mv1) && !isEQ(v2, mv2))
    {
      sum0 += w * v1;
      sum1 += w * v2;
      sum01 += w * v1 * v2;
      wsum0 += w;
    }
};

template <typename T1, typename T2>
static double
covariance(const Varray<T1> &v1, const Varray<T2> &v2, const Varray<double> &weight, double missval1, double missval2,
           size_t gridsize)
{
  double sum0 = 0.0, sum1 = 0.0, sum01 = 0.0, wsum0 = 0.0;

  if (std::isnan(missval1) || std::isnan(missval2))
    {
      for (size_t i = 0; i < gridsize; ++i)
        covariance_kernel(v1[i], missval1, v2[i], missval2, weight[i], sum0, sum1, sum01, wsum0, dbl_is_equal);
    }
  else
    {
      for (size_t i = 0; i < gridsize; ++i)
        covariance_kernel(v1[i], missval1, v2[i], missval2, weight[i], sum0, sum1, sum01, wsum0, is_equal);
    }

  const auto out = is_not_equal(wsum0, 0.0) ? (sum01 * wsum0 - sum0 * sum1) / (wsum0 * wsum0) : missval1;

  return out;
}

static double
covariance(const Field &field1, const Field &field2, const Varray<double> &weight)
{
  if (field1.memType == MemType::Float && field2.memType == MemType::Float)
    return covariance(field1.vec_f, field2.vec_f, weight, field1.missval, field2.missval, field1.size);
  else if (field1.memType == MemType::Float && field2.memType == MemType::Double)
    return covariance(field1.vec_f, field2.vec_d, weight, field1.missval, field2.missval, field1.size);
  else if (field1.memType == MemType::Double && field2.memType == MemType::Float)
    return covariance(field1.vec_d, field2.vec_f, weight, field1.missval, field2.missval, field1.size);
  else
    return covariance(field1.vec_d, field2.vec_d, weight, field1.missval, field2.missval, field1.size);
}

void *
Fldstat2(void *process)
{
  auto wstatus = false;
  auto needWeights = true;

  cdo_initialize(process);

  // clang-format off
  cdo_operator_add("fldcor",   FieldFunc_Cor,   0, nullptr);
  cdo_operator_add("fldcovar", FieldFunc_Covar, 0, nullptr);
  // clang-format on

  const auto operatorID = cdo_operator_id();
  const auto operfunc = cdo_operator_f1(operatorID);

  const auto streamID1 = cdo_open_read(0);
  const auto streamID2 = cdo_open_read(1);

  const auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  const auto vlistID2 = cdo_stream_inq_vlist(streamID2);
  const auto vlistID3 = vlistDuplicate(vlistID1);

  vlist_compare(vlistID1, vlistID2, CMP_ALL);

  const auto taxisID1 = vlistInqTaxis(vlistID1);
  const auto taxisID3 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID3, taxisID3);

  VarList varList1, varList2;
  varListInit(varList1, vlistID1);
  varListInit(varList2, vlistID2);

  double slon = 0.0, slat = 0.0;
  const auto gridID3 = gridCreate(GRID_LONLAT, 1);
  gridDefXsize(gridID3, 1);
  gridDefYsize(gridID3, 1);
  gridDefXvals(gridID3, &slon);
  gridDefYvals(gridID3, &slat);

  const auto ngrids = vlistNgrids(vlistID1);

  for (int index = 0; index < ngrids; ++index) vlistChangeGridIndex(vlistID3, index, gridID3);

  Field field1, field2;

  const auto streamID3 = cdo_open_write(2);
  cdo_def_vlist(streamID3, vlistID3);

  const auto gridsizemax = vlistGridsizeMax(vlistID1);
  Varray<double> weight;
  if (needWeights) weight.resize(gridsizemax);

  int lastgridID = -1;
  int tsID = 0;
  while (true)
    {
      const auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
      if (nrecs == 0) break;

      const auto nrecs2 = cdo_stream_inq_timestep(streamID2, tsID);
      if (nrecs2 == 0)
        {
          cdo_warning("Input streams have different number of time steps!");
          break;
        }

      cdo_taxis_copy_timestep(taxisID3, taxisID1);
      cdo_def_timestep(streamID3, tsID);

      for (int recID = 0; recID < nrecs; ++recID)
        {
          int varID, levelID;
          cdo_inq_record(streamID1, &varID, &levelID);
          field1.init(varList1[varID]);
          cdo_inq_record(streamID2, &varID, &levelID);
          field2.init(varList2[varID]);
          cdo_read_record(streamID1, field1);
          cdo_read_record(streamID2, field2);

          const auto gridID = varList1[varID].gridID;
          if (needWeights && gridID != lastgridID)
            {
              lastgridID = gridID;
              wstatus = (gridcell_weights(gridID, weight) != 0);
            }
          if (wstatus && tsID == 0 && levelID == 0)
            cdo_warning("Using constant grid cell area weights for variable %s!", varList1[varID].name);

          double sglval = 0.0;
          if (operfunc == FieldFunc_Cor)
            sglval = correlation(field1, field2, weight);
          else if (operfunc == FieldFunc_Covar)
            sglval = covariance(field1, field2, weight);

          const auto nmiss3 = DBL_IS_EQUAL(sglval, varList1[varID].missval) ? 1 : 0;

          cdo_def_record(streamID3, varID, levelID);
          cdo_write_record(streamID3, &sglval, nmiss3);
        }

      tsID++;
    }

  cdo_stream_close(streamID3);
  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  cdo_finish();

  return nullptr;
}
