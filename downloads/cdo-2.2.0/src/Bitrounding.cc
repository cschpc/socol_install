/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include "process_int.h"
#include "cdo_vlist.h"
#include "cdo_options.h"
#include "pmlist.h"
#include "param_conversion.h"
#include "progress.h"
#include "bitinformation.h"

struct BitroundParams
{
  double infLevel = 0.9999;
  int minBits = 1;
  int maxBits = 23;
  int addBits = 0;
  int numBits = -1;
  int numSteps = -1;
  std::string filename;
  bool printBits = false;
};

struct VarStat
{
  int nsbMin = 1000;
  int nsbMax = -1000;
};

static int
get_keepbits(const MutualInformation &bitInfo, double inflevel)
{
  // xbitinfo::get_keepbits v0.0.1 (https://github.com/observingClouds/xbitinfo)
  // Converted from Python to C++ by Uwe Schulzweida

  constexpr int floatNMBITS = 9;  // number of non mantissa bits for float
  int keepMantissaBits = 23;

  double bitInfoMax = -9.e33;
  for (int i = 0; i < NBITS; ++i) bitInfoMax = std::max(bitInfoMax, bitInfo.M[i]);
  // printf("bitInfoMax %g\n", bitInfoMax);

  double bitInfoMaxLast4 = -9.e33;
  for (int i = NBITS - 4; i < NBITS; ++i) bitInfoMaxLast4 = std::max(bitInfoMaxLast4, bitInfo.M[i]);
  // printf("bitInfoMax/bitInfoMaxLast4 %g\n", bitInfoMax/bitInfoMaxLast4);
  bitInfoMaxLast4 *= 1.5;
  // printf("bitInfoMaxLast4 %g\n", bitInfoMaxLast4);

  MutualInformation infoPerBitCleaned;
  for (int i = 0; i < NBITS; ++i) infoPerBitCleaned.M[i] = (bitInfo.M[i] > bitInfoMaxLast4) ? bitInfo.M[i] : 0.0;
  // for (int i = 0; i < NBITS; ++i) printf("cleaned[%d] %g\n", i + 1, infoPerBitCleaned.M[i]);

  for (int i = 1; i < NBITS; ++i) infoPerBitCleaned.M[i] += infoPerBitCleaned.M[i - 1];
  // for (int i = 0; i < NBITS; ++i) printf("cumsum[%d] %g\n", i + 1, infoPerBitCleaned.M[i]);

  auto lastBit = infoPerBitCleaned.M[NBITS - 1];
  if (lastBit > 0.0)
    {
      MutualInformation cdf;
      for (int i = 0; i < NBITS; ++i) cdf.M[i] = infoPerBitCleaned.M[i] / lastBit;
      // for (int i = 0; i < NBITS; ++i) printf("cdf[%d] %g\n", i + 1, infoPerBitCleaned.M[i]);

      constexpr int nonMantissaBits = floatNMBITS;

      for (int i = 0; i < NBITS; ++i)
        if (cdf.M[i] > inflevel)
          {
            keepMantissaBits = i + 1 - nonMantissaBits;
            break;
          }
    }

  // printf("keepMantissaBits: %d\n", keepMantissaBits);

  int nsb = std::min(std::max(keepMantissaBits, 1), 23);

  return nsb;
}

static int
bit_rounding(size_t len, Varray<float> v, double infLevel)  // copy v!
{
  signed_exponent(v.data(), len);

  auto bitInfo = bitinformation(v.data(), len);
  // if (Options::cdoVerbose) for (int i = 0; i < NBITS; ++i) fprintf(stderr, "bitInfo[%d] %.8e %g\n", i+1, bitInfo.M[i],
  // bitInfo.M[i]);

  return get_keepbits(bitInfo, infLevel);
}

static void
bitround(int nsb, size_t len, Varray<float> &v, float missval)
{
  // BitRound from NetCDF 4.9.0; routine nv4var.c

  constexpr uint32_t BIT_XPL_NBR_SGN_FLT = 23;

  // BitRound interprets nsd as number of significant binary digits (bits)
  uint32_t prc_bnr_xpl_rqr = nsb;

  uint32_t bit_xpl_nbr_zro = BIT_XPL_NBR_SGN_FLT - prc_bnr_xpl_rqr;

  // Create mask
  uint32_t msk_f32_u32_zro = 0u;       // Zero all bits
  msk_f32_u32_zro = ~msk_f32_u32_zro;  // Turn all bits to ones

  // BitShave mask for AND: Left shift zeros into bits to be rounded, leave ones in untouched bits.
  msk_f32_u32_zro <<= bit_xpl_nbr_zro;

  // BitSet mask for OR: Put ones into bits to be set, zeros in untouched bits.
  uint32_t msk_f32_u32_one = ~msk_f32_u32_zro;

  // BitRound mask for ADD: Set one bit: the MSB of LSBs
  uint32_t msk_f32_u32_hshv = msk_f32_u32_one & (msk_f32_u32_zro >> 1);

  // BitRound: Quantize to user-specified NSB with IEEE-rounding
  uint32_t *u32_ptr = (uint32_t *) v.data();

  for (size_t idx = 0; idx < len; idx++)
    {
      if (is_not_equal(v[idx], missval))
        {
          u32_ptr[idx] += msk_f32_u32_hshv;  // Add 1 to the MSB of LSBs, carry 1 to mantissa or even exponent
          u32_ptr[idx] &= msk_f32_u32_zro;   // Shave it
        }
    }
}

static void
check_range(double value, double minVal, double maxVal, const std::string &key)
{
  if (value < minVal || value > maxVal) cdo_abort("Parameter %s=%g out of range (min=%g/max=%g)!", key, value, minVal, maxVal);
}

static void
check_range(int value, int minVal, int maxVal, const std::string &key)
{
  if (value < minVal || value > maxVal) cdo_abort("Parameter %s=%d out of range (min=%d/max=%d)!", key, value, minVal, maxVal);
}

static BitroundParams
get_parameter()
{
  BitroundParams params;

  auto pargc = cdo_operator_argc();
  if (pargc)
    {
      const auto &pargv = cdo_get_oper_argv();

      KVList kvlist;
      kvlist.name = cdo_module_name();
      if (kvlist.parse_arguments(pargc, pargv) != 0) cdo_abort("Parse error!");
      if (Options::cdoVerbose) kvlist.print();

      for (const auto &kv : kvlist)
        {
          const auto &key = kv.key;
          if (kv.nvalues > 1) cdo_abort("Too many values for parameter key >%s<!", key);
          if (kv.nvalues < 1) cdo_abort("Missing value for parameter key >%s<!", key);
          const auto &value = kv.values[0];

          // clang-format off
          if      (key == "inflevel")  check_range(params.infLevel = parameter_to_double(value), 0.0, 1.0, key);
          else if (key == "minbits")   check_range(params.minBits = parameter_to_int(value), 1, 23, key);
          else if (key == "maxbits")   check_range(params.maxBits = parameter_to_int(value), 1, 23, key);
          else if (key == "addbits")   check_range(params.addBits = parameter_to_int(value), 0, 22, key);
          else if (key == "numbits")   check_range(params.numBits = parameter_to_int(value), 1, 23, key);
          else if (key == "numsteps")  check_range(params.numSteps = parameter_to_int(value), 1, 1, key);
          else if (key == "printbits") params.printBits = parameter_to_bool(value);
          else if (key == "filename")  params.filename = parameter_to_word(value);
          else cdo_abort("Invalid parameter key >%s<!", key);
          // clang-format on
        }
    }

  return params;
}

static void
print_parameter(const BitroundParams &params)
{
  std::stringstream outbuffer;

  outbuffer << "inflevel=" << params.infLevel;
  outbuffer << ", minbits=" << params.minBits;
  outbuffer << ", maxbits=" << params.maxBits;
  outbuffer << ", addbits=" << params.addBits;
  outbuffer << ", numbits=" << params.numBits;
  outbuffer << ", numsteps=" << params.numSteps;
  outbuffer << ", printbits=" << params.printBits;
  outbuffer << ", filename=" << params.filename;

  cdo_verbose("%s", outbuffer.str());
}

static void
check_attributes(int vlistID)
{
  int numBits = -1;
  auto status1 = cdiInqAttInt(vlistID, CDI_GLOBAL, "cdo_bitrounding_numbits", 1, &numBits);
  double infLevel = -1.0;
  auto status2 = cdiInqAttFlt(vlistID, CDI_GLOBAL, "cdo_bitrounding_inflevel", 1, &infLevel);
  char filename[2];
  auto status3 = cdiInqAttTxt(vlistID, CDI_GLOBAL, "cdo_bitrounding_filename", 1, filename);

  if ((status1 == 0 && numBits != -1) || (status2 == 0 && infLevel > 0.0) || status3 == 0)
    cdo_warning("It looks like CDO bitrounding has been applied to the input data before!");
}

static void
set_local_attributes(int vlistID, int varID, int numBits)
{
  cdiDefAttInt(vlistID, varID, "_QuantizeBitRoundNumberOfSignificantBits", CDI_DATATYPE_INT32, 1, &numBits);
}

static void
set_global_attributes(int vlistID, const BitroundParams &params, int numVarsHaveNumbits)
{
  if (params.filename.size() && numVarsHaveNumbits > 0)
    cdiDefAttTxt(vlistID, CDI_GLOBAL, "cdo_bitrounding_filename", (int) params.filename.size(), params.filename.c_str());

  if (numVarsHaveNumbits == vlistNvars(vlistID)) return;

  if (params.numBits != -1)
    {
      cdiDefAttInt(vlistID, CDI_GLOBAL, "cdo_bitrounding_numbits", CDI_DATATYPE_INT32, 1, &params.numBits);
    }
  else
    {
      cdiDefAttFlt(vlistID, CDI_GLOBAL, "cdo_bitrounding_inflevel", CDI_DATATYPE_FLT64, 1, &params.infLevel);
      if (params.addBits) cdiDefAttInt(vlistID, CDI_GLOBAL, "cdo_bitrounding_addbits", CDI_DATATYPE_INT32, 1, &params.addBits);
      if (params.minBits > 1) cdiDefAttInt(vlistID, CDI_GLOBAL, "cdo_bitrounding_minbits", CDI_DATATYPE_INT32, 1, &params.minBits);
      if (params.maxBits < 23) cdiDefAttInt(vlistID, CDI_GLOBAL, "cdo_bitrounding_maxbits", CDI_DATATYPE_INT32, 1, &params.maxBits);
      if (params.numSteps != -1)
        cdiDefAttInt(vlistID, CDI_GLOBAL, "cdo_bitrounding_numsteps", CDI_DATATYPE_INT32, 1, &params.numSteps);
    }
}

static std::vector<int>
get_vars_numbits(const VarList &varList, const std::string &filename)
{
  auto numVars = varList.size();
  std::vector<int> varsNumbits(numVars, -1);

  if (filename.size())
    {
      auto fp = std::fopen(filename.c_str(), "r");
      if (fp == nullptr) cdo_abort("Open failed on: %s\n", filename);
      PMList pmlist;
      pmlist.read_namelist(fp, filename.c_str());
      auto &kvlist = pmlist.front();
      std::fclose(fp);
      if (Options::cdoVerbose) kvlist.print();

      for (const auto &kv : kvlist)
        {
          const auto &key = kv.key;
          if (kv.nvalues > 1) cdo_abort("Too many values for parameter key >%s<!", key);
          if (kv.nvalues < 1) cdo_abort("Missing value for parameter key >%s<!", key);

          for (size_t varID = 0; varID < numVars; ++varID)
            {
              if (key == varList[varID].name)
                {
                  const auto &value = kv.values[0];
                  auto numBits = parameter_to_int(value);
                  check_range(numBits, 1, 23, key);
                  varsNumbits[varID] = numBits;
                }
            }
        }
    }

  return varsNumbits;
}

static int
num_vars_have_numbits(const std::vector<int> &varsNumbits)
{
  int numVarsHaveNumbits = 0;
  for (size_t i = 0; i < varsNumbits.size(); ++i)
    if (varsNumbits[i] != -1) numVarsHaveNumbits++;

  return numVarsHaveNumbits;
}

void *
Bitrounding(void *process)
{
  cdo_initialize(process);

  cdo_operator_add("bitrounding", 0, 0, nullptr);

  auto params = get_parameter();
  if (Options::cdoVerbose) print_parameter(params);

  auto streamID1 = cdo_open_read(0);
  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto taxisID1 = vlistInqTaxis(vlistID1);

  VarList varList1;
  varListInit(varList1, vlistID1);

  auto varsNumbits = get_vars_numbits(varList1, params.filename);
  auto numVarsHaveNumbits = num_vars_have_numbits(varsNumbits);

  auto vlistID2 = vlistDuplicate(vlistID1);
  auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  check_attributes(vlistID1);
  set_global_attributes(vlistID2, params, numVarsHaveNumbits);

  auto nvars = vlistNvars(vlistID1);
  for (int varID = 0; varID < nvars; ++varID)
    {
      if (varList1[varID].memType == MemType::Float)
        {
          int nsb = (varsNumbits[varID] != -1) ? varsNumbits[varID] : params.numBits;
          if (nsb >= 1 && nsb <= 23) set_local_attributes(vlistID2, varID, nsb);
        }
    }

  auto streamID2 = cdo_open_write(1);
  cdo_def_vlist(streamID2, vlistID2);

  std::vector<VarStat> varsStatGlob(nvars);
  std::vector<bool> varsCheckMiss(nvars, true);
  std::vector<bool> varsCheckFloat(nvars, true);

  std::vector<std::vector<int>> nsbVarLevels(nvars);
  for (int varID = 0; varID < nvars; ++varID) nsbVarLevels[varID].resize(varList1[varID].nlevels, 0);

  auto ntsteps = vlistNtsteps(vlistID1);

  Field field;

  progress::init();

  int tsID = 0;
  while (true)
    {
      auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
      if (nrecs == 0) break;

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_def_timestep(streamID2, tsID);

      std::vector<VarStat> varsStat(nvars);

      for (int recID = 0; recID < nrecs; ++recID)
        {
          auto fstatus = (ntsteps > 1) ? (tsID + (recID + 1.0) / nrecs) / ntsteps : 1.0;
          if (!Options::cdoVerbose) progress::update(0, 1, fstatus);

          int varID, levelID;
          cdo_inq_record(streamID1, &varID, &levelID);
          cdo_def_record(streamID2, varID, levelID);

          const auto &var = varList1[varID];
          field.init(var);
          cdo_read_record(streamID1, field);

          if (field.memType == MemType::Double)
            {
              if (varsCheckFloat[varID])
                {
                  varsCheckFloat[varID] = false;
                  cdo_warning("64-bit floats unsupported, bitrounding disabled for %s!", var.name);
                }
            }
          else if (field.memType == MemType::Float)
            {
              int nsb = (varsNumbits[varID] != -1) ? varsNumbits[varID] : params.numBits;

              if (field.nmiss == 0)
                {
                  if (nsb == -1 && (tsID == 0 || params.numSteps == -1))
                    {
                      nsb = bit_rounding(field.size, field.vec_f, params.infLevel);
                      // printf("nsb=%d\n", nsb);
                      if (params.addBits) nsb += params.addBits;
                      nsb = std::min(std::max(nsb, params.minBits), params.maxBits);
                    }

                  if (tsID == 0) { nsbVarLevels[varID][levelID] = nsb; }
                  else if (params.numSteps == 1) { nsb = nsbVarLevels[varID][levelID]; }
                }

              varsStat[varID].nsbMin = std::min(varsStat[varID].nsbMin, nsb);
              varsStat[varID].nsbMax = std::max(varsStat[varID].nsbMax, nsb);

              if (nsb >= 1 && nsb <= 23) bitround(nsb, field.size, field.vec_f, var.missval);

              if (nsb == -1 && field.nmiss > 0 && varsCheckMiss[varID])
                {
                  varsCheckMiss[varID] = false;
                  cdo_warning("Missing values unsupported, bitrounding disabled for %s!", var.name);
                }
            }

          cdo_write_record(streamID2, field);
        }

      if (Options::cdoVerbose && params.numBits == -1)
        {
          fprintf(stderr, "NSB: step=%d:", tsID + 1);
          for (int varID = 0; varID < nvars; ++varID)
            if (varsStat[varID].nsbMin >= 1 && varsStat[varID].nsbMin <= 23)
              {
                fprintf(stderr, " %s=%d", varList1[varID].name.c_str(), varsStat[varID].nsbMin);
                if (varList1[varID].nlevels > 1) fprintf(stderr, "-%d ", varsStat[varID].nsbMax);
              }
          fprintf(stderr, "\n");
        }

      for (int varID = 0; varID < nvars; ++varID)
        {
          varsStatGlob[varID].nsbMin = std::min(varsStatGlob[varID].nsbMin, varsStat[varID].nsbMin);
          varsStatGlob[varID].nsbMax = std::max(varsStatGlob[varID].nsbMax, varsStat[varID].nsbMax);
        }

      if (params.printBits) break;

      tsID++;
    }

  if (params.printBits)
    {
      for (int varID = 0; varID < nvars; ++varID)
        if (varsStatGlob[varID].nsbMin >= 1 && varsStatGlob[varID].nsbMin <= 23)
          fprintf(stdout, "%s=%d\n", varList1[varID].name.c_str(), varsStatGlob[varID].nsbMax);
    }
  else if (Options::cdoVerbose && params.numBits == -1)
    {
      fprintf(stderr, "NSB: step=all:");
      for (int varID = 0; varID < nvars; ++varID)
        if (varsStatGlob[varID].nsbMin >= 1 && varsStatGlob[varID].nsbMin <= 23)
          {
            fprintf(stderr, " %s=%d", varList1[varID].name.c_str(), varsStatGlob[varID].nsbMin);
            if (varsStatGlob[varID].nsbMin != varsStatGlob[varID].nsbMax) fprintf(stderr, "-%d", varsStatGlob[varID].nsbMax);
          }
      fprintf(stderr, "\n");
    }

  cdo_stream_close(streamID1);
  cdo_stream_close(streamID2);

  vlistDestroy(vlistID2);

  cdo_finish();

  return nullptr;
}
