/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Fabian Wachsmann

*/

/*
   This module contains the following operators:

      EcaEtccdi    eca_etccdi         Etccdi conform indices
*/

#include <cdi.h>

#include "cdo_options.h"
#include "cdo_vlist.h"
#include "datetime.h"
#include "process_int.h"
#include "param_conversion.h"
#include "percentiles_hist.h"
#include "percentiles.h"
#include "util_date.h"
#include "ecautil.h"
#include "ecacore.h"
#include "field_functions.h"

enum functions_select
{
  func_selle = 1,
  func_selge = 2
};

static const char TX90P_UNITS[] = "%";
static const char TX90P_NAME[] = "tx90pETCCDI";
static const char TX90P_LONGNAME[] = "Percentage of Days when Daily Maximum Temperature is Above the 90th Percentile";

static const char TX10P_UNITS[] = "%";
static const char TX10P_NAME[] = "tx10pETCCDI";
static const char TX10P_LONGNAME[] = "Percentage of Days when Daily Maximum Temperature is Below the 10th Percentile";

static const char TN90P_UNITS[] = "%";
static const char TN90P_NAME[] = "tn90pETCCDI";
static const char TN90P_LONGNAME[] = "Percentage of Days when Daily Minimum Temperature is Above the 90th Percentile";

static const char TN10P_UNITS[] = "%";
static const char TN10P_NAME[] = "tn10pETCCDI";
static const char TN10P_LONGNAME[] = "Percentage of Days when Daily Minimum Temperature is Below the 10th Percentile";

static const char R99P_UNITS[] = "mm";
static const char R99P_NAME[] = "r99pETCCDI";
static const char R99P_LONGNAME[]
    = "Annual Total Precipitation when Daily Precipitation Exceeds the 99th Percentile of Wet Day Precipitation";

static const char R95P_UNITS[] = "mm";
static const char R95P_NAME[] = "r95pETCCDI";
static const char R95P_LONGNAME[]
    = "Annual Total Precipitation when Daily Precipitation Exceeds the 95th Percentile of Wet Day Precipitation";

/* windowDays()
   Saves for each bootstrap year
   for each day of year
   the day of each window day */
static void
windowDays(int dayOfYear, std::vector<int> &wdays, const std::vector<bool> &wdaysRead, int MaxDays, int ndates, int sumboot)
{
  int wdayOfYear = dayOfYear;
  for (int gobackDays = ceil(ndates / 2. - 1); gobackDays != 0; gobackDays--)
    {
      wdayOfYear--;
      if (wdayOfYear < 1) wdayOfYear += (MaxDays - 1) * sumboot;
      while (wdaysRead[wdayOfYear] == false)
        {
          wdayOfYear--;
          if (wdayOfYear == dayOfYear) cdo_abort("Too less timesteps!");
          if (wdayOfYear < 1) wdayOfYear += (MaxDays - 1) * sumboot;
        }
    }
  int base = (dayOfYear - 1) * ndates + 1;
  wdays[base] = wdayOfYear;
  int nndates = 1;
  while (nndates != ndates)
    {
      wdayOfYear++;
      if (wdayOfYear > sumboot * (MaxDays - 1)) wdayOfYear -= sumboot * (MaxDays - 1);
      if (wdaysRead[wdayOfYear] != false)
        {
          wdays[base + nndates] = wdayOfYear;
          nndates++;
        }
    }
}

static void
writeTimesteps(int MaxMonths, int recentYear, FieldVector3D &cei, int frequency, int taxisID4, const CdoStreamID &streamID4,
               int *otsID, const VarList &varList1, const std::vector<RecordInfo> &recList, const std::vector<int> &tempdpm,
               int tempdpy, int func2)
{
  const int maxRecs = recList.size();

  if (frequency == 8)
    {
      for (int loopmonth = 2; loopmonth < MaxMonths + 2; loopmonth++)
        {
          define_mid_of_time(frequency, taxisID4, recentYear, loopmonth, MaxMonths);
          cdo_def_timestep(streamID4, *otsID);
          for (int recID = 0; recID < maxRecs; ++recID)
            {
              auto [varIDo, levelIDo] = recList[recID].get();
              if (*otsID && varList1[varIDo].isConstant) continue;

              fieldc_div(cei[loopmonth - 2][0][levelIDo], (double) (tempdpm[loopmonth - 2] / 100.0));
              cdo_def_record(streamID4, varIDo, levelIDo);
              cdo_write_record(streamID4, cei[loopmonth - 2][0][levelIDo].vec_d.data(), cei[loopmonth - 2][0][levelIDo].nmiss);
            }
          (*otsID)++;
        }
    }
  else
    {
      define_mid_of_time(frequency, taxisID4, recentYear, 0, MaxMonths);
      cdo_def_timestep(streamID4, *otsID);
      for (int recID = 0; recID < maxRecs; ++recID)
        {
          auto [varIDo, levelIDo] = recList[recID].get();
          if (*otsID && varList1[varIDo].isConstant) continue;

          if (func2 == FieldFunc_Avg) fieldc_div(cei[0][0][levelIDo], (double) (tempdpy / 100.0));
          cdo_def_record(streamID4, varIDo, levelIDo);
          cdo_write_record(streamID4, cei[0][0][levelIDo].vec_d.data(), cei[0][0][levelIDo].nmiss);
        }
      (*otsID)++;
    }
}

static void
calculateOuterPeriod(Field &field, int MaxMonths, int recentYear, int endOfCalc, FieldVector3D &cei, FieldVector3D &varsPtemp,
                     int frequency, int taxisID4, const CdoStreamID &streamID4, int *otsID, const VarList &varList1,
                     const std::vector<RecordInfo> &recList, int selection, int func2)
{
  if (cdo_assert_files_only() == false) cdo_abort("This operator can't be combined with other operators!");

  auto cdiStream = streamOpenRead(cdo_get_stream_name(0));

  const auto cdiVlistID = streamInqVlist(cdiStream);
  const auto cdiTaxisID = vlistInqTaxis(cdiVlistID);
  std::vector<int> tempdpm(MaxMonths);
  int tempdpy = 0;
  for (int i = 0; i < MaxMonths; ++i) tempdpm[i] = 0;
  int year, month, day, tsID = 0, varID, levelID;
  bool lHasStarted = false;

  if (Options::cdoVerbose) cdo_print("Start to process variables");

  while (true)
    {
      const auto nrecs = streamInqTimestep(cdiStream, tsID++);
      if (nrecs == 0) break;

      const auto vDateTime = taxisInqVdatetime(cdiTaxisID);
      cdiDate_decode(vDateTime.date, &year, &month, &day);
      if (!lHasStarted && year != recentYear)
        continue;
      else if (!lHasStarted)
        lHasStarted = true;

      if (year != recentYear)
        {
          writeTimesteps(MaxMonths, recentYear, cei, frequency, taxisID4, streamID4, otsID, varList1, recList, tempdpm, tempdpy,
                         func2);
          recentYear = year;
          tempdpy = 0;
          tempdpm[0] = 0;
          field_fill(cei[0][0][0], 0.);
          if (frequency == 8)
            for (int loopmonth = 1; loopmonth < MaxMonths; loopmonth++)
              {
                tempdpm[loopmonth] = 0;
                field_fill(cei[loopmonth][0][0], 0.);
              }
        }
      if (year == endOfCalc && func2 == FieldFunc_Avg) break;
      tempdpy++;
      auto dayOfYear = decode_day_of_year(vDateTime.date);
      tempdpm[month - 1]++;

      if (func2 == FieldFunc_Sum) dayOfYear = 1;

      for (int recID = 0; recID < nrecs; ++recID)
        {
          streamInqRecord(cdiStream, &varID, &levelID);
          streamReadRecord(cdiStream, field.vec_d.data(), &field.nmiss);

          Field &pctls = varsPtemp[dayOfYear][0][levelID];
          if (selection == func_selle)
            vfarselle(field, pctls);
          else if (selection == func_selge)
            vfarselge(field, pctls);

          auto &array = field.vec_d;
          if (func2 == FieldFunc_Avg)
            for (size_t i = 0; i < field.size; ++i) array[i] = (DBL_IS_EQUAL(array[i], field.missval)) ? 0.0 : 1.0;
          else
            for (size_t i = 0; i < field.size; ++i) array[i] = (DBL_IS_EQUAL(array[i], field.missval)) ? 0.0 : array[i];
          if (frequency == 8)
            field2_add(cei[(int) ((dayOfYear - 1) / 31.)][0][levelID], field);
          else
            field2_add(cei[0][0][levelID], field);
        }
    }
  if (Options::cdoVerbose) cdo_print("Finished Processing variables");
  if (year != endOfCalc)
    writeTimesteps(MaxMonths, year, cei, frequency, taxisID4, streamID4, otsID, varList1, recList, tempdpm, tempdpy, func2);

  field_fill(cei[0][0][0], 0.);

  if (frequency == 8)
    for (int loopmonth = 1; loopmonth < MaxMonths; loopmonth++)
      {
        tempdpm[loopmonth] = 0;
        field_fill(cei[loopmonth][0][0], 0.);
      }

  streamClose(cdiStream);
}

void
etccdi_op(ETCCDI_REQUEST &request)
{
  constexpr int MaxDays = 373;
  constexpr int MaxMonths = 12;
  FieldVector2D vars2[MaxDays];
  HistogramSet hsets[MaxDays];

  const int operatorID = cdo_operator_id();
  auto selection = cdo_operator_f1(operatorID);
  auto frequency = cdo_operator_f2(operatorID);

  percentile_set_method("rtype8");

  int FIELD_MEMTYPE = (Options::CDO_Memtype == MemType::Float) ? FIELD_FLT : 0;

  bool wdaysSrc[MaxDays];

  if (request.endboot < request.startboot)
    {
      cdo_warning("Your interval end '%d' is before the interval start '%d'. Switched interval years.", request.endboot,
                  request.startboot);
      request.startboot = request.endboot;
      request.endboot = request.startboot;
    }
  int sumboot = request.endboot - request.startboot + 1;
  std::vector<bool> wdaysRead((MaxDays - 1) * sumboot + 1);
  std::vector<int> wdays(request.ndates * (MaxDays - 1) * sumboot + 1);
  std::vector<int> dpy(sumboot), dpm(MaxMonths * sumboot);

  for (int year = 0; year < sumboot; year++)
    {
      dpy[year] = 0;
      for (int month = 0; month < MaxMonths; ++month) dpm[month + year * MaxMonths] = 0;
    }

  const auto streamID1 = cdo_open_read(0);
  const auto streamID2 = cdo_open_read(1);
  const auto streamID3 = cdo_open_read(2);

  const auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  const auto vlistID2 = cdo_stream_inq_vlist(streamID2);
  const auto vlistID3 = cdo_stream_inq_vlist(streamID3);
  const auto vlistID4 = vlistDuplicate(vlistID1);

  vlist_compare(vlistID1, vlistID2, CMP_ALL);
  vlist_compare(vlistID1, vlistID3, CMP_ALL);

  const auto taxisID1 = vlistInqTaxis(vlistID1);
  const auto taxisID2 = vlistInqTaxis(vlistID2);
  const auto taxisID3 = vlistInqTaxis(vlistID3);
  // TODO - check that time axes 2 and 3 are equal

  const auto taxisID4 = taxisDuplicate(taxisID1);
  if (taxisHasBounds(taxisID4)) taxisDeleteBounds(taxisID4);
  vlistDefTaxis(vlistID4, taxisID4);

  const auto streamID4 = cdo_open_write(3);

  const auto nvars = vlistNvars(vlistID1);
  bool lOnlyOneVar = true;

  if (nvars == 1)
    {
      cdiDefKeyString(vlistID4, 0, CDI_KEY_NAME, request.name);
      cdiDefKeyString(vlistID4, 0, CDI_KEY_LONGNAME, request.longname);
      cdiDefKeyString(vlistID4, 0, CDI_KEY_UNITS, request.units);
      cdiDefAttTxt(vlistID4, 0, "cell_methods", (int) strlen("time: maximum"), "time: maximum");
    }
  else
    {
      lOnlyOneVar = false;
      cdo_warning("Your input file has more than one variable. No attributes can be set.");
    }

  cdo_def_vlist(streamID4, vlistID4);

  const auto maxrecs = vlistNrecs(vlistID1);
  std::vector<RecordInfo> recList(maxrecs);

  const auto gridsizemax = vlistGridsizeMax(vlistID1);

  Field field;

  VarList varList1;
  varListInit(varList1, vlistID1);

  FieldVector3D vars1(sumboot * (MaxDays - 1) + 1), cei(sumboot * MaxMonths), varsPtemp(MaxDays);
  for (int dayOfYear = 0; dayOfYear < MaxDays; dayOfYear++)
    {
      fields_from_vlist(vlistID1, vars1[dayOfYear], FIELD_VEC | FIELD_MEMTYPE);
      wdaysSrc[dayOfYear] = false;
      for (int year = 0; year < sumboot; year++) wdaysRead[dayOfYear + year * (MaxDays - 1)] = false;
    }

  for (int dayOfYear = MaxDays; dayOfYear < sumboot * (MaxDays - 1) + 1; dayOfYear++)
    fields_from_vlist(vlistID1, vars1[dayOfYear], FIELD_VEC | FIELD_MEMTYPE);

  int tsID = 0;

  while (true)
    {
      const auto nrecs = cdo_stream_inq_timestep(streamID2, tsID);
      if (nrecs == 0) break;

      if (nrecs != cdo_stream_inq_timestep(streamID3, tsID))
        cdo_abort("Number of records at time step %d of %s and %s differ!", tsID + 1, cdo_get_stream_name(1),
                  cdo_get_stream_name(2));

      const auto vDateTime2 = taxisInqVdatetime(taxisID2);
      const auto vDateTime3 = taxisInqVdatetime(taxisID3);

      if (cdiDate_get(vDateTime2.date) != cdiDate_get(vDateTime3.date))
        cdo_abort("Verification dates at time step %d of %s and %s differ!", tsID + 1, cdo_get_stream_name(1),
                  cdo_get_stream_name(2));

      auto dayOfYear = decode_day_of_year(vDateTime2.date);

      if (request.func2 == FieldFunc_Sum) dayOfYear = 1;

      if (dayOfYear < 0 || dayOfYear >= MaxDays) cdo_abort("Day %d out of range!", dayOfYear);

      if (!vars2[dayOfYear].size())
        {
          wdaysSrc[dayOfYear] = true;
          fields_from_vlist(vlistID2, vars2[dayOfYear], FIELD_VEC | FIELD_MEMTYPE);
          fields_from_vlist(vlistID2, varsPtemp[dayOfYear], FIELD_VEC);
          hsets[dayOfYear].create(nvars);

          for (int varID = 0; varID < nvars; ++varID)
            hsets[dayOfYear].createVarLevels(varID, varList1[varID].nlevels, varList1[varID].gridsize);
        }

      for (int recID = 0; recID < nrecs; ++recID)
        {
          int varID, levelID;
          cdo_inq_record(streamID2, &varID, &levelID);
          cdo_read_record(streamID2, vars2[dayOfYear][varID][levelID]);
          varsPtemp[dayOfYear][varID][levelID].nmiss = vars2[dayOfYear][varID][levelID].nmiss;
        }

      for (int recID = 0; recID < nrecs; ++recID)
        {
          int varID, levelID;
          cdo_inq_record(streamID3, &varID, &levelID);
          field.init(varList1[varID]);
          cdo_read_record(streamID3, field);

          hsets[dayOfYear].defVarLevelBounds(varID, levelID, vars2[dayOfYear][varID][levelID], field);
        }
      // fieldsFree(vlistID2, vars2[dayOfYear]);
      // fields_from_vlist(vlistID2, vars2[dayOfYear], FIELD_VEC);

      tsID++;
    }

  if (Options::cdoVerbose) cdo_print("Defined the boundaries for the histograms");
  tsID = 0;
  bool lOnlyRefPeriod = true;
  int firstYear = 0, lastYear = 0;
  int year;
  while (true)
    {
      const auto nrecs = cdo_stream_inq_timestep(streamID1, tsID++);
      if (nrecs == 0) break;

      const auto vDateTime = taxisInqVdatetime(taxisID1);

      int month, day;
      cdiDate_decode(vDateTime.date, &year, &month, &day);
      if (tsID == 1)
        {
          if (year > request.startboot)
            cdo_abort("The interval start year '%d' is before infile start year '%d'.", request.startboot, year);
          firstYear = year;
        }
      lastYear = year;

      if (year >= request.startboot && year <= request.endboot)
        {
          auto dayOfYear = decode_day_of_year(vDateTime.date);
          if (dayOfYear < 0 || dayOfYear >= MaxDays) cdo_abort("Day %d out of range!", dayOfYear);
          if (wdaysSrc[dayOfYear] || request.func2 == FieldFunc_Sum)
            {
              // Variable independent ?
              wdaysRead[dayOfYear + (year - request.startboot) * (MaxDays - 1)]
                  = dayOfYear + (year - request.startboot) * (MaxDays - 1);
              dpy[year - request.startboot]++;
              dpm[(year - request.startboot) * MaxMonths + (int) ((dayOfYear - 1) / 31.)]++;
              for (int recID = 0; recID < nrecs; ++recID)
                {
                  int varID, levelID;
                  cdo_inq_record(streamID1, &varID, &levelID);
                  cdo_read_record(streamID1, vars1[dayOfYear + (year - request.startboot) * (MaxDays - 1)][varID][levelID]);

                  if (tsID == 0) recList[recID].set(varID, levelID);
                }
            }
          else
            cdo_warning("Could not find histogram minimum or maximum for day of year: '%d'", dayOfYear);
        }
      else
        lOnlyRefPeriod = false;
    }
  if (Options::cdoVerbose) cdo_print("Read in variables");

  if (year < request.endboot) cdo_abort("The interval end year '%d' is after infile end year '%d'.", request.endboot, year);

  for (year = 0; year < sumboot; year++)
    {
      fields_from_vlist(vlistID1, cei[year * MaxMonths], FIELD_VEC);
      if (frequency == 8)
        for (int month = 1; month < MaxMonths; ++month) fields_from_vlist(vlistID1, cei[year * MaxMonths + month], FIELD_VEC);
    }

  // printf("Wir beginnen nun mit der Schleife.\n");
  int bootsyear = 0;
  int subyear = 0;
  int otsID = 0;

  for (int loopdoy = 1; loopdoy < MaxDays; loopdoy++)
    {
      for (int ytoadd = 0; ytoadd < sumboot; ytoadd++)
        {
          if (wdaysRead[loopdoy + ytoadd * (MaxDays - 1)])
            {
              windowDays(loopdoy + ytoadd * (MaxDays - 1), wdays, wdaysRead, MaxDays, request.ndates, sumboot);
            }
        }
    }
  if (Options::cdoVerbose) cdo_print("Calculated window days");

  for (int varID = 0; varID < nvars; ++varID)
    {
      if (varList1[varID].isConstant) continue;
      for (int levelID = 0; levelID < varList1[varID].nlevels; ++levelID)
        {
          if (request.func2 == FieldFunc_Sum)
            {
              for (int loopdoy = 1; loopdoy < MaxDays; loopdoy++)
                {
                  for (int ytoadd = 0; ytoadd < sumboot; ytoadd++)
                    {
                      if (wdaysRead[loopdoy + ytoadd * (MaxDays - 1)])
                        {
                          auto &source = vars1[loopdoy + ytoadd * (MaxDays - 1)][varID][levelID];
                          auto &hset = hsets[1];
                          hset.addVarLevelValues(varID, levelID, source);
                        }
                    }
                }
            }
          else
            {
#ifdef _OPENMP
#pragma omp parallel for shared(sumboot, vars1, request, varID, levelID, hsets, wdays, wdaysRead) schedule(dynamic)
#endif
              for (int loopdoy = 1; loopdoy < MaxDays; loopdoy++)
                {
                  for (int ytoadd = 0; ytoadd < sumboot; ytoadd++)
                    {
                      if (wdaysRead[loopdoy + ytoadd * (MaxDays - 1)])
                        {
                          for (int ano = 0; ano < request.ndates; ano++)
                            {
                              auto &source
                                  = vars1[wdays[ytoadd * request.ndates * (MaxDays - 1) + (loopdoy - 1) * request.ndates + ano + 1]]
                                         [varID][levelID];
                              auto &hset = hsets[loopdoy];
                              hset.addVarLevelValues(varID, levelID, source);
                            }
                        }
                    }
                }
            }
        }
    }

  tsID = 0;
  if (Options::cdoVerbose) cdo_print("Added 30 years to histograms");
  if (lOnlyOneVar && ((!lOnlyRefPeriod && firstYear != request.startboot) || request.func2 == FieldFunc_Sum))
    {
      for (int varID = 0; varID < nvars; ++varID)
        {
          if (varList1[varID].isConstant) continue;
          for (int levelID = 0; levelID < varList1[varID].nlevels; ++levelID)
            {
              if (request.func2 == FieldFunc_Sum)
                {
                  auto &pctls = varsPtemp[1][varID][levelID];
                  hsets[1].getVarLevelPercentiles(pctls, varID, levelID, request.pn);
                }
              else
                {
#ifdef _OPENMP
#pragma omp parallel for shared(request, wdaysSrc, varID, levelID, hsets, varsPtemp) schedule(dynamic)
#endif
                  for (int loopdoy = 1; loopdoy < MaxDays; loopdoy++)
                    {
                      if (wdaysSrc[loopdoy])
                        {
                          auto &pctls = varsPtemp[loopdoy][varID][levelID];
                          hsets[loopdoy].getVarLevelPercentiles(pctls, varID, levelID, request.pn);
                        }
                    }
                }
            }
        }
      field.resize(gridsizemax);
      calculateOuterPeriod(field, MaxMonths, firstYear, request.startboot, cei, varsPtemp, frequency, taxisID4, streamID4, &otsID,
                           varList1, recList, selection, request.func2);
    }
  else if (!lOnlyRefPeriod && firstYear != request.startboot)
    cdo_warning("Since you have more than one variable in the input file, only the bootstrapping period can be calculated");

  if (request.func2 == FieldFunc_Avg)
    {
      for (bootsyear = request.startboot; bootsyear < request.endboot + 1; bootsyear++)
        {
          if (Options::cdoVerbose) cdo_print("Bootsyear: %d", bootsyear);
          for (int varID = 0; varID < nvars; ++varID)
            {
              if (varList1[varID].isConstant) continue;
              for (int levelID = 0; levelID < varList1[varID].nlevels; ++levelID)
                {
#ifdef _OPENMP
#pragma omp parallel for shared(sumboot, wdaysRead, request, vars1, varID, levelID, hsets, wdays, cei) schedule(dynamic)
#endif
                  for (int loopdoy = 1; loopdoy < MaxDays; loopdoy++)
                    {
                      for (int ytoadd = 0; ytoadd < sumboot; ytoadd++)
                        {
                          if (wdaysRead[loopdoy + ytoadd * (MaxDays - 1)])
                            {
                              for (int ano = 0; ano < request.ndates; ano++)
                                {
                                  int recentWday
                                      = ytoadd * request.ndates * (MaxDays - 1) + (loopdoy - 1) * request.ndates + ano + 1;
                                  if ((int((wdays[recentWday] - 1) / (MaxDays - 1)) + request.startboot) == bootsyear)
                                    {
                                      auto &source = vars1[wdays[recentWday]][varID][levelID];
                                      auto &hset = hsets[loopdoy];
                                      hset.subVarLevelValues(varID, levelID, source);
                                    }
                                  // percyear cannot be smaller than request.startboot
                                  if ((int((wdays[recentWday] - 1) / (MaxDays - 1)) + request.startboot) == bootsyear - 1)
                                    {
                                      auto &source = vars1[wdays[recentWday]][varID][levelID];
                                      auto &hset = hsets[loopdoy];
                                      hset.addVarLevelValues(varID, levelID, source);
                                    }
                                }
                            }
                        }
                    }
                  for (subyear = request.startboot; subyear < request.endboot + 1; subyear++)
                    {
                      if (Options::cdoVerbose) cdo_print("Subyear: %d", subyear);
                      if (subyear != bootsyear)
                        {
#ifdef _OPENMP
#pragma omp parallel for shared(sumboot, request, vars1, varID, levelID, hsets, wdaysRead, varsPtemp, vars2, cei, subyear, \
                                bootsyear, wdays, frequency) schedule(dynamic)
#endif
                          for (int loopdoy = 1; loopdoy < MaxDays; loopdoy++)
                            {
                              for (int ytoadd = 0; ytoadd < sumboot; ytoadd++)
                                {
                                  if (wdaysRead[loopdoy + ytoadd * (MaxDays - 1)])
                                    {
                                      for (int ano = 0; ano < request.ndates; ano++)
                                        {
                                          int recentWday
                                              = ytoadd * request.ndates * (MaxDays - 1) + (loopdoy - 1) * request.ndates + ano + 1;
                                          if ((int((wdays[recentWday] - 1) / (MaxDays - 1)) + request.startboot) == subyear)
                                            {
                                              auto &source = vars1[wdays[recentWday]][varID][levelID];
                                              auto &hset = hsets[loopdoy];
                                              if (hset.addVarLevelValues(varID, levelID, source) == 1)
                                                cdo_print("'%d', '%d", loopdoy, wdays[recentWday]);
                                            }
                                        }
                                    }
                                }
                              // printf("Haben es zum temp array addiert.\n");

                              /*** Calculate percentile  ***/
                              if (wdaysRead[loopdoy + (bootsyear - request.startboot) * (MaxDays - 1)])
                                {
                                  auto &pctls = varsPtemp[loopdoy][varID][levelID];
                                  hsets[loopdoy].getVarLevelPercentiles(pctls, varID, levelID, request.pn);
                                  /*** Compare data with percentile ***/
                                  auto &source = vars1[loopdoy + (bootsyear - request.startboot) * (MaxDays - 1)][varID][levelID];
                                  auto &toCompare = vars2[loopdoy][varID][levelID];
                                  field_copy(source, toCompare);
                                  if (selection == func_selle)
                                    vfarselle(toCompare, pctls);
                                  else if (selection == func_selge)
                                    vfarselge(toCompare, pctls);
                                  if (request.func2 == FieldFunc_Avg)
                                    {
                                      auto &array = toCompare.vec_d;
                                      for (size_t i = 0; i < toCompare.size; ++i)
                                        array[i] = (DBL_IS_EQUAL(array[i], toCompare.missval)) ? 0.0 : 1.0;
                                    }
                                  else
                                    {
                                      auto &array = toCompare.vec_d;
                                      for (size_t i = 0; i < toCompare.size; ++i)
                                        array[i] = (DBL_IS_EQUAL(array[i], toCompare.missval)) ? 0.0 : array[i];
                                    }
                                  // printf("Haben ein Percentil berechnet.\n");
                                  // Year sum
                                  if (frequency == 8)
#ifdef _OPENMP
#pragma omp critical
#endif
                                    field2_add(cei[(bootsyear - request.startboot) * MaxMonths + (int) ((loopdoy - 1) / 31.)][varID]
                                                  [levelID],
                                               toCompare);
                                  else
#ifdef _OPENMP
#pragma omp critical
#endif
                                    field2_add(cei[(bootsyear - request.startboot) * MaxMonths][varID][levelID], toCompare);
                                }
                              for (int ytoadd = 0; ytoadd < sumboot; ytoadd++)
                                {
                                  if (wdaysRead[loopdoy + ytoadd * (MaxDays - 1)])
                                    {
                                      for (int ano = 0; ano < request.ndates; ano++)
                                        {
                                          int recentWday
                                              = ytoadd * request.ndates * (MaxDays - 1) + (loopdoy - 1) * request.ndates + ano + 1;
                                          if ((int((wdays[recentWday] - 1) / (MaxDays - 1)) + request.startboot) == subyear)
                                            {
                                              auto &source = vars1[wdays[recentWday]][varID][levelID];
                                              auto &hset = hsets[loopdoy];
                                              if (hset.subVarLevelValues(varID, levelID, source) == 1)
                                                cdo_print("'%d', '%d", loopdoy, wdays[recentWday]);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                  if (frequency == 8)
                    {
                      for (int month = 0; month < MaxMonths; ++month)
                        if (cei[(bootsyear - request.startboot) * MaxMonths + month][varID][levelID].vec_d.data())
                          {
                            // Divide vars2 to receive average
                            fieldc_div(cei[(bootsyear - request.startboot) * MaxMonths + month][varID][levelID],
                                       (double) ((sumboot - 1) * dpm[(bootsyear - request.startboot) * MaxMonths + month] / 100.0));
                          }
                    }
                  else if (cei[(bootsyear - request.startboot) * MaxMonths][varID][levelID].vec_d.data())
                    {
                      fieldc_div(cei[(bootsyear - request.startboot) * MaxMonths][varID][levelID],
                                 (double) ((sumboot - 1) * dpy[bootsyear - request.startboot] / 100.0));
                    }
                }
            }
          if (frequency == 8)
            {
              for (int month = 2; month < MaxMonths + 2; ++month)
                {
                  define_mid_of_time(frequency, taxisID4, bootsyear, month, MaxMonths);
                  cdo_def_timestep(streamID4, otsID);

                  for (int recID = 0; recID < maxrecs; ++recID)
                    {
                      auto [varIDo, levelIDo] = recList[recID].get();
                      if (otsID && varList1[varIDo].isConstant) continue;

                      cdo_def_record(streamID4, varIDo, levelIDo);
                      cdo_write_record(
                          streamID4, cei[(bootsyear - request.startboot) * MaxMonths + (month - 2)][varIDo][levelIDo].vec_d.data(),
                          cei[(bootsyear - request.startboot) * MaxMonths + (month - 2)][varIDo][levelIDo].nmiss);
                    }
                  otsID++;
                }
            }
          else
            {
              define_mid_of_time(frequency, taxisID4, bootsyear, 0, MaxMonths);
              cdo_def_timestep(streamID4, otsID);

              for (int recID = 0; recID < maxrecs; ++recID)
                {
                  auto [varIDo, levelIDo] = recList[recID].get();
                  if (otsID && varList1[varIDo].isConstant) continue;

                  cdo_def_record(streamID4, varIDo, levelIDo);
                  cdo_write_record(streamID4, cei[(bootsyear - request.startboot) * MaxMonths][varIDo][levelIDo].vec_d.data(),
                                   cei[(bootsyear - request.startboot) * MaxMonths][varIDo][levelIDo].nmiss);
                }
              otsID++;
            }
          // printf("Haben ein Mittel für Jahr '%d' berechnet.\n", bootsyear);
        }
    }
  if (!lOnlyRefPeriod && lOnlyOneVar && lastYear != request.endboot && request.func2 == FieldFunc_Avg)
    {
      field_fill(cei[0][0][0], 0.);
      if (frequency == 8)
        for (int loopmonth = 1; loopmonth < MaxMonths; loopmonth++) field_fill(cei[loopmonth][0][0], 0.);

      for (int varID = 0; varID < nvars; ++varID)
        {
          if (varList1[varID].isConstant) continue;
          for (int levelID = 0; levelID < varList1[varID].nlevels; ++levelID)
            {
#ifdef _OPENMP
#pragma omp parallel for shared(request, wdaysRead, varID, levelID, hsets, varsPtemp) schedule(dynamic)
#endif
              for (int loopdoy = 1; loopdoy < MaxDays; loopdoy++)
                {
                  for (int ytoadd = 0; ytoadd < sumboot; ytoadd++)
                    {
                      if (wdaysRead[loopdoy + ytoadd * (MaxDays - 1)])
                        {
                          for (int ano = 0; ano < request.ndates; ano++)
                            {
                              int recentWday = ytoadd * request.ndates * (MaxDays - 1) + (loopdoy - 1) * request.ndates + ano + 1;
                              // percyear cannot be smaller than request.startboot
                              if ((int((wdays[recentWday] - 1) / (MaxDays - 1)) + request.startboot) == bootsyear - 1)
                                {
                                  auto &source = vars1[wdays[recentWday]][varID][levelID];
                                  auto &hset = hsets[loopdoy];
                                  hset.addVarLevelValues(varID, levelID, source);
                                }
                            }
                        }
                    }
                  if (wdaysSrc[loopdoy])
                    {
                      auto &pctls = varsPtemp[loopdoy][varID][levelID];
                      hsets[loopdoy].getVarLevelPercentiles(pctls, varID, levelID, request.pn);
                    }
                }
            }
        }
      field.resize(gridsizemax);
      field.missval = vars1[1][0][0].missval;
      field.size = vars1[1][0][0].size;
      field.grid = vars1[1][0][0].grid;
      calculateOuterPeriod(field, MaxMonths, request.endboot + 1, lastYear + 1, cei, varsPtemp, frequency, taxisID4, streamID4,
                           &otsID, varList1, recList, selection, request.func2);
    }

  cdo_stream_close(streamID4);
  cdo_stream_close(streamID3);
  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);
}

void *
EcaEtccdi(void *process)
{
  cdo_initialize(process);

  if (cdo_operator_argc() < 3) cdo_abort("Too few arguments!");
  if (cdo_operator_argc() > 4) cdo_abort("Too many arguments!");

  int TX90P, TX10P, TN90P, TN10P, ALLX, R99P, R95P;
  if (cdo_operator_argc() == 4 && 'm' == cdo_operator_argv(3)[0])
    {
      TX90P = cdo_operator_add("etccdi_tx90p", func_selge, CMP_MONTH, nullptr);  // monthly mode
      R99P = cdo_operator_add("etccdi_r99p", func_selge, CMP_MONTH, nullptr);    // monthly mode
      R95P = cdo_operator_add("etccdi_r95p", func_selge, CMP_MONTH, nullptr);    // monthly mode
      TX10P = cdo_operator_add("etccdi_tx10p", func_selle, CMP_MONTH, nullptr);  // monthly mode
      TN90P = cdo_operator_add("etccdi_tn90p", func_selge, CMP_MONTH, nullptr);  // monthly mode
      TN10P = cdo_operator_add("etccdi_tn10p", func_selle, CMP_MONTH, nullptr);  // monthly mode
      ALLX = cdo_operator_add("etccdi", 0, CMP_MONTH, nullptr);                  // monthly mode
    }
  else
    {
      TX90P = cdo_operator_add("etccdi_tx90p", func_selge, CMP_DATE, nullptr);
      R99P = cdo_operator_add("etccdi_r99p", func_selge, CMP_DATE, nullptr);
      R95P = cdo_operator_add("etccdi_r95p", func_selge, CMP_DATE, nullptr);
      TX10P = cdo_operator_add("etccdi_tx10p", func_selle, CMP_DATE, nullptr);
      TN90P = cdo_operator_add("etccdi_tn90p", func_selge, CMP_DATE, nullptr);
      TN10P = cdo_operator_add("etccdi_tn10p", func_selle, CMP_DATE, nullptr);
      ALLX = cdo_operator_add("etccdi", 0, CMP_DATE, nullptr);
    }

  ETCCDI_REQUEST request;

  request.ndates = parameter_to_int(cdo_operator_argv(0));
  request.startboot = parameter_to_int(cdo_operator_argv(1));

  const auto operatorID = cdo_operator_id();
  if (operatorID == TX90P || operatorID == TN90P || operatorID == R95P || operatorID == R99P)
    {
      if (operatorID == TX90P || operatorID == TN90P)
        {
          if (cdo_operator_argc() < 3)
            cdo_abort("Operator requires at least 3 parameter values, you provided '%d'!", cdo_operator_argc());
          request.endboot = parameter_to_int(cdo_operator_argv(2));
          if (operatorID == TX90P)
            {
              request.name = TX90P_NAME;
              request.longname = TX90P_LONGNAME;
              request.units = TX90P_UNITS;
              request.func2 = FieldFunc_Avg;
            }
          else if (operatorID == TN90P)
            {
              request.name = TN90P_NAME;
              request.longname = TN90P_LONGNAME;
              request.units = TN90P_UNITS;
              request.func2 = FieldFunc_Avg;
            }
          request.pn = 90;
        }
      else
        {
          if (cdo_operator_argc() < 2)
            cdo_abort("Operator requires at least 2 parameter values, you provided '%d'!", cdo_operator_argc());
          request.ndates = 1;
          request.startboot = parameter_to_int(cdo_operator_argv(0));
          request.endboot = parameter_to_int(cdo_operator_argv(1));
          if (operatorID == R95P)
            {
              request.name = R95P_NAME;
              request.longname = R95P_LONGNAME;
              request.units = R95P_UNITS;
              request.pn = 95;
              request.func2 = FieldFunc_Sum;
            }
          else if (operatorID == R99P)
            {
              request.name = R99P_NAME;
              request.longname = R99P_LONGNAME;
              request.units = R99P_UNITS;
              request.pn = 99;
              request.func2 = FieldFunc_Sum;
            }
        }
    }
  else if (operatorID == TX10P || operatorID == TN10P)
    {
      if (cdo_operator_argc() < 3)
        cdo_abort("Operator requires at least 3 parameter values, you provided '%d'!", cdo_operator_argc());
      request.endboot = parameter_to_int(cdo_operator_argv(2));
      if (operatorID == TX10P)
        {
          request.name = TX10P_NAME;
          request.longname = TX10P_LONGNAME;
          request.units = TX10P_UNITS;
          request.func2 = FieldFunc_Avg;
        }
      else
        {
          request.name = TN10P_NAME;
          request.longname = TN10P_LONGNAME;
          request.units = TN10P_UNITS;
          request.func2 = FieldFunc_Avg;
        }
      request.pn = 10;
    }

  etccdi_op(request);
  /*  else
      EcaEtccdi(-1, ndates, startboot, endboot); */

  cdo_finish();

  return nullptr;
}
