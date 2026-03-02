/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Copyright (C) 2006 Brockmann Consult

  Author: Ralf Quast

*/

// clang-format off
/*
      MODULE      OPERATOR     INDEX    DESCRIPTION

      EcaCfd      eca_cfd      CFD      maximum number of consecutive frost days
      EcaCsu      eca_csu      CSU      maximum number of consecutive summer days
      EcaCwdi     eca_cwdi     CWDI     cold wave duration index 
      EcaCwfi     eca_cwfi     CWFI     number of cold-spell days
      EcaEtr      eca_etr      ETR      intra-period extreme temperature range
      EcaFd       eca_fd       FD       number of frost days
      EcaGsl      eca_gsl      GSL      growing season length
      EcaHd       eca_hd       HD       heating degree days
      EcaHwdi     eca_hwdi     HWDI     heat wave duration index
      EcaHwfi     eca_hwfi     HWFI     number of warm-spell days
      EcaId       eca_id       ID       number of ice days
      EcaSu       eca_su       SU       number of summer days
      EcaTg10p    eca_tg10p    TG10p    percent of time TX < 10th percentile of daily mean temperature
      EcaTg90p    eca_tg90p    TG90p    percent of time TX > 90th percentile of daily mean temperature
      EcaTn10p    eca_tn10p    TN10p    percent of time TX < 10th percentile of daily minimum temperature
      EcaTn90p    eca_tn90p    TN90p    percent of time TX > 90th percentile of daily minimum temperature
      EcaTr       eca_tr       TR       number of tropical nights
      EcaTx10p    eca_tx10p    TX10p    percent of time TX < 10th percentile of daily maximum temperature
      EcaTx90p    eca_tx90p    TX90p    percent of time TX > 90th percentile of daily maximum temperature

      EcaCdd      eca_cdd      CDD      maximum number of consecutive dry days
      EcaCwd      eca_cwd      CWD      maximum number of consecutive wet days
      EcaR10mm    eca_r10mm    R10mm    number of days with precipitation >= 10 mm
      EcaR20mm    eca_r20mm    R20mm    number of days with precipitation >= 20 mm
      EcaR75p     eca_r75p     R75p     Percent of time RR > 75th percentile of daily precipitation amount
      EcaR75ptot  eca_r75ptot  R75pTOT  Percentage of annual total precipitation due to events with RR > 75th percentile of daily precipitation amount
      EcaR90p     eca_r90p     R90p     Percent of time RR > 90th percentile of daily precipitation amount
      EcaR90ptot  eca_r90ptot  R90pTOT  Percentage of annual total precipitation due to events with RR > 90th percentile of daily precipitation amount
      EcaR95p     eca_r95p     R95p     Percent of time RR > 95th percentile of daily precipitation amount
      EcaR95ptot  eca_r95ptot  R95pTOT  Percentage of annual total precipitation due to events with RR > 95th percentile of daily precipitation amount
      EcaR99p     eca_r99p     R99p     Percent of time RR > 75th percentile of daily precipitation amount
      EcaR99ptot  eca_r99ptot  R99pTOT  Percentage of annual total precipitation due to events with RR > 99th percentile of daily precipitation amount
      EcaRr1      eca_rr1      RR1      number of wet days
      EcaSdii     eca_sdii     SDII     simple daily intensity index

      Fdns        fdns                  frost days without surface snow

      Strwin      strwin                number of strong-wind days
      Strbre      strbre                number of strong-breeze days
      Strgal      strgal                number of strong-gale days
      Hurr        hurr                  number of hurricane days
*/
// clang-format on

#include "process_int.h"
#include "cdo_options.h"
#include "param_conversion.h"
#include "ecacore.h"
#include "ecautil.h"
#include "util_date.h"
#include "pmlist.h"
#include "field_functions.h"

#define TO_DEG_CELSIUS(x) ((x) -273.15)
#define TO_KELVIN(x) ((x) + 273.15)

constexpr int ECA_refdate = 19550101;
constexpr int ETC_refdate = 18500101;

// clang-format off

static const char CFD_NAME[]         = "consecutive_frost_days_index_per_time_period";
static const char CFD_LONGNAME[]     = "Consecutive frost days index is the greatest number of consecutive frost days in a given time period. Frost days is the number of days where minimum of temperature is below 0 degree Celsius. The time period should be defined by the bounds of the time coordinate.";
//static const char CFD_UNITS[]        = "No.";
static const char CFD_NAME2[]        = "number_of_cfd_periods_with_more_than_%ddays_per_time_period";
static const char CFD_LONGNAME2[]    = "Number of cfd periods in given time period with more than %d days. The time period should be defined by the bounds of the time coordinate.";
static const char CFD_UNITS2[]       = "No.";

static const char CSU_NAME[]         = "consecutive_summer_days_index_per_time_period";
static const char CSU_LONGNAME[]     = "Consecutive summer days index is the greatest number of consecutive summer days in a given time period. Summer days is the number of days where maximum of temperature is above 25 degree Celsius. The time period should be defined by the bounds of the time coordinate.";
//static const char CSU_UNITS[]        = "No.";
static const char CSU_NAME2[]        = "number_of_csu_periods_with_more_than_%ddays_per_time_period";
static const char CSU_LONGNAME2[]    = "Number of csu periods in given time period with more than %d days. The time period should be defined by the bounds of the time coordinate.";
static const char CSU_UNITS2[]       = "No.";

static const char CWDI_NAME[]        = "cold_wave_duration_index_wrt_mean_of_reference_period";
static const char CWDI_LONGNAME[]    = "This is the number of days per time period where in intervals of at least %d consecutive days the daily minimum temperature is more than %1.0f degrees below a reference value. The reference value is calculated  as the mean of minimum temperatures of a five day window centred on each calendar day of a given 30 year climate reference period. The time period should be defined by the bounds of the time coordinate.";
static const char CWDI_UNITS[]       = "No.";
static const char CWDI_NAME2[]       = "cold_waves_per_time_period";
static const char CWDI_LONGNAME2[]   = "Number of cold waves per time period. The time period should be defined by the bounds of the time coordinate.";
static const char CWDI_UNITS2[]      = "No.";

static const char CWFI_NAME[]        = "cold_spell_days_index_wrt_10th_percentile_of_reference_period";
static const char CWFI_NAME_ET[]     = "csdiETCCDI";
static const char CWFI_LONGNAME[]    = "This is the number of days per time period where in intervals of at least %d consecutive days the daily mean temperature is below a reference value. The reference value is calculated  as the 10th percentile of daily mean temperatures of a five day window centred on each calendar day of a given 30 year climate reference period. The time period should be defined by the bounds of the time coordinate.";
static const char CWFI_LONGNAME_ET[] = "Cold Spell Duration Index";
static const char CWFI_UNITS[]       = "No.";
static const char CWFI_UNITS_ET[]    = "days";
static const char CWFI_NAME2[]       = "cold_spell_periods_per_time_period";
static const char CWFI_LONGNAME2[]   = "Number of cold spell periods per time period. The time period should be defined by the bounds of the time coordinate.";
static const char CWFI_UNITS2[]      = "No.";

static const char ETR_NAME[]         = "intra_period_extreme_temperature_range";
static const char ETR_LONGNAME[]     = "Difference between the absolute extreme temperatures in observation period. The time period should be defined by the bounds of the time coordinate.";
//static const char ETR_UNITS[]        = "K";

static const char FD_NAME[]          = "frost_days_index_per_time_period";
static const char FD_NAME_ET[]       = "fdETCCDI";
static const char FD_LONGNAME[]      = "Frost days index is the number of days where minimum of temperature is below 0 degree Celsius. The time period should be defined by the bounds of the time coordinate.";
static const char FD_LONGNAME_ET[]   = "Number of Frost Days";
//static const char FD_UNITS[]         = "No.";
static const char FD_UNITS_ET[]      = "days";

static const char GSL_NAME[]         = "thermal_growing_season_length";
static const char GSL_LONGNAME[]     = "Counted are the number of days per calendar year between the first occurrence of at least %d consecutive days where the daily mean temperature is above %1.0f degree Celsius and the first occurrence of at least %d consecutive days after 1st of July where the daily mean temperature is below %1.0f degree Celsius. The time period should be defined by the bounds of the time coordinate.";
static const char GSL_UNITS[]        = "No.";
static const char GSL_NAME2[]        = "day_of_year_of_growing_season_start";
static const char GSL_LONGNAME2[]    = "Day of year of growing season start. The time period should be defined by the bounds of the time coordinate.";
static const char GSL_UNITS2[]       = "No.";

static const char HD_NAME[]          = "heating_degree_days_per_time_period";
static const char HD_LONGNAME[]      = "Heating degree days relates the outside temperature with the room temperature during the heating period. It is the sum of the difference between room temperature X and daily mean temperature Y on days where Y is below a given constant A. X is 20 degree Celsius and A is 15 degree Celsius according to VDI guidelines. According to ECAD both X and A are 17 degree Celsius. The time period should be defined by the bounds of the time coordinate.";
static const char HD_UNITS[]         = "No.";

static const char HWDI_NAME[]        = "heat_wave_duration_index_wrt_mean_of_reference_period";
static const char HWDI_LONGNAME[]    = "This is the number of days per time period where in intervals of at least %d consecutive days the daily maximum temperature is more than %1.0f degrees above a reference value. The reference value is calculated  as the mean of maximum temperatures of a five day window centred on each calendar day of a given 30 year climate reference period. The time period should be defined by the bounds of the time coordinate.";
static const char HWDI_UNITS[]       = "No.";
static const char HWDI_NAME2[]       = "heat_waves_per_time_period";
static const char HWDI_LONGNAME2[]   = "Number of heat waves per time period. The time period should be defined by the bounds of the time coordinate.";
static const char HWDI_UNITS2[]      = "No.";

static const char HWFI_NAME[]        = "warm_spell_days_index_wrt_90th_percentile_of_reference_period";
static const char HWFI_NAME_ET[]     = "wsdiETCCDI";
static const char HWFI_LONGNAME[]    = "This is the number of days per time period where in intervals of at least %d consecutive days the daily mean temperature is above a reference value. The reference value is calculated  as the 90th percentile of daily mean temperatures of a five day window centred on each calendar day of a given 30 year climate reference period. The time period should be defined by the bounds of the time coordinate.";
static const char HWFI_LONGNAME_ET[] = "Warm Spell Duration Index";
static const char HWFI_UNITS[]       = "No.";
static const char HWFI_UNITS_ET[]    = "days";
static const char HWFI_NAME2[]       = "warm_spell_periods_per_time_period";
static const char HWFI_LONGNAME2[]   = "Number of warm spell periods per time period. The time period should be defined by the bounds of the time coordinate.";
static const char HWFI_UNITS2[]      = "No.";

static const char ID_NAME[]          = "ice_days_index_per_time_period";
static const char ID_NAME_ET[]       = "idETCCDI";
static const char ID_LONGNAME[]      = "Ice days index is the number of days where maximum of temperature is below 0 degree Celsius. The time period should be defined by the bounds of the time coordinate.";
static const char ID_LONGNAME_ET[]   = "Number of Icing Days";
static const char ID_UNITS[]         = "No.";
static const char ID_UNITS_ET[]      = "days";

static const char SU_NAME[]          = "summer_days_index_per_time_period";
static const char SU_NAME_ET[]       = "suETCCDI";
static const char SU_LONGNAME[]      = "Summer days index is the number of days where maximum of temperature is above %1.0f degree Celsius. The time period should be defined by the bounds of the time coordinate.";
static const char SU_LONGNAME_ET[]   = "Number of Summer Days";
//static const char SU_UNITS[]         = "No.";
static const char SU_UNITS_ET[]      = "days";

static const char TG10P_NAME[]       = "cold_days_percent_wrt_10th_percentile_of_reference_period";
static const char TG10P_LONGNAME[]   = "This is the percent of time per time period where daily mean temperature is below a reference value. The reference value is calculated as the 10th percentile of daily mean temperatures of a five day window centred on each calendar day of a given 30 year climate reference period. The time period should be defined by the bounds of the time coordinate.";
static const char TG10P_UNITS[]      = "Percent";

static const char TG90P_NAME[]       = "warm_days_percent_wrt_90th_percentile_of_reference_period";
static const char TG90P_LONGNAME[]   = "This is the percent of time per time period where daily mean temperature is above a reference value. The reference value is calculated as the 90th percentile of daily mean temperatures of a five day window centred on each calendar day of a given 30 year climate reference period. The time period should be defined by the bounds of the time coordinate.";
static const char TG90P_UNITS[]      = "Percent";

static const char TN10P_NAME[]       = "cold_nights_percent_wrt_10th_percentile_of_reference_period";
static const char TN10P_LONGNAME[]   = "This is the percent of time per time period where daily minimum temperature is below a reference value. The reference value is calculated as the 10th percentile of daily minimum temperatures of a five day window centred on each calendar day of a given 30 year climate reference period. The time period should be defined by the bounds of the time coordinate.";
static const char TN10P_UNITS[]      = "Percent";

static const char TN90P_NAME[]       = "warm_nights_percent_wrt_90th_percentile_of_reference_period";
static const char TN90P_LONGNAME[]   = "This is the percent of time per time period where daily minimum temperature is above a reference value. The reference value is calculated as the 90th percentile of daily minimum temperatures of a five day window centred on each calendar day of a given 30 year climate reference period. The time period should be defined by the bounds of the time coordinate.";
static const char TN90P_UNITS[]      = "Percent";

static const char TR_NAME[]          = "tropical_nights_index_per_time_period";
static const char TR_NAME_ET[]       = "trETCCDI";
static const char TR_LONGNAME[]      = "Tropical nights index is the number of days where minimum of temperature is above %1.0f degree Celsius. The time period should be defined by the bounds of the time coordinate.";
static const char TR_LONGNAME_ET[]   = "Number of Tropical Nights";
static const char TR_UNITS[]         = "No.";
static const char TR_UNITS_ET[]      = "days";

static const char TX10P_NAME[]       = "very_cold_days_percent_wrt_10th_percentile_of_reference_period";
static const char TX10P_LONGNAME[]   = "This is the percent of time per time period where daily maximum temperature is below a reference value. The reference value is calculated as the 10th percentile of daily maximum temperatures of a five day window centred on each calendar day of a given 30 year climate reference period. The time period should be defined by the bounds of the time coordinate.";
static const char TX10P_UNITS[]      = "Percent";

static const char TX90P_NAME[]       = "very_warm_days_percent_wrt_90th_percentile_of_reference_period";
static const char TX90P_LONGNAME[]   = "This is the percent of time per time period where daily maximum temperature is above a reference value. The reference value is calculated as the 90th percentile of daily maximum temperatures of a five day window centred on each calendar day of a given 30 year climate reference period. The time period should be defined by the bounds of the time coordinate.";
static const char TX90P_UNITS[]      = "Percent";

static const char CDD_NAME[]         = "consecutive_dry_days_index_per_time_period";
static const char CDD_NAME_ET[]      = "cddETCCDI";
static const char CDD_LONGNAME[]     = "Consecutive dry days is the greatest number of consecutive days per time period with daily precipitation amount below %g mm. The time period should be defined by the bounds of the time coordinate.";
static const char CDD_LONGNAME_ET[]  = "Maximum Number of Consecutive Days with Less Than 1mm of Precipitation [days]";
static const char CDD_UNITS[]        = "No.";
static const char CDD_UNITS_ET[]      = "days";
static const char CDD_NAME2[]        = "number_of_cdd_periods_with_more_than_%ddays_per_time_period";
static const char CDD_LONGNAME2[]    = "Number of cdd periods in given time period with more than %d days. The time period should be defined by the bounds of the time coordinate.";
static const char CDD_UNITS2[]       = "No.";

static const char CWD_NAME[]         = "consecutive_wet_days_index_per_time_period";
static const char CWD_NAME_ET[]      = "cwdETCCDI";
static const char CWD_LONGNAME[]     = "Consecutive wet days is the greatest number of consecutive days per time period with daily precipitation above %g mm. The time period should be defined by the bounds of the time coordinate.";
static const char CWD_LONGNAME_ET[]  = "Maximum Number of Consecutive Days with At Least 1mm of Precipitation";
static const char CWD_UNITS[]        = "No.";
static const char CWD_UNITS_ET[]     = "days";
static const char CWD_NAME2[]        = "number_of_cwd_periods_with_more_than_%ddays_per_time_period";
static const char CWD_LONGNAME2[]    = "Number of cwd periods in given time period with more than %d days. The time period should be defined by the bounds of the time coordinate.";
static const char CWD_UNITS2[]       = "No.";

static const char PD_NAME[]          = "precipitation_days_index_per_time_period";
static const char PD_NAME_ET[]       = "r1mmETCCDI";
static const char PD_LONGNAME[]      = "precipitation days is the number of days per time period with daily precipitation sum exceeding %g mm. The time period should be defined by the bounds of the time coordinate.";
static const char PD_LONGNAME_ET[]   = "Count of Days with At Least 1mm of Precipitation";
static const char PD_UNITS[]         = "No.";
static const char PD_UNITS_ET[]      = "days";

static const char R10MM_NAME[]       = "heavy_precipitation_days_index_per_time_period";
static const char R10MM_NAME_ET[]    = "r10mmETCCDI";
static const char R10MM_LONGNAME[]   = "Heavy precipitation days is the number of days per time period with daily precipitation sum exceeding 10 mm. The time period should be defined by the bounds of the time coordinate.";
static const char R10MM_LONGNAME_ET[]= "Count of Days with At Least 10mm of Precipitation";
static const char R10MM_UNITS[]      = "No.";
static const char R10MM_UNITS_ET[]   = "days";

static const char R20MM_NAME[]       = "very_heavy_precipitation_days_index_per_time_period";
static const char R20MM_NAME_ET[]    = "r20mmETCCDI";
static const char R20MM_LONGNAME[]   = "Very heavy precipitation days is the number of days with daily precipitation sum exceeding 20 mm. The time period should be defined by the bounds of the time coordinate.";
static const char R20MM_LONGNAME_ET[]= "Count of Days with At Least 20mm of Precipitation";
static const char R20MM_UNITS[]      = "No.";
static const char R20MM_UNITS_ET[]   = "days";

static const char R75P_NAME[]        = "moderate_wet_days_wrt_75th_percentile_of_reference_period";
static const char R75P_LONGNAME[]    = "This is the percent of time per time period of wet days (daily sum at least 1 mm / day) where daily precipitation amount of a wet day is above a reference value. The reference value is calculated as the 75th percentile of all wet days of a given 30 year climate reference period. The time period should be defined by the bounds of the time coordinate.";
static const char R75P_UNITS[]       = "Percent";

static const char R75PTOT_NAME[]     = "precipitation_percent_due_to_R75p_days";
static const char R75PTOT_LONGNAME[] = "Percentage of total precipitation amount per time period due to moderate_wet_days_wrt_75th_percentile_of_reference_period. The time period should be defined by the bounds of the time coordinate.";
static const char R75PTOT_UNITS[]    = "Percent";

static const char R90P_NAME[]        = "wet_days_wrt_90th_percentile_of_reference_period";
static const char R90P_LONGNAME[]    = "This is the percent of time per time period of wet days (daily sum at least 1 mm / day) where daily precipitation amount of a wet day is above a reference value. The reference value is calculated as the 90th percentile of all wet days of a given 30 year climate reference period. The time period should be defined by the bounds of the time coordinate.";
static const char R90P_UNITS[]       = "Percent";

static const char R90PTOT_NAME[]     = "precipitation_percent_due_to_R90p_days";
static const char R90PTOT_LONGNAME[] = "Percentage of total precipitation amount per time period due towet_days_wrt_90th_percentile_of_reference_period. The time period should be defined by the bounds of the time coordinate.";
static const char R90PTOT_UNITS[]    = "Percent";

static const char R95P_NAME[]        = "very_wet_days_wrt_95th_percentile_of_reference_period";
static const char R95P_LONGNAME[]    = "This is the percent of time per time period of wet days (daily sum at least 1 mm / day) where daily precipitation amount of a wet day is above a reference value. The reference value is calculated as the 95th percentile of all wet days of a given 30 year climate reference period. The time period should be defined by the bounds of the time coordinate.";
static const char R95P_UNITS[]       = "Percent";

static const char R95PTOT_NAME[]     = "precipitation_percent_due_to_R95p_days";
static const char R95PTOT_LONGNAME[] = "Percentage of total precipitation amount per time period due to very_wet_days_wrt_95th_percentile_of_reference_period. The time period should be defined by the bounds of the time coordinate.";
static const char R95PTOT_UNITS[]    = "Percent";

static const char R99P_NAME[]        = "extremely_wet_days_wrt_99th_percentile_of_reference_period";
static const char R99P_LONGNAME[]    = "This is the percent of time per time period of wet days (daily sum at least 1 mm / day) where daily precipitation amount of a wet day is above a reference value. The reference value is calculated as the 99th percentile of all wet days of a given 30 year climate reference period. The time period should be defined by the bounds of the time coordinate.";
static const char R99P_UNITS[]       = "Percent";

static const char R99PTOT_NAME[]     = "precipitation_percent_due_to_R99p_days";
static const char R99PTOT_LONGNAME[] = "percentage of total  precipitation amount per time period due to extremely_wet_days_wrt_99th_percentile_of_reference_period. The time period should be defined by the bounds of the time coordinate.";
//static const char R99PTOT_UNITS[]    = "Percent";

static const char RR1_NAME[]         = "wet_days_index_per_time_period";
static const char RR1_LONGNAME[]     = "Wet days index is the number of days per time period with daily precipitation of at least %g mm. The time period should be defined by the bounds of the time coordinate.";
static const char RR1_UNITS[]        = "No.";

static const char RX1DAY_NAME[]      = "highest_one_day_precipitation_amount_per_time_period";
static const char RX1DAY_NAME_ET[]   = "rx1dayETCCDI";
static const char RX1DAY_LONGNAME[]  = "Highest one day precipitation is the maximum of one day precipitation amount in a given time period. The time period should be defined by the bounds of the time coordinate.";
static const char RX1DAY_LONGNAME_ET[]= "Maximum 1-day Precipitation";
static const char RX1DAY_UNITS[]     = "mm per day";
static const char RX1DAY_UNITS_ET[]  = "mm";

static const char RX5DAY_NAME[]      = "highest_five_day_precipitation_amount_per_time_period";
static const char RX5DAY_NAME_ET[]   = "rx5dayETCCDI";
static const char RX5DAY_LONGNAME[]  = "Highest precipitation amount for five day interval (including the calendar day as the last day). The time period should be defined by the bounds of the time coordinate.";
static const char RX5DAY_LONGNAME_ET[]= "Maximum Consecutive 5-day Precipitation";
static const char RX5DAY_UNITS[]     = "mm per 5 day";
static const char RX5DAY_UNITS_ET[]  = "mm";
static const char RX5DAY_NAME2[]     = "number_of_5day_heavy_precipitation_periods_per_time_period";
static const char RX5DAY_LONGNAME2[] = "Number of 5day periods in given time period with precipitation amount exceeding %1.0f mm / 5 days. The time period should be defined by the bounds of the time coordinate.";
static const char RX5DAY_UNITS2[]    = "No.";

static const char SDII_NAME[]        = "simple_daily_intensity_index_per_time_period";
static const char SDII_NAME_ET[]     = "sdiiETCCDI";
static const char SDII_LONGNAME[]    = "Simple daily intensity index is the mean of precipitation amount on wet days. A wet day is a day with precipitation sum of at least %g mm. The time period should be defined by the bounds of the time coordinate.";
static const char SDII_LONGNAME_ET[] = "Simple Precipitation Intensity Index";
static const char SDII_UNITS[]       = "mm";
static const char SDII_UNITS_ET[]    = "mm d-1";

static const char FDNS_NAME[]        = "frost_days_where_no_snow_index_per_time_period";
static const char FDNS_LONGNAME[]    = "Frost days where no snow index is the number of days without snowcover and where the minimum of temperature is below 0 degree Celsius. The time period should be defined by the bounds of the time coordinate.";
static const char FDNS_UNITS[]       = "No.";

static const char STRWIN_NAME[]      = "strong_wind_days_index_per_time_period";
static const char STRWIN_LONGNAME[]  = "Strong wind days index is the number of days per time period where maximum wind speed is above %1.0f m/s. The time period should be defined by the bounds of the time coordinate.";
static const char STRWIN_UNITS[]     = "No.";
static const char STRWIN_NAME2[]     = "consecutive_strong_wind_days_index_per_time_period";
static const char STRWIN_LONGNAME2[] = "Greatest number of consecutive strong wind days per time period. The time period should be defined by the bounds of the time coordinate.";
static const char STRWIN_UNITS2[]    = "No.";

static const char STRBRE_NAME[]      = "strong_breeze_days_index_per_time_period";
static const char STRBRE_LONGNAME[]  = "Strong breeze days index is the number of days per time period where maximum wind speed is above 10.5 m/s. The time period should be defined by the bounds of the time coordinate.";
static const char STRBRE_NAME2[]     = "consecutive_strong_breeze_days_index_per_time_period";
static const char STRBRE_LONGNAME2[] = "Greatest number of consecutive strong breeze days per time period. The time period should be defined by the bounds of the time coordinate.";

//static const char STRGAL_NAME[]      = "strong_gale_days_index_per_time_period";
//static const char STRGAL_LONGNAME[]  = "Strong gale days index is the number of days per time period where maximum wind speed is above 20.5 m/s. The time period should be defined by the bounds of the time coordinate.";
//static const char STRGAL_NAME2[]     = "consecutive_strong_gale_days_index_per_time_period";
//static const char STRGAL_LONGNAME2[] = "Greatest number of consecutive strong gale days per time period. The time period should be defined by the bounds of the time coordinate.";

static const char HURR_NAME[]        = "hurricane_days_index_per_time_period";
static const char HURR_LONGNAME[]    = "Hurricane days index is the number of days per time period where maximum wind speed is above 32.5 m/s. The time period should be defined by the bounds of the time coordinate.";
static const char HURR_NAME2[]       = "consecutive_hurricane_days_index_per_time_period";
static const char HURR_LONGNAME2[]   = "Greatest number of consecutive hurricane days per time period. The time period should be defined by the bounds of the time coordinate.";

// clang-format on

/* ECA temperature indices */

static int
addWithFrequency(const std::vector<std::string> &params, const char *operatorName, size_t defaultDays)
{
  int opID = 0;

  KVList kvlist;
  if (kvlist.parse_arguments(1, params) != 0) cdo_abort("Argument parse error!");
  auto kv = kvlist.search("freq");
  if (kv && kv->nvalues > 0)
    {
      // clang-format off
      if      (kv->values[0] == "month") opID = cdo_operator_add(operatorName, 0, CMP_MONTH, nullptr);
      else if (kv->values[0] == "year")  opID = cdo_operator_add(operatorName, 0, CMP_YEAR, nullptr);
      else cdo_abort("Frequency '%s' unknown.", kv->values[0]);
      // clang-format on
    }
  else
    opID = cdo_operator_add(operatorName, 0, defaultDays, nullptr);

  return opID;
}

void *
EcaCfd(void *process)
{
  int ndays = 5;

  cdo_initialize(process);

  if (cdo_operator_argc() > 2) cdo_abort("Too many arguments!");
  if (cdo_operator_argc() > 1)
    {
      auto params = cdo_get_oper_argv();
      params = std::vector<std::string>(params.begin() + 1, params.end());
      addWithFrequency(params, "eca_cfd", CMP_DATE);
    }
  else
    {
      if (cdo_operator_argc() > 0) ndays = parameter_to_int(cdo_operator_argv(0));
      cdo_operator_add("eca_cfd", 0, CMP_DATE, nullptr);
    }

  char cfd_longname2[1024];
  char cfd_name2[1024];
  sprintf(cfd_longname2, CFD_LONGNAME2, ndays);
  sprintf(cfd_name2, CFD_NAME2, ndays);

  ECA_REQUEST_1 request;

  request.var1.name = CFD_NAME;
  request.var1.longname = CFD_LONGNAME;
  request.var1.refdate = ECA_refdate;
  request.var1.f1 = vfarselltc;
  request.var1.f1arg = TO_KELVIN(0.0);
  request.var1.f2 = vfarnum2;
  request.var1.f3 = field2_max;
  request.var2.name = cfd_name2;
  request.var2.longname = cfd_longname2;
  request.var2.units = CFD_UNITS2;
  request.var2.h1 = vfarseleqc;
  request.var2.h1arg = ndays + 1;
  request.var2.h3 = vfarnum;

  eca1(request);

  cdo_finish();

  return 0;
}

void *
EcaCsu(void *process)
{
  double argT = 25.0;
  int ndays = 5;

  cdo_initialize(process);

  if (cdo_operator_argc() > 3) cdo_abort("Too many arguments!");
  if (cdo_operator_argc() > 2)
    {
      auto params = cdo_get_oper_argv();
      params = std::vector<std::string>(params.begin() + 2, params.end());
      addWithFrequency(params, "eca_csu", CMP_DATE);
    }
  else if (cdo_operator_argc() > 0)
    {
      cdo_operator_add("eca_csu", 0, CMP_DATE, nullptr);
      argT = parameter_to_double(cdo_operator_argv(0));
      if (cdo_operator_argc() == 2) ndays = parameter_to_int(cdo_operator_argv(1));
    }
  else
    cdo_operator_add("eca_csu", 0, CMP_DATE, nullptr);

  char csu_longname2[1024];
  char csu_name2[1024];
  sprintf(csu_longname2, CSU_LONGNAME2, ndays);
  sprintf(csu_name2, CSU_NAME2, ndays);

  ECA_REQUEST_1 request;

  request.var1.name = CSU_NAME;
  request.var1.longname = CSU_LONGNAME;
  request.var1.refdate = ECA_refdate;
  request.var1.f1 = vfarselgtc;
  request.var1.f1arg = TO_KELVIN(argT);
  request.var1.f2 = vfarnum2;
  request.var1.f3 = field2_max;
  request.var2.name = csu_name2;
  request.var2.longname = csu_longname2;
  request.var2.units = CSU_UNITS2;
  request.var2.h1 = vfarseleqc;
  request.var2.h1arg = ndays + 1;
  request.var2.h3 = vfarnum;

  eca1(request);

  cdo_finish();

  return 0;
}

void *
EcaCwdi(void *process)
{
  int argN = 6;
  double argT = 5.0;

  cdo_initialize(process);

  if (cdo_operator_argc() > 2)
    {
      auto params = cdo_get_oper_argv();
      params = std::vector<std::string>(params.begin() + 2, params.end());
      addWithFrequency(params, "eca_cwdi", CMP_DATE);
      argT = parameter_to_double(cdo_operator_argv(1));
      argN = parameter_to_int(cdo_operator_argv(0));
    }
  else
    {
      if (cdo_operator_argc() > 1)
        argT = parameter_to_double(cdo_operator_argv(1));
      else if (cdo_operator_argc() > 0)
        argN = parameter_to_int(cdo_operator_argv(0));
      cdo_operator_add("eca_cwdi", 0, CMP_DATE, nullptr);
    }

  char longname[sizeof(CWDI_LONGNAME) + 80];
  sprintf(longname, CWDI_LONGNAME, argN, argT);

  ECA_REQUEST_2 request;

  request.var1.name = CWDI_NAME;
  request.var1.longname = longname;
  request.var1.refdate = ECA_refdate;
  request.var1.units = CWDI_UNITS;
  request.var1.f2 = fieldc_sub;
  request.var1.f2arg = argT;
  request.var1.f3 = vfarsellt;
  request.var1.f4 = vfarnum2;
  request.var1.f5 = vfarnum3;
  request.var1.f5arg = argN;
  request.var2.name = CWDI_NAME2;
  request.var2.longname = CWDI_LONGNAME2;
  request.var2.units = CWDI_UNITS2;
  request.var2.h1 = vfarseleqc;
  request.var2.h1arg = argN;
  request.var2.h2 = vfarnum;

  eca2(request);

  cdo_finish();

  return 0;
}

void *
EcaCwfi(void *process)
{
  int argN = 6;

  cdo_initialize(process);

  int OPID_ECA = 0, OPID_ETC = 0;
  if (cdo_operator_argc() > 1)
    {
      auto params = cdo_get_oper_argv();
      params = std::vector<std::string>(params.begin() + 1, params.end());
      OPID_ECA = addWithFrequency(params, "eca_cwfi", CMP_DATE);
      OPID_ETC = addWithFrequency(params, "etccdi_csdi", CMP_YEAR);
      argN = parameter_to_int(cdo_operator_argv(0));
    }
  else
    {
      if (cdo_operator_argc() > 0) argN = parameter_to_int(cdo_operator_argv(0));
      OPID_ECA = cdo_operator_add("eca_cwfi", 0, CMP_DATE, nullptr);
      OPID_ETC = cdo_operator_add("etccdi_csdi", 0, CMP_YEAR, nullptr);
    }

  char longname[sizeof(CWFI_LONGNAME) + 40];
  sprintf(longname, CWFI_LONGNAME, argN);

  ECA_REQUEST_2 request;

  if (OPID_ECA == cdo_operator_id())
    {
      request.var1.name = CWFI_NAME;
      request.var1.longname = longname;
      request.var1.units = CWFI_UNITS;
      request.var1.refdate = ECA_refdate;
    }
  else if (OPID_ETC == cdo_operator_id())
    {
      request.var1.name = CWFI_NAME_ET;
      request.var1.longname = CWFI_LONGNAME_ET;
      request.var1.units = CWFI_UNITS_ET;
      request.var1.refdate = ETC_refdate;
    }

  request.var1.f3 = vfarsellt;
  request.var1.f4 = vfarnum2;
  request.var1.f5 = vfarnum3;
  request.var1.f5arg = argN;
  request.var2.name = CWFI_NAME2;
  request.var2.longname = CWFI_LONGNAME2;
  request.var2.units = CWFI_UNITS2;
  request.var2.h1 = vfarseleqc;
  request.var2.h1arg = argN;
  request.var2.h2 = vfarnum;

  eca2(request);

  cdo_finish();

  return 0;
}

void *
EcaEtr(void *process)
{
  cdo_initialize(process);

  cdo_operator_add("eca_etr", 0, CMP_DATE, nullptr);

  ECA_REQUEST_3 request;

  request.name = ETR_NAME;
  request.longname = ETR_LONGNAME;
  request.refdate = ECA_refdate;
  request.f1 = field2_max;
  request.f2 = field2_min;
  request.f3 = field2_sub;

  eca3(request);
  cdo_finish();

  return 0;
}

void *
EcaFd(void *process)
{
  cdo_initialize(process);

  int OPID_ECA = 0, OPID_ETC = 0;
  if (cdo_operator_argc() > 0)
    {
      const auto &params = cdo_get_oper_argv();
      OPID_ECA = addWithFrequency(params, "eca_fd", CMP_DATE);
      OPID_ETC = addWithFrequency(params, "etccdi_fd", CMP_YEAR);
    }
  else
    {
      OPID_ECA = cdo_operator_add("eca_fd", 0, CMP_DATE, nullptr);
      OPID_ETC = cdo_operator_add("etccdi_fd", 0, CMP_YEAR, nullptr);
    }

  ECA_REQUEST_1 request;

  if (OPID_ECA == cdo_operator_id())
    {
      request.var1.name = FD_NAME;
      request.var1.longname = FD_LONGNAME;
      request.var1.refdate = ECA_refdate;
    }
  else if (OPID_ETC == cdo_operator_id())
    {
      request.var1.name = FD_NAME_ET;
      request.var1.longname = FD_LONGNAME_ET;
      request.var1.units = FD_UNITS_ET;
      request.var1.refdate = ETC_refdate;
    }

  request.var1.f1 = vfarselltc;
  request.var1.f1arg = TO_KELVIN(0.0);
  request.var1.f2 = vfarnum;

  eca1(request);

  cdo_finish();

  return 0;
}

/*
 * Definition of GSL: (Thermal) Growing Season Length start at the first span
 * of at least 6 (argN) days with T > 5.0°C (argT) in first half of the year
 * and ends at the first span of ar least 6 (argN) days with T < 5.0°C (argT)
 * in the second half.
 * ATTENTION: Year of the northern hemisphere starts in january to
 * december, whereas for the southern hemisphere is goes from july to june!
 * Hence, at least 18 Month of data is needed for computing the gsl of the whole earth.
 */
void *
EcaGsl(void *process)
{
  int argN = 6;
  double argT = 5.0;
  double minLandFraction = 0.5;

  cdo_initialize(process);
  cdo_operator_add("eca_gsl", 0, CMP_YEAR, nullptr);

  if (cdo_operator_argc() > 0) argN = parameter_to_int(cdo_operator_argv(0));
  if (cdo_operator_argc() > 1) argT = parameter_to_double(cdo_operator_argv(1));
  if (cdo_operator_argc() > 2) minLandFraction = parameter_to_double(cdo_operator_argv(2));

  char longname[sizeof(GSL_LONGNAME) + 160];
  sprintf(longname, GSL_LONGNAME, argN, argT, argN, argT);

  ECA_REQUEST_4 request;

  request.name = GSL_NAME;
  request.longname = longname;
  request.units = GSL_UNITS;
  request.name2 = GSL_NAME2;
  request.longname2 = GSL_LONGNAME2;
  request.units2 = GSL_UNITS2;
  request.s1 = vfarselgtc;
  request.s1arg = TO_KELVIN(argT);
  request.s2 = vfarselltc;
  request.s2arg = TO_KELVIN(argT);
  request.s3 = vfarselgec;
  request.s3arg = minLandFraction;
  request.consecutiveDays = argN;

  eca4(request);

  cdo_finish();

  return 0;
}

void *
EcaHd(void *process)
{
  double argX = 17.0;
  double argA = 17.0;

  cdo_initialize(process);

  cdo_operator_add("eca_hd", 0, CMP_DATE, nullptr);

  if (cdo_operator_argc() > 0)
    {
      argX = parameter_to_double(cdo_operator_argv(0));
      argA = argX;
    }
  if (cdo_operator_argc() > 1) argA = parameter_to_double(cdo_operator_argv(1));

  ECA_REQUEST_1 request;

  request.var1.name = HD_NAME;
  request.var1.longname = HD_LONGNAME;
  request.var1.refdate = ECA_refdate;
  request.var1.units = HD_UNITS;
  request.var1.f1 = vfarselltc;
  request.var1.f1arg = TO_KELVIN(argA);
  request.var1.f2 = field2_sum;
  request.var1.mulc = -1.0;
  request.var1.addc = TO_KELVIN(argX);

  eca1(request);

  cdo_finish();

  return 0;
}

void *
EcaHwdi(void *process)
{
  int argN = 6;
  double argT = 5.0;

  cdo_initialize(process);

  if (cdo_operator_argc() > 2)
    {
      auto params = cdo_get_oper_argv();
      params = std::vector<std::string>(params.begin() + 2, params.end());
      addWithFrequency(params, "eca_hwdi", CMP_DATE);
      argN = parameter_to_int(cdo_operator_argv(0));
      argT = parameter_to_double(cdo_operator_argv(1));
    }
  else
    {
      if (cdo_operator_argc() > 0) argN = parameter_to_int(cdo_operator_argv(0));
      if (cdo_operator_argc() > 1) argT = parameter_to_double(cdo_operator_argv(1));
      cdo_operator_add("eca_hwdi", 0, CMP_DATE, nullptr);
    }

  char longname[sizeof(HWDI_LONGNAME) + 80];
  sprintf(longname, HWDI_LONGNAME, argN, argT);

  ECA_REQUEST_2 request;

  request.var1.name = HWDI_NAME;
  request.var1.longname = longname;
  request.var1.refdate = ECA_refdate;
  request.var1.units = HWDI_UNITS;
  request.var1.f2 = fieldc_add;
  request.var1.f2arg = argT;
  request.var1.f3 = vfarselgt;
  request.var1.f4 = vfarnum2;
  request.var1.f5 = vfarnum3;
  request.var1.f5arg = argN;
  request.var2.name = HWDI_NAME2;
  request.var2.longname = HWDI_LONGNAME2;
  request.var2.units = HWDI_UNITS2;
  request.var2.h1 = vfarseleqc;
  request.var2.h1arg = argN;
  request.var2.h2 = vfarnum;

  eca2(request);

  cdo_finish();

  return 0;
}

void *
EcaHwfi(void *process)
{
  int argN = 6;

  cdo_initialize(process);

  int OPID_ECA = 0, OPID_ETC = 0;
  if (cdo_operator_argc() > 1)
    {
      auto params = cdo_get_oper_argv();
      params = std::vector<std::string>(params.begin() + 1, params.end());
      OPID_ECA = addWithFrequency(params, "eca_hwfi", CMP_DATE);
      OPID_ETC = addWithFrequency(params, "etccdi_wsdi", CMP_YEAR);
      argN = parameter_to_int(cdo_operator_argv(0));
    }
  else
    {
      if (cdo_operator_argc() > 0) argN = parameter_to_int(cdo_operator_argv(0));
      OPID_ECA = cdo_operator_add("eca_hwfi", 0, CMP_DATE, nullptr);
      OPID_ETC = cdo_operator_add("etccdi_wsdi", 0, CMP_YEAR, nullptr);
    }

  char longname[sizeof(HWFI_LONGNAME) + 40];
  sprintf(longname, HWFI_LONGNAME, argN);

  ECA_REQUEST_2 request;

  if (OPID_ECA == cdo_operator_id())
    {
      request.var1.name = HWFI_NAME;
      request.var1.longname = longname;
      request.var1.units = HWFI_UNITS;
      request.var1.refdate = ECA_refdate;
    }
  else if (OPID_ETC == cdo_operator_id())
    {
      request.var1.name = HWFI_NAME_ET;
      request.var1.longname = HWFI_LONGNAME_ET;
      request.var1.units = HWFI_UNITS_ET;
      request.var1.refdate = ETC_refdate;
    }

  request.var1.f3 = vfarselgt;
  request.var1.f4 = vfarnum2;
  request.var1.f5 = vfarnum3;
  request.var1.f5arg = argN;
  request.var2.name = HWFI_NAME2;
  request.var2.longname = HWFI_LONGNAME2;
  request.var2.units = HWFI_UNITS2;
  request.var2.h1 = vfarseleqc;
  request.var2.h1arg = argN;
  request.var2.h2 = vfarnum;

  eca2(request);

  cdo_finish();

  return 0;
}

void *
EcaId(void *process)
{
  cdo_initialize(process);

  int OPID_ECA = 0, OPID_ETC = 0;
  if (cdo_operator_argc() > 0)
    {
      const auto &params = cdo_get_oper_argv();
      OPID_ECA = addWithFrequency(params, "eca_id", CMP_DATE);
      OPID_ETC = addWithFrequency(params, "etccdi_id", CMP_YEAR);
    }
  else
    {
      OPID_ECA = cdo_operator_add("eca_id", 0, CMP_DATE, nullptr);
      OPID_ETC = cdo_operator_add("etccdi_id", 0, CMP_YEAR, nullptr);
    }

  ECA_REQUEST_1 request;

  if (OPID_ECA == cdo_operator_id())
    {
      request.var1.name = ID_NAME;
      request.var1.longname = ID_LONGNAME;
      request.var1.units = ID_UNITS;
    }
  else if (OPID_ETC == cdo_operator_id())
    {
      request.var1.name = ID_NAME_ET;
      request.var1.longname = ID_LONGNAME_ET;
      request.var1.units = ID_UNITS_ET;
    }

  request.var1.f1 = vfarselltc;
  request.var1.f1arg = TO_KELVIN(0.0);
  request.var1.f2 = vfarnum;

  eca1(request);

  cdo_finish();

  return 0;
}

void *
EcaSu(void *process)
{
  double argT = 25.0;

  cdo_initialize(process);

  int OPID_ECA = 0, OPID_ETC = 0;
  if (cdo_operator_argc() > 0) argT = parameter_to_double(cdo_operator_argv(0));
  if (cdo_operator_argc() > 1)
    {
      auto params = cdo_get_oper_argv();
      params = std::vector<std::string>(params.begin() + 1, params.end());
      OPID_ECA = addWithFrequency(params, "eca_su", CMP_DATE);
      OPID_ETC = addWithFrequency(params, "etccdi_su", CMP_YEAR);
    }
  else
    {
      OPID_ECA = cdo_operator_add("eca_su", 0, CMP_DATE, nullptr);
      OPID_ETC = cdo_operator_add("etccdi_su", 0, CMP_YEAR, nullptr);
    }

  char longname[sizeof(SU_LONGNAME) + 40];
  sprintf(longname, SU_LONGNAME, argT);

  ECA_REQUEST_1 request;

  if (OPID_ECA == cdo_operator_id())
    {
      request.var1.name = SU_NAME;
      request.var1.longname = longname;
      request.var1.refdate = ECA_refdate;
    }
  else if (OPID_ETC == cdo_operator_id())
    {
      request.var1.name = SU_NAME_ET;
      request.var1.longname = SU_LONGNAME_ET;
      request.var1.units = SU_UNITS_ET;
      request.var1.refdate = ETC_refdate;
    }

  request.var1.f1 = vfarselgtc;
  request.var1.f1arg = TO_KELVIN(argT);
  request.var1.f2 = vfarnum;

  eca1(request);

  cdo_finish();

  return 0;
}

void *
EcaTg10p(void *process)
{
  cdo_initialize(process);

  cdo_operator_add("eca_tg10p", 0, CMP_DATE, nullptr);

  ECA_REQUEST_2 request;

  request.var1.name = TG10P_NAME;
  request.var1.longname = TG10P_LONGNAME;
  request.var1.refdate = ECA_refdate;
  request.var1.units = TG10P_UNITS;
  request.var1.f3 = vfarsellt;
  request.var1.f4 = vfarnum;
  request.var1.epilog = PERCENT_OF_TIME;

  eca2(request);

  cdo_finish();

  return 0;
}

void *
EcaTg90p(void *process)
{
  cdo_initialize(process);

  cdo_operator_add("eca_tg90p", 0, CMP_DATE, nullptr);

  ECA_REQUEST_2 request;

  request.var1.name = TG90P_NAME;
  request.var1.longname = TG90P_LONGNAME;
  request.var1.refdate = ECA_refdate;
  request.var1.units = TG90P_UNITS;
  request.var1.f3 = vfarselgt;
  request.var1.f4 = vfarnum;
  request.var1.epilog = PERCENT_OF_TIME;

  eca2(request);

  cdo_finish();

  return 0;
}

void *
EcaTn10p(void *process)
{
  cdo_initialize(process);

  cdo_operator_add("eca_tn10p", 0, CMP_DATE, nullptr);

  ECA_REQUEST_2 request;

  request.var1.name = TN10P_NAME;
  request.var1.longname = TN10P_LONGNAME;
  request.var1.refdate = ECA_refdate;
  request.var1.units = TN10P_UNITS;
  request.var1.f3 = vfarsellt;
  request.var1.f4 = vfarnum;
  request.var1.epilog = PERCENT_OF_TIME;

  eca2(request);

  cdo_finish();

  return 0;
}

void *
EcaTn90p(void *process)
{
  cdo_initialize(process);

  cdo_operator_add("eca_tn90p", 0, CMP_DATE, nullptr);

  ECA_REQUEST_2 request;

  request.var1.name = TN90P_NAME;
  request.var1.longname = TN90P_LONGNAME;
  request.var1.refdate = ECA_refdate;
  request.var1.units = TN90P_UNITS;
  request.var1.f3 = vfarselgt;
  request.var1.f4 = vfarnum;
  request.var1.epilog = PERCENT_OF_TIME;

  eca2(request);

  cdo_finish();

  return 0;
}

void *
EcaTr(void *process)
{
  double argT = 20.0;

  cdo_initialize(process);

  int OPID_ECA = 0, OPID_ETC = 0;
  if (cdo_operator_argc() > 0) argT = parameter_to_double(cdo_operator_argv(0));
  if (cdo_operator_argc() > 1)
    {
      auto params = cdo_get_oper_argv();
      params = std::vector<std::string>(params.begin() + 1, params.end());
      OPID_ECA = addWithFrequency(params, "eca_tr", CMP_DATE);
      OPID_ETC = addWithFrequency(params, "etccdi_tr", CMP_YEAR);
    }
  else
    {
      OPID_ECA = cdo_operator_add("eca_tr", 0, CMP_DATE, nullptr);
      OPID_ETC = cdo_operator_add("etccdi_tr", 0, CMP_YEAR, nullptr);
    }

  char tr_longname[1024];
  sprintf(tr_longname, TR_LONGNAME, argT);

  ECA_REQUEST_1 request;

  if (OPID_ECA == cdo_operator_id())
    {
      request.var1.name = TR_NAME;
      request.var1.longname = tr_longname;
      request.var1.units = TR_UNITS;
      request.var1.refdate = ECA_refdate;
    }
  else if (OPID_ETC == cdo_operator_id())
    {
      request.var1.name = TR_NAME_ET;
      request.var1.longname = TR_LONGNAME_ET;
      request.var1.units = TR_UNITS_ET;
      request.var1.refdate = ETC_refdate;
    }

  request.var1.f1 = vfarselgtc;
  request.var1.f1arg = TO_KELVIN(argT);
  request.var1.f2 = vfarnum;

  eca1(request);

  cdo_finish();

  return 0;
}

void *
EcaTx10p(void *process)
{
  cdo_initialize(process);

  cdo_operator_add("eca_tx10p", 0, CMP_DATE, nullptr);

  ECA_REQUEST_2 request;

  request.var1.name = TX10P_NAME;
  request.var1.longname = TX10P_LONGNAME;
  request.var1.refdate = ECA_refdate;
  request.var1.units = TX10P_UNITS;
  request.var1.f3 = vfarsellt;
  request.var1.f4 = vfarnum;
  request.var1.epilog = PERCENT_OF_TIME;

  eca2(request);

  cdo_finish();

  return 0;
}

void *
EcaTx90p(void *process)
{
  cdo_initialize(process);
  if (cdo_operator_argc() > 0)
    {
      if ('m' == cdo_operator_argv(0)[0])
        cdo_operator_add("eca_tx90p", 0, CMP_MONTH, nullptr);  // monthly mode
      else
        cdo_warning("Parameter value '%s' is invalid. The only valid value is "
                    "'m' indicating monthly mode. Operating in yearly mode now.",
                    cdo_operator_argv(0));
    }
  else
    cdo_operator_add("eca_tx90p", 0, CMP_DATE, nullptr);

  ECA_REQUEST_2 request;

  request.var1.name = TX90P_NAME;
  request.var1.longname = TX90P_LONGNAME;
  request.var1.refdate = ECA_refdate;
  request.var1.units = TX90P_UNITS;
  request.var1.f3 = vfarselgt;
  request.var1.f4 = vfarnum;
  request.var1.epilog = PERCENT_OF_TIME;

  eca2(request);

  cdo_finish();

  return 0;
}

// ECA precipitation indices

void *
EcaCdd(void *process)
{
  double threshold = 1;
  int ndays = 5;

  cdo_initialize(process);

  int OPID_ECA = 0, OPID_ETC = 0;
  if (cdo_operator_argc() > 3)
    cdo_abort("Too many arguments!");
  else if (cdo_operator_argc() > 2)
    {
      auto params = cdo_get_oper_argv();
      params = std::vector<std::string>(params.begin() + 2, params.end());
      OPID_ECA = addWithFrequency(params, "eca_cdd", CMP_DATE);
      OPID_ETC = addWithFrequency(params, "etccdi_cdd", CMP_YEAR);
    }
  else if (cdo_operator_argc() > 0)
    {
      threshold = parameter_to_double(cdo_operator_argv(0));
      if (cdo_operator_argc() == 2) ndays = parameter_to_int(cdo_operator_argv(1));
      OPID_ECA = cdo_operator_add("eca_cdd", 0, CMP_DATE, nullptr);
      OPID_ETC = cdo_operator_add("etccdi_cdd", 0, CMP_YEAR, nullptr);
    }
  else
    {
      OPID_ECA = cdo_operator_add("eca_cdd", 0, CMP_DATE, nullptr);
      OPID_ETC = cdo_operator_add("etccdi_cdd", 0, CMP_YEAR, nullptr);
    }

  char cdd_longname[1024];
  char cdd_longname2[1024];
  char cdd_name2[1024];
  sprintf(cdd_longname, CDD_LONGNAME, threshold);
  sprintf(cdd_longname2, CDD_LONGNAME2, ndays);
  sprintf(cdd_name2, CDD_NAME2, ndays);

  ECA_REQUEST_1 request;

  if (OPID_ECA == cdo_operator_id())
    {
      request.var1.name = CDD_NAME;
      request.var1.longname = cdd_longname;
      request.var1.units = CDD_UNITS;
      request.var1.refdate = ECA_refdate;
    }
  else if (OPID_ETC == cdo_operator_id())
    {
      request.var1.name = CDD_NAME_ET;
      request.var1.longname = CDD_LONGNAME_ET;
      request.var1.units = CDD_UNITS_ET;
      request.var1.refdate = ETC_refdate;
    }

  request.var1.f1 = vfarselltc;
  request.var1.f1arg = threshold;
  request.var1.f2 = vfarnum2;
  request.var1.f3 = field2_max;
  request.var2.name = cdd_name2;
  request.var2.longname = cdd_longname2;
  request.var2.units = CDD_UNITS2;
  request.var2.h1 = vfarseleqc;
  request.var2.h1arg = ndays + 1;
  request.var2.h3 = vfarnum;

  eca1(request);

  cdo_finish();

  return 0;
}

void *
EcaCwd(void *process)
{
  double threshold = 1;
  int ndays = 5;

  cdo_initialize(process);

  int OPID_ECA = 0, OPID_ETC = 0;

  if (cdo_operator_argc() > 3)
    cdo_abort("Too many arguments!");
  else if (cdo_operator_argc() > 2)
    {
      auto params = cdo_get_oper_argv();
      params = std::vector<std::string>(params.begin() + 2, params.end());
      OPID_ECA = addWithFrequency(params, "eca_cwd", CMP_DATE);
      OPID_ETC = addWithFrequency(params, "etccdi_cwd", CMP_YEAR);
    }
  else if (cdo_operator_argc() > 0)
    {
      threshold = parameter_to_double(cdo_operator_argv(0));
      if (cdo_operator_argc() == 2) ndays = parameter_to_int(cdo_operator_argv(1));
      OPID_ECA = cdo_operator_add("eca_cwd", 0, CMP_DATE, nullptr);
      OPID_ETC = cdo_operator_add("etccdi_cwd", 0, CMP_YEAR, nullptr);
    }
  else
    {
      OPID_ECA = cdo_operator_add("eca_cwd", 0, CMP_DATE, nullptr);
      OPID_ETC = cdo_operator_add("etccdi_cwd", 0, CMP_YEAR, nullptr);
    }

  char cwd_longname[1024];
  char cwd_longname2[1024];
  char cwd_name2[1024];
  sprintf(cwd_longname, CWD_LONGNAME, threshold);
  sprintf(cwd_longname2, CWD_LONGNAME2, ndays);
  sprintf(cwd_name2, CWD_NAME2, ndays);

  ECA_REQUEST_1 request;

  if (OPID_ECA == cdo_operator_id())
    {
      request.var1.name = CWD_NAME;
      request.var1.longname = cwd_longname;
      request.var1.units = CWD_UNITS;
      request.var1.refdate = ECA_refdate;
    }
  else if (OPID_ETC == cdo_operator_id())
    {
      request.var1.name = CWD_NAME_ET;
      request.var1.longname = CWD_LONGNAME_ET;
      request.var1.units = CWD_UNITS_ET;
      request.var1.refdate = ETC_refdate;
    }
  request.var1.f1 = vfarselgec;
  request.var1.f1arg = threshold;
  request.var1.f2 = vfarnum2;
  request.var1.f3 = field2_max;
  request.var2.name = cwd_name2;
  request.var2.longname = cwd_longname2;
  request.var2.units = CWD_UNITS2;
  request.var2.h1 = vfarseleqc;
  request.var2.h1arg = ndays + 1;
  request.var2.h3 = vfarnum;

  eca1(request);

  cdo_finish();

  return 0;
}

void *
EcaPd(void *process)
{
  char lnamebuffer[1024];
  double threshold = 0;
  ECA_REQUEST_1 request;

  cdo_initialize(process);

  int datelenOp = CMP_DATE;
  if (cdo_operator_argc() > 0)
    {
      auto params = cdo_get_oper_argv();
      KVList kvlist;
      if (strstr(cdo_operator_argv(0).c_str(), "=") || cdo_operator_argc() > 1)
        {
          if (cdo_operator_argc() > 1) params = std::vector<std::string>(params.begin() + 1, params.end());
          if (kvlist.parse_arguments(1, params) != 0) cdo_abort("Argument parse error!");
          auto kv = kvlist.search("freq");
          if (kv && kv->nvalues > 0)
            {
              // clang-format off
              if      (kv->values[0] == "month") datelenOp = CMP_MONTH;
              else if (kv->values[0] == "year")  datelenOp = CMP_YEAR;
              // clang-format on
            }
        }
    }

  // clang-format off
  const auto ECA_PD       = cdo_operator_add("eca_pd",       0, datelenOp, nullptr);
  const auto ETCCDI_PD    = cdo_operator_add("etccdi_r1mm",  0, datelenOp, nullptr);
  const auto ECA_R10MM    = cdo_operator_add("eca_r10mm",    0, datelenOp, nullptr);
  const auto ETCCDI_R10MM = cdo_operator_add("etccdi_r10mm", 0, datelenOp, nullptr);
  const auto ECA_R20MM    = cdo_operator_add("eca_r20mm",    0, datelenOp, nullptr);
  const auto ETCCDI_R20MM = cdo_operator_add("etccdi_r20mm", 0, datelenOp, nullptr);
  // clang-format on

  const auto operatorID = cdo_operator_id();

  if (operatorID == ECA_PD || operatorID == ETCCDI_PD)
    {
      if (operatorID == ECA_PD)
        {
          operator_input_arg("daily precipitation amount threshold in [mm]");

          if (cdo_operator_argc() < 1) cdo_abort("Too few arguments!");
          if (cdo_operator_argc() > 2) cdo_abort("Too many arguments!");
          threshold = parameter_to_double(cdo_operator_argv(0));
          sprintf(lnamebuffer, PD_LONGNAME, threshold);
          request.var1.name = PD_NAME;
          request.var1.longname = lnamebuffer;
          request.var1.units = PD_UNITS;
        }
      else
        {
          threshold = 1;
          request.var1.name = PD_NAME_ET;
          request.var1.longname = PD_LONGNAME_ET;
          request.var1.units = PD_UNITS_ET;
        }
      if (threshold < 0) cdo_abort("Parameter out of range: threshold = %g", threshold);
    }
  else if (operatorID == ECA_R10MM || operatorID == ETCCDI_R10MM)
    {
      threshold = 10;
      if (operatorID == ECA_R10MM)
        {
          request.var1.name = R10MM_NAME;
          request.var1.longname = R10MM_LONGNAME;
          request.var1.units = R10MM_UNITS;
          request.var1.refdate = ECA_refdate;
        }
      else
        {
          request.var1.name = R10MM_NAME_ET;
          request.var1.longname = R10MM_LONGNAME_ET;
          request.var1.units = R10MM_UNITS_ET;
          request.var1.refdate = ETC_refdate;
        }
    }
  else if (operatorID == ECA_R20MM || operatorID == ETCCDI_R20MM)
    {
      threshold = 20;
      if (operatorID == ECA_R20MM)
        {
          request.var1.name = R20MM_NAME;
          request.var1.longname = R20MM_LONGNAME;
          request.var1.units = R20MM_UNITS;
          request.var1.refdate = ECA_refdate;
        }
      else
        {
          request.var1.name = R20MM_NAME_ET;
          request.var1.longname = R20MM_LONGNAME_ET;
          request.var1.units = R20MM_UNITS_ET;
          request.var1.refdate = ETC_refdate;
        }
    }

  if (Options::cdoVerbose) cdo_print("threshold = %g", threshold);

  request.var1.f1 = vfarselgec;
  request.var1.f1arg = threshold;
  request.var1.f2 = vfarnum;

  eca1(request);

  cdo_finish();

  return 0;
}

void *
EcaR75p(void *process)
{
  cdo_initialize(process);

  cdo_operator_add("eca_r75p", 0, CMP_DATE, nullptr);

  ECA_REQUEST_2 request;

  request.var1.name = R75P_NAME;
  request.var1.longname = R75P_LONGNAME;
  request.var1.refdate = ECA_refdate;
  request.var1.units = R75P_UNITS;
  request.var1.f1 = vfarselgec;
  request.var1.f3 = vfarselgt;
  request.var1.f4 = vfarnum;
  request.var1.epilog = PERCENT_OF_TIME;

  eca2(request);

  cdo_finish();

  return 0;
}

void *
EcaR75ptot(void *process)
{
  cdo_initialize(process);

  cdo_operator_add("eca_r75ptot", 0, CMP_DATE, nullptr);

  ECA_REQUEST_2 request;

  request.var1.name = R75PTOT_NAME;
  request.var1.longname = R75PTOT_LONGNAME;
  request.var1.refdate = ECA_refdate;
  request.var1.units = R75PTOT_UNITS;
  request.var1.f1 = vfarselgec;
  request.var1.f3 = vfarselgt;
  request.var1.f4 = field2_sum;
  request.var1.epilog = PERCENT_OF_TOTAL_AMOUNT;

  eca2(request);

  cdo_finish();

  return 0;
}

void *
EcaR90p(void *process)
{
  cdo_initialize(process);

  cdo_operator_add("eca_r90p", 0, CMP_DATE, nullptr);

  ECA_REQUEST_2 request;

  request.var1.name = R90P_NAME;
  request.var1.longname = R90P_LONGNAME;
  request.var1.refdate = ECA_refdate;
  request.var1.units = R90P_UNITS;
  request.var1.f1 = vfarselgec;
  request.var1.f3 = vfarselgt;
  request.var1.f4 = vfarnum;
  request.var1.epilog = PERCENT_OF_TIME;

  eca2(request);

  cdo_finish();

  return 0;
}

void *
EcaR90ptot(void *process)
{
  cdo_initialize(process);

  cdo_operator_add("eca_r90ptot", 0, CMP_DATE, nullptr);

  ECA_REQUEST_2 request;

  request.var1.name = R90PTOT_NAME;
  request.var1.longname = R90PTOT_LONGNAME;
  request.var1.refdate = ECA_refdate;
  request.var1.units = R90PTOT_UNITS;
  request.var1.f1 = vfarselgec;
  request.var1.f3 = vfarselgt;
  request.var1.f4 = field2_sum;
  request.var1.epilog = PERCENT_OF_TOTAL_AMOUNT;

  eca2(request);

  cdo_finish();

  return 0;
}

void *
EcaR95p(void *process)
{
  cdo_initialize(process);

  cdo_operator_add("eca_r95p", 0, CMP_DATE, nullptr);

  ECA_REQUEST_2 request;

  request.var1.name = R95P_NAME;
  request.var1.longname = R95P_LONGNAME;
  request.var1.refdate = ECA_refdate;
  request.var1.units = R95P_UNITS;
  request.var1.f1 = vfarselgec;
  request.var1.f3 = vfarselgt;
  request.var1.f4 = vfarnum;
  request.var1.epilog = PERCENT_OF_TIME;

  eca2(request);

  cdo_finish();

  return 0;
}

void *
EcaR95ptot(void *process)
{
  cdo_initialize(process);

  cdo_operator_add("eca_r95ptot", 0, CMP_DATE, nullptr);

  ECA_REQUEST_2 request;

  request.var1.name = R95PTOT_NAME;
  request.var1.longname = R95PTOT_LONGNAME;
  request.var1.refdate = ECA_refdate;
  request.var1.units = R95PTOT_UNITS;
  request.var1.f1 = vfarselgec;
  request.var1.f3 = vfarselgt;
  request.var1.f4 = field2_sum;
  request.var1.epilog = PERCENT_OF_TOTAL_AMOUNT;

  eca2(request);

  cdo_finish();

  return 0;
}

void *
EcaR99p(void *process)
{
  cdo_initialize(process);

  cdo_operator_add("eca_r99p", 0, CMP_DATE, nullptr);

  ECA_REQUEST_2 request;

  request.var1.name = R99P_NAME;
  request.var1.longname = R99P_LONGNAME;
  request.var1.refdate = ECA_refdate;
  request.var1.units = R99P_UNITS;
  request.var1.f1 = vfarselgec;
  request.var1.f3 = vfarselgt;
  request.var1.f4 = vfarnum;
  request.var1.epilog = PERCENT_OF_TIME;

  eca2(request);

  cdo_finish();

  return 0;
}

void *
EcaR99ptot(void *process)
{
  cdo_initialize(process);

  cdo_operator_add("eca_r99ptot", 0, CMP_DATE, nullptr);

  ECA_REQUEST_2 request;

  request.var1.name = R99PTOT_NAME;
  request.var1.longname = R99PTOT_LONGNAME;
  request.var1.refdate = ECA_refdate;
  request.var1.f1 = vfarselgec;
  request.var1.f3 = vfarselgt;
  request.var1.f4 = field2_sum;
  request.var1.epilog = PERCENT_OF_TOTAL_AMOUNT;

  eca2(request);

  cdo_finish();

  return 0;
}

void *
EcaRr1(void *process)
{
  double threshold = 1;

  cdo_initialize(process);

  if (cdo_operator_argc() > 2)
    cdo_abort("Too many arguments!");
  else if (cdo_operator_argc() > 1)
    {
      auto params = cdo_get_oper_argv();
      params = std::vector<std::string>(params.begin() + 1, params.end());
      addWithFrequency(params, "eca_rr1", CMP_DATE);
    }
  else
    {
      if (cdo_operator_argc() == 1) threshold = parameter_to_double(cdo_operator_argv(0));
      cdo_operator_add("eca_rr1", 0, CMP_DATE, nullptr);
    }

  char lnamebuffer[1024];
  sprintf(lnamebuffer, RR1_LONGNAME, threshold);

  ECA_REQUEST_1 request;

  request.var1.name = RR1_NAME;
  request.var1.longname = lnamebuffer;
  request.var1.units = RR1_UNITS;
  request.var1.f1 = vfarselgec;
  request.var1.f1arg = threshold;
  request.var1.f2 = vfarnum;

  eca1(request);

  cdo_finish();

  return 0;
}

void *
EcaRx1day(void *process)
{
  cdo_initialize(process);

  int OPID_ECA = 0, OPID_ETC = 0;

  if (cdo_operator_argc() > 0)
    {
      const auto &params = cdo_get_oper_argv();
      OPID_ECA = addWithFrequency(params, "eca_rx1day", CMP_DATE);
      OPID_ETC = addWithFrequency(params, "etccdi_rx1day", CMP_YEAR);
    }
  else
    {
      OPID_ECA = cdo_operator_add("eca_rx1day", 0, CMP_DATE, nullptr);
      OPID_ETC = cdo_operator_add("etccdi_rx1day", 0, CMP_YEAR, nullptr);
    }

  ECA_REQUEST_1 request;

  if (OPID_ECA == cdo_operator_id())
    {
      request.var1.name = RX1DAY_NAME;
      request.var1.longname = RX1DAY_LONGNAME;
      request.var1.units = RX1DAY_UNITS;
      request.var1.refdate = ECA_refdate;
    }
  else if (OPID_ETC == cdo_operator_id())
    {
      request.var1.name = RX1DAY_NAME_ET;
      request.var1.longname = RX1DAY_LONGNAME_ET;
      request.var1.units = RX1DAY_UNITS_ET;
      request.var1.refdate = ETC_refdate;
    }
  request.var1.f2 = field2_max;

  eca1(request);

  cdo_finish();

  return 0;
}

void *
EcaRx5day(void *process)
{
  double argX = 50.0;

  cdo_initialize(process);

  int OPID_ECA = 0, OPID_ETC = 0;
  if (cdo_operator_argc() > 0)
    {
      argX = parameter_to_double(cdo_operator_argv(0));
      if (cdo_operator_argc() > 1)
        {
          auto params = cdo_get_oper_argv();
          params = std::vector<std::string>(params.begin() + 1, params.end());
          OPID_ECA = addWithFrequency(params, "eca_rx5day", CMP_DATE);
          OPID_ETC = addWithFrequency(params, "etccdi_rx5day", CMP_YEAR);
        }
      else
        {
          OPID_ECA = cdo_operator_add("eca_rx5day", 0, CMP_DATE, nullptr);
          OPID_ETC = cdo_operator_add("etccdi_rx5day", 0, CMP_YEAR, nullptr);
        }
    }
  else
    {
      OPID_ECA = cdo_operator_add("eca_rx5day", 0, CMP_DATE, nullptr);
      OPID_ETC = cdo_operator_add("etccdi_rx5day", 0, CMP_YEAR, nullptr);
    }

  char longname[sizeof(RX5DAY_LONGNAME2) + 40];
  sprintf(longname, RX5DAY_LONGNAME2, argX);

  ECA_REQUEST_1 request;

  if (OPID_ECA == cdo_operator_id())
    {
      request.var1.name = RX5DAY_NAME;
      request.var1.longname = RX5DAY_LONGNAME;
      request.var1.units = RX5DAY_UNITS;
      request.var1.refdate = ECA_refdate;
    }
  else if (OPID_ETC == cdo_operator_id())
    {
      request.var1.name = RX5DAY_NAME_ET;
      request.var1.longname = RX5DAY_LONGNAME_ET;
      request.var1.units = RX5DAY_UNITS_ET;
      request.var1.refdate = ETC_refdate;
    }
  request.var1.f2 = field2_max;
  request.var2.name = RX5DAY_NAME2;
  request.var2.longname = longname;
  request.var2.units = RX5DAY_UNITS2;
  request.var2.h1 = vfarselgec;
  request.var2.h1arg = argX;
  request.var2.h2 = vfarnum;

  eca1(request);

  cdo_finish();

  return 0;
}

void *
EcaSdii(void *process)
{
  ECA_REQUEST_1 request;
  char lnamebuffer[1024];
  double threshold = 1;

  cdo_initialize(process);
  int OPID_ECA = 0, OPID_ETC = 0;

  if (cdo_operator_argc() > 2)
    cdo_abort("Too many arguments!");
  else if (cdo_operator_argc() > 1)
    {
      auto params = cdo_get_oper_argv();
      params = std::vector<std::string>(params.begin() + 1, params.end());
      OPID_ECA = addWithFrequency(params, "eca_sdii", CMP_DATE);
      OPID_ETC = addWithFrequency(params, "etccdi_sdii", CMP_YEAR);
    }
  else
    {
      OPID_ECA = cdo_operator_add("eca_sdii", 0, CMP_DATE, nullptr);
      OPID_ETC = cdo_operator_add("etccdi_sdii", 0, CMP_DATE, nullptr);
      if (cdo_operator_argc() == 1) threshold = parameter_to_double(cdo_operator_argv(0));
    }

  sprintf(lnamebuffer, SDII_LONGNAME, threshold);
  if (OPID_ECA == cdo_operator_id())
    {
      request.var1.name = SDII_NAME;
      request.var1.longname = lnamebuffer;
      request.var1.units = SDII_UNITS;
    }
  else if (OPID_ETC == cdo_operator_id())
    {
      request.var1.name = SDII_NAME_ET;
      request.var1.longname = SDII_LONGNAME_ET;
      request.var1.units = SDII_UNITS_ET;
    }

  request.var1.f1 = vfarselgec;
  request.var1.f1arg = threshold;
  request.var1.f2 = field2_sum;
  request.var1.epilog = MEAN;

  eca1(request);
  cdo_finish();

  return 0;
}

void *
Fdns(void *process)
{
  ECA_REQUEST_2 request;

  cdo_initialize(process);
  cdo_operator_add("fdns", 0, CMP_DATE, nullptr);

  request.var1.name = FDNS_NAME;
  request.var1.longname = FDNS_LONGNAME;
  request.var1.refdate = ECA_refdate;
  request.var1.units = FDNS_UNITS;
  request.var1.f1 = vfarsellec;
  request.var1.f1arg = TO_KELVIN(0.0);
  request.var1.f2 = vfarsellec;
  request.var1.f2arg = 0.01;
  request.var1.f3 = field2_add;  // any f with f(a, b) = miss, if a = miss or b = miss will do here
  request.var1.f4 = vfarnum;

  eca2(request);

  cdo_finish();

  return 0;
}

void *
Strwin(void *process)
{
  double maxWind = 10.5;

  cdo_initialize(process);

  if (cdo_operator_argc() > 2)
    cdo_abort("Too many arguments!");
  else if (cdo_operator_argc() > 1)
    {
      auto params = cdo_get_oper_argv();
      params = std::vector<std::string>(params.begin() + 1, params.end());
      addWithFrequency(params, "strwin", CMP_DATE);
    }
  else
    {
      if (cdo_operator_argc() > 0) maxWind = parameter_to_double(cdo_operator_argv(0));
      cdo_operator_add("strwin", 0, CMP_DATE, nullptr);
    }

  char longname[sizeof(STRWIN_LONGNAME) + 40];
  sprintf(longname, STRWIN_LONGNAME, maxWind);

  ECA_REQUEST_1 request;

  request.var1.name = STRWIN_NAME;
  request.var1.longname = longname;
  request.var1.refdate = ECA_refdate;
  request.var1.units = STRWIN_UNITS;
  request.var1.f1 = vfarselgec;
  request.var1.f1arg = maxWind;
  request.var1.f2 = vfarnum;
  request.var2.name = STRWIN_NAME2;
  request.var2.longname = STRWIN_LONGNAME2;
  request.var2.units = STRWIN_UNITS2;
  request.var2.h1 = vfarselgec;
  request.var2.h1arg = maxWind;
  request.var2.h2 = vfarnum2;
  request.var2.h3 = field2_max;

  eca1(request);

  cdo_finish();

  return 0;
}

void *
Strbre(void *process)
{
  static const double maxWind = 10.5;
  ECA_REQUEST_1 request;

  cdo_initialize(process);
  if (cdo_operator_argc() > 0)
    {
      const auto &params = cdo_get_oper_argv();
      addWithFrequency(params, "strbre", CMP_DATE);
    }
  else
    cdo_operator_add("strbre", 0, CMP_DATE, nullptr);

  request.var1.name = STRBRE_NAME;
  request.var1.longname = STRBRE_LONGNAME;
  request.var1.refdate = ECA_refdate;
  request.var1.units = STRWIN_UNITS;
  request.var1.f1 = vfarselgec;
  request.var1.f1arg = maxWind;
  request.var1.f2 = vfarnum;
  request.var2.name = STRBRE_NAME2;
  request.var2.longname = STRBRE_LONGNAME2;
  request.var2.units = STRWIN_UNITS2;
  request.var2.h1 = vfarselgec;
  request.var2.h1arg = maxWind;
  request.var2.h2 = vfarnum2;
  request.var2.h3 = field2_max;

  eca1(request);
  cdo_finish();

  return 0;
}

void *
Strgal(void *process)
{
  static const double maxWind = 20.5;
  ECA_REQUEST_1 request;

  cdo_initialize(process);
  if (cdo_operator_argc() > 0)
    {
      const auto &params = cdo_get_oper_argv();
      addWithFrequency(params, "strgal", CMP_DATE);
    }
  else
    cdo_operator_add("strgal", 0, CMP_DATE, nullptr);

  request.var1.name = STRBRE_NAME;
  request.var1.longname = STRBRE_LONGNAME;
  request.var1.refdate = ECA_refdate;
  request.var1.units = STRWIN_UNITS;
  request.var1.f1 = vfarselgec;
  request.var1.f1arg = maxWind;
  request.var1.f2 = vfarnum;
  request.var2.name = STRBRE_NAME2;
  request.var2.longname = STRBRE_LONGNAME2;
  request.var2.units = STRWIN_UNITS2;
  request.var2.h1 = vfarselgec;
  request.var2.h1arg = maxWind;
  request.var2.h2 = vfarnum2;
  request.var2.h3 = field2_max;

  eca1(request);
  cdo_finish();

  return 0;
}

void *
Hurr(void *process)
{
  static const double maxWind = 32.5;
  ECA_REQUEST_1 request;

  cdo_initialize(process);
  if (cdo_operator_argc() > 0)
    {
      const auto &params = cdo_get_oper_argv();
      addWithFrequency(params, "hurr", CMP_DATE);
    }
  else
    cdo_operator_add("hurr", 0, CMP_DATE, nullptr);

  request.var1.name = HURR_NAME;
  request.var1.longname = HURR_LONGNAME;
  request.var1.refdate = ECA_refdate;
  request.var1.units = STRWIN_UNITS;
  request.var1.f1 = vfarselgec;
  request.var1.f1arg = maxWind;
  request.var1.f2 = vfarnum;
  request.var2.name = HURR_NAME2;
  request.var2.longname = HURR_LONGNAME2;
  request.var2.units = STRWIN_UNITS2;
  request.var2.h1 = vfarselgec;
  request.var2.h1arg = maxWind;
  request.var2.h2 = vfarnum2;
  request.var2.h3 = field2_max;

  eca1(request);

  cdo_finish();

  return 0;
}
