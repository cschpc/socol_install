/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Modali Kameswarrao

*/

#ifdef HAVE_CONFIG_H
#include "config.h" /* HAVE_LIBMAGICS */
#endif

#include <climits>
#include <cdi.h>
#include <cctype>

#include "process_int.h"
#include <mpim_grid.h>

#ifdef HAVE_LIBMAGICS

#include "magics_api.h"
#include "magics_template_parser.h"
#include "util_string.h"
#include "cdo_vlist.h"
#include "printinfo.h"

#define DBG 0

const char *line_colours[] = {
  "red",
  "green",
  "blue",
  "yellow",
  "cyan",
  "magenta",
  "avocado",
  "beige",
  "brick",
  "brown",
  "burgundy",
  "charcoal",
  "chestnut",
  "coral",
  "cream",
  "evergreen",
  "gold",
  "khaki",
  "kellygreen",
  "lavender",
  "mustard",
  "navy",
  "ochre",
  "olive",
  "peach",
  "pink",
  "rose",
  "rust",
  "sky",
  "tan",
  "tangerine",
  "turquoise",
  "violet",
  "reddishpurple",
  "purplered",
  "purplishred",
  "orangishred",
  "redorange",
  "reddishorange",
  "orange",
  "yellowishorange",
  "orangeyellow",
  "orangishyellow",
  "greenishyellow",
  "yellowgreen",
  "yellowishgreen",
  "bluishgreen",
  "bluegreen",
  "greenishblue",
  "purplishblue",
  "bluepurple",
  "bluishpurple",
  "purple",
};

const char *graph_params[] = { "ymin", "ymax", "sigma", "stat", "obsv", "device", "linewidth" };

int graph_param_count = sizeof(graph_params) / sizeof(char *);
int num_colours = sizeof(line_colours) / sizeof(char *);

int checkdevice(char *device_in);

extern const char *DEVICE;

static int
compareDate(CdiDate date1, CdiDate date2)
{
  int c1[3], c2[3];
  cdiDate_decode(date1, c1, c1 + 1, c1 + 2);
  cdiDate_decode(date2, c2, c2 + 1, c2 + 2);

  for (int i = 0; i < 3; ++i)
    {
      auto flag = c1[i] - c2[i];
      if (flag > 0)
        return 1;
      else if (flag < 0)
        return -1;
    }

  return 0;
}

static int
compareTime(CdiTime time1, CdiTime time2)
{
  int c1[4], c2[4];
  cdiTime_decode(time1, c1, c1 + 1, c1 + 2, c1 + 3);
  cdiTime_decode(time2, c2, c2 + 1, c2 + 2, c1 + 3);

  for (int i = 0; i < 3; ++i)
    {
      auto flag = c1[i] - c2[i];
      if (flag > 0)
        return 1;
      else if (flag < 0)
        return -1;
    }

  return 0;
}

static void
maggraph(const char *plotfile, const std::string &varname, const std::string &varunits, long nfiles, std::vector<long> nts,
         std::vector<std::vector<CdiDateTime>> vDateTimes, std::vector<std::vector<double>> datatab, int nparam,
         std::vector<std::string> &params)
{
  char min_date_time_str[1024], max_date_time_str[1024];
  int min_index = 0, max_index = 0;
  char legend_text_data[256];
  int num_sigma = 2;
  bool stat = false, obsv = false;
  int file_begin = 0;
  int count;
  int ret;
  long tsID, fileID, i, ntime_steps = 0;
  constexpr double MinVal = -1.0e+200;
  constexpr double MaxVal = 1.0e+200;
  double min_val = MaxVal, max_val = MinVal;
  double y_min_val = MaxVal, y_max_val = MinVal;
  int linewidth_val = 8;

  if (DBG)
    {
      fprintf(stderr, "Num params %d\n", nparam);
      for (i = 0; i < nparam; ++i) fprintf(stderr, "Param %s\n", params[i].c_str());
    }

  std::string temp_str;
  for (i = 0; i < nparam; ++i)
    {
      auto splitStrings = split_with_seperator(params[i], '=');
      const auto &key = splitStrings[0];
      const auto &value = splitStrings[1];

      if (key == "obsv")
        {
          temp_str = string_to_lower(value);
          if (temp_str == "true")
            {
              obsv = true;
              file_begin = 1;
              if (DBG) fprintf(stderr, "OBSV true\n");
            }
        }

      if (key == "stat")
        {
          temp_str = string_to_lower(value);
          if (temp_str == "true")
            {
              stat = true;
              if (DBG) fprintf(stderr, "STAT true\n");
            }
        }

      if (key == "ymin")
        {
          y_min_val = std::stod(value);
          if (DBG) fprintf(stderr, "Y min Val %g\n", y_min_val);
        }

      if (key == "ymax")
        {
          y_max_val = std::stod(value);
          if (DBG) fprintf(stderr, "Y max Val %g\n", y_max_val);
        }

      if (key == "linewidth")
        {
          linewidth_val = std::stoi(value);
          if (DBG) fprintf(stderr, "linewidth Val %d\n", linewidth_val);
        }

      if (key == "sigma")
        {
          num_sigma = std::stod(value);
          if (DBG) fprintf(stderr, "SIGMA %d\n", num_sigma);
        }

      if (key == "device")
        {
          temp_str = string_to_upper(value);
          DEVICE = strdup(temp_str.c_str());
          if (DBG) fprintf(stderr, "DEVICE %s\n", DEVICE);

          mag_setc("output_format", DEVICE);
        }
    }

  if (DBG)
    {
      ntime_steps = nts[0];
      fprintf(stderr, " nfiles=%ld  ntime_steps=%ld\n", nfiles, ntime_steps);
      fprintf(stderr, "STAT  %d\n", (int) stat);
    }

  if (stat)
    {
      ntime_steps = nts[0];

      for (fileID = 1; fileID < nfiles; ++fileID)
        {
          if (nts[fileID] != ntime_steps)
            {
              cdo_warning("  Unequal number of time steps! Statistics disabled.");
              stat = false;
              break;
            }

          // First date & time of the present file
          if (compareDate(vDateTimes[0][0].date, vDateTimes[fileID][0].date))
            {
              cdo_warning("  Incosistent start date! Statistics disabled.");
              stat = false;
              break;
            }

          // First time of the present file
          if (compareTime(vDateTimes[0][0].time, vDateTimes[fileID][0].time))
            {
              cdo_warning("  Incosistent start time! Statistics disabled.");
              stat = false;
              break;
            }

          // Last date of the present file
          if (compareDate(vDateTimes[fileID][nts[fileID] - 1].date, vDateTimes[0][nts[0] - 1].date))
            {
              cdo_warning("  Incosistent end date! Statistics disabled.");
              stat = false;
              break;
            }

          // Last time of the present file
          if (compareTime(vDateTimes[fileID][nts[fileID] - 1].time, vDateTimes[0][nts[0] - 1].time))
            {
              cdo_warning("  Incosistent end time! Statistics disabled.");
              stat = false;
              break;
            }
        }
    }

  if (DBG) fprintf(stderr, "STAT  %d\n", (int) stat);

  char ***date_time_str = (char ***) malloc(nfiles * sizeof(char **));

  std::vector<double> date_time;
  std::vector<double> mean_val, std_dev_val;
  std::vector<double> spread_min, spread_max;

  if (stat)
    {
      // if all files are of same number of steps, only one date_time_str array is being used
      date_time_str[0] = (char **) malloc(ntime_steps * sizeof(char *));

      date_time.resize(ntime_steps);
      mean_val.resize(ntime_steps);
      std_dev_val.resize(ntime_steps);
      spread_min.resize(ntime_steps);
      spread_max.resize(ntime_steps);

      for (tsID = 0; tsID < ntime_steps; ++tsID)
        {
          date_time[tsID] = tsID + 1;
          date_time_str[0][tsID] = (char *) malloc(256);
          sprintf(date_time_str[0][tsID], "%s", datetime_to_string(vDateTimes[0][tsID]).c_str());
          mean_val[tsID] = 0.;
          std_dev_val[tsID] = 0.;

          if (DBG) fprintf(stderr, "tsID=%ld: %s\n", tsID, date_time_str[0][tsID]);

          for (fileID = 0; fileID < nfiles; ++fileID)
            {
              if (DBG) fprintf(stderr, "fileID=%ld\n", fileID);

              if (datatab[fileID][tsID] < min_val) min_val = datatab[fileID][tsID];
              if (datatab[fileID][tsID] > max_val) max_val = datatab[fileID][tsID];

              mean_val[tsID] += datatab[fileID][tsID];
              std_dev_val[tsID] = 0.;
              spread_min[tsID] = 0.;
              spread_max[tsID] = 0.;

              if (DBG)
                {
                  fprintf(stderr, " %6g", datatab[fileID][tsID]);
                  fprintf(stderr, "\n");
                }
            }
        }

      for (tsID = 0; tsID < ntime_steps; ++tsID)
        {
          mean_val[tsID] /= (double) nfiles;
          spread_min[tsID] = mean_val[tsID];
          spread_max[tsID] = mean_val[tsID];

          for (fileID = 0; fileID < nfiles; ++fileID)
            {
              std_dev_val[tsID] += (datatab[fileID][tsID] - mean_val[tsID]) * (datatab[fileID][tsID] - mean_val[tsID]);
            }
          std_dev_val[tsID] /= (double) nfiles;
          std_dev_val[tsID] = std::pow(std_dev_val[tsID], 0.5);

          if (DBG) fprintf(stderr, " Mean : %g Std Dev: %g\n", mean_val[tsID], std_dev_val[tsID]);

          spread_min[tsID] = mean_val[tsID] - num_sigma * std_dev_val[tsID];
          spread_max[tsID] = mean_val[tsID] + num_sigma * std_dev_val[tsID];

          if (DBG) fprintf(stderr, " Min : %g Max: %g\n", spread_min[tsID], spread_max[tsID]);
        }

      for (tsID = 0; tsID < ntime_steps; ++tsID)
        {
          if (spread_min[tsID] < min_val) min_val = spread_min[tsID];
          if (spread_max[tsID] > max_val) max_val = spread_max[tsID];
        }

      if (DBG)
        {
          fprintf(stderr, " %6g %6g\n", min_val, max_val);
          fprintf(stderr, " %s %s\n", date_time_str[0][0], date_time_str[0][ntime_steps - 1]);
          fprintf(stderr, "\n");
        }

      strcpy(min_date_time_str, date_time_str[0][0]);
      strcpy(max_date_time_str, date_time_str[0][ntime_steps - 1]);
    }
  else
    {
      /* Find the min_date_time_str from the min's of nfiles
         Find the max_date_time_str from the max's of nfiles
         Construct the date_time_str array
      */

      if (DBG) fprintf(stderr, "STAT  %d\n", (int) stat);

      for (fileID = 0; fileID < nfiles; ++fileID)
        {
          if (DBG) fprintf(stderr, "FILE  %ld\n", fileID);
          date_time.resize(nts[fileID]);
          date_time_str[fileID] = (char **) malloc(nts[fileID] * sizeof(char *));

          for (tsID = 0; tsID < nts[fileID]; ++tsID)
            {
              date_time[tsID] = tsID + 1;

              date_time_str[fileID][tsID] = (char *) malloc(256);
              sprintf(date_time_str[fileID][tsID], "%s", datetime_to_string(vDateTimes[fileID][tsID]).c_str());
              if (DBG && (tsID == 0 || tsID == nts[fileID] - 1)) fprintf(stderr, "%s\n", date_time_str[fileID][tsID]);

              if (datatab[fileID][tsID] < min_val) min_val = datatab[fileID][tsID];
              if (datatab[fileID][tsID] > max_val) max_val = datatab[fileID][tsID];
            }

          if (fileID == 0)
            {
              if (DBG) fprintf(stderr, "\n %s %s\n", date_time_str[fileID][0], date_time_str[fileID][nts[0] - 1]);
              min_index = 0;
              max_index = 0;
            }
          else
            {
              ret = compareDate(vDateTimes[min_index][0].date, vDateTimes[fileID][0].date);
              if (ret == 1)
                min_index = fileID;
              else if (!ret)
                {
                  ret = compareTime(vDateTimes[min_index][0].time, vDateTimes[fileID][0].time);
                  if (ret == -999)
                    cdo_abort("Error in input Date Time");
                  else if (ret == 1)
                    min_index = fileID;
                }
              if (DBG) fprintf(stderr, "Min File ID %d\n", min_index);

              if (DBG) fprintf(stderr, "compareDateOrTime  %s\n", date_time_str[fileID][nts[fileID] - 1]);

              ret = compareDate(vDateTimes[max_index][nts[max_index] - 1].date, vDateTimes[fileID][nts[fileID] - 1].date);
              if (ret == -1)
                max_index = fileID;
              else if (!ret)
                {
                  ret = compareTime(vDateTimes[max_index][nts[max_index] - 1].time, vDateTimes[fileID][nts[fileID] - 1].time);
                  if (ret == -999)
                    cdo_abort("Error in input Date Time");
                  else if (ret == -1)
                    max_index = fileID;
                }

              if (DBG) fprintf(stderr, "Max File ID %d\n", max_index);
            }
        }

      strcpy(min_date_time_str, date_time_str[min_index][0]);
      strcpy(max_date_time_str, date_time_str[max_index][nts[max_index] - 1]);
      if (DBG) fprintf(stderr, "%s %s\n", min_date_time_str, max_date_time_str);
    }

  if (DBG) fprintf(stderr, "%s %s\n", min_date_time_str, max_date_time_str);

  auto splitStrings = split_with_seperator(max_date_time_str, '-');

  auto num_years = std::stoi(splitStrings[0]);
  auto num_months = std::stoi(splitStrings[1]);
  auto num_days = std::stoi(splitStrings[2]);

  splitStrings = split_with_seperator(min_date_time_str, '-');
  num_years -= std::stoi(splitStrings[0]);

  if (num_years <= 1)
    {
      if (num_years == 1)
        num_months += (12 - std::stoi(splitStrings[1]));
      else
        num_months -= (std::stoi(splitStrings[1]));

      if (!num_months)
        num_days -= std::stoi(splitStrings[2]);
      else if (num_months == 1)
        num_days += (31 - std::stoi(splitStrings[2]));
    }

  if (DBG) fprintf(stderr, " num_years=%d  num_months=%d  num_days=%d\n", num_years, num_months, num_days);

  /*
    1. Loop over the Files
    2. Loop over the number of time steps
    3. Set the attributes for the magics data and plot
  */

  // magics_template_parser( magics_node );

  mag_setc("output_name", plotfile);
  mag_setc("subpage_map_projection", "cartesian");
  mag_setr("subpage_y_length", 14.);
  mag_setr("subpage_y_position", 1.5);

  // Horizontal Axis attributes
  mag_setc("axis_orientation", "horizontal");
  mag_setc("axis_grid", "on");
  mag_setc("axis_grid_colour", "grey");
  mag_seti("axis_grid_thickness", 1);
  mag_setc("axis_grid_line_style", "dot");
  mag_setc("axis_type", "date");

  const char *dateType = (num_years > 1) ? "years" : (num_months > 1) ? "months" : (num_months == 1 || num_days) ? "days" : "hours";
  mag_setc("axis_date_type", dateType);

  mag_setc("axis_date_min_value", min_date_time_str);
  mag_setc("axis_date_max_value", max_date_time_str);
  mag_setc("axis_title_text", "Time");
  mag_setc("axis_title_orientation", "horizontal");

  mag_seti("axis_tick_label_frequency", 2);
  mag_setr("axis_years_label_height", 0.4);

  mag_axis();

  // Vertical Axis attributes
  mag_setc("axis_orientation", "vertical");
  mag_setc("axis_grid", "on");
  mag_setc("axis_type", "regular");
  mag_setc("axis_grid_colour", "grey");
  mag_seti("axis_grid_thickness", 1);
  mag_setc("axis_grid_line_style", "dot");

  // To redefine the y- axis scale based on user input in .xml file

  // min & max values from the input data files
  mag_setr("axis_min_value", min_val);
  mag_setr("axis_max_value", max_val);

  // min & max values specified by the user in the command line args
  if (y_min_val < MaxVal) mag_setr("axis_min_value", y_min_val);
  if (y_max_val > MinVal) mag_setr("axis_max_value", y_max_val);

  mag_setc("axis_title_text", varname.c_str());

  mag_setc("axis_title_orientation", "vertical");

  mag_seti("axis_tick_label_frequency", 2);
  mag_setr("axis_tick_label_height", 0.5);

  mag_axis();

  // Legend
  mag_setc("legend", "on");
  mag_setc("legend_text_colour", "black");

  mag_setc("graph_symbol", "off");
  mag_seti("graph_line_thickness", linewidth_val);

  if (DBG) fprintf(stderr, "FILE BEGIN %d\n", file_begin);

  for (i = file_begin; i < nfiles; ++i)
    {
      count = obsv ? i - 1 : i;
      if (DBG) fprintf(stderr, "Current File %ld\n", i);
      // sprintf(legend_text_data, "ens_%d", count + 1);
      sprintf(legend_text_data, "data_%d", count + 1);
      mag_setc("graph_line_colour", line_colours[count % num_colours]);
      mag_setc("legend_user_text", legend_text_data);
      if (stat)
        mag_set1c("graph_curve_date_x_values", (const char **) date_time_str[0], ntime_steps);
      else
        mag_set1c("graph_curve_date_x_values", (const char **) date_time_str[i], nts[i]);

      // TEMPORARY FIX, UNITL NEW MAGICS LIBRARY RELEASE *  begin
      mag_setr("graph_x_suppress_below", (double) LLONG_MIN);
      mag_setr("graph_x_suppress_above", (double) LLONG_MAX);
      // TEMPORARY FIX, UNITL NEW MAGICS LIBRARY RELEASE *  end

      mag_set1r("graph_curve_y_values", datatab[i].data(), nts[i]);
      mag_graph();
    }

  if (obsv)
    {
      mag_setc("graph_line_colour", "black");
      sprintf(legend_text_data, "%s", "Obsv");
      mag_setc("legend_user_text", legend_text_data);
      mag_set1c("graph_curve_date_x_values", (const char **) date_time_str[0], nts[0]);

      // TEMPORARY FIX, UNITL NEW MAGICS LIBRARY RELEASE *  begin
      mag_setr("graph_x_suppress_below", (double) LLONG_MIN);
      mag_setr("graph_x_suppress_above", (double) LLONG_MAX);
      // TEMPORARY FIX, UNITL NEW MAGICS LIBRARY RELEASE *  end

      mag_set1r("graph_curve_y_values", datatab[0].data(), nts[0]);
      mag_setc("graph_line_style", "dot");
      mag_seti("graph_line_thickness", linewidth_val + 2);
      mag_graph();
    }

  if (DBG) fprintf(stderr, "NTIME STEPS %ld\n", ntime_steps);

  if (stat)
    {
      if (DBG) fprintf(stderr, "NTIME STEPS %ld\n", ntime_steps);

      mag_seti("graph_line_thickness", linewidth_val);
      mag_setc("graph_line_colour", "grey");
      mag_setc("graph_line_style", "dash");
      mag_set1c("graph_curve_date_x_values", (const char **) date_time_str[0], ntime_steps);

      // TEMPORARY FIX, UNITL NEW MAGICS LIBRARY RELEASE *  begin
      mag_setr("graph_x_suppress_below", (double) LLONG_MIN);
      mag_setr("graph_x_suppress_above", (double) LLONG_MAX);
      // TEMPORARY FIX, UNITL NEW MAGICS LIBRARY RELEASE *  end

      mag_set1r("graph_curve_y_values", mean_val.data(), ntime_steps);
      sprintf(legend_text_data, "Mean");
      mag_setc("legend_user_text", legend_text_data);
      mag_graph();

      mag_reset("graph_type");
      mag_setc("graph_type", "area");
      mag_seti("graph_line_thickness", 1);
      mag_setc("graph_shade_style", "dot");
      mag_setr("graph_shade_dot_size", 1.);
      mag_set1c("graph_curve2_date_x_values", (const char **) date_time_str[0], ntime_steps);
      mag_set1r("graph_curve2_y_values", spread_max.data(), ntime_steps);
      mag_set1c("graph_curve_date_x_values", (const char **) date_time_str[0], ntime_steps);
      mag_set1r("graph_curve_y_values", spread_min.data(), ntime_steps);
      mag_setc("graph_shade_colour", "grey");
      sprintf(legend_text_data, "%dSigma", num_sigma);
      mag_setc("legend_user_text", legend_text_data);

      // TEMPORARY FIX, UNITL NEW MAGICS LIBRARY RELEASE *  begin
      mag_setr("graph_x_suppress_below", (double) LLONG_MIN);
      mag_setr("graph_x_suppress_above", (double) LLONG_MAX);
      // TEMPORARY FIX, UNITL NEW MAGICS LIBRARY RELEASE *  end

      mag_graph();
    }

  char *lines[1];
  lines[0] = (char *) malloc(1024);
  // To be obtained from Meta Data
  // sprintf( lines[0],"%s","ExpID : " );
  // sprintf( lines[0],"%sxxxx  Variable : %s[%s]",lines[0], varname.c_str(), varunits.c_str() );
  // sprintf( lines[0],"Variable : %s[%s]",varname.c_str(), varunits.c_str() );
  // sprintf( lines[0],"%s  Date : %s --%s",lines[0], min_date_time_str, max_date_time_str );
  sprintf(lines[0], "Variable : %s[%s]  Date : %s --%s", varname.c_str(), varunits.c_str(), min_date_time_str, max_date_time_str);
  mag_set1c("text_lines", (const char **) lines, 1);

  mag_setc("text_html", "true");
  mag_setc("text_colour", "black");
  mag_setr("text_font_size", 0.6);
  mag_setc("text_mode", "positional");
  mag_setr("text_box_x_position", 1.5);
  mag_setr("text_box_y_position", 16.5);
  mag_setr("text_box_x_length", 20.);
  mag_setr("text_box_y_length", 2.5);
  mag_setc("text_border", "off");
  mag_setc("text_justification", "left");
  mag_text();

  free(date_time_str);

  if (DBG) fprintf(stderr, "lines=%s\n", lines[0]);

  free(lines[0]);
}

static void
init_MAGICS()
{
  setenv("MAGPLUS_QUIET", "1", 1); /* To suppress magics messages */
  mag_open();

  // Some standard parameters affectng the magics environment, moved from the xml file  ** begin **
  mag_setc("page_id_line", "off");
  // Some standard parameters affectng the magics environment, moved from the xml file  ** end **
}

static void
quit_MAGICS()
{
  mag_close();
  if (DBG) fprintf(stdout, "Exiting From MAGICS\n");
}

static void
VerifyGraphParameters(int num_param, std::vector<std::string> &param_names)
{
  int i, j;
  bool found = false, syntax = true, halt_flag = false;
  char *temp_str;

  for (i = 0; i < num_param; ++i)
    {
      found = false;
      syntax = true;
      auto splitStrings = split_with_seperator(param_names[i], '=');

      if (splitStrings.size() > 1)
        {
          const auto &key = splitStrings[0];
          const auto &value = splitStrings[1];

          for (j = 0; j < graph_param_count; ++j)
            {
              if (key == graph_params[j])
                {
                  found = true;
                  if (key == "obsv" || key == "stat")
                    {
                      if (string_is_float(value))
                        syntax = false;
                      else
                        {
                          temp_str = strdup(value.c_str());
                          cstr_to_lower(temp_str);
                          if (strcmp(temp_str, "true") && strcmp(temp_str, "false")) syntax = false;
                        }
                    }

                  if (key == "ymin" || key == "ymax" || key == "sigma" || key == "linewidth")
                    {
                      if (!string_is_float(value)) syntax = false;
                    }

                  if (key == "device")
                    {
                      if (string_is_float(value))
                        syntax = false;
                      else
                        {
                          if (DBG) fprintf(stderr, "Parameter value '%s'\n", value.c_str());
                          char *deviceCstr = strdup(value.c_str());
                          if (checkdevice(deviceCstr)) syntax = false;

                          // Graph not supported in google earth format
                          if (value == "GIF_ANIMATION" || value == "gif_animation")
                            {
                              syntax = false;
                              fprintf(stderr, "Animation not supported for Graph!\n");
                              if (DBG) fprintf(stderr, "Parameter value '%s'\n", value.c_str());
                            }
                          if (value == "KML" || value == "kml")
                            {
                              syntax = false;
                              fprintf(stderr, " 'kml' format not supported for  Graph!\n");
                              if (DBG) fprintf(stderr, "Parameter value '%s'\n", value.c_str());
                            }
                        }
                    }

                  /*
                    if(key == "xml")
                      {
                        if( ( fp = std::fopen(value.c_str(),"r") ) == nullptr )
                          {
                            fprintf( stderr,"Input XML File not found in specified path '%s'\n", value.c_str() );
                            halt_flag = true;
                          }
                        else
                          {
                            // HARDCODED THE FILE NAME .. TO BE SENT AS COMMAND LINE ARGUMENT FOR THE MAGICS OPERATOR
                            std::fclose(fp);
                            init_XML_template_parser( value.c_str() ); updatemagics_and_results_nodes( );
                          }
                      }
                  */
                }
            }
        }
      else { syntax = false; }

      if (!found)
        {
          halt_flag = true;
          fprintf(stderr, "Unknown parameter  '%s'!\n", param_names[i].c_str());
        }
      if (found && !syntax)
        {
          halt_flag = true;
          fprintf(stderr, "Invalid parameter specification  '%s'!\n", param_names[i].c_str());
        }
    }

  if (halt_flag) exit(0);
}
#endif

void *
Maggraph(void *process)
{
  cdo_initialize(process);

#ifdef HAVE_LIBMAGICS
  std::string varname, units;
  int gridID;
  int vlistID0 = -1;

  int nparam = cdo_operator_argc();
  auto pnames = cdo_get_oper_argv();

  if (nparam) VerifyGraphParameters(nparam, pnames);

  int nfiles = cdo_stream_cnt() - 1;
  const char *ofilename = cdo_get_stream_name(nfiles);

  if (DBG)
    {
      fprintf(stderr, " Num of files %d\n", nfiles);
      fprintf(stderr, " files %s\n", ofilename);
    }

  std::vector<std::vector<double>> datatab(nfiles);
  std::vector<std::vector<CdiDateTime>> vDateTimes(nfiles);
  std::vector<long> nts(nfiles, 0);

  for (int fileID = 0; fileID < nfiles; ++fileID)
    {
      if (DBG) fprintf(stderr, " file %d is %s\n", fileID, cdo_get_stream_name(fileID));

      auto streamID = cdo_open_read(fileID);
      auto vlistID = cdo_stream_inq_vlist(streamID);
      auto taxisID = vlistInqTaxis(vlistID);

      units = cdo::inq_var_units(vlistID, 0);
      if (DBG) fprintf(stderr, " units=%s\n", units.c_str());
      if (fileID == 0)
        {
          varname = cdo::inq_var_name(vlistID, 0);

          gridID = vlistInqVarGrid(vlistID, 0);
          auto zaxisID = vlistInqVarZaxis(vlistID, 0);
          auto nvars = vlistNvars(vlistID);

          if (nvars > 1) cdo_abort("Input stream has more than on variable!");
          if (gridInqSize(gridID) != 1) cdo_abort("Variable has more than one grid point!");
          if (zaxisInqSize(zaxisID) != 1) cdo_abort("Variable has more than one level!");

          vlistID0 = vlistDuplicate(vlistID);
        }
      else { vlist_compare(vlistID0, vlistID, CMP_ALL); }

      int tsID = 0;
      size_t numTsAlloc = 0;
      while (true)
        {
          auto nrecs = cdo_stream_inq_timestep(streamID, tsID);
          if (nrecs == 0) break;

          if (nrecs != 1) cdo_abort("Input stream has more than one point in time!");

          if ((size_t) tsID >= numTsAlloc)
            {
              constexpr size_t NALLOC_INC = 1024;
              numTsAlloc += NALLOC_INC;
              datatab[fileID].resize(numTsAlloc);
              vDateTimes[fileID].resize(numTsAlloc);
            }

          nts[fileID]++;

          int varID, levelID;
          cdo_inq_record(streamID, &varID, &levelID);
          size_t nmiss;
          double val;
          cdo_read_record(streamID, &val, &nmiss);
          datatab[fileID][tsID] = val;
          vDateTimes[fileID][tsID] = taxisInqVdatetime(taxisID);

          tsID++;
        }

      cdo_stream_close(streamID);
    }

  // HARDCODED THE FILE NAME .. TO BE SENT AS COMMAND LINE ARGUMENT FOR THE MAGICS OPERATOR
  // init_XML_template_parser( Filename );
  // updatemagics_and_results_nodes( );

  init_MAGICS();

  cdo_print("Creating PLOT for %s", varname);
  if (DBG)
    {
      fprintf(stderr, "Num params %d\n", nparam);
      for (int i = 0; i < nparam; ++i) fprintf(stderr, "Param %s\n", pnames[i].c_str());
    }

  maggraph(ofilename, varname, units, nfiles, nts, vDateTimes, datatab, nparam, pnames);

  // quit_XML_template_parser( );

  quit_MAGICS();

  if (vlistID0 != -1) vlistDestroy(vlistID0);

#else

  cdo_abort("MAGICS support not compiled in!");

#endif

  cdo_finish();

  return nullptr;
}
