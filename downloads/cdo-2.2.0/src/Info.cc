/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Info       info            Dataset information
      Info       map             Dataset information and simple map
*/

#include <cfloat>

#include <cdi.h>

#include "cdo_options.h"
#include "cdo_vlist.h"
#include "process_int.h"
#include "mpmo_color.h"
#include "varray.h"
#include "datetime.h"
#include "printinfo.h"
#include "cdo_zaxis.h"

struct Infostat
{
  double min = DBL_MAX;
  double max = -DBL_MAX;
  double sum = 0.0;
  double sumi = 0.0;
  size_t nvals = 0;
  size_t nmiss = 0;
  int nlevels = 0;
};

static void
field_min_max_sum(const Field &field, double &min, double &max, double &sum)
{
  auto mms = MinMaxSum(min, max, sum);
  if (field.memType == MemType::Float)
    mms = varray_min_max_sum(field.size, field.vec_f, mms);
  else
    mms = varray_min_max_sum(field.size, field.vec_d, mms);

  min = mms.min;
  max = mms.max;
  sum = mms.sum;
}

static size_t
field_min_max_sum_mv(const Field &field, double &min, double &max, double &sum)
{
  auto mms = MinMaxSum(min, max, sum);
  if (field.memType == MemType::Float)
    mms = varray_min_max_sum_mv(field.size, field.vec_f, static_cast<float>(field.missval), mms);
  else
    mms = varray_min_max_sum_mv(field.size, field.vec_d, field.missval, mms);

  min = mms.min;
  max = mms.max;
  sum = mms.sum;
  return mms.n;
}

static void
print_grid_index(int nlon, int nlat, int i)
{
  int index = (nlat < 10) ? 2 : (nlat < 100) ? 3 : (nlat < i) ? 4 : 5;

  std::stringstream s;
  s << std::string(index, ' ');
  for (int ilon = 0; ilon < nlon; ilon++) s << ((ilon + 1) / i) % 10;

  printf("%s\n", s.str().c_str());
}

template <typename T>
static void
print_map(int nlon, int nlat, const Varray<T> &varray, double missval, double min, double max)
{
  // source code from PINGO
  double level[10] = {};
  int bmin = 1, bmax = 1;
  unsigned char c;

  auto step = (max - min) / 10;

  if (is_not_equal(step, 0))
    {
      auto a = std::pow(10, std::floor(std::log(step) / M_LN10));
      auto b = step / a;

      // clang-format off
      if      (b > 5) b = 0.5 * std::ceil(b / 0.5);
      else if (b > 2) b = 0.2 * std::ceil(b / 0.2);
      else if (b > 1) b = 0.1 * std::ceil(b / 0.1);
      else            b = 1;
      // clang-format on

      step = b * a;

      if (min < 0 && max > 0)
        {
          int min_n = (int) std::floor(10 * (-min) / (max - min) - 0.5);
          int max_n = (int) std::ceil(10 * (-min) / (max - min) - 0.5);
          level[min_n] = 0;
          for (int i = min_n - 1; i >= 0; i--) level[i] = level[i + 1] - step;
          for (int i = max_n; i < 9; ++i) level[i] = level[i - 1] + step;
        }
      else
        {
          level[0] = step * std::ceil(min / step + 0.5);
          for (int i = 1; i < 9; ++i) level[i] = level[i - 1] + step;
        }
    }
  else
    for (int i = 0; i < 9; ++i) level[i] = min;

  printf("\n");

  for (int i = 1; i <= 4; ++i)
    {
      int current = 10000 / std::pow(10, i);
      if (nlon >= current) print_grid_index(nlon, nlat, current);
    }
  printf("\n");

  for (int ilat = 0; ilat < nlat; ilat++)
    {
      printf("%0*d ", (nlat < 10) ? 1 : (nlat < 100) ? 2 : (nlat < 1000) ? 3 : 4, ilat + 1);
      for (int ilon = 0; ilon < nlon; ilon++)
        {
          auto x = varray[ilat * nlon + ilon];
          if (dbl_is_equal(x, missval))
            c = '.';
          else if (dbl_is_equal(x, min) && !dbl_is_equal(min, max))
            c = 'm';
          else if (dbl_is_equal(x, max) && !dbl_is_equal(min, max))
            c = 'M';
          else if (dbl_is_equal(x, 0.0))
            c = '*';
          else if (x < 0)
            {
              c = '9';
              for (int i = 0; i < 9; ++i)
                if (level[i] > x)
                  {
                    c = i + '0';
                    break;
                  }
            }
          else
            {
              c = '0';
              for (int i = 8; i >= 0; i--)
                if (level[i] < x)
                  {
                    c = i + 1 + '0';
                    break;
                  }
            }

          TextMode mode(MODELESS);
          TextColor color(BLACK);
          switch (c)
            {
            // clang-format off
            case '0': mode = BRIGHT   ; color = BLUE    ; break ;
            case '1': mode = MODELESS ; color = BLUE    ; break ;
            case '2': mode = BRIGHT   ; color = CYAN    ; break ;
            case '3': mode = MODELESS ; color = CYAN    ; break ;
            case '4': mode = MODELESS ; color = GREEN   ; break ;
            case '5': mode = MODELESS ; color = YELLOW  ; break ;
            case '6': mode = MODELESS ; color = RED     ; break ;
            case '7': mode = BRIGHT   ; color = RED     ; break ;
            case '8': mode = MODELESS ; color = MAGENTA ; break ;
            case '9': mode = BRIGHT   ; color = MAGENTA ; break ;
            // clang-format on
            case 'm':
              (bmax == 1) ? mode = BLINK : mode = MODELESS, color = BLACK;
              if (bmax) bmax = 0;
              break;
            case 'M':
              (bmin == 1) ? mode = BLINK : mode = MODELESS, color = BLACK;
              if (bmin) bmin = 0;
              break;
            }

          set_text_color(stdout, mode, color);
          putchar(c);
          reset_text_color(stdout);
        }
      printf(" %0*d\n", (nlat < 10) ? 1 : (nlat < 100) ? 2 : (nlat < 1000) ? 3 : 4, ilat + 1);
    }
  printf("\n");

  for (int i = 1; i <= 4; ++i)
    {
      int current = 10000 / std::pow(10, i);
      if (nlon >= current) print_grid_index(nlon, nlat, current);
    }
  printf("\n");

  for (int i = 0; i < 10; ++i)
    {
      printf("%d=%c%+9.3e,%+9.3e%c%s", i, '[', (i == 0) ? min : level[i - 1], (i == 9) ? max : level[i],
             ']', (i != 2 && i != 5 && i != 8) ? "  " : "");

      if (i == 2 || i == 5 || i == 8) printf("\n");
    }

  printf("*=0  .=miss  m=min=%+9.3e  M=max=%+9.3e\n", min, max);
  printf("\n");
}

static void
print_map(int nlon, int nlat, const Field &field, const Infostat &infostat)
{
  if (field.memType == MemType::Float)
    print_map(nlon, nlat, field.vec_f, static_cast<float>(field.missval), infostat.min, infostat.max);
  else
    print_map(nlon, nlat, field.vec_d, field.missval, infostat.min, infostat.max);
}

template <typename T>
static size_t
complex_sum(size_t gridsize, const Varray<T> &varray, double missval, double &sum, double &sumi)
{
  size_t n = 0;
  for (size_t i = 0; i < gridsize; ++i)
    {
      if (!dbl_is_equal(varray[i * 2], missval) && !dbl_is_equal(varray[i * 2 + 1], missval))
        {
          sum += varray[i * 2];
          sumi += varray[i * 2 + 1];
          n++;
        }
    }

  return n;
}

static size_t
field_complex_sum(const Field &field, double &sum, double &sumi)
{
  if (field.memType == MemType::Float)
    return complex_sum(field.gridsize, field.vec_f, static_cast<float>(field.missval), sum, sumi);
  else
    return complex_sum(field.gridsize, field.vec_d, field.missval, sum, sumi);
}

static void
infostat_init(Infostat &infostat)
{
  infostat.nvals = 0;
  infostat.nmiss = 0;
  infostat.nlevels = 0;
  infostat.min = DBL_MAX;
  infostat.max = -DBL_MAX;
  infostat.sum = 0.0;
  infostat.sumi = 0.0;
}

static void
print_header(int fileIndex, bool lvinfo, const std::string &e, const std::string &v)
{
  set_text_color(stdout, BRIGHT);
  if (fileIndex)
    fprintf(stdout, "%6d :       Date     Time   %s Gridsize    Miss :     Minimum        Mean     Maximum : %s%s\n", fileIndex,
            lvinfo ? "Nlevs" : "Level", e.c_str(), v.c_str());
  else
    fprintf(stdout, "       :       Date     Time   %s Gridsize    Miss :     Minimum        Mean     Maximum : %s%s\n",
            lvinfo ? "Nlevs" : "Level", e.c_str(), v.c_str());
  reset_text_color(stdout);
}

static void
compute_stat_real(const Field &field, Infostat &infostat, size_t &imiss, size_t gridsize)
{
  if (infostat.nmiss)
    {
      auto nvals = field_min_max_sum_mv(field, infostat.min, infostat.max, infostat.sum);
      imiss = gridsize - nvals;
      infostat.nvals += nvals;
    }
  else if (gridsize == 1)
    {
      double val = (field.memType == MemType::Float) ? field.vec_f[0] : field.vec_d[0];
      infostat.sum = (infostat.nvals == 0) ? val : infostat.sum + val;
      infostat.nvals += 1;
    }
  else
    {
      field_min_max_sum(field, infostat.min, infostat.max, infostat.sum);
      infostat.nvals += gridsize;
    }
}

static void
compute_stat_comp(const Field &field, Infostat &infostat, size_t &imiss, size_t gridsize)
{
  auto nvals = field_complex_sum(field, infostat.sum, infostat.sumi);

  imiss = gridsize - nvals;
  infostat.nvals += nvals;
}

static void
print_stat_real(const Infostat &infostat)
{
  if (infostat.nvals == 0)
    fprintf(stdout, "                     nan            ");
  else if (infostat.nvals == 1)
    fprintf(stdout, "            %#12.5g            ", infostat.sum);
  else
    fprintf(stdout, "%#12.5g%#12.5g%#12.5g", infostat.min, infostat.sum / infostat.nvals, infostat.max);
}

static void
print_stat_comp(const Infostat &infostat)
{
  auto arrmean_r = (infostat.nvals > 0) ? infostat.sum / infostat.nvals : 0.0;
  auto arrmean_i = (infostat.nvals > 0) ? infostat.sumi / infostat.nvals : 0.0;
  fprintf(stdout, "   -  (%#12.5g,%#12.5g)  -", arrmean_r, arrmean_i);
}

void *
Info(void *process)
{
  enum
  {
    E_NAME,
    E_CODE,
    E_PARAM
  };
  size_t imiss = 0;
  char paramstr[32];

  cdo_initialize(process);

  // clang-format off
                cdo_operator_add("info",   E_PARAM,  0, nullptr);
                cdo_operator_add("infop",  E_PARAM,  0, nullptr);
                cdo_operator_add("infon",  E_NAME,   0, nullptr);
                cdo_operator_add("infoc",  E_CODE,   0, nullptr);
  auto VINFON = cdo_operator_add("vinfon", E_NAME,   0, nullptr);
  auto XINFON = cdo_operator_add("xinfon", E_NAME,   0, nullptr);
  auto MAP    = cdo_operator_add("map",    E_PARAM,  0, nullptr);
  // clang-format on

  auto operatorID = cdo_operator_id();
  auto operfunc = cdo_operator_f1(operatorID);

  operator_check_argc(0);

  auto lvinfo = (operatorID == VINFON || operatorID == XINFON);

  DateTimeList dtlist;

  std::string e = (operfunc == E_NAME) ? "Parameter name" : ((operfunc == E_CODE) ? "Code number" : "Parameter ID");

  std::string v = (Options::cdoVerbose) ? " : Extra" : "";
  int indg = 0;

  for (int indf = 0; indf < cdo_stream_cnt(); indf++)
    {
      auto streamID = cdo_open_read(indf);
      auto vlistID = cdo_stream_inq_vlist(streamID);
      auto taxisID = vlistInqTaxis(vlistID);

      auto nvars = vlistNvars(vlistID);
      if (nvars == 0) continue;

      std::vector<Infostat> infostat(nvars);

      VarList varList;
      varListInit(varList, vlistID);

      Field field;

      indg = 0;
      int tsID = 0;
      while (true)
        {
          auto nrecs = cdo_stream_inq_timestep(streamID, tsID);
          if (nrecs == 0) break;

          dtlist.taxis_inq_timestep(taxisID, 0);
          auto vdateString = date_to_string(dtlist.get_vDateTime(0).date);
          auto vtimeString = time_to_string(dtlist.get_vDateTime(0).time);

          for (int varID = 0; varID < nvars; ++varID) infostat_init(infostat[varID]);

          for (int recID = 0; recID < nrecs; ++recID)
            {
              if ((tsID == 0 && recID == 0) || operatorID == MAP) print_header(-(indf + 1), lvinfo, e, v);

              int varID, levelID;
              cdo_inq_record(streamID, &varID, &levelID);
              const auto &var = varList[varID];
              field.init(var);
              cdo_read_record(streamID, field);
              auto nmiss = field.nmiss;

              indg = lvinfo ? varID + 1 : indg + 1;

              auto gridsize = var.gridsize;

              auto loutput = !lvinfo;

              if (loutput) infostat_init(infostat[varID]);

              auto &infostatr = infostat[varID];
              infostatr.nlevels += 1;
              infostatr.nmiss += nmiss;

              if (var.nlevels == infostatr.nlevels) loutput = true;

              if (loutput)
                {
                  cdiParamToString(var.param, paramstr, sizeof(paramstr));

                  fprintf(stdout, "%6d ", indg);
                  fprintf(stdout, ":");

                  set_text_color(stdout, MAGENTA);
                  fprintf(stdout, "%s %s ", vdateString.c_str(), vtimeString.c_str());
                  reset_text_color(stdout);

                  set_text_color(stdout, GREEN);
                  if (lvinfo)
                    fprintf(stdout, "%7d ", var.nlevels);
                  else
                    fprintf(stdout, "%7g ", cdo_zaxis_inq_level(var.zaxisID, levelID));

                  fprintf(stdout, "%8zu %7zu ", gridsize, infostatr.nmiss);
                  reset_text_color(stdout);

                  fprintf(stdout, ":");

                  set_text_color(stdout, BLUE);
                }

              // clang-format off
              if (var.nwpv == CDI_REAL) compute_stat_real(field, infostatr, imiss, gridsize);
              else                      compute_stat_comp(field, infostatr, imiss, gridsize);
              // clang-format on

              if (loutput)
                {
                  // clang-format off
                  if (var.nwpv == CDI_REAL) print_stat_real(infostatr);
                  else                      print_stat_comp(infostatr);
                  // clang-format on

                  reset_text_color(stdout);

                  fprintf(stdout, " : ");

                  // set_text_color(stdout, GREEN);
                  // clang-format off
                  if      (operfunc == E_NAME) fprintf(stdout, "%-14s", var.name.c_str());
                  else if (operfunc == E_CODE) fprintf(stdout, "%4d   ", var.code);
                  else                         fprintf(stdout, "%-14s", paramstr);
                  // clang-format on
                  // reset_text_color(stdout);

                  fprintf(stdout, "\n");
                }

              if (imiss != nmiss && nmiss) cdo_warning("Found %zu of %zu missing values!", imiss, nmiss);

              if (operatorID == MAP)
                {
                  auto gridID = var.gridID;
                  auto gridtype = gridInqType(gridID);
                  auto nlon = gridInqXsize(gridID);
                  auto nlat = gridInqYsize(gridID);

                  if (gridtype == GRID_GAUSSIAN || gridtype == GRID_LONLAT || gridtype == GRID_CURVILINEAR
                      || (gridtype == GRID_GENERIC && nlon * nlat == gridInqSize(gridID) && nlon < 1024))
                    {
                      print_map(nlon, nlat, field, infostatr);
                    }
                }
            }

          tsID++;
        }

      cdo_stream_close(streamID);
    }

  if (indg > 36 && operatorID != MAP) print_header(0, lvinfo, e, v);

  cdo_finish();

  return nullptr;
}
