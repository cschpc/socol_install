/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include "cdo_options.h"
#include "cdo_vlist.h"
#include "param_conversion.h"
#include "process_int.h"
#include "timer.h"
#include "util_files.h"

static void
print_stat(const char *sinfo, MemType memtype, int datatype, int filetype, off_t nvalues, double data_size, double file_size,
           double tw)
{
  nvalues /= 1000000;
  data_size /= 1024. * 1024. * 1024.;

  cdo_print("%s Read %.1f GB of %d bit floats from %s %s, %.1f MVal/s", sinfo, data_size, (memtype == MemType::Float) ? 32 : 64,
            cdo::datatype_to_cstr(datatype), cdo::filetype_to_cstr(filetype), (tw > 0) ? nvalues / tw : -1);

  file_size /= 1024. * 1024. * 1024.;
  cdo_print("%s Read %.1f GB in %.1f seconds, total %.1f MB/s", sinfo, file_size, tw, (tw > 0) ? 1024 * file_size / tw : -1);
}

void *
CDIread(void *process)
{
  const auto memtype = Options::CDO_Memtype;
  int filetype = -1, datatype = -1;
  char sinfo[64];
  off_t nvalues = 0;
  double file_size = 0, data_size = 0;
  double twsum = 0;

  sinfo[0] = 0;

  cdo_initialize(process);

  if (Options::cdoVerbose) cdo_print("parameter: <nruns>");

  if (cdo_operator_argc() > 1) cdo_abort("Too many arguments!");

  auto nruns = (cdo_operator_argc() == 1) ? parameter_to_int(cdo_operator_argv(0)) : 1;
  if (nruns < 0) nruns = 0;
  if (nruns > 99) nruns = 99;

  if (Options::cdoVerbose) cdo_print("nruns      : %d", nruns);

  // vlistDefNtsteps(vlistID, 1);

  for (int irun = 0; irun < nruns; ++irun)
    {
      const auto tw0 = timer_val(timer_read);
      data_size = 0;
      nvalues = 0;

      const auto streamID = cdo_open_read(0);
      const auto vlistID = cdo_stream_inq_vlist(streamID);

      VarList varList;
      varListInit(varList, vlistID);

      filetype = cdo_inq_filetype(streamID);
      datatype = vlistInqVarDatatype(vlistID, 0);

      const auto gridsizemax = vlistGridsizeMax(vlistID);

      Varray<float> farray;
      Varray<double> darray;
      if (memtype == MemType::Float)
        farray.resize(gridsizemax);
      else
        darray.resize(gridsizemax);

      auto t0 = timer_val(timer_read);

      int tsID = 0;
      while (true)
        {
          const auto nrecs = cdo_stream_inq_timestep(streamID, tsID);
          if (nrecs == 0) break;

          for (int recID = 0; recID < nrecs; ++recID)
            {
              int varID, levelID;
              cdo_inq_record(streamID, &varID, &levelID);

              const auto gridsize = varList[varID].gridsize;
              nvalues += gridsize;

              size_t nmiss;
              if (memtype == MemType::Float)
                {
                  cdo_read_record_f(streamID, farray.data(), &nmiss);
                  data_size += gridsize * 4;
                }
              else
                {
                  cdo_read_record(streamID, darray.data(), &nmiss);
                  data_size += gridsize * 8;
                }
            }

          if (Options::cdoVerbose)
            {
              const auto tw = timer_val(timer_read) - t0;
              t0 = timer_val(timer_read);
              cdo_print("Timestep %d: %.2f seconds", tsID + 1, tw);
            }

          tsID++;
        }

      cdo_stream_close(streamID);

      const auto tw = timer_val(timer_read) - tw0;
      twsum += tw;

      file_size = (double) FileUtils::size(cdo_get_stream_name(0));

      if (nruns > 1) std::snprintf(sinfo, sizeof(sinfo), "(run %d)", irun + 1);

      print_stat(sinfo, memtype, datatype, filetype, nvalues, data_size, file_size, tw);
    }

  if (nruns > 1) print_stat("(mean)", memtype, datatype, filetype, nvalues, data_size, file_size, twsum / nruns);

  cdo_finish();

  return nullptr;
}
