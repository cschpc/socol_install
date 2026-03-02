/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

*/

#include <cdi.h>

#include "cdo_options.h"
#include "cdo_vlist.h"
#include "process_int.h"
#include "param_conversion.h"

class HistoGram
{
private:
  CdoStreamID streamID1;
  CdoStreamID streamID2;

  int taxisID1;
  int taxisID2;

  Varray2D<double> vardata;
  Varray2D<double> varcount;
  Varray2D<double> vartcount;

  Varray<double> array;
  VarList varList1, varList2;

  int nbins;

  std::vector<double> fltarr;
  int nvars;
  int operatorID;

  int HISTCOUNT;
  int HISTSUM;
  int HISTMEAN;
  int HISTFREQ;

public:
  void
  init(void *process)
  {
    cdo_initialize(process);

    // clang-format off
    HISTCOUNT = cdo_operator_add("histcount", 0, 0, nullptr);
    HISTSUM   = cdo_operator_add("histsum",   0, 0, nullptr);
    HISTMEAN  = cdo_operator_add("histmean",  0, 0, nullptr);
    HISTFREQ  = cdo_operator_add("histfreq",  0, 0, nullptr);
    // clang-format on

    (void) (HISTSUM);  // unused

    operatorID = cdo_operator_id();

    operator_input_arg("bins");

    fltarr = cdo_argv_to_flt(cdo_get_oper_argv());
    nbins = fltarr.size() - 1;
    if (nbins < 1) cdo_abort("Too few arguments!");

    if (Options::cdoVerbose)
      {
        printf("nbins = %d\n", nbins);
        for (int i = 0; i < nbins; ++i) printf("flt %d = %g\n", i + 1, fltarr[i]);
      }

    streamID1 = cdo_open_read(0);
    const auto vlistID1 = cdo_stream_inq_vlist(streamID1);

    taxisID1 = vlistInqTaxis(vlistID1);

    const auto vlistID2 = vlistDuplicate(vlistID1);

    /* create zaxis for output bins */
    const auto zaxisID2 = zaxisCreate(ZAXIS_GENERIC, nbins);

    zaxisDefLevels(zaxisID2, &fltarr[0]);
    zaxisDefLbounds(zaxisID2, &fltarr[0]);
    zaxisDefUbounds(zaxisID2, &fltarr[1]);

    cdiDefKeyString(zaxisID2, CDI_GLOBAL, CDI_KEY_NAME, "bin");
    cdiDefKeyString(zaxisID2, CDI_GLOBAL, CDI_KEY_LONGNAME, "histogram bins");
    cdiDefKeyString(zaxisID2, CDI_GLOBAL, CDI_KEY_UNITS, "level");

    /* check zaxis: only 2D fields allowed */
    const auto nzaxis = vlistNzaxis(vlistID1);
    for (int index = 0; index < nzaxis; ++index)
      {
        const auto zaxisID = vlistZaxis(vlistID1, index);
        const auto nlevel = zaxisInqSize(zaxisID);
        if (nlevel > 1) cdo_abort("Found 3D field with %d levels. Only 2D fields allowed!", nlevel);
        vlistChangeZaxisIndex(vlistID2, index, zaxisID2);
      }

    streamID2 = cdo_open_write(1);

    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    cdo_def_vlist(streamID2, vlistID2);

    varListInit(varList1, vlistID1);
    varListInit(varList2, vlistID2);

    nvars = vlistNvars(vlistID2);

    vardata = Varray2D<double>(nvars);
    varcount = Varray2D<double>(nvars);
    vartcount = Varray2D<double>(nvars);

    for (int varID = 0; varID < nvars; ++varID)
      {
        const auto gridsize = varList1[varID].gridsize;
        vardata[varID].resize(nbins * gridsize, 0);
        varcount[varID].resize(nbins * gridsize, 0);
        vartcount[varID].resize(gridsize, 0);
      }

    const auto gridsizemax = vlistGridsizeMax(vlistID1);
    array = Varray<double>(gridsizemax);
  }
  void
  run()
  {
    int tsID = 0;
    while (true)
      {
        const auto nrecs = cdo_stream_inq_timestep(streamID1, tsID++);
        if (nrecs == 0) break;

        cdo_taxis_copy_timestep(taxisID2, taxisID1);

        for (int recID = 0; recID < nrecs; ++recID)
          {
            int varID, levelID;
            cdo_inq_record(streamID1, &varID, &levelID);
            size_t nmiss;
            cdo_read_record(streamID1, array.data(), &nmiss);

            const auto gridsize = varList1[varID].gridsize;
            const auto missval = varList1[varID].missval;

            nmiss = 0;
            for (size_t i = 0; i < gridsize; ++i)
              {
                if (!DBL_IS_EQUAL(array[i], missval))
                  {
                    vartcount[varID][i] += 1;
                    int index = 0;
                    while (index < nbins)
                      {
                        const auto offset = gridsize * index;
                        if (!DBL_IS_EQUAL(vardata[varID][offset + i], missval) && array[i] >= fltarr[index]
                            && array[i] < fltarr[index + 1])
                          {
                            vardata[varID][offset + i] += array[i];
                            varcount[varID][offset + i] += 1;
                            break;
                          }
                        index++;
                      }
                  }
                else
                  {
                    nmiss++;  // missing value
                  }
              }
          }
      }

    cdo_def_timestep(streamID2, 0);

    for (int varID = 0; varID < nvars; ++varID)
      {
        const auto gridsize = varList2[varID].gridsize;
        const auto missval = varList2[varID].missval;

        // fix missing values

        for (int index = 0; index < nbins; ++index)
          {
            size_t nmiss = 0;
            const auto offset = gridsize * index;

            for (size_t i = 0; i < gridsize; ++i)
              {
                if (vartcount[varID][i] > 0)
                  {
                    if (operatorID == HISTMEAN || operatorID == HISTFREQ)
                      {
                        if (varcount[varID][offset + i] > 0)
                          {
                            if (operatorID == HISTMEAN)
                              vardata[varID][offset + i] /= varcount[varID][offset + i];
                            else
                              vardata[varID][offset + i] = varcount[varID][offset + i] / vartcount[varID][i];
                          }
                      }
                  }
                else
                  {
                    nmiss++;
                    varcount[varID][offset + i] = missval;
                    vardata[varID][offset + i] = missval;
                  }
              }

            cdo_def_record(streamID2, varID, index);

            if (operatorID == HISTCOUNT)
              cdo_write_record(streamID2, &varcount[varID][offset], nmiss);
            else
              cdo_write_record(streamID2, &vardata[varID][offset], nmiss);
          }
      }
  }

  void
  close()
  {
    cdo_stream_close(streamID1);
    cdo_stream_close(streamID2);

    cdo_finish();
  }
};

void *
Histogram(void *process)
{
  HistoGram histogram;
  histogram.init(process);
  histogram.run();
  histogram.close();
  return nullptr;
}
