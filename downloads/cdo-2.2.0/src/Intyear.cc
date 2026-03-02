/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Intyear    intyear         Year interpolation
*/

#include <cdi.h>

#include "cdo_rlimit.h"
#include "process_int.h"
#include "cdo_vlist.h"
#include "util_files.h"
#include "param_conversion.h"

static size_t
intlin_year(double fac1, double fac2, size_t gridsize, const Varray<double> &array1, const Varray<double> &array2,
            Varray<double> &array3, bool withMissval, double missval1, double missval2)
{
  size_t nmiss3 = 0;

  if (withMissval)
    {
      for (size_t i = 0; i < gridsize; ++i)
        {
          if (!dbl_is_equal(array1[i], missval1) && !dbl_is_equal(array2[i], missval2))
            {
              array3[i] = array1[i] * fac1 + array2[i] * fac2;
            }
          else
            {
              array3[i] = missval1;
              nmiss3++;
            }
        }
    }
  else
    {
      for (size_t i = 0; i < gridsize; ++i) array3[i] = array1[i] * fac1 + array2[i] * fac2;
    }

  return nmiss3;
}

class ModuleIntYear
{
  CdoStreamID streamID1;
  CdoStreamID streamID2;

  int taxisID1;
  int taxisID2;
  int taxisID3;

  std::vector<int> iyears;
  int nyears;

  VarList varList1, varList2;

  std::vector<CdoStreamID> streamIDs;

  Varray<double> array1;
  Varray<double> array2;
  Varray<double> array3;

public:
  void
  init(void *process)
  {
    cdo_initialize(process);

    operator_input_arg("years");

    iyears = cdo_argv_to_int(cdo_get_oper_argv());
    nyears = iyears.size();

    cdo::set_numfiles(nyears + 8);

    streamIDs = std::vector<CdoStreamID>(nyears);

    streamID1 = cdo_open_read(0);
    streamID2 = cdo_open_read(1);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = cdo_stream_inq_vlist(streamID2);
    auto vlistID3 = vlistDuplicate(vlistID1);

    vlist_compare(vlistID1, vlistID2, CMP_ALL);

    auto gridsizemax = vlistGridsizeMax(vlistID1);

    array1 = Varray<double>(gridsizemax);
    array2 = Varray<double>(gridsizemax);
    array3 = Varray<double>(gridsizemax);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = vlistInqTaxis(vlistID2);
    taxisID3 = taxisDuplicate(taxisID1);
    if (taxisHasBounds(taxisID3)) taxisDeleteBounds(taxisID3);
    vlistDefTaxis(vlistID3, taxisID3);

    char filename[8192];
    strcpy(filename, cdo_get_obase().c_str());
    const int nchars = strlen(filename);

    auto refname = cdo_get_stream_name(0);
    char filesuffix[32] = { 0 };
    FileUtils::gen_suffix(filesuffix, sizeof(filesuffix), cdo_inq_filetype(streamID1), vlistID1, refname);

    for (int iy = 0; iy < nyears; iy++)
      {
        sprintf(filename + nchars, "%04d", iyears[iy]);
        if (filesuffix[0]) sprintf(filename + nchars + 4, "%s", filesuffix);

        streamIDs[iy] = cdo_open_write(filename);
        cdo_def_vlist(streamIDs[iy], vlistID3);
      }

    varListInit(varList1, vlistID1);
    varListInit(varList2, vlistID2);
  }

  void
  run()
  {
    int tsID = 0;
    while (true)
      {
        auto nrecs = cdo_stream_inq_timestep(streamID1, tsID);
        if (nrecs == 0) break;
        nrecs = cdo_stream_inq_timestep(streamID2, tsID);
        if (nrecs == 0) cdo_abort("Too few timesteps in second inputfile!");

        auto vDateTime1 = taxisInqVdatetime(taxisID1);
        auto vDateTime2 = taxisInqVdatetime(taxisID2);
        auto year1 = vDateTime1.date.year;
        auto year2 = vDateTime2.date.year;

        for (int iy = 0; iy < nyears; iy++)
          {
            if (iyears[iy] < year1 || iyears[iy] > year2)
              cdo_abort("Year %d out of bounds (first year %d; last year %d)!", iyears[iy], year1, year2);

            auto vDateTime3 = vDateTime1;
            vDateTime3.date.year = iyears[iy];
            taxisDefVdatetime(taxisID3, vDateTime3);
            cdo_def_timestep(streamIDs[iy], tsID);
          }

        for (int recID = 0; recID < nrecs; ++recID)
          {
            int varID, levelID;
            cdo_inq_record(streamID1, &varID, &levelID);
            cdo_inq_record(streamID2, &varID, &levelID);

            size_t nmiss1, nmiss2;
            cdo_read_record(streamID1, array1.data(), &nmiss1);
            cdo_read_record(streamID2, array2.data(), &nmiss2);

            auto gridsize = varList1[varID].gridsize;
            auto missval1 = varList1[varID].missval;
            auto missval2 = varList2[varID].missval;
            auto withMissval = (nmiss1 || nmiss2);

            for (int iy = 0; iy < nyears; iy++)
              {
                auto fac1 = ((double) year2 - iyears[iy]) / (year2 - year1);
                auto fac2 = ((double) iyears[iy] - year1) / (year2 - year1);

                auto nmiss3 = intlin_year(fac1, fac2, gridsize, array1, array2, array3, withMissval, missval1, missval2);

                cdo_def_record(streamIDs[iy], varID, levelID);
                cdo_write_record(streamIDs[iy], array3.data(), nmiss3);
              }
          }

        tsID++;
      }
  }

  void
  close()
  {
    for (int iy = 0; iy < nyears; iy++) cdo_stream_close(streamIDs[iy]);

    cdo_stream_close(streamID2);
    cdo_stream_close(streamID1);

    cdo_finish();
  }
};

void *
Intyear(void *process)
{
  ModuleIntYear intyear;
  intyear.init(process);
  intyear.run();
  intyear.close();
  return nullptr;
}
