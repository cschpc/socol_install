/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include "cdo_options.h"
#include "cdo_output.h"
#include "compare.h"
#include "cdi_lockedIO.h"
#include "varray.h"

int
cdo_read_timestepmask(const char *maskfile, std::vector<bool> &imask)
{
  auto streamID = stream_open_read_locked(maskfile);
  const auto vlistID = streamInqVlist(streamID);

  const auto nvars = vlistNvars(vlistID);
  if (nvars > 1) cdo_abort("timestepmask %s contains more than one variable!", maskfile);

  const auto nlev = zaxisInqSize(vlistInqVarZaxis(vlistID, 0));
  if (nlev > 1) cdo_abort("timestepmask %s has more than one level!", maskfile);

  const auto gridsize = gridInqSize(vlistInqVarGrid(vlistID, 0));
  if (gridsize > 1) cdo_abort("timestepmask %s has more than one gridpoint!", maskfile);

  auto nts = vlistNtsteps(vlistID);
  if (nts == -1)
    {
      nts = 0;
      while (streamInqTimestep(streamID, nts)) nts++;

      if (Options::cdoVerbose) cdo_print("%s: counted %i timeSteps in %s", __func__, nts, maskfile);

      streamClose(streamID);
      streamID = stream_open_read_locked(maskfile);
    }
  else if (Options::cdoVerbose) { cdo_print("%s: found %i timeSteps in %s", __func__, nts, maskfile); }

  int n = nts;
  imask.resize(nts);

  int tsID = 0;
  while (true)
    {
      auto nrecs = streamInqTimestep(streamID, tsID);
      if (nrecs == 0) break;

      if (nrecs != 1) cdo_abort("Internal error; unexprected number of records!");

      int varID, levelID;
      size_t nmiss;
      double value;
      streamInqRecord(streamID, &varID, &levelID);
      streamReadRecord(streamID, &value, &nmiss);

      imask[tsID] = !(nmiss || IS_EQUAL(value, 0));

      tsID++;
    }

  streamClose(streamID);

  return n;
}

static void
read_one_field(const char *text, const char *filename, Varray<double> &array)
{
  const auto streamID = stream_open_read_locked(filename);
  const auto vlistID = streamInqVlist(streamID);

  const auto nvars = vlistNvars(vlistID);
  if (nvars > 1) cdo_abort("%s file %s contains more than one variable!", text, filename);

  const auto nlev = zaxisInqSize(vlistInqVarZaxis(vlistID, 0));
  if (nlev > 1) cdo_abort("%s file %s has more than one level!", text, filename);

  const auto nrecs = streamInqTimestep(streamID, 0);
  if (nrecs != 1) cdo_abort("%s file %s contains more than one field!", text, filename);

  const auto gridsize = gridInqSize(vlistInqVarGrid(vlistID, 0));
  array.resize(gridsize);

  int varID, levelID;
  size_t nmiss;
  streamInqRecord(streamID, &varID, &levelID);
  streamReadRecord(streamID, array.data(), &nmiss);
  streamClose(streamID);
}

size_t
cdo_read_mask(const char *maskfile, std::vector<bool> &imask)
{
  Varray<double> array;
  read_one_field("Mask", maskfile, array);

  const auto gridsize = array.size();
  imask.resize(gridsize);

  for (size_t i = 0; i < gridsize; ++i) imask[i] = IS_NOT_EQUAL(array[i], 0);

  return gridsize;
}

size_t
cdo_read_index(const char *indexfile, std::vector<int> &index)
{
  Varray<double> array;
  read_one_field("Index", indexfile, array);

  const auto gridsize = array.size();
  index.resize(gridsize);

  for (size_t i = 0; i < gridsize; ++i) index[i] = (int) std::lround(array[i]) - 1;

  return gridsize;
}
