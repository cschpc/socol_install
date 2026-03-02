/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include "config.h"

#include <cdi.h>
#include "cdo_options.h"
#include "cdo_output.h"

#include <cstring>

namespace cdo
{
const char *progname;
char File_Suffix[32];
const char *Version = "Climate Data Operators version " VERSION " (https://mpimet.mpg.de/cdo)";
std::string DownloadPath;
std::string IconGrids;
bool stdinIsTerminal = false;
bool stdoutIsTerminal = false;
bool stderrIsTerminal = false;
}  // namespace cdo

namespace Options
{
long coresize = 0;
int numStreamWorker = 0;
int nsb = 0;  // Number of significant bits
bool benchmark = false;
bool silentMode = false;
bool test = false;

// NetCDF4/HDF5 filter
int filterId = 0;
std::vector<int> filterParams;

bool cdoCompress = false;
int cdoCompType = CDI_COMPRESS_NONE;
int cdoCompLevel = 0;
bool cdoInteractive = false;
bool cdoVerbose = false;
int cdoExitStatus = 0;
bool Timer = false;

bool CheckDatarange = false;

int CDO_flt_digits = 7;   // TODO:rename
int CDO_dbl_digits = 15;  // TODO:rename

bool Use_FFTW = true;
bool VersionInfo = true;
int CMOR_Mode = false;

bool cdoDiag = false;

MemType CDO_Memtype(MemType::Native);
bool CDO_Parallel_Read = false;

int CDO_Reduce_Dim = false;
int CDO_Append_History = true;
bool CDO_Reset_History = false;
bool CDO_task = false;

unsigned Random_Seed = 1;

int cdoChunkType = CDI_UNDEFID;
int cdoChunkSize = CDI_UNDEFID;
bool cdoOverwriteMode = false;
bool cdoParIO = false;
bool cdoRegulargrid = false;
std::vector<std::string> cdoVarnames;
size_t
cdo_num_varnames()
{
  return cdoVarnames.size();
}

bool REMAP_genweights = true;

const char *cdoExpName = nullptr;
}  // namespace Options

namespace Threading
{
int ompNumThreads = 1;
bool cdoLockIO = false;
}  // namespace Threading

const char *
cdo_comment(void)
{
  return cdo::Version;
}

static bool
filetype_has_szip(int filetype)
{
  return (filetype == CDI_FILETYPE_GRB || filetype == CDI_FILETYPE_GRB2 || filetype == CDI_FILETYPE_NC4
          || filetype != CDI_FILETYPE_NC4C);
}

static bool
filetype_has_zip(int filetype)
{
  return (filetype == CDI_FILETYPE_NC4 || filetype == CDI_FILETYPE_NC4C || filetype == CDI_FILETYPE_NCZARR);
}

void
set_compression(int streamID, int filetype)
{
  if (Options::cdoCompress)
    {
      if (filetype == CDI_FILETYPE_GRB || filetype == CDI_FILETYPE_GRB2)
        {
          Options::cdoCompType = CDI_COMPRESS_SZIP;
          Options::cdoCompLevel = 0;
        }
      else if (filetype == CDI_FILETYPE_NC4 || filetype == CDI_FILETYPE_NC4C || filetype == CDI_FILETYPE_NCZARR)
        {
          Options::cdoCompType = CDI_COMPRESS_ZIP;
          Options::cdoCompLevel = 1;
        }
    }

  if (Options::cdoCompType != CDI_COMPRESS_NONE)
    {
      streamDefCompType(streamID, Options::cdoCompType);
      streamDefCompLevel(streamID, Options::cdoCompLevel);

      if (Options::cdoCompType == CDI_COMPRESS_SZIP && !filetype_has_szip(filetype))
        cdo_warning("SZIP compression not available for non GRIB/NetCDF4 data!");

      if (Options::cdoCompType == CDI_COMPRESS_JPEG && filetype != CDI_FILETYPE_GRB2)
        cdo_warning("JPEG compression not available for non GRIB2 data!");

      if (Options::cdoCompType == CDI_COMPRESS_ZIP && !filetype_has_zip(filetype))
        cdo_warning("Deflate compression not available for non NetCDF4 data!");
    }

  if (Options::filterId != 0)
    {
      streamDefFilter(streamID, Options::filterId, (int)Options::filterParams.size(), Options::filterParams.data());
    }
}
