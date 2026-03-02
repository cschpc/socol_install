#ifndef CDO_OPTIONS_H
#define CDO_OPTIONS_H

#include <vector>
#include <string>

#include "config.h"
#ifndef VERSION
#define VERSION "0.0.1"
#endif

namespace cdo
{
extern const char *progname;
extern char File_Suffix[32];
extern const char *Version;
extern std::string DownloadPath;
extern std::string IconGrids;
extern bool stdinIsTerminal;
extern bool stdoutIsTerminal;
extern bool stderrIsTerminal;
}  // namespace cdo

enum class MemType
{
  Native,
  Float,
  Double
};

namespace Options
{
extern long coresize;
extern int numStreamWorker;
extern int nsb;  // Number of significant bits
extern bool benchmark;
extern bool silentMode;
extern bool test;

extern int filterId;
extern std::vector<int> filterParams;

extern bool cdoCompress;
extern int cdoCompType;
extern int cdoCompLevel;
extern bool cdoInteractive;
extern bool cdoVerbose;
extern bool cdoProcessInfo;
extern int cdoExitStatus;
extern bool Timer;

extern bool CheckDatarange;

extern int CDO_flt_digits;
extern int CDO_dbl_digits;

extern bool Use_FFTW;
extern bool VersionInfo;
extern int CMOR_Mode;

extern bool cdoDiag;

extern MemType CDO_Memtype;

extern bool CDO_Parallel_Read;
extern bool CDO_task;

extern int CDO_Reduce_Dim;
extern int CDO_Append_History;
extern bool CDO_Reset_History;

extern unsigned Random_Seed;

extern int cdoChunkType;
extern int cdoChunkSize;
extern bool cdoOverwriteMode;
extern bool cdoParIO;
extern bool cdoRegulargrid;
size_t cdo_num_varnames();
extern std::vector<std::string> cdoVarnames;

extern bool REMAP_genweights;

}  // namespace Options

namespace Threading
{
extern int ompNumThreads;
extern bool cdoLockIO;
}  // namespace Threading

const char *cdo_comment(void);

void set_compression(int fileID, int filetype);

#endif
