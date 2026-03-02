/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#include <iostream>
#include <algorithm>
#include <vector>
#include <limits>

#ifdef HAVE_EXECINFO_H
#include <execinfo.h>
#endif

#ifdef HAVE_FEENABLEEXCEPT
#ifndef __USE_GNU
#define __USE_GNU  // gives us feenableexcept()
#endif
#endif
#include <cfenv>
#include <sys/stat.h>
#include <unistd.h> /* sysconf */
#include <cstring>
#include <csignal>

#include <cdi.h>

#include "cdo_getopt.h"
#include "cdo_rlimit.h"
#include <mpim_grid.h>
#include <griddes.h>
#include "cdo_default_values.h"
#include "param_conversion.h"
#include "progress.h"

#include "module_list.h"
#include "module_info.h"
#include "percentiles.h"
#include "util_wildcards.h"
#include "util_string.h"
#include "process_int.h"
#include "cdo_options.h"
#include "timer.h"
#include "commandline.h"
#include "mpmo_color.h"
#include "cdo_output.h"
#include "cdo_features.h"
#include "cdo_zaxis.h"
#include "compare.h"
#include "dmemory.h"
#include "table.h"
#include "datetime.h"
#include "remap_grid_cell_search.h"
#include "cdo_pthread.h"
#include "institution.h"
#include "cdo_apply.h"
#include "parser.h"

static ProcessManager g_processManager;

void
cdo_exit()
{
  g_processManager.kill_processes();
  exit(EXIT_FAILURE);
}

static int CDO_numThreads = 0;
static int timer_total;
static int CDO_netcdf_hdr_pad = 0;
static int CDO_Rusage = 0;
static bool applyDryRun = false;

#ifdef HIRLAM_EXTENSIONS
extern "C" void streamGrbDefDataScanningMode(int scanmode);
#endif

void set_point_search_method(const std::string &methodstr);

static void
cdo_stackframe()
{
#if defined HAVE_EXECINFO_H && defined HAVE_BACKTRACE
  void *callstack[32];
  auto frames = backtrace(callstack, 32);
  auto messages = backtrace_symbols(callstack, frames);

  fprintf(stderr, "[bt] Execution path:\n");
  if (messages)
    {
      for (int i = 0; i < frames; ++i) fprintf(stderr, "[bt] %s\n", messages[i]);
      free(messages);
    }
#endif
}

#ifdef HAVE_FEENABLEEXCEPT
static int
cdo_feenableexcept(int excepts)
{
  return feenableexcept(excepts);
}
#else
static int
cdo_feenableexcept(int excepts)
{
  static fenv_t fenv;
  if (std::fegetenv(&fenv)) return -1;

  (void) excepts;
  int oldExcepts = -1;  // previous masks
#if defined HAVE_FENV_T___CONTROL && defined HAVE_FENV_T___MXCSR
  unsigned newExcepts = ((unsigned) excepts) & FE_ALL_EXCEPT;
  oldExcepts = (int) (fenv.__control & FE_ALL_EXCEPT);

  // unmask
  fenv.__control &= ~newExcepts;
  fenv.__mxcsr &= ~(newExcepts << 7);
#endif

  return (std::fesetenv(&fenv) ? -1 : oldExcepts);
}
#endif

static void
cdo_signal_handler(int signo)
{
  if (signo == SIGFPE)
    {
      cdo_stackframe();
      cdo_abort("floating-point exception!");
    }
}

static void
cdo_set_digits(const char *const arg)
{
  char *ptr1 = 0;
  if (arg != 0 && (int) strlen(arg) > 0 && arg[0] != ',') Options::CDO_flt_digits = (int) strtol(arg, &ptr1, 10);

  if (Options::CDO_flt_digits < 1 || Options::CDO_flt_digits > 20)
    cdo_abort("Unreasonable value for float significant digits: %d", Options::CDO_flt_digits);

  if (ptr1 && *ptr1 == ',')
    {
      char *ptr2 = 0;
      Options::CDO_dbl_digits = (int) strtol(ptr1 + 1, &ptr2, 10);
      if (ptr2 == ptr1 + 1 || Options::CDO_dbl_digits < 1 || Options::CDO_dbl_digits > 20)
        cdo_abort("Unreasonable value for double significant digits: %d", Options::CDO_dbl_digits);
    }
}

static void
cdo_version()
{
  int filetypes[] = { CDI_FILETYPE_SRV, CDI_FILETYPE_EXT, CDI_FILETYPE_IEG,  CDI_FILETYPE_GRB, CDI_FILETYPE_GRB2,  CDI_FILETYPE_NC,
                      CDI_FILETYPE_NC2, CDI_FILETYPE_NC4, CDI_FILETYPE_NC4C, CDI_FILETYPE_NC5, CDI_FILETYPE_NCZARR };
  const char *typenames[] = { "srv", "ext", "ieg", "grb1", "grb2", "nc1", "nc2", "nc4", "nc4c", "nc5", "nczarr" };

  auto fp = stdout;
  fprintf(fp, "%s\n", cdo::Version);
#ifdef SYSTEM_TYPE
  fprintf(fp, "System: %s\n", SYSTEM_TYPE);
#endif
#ifdef CXX_COMPILER
  fprintf(fp, "CXX Compiler: %s\n", CXX_COMPILER);
#ifdef CXX_VERSION
  fprintf(fp, "CXX version : %s\n", CXX_VERSION);
#endif
#endif
#ifdef C_COMPILER
  fprintf(fp, "C Compiler: %s\n", C_COMPILER);
#ifdef C_VERSION
  fprintf(fp, "C version : %s\n", C_VERSION);
#endif
#endif
#ifdef F77_COMPILER
  fprintf(fp, "F77 Compiler: %s\n", F77_COMPILER);
#ifdef F77_VERSION
  fprintf(fp, "F77 version : %s\n", F77_VERSION);
#endif
#endif

  cdo_print_features();
  cdo_print_libraries();

#ifdef CDI_SIZE_TYPE
#define CDO_STRINGIFY(x) #x
#define CDO_TOSTRING(x) CDO_STRINGIFY(x)
  fprintf(fp, "CDI data types: SizeType=%s\n", CDO_TOSTRING(CDI_SIZE_TYPE));
#endif
  fprintf(fp, "CDI file types: ");
  set_text_color(fp, BRIGHT, GREEN);
  for (size_t i = 0; i < sizeof(filetypes) / sizeof(int); ++i)
    if (cdiHaveFiletype(filetypes[i])) fprintf(fp, "%s ", typenames[i]);
  reset_text_color(fp);
  fprintf(fp, "\n");

  cdiPrintVersion();
  fprintf(fp, "\n");
}

static void
cdo_variableInputs()
{
  set_text_color(stderr, BRIGHT, BLUE);
  fprintf(stderr, "#==============================================================================#\n");
  reset_text_color(stderr);
  fprintf(stderr, "    For operators with variable number of inputs:\n");
  fprintf(stderr, "    Brackets can be used for grouping input to the right operator.\n");
  reset_text_color(stderr);
  fprintf(stderr, "    example:\n");
  fprintf(stderr, "       -add -select,x=0 [ file1 -add -topo -file2 ] -merge [ file3 file4 ] out\n");
  set_text_color(stderr, BRIGHT, BLUE);
  fprintf(stderr, "#==============================================================================#\n");
  reset_text_color(stderr);
}

static void
cdo_init_is_tty()
{
  struct stat statbuf;
  fstat(0, &statbuf);
  if (S_ISCHR(statbuf.st_mode)) cdo::stdinIsTerminal = true;
  fstat(1, &statbuf);
  if (S_ISCHR(statbuf.st_mode))
    {
      cdo::stdoutIsTerminal = true;
      progress::stdoutIsTerminal = true;
    }
  fstat(2, &statbuf);
  if (S_ISCHR(statbuf.st_mode)) cdo::stderrIsTerminal = true;
}

static void
cdo_print_help(const char **help)
{
  if (!help)
    fprintf(stderr, "No help available for this operator!\n");
  else
    {
      size_t help_size = 0;
      while (help[help_size]) help_size++;
      for (size_t i = 0; i < help_size; ++i)
        {
          auto doPrint = !(help[i][0] == '\0' && help[i + 1][0] == ' ');
          if (doPrint)
            {
              if (color_enabled())
                {
                  if (cdo_cmpstr(help[i], "NAME") || cdo_cmpstr(help[i], "SYNOPSIS") || cdo_cmpstr(help[i], "DESCRIPTION")
                      || cdo_cmpstr(help[i], "OPERATORS") || cdo_cmpstr(help[i], "NAMELIST") || cdo_cmpstr(help[i], "PARAMETER")
                      || cdo_cmpstr(help[i], "ENVIRONMENT") || cdo_cmpstr(help[i], "NOTE") || cdo_cmpstr(help[i], "EXAMPLES"))
                    {
                      set_text_color(stdout, BRIGHT);
                      fprintf(stdout, "%s", help[i]);
                      reset_text_color(stdout);
                      fprintf(stdout, "\n");
                    }
                  else
                    fprintf(stdout, "%s\n", help[i]);
                }
              else { fprintf(stdout, "%s\n", help[i]); }
            }
        }
    }
}

#undef IsBigendian
#define IsBigendian() (u_byteorder.c[sizeof(long) - 1])

static void
set_default_datatype(const char *datatypestr)
{
  static const union
  {
    unsigned long l;
    unsigned char c[sizeof(long)];
  } u_byteorder = { 1 };
  enum
  {
    D_UINT,
    D_INT,
    D_FLT,
    D_CPX
  };
  int dtype = -1;

  auto datatype = tolower(*datatypestr);
  // clang-format off
  if      (datatype == 'i') { dtype = D_INT;  datatypestr++; }
  else if (datatype == 'u') { dtype = D_UINT; datatypestr++; }
  else if (datatype == 'f') { dtype = D_FLT;  datatypestr++; }
  else if (datatype == 'c') { dtype = D_CPX;  datatypestr++; }
  else if (datatype == 'p') {                 datatypestr++; }
  // clang-format on

  if (isdigit((int) *datatypestr))
    {
      auto nbits = atoi(datatypestr);
      datatypestr += 1;
      if (nbits >= 10) datatypestr += 1;

      if (dtype == -1)
        {
          if (nbits > 0 && nbits < 32)
            CdoDefault::DataType = nbits;
          else if (nbits == 32)
            CdoDefault::DataType = (CdoDefault::FileType == CDI_FILETYPE_GRB) ? CDI_DATATYPE_PACK32 : CDI_DATATYPE_FLT32;
          else if (nbits == 64)
            CdoDefault::DataType = CDI_DATATYPE_FLT64;
          else
            {
              cdo_warning("Unsupported number of bits %d!", nbits);
              cdo_warning("Use I8/I16/I32/F32/F64 for nc1/nc2/nc4/nc4c/nc5/nczarr; U8/U16/U32 for nc4/nc4c/nc5/nczarr; F32/F64 for "
                          "grb2/srv/ext/ieg; P1 - P24 for grb1/grb2.");
              cdo_abort("Unsupported number of bits!");
            }
        }
      else
        {
          // clang-format off
          if (dtype == D_INT)
            {
              if      (nbits ==  8) CdoDefault::DataType = CDI_DATATYPE_INT8;
              else if (nbits == 16) CdoDefault::DataType = CDI_DATATYPE_INT16;
              else if (nbits == 32) CdoDefault::DataType = CDI_DATATYPE_INT32;
              else cdo_abort("Unsupported number of bits = %d for datatype INT!", nbits);
            }
          else if (dtype == D_UINT)
            {
              if      (nbits ==  8) CdoDefault::DataType = CDI_DATATYPE_UINT8;
              else if (nbits == 16) CdoDefault::DataType = CDI_DATATYPE_UINT16;
              else if (nbits == 32) CdoDefault::DataType = CDI_DATATYPE_UINT32;
              else cdo_abort("Unsupported number of bits = %d for datatype UINT!", nbits);
            }
          else if (dtype == D_FLT)
            {
              if      (nbits == 32) CdoDefault::DataType = CDI_DATATYPE_FLT32;
              else if (nbits == 64) CdoDefault::DataType = CDI_DATATYPE_FLT64;
              else cdo_abort("Unsupported number of bits = %d for datatype FLT!", nbits);
            }
          else if (dtype == D_CPX)
            {
              if      (nbits == 32) CdoDefault::DataType = CDI_DATATYPE_CPX32;
              else if (nbits == 64) CdoDefault::DataType = CDI_DATATYPE_CPX64;
              else cdo_abort("Unsupported number of bits = %d for datatype CPX!", nbits);
            }
          // clang-format on
        }
    }

  if (*datatypestr != 0)
    {
      if (*datatypestr == 'l' || *datatypestr == 'L')
        {
          if (IsBigendian()) CdoDefault::Byteorder = CDI_LITTLEENDIAN;
        }
      else if (*datatypestr == 'b' || *datatypestr == 'B')
        {
          if (!IsBigendian()) CdoDefault::Byteorder = CDI_BIGENDIAN;
        }
      else { cdo_abort("Unsupported character in number of bytes: >%s< !", datatypestr); }
    }
}

static void
set_default_filetype(const std::string &filetypeString)
{

  if (filetypeString.size() > 0)
    {
      size_t len = 0;

      // clang-format off
      if      (cdo_cmpstrLenRhs(filetypeString, "grb2",   len)) CdoDefault::FileType = CDI_FILETYPE_GRB2;
      else if (cdo_cmpstrLenRhs(filetypeString, "grb1",   len)) CdoDefault::FileType = CDI_FILETYPE_GRB;
      else if (cdo_cmpstrLenRhs(filetypeString, "grb",    len)) CdoDefault::FileType = CDI_FILETYPE_GRB;
      else if (cdo_cmpstrLenRhs(filetypeString, "nc2",    len)) CdoDefault::FileType = CDI_FILETYPE_NC2;
      else if (cdo_cmpstrLenRhs(filetypeString, "nc4c",   len)) CdoDefault::FileType = CDI_FILETYPE_NC4C;
      else if (cdo_cmpstrLenRhs(filetypeString, "nc4",    len)) CdoDefault::FileType = CDI_FILETYPE_NC4;
      else if (cdo_cmpstrLenRhs(filetypeString, "nc5",    len)) CdoDefault::FileType = CDI_FILETYPE_NC5;
      else if (cdo_cmpstrLenRhs(filetypeString, "nc1",    len)) CdoDefault::FileType = CDI_FILETYPE_NC;
      else if (cdo_cmpstrLenRhs(filetypeString, "nczarr", len)) CdoDefault::FileType = CDI_FILETYPE_NCZARR;
      else if (cdo_cmpstrLenRhs(filetypeString, "nc",     len)) CdoDefault::FileType = CDI_FILETYPE_NC2;
      else if (cdo_cmpstrLenRhs(filetypeString, "srv",    len)) CdoDefault::FileType = CDI_FILETYPE_SRV;
      else if (cdo_cmpstrLenRhs(filetypeString, "ext",    len)) CdoDefault::FileType = CDI_FILETYPE_EXT;
      else if (cdo_cmpstrLenRhs(filetypeString, "ieg",    len)) CdoDefault::FileType = CDI_FILETYPE_IEG;
      else
        {
          cdo_warning("Unsupported filetype %s!", filetypeString);
          cdo_warning("Available filetypes: grb1/grb2/nc1/nc2/nc4/nc4c/nc5/nczarr/srv/ext/ieg");
          cdo_abort("Unsupported filetype %s!", filetypeString);
        }
      // clang-format on

      const char *ftstr = filetypeString.c_str() + len;

      if (CdoDefault::FileType != CDI_UNDEFID && *ftstr != 0)
        {
          if (*ftstr == '_') { set_default_datatype(++ftstr); }
          else
            {
              cdo_warning("Unexpected character >%c< in file type >%s<!", *ftstr, filetypeString);
              cdo_warning("Use format[_nbits] with:");
              cdo_warning("    format = grb1, grb2, nc1, nc2, nc4, nc4c, nc5, nczarr, srv, ext or ieg");
              cdo_warning("    nbits  = 32/64 for grb2/nc1/nc2/nc4/nc4c/nc5/nczarr/srv/ext/ieg; 1 - 24 for grb1/grb2");
              cdo_abort("Unexpected character in file type option!");
            }
        }
    }
}

#include <inttypes.h>

static auto
alignof_address(void *ptr) -> int
{
  auto n = reinterpret_cast<int64_t>(ptr);
  return (int) (n & (-n));
}

static auto
alignof_malloc_data(const std::vector<int> &tsize) -> int
{
  int align = (1 << 30);
  auto n = tsize.size();

  std::vector<double *> ptr(n);

  for (size_t i = 0; i < n; ++i)
    {
      ptr[i] = (double *) malloc(tsize[i] * sizeof(double));
      align = std::min(align, alignof_address(ptr[i]));
    }
  for (auto &p : ptr) free(p);

  return align;
}

static auto
alignof_vector_data(const std::vector<int> &tsize) -> int
{
  int align = 1 << 30;
  auto n = tsize.size();

  std::vector<std::vector<double>> ptr(n);

  for (size_t i = 0; i < n; ++i)
    {
      ptr[i].resize(tsize[i]);
      align = std::min(align, alignof_address(ptr[i].data()));
    }

  return align;
}

static void
define_filter(const std::string &argString)
{
  if (argString.size() > 0)
    {
      if (Options::filterId != 0) { cdo_abort("Filter already defined! Only one filter is allowed."); }

      auto filterArgs = split_with_seperator(argString, ',');
      Options::filterId = parameter_to_int(filterArgs[0]);
      if (Options::filterId <= 0) { cdo_abort("Undefined filter id: %d", Options::filterId); }
      for (size_t i = 1; i < filterArgs.size(); ++i)
        {
          if (filterArgs[i].size() > 0) Options::filterParams.push_back(parameter_to_int(filterArgs[i]));
        }
    }
  else { cdo_abort("Filter id missing!"); }
}

static void
define_compress(const std::string &argString)
{
  const char *arg = argString.c_str();
  size_t len = argString.size();

  if (argString == "szip")
    {
      Options::cdoCompType = CDI_COMPRESS_SZIP;
      Options::cdoCompLevel = 0;
    }
  else if (argString == "aec" || argString == "ccsds")
    {
      Options::cdoCompType = CDI_COMPRESS_AEC;
      Options::cdoCompLevel = 0;
    }
  else if (argString == "jpeg")
    {
      Options::cdoCompType = CDI_COMPRESS_JPEG;
      Options::cdoCompLevel = 0;
    }
  else if (strncmp(arg, "zip", 3) == 0)
    {
      Options::cdoCompType = CDI_COMPRESS_ZIP;
      Options::cdoCompLevel = (len == 5 && arg[3] == '_' && isdigit(arg[4])) ? atoi(&arg[4]) : 1;
    }
  else if (strncmp(arg, "zstd", 4) == 0)
    {
      if (Options::filterId != 0) { cdo_abort("Filter already defined! Only one filter is allowed."); }
      Options::filterId = 32015;
      int zstdLevel = (len >= 6 && len <= 7 && arg[4] == '_' && isdigit(arg[5])) ? atoi(&arg[5]) : 1;
      Options::filterParams.push_back(zstdLevel);
    }
  else { cdo_abort("Compression type '%s' unsupported!", arg); }
}

static void
define_chunktype(const std::string &arg)
{
  // clang-format off
  if      ("auto"  == arg) Options::cdoChunkType = CDI_CHUNK_AUTO;
  else if ("grid"  == arg) Options::cdoChunkType = CDI_CHUNK_GRID;
  else if ("lines" == arg) Options::cdoChunkType = CDI_CHUNK_LINES;
  else cdo_abort("Chunk type '%s' unsupported!", arg);
  // clang-format on
}

std::vector<std::string>
define_varnames(const char *const arg)
{

  std::string strArgs = std::string(arg);
  std::vector<std::string> newVarnames;

  char delim = ',';
  size_t previous = 0;
  size_t current = strArgs.find(delim);

  while (current != std::string::npos)
    {
      newVarnames.push_back(strArgs.substr(previous, current - previous));
      previous = current + 1;
      current = strArgs.find(delim, previous);
    }
  newVarnames.push_back(strArgs.substr(previous, current - previous));

  return newVarnames;
}

static void
get_env_vars()
{
  CLIOptions::envvar("CDO_TEST")
      ->add_effect([&](const std::string &envstr) { Options::test = parameter_to_bool(envstr); })
      ->describe_argument("true|false")
      ->add_default("false")
      ->add_help("'true' test new features [default: false].");

  CLIOptions::envvar("CDO_CORESIZE")
      ->add_effect([&](const std::string &envstr) { Options::coresize = parameter_to_long(envstr); })
      ->describe_argument("max. core dump size")
      ->add_help("The largest size (in bytes) core file that may be created.");

  CLIOptions::envvar("CDO_DOWNLOAD_PATH")
      ->add_effect([&](const std::string &downloadPath) { cdo::DownloadPath = downloadPath; })
      ->describe_argument("path")
      ->add_help("Path where CDO can store downloads.");

  CLIOptions::envvar("CDO_ICON_GRIDS")
      ->add_effect([&](const std::string &iconGrid) { cdo::IconGrids = iconGrid; })
      ->describe_argument("path")
      ->add_help("Root directory of the installed ICON grids (e.g. /pool/data/ICON).");

  CLIOptions::envvar("CDO_DISABLE_HISTORY")
      ->add_effect([&](const std::string &envstr) {
        if (parameter_to_bool(envstr) == true)
          {
            Options::CDO_Reset_History = true;
            Options::CDO_Append_History = false;
          }
      })
      ->describe_argument("true|false")
      ->add_help("MISSING HELP");

  CLIOptions::envvar("CDO_RESET_HISTORY")
      ->add_effect([&](const std::string &envstr) { Options::CDO_Reset_History = parameter_to_bool(envstr); })
      ->describe_argument("true|false")
      ->add_default("false")
      ->add_help("'true' resets the global history attribute [default: false].");

  CLIOptions::envvar("CDO_HISTORY_INFO")
      ->add_effect([&](const std::string &envstr) { Options::CDO_Append_History = parameter_to_bool(envstr); })
      ->describe_argument("true|false")
      ->add_default("true")
      ->add_help("'false' don't write information to the global history attribute [default: true].");

  cdo::File_Suffix[0] = 0;
  CLIOptions::envvar("CDO_FILE_SUFFIX")
      ->add_effect([&](const std::string &envstr) { strncat(cdo::File_Suffix, envstr.c_str(), sizeof(cdo::File_Suffix) - 1); })
      ->describe_argument("suffix")
      ->add_help("Default filename suffix.");

  CLIOptions::envvar("CDO_DISABLE_FILESUFFIX")
      ->add_effect([&]() { strcat(cdo::File_Suffix, "nullptr"); })
      ->describe_argument("true|false")
      ->add_help("MISSING HELP");

  CLIOptions::envvar("CDO_VERSION_INFO")
      ->add_effect([&](const std::string &envstr) { Options::VersionInfo = parameter_to_bool(envstr); })
      ->describe_argument("true|false")
      ->add_default("true")
      ->add_help("'false' disables the global NetCDF attribute CDO [default: true].");
}

static void
print_system_info()
{
  fprintf(stderr, "\n");
  fprintf(stderr, "CDO_Color           = %d\n", mpmo_get_color_mode());
  fprintf(stderr, "Options::CDO_Reset_History   = %d\n", Options::CDO_Reset_History);
  fprintf(stderr, "CDO_File_Suffix     = %s\n", cdo::File_Suffix);
  fprintf(stderr, "CdoDefault::FileType  = %d\n", CdoDefault::FileType);
  fprintf(stderr, "CdoDefault::DataType  = %d\n", CdoDefault::DataType);
  fprintf(stderr, "CdoDefault::Byteorder = %d\n", CdoDefault::Byteorder);
  fprintf(stderr, "CdoDefault::TableID   = %d\n", CdoDefault::TableID);
  fprintf(stderr, "\n");

  const char *envstr;
  envstr = getenv("HOSTTYPE");
  if (envstr) fprintf(stderr, "HOSTTYPE            = %s\n", envstr);
  envstr = getenv("VENDOR");
  if (envstr) fprintf(stderr, "VENDOR              = %s\n", envstr);
  envstr = getenv("OSTYPE");
  if (envstr) fprintf(stderr, "OSTYPE              = %s\n", envstr);
  envstr = getenv("MACHTYPE");
  if (envstr) fprintf(stderr, "MACHTYPE            = %s\n", envstr);
  fprintf(stderr, "\n");

#if defined(_ARCH_PWR6)
  fprintf(stderr, "Predefined: _ARCH_PWR6\n");
#elif defined(_ARCH_PWR7)
  fprintf(stderr, "Predefined: _ARCH_PWR7\n");
#endif

#if defined(__AVX2__)
  fprintf(stderr, "Predefined: __AVX2__\n");
#elif defined(__AVX__)
  fprintf(stderr, "Predefined: __AVX__\n");
#elif defined(__SSE4_2__)
  fprintf(stderr, "Predefined: __SSE4_2__\n");
#elif defined(__SSE4_1__)
  fprintf(stderr, "Predefined: __SSE4_1__\n");
#elif defined(__SSE3__)
  fprintf(stderr, "Predefined: __SSE3__\n");
#elif defined(__SSE2__)
  fprintf(stderr, "Predefined: __SSE2__\n");
#endif
  fprintf(stderr, "\n");

  fprintf(stderr, "sizeof(size_t)      = %zu\n", sizeof(size_t));
  {
    constexpr size_t megaByte = 1024 * 1024;
    std::vector<int> numElements = { 1, 3, 5, 9, 17, 33, 69, 121, 251, 510, 1025, 1 * megaByte };
    fprintf(stderr, "alignof malloc data = %d\n", alignof_malloc_data(numElements));
    fprintf(stderr, "alignof malloc big  = %d\n", alignof_malloc_data({ 8 * megaByte, 16 * megaByte, 32 * megaByte }));
    fprintf(stderr, "alignof vector data = %d\n", alignof_vector_data(numElements));
    fprintf(stderr, "alignof vector big  = %d\n", alignof_vector_data({ 8 * megaByte, 16 * megaByte, 32 * megaByte }));
  }
  fprintf(stderr, "\n");

#ifdef HAVE_MMAP
  fprintf(stderr, "HAVE_MMAP\n");
#endif
#ifdef HAVE_MEMORY_H
  fprintf(stderr, "HAVE_MEMORY_H\n");
#endif
  fprintf(stderr, "\n");

#ifdef _OPENACC
  fprintf(stderr, "OPENACC VERSION     = %d\n", _OPENACC);
#endif
  // OPENMP3:   201107
  // OPENMP4:   201307 gcc 4.9
  // OPENMP45:  201511
#ifdef _OPENMP
  fprintf(stderr, "OPENMP VERSION      = %d\n", _OPENMP);
#endif
  fprintf(stderr, "__cplusplus         = %ld\n", (long) __cplusplus);
#ifdef __GNUC__
  fprintf(stderr, "GNUC VERSION        = %d\n", __GNUC__);
#endif
#ifdef __GNUC_MINOR__
  fprintf(stderr, "GNUC MINOR          = %d\n", __GNUC_MINOR__);
#endif
#ifdef __ICC
  fprintf(stderr, "ICC VERSION         = %d\n", __ICC);
#endif
#ifdef __STDC__
  fprintf(stderr, "STD ANSI C          = %d\n", __STDC__);
#endif
#ifdef __STD_VERSION__
  fprintf(stderr, "STD VERSION         = %ld\n", __STD_VERSION__);
#endif
#ifdef __STDC_VERSION__
  fprintf(stderr, "STDC VERSION        = %ld\n", __STDC_VERSION__);
#endif
#ifdef __STD_HOSTED__
  fprintf(stderr, "STD HOSTED          = %d\n", __STD_HOSTED__);
#endif
#ifdef FLT_EVAL_METHOD
  fprintf(stderr, "FLT_EVAL_METHOD     = %d\n", FLT_EVAL_METHOD);
#endif
#ifdef FP_FAST_FMA
  fprintf(stderr, "FP_FAST_FMA         = defined\n");
#endif
#ifdef __FAST_MATH__
  fprintf(stderr, "__FAST_MATH__       = defined\n");
#endif
  fprintf(stderr, "\n");

#ifdef _SC_VERSION
  fprintf(stderr, "POSIX.1 VERSION     = %ld\n", sysconf(_SC_VERSION));
#endif
#ifdef _SC_ARG_MAX
  fprintf(stderr, "POSIX.1 ARG_MAX     = %ld\n", sysconf(_SC_ARG_MAX));
#endif
#ifdef _SC_CHILD_MAX
  fprintf(stderr, "POSIX.1 CHILD_MAX   = %ld\n", sysconf(_SC_CHILD_MAX));
#endif
#ifdef _SC_STREAM_MAX
  fprintf(stderr, "POSIX.1 STREAM_MAX  = %ld\n", sysconf(_SC_STREAM_MAX));
#endif
#ifdef _SC_OPEN_MAX
  fprintf(stderr, "POSIX.1 OPEN_MAX    = %ld\n", sysconf(_SC_OPEN_MAX));
#endif
#ifdef _SC_PAGESIZE
  fprintf(stderr, "POSIX.1 PAGESIZE    = %ld\n", sysconf(_SC_PAGESIZE));
#endif

  fprintf(stderr, "\n");

  cdo::print_rlimits();

  fprintf(stderr, "\n");
}

static void
cdo_set_options()
{
  if (cdo::dbg())
    {
      fprintf(stderr, "CMOR_Mode           = %d\n", Options::CMOR_Mode);
      fprintf(stderr, "CDO_netcdf_hdr_pad  = %d\n", CDO_netcdf_hdr_pad);
      fprintf(stderr, "\n");
    }

  if (Options::CMOR_Mode) cdiDefGlobal("CMOR_MODE", Options::CMOR_Mode);  // TODO maybe reposition into effect of "cmor"
  if (Options::CDO_Reduce_Dim) cdiDefGlobal("REDUCE_DIM", Options::CDO_Reduce_Dim);
  if (CDO_netcdf_hdr_pad > 0) cdiDefGlobal("NETCDF_HDR_PAD", CDO_netcdf_hdr_pad);
}

void
evaluate_color_options(const std::string &arg)
{
  // clang-format off
  if      ("all"  == arg) mpmo_color_set(All);
  else if ("auto" == arg) mpmo_color_set(Auto);
  else if ("no"   == arg) mpmo_color_set(No);
  else cdo_abort("Color option <%s> unknown. Known options: auto, all, no", Yellow(arg));
  // clang-format on
}

int
evaluate_except_options(const std::string &arg)
{
  int except = -1;
  // clang-format off
  if      (arg == "DIVBYZERO")  except = FE_DIVBYZERO;
  else if (arg == "INEXACT")    except = FE_INEXACT;
  else if (arg == "INVALID")    except = FE_INVALID;
  else if (arg == "OVERFLOW")   except = FE_OVERFLOW;
  else if (arg == "UNDERFLOW")  except = FE_UNDERFLOW;
  else if (arg == "ALL_EXCEPT") except = FE_ALL_EXCEPT;
  // clang-format on
  return except;
}

static void
cdo_rusage(void)
{
#if defined HAVE_SYS_RESOURCE_H && defined RUSAGE_SELF
  struct rusage ru;
  auto status = getrusage(RUSAGE_SELF, &ru);
  if (status == 0)
    {
      double ut = ru.ru_utime.tv_sec + 0.000001 * ru.ru_utime.tv_usec;
      double st = ru.ru_stime.tv_sec + 0.000001 * ru.ru_stime.tv_usec;

      fprintf(stderr, "  User time:     %.3f seconds\n", ut);
      fprintf(stderr, "  System time:   %.3f seconds\n", st);
      fprintf(stderr, "  Total time:    %.3f seconds\n", ut + st);
      fprintf(stderr, "  Memory usage:  %.2f MBytes\n", ru.ru_maxrss / (1024.0 * 1024.0));
      fprintf(stderr, "  Page reclaims: %5ld page%s\n", ru.ru_minflt, ADD_PLURAL(ru.ru_minflt));
      fprintf(stderr, "  Page faults:   %5ld page%s\n", ru.ru_majflt, ADD_PLURAL(ru.ru_majflt));
      fprintf(stderr, "  Swaps:         %5ld\n", ru.ru_nswap);
      fprintf(stderr, "  Disk read:     %5ld block%s\n", ru.ru_inblock, ADD_PLURAL(ru.ru_inblock));
      fprintf(stderr, "  Disk Write:    %5ld block%s\n", ru.ru_oublock, ADD_PLURAL(ru.ru_oublock));
    }
#endif
}
// clang-format on

#ifdef _OPENMP
static void
print_openmp_info()
{
  fprintf(stderr, "OMP num procs       = %d\n", omp_get_num_procs());
  fprintf(stderr, "OMP max threads     = %d\n", omp_get_max_threads());
  fprintf(stderr, "OMP num threads     = %d\n", omp_get_num_threads());
#ifndef HAVE_OPENMP3
  fprintf(stderr, "OMP thread limit    = %d\n", omp_get_thread_limit());
  omp_sched_t kind;
  int modifer;
  omp_get_schedule(&kind, &modifer);
  fprintf(stderr, "OMP schedule        = %d (1:static; 2:dynamic; 3:guided; 4:auto)\n", (int) kind);
#endif
#ifdef HAVE_OPENMP4
  fprintf(stderr, "OMP proc bind       = %d (0:false; 1:true; 2:master; 3:close; 4:spread)\n", (int) omp_get_proc_bind());
#ifndef __ICC
  fprintf(stderr, "OMP num devices     = %d\n", omp_get_num_devices());
#endif
#endif
}
#endif

static void
set_external_proj_func(void)
{
#ifdef HAVE_CDI_PROJ_FUNCS
  proj_lonlat_to_lcc_func = proj_lonlat_to_lcc;
  proj_lcc_to_lonlat_func = proj_lcc_to_lonlat;

  proj_lonlat_to_stere_func = proj_lonlat_to_stere;
  proj_stere_to_lonlat_func = proj_stere_to_lonlat;
#endif
}

static const char *
get_progname(char *string)
{
#ifdef _WIN32
  //  progname = strrchr(string, '\\');
  char *progname = " cdo";
#else
  char *progname = strrchr(string, '/');
#endif

  return (progname == nullptr) ? string : ++progname;
}

#ifdef HAVE_H5DONT_ATEXIT
extern "C" void H5dont_atexit(void);
#endif

static void
print_operator_attributes(const std::string &argument)
{
  ModListOptions local_modListOpt;
  local_modListOpt.parse_request(argument);
  operator_print_list(local_modListOpt);
}

static void
setup_openMP()
{
#ifdef _OPENMP
  if (CDO_numThreads <= 0) CDO_numThreads = 1;
  omp_set_num_threads(CDO_numThreads);

  Threading::ompNumThreads = omp_get_max_threads();
  if (omp_get_max_threads() > omp_get_num_procs())
    fprintf(stderr, "Warning: Number of OMP threads=%d is greater than number of Cores=%d!\n", omp_get_max_threads(),
            omp_get_num_procs());

  if (Threading::ompNumThreads < CDO_numThreads)
    fprintf(stderr, "Warning: omp_get_max_threads() returns %d!\n", Threading::ompNumThreads);

  if (cdo::dbg()) print_openmp_info();

  if (Options::cdoVerbose)
    {
      fprintf(stderr, " OpenMP:  num_procs=%d  max_threads=%d", omp_get_num_procs(), omp_get_max_threads());
#ifdef HAVE_OPENMP4
#ifndef __ICC
      fprintf(stderr, "  num_devices=%d", omp_get_num_devices());
#endif
#endif
      fprintf(stderr, "\n");
    }
#else
  if (CDO_numThreads > 1) fprintf(stderr, "Warning: Option -P failed, OpenMP support not compiled in!\n");
#endif
}

static void
cdo_print_debug_info()
{
  fprintf(stderr, "stdinIsTerminal:   %d\n", cdo::stdinIsTerminal);
  fprintf(stderr, "stdoutIsTerminal:  %d\n", cdo::stdoutIsTerminal);
  fprintf(stderr, "stderrIsTerminal:  %d\n", cdo::stderrIsTerminal);
  print_system_info();
  print_pthread_info();
}

static std::string
predefined_tables(int p_padding)
{
  const char *name;
  constexpr int id_padding = 4;
  int padding = p_padding + id_padding;
  int numTables = tableInqNumber();
  std::string tables = std::string("Predefined tables: ");
  for (int id = 0; id < numTables; id++)
    {
      if (id % 7 == 6) tables += "\n" + std::string(padding, ' ');
      if ((name = tableInqNamePtr(id))) tables += std::string(name);
      if (id < numTables - 1) tables += ",";
    }
  return tables;
}

static void
create_options_from_envvars()
{
  CLIOptions::option_from_envvar("CDO_DISABLE_FILESUFFIX");
  CLIOptions::option_from_envvar("CDO_DISABLE_HISTORY");
  CLIOptions::option_from_envvar("CDO_DOWNLOAD_PATH");
  CLIOptions::option_from_envvar("CDO_FILE_SUFFIX");
  CLIOptions::option_from_envvar("CDO_HISTORY_INFO");
  CLIOptions::option_from_envvar("CDO_ICON_GRIDS");
  CLIOptions::option_from_envvar("CDO_RESET_HISTORY");
  CLIOptions::option_from_envvar("CDO_TEST");
  CLIOptions::option_from_envvar("CDO_VERSION_INFO");
}

static void
setup_cli_options()
{
  CLIOptions::option("envvars")
      ->add_effect([&]() { CLIOptions::print_envvars = true; })
      ->aborts_program(true)
      ->add_help("Prints the environment variables of CDO.");

  CLIOptions::option("settings")
      ->add_effect([&]() { CLIOptions::print_settings = true; })
      ->aborts_program(true)
      ->add_help("Prints the settings of CDO.");

  CLIOptions::option("debug", "d")
      ->add_effect([&]() {
        unsigned cdoDebugLevel = 0;
        unsigned cdiDebugLevel = 0;
        cdo::parse_debug_arguments({ "1" }, cdoDebugLevel, cdiDebugLevel);

        cdiDebug(cdiDebugLevel);
        cdo::set_debug(cdoDebugLevel);
        cdo_version();
      })
      ->add_help("Pring all available debug messages");

  CLIOptions::option("scoped_debug", "D")
      ->describe_argument("comma seperated scopes")
      ->set_argument_optional(true)
      ->add_effect([&](const std::string &argument) {
        auto [success, tokens] = tokenize_comma_seperated_int_list(argument);
        if (argument.empty() || success == false)
          {
            std::cout << "No debug level given please choose: " << std::endl;
            print_debug_options();
            exit(EXIT_SUCCESS);
          }
        else
          {

            unsigned cdoDebugLevel = 0;
            unsigned cdiDebugLevel = 0;
            cdo::parse_debug_arguments(tokens, cdoDebugLevel, cdiDebugLevel);

            cdiDebug(cdiDebugLevel);
            cdo::set_debug(cdoDebugLevel);

            cdo_version();
          }
      })
      ->add_help("Multiple scopes suimultaneusly possible. Use option without arguments to get a list of possible scopes");

  CLIOptions::option("worker")
      ->describe_argument("num")
      ->add_effect([&](const std::string &argument) { Options::numStreamWorker = parameter_to_int(argument); })
      ->add_help("Number of worker to decode/decompress GRIB records.");

  CLIOptions::option("precision")
      ->describe_argument("float_digits[,double_digits]")
      ->add_effect([&](const std::string &argument) { cdo_set_digits(argument.c_str()); })
      ->add_help("Precision to use in displaying floating-point data (default: 7,15).");

  CLIOptions::option("percentile")
      ->describe_argument("method")
      ->add_effect([&](const std::string &argument) { percentile_set_method(argument); })
      ->add_help("Methods: nrank, nist, rtype8, <NumPy method (linear|lower|higher|nearest|...)>");

  CLIOptions::option("netcdf_hdr_pad")
      ->describe_argument("nbr")
      ->add_effect([&](const std::string &argument) {
        int netcdf_hdr_pad = parameter_to_bytes(argument);
        if (netcdf_hdr_pad >= 0) CDO_netcdf_hdr_pad = netcdf_hdr_pad;
      })
      ->add_help("Pad NetCDF output header with nbr bytes.");

  CLIOptions::option("use_fftw")
      ->describe_argument("true|false")
      ->add_effect([&](const std::string &argument) { Options::Use_FFTW = (int) parameter_to_bool(argument); })
      ->add_help("Sets fftw usage.");

  CLIOptions::option("cellsearchmethod")
      ->describe_argument("spherepart|latbins")
      ->add_effect([&](const std::string &argument) { set_cell_search_method(argument); })
      ->add_help("Sets the cell search method.");

  CLIOptions::option("config")
      ->describe_argument("all|all-json|<specific_feature_name>")
      ->add_effect([&](const std::string &argument) { cdo_print_config(argument); })
      ->aborts_program(true)
      ->add_help("Prints all features and the enabled status.", "Use option <all> to see explicit feature names.");

  CLIOptions::option("pointsearchmethod")
      ->set_internal(true)
      ->describe_argument("<full|kdtree|nanoflann|spherepart|latbins>")
      ->add_effect([&](const std::string &argument) { set_point_search_method(argument); })
      ->add_help("Sets the point search method.");

  CLIOptions::option("gridsearchradius")
      ->describe_argument("degrees[0..180]")
      ->add_effect([&](const std::string &argument) {
        extern double pointSearchRadius;
        auto fval = radius_str_to_deg(argument);
        if (fval < 0 || fval > 180) cdo_abort("%s=%g out of bounds (0-180 deg)!", "gridsearchradius", fval);
        pointSearchRadius = fval;
      })
      ->add_help("Sets the grid search radius (0-180 deg).");

  CLIOptions::option("remap_weights")
      ->describe_argument("0|1")
      ->add_effect([&](const std::string &argument) {
        auto intarg = parameter_to_int(argument);
        if (intarg != 0 && intarg != 1) cdo_abort("Unsupported value for option --remap_weights %d [0/1]", intarg);
        Options::REMAP_genweights = intarg;
      })
      ->add_help("Generate remap weights (default: 1).");

  CLIOptions::option("no_remap_weights")
      ->add_effect([&]() { Options::REMAP_genweights = 0; })
      ->add_help("Switch off generation of remap weights.");

  CLIOptions::option("enableexcept")
      ->describe_argument("except")
      ->add_effect([&](const std::string &argument) {
        auto except = evaluate_except_options(argument);
        if (except < 0) cdo_abort("option --%s: unsupported argument: %s", "enableexcept", argument);
        cdo_feenableexcept(except);
        if (signal(SIGFPE, cdo_signal_handler) == SIG_ERR) cdo_warning("can't catch SIGFPE!");
      })
      ->add_help("Set individual floating-point traps ", "(DIVBYZERO, INEXACT, INVALID, OVERFLOW, UNDERFLOW, ALL_EXCEPT)");

  CLIOptions::option("timestat_date")
      ->describe_argument("srcdate")
      ->add_effect([&](const std::string &argument) { set_timestat_date(argument); })
      ->add_help("Target timestamp (temporal statistics): ", "first, middle, midhigh or last source timestep.");

  CLIOptions::option("ignore_time_bounds")
      ->add_effect([&]() {
        extern bool CDO_Ignore_Time_Bounds;
        CDO_Ignore_Time_Bounds = true;
      })
      ->add_help("Ignores time bounds for time range statistics.");

  CLIOptions::option("use_time_bounds")
      ->add_effect([&]() {
        extern bool CDO_Use_Time_Bounds;
        CDO_Use_Time_Bounds = true;
      })
      ->add_help("Enables use of timebounds.");

  CLIOptions::option("cmor")->add_effect([&]() { Options::CMOR_Mode = 1; })->add_help("CMOR conform NetCDF output.");

  CLIOptions::option("reduce_dim")->add_effect([&]() { Options::CDO_Reduce_Dim = 1; })->add_help("Reduce NetCDF dimensions.");

  CLIOptions::option("float")
      ->add_effect([&]() { Options::CDO_Memtype = MemType::Float; })
      ->add_help("Using single precision floats for data in memory.");

  CLIOptions::option("single")
      ->add_effect([&]() { Options::CDO_Memtype = MemType::Float; })
      ->add_help("Using single precision floats for data in memory.");

  CLIOptions::option("double")
      ->add_effect([&]() { Options::CDO_Memtype = MemType::Double; })
      ->add_help("Using double precision floats for data in memory.");

  CLIOptions::option("rusage")->add_effect([&]() { CDO_Rusage = 1; })->add_help("Print information about resource utilization.");

  CLIOptions::option("attribs")
      ->describe_argument("arbitrary|filesOnly|onlyFirst|noOutput|obase")
      ->aborts_program(true)
      ->add_effect([&](const std::string &argument) { print_operator_attributes(argument); })
      ->add_help("Lists all operators with choosen features or the attributes of given operator(s)",
                 "operator name or a combination of [arbitrary,filesOnly,onlyFirst,noOutput,obase].");

  CLIOptions::option("operators")
      ->aborts_program(true)
      ->add_effect([&]() { print_operator_attributes(std::string()); })
      ->add_help("Prints list of operators.");

  CLIOptions::option("module_info")
      ->aborts_program(true)
      ->describe_argument("module name")
      ->add_effect([&](const std::string &argument) {
        auto names = get_module_operator_names(argument);
        if (names.empty())
          {
            std::string errstr = "Module " + argument + " not found\n";
            std::cerr << errstr;
          }
        else
          {

            std::string info_string = "\n" + argument + ":\n";
            for (const auto &name : names) { info_string += std::string(4, ' ') + name + "\n"; }
            std::cerr << info_string + "\n";
          }
      })
      ->add_help("Prints list of operators.");

  CLIOptions::option("operators_no_output")
      ->aborts_program(true)
      ->add_effect([&]() { print_operator_attributes("noOutput"); })
      ->add_help("Prints all operators which produce no output.");

  CLIOptions::option("pedantic")->add_effect([&]() { MpMO::enable_pedantic(true); })->add_help("Warnings count as errors.");

  CLIOptions::option("color", "C")
      ->describe_argument("auto|no|all")
      ->add_effect([&](const std::string &argument) { evaluate_color_options(argument); })
      ->add_help("Set behaviour of colorized output messages.");

  CLIOptions::option("eccodes")
      ->add_effect([&]() { cdiDefGlobal("ECCODES_GRIB1", true); })
      ->add_help("Use ecCodes to decode/encode GRIB1 messages.");

  CLIOptions::option("format", "f")
      ->describe_argument("grb1|grb2|nc1|nc2|nc4|nc4c|nc5|nczarr|srv|ext|ieg")
      ->add_effect([&](const std::string &argument) { set_default_filetype(argument); })
      ->add_help("Format of the output file.");

  CLIOptions::option("help", "h")
      ->describe_argument("operator")
      ->set_argument_optional(true)
      ->add_effect([&](const std::string &argument) {
        if (argument.empty())
          CLIOptions::usage();
        else
          cdo_print_help(operator_help(argument));
      })
      ->aborts_program(true)
      ->add_help("Shows either help information for the given operator or the usage of CDO.");

  CLIOptions::option("history")
      ->add_effect([&]() { Options::CDO_Append_History = 1; })
      ->add_help("Do append to NetCDF \"history\" global attribute.");

  CLIOptions::option("no_history")
      ->add_effect([&]() { Options::CDO_Append_History = 0; })
      ->add_help("Do not append to NetCDF \"history\" global attribute.");

  CLIOptions::option("version", "V")
      ->add_effect([&]() { cdo_version(); })
      ->aborts_program(true)
      ->add_help("Print the version number.");

  CLIOptions::option("dryrun", "A")->add_effect([&]() { applyDryRun = true; })->add_help("Dry run that shows processed CDO call.");

  CLIOptions::option("absolute_taxis", "a")
      ->add_effect([&]() {
        if (CdoDefault::TaxisType == TAXIS_RELATIVE)
          cdo_abort("option --%s: can't be combined with option --%s", "absolute_taxis (-a)", "relative_taxis (-r)");
        CdoDefault::TaxisType = TAXIS_ABSOLUTE;
      })
      ->add_help("Generate an absolute time axis.");

  // clang-format off
  CLIOptions::option("default_datatype", "b")
      ->describe_argument("nbits")
      ->add_effect([&](const std::string &argument) { set_default_datatype(argument.c_str()); })
      ->add_help("Set the number of bits for the output precision",
                 "    I8|I16|I32|F32|F64     for nc1,nc2,nc4,nc4c,nc5,nczarr;",
                 "    U8|U16|U32             for nc4,nc4c,nc5;",
                 "    F32|F64                for grb2,srv,ext,ieg;",
                 "    P1 - P24               for grb1,grb2");
  // clang-format on

  CLIOptions::option("check_data_range", "c")
      ->add_effect([&]() { Options::CheckDatarange = true; })
      ->add_help("Enables checks for data overflow.");

  CLIOptions::option("grid", "g")
      ->describe_argument("grid")
      ->add_effect([&](const std::string &argument) { cdo_set_grids(argument.c_str()); })
      ->add_help("Set default grid name or file. Available grids: ",
                 "F<XXX>, t<RES>, tl<RES>, r<NX>x<NY>, global_<DXY>, zonal_<DY>, gme<NI>, lon=<LON>/lat=<LAT>");

  CLIOptions::option("institution", "i")
      ->describe_argument("institute_name")
      ->add_effect([&](const std::string &argument) { define_institution(argument.c_str()); })
      ->add_help("Sets institution name.");

  CLIOptions::option("chunktype", "k")
      ->describe_argument("auto|grid|lines")
      ->add_effect([&](const std::string &argument) { define_chunktype(argument); })
      ->add_help("NetCDF4 chunk type: auto, grid or lines.");

  CLIOptions::option("chunksize")
      ->describe_argument("size")
      ->add_effect([&](const std::string &argument) {
        int chunkSize = parameter_to_bytes(argument);
        if (chunkSize >= 0) Options::cdoChunkSize = chunkSize;
      })
      ->add_help("NetCDF4 chunk size.");

  CLIOptions::option("lock_io", "L")->add_effect([&]() { Threading::cdoLockIO = true; })->add_help("Lock IO (sequential access).");

  CLIOptions::option("zaxis", "l")
      ->describe_argument("zaxis")
      ->add_effect([&](const std::string &argument) { define_zaxis(argument.c_str()); })
      ->add_help("Set default zaxis name or file.");

  CLIOptions::option("set_missval", "m")
      ->describe_argument("missval")
      ->add_effect([&](const std::string &argument) { cdiDefMissval(atof(argument.c_str())); })
      ->add_help("Set the missing value of non NetCDF files (default: " + get_scientific(cdiInqMissval()) + ").");

  CLIOptions::option("has_missval", "M")
      ->add_effect([&]() { cdiDefGlobal("HAVE_MISSVAL", true); })
      ->add_help("Set HAS_MISSVAL to true.");

  CLIOptions::option("varnames", "n")
      ->set_internal(true)
      ->describe_argument("<varname| file>")
      ->add_effect([&](const std::string &argument) { Options::cdoVarnames = define_varnames(argument.c_str()); })
      ->add_help("Set default varnames or file.");

  CLIOptions::option("overwrite", "O")
      ->add_effect([&]() { Options::cdoOverwriteMode = true; })
      ->add_help("Overwrite existing output file, if checked.");

  CLIOptions::option("num_threads", "P")
      ->describe_argument("nthreads")
      ->add_effect([&](const std::string &argument) { CDO_numThreads = parameter_to_int(argument); })
      ->add_help("Set number of OpenMP threads.");

  CLIOptions::option("parrallel_read", "p")
      ->set_internal(true)
      ->add_effect([&]() {
        Options::CDO_Parallel_Read = true;
        Options::CDO_task = true;
      })
      ->add_help("Enables parallel read.");

  CLIOptions::option("sortname", "Q")
      ->add_effect([&]() { cdiDefGlobal("SORTNAME", true); })
      ->add_help("Alphanumeric sorting of NetCDF parameter names.");

  CLIOptions::option("seed")
      ->describe_argument("seed")
      ->add_effect([&](const std::string &argument) {
        int intarg = parameter_to_int(argument);
        if (intarg < 0) cdo_abort("Unsupported value for option --seed %d [>=0]", intarg);
        Options::Random_Seed = intarg;
      })
      ->add_help("Seed for a new sequence of pseudo-random numbers. <seed> must be >= 0");

  CLIOptions::option("regular", "R")
      ->add_effect([&]() {
        Options::cdoRegulargrid = true;
        cdiDefGlobal("REGULARGRID", true);
      })
      ->add_help("Convert GRIB1 data from global reduced to regular Gaussian grid (cgribex only).");

  CLIOptions::option("relative_taxis", "r")
      ->add_effect([&]() {
        if (CdoDefault::TaxisType == TAXIS_ABSOLUTE)
          cdo_abort("option --%s: can't be combined with option --%s", "relative_taxis (-r)", "absolute_taxis (-a)");
        CdoDefault::TaxisType = TAXIS_RELATIVE;
      })
      ->add_help("Generate a relative time axis.");

  CLIOptions::option("cdo_diagnostic", "S")
      ->add_effect([&]() { Options::cdoDiag = true; })
      ->add_help("Create an extra output stream for the module TIMSTAT. This stream",
                 "contains the number of non missing values for each output period.");

  CLIOptions::option("silent", "s")
      ->add_effect([&]() {
        Options::silentMode = true;
        MpMO::enable_silent_mode(Options::silentMode);
        progress::silentMode = true;
      })
      ->add_help("Silent mode.");

  CLIOptions::option("timer", "T")->add_effect([&]() { Options::Timer = true; })->add_help("Enable timer.");

  CLIOptions::option("table", "t")
      ->describe_argument("codetab")
      ->add_effect([&](const std::string &argument) { CdoDefault::TableID = cdo::define_table(argument); })
      ->add_help("Set GRIB1 default parameter code table name or file (cgribex only).", predefined_tables(CLIOptions::padding));

  CLIOptions::option("interactive", "u")
      ->add_effect([&]() { Options::cdoInteractive = true; })
      ->add_help("Enable CDO interactive mode.");

  CLIOptions::option("verbose", "v")
      ->add_effect([&]() {
        Options::cdoVerbose = true;
        MpMO::enable_verbose(true);
        CLIOptions::print_envvars = true;
        gridEnableVerbose(Options::cdoVerbose);
      })
      ->add_help("Print extra details for some operators.");

  CLIOptions::option("disable_warnings", "w")
      ->add_effect([&]() {  // disable warning messages
        MpMO::enable_warnings(false);
        extern int _Verbose;  // CDI Warnings
        _Verbose = 0;
      })
      ->add_help("Disable warning messages.");

  CLIOptions::option("par_io", "X")
      ->set_internal(true)
      ->add_effect([&]() {
        Options::cdoParIO = true;  // multi threaded I/O
      })
      ->add_help("Enables multithreaded I/O.");

  CLIOptions::option("compress", "Z")
      ->add_effect([&]() { Options::cdoCompress = true; })
      ->add_help("Enables compression. Default = SZIP");

  CLIOptions::option("filter", "F")
      ->describe_argument("filterId,params")
      ->add_effect([&](const std::string &argument) { define_filter(argument); })
      ->add_help("NetCDF4/HDF5 filter description.");

  CLIOptions::option("compression_type", "z")
      ->describe_argument("aec|jpeg|zip[_1-9]|zstd[1-19]")
      ->add_effect([&](const std::string &argument) { define_compress(argument); })
      ->add_help("aec         AEC compression of GRIB2 records", "jpeg        JPEG compression of GRIB2 records",
                 "zip[_1-9]   Deflate compression of NetCDF4 variables", "zstd[_1-19] Zstandard compression of NetCDF4 variables");

  CLIOptions::option("nsb")
      ->set_internal(true)
      ->describe_argument("1-23")
      ->add_effect([&](const std::string &argument) { Options::nsb = parameter_to_int(argument); })
      ->add_help("Number of significant bits used for bit-rounding.");

  CLIOptions::option("show_available_options")
      ->set_internal(true)
      ->aborts_program(true)
      ->add_effect([&]() { CLIOptions::print_available_options(); })
      ->add_help("Shows all available optins and prints all shortforms, only internal use for testing.");

  CLIOptions::option("argument_groups")->aborts_program(true)->add_effect([&]() { cdo_variableInputs(); });
  CLIOptions::option("sortparam")->add_effect([]() { cdiDefGlobal("SORTPARAM", true); });

#ifdef HIRLAM_EXTENSIONS
  CLIOptions::option("Dkext")
      ->describe_argument("debLev")
      ->add_effect([&](const std::string &argument) {
        auto extDebugVal = parameter_to_int(argument);
        if (extDebugVal > 0)
          {
            extern int cdiDebugExt;
            cdoDebugExt = extDebugVal;
            cdiDebugExt = extDebugVal;
          }
      })
      ->add_help("Setting debugLevel for extensions.");

  CLIOptions::option("outputGribDataScanningMode")
      ->describe_argument("mode")
      ->add_effect([&](const std::string &argument) {
        auto scanningModeValue = parameter_to_int(argument);
        if (cdoDebugExt) printf("scanningModeValue=%d\n", scanningModeValue);

        if ((scanningModeValue == 0) || (scanningModeValue == 64) || (scanningModeValue == 96))
          {
            streamGrbDefDataScanningMode(scanningModeValue);  // -1: not used; allowed modes: <0,
                                                              // 64, 96>; Default is 64
          }
        else
          {
            cdo_warning("Warning: %d not in allowed modes: <0, 64, 96>; Using default: 64\n", scanningModeValue);
            streamGrbDefDataScanningMode(64);
          }
      })
      ->add_help("Setting grib scanning mode for data in output file <0, 64, 96>.", "Default is 64");
#endif  // HIRLAM_EXTENSIONS
}



int
main(int argc, char *argv[])
{
  cdo::set_exit_function(cdo_exit);
  cdo::set_context_function(process_inq_prompt);
  progress::set_context_function(process_inq_prompt);

  mpmo_color_set(Auto);

  cdo_init_is_tty();
  CLIOptions::set_tty_status(cdo::stderrIsTerminal);

  memExitOnError();

  Options::CDO_Reduce_Dim = 0;

  // mallopt(M_MMAP_MAX, 0);

  set_command_line(argc, argv);

  cdo::progname = get_progname(argv[0]);

  get_env_vars();
  create_options_from_envvars();
  CLIOptions::get_env_vars();

  setup_cli_options();

  auto CDO_optind = CLIOptions::parse(std::vector<std::string>(argv, argv + argc));

  if (CDO_optind == CLIOptions::ABORT_REQUESTED)
    exit(EXIT_FAILURE);
  else if (CDO_optind == CLIOptions::EXIT_REQUESTED)
    exit(EXIT_SUCCESS);

  if (CDO_optind >= argc)
    {
      fprintf(stderr, "\nNo operator given!\n\n");
      CLIOptions::usage();
      exit(EXIT_FAILURE);
    }
  else
    {
      cdo_set_options();
      set_external_proj_func();
      cdo::set_stacksize(67108864);  // 64MB
      cdo::set_coresize(Options::coresize);
      setup_openMP();

      if (cdo::dbg()) cdo_print_debug_info();

      std::vector<std::string> new_argv(&argv[CDO_optind], argv + argc);

      new_argv = expand_wild_cards(new_argv);
      if (applyDryRun == true)
        {
          std::cerr << argv_to_string(new_argv) << std::endl;
          exit(applyDryRun ? 0 : -1);
        }

      if (CdoDefault::TableID != CDI_UNDEFID) cdo_def_table_id(CdoDefault::TableID);

      timer_total = timer_new("total");
      timer_read = timer_new("read");
      timer_write = timer_new("write");

#ifdef HAVE_H5DONT_ATEXIT
      H5dont_atexit();  // don't call H5close on exit
#endif
#ifdef CUSTOM_MODULES
      load_custom_modules("custom_modules");
      close_library_handles();
#endif

      auto processStructure = Parser::parse(new_argv, process_inq_prompt);
      g_processManager.buildProcessTree(processStructure);
      timer_start(timer_total);
      g_processManager.run_processes();
      timer_stop(timer_total);
      g_processManager.clear_processes();

      if (Options::Timer) timer_report();
    }

  if (CDO_Rusage) cdo_rusage();

  return Options::cdoExitStatus;
}
