/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef _OPENMP
#include <omp.h>  // omp_get_num_procs()
#endif

#include <cdi.h>

#ifdef HAVE_NETCDF_META_H
#include <netcdf_meta.h>
#endif

#ifdef HAVE_HDF5_H
#include <hdf5.h>
#endif

#ifdef HAVE_ZLIB_H
#include <zlib.h>
#endif

#ifdef HAVE_LIBXML2
#include <libxml/xmlversion.h>
#endif

#ifdef HAVE_CURL_CURL_H
#include <curl/curl.h>
#endif

#ifdef HAVE_PROJ_H
#include <proj.h>
#endif

#ifdef HAVE_LIBCMOR
extern "C"
{
#include "cmor.h"
}
#endif

#include <stdio.h>

#include "process_int.h"
#include "cdo_features.h"
#include "cdo_rlimit.h"
#include "cimdOmp.h"
#include "lib/yac/yac_version.h"

#include <thread>  // std::thread::hardware_concurrency()

extern "C" size_t getMemorySize(void);

void
cdo_print_features(void)
{
  constexpr size_t gigaByte = 1024 * 1024 * 1024;
  auto fp = stdout;
  fprintf(fp, "Features: ");
  size_t rssCur = (size_t)cdo::get_rss_cur() / gigaByte;
  size_t memorySize = getMemorySize() / gigaByte;
  if (rssCur > 0 && rssCur < memorySize) fprintf(fp, "%zu/", rssCur);
  if (memorySize > 0) fprintf(fp, "%zuGB ", memorySize);
  auto concurrentThreads = std::thread::hardware_concurrency();
#ifdef _OPENMP
  unsigned numProcs = omp_get_num_procs();
  if (numProcs < concurrentThreads) fprintf(fp, "%u/", numProcs);
#endif
  fprintf(fp, "%uthreads", concurrentThreads);
  fprintf(fp, " c++%d", (int) ((__cplusplus - 200000) / 100));
#ifdef _OPENMP
  fprintf(fp, " OpenMP");
#if defined(HAVE_OPENMP45)
  fprintf(fp, "45");
#elif defined(HAVE_OPENMP4)
  fprintf(fp, "4");
#elif defined(HAVE_OPENMP3)
  fprintf(fp, "3");
#endif
#endif
#ifdef HAVE_CF_INTERFACE
  fprintf(fp, " Fortran");
#endif
#ifdef HAVE_LIBPTHREAD
  fprintf(fp, " pthreads");
#endif
#ifdef HAVE_LIBHDF5
  fprintf(fp, " HDF5");
#endif
#ifdef HAVE_NETCDF4
  fprintf(fp, " NC4");
#ifdef HAVE_NC4HDF5
  fprintf(fp, "/HDF5");
#ifdef HAVE_NC4HDF5_THREADSAFE
  fprintf(fp, "/threadsafe");
#endif
#endif
#endif
#ifdef HAVE_LIBNC_DAP
  fprintf(fp, " OPeNDAP");
#endif
#ifdef HAVE_LIBSZ
  fprintf(fp, " sz");
#endif
/*
#ifdef HAVE_LIBZ
fprintf(fp, " z");
#endif
*/
#ifdef HAVE_LIBUDUNITS2
  fprintf(fp, " udunits2");
#endif
#ifdef HAVE_LIBPROJ
  fprintf(fp, " proj");
#endif
#ifdef HAVE_LIBXML2
  fprintf(fp, " xml2");
#endif
#ifdef HAVE_LIBMAGICS
  fprintf(fp, " magics");
#endif
#ifdef HAVE_LIBDRMAA
  fprintf(fp, " drmaa");
#endif
#ifdef HAVE_LIBCURL
  fprintf(fp, " curl");
#endif
#ifdef HAVE_LIBFFTW3
  fprintf(fp, " fftw3");
#endif
#ifdef HAVE_LIBCMOR
  fprintf(fp, " cmor");
#endif
#ifdef HIRLAM_EXTENSIONS
  fprintf(fp, " hirlam_extensions");
#endif
#if defined(__AVX2__)
  fprintf(fp, " avx2");
#elif defined(__AVX__)
  fprintf(fp, " avx");
#elif defined(__SSE4_2__)
  fprintf(fp, " sse4_2");
#elif defined(__SSE4_1__)
  fprintf(fp, " sse4_1");
#elif defined(__SSE3__)
  fprintf(fp, " sse3");
#elif defined(__SSE2__)
  fprintf(fp, " sse2");
#endif
  fprintf(fp, "\n");
}

void
cdo_print_libraries(void)
{
  auto fp = stdout;
  fprintf(fp, "Libraries:");
  fprintf(fp, " yac/%s", YAC_VERSION);
#ifdef HAVE_LIBNETCDF
  fprintf(fp, " NetCDF");
#ifdef NC_VERSION
  fprintf(fp, "/%s", NC_VERSION);
#endif
#endif
#ifdef HAVE_LIBHDF5
  fprintf(fp, " HDF5");
#ifdef H5_VERS_MAJOR
  unsigned h5l_majnum, h5l_minnum, h5l_relnum;
  H5get_libversion(&h5l_majnum, &h5l_minnum, &h5l_relnum);
  fprintf(fp, "/%u.%u.%u", h5l_majnum, h5l_minnum, h5l_relnum);

  unsigned h5h_majnum = H5_VERS_MAJOR, h5h_minnum = H5_VERS_MINOR, h5h_relnum = H5_VERS_RELEASE;
  if ((h5h_majnum != h5l_majnum) || (h5h_minnum != h5l_minnum) || (h5h_relnum != h5l_relnum))
    fprintf(fp, "(h%u.%u.%u)", h5h_majnum, h5h_minnum, h5h_relnum);
#endif
#endif
/*
#ifdef HAVE_LIBZ
{
  fprintf(fp, " zlib/%s", zlibVersion());
#ifdef ZLIB_VERSION
  if ( strcmp(ZLIB_VERSION, zlibVersion()) != 0 )
    fprintf(fp, "(h%s)", ZLIB_VERSION);
#else
  fprintf(fp, "(header not found)");
#endif
}
#endif
*/
#ifdef HAVE_LIBPROJ
  fprintf(fp, " proj");
#ifdef PROJ_VERSION_MAJOR
  fprintf(fp, "/%u.%u.%u", PROJ_VERSION_MAJOR, PROJ_VERSION_MINOR, PROJ_VERSION_PATCH);
#endif
#endif

#ifdef HAVE_LIBCMOR
  fprintf(fp, " cmor");
#ifdef CMOR_VERSION_MAJOR
  fprintf(fp, "/%u.%u.%u", CMOR_VERSION_MAJOR, CMOR_VERSION_MINOR, CMOR_VERSION_PATCH);
#endif
#endif

#ifdef HAVE_LIBXML2
  fprintf(fp, " xml2");
#ifdef LIBXML_DOTTED_VERSION
  fprintf(fp, "/%s", LIBXML_DOTTED_VERSION);
#endif
#endif

#ifdef HAVE_LIBCURL
  {
    auto version_data = curl_version_info(CURLVERSION_NOW);
    fprintf(fp, " curl/%s", version_data->version);
#ifdef LIBCURL_VERSION
    if (strcmp(LIBCURL_VERSION, version_data->version) != 0) fprintf(fp, "(h%s)", LIBCURL_VERSION);
#else
    fprintf(fp, "(header not found)");
#endif
  }
#endif

#ifdef HAVE_LIBMAGICS
  {
#ifdef HAVE_STDINT_H
#undef HAVE_STDINT_H
#endif
#ifdef HAVE_SYS_TYPES_H
#undef HAVE_SYS_TYPES_H
#endif
#include <magics_config.h>
#ifdef MAGICS_VERSION
    fprintf(fp, " magics/%s", MAGICS_VERSION);
#endif
  }
#endif

  fprintf(fp, "\n");
}

int
cdo_print_config(const std::string &option)
{
  int status = EXIT_SUCCESS;

  std::map<std::string, std::pair<std::string, bool>> configMap;

  // clang-format off
  configMap["has-srv"]               = {"SERVICE",           cdiHaveFiletype(CDI_FILETYPE_SRV)};
  configMap["has-ext"]               = {"EXTRA",             cdiHaveFiletype(CDI_FILETYPE_EXT)};
  configMap["has-ieg"]               = {"IEG",               cdiHaveFiletype(CDI_FILETYPE_IEG)};
  configMap["has-grb"]               = {"GRIB 1",            cdiHaveFiletype(CDI_FILETYPE_GRB)};
  configMap["has-grb1"]              = {"GRIB 1",            cdiHaveFiletype(CDI_FILETYPE_GRB)};
  configMap["has-grb2"]              = {"GRIB 2",            cdiHaveFiletype(CDI_FILETYPE_GRB2)};
  configMap["has-nc"]                = {"NetCDF",            cdiHaveFiletype(CDI_FILETYPE_NC)};
  configMap["has-nc2"]               = {"NetCDF 2",          cdiHaveFiletype(CDI_FILETYPE_NC2)};
  configMap["has-nc4"]               = {"NetCDF 4",          cdiHaveFiletype(CDI_FILETYPE_NC4)};
  configMap["has-nc4c"]              = {"NetCDF 4 classic",  cdiHaveFiletype(CDI_FILETYPE_NC4C)};
  configMap["has-nc5"]               = {"NetCDF 5",          cdiHaveFiletype(CDI_FILETYPE_NC5)};
  configMap["has-nczarr"]            = {"NetCDF 4 zarr",     cdiHaveFiletype(CDI_FILETYPE_NCZARR)};
  configMap["has-hdf5"]              = {"HDF5",              false};
  configMap["has-cgribex"]           = {"CGRIBEX",           false};
  configMap["has-cmor"]              = {"CMOR",              false};
  configMap["has-magics"]            = {"MAGICS",            false};
  configMap["has-openmp"]            = {"OPENMP",            false};
  configMap["has-proj"]              = {"PROJ",              false};
  configMap["has-threads"]           = {"PTHREADS",          false};
  configMap["has-wordexp"]           = {"WORDEXP",           false};
  configMap["has-hirlam_extensions"] = {"HIRLAM_EXTENSIONS", false};
  // clang-format on

#ifdef HAVE_LIBHDF5
  configMap["has-hdf5"].second = true;
#endif

#ifdef HAVE_LIBCGRIBEX
  configMap["has-cgribex"].second = true;
#endif

#ifdef HAVE_LIBCMOR
  configMap["has-cmor"].second = true;
#endif

#ifdef HAVE_LIBMAGICS
  configMap["has-magics"].second = true;
#endif

#ifdef _OPENMP
  configMap["has-openmp"].second = true;
#endif

#ifdef HAVE_LIBPROJ
  configMap["has-proj"].second = true;
#endif

#ifdef HAVE_LIBPTHREAD
  configMap["has-threads"].second = true;
#endif

#ifdef HAVE_WORDEXP_H
  configMap["has-wordexp"].second = true;
#endif

#ifdef HIRLAM_EXTENSIONS
  configMap["has-hirlam_extensions"].second = true;
#endif

  if ("all-json" == option || "all" == option)
    {
      std::cout << "{";
      int i = 0;
      for (const auto &entry : configMap)
        {
          if (i++) fprintf(stdout, ",");
          std::cout << "\"" << entry.first << "\":\"" << (entry.second.second ? "yes" : "no") << "\"";
        }
      std::cout << "}\n";
    }
  else
    {
      auto foundOption = false;
      for (const auto &entry : configMap)
        {
          if (entry.first == option)
            {
              foundOption = true;
              std::cout << (entry.second.second ? "yes" : "no") << "\n";
            }
        }

      if (!foundOption)
        {
          fprintf(stdout, "unknown config option: %s\n", option.c_str());
          fprintf(stdout, "\n");
          fprintf(stdout, "Available config option:\n");
          fprintf(stdout, "\n");
          for (const auto &entry : configMap)
            fprintf(stdout, "  %-12s  whether %s is enabled\n", entry.first.c_str(), entry.second.first.c_str());

          status = EXIT_FAILURE;
        }
    }

  return status;
}
