AC_INIT([oasis3-mct], [340], [rene.redler@zmaw.de])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_SRCDIR([libpsmile/oasis_os.h])
AC_CONFIG_HEADERS([config/config.h])

AM_INIT_AUTOMAKE

# Checks for programs.
AC_PROG_CC
AC_PROG_FC
AC_PROG_F77
AC_PROG_MKDIR_P
AC_PROG_INSTALL

# Checks for libraries.

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
AX_LIB_ZLIB
AX_LIB_SZIP
AX_LIB_HDF5

# AX_LIB_NETCDF4 # This script is broken


## Edit: manually supply netcdf libs etc {{{
NETCDF_INCLUDE=${NETCDF_INCLUDE}
NETCDFROOT=${NETCDFROOT}
NETCDFFROOT=${NETCDFFROOT}
NETCDF_FINCLUDE=${NETCDF_FINCLUDE}
NETCDF_FLIBS=${NETCDF_FLIBS}
NETCDF_LIBS=${NETCDF_LIBS}
NETCDF_VERSION=${NETCDF_VERSION}


AC_SUBST([NETCDF_INCLUDE])
AC_SUBST([NETCDFROOT])
AC_SUBST([NETCDFFROOT])
AC_SUBST([NETCDF_FINCLUDE])
AC_SUBST([NETCDF_FLIBS])
AC_SUBST([NETCDF_LIBS])
AC_SUBST([NETCDF_VERSION])
## Edit }}}

# Create output files
MPIFC=${MPIFC:-$FC}
AC_SUBST([MPIFC])

AC_CONFIG_FILES([Makefile libscrip/Makefile libpsmile/Makefile])

AC_CONFIG_SUBDIRS([mct])

AC_OUTPUT
