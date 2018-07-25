#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([f5c], [1.0], [hasindu@unsw.edu.au])
AM_INIT_AUTOMAKE([foreign])
AC_CONFIG_SRCDIR([f5c.c])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CC([g++])
AC_PROG_CXX([g++])

AC_SYS_LARGEFILE
AC_FUNC_FSEEKO

# Checks for header files.
AC_CHECK_HEADERS([float.h inttypes.h stdint.h stdlib.h string.h sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_INT32_T
AC_TYPE_SIZE_T
AC_TYPE_SSIZE_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_TYPE_UINT8_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([gettimeofday memset sqrt strerror])



# Checks for libraries.
libhts=yes
AC_CHECK_HEADER([htslib/hts.h], [], [libhts=no], [;])
AC_CHECK_LIB([hts], [faidx_fetch_seq], [], [libhts=no])

if test $libhts != yes; then
  AC_MSG_ERROR([HTS library development files not found! 
  On Debian based systems such as Ubuntu you may use sudo apt-get install libhts-dev.
  Or you download and install the latest release of htslib from http://www.htslib.org/download/."
  ])  
fi



AC_CHECK_HEADER([hdf5.h], [libhdf=yes], [], [;])

AC_CHECK_LIB([hdf5_serial], [H5Fopen], [], [libhdf=no])

if test $libhdf != yes; then
  AC_MSG_ERROR([HDF5 library development files not found! 
  On Debian based systems such as Ubuntu you may use sudo apt-get install libhdf5-dev.
  ])  
fi



AC_CONFIG_FILES([Makefile])
AC_OUTPUT([])
