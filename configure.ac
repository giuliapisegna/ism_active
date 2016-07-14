#
# configure.ac -- autoconf script for active simulations
#
# copyright (C) 2015, 2016 by Tomas S. Grigera
#
# Copying / redistribution not permitted
#

###############################################################################
#
# Initialization
#
AC_PREREQ([2.67])
AC_INIT([active],[0.0.1],[tgrigera@iflysib.unlp.edu.ar])
AC_CONFIG_SRCDIR([ISM/social.hh])  # Check just one file as safety check
AC_CONFIG_MACRO_DIR([m4])
AC_LANG(C++)
AM_INIT_AUTOMAKE
AC_OPENMP

###############################################################################
#
# Checks for programs needed to build
#
AC_PROG_MAKE_SET
AC_PROG_CPP
AC_PROG_CC
AC_PROG_CXX
AX_CXX_COMPILE_STDCXX(11,noext,mandatory)
AC_PROG_LN_S
AC_PROG_INSTALL

###############################################################################
#
# Checks for header files.
#
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h unistd.h])
AC_CHECK_HEADERS([glsim/blib.hh],,exit)


###############################################################################
#
# Checks for libraries
#

# GSL
AC_CHECK_LIB([gslcblas], [cblas_dsdot],,exit)
AC_CHECK_LIB([gsl], [gsl_rng_alloc],,exit)

# Boost::system
AC_MSG_CHECKING([for working -lboost_system])
SAVED_LIBS=$LIBS
LIBS="-lboost_system $LIBS"
AC_LINK_IFELSE(
  [AC_LANG_PROGRAM([#include <boost/system/error_code.hpp>],
  	           [boost::system::error_code e;])],
  AC_MSG_RESULT(yes) ; AC_DEFINE([HAVE_LIBBOOST_SYSTEM],[1],[Define if Boost::system is available]),
  AC_MSG_RESULT(no) ; [LIBS=$SAVED_LIBS] ; [AC_MSG_FAILURE([cannot build without boost::system])])

# Boost::program_options
AC_MSG_CHECKING([for working -lboost_program_options])
SAVED_LIBS=$LIBS
LIBS="-lboost_program_options $LIBS"
AC_LINK_IFELSE(
  [AC_LANG_PROGRAM([#include <boost/program_options.hpp>],
    [boost::program_options::options_description opt])],
  AC_MSG_RESULT(yes) ; AC_DEFINE([HAVE_LIBBOOST_PROGRAM_OPTIONS],[1],[Define if Boost::program_options available]) ,
  AC_MSG_RESULT(no) ; [LIBS=$SAVED_LIBS] ; [AC_MSG_FAILURE([cannot build without boost::program_options])])

# Boost::serialization
AC_MSG_CHECKING([for working -lboost_serialization])
SAVED_LIBS=$LIBS
LIBS="-lboost_serialization $LIBS"
AC_LINK_IFELSE(
  [AC_LANG_PROGRAM([#include <boost/archive/binary_oarchive.hpp>
                    #include <fstream>],
    [ std::ofstream os("tst.dat",std::ios::binary);
      boost::archive::binary_oarchive oa(os);])],
  AC_MSG_RESULT(yes) ; AC_DEFINE([HAVE_LIBBOOST_SERIALIZATION],[1],[Define if Boost::serialization is available]),
  AC_MSG_RESULT(no) ; [LIBS=$SAVED_LIBS] ; [AC_MSG_FAILURE([cannot build without boost::serialization])])

# Boost::date_time
AC_MSG_CHECKING([for working -lboost_date_time])
SAVED_LIBS=$LIBS
LIBS="-lboost_date_time $LIBS"
AC_LINK_IFELSE(
  [AC_LANG_PROGRAM([#include <boost/date_time/posix_time/posix_time.hpp>],
  	           [boost::posix_time::ptime now = boost::posix_time::microsec_clock::local_time();])],
  AC_MSG_RESULT(yes) ; AC_DEFINE([HAVE_LIBBOOST_DATE_TIME],[1],[Define if Boost::date_time is available]),
  AC_MSG_RESULT(no) ; [LIBS=$SAVED_LIBS] ; [AC_MSG_FAILURE([cannot build without boost::date_time])])

# Boost::timer
AC_MSG_CHECKING([for working -lboost_timer])
SAVED_LIBS=$LIBS
LIBS="-lboost_timer $LIBS"
AC_LINK_IFELSE(
  [AC_LANG_PROGRAM([#include <boost/timer/timer.hpp>],
  	           [boost::timer::auto_cpu_timer t;])],
  AC_MSG_RESULT(yes) ; AC_DEFINE([HAVE_LIBBOOST_TIMER],[1],[Define if Boost::timer is available]),
  AC_MSG_RESULT(no) ; [LIBS=$SAVED_LIBS] ; [AC_MSG_FAILURE([cannot build without boost::timer])])

# HDF5
AC_CHECK_LIB([hdf5], H5Fcreate,,exit)
AC_CHECK_LIB([hdf5_hl], H5LTset_attribute_short,,exit)

# FFTW3
AC_ARG_WITH([fftw3],AS_HELP_STRING([--without-fftw3],[do not use the FFTW3 library even if installed]),,with_fftw3=yes)
AS_IF([test x"$with_fftw3" = xyes],
   [ AC_CHECK_LIB([fftw3], fftw_execute,,
     [echo "**** WARNING"
      echo "**** FFTW3 library not found"
      echo "**** FFT calls will will use the GSL FFT routines, but"
      echo "**** will fail unless the number of data points is a power of 2"]) ]
)

#
# glsim, glsim_ol, glsim_an
#
AC_MSG_CHECKING([for working -lglsim])
SAVED_LIBS=$LIBS
LIBS="-lglsim $LIBS"
AC_LINK_IFELSE(
  [AC_LANG_PROGRAM([#include <glsim/environment.hh>],
  	           [glsim::Environment e;])],
  AC_MSG_RESULT(yes) ; AC_DEFINE([HAVE_GLSIM],[1],[Define if libglsim is available]),
  AC_MSG_RESULT(no) ; [LIBS=$SAVED_LIBS] ; [AC_MSG_FAILURE([cannot build without libglsim])])
AC_MSG_CHECKING([for working -lglsim_ol])
SAVED_LIBS=$LIBS
LIBS="-lglsim_ol $LIBS"
AC_LINK_IFELSE(
  [AC_LANG_PROGRAM([#include <glsim/olconfiguration.hh>],
  	           [glsim::OLconfiguration conf;])],
  AC_MSG_RESULT(yes) ; AC_DEFINE([HAVE_GLSIM_OL],[1],[Define if libglsim_ol is available]),
  AC_MSG_RESULT(no) ; [LIBS=$SAVED_LIBS] ; [AC_MSG_FAILURE([cannot build without libglsim_ol])])


###############################################################################
#
# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T
AC_STRUCT_TM
AC_C_VOLATILE

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_MALLOC
AC_FUNC_STAT
AC_FUNC_STRFTIME

###############################################################################
#
# Generates files
#
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile
		 ISM/Makefile])
AC_OUTPUT
