#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.68])
AC_INIT([fcc_tools],[dev],[fabrizio.santoro@iccom.cnr.it])
AM_INIT_AUTOMAKE
# AM_EXTRA_RECURSIVE_TARGETS([python-update-version python-standalone])
AC_LANG([Fortran])

###################################
#
# COMPILER CHECKS
#
###################################
AC_LANG([Fortran])

#--------------
# Flags f77
#--------------
# Prevent autoconf from adding flags automatically  (e.g. -g)
if test -z $FFLAGS; then
#    FFLAGS='-O2 -pg'
#    FFLAGS='-g -fbounds-check -pg'
    FFLAGS='-O3' 
#    FFLAGS='-g -fcheck=all -fbounds-check' 
fi

#--------------
# Flags f90
#--------------
if test -z $FCFLAGS; then
#    FCFLAGS='-O2 -pg'
#    FCFLAGS='-g -fbounds-check -pg'
    FCFLAGS='-O3' 
#    FCFLAGS='-g -fcheck=all -fbounds-check' # -ffpe-trap=zero,denormal' 
fi

#--------------
# Set compiler
#--------------
# FCFLAGS="$FCFLAGS -fopenmp"
# FFLAGS="$FCFLAGS -fopenmp"
AC_PROG_F77([ifort,gfortran])
 # ifort compiler works, but gfortran does not!
AC_PROG_FC([ifort,gfortran])


#########################
#
# LIBRARY CHECKS 
#
#########################

#--------------------------------
# MKL support (LAPACK and BLAS)
#--------------------------------
AC_ARG_WITH([mkl],
    AS_HELP_STRING([--with-mkl], [Use Intel MKL library for lapak subroutines]),
    ,
    [with_mkl=no])

AS_IF([test "x$with_mkl" = "xyes"],
      # Using macro:
      # AC_CHECK_LIB (library, function, [action-if-found], [action-if-not-found], [other-libraries])
      # where we use the last (optional) argument [other-libraries] to enter the additional libs
      # needed to use mkl lapack. (Previously this was done post-processing the configure: not nice)
      [AC_CHECK_LIB([mkl_intel_lp64],[zsytrf],
                    [LIBS="$LIBS -lmkl_intel_lp64 -lmkl_sequential -lmkl_core"], # default: add libraries to LIBS
                    [AC_MSG_WARN([MKL library requested but not found (LAPACK)]);use_lapack=yes],
                    [-lmkl_sequential -lmkl_core])],
      [use_lapack=yes])

AS_IF([test "x$with_mkl" = "xyes"],
      # Using macro:
      # AC_CHECK_LIB (library, function, [action-if-found], [action-if-not-found], [other-libraries])
      # where we use the last (optional) argument [other-libraries] to enter the additional libs
      # needed to use mkl lapack. (Previously this was done post-processing the configure: not nice)
      [AC_CHECK_LIB([mkl_intel_lp64],[dgemm],
                    [LIBS="$LIBS -lmkl_intel_lp64 -lmkl_sequential -lmkl_core"], # default: add libraries to LIBS
                    [AC_MSG_WARN([MKL library requested but not found (BLAS)]);use_blas=yes],
                    [-lmkl_sequential -lmkl_core])],
      [use_blas=yes])
      
      
#-------------
# LAPACK
#-------------
AS_IF([test "x$use_lapack" = "xyes"],
      [AC_CHECK_LIB([lapack],[zsytrf],
                    [have_lapack_funct=yes;LIBS="$LIBS -llapack"],
                    [have_lapack_funct=no])],
      [have_lapack_funct=yes])

AS_IF([test "x$have_lapack_funct" != "xyes"],
      [AC_MSG_ERROR([No lapack library found])],
      [have_lapack_funct=no])
      
      
#-------------
# BLAS
#-------------
AS_IF([test "x$use_blas" = "xyes"],
      [AC_CHECK_LIB([blas],[dgemm],
                    [have_blas_funct=yes;LIBS="$LIBS -lblas"],
                    [have_blas_funct=no])],
      [have_blas_funct=yes])

AS_IF([test "x$have_blas_funct" != "xyes"],
      [AC_MSG_ERROR([No blas library found])],
      [have_blas_funct=no])
          
      
#-------------
# FFTW3
#-------------
AC_ARG_WITH([fftw],
    AS_HELP_STRING([--with-fftw], [Use the Fastes Fourier Transform in the West (FFTW)]),
    ,
    [with_fftw=no;use_fftw=no])
# FFTW is mandatory since FFT is broken for the corrected DFT -> FT manipulation (to be resolved..)
# with_fftw=yes

AS_IF([test "x$with_fftw" = "xyes"],
      [AC_CHECK_LIB([fftw3],[dfftw_plan_dft_1d],
                    [use_fftw=yes;LIBS="$LIBS -lfftw3"],
                    [AC_MSG_ERROR([No fftw library found])])],
      [use_fftw=no])
#The WARNING instead of the error when FFT is fixed will be:
#                    [AC_MSG_WARN([FFTW functionalities cannot be used]);use_fftw=no])],

# Set conditional to be used by automake 
AM_CONDITIONAL([FFTW], [test x$use_fftw = xyes])

# # Checks for files in the distribution
# AC_CHECK_FILE ([$srcdir/python/dist/fcc_analyzer]
#                [install_pythonbin="yes"]
#                [install_pythonbin="no"])
# AM_CONDITIONAL([INSTALL_PYTHONBIN], [test "$install_pythonbin" = "yes"])



AC_CONFIG_FILES([Makefile
                 generators/Makefile
                 postprocessing/Makefile
                 python/Makefile
                 tests/Makefile])
AC_OUTPUT

