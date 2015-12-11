#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.68])
AC_INIT([fcc_tools],[dev],[fabrizio.santoro@iccom.cnr.it])
AM_INIT_AUTOMAKE
AM_EXTRA_RECURSIVE_TARGETS([python-update-version python-standalone])
AC_LANG([Fortran])

# Checks for programs.
# Prevent autoconf from adding flags automatically  (e.g. -g)
if test -z $FFLAGS; then
    FFLAGS=''
fi
if test -z $FCFLAGS; then
    FCFLAGS=''
fi
# AC_PROG_F77([gfortran,ifort])
# AC_PROG_FC([gfortran,ifort])


# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
# Checks for libraries: LAPACK
AC_ARG_WITH([mkl],
    AS_HELP_STRING([--with-mkl], [Use Intel MKL library for lapak subroutines]),
    [with_mkl=yes],
    [with_mkl=no])

AS_IF([test "x$with_mkl" = "xyes"],
      [AC_CHECK_LIB([mkl_intel],[dsyev],
                    [use_lapack=no;LIBS="-lmkl_intel"],
                    [AC_MSG_WARN([MKL library requested but not found]);use_lapack=yes])],
      # MKL test needs some changes in the configure file:
      # sed -i "s/-lmkl_intel/-lmkl_intel_lp64 -lmkl_sequential -lmkl_core/g" configure
      [use_lapack=yes])

AS_IF([test "x$use_lapack" = "xyes"],
      [AC_CHECK_LIB([lapack],[dsyev],
                    [have_lapack_funct=yes;LIBS=-llapack],
                    [have_lapack_funct=no])],
      [have_lapack_funct=yes])

AS_IF([test "x$have_lapack_funct" != "xyes"],
      [AC_MSG_ERROR([No lapack library found])],
      [have_lapack_funct=yes])

# # Checks for files in the distribution
# AC_CHECK_FILE ([$srcdir/python/dist/fcc_analyzer]
#                [install_pythonbin="yes"]
#                [install_pythonbin="no"])
# AM_CONDITIONAL([INSTALL_PYTHONBIN], [test "$install_pythonbin" = "yes"])



AC_CONFIG_FILES([Makefile
                 src/generators/Makefile
                 python/Makefile
                 tests/Makefile])
AC_OUTPUT

