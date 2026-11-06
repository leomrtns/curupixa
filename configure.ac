dnl Retained for downstream Autotools projects; CMake is the recommended standalone build.
AC_INIT([curupixa],[0.1.0],[leomrtns@gmail.com],[],[https://github.com/leomrtns/curupixa])
AC_CONFIG_AUX_DIR([.])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS([lib/global/config.h])
AM_INIT_AUTOMAKE([-Wall foreign subdir-objects])
AC_USE_SYSTEM_EXTENSIONS
AM_MAINTAINER_MODE
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])
AC_LANG([C])
AC_PROG_CC
AM_PROG_CC_C_O
AC_PROG_INSTALL
LT_INIT

AC_CHECK_HEADERS([unistd.h sys/time.h sys/times.h sys/syscall.h threads.h])
AC_SEARCH_LIBS([clock_gettime], [rt])
AC_SEARCH_LIBS([sqrt], [m])
AC_SEARCH_LIBS([pthread_create], [pthread])

dnl A parent package controls optimisation, architecture and instrumentation flags.
AX_APPEND_COMPILE_FLAGS([-std=c17], [AM_CFLAGS])
AC_ARG_ENABLE([openmp], AS_HELP_STRING([--enable-openmp], [enable optional OpenMP compilation]), [],
  [enable_openmp=no])
AS_IF([test "x$enable_openmp" = xyes], [AC_OPENMP])
AC_SUBST([OPENMP_CFLAGS])

AC_ARG_ENABLE([sse], AS_HELP_STRING([--enable-sse], [require SSE4.2 on the target CPU]), [], [enable_sse=no])
AS_IF([test "x$enable_sse" = xyes], [
  AC_DEFINE([HAVE_SSE], [1], [Define if SSE4.2 instructions are requested])
  AX_APPEND_COMPILE_FLAGS([-msse4.2], [AM_CFLAGS])
])
AC_ARG_ENABLE([avx2], AS_HELP_STRING([--enable-avx2], [require AVX2 and FMA on the target CPU]), [],
  [enable_avx2=no])
AS_IF([test "x$enable_avx2" = xyes], [
  AC_DEFINE([HAVE_AVX], [1], [Define if AVX2 and FMA instructions are requested])
  AX_APPEND_COMPILE_FLAGS([-mavx2 -mfma], [AM_CFLAGS])
])
AC_ARG_ENABLE([rdrnd], AS_HELP_STRING([--enable-rdrnd], [require RDRAND on the target CPU]), [],
  [enable_rdrnd=no])
AS_IF([test "x$enable_rdrnd" = xyes], [
  AC_DEFINE([HAVE_RDRND], [1], [Define if RDRAND instructions are requested])
  AX_APPEND_COMPILE_FLAGS([-mrdrnd], [AM_CFLAGS])
])
AC_ARG_ENABLE([debug], AS_HELP_STRING([--enable-debug], [enable compiler diagnostics]), [], [enable_debug=no])
AS_IF([test "x$enable_debug" = xyes], [
  AC_DEFINE([CURUPIXA_DEBUG], [1], [Enable debugging diagnostics])
  AX_APPEND_COMPILE_FLAGS([-g -Wall -Wextra], [AM_CFLAGS])
])

AC_ARG_ENABLE([zlib], AS_HELP_STRING([--enable-zlib], [enable gzip file streams]), [], [enable_zlib=no])
AS_IF([test "x$enable_zlib" = xyes], [
  AC_CHECK_HEADER([zlib.h], [], [AC_MSG_ERROR([zlib headers are required for --enable-zlib])])
  AC_CHECK_LIB([z], [gzopen], [], [AC_MSG_ERROR([zlib is required for --enable-zlib])])
  AC_DEFINE([HAVE_ZLIB], [1], [Enable gzip file streams])
])
AC_ARG_ENABLE([bzip2], AS_HELP_STRING([--enable-bzip2], [enable bzip2 file streams]), [], [enable_bzip2=no])
AS_IF([test "x$enable_bzip2" = xyes], [
  AC_CHECK_HEADER([bzlib.h], [], [AC_MSG_ERROR([bzip2 headers are required for --enable-bzip2])])
  AC_CHECK_LIB([bz2], [BZ2_bzopen], [], [AC_MSG_ERROR([bzip2 is required for --enable-bzip2])])
  AC_DEFINE([HAVE_BZIP2], [1], [Enable bzip2 file streams])
])
AC_ARG_ENABLE([lzma], AS_HELP_STRING([--enable-lzma], [enable XZ file streams]), [], [enable_lzma=no])
AS_IF([test "x$enable_lzma" = xyes], [
  AC_CHECK_HEADER([lzma.h], [], [AC_MSG_ERROR([liblzma headers are required for --enable-lzma])])
  AC_CHECK_LIB([lzma], [lzma_stream_decoder], [], [AC_MSG_ERROR([liblzma is required for --enable-lzma])])
  AC_DEFINE([HAVE_LZMA], [1], [Enable XZ file streams])
])

CURUPIXA_PC_LIBS_PRIVATE="$LIBS $OPENMP_CFLAGS"
CURUPIXA_PC_CFLAGS="$OPENMP_CFLAGS"
AC_SUBST([CURUPIXA_PC_LIBS_PRIVATE])
AC_SUBST([CURUPIXA_PC_CFLAGS])
AC_SUBST([AM_CFLAGS])
AC_CONFIG_FILES([Makefile lib/Makefile lib/global/Makefile tests/Makefile curupixa.pc])
AC_OUTPUT
