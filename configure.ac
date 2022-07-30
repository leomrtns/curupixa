dnl Autoconf init, autoheader output location, config and prefix directories
AC_INIT([curupixa],[0.1],[leomrtns@gmail.com],[],[https://github.com/leomrtns/curupixa])
AC_CONFIG_AUX_DIR([.])

dnl see AX_CHECK_X86_FEATURES https://www.gnu.org/software/autoconf-archive/ax_check_x86_features.html 
AC_CONFIG_HEADERS([lib/config.h])
AC_PREFIX_DEFAULT(${HOME})
dnl libtoolize recomendation
AC_CONFIG_MACRO_DIR([m4])
dnl Override default O2
CFLAGS=${CFLAGS-""}

dnl automake initialization (completely unrelated to compiler arguments; subdir-objects generates .o in subdirs so avoid name clashes)
AM_INIT_AUTOMAKE([-Wall -Werror foreign color-tests subdir-objects])
dnl AM_SILENT_RULES([yes])
AC_USE_SYSTEM_EXTENSIONS
AM_MAINTAINER_MODE
dnl must be called after AC_USE_SYSTEM_EXTENSIONS (new restriction in autoconf2.68) 
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])
dnl Basic compiler and related tools (CC_C_O is a prerequisite for automake's subdir-objects)
AC_LANG([C])
AC_PROG_CC
AM_PROG_CC_C_O
AC_PROG_INSTALL
AC_PROG_LN_S
dnl directs the libc header files to provide the standard GNU system interfaces including all GNU extensions
AC_USE_SYSTEM_EXTENSIONS
dnl Libtool (for library): in lib/Makefile.am, we include an "abstract" libfoo.la
LT_INIT
dnl will add lines #define SIZEOF_INT to config.h etc.
AC_CHECK_SIZEOF(int)
AC_CHECK_SIZEOF(double)
AC_CHECK_SIZEOF(long long int)
dnl check for the clock_gettime function in rt library (linking with -lrt if found)
dnl (preferred way of checking for libraries, over AC_CHECK_LIB...)
AC_SEARCH_LIBS([clock_gettime],[rt])                                               

dnl Check for additional headers 
AC_CHECK_HEADERS( unistd.h sys/time.h sys/times.h)

AC_MSG_RESULT([                ===    configuration options specific to curupixa])

dnl openMP: sets $OPENMP_CFLAGS which should be passed to CFLAGS, CPPFLAGS; creates preprocessor macro _OPENMP 
dnl         (checked with "ifdef _OPENMP"); user can disable it through "--disable-openmp"
AC_OPENMP
AC_SUBST(OPENMP_CFLAGS)
AC_SUBST(OPENMP_CPPFLAGS)
dnl notice that default if-not-found action is to exist with error
PKG_CHECK_MODULES([CHECK], [check >= 0.9.10],
                  AC_DEFINE([HAVE_CHECK],[1],[set to one if check library is present]),
                  [AC_MSG_NOTICE([optional check library not found, not possible to run tests; please install it with "apt-get install check" or "conda check"])] )
PKG_CHECK_MODULES([ZLIB], [zlib],
                  AC_DEFINE([HAVE_ZLIB],[1],[set to one if gz compression library found]),
                  [AC_MSG_NOTICE([optional zlib headers not found; please install them with "apt-get install zlib1g-dev" or "conda zlib"])] )
PKG_CHECK_MODULES([LZMA], [liblzma],
                  AC_DEFINE([HAVE_LZMA],[1],[set to one if xz compression library found]),
                  [AC_MSG_NOTICE([optional lzma headers not found; please install them with "apt-get install liblzma-dev" or "conda xz"])] )
AC_CHECK_HEADER(bzlib.h,
                  AC_DEFINE([HAVE_BZIP2],[1],[set to one if bz2 compression library found]),
                  [AC_MSG_NOTICE([optional bzip2 headers not found; please install them with "apt-get install libbz2-dev" or "conda bzip2"])] )
AC_SEARCH_LIBS([BZ2_bzlibVersion], [bz2]) # AC_CHECK_LIB([bz2],[BZ2_bzlibVersion])
dnl check for the math library
AC_SEARCH_LIBS([sqrt], [m]) # modern equivalent of AC_CHECK_LIB([m],[main])

dnl see if SSE, AVX should be disabled (do not check against build machine since conda et al cross-compile)
dnl original configure.ac had AC_LINK_ELSE 

dnl if --enable-FEATURE is missing, it assumes --disable-FEATURE (last argument) by default; we invert it below
AC_ARG_ENABLE(sse, AS_HELP_STRING([--disable-sse], [Build without SSE support]), , enable_sse=yes)
AC_MSG_CHECKING([SSE support])
AS_IF([test "x$enable_sse" != "xno"], [
    has_sse=yes
    AC_DEFINE([HAVE_SSE], [1], [Define if SSE support is requested])
    AX_APPEND_COMPILE_FLAGS([-msse4.2], [AM_CFLAGS]) # alternative is to add to Makefile.am only in files where needed (using MK_HAVE_SSE)
  ], [
    AX_APPEND_COMPILE_FLAGS([-mno-sse4.2], [AM_CFLAGS])
  ])

dnl this follows usual pattern: disabled by default 
AC_ARG_ENABLE(avx2, AS_HELP_STRING([--enable-avx2],[enable AVX2 instructions (makes executable less portable)]),enable_avx2=yes, enable_avx2=no)
AC_MSG_CHECKING([AVX2 support])
AS_IF([test "x$enable_avx2" != "xno"], [
    has_avx2=yes
    AC_DEFINE([HAVE_AVX2], [1], [Define if AVX2 support is requested])
    AX_APPEND_COMPILE_FLAGS([-mavx2 -mfma], [AM_CFLAGS]) 
  ], [
    AX_APPEND_COMPILE_FLAGS([-mno-avx2 -mno-mfma], [AM_CFLAGS])
    ## AC_MSG_RESULT([disabled (to enable AVX2, add --enable-avx2 to configure command line)]) # append_compile_flags already announces result
  ])

AM_CONDITIONAL(MK_HAVE_SSE,  test "x${has_sse}" = "xyes") # these are automake conditionals AM, not the preprocessor symbol HAVE_SSE 
AM_CONDITIONAL(MK_HAVE_AVX2, test "x${has_avx2}" = "xyes") 

AC_ARG_ENABLE(rdrnd, AS_HELP_STRING([--disable-rdrnd], [Build without Intel RDRAND instruction support]), , enable_rdrnd=yes) 
AC_MSG_CHECKING([RDRAND support])
AS_IF([test "x$enable_rdrnd" != "xno"], [
    has_rdrnd=yes
    AC_DEFINE([HAVE_RDRND], [1], [Define if Intel RDRAND instruction support is requested])
    AX_APPEND_COMPILE_FLAGS([-mrdrnd], [AM_CFLAGS])
  ], []) # no need to create MK_HAVE_RDRND since it is used only by automake

AC_MSG_CHECKING([whether to build debug code (development version only, not for final user)])
AC_ARG_ENABLE([debug],
    [AS_HELP_STRING([--enable-debug], [enable debugging with gdb and friends (default=no)])],
    [debugit="$enableval"], [debugit=no])
AC_MSG_RESULT([$debugit])
if test x"$debugit" = x"yes"; then
    AC_DEFINE([CURUPIXA_DEBUG],[],[Debug Mode, with assert()-like error checking])
    AX_APPEND_COMPILE_FLAGS([-g -pg -W -Wall -Werror -Wno-uninitialized -O0 -fsanitize=undefined -fsanitize=address -fno-omit-frame-pointer], [AM_CFLAGS])
else
    AC_DEFINE([NDEBUG],[],[not in debug mode, disabling all assert() macros])
    AX_APPEND_COMPILE_FLAGS([-funroll-loops -fomit-frame-pointer -finline-functions -O3], [AM_CFLAGS])
fi

AC_MSG_CHECKING([whether you want static binaries (not the library, but the executables)])
AC_ARG_ENABLE(static-binary,
    [  --enable-static-binary  static binaries, that run on same arch without the libraries [[default=no]]],
    [ statbin_use="yes" ], [ statbin_use=""    ])
if test -n "${statbin_use}"; then
    AC_MSG_RESULT([yes])
    AM_LDFLAGS="-static ${AM_LDFLAGS}";
else
    AC_MSG_RESULT([no])
fi

AC_MSG_RESULT([                ===    end of specific configuration options])

dnl propagate changed vars among final makefiles
AC_SUBST([AM_CFLAGS])
AC_SUBST([AM_LDFLAGS])
AC_SUBST([MPI_CXXLIBS])
AC_SUBST([MPI_CXXFLAGS])

dnl generate makefiles (last recipes, after defining CFLAGS etc.)
AC_CONFIG_FILES([Makefile lib/Makefile tests/Makefile])
AC_OUTPUT

echo \
"---- 
Configuration parameters for curupixa-lib [$PACKAGE $VERSION]:
 Target:          $host_os $host_cpu
 Source location: ${srcdir}
 Compiler:        ${CC}
 Compiler flags:  ${AM_CFLAGS} 
 Linker flags:    ${AM_LDFLAGS}
 LIBS:            ${LIBS} ${LDFLAGS}
 Install path:    ${prefix}

----"
