dnl Autoconf init, autoheader output location, config and prefix directories
AC_INIT([curupixa], [0.1], [leomrtns@gmail.com],,[https://github.com/leomrtns/curupixa])
AC_CONFIG_AUX_DIR([.])

dnl see AX_CHECK_X86_FEATURES https://www.gnu.org/software/autoconf-archive/ax_check_x86_features.html 
AC_CONFIG_HEADER([lib/config.h])
AC_PREFIX_DEFAULT(${HOME})
dnl libtoolize recomendation
AC_CONFIG_MACRO_DIR([m4])
dnl Override default O2
CFLAGS=${CFLAGS-""}

dnl automake initialization (completely unrelated to compiler arguments; subdir-objects generates .o in subdirs so avoid name clashes)
AM_INIT_AUTOMAKE([-Wall -Werror foreign color-tests subdir-objects])
AM_SILENT_RULES([yes])
AC_USE_SYSTEM_EXTENSIONS
AM_MAINTAINER_MODE
dnl must be called after AC_USE_SYSTEM_EXTENSIONS (new restriction in autoconf2.68) 
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])
dnl Basic compiler and related tools (CC_C_O is a prerequisite for automake's subdir-objects)
AC_LANG_C
AC_PROG_CC
AM_PROG_CC_C_O
AC_PROG_INSTALL
AC_PROG_LN_S
dnl directs the libc header files to provide the standard GNU system interfaces including all GNU extensions
AC_GNU_SOURCE
dnl Libtool (for library): in lib/Makefile.am, we include an "abstract" libfoo.la
LT_INIT
AC_PROG_LIBTOOL
dnl will add lines #define SIZEOF_INT to config.h etc.
AC_CHECK_SIZEOF(int)
AC_CHECK_SIZEOF(double)
AC_CHECK_SIZEOF(long long int)
dnl check for the clock_gettime function in rt library (linking with -lrt if found)
dnl (preferred way of checking for libraries, over AC_CHECK_LIB...)
AC_SEARCH_LIBS([clock_gettime],[rt])                                               
dnl Check for the 15 Standard C (ANSI-C/ISO-C) headers
AC_HEADER_STDC
dnl Check for additional headers 
AC_CHECK_HEADERS( \
    unistd.h \
    sys/time.h \
    sys/times.h \
)

AC_MSG_RESULT([                ===    configuration options specific to curupixa])

dnl openMP: sets $OPENMP_CFLAGS which should be passed to CFLAGS, CPPFLAGS; creates preprocessor macro _OPENMP 
dnl         (checked with "ifdef _OPENMP"); user can disable it through "--disable-openmp"
AC_OPENMP
AC_SUBST(OPENMP_CFLAGS)
AC_SUBST(OPENMP_CPPFLAGS)
dnl notice that default if-not-found action is to exist with error
PKG_CHECK_MODULES([CHECK], [check >= 0.9.10]) # unit tests
PKG_CHECK_MODULES([ZLIB], [zlib],
                  AC_DEFINE([HAVE_ZLIB],[1],[set to one if gz compression library found]),
                  [AC_MSG_NOTICE([optional zlib headers not found])]  )
PKG_CHECK_MODULES([LZMA], [liblzma],
                  AC_DEFINE([HAVE_LZMA],[1],[set to one if xz compression library found]),
                  [AC_MSG_NOTICE([optional lzma headers not found])]  )
AC_CHECK_HEADER(bzlib.h,
                  AC_DEFINE([HAVE_BZIP2],[1],[set to one if bz2 compression library found]),
                  [AC_MSG_NOTICE([optional bzip2 headers not found])] )
AC_SEARCH_LIBS([BZ2_bzlibVersion], [bz2]) # AC_CHECK_LIB([bz2],[BZ2_bzlibVersion])
dnl check for the math library
AC_SEARCH_LIBS([sqrt], [m]) # modern equivalent of AC_CHECK_LIB([m],[main])

dnl if --enable-FEATURE is missing, it assumes --disable-FEATURE (last argument) by default; we invert it below
AC_ARG_ENABLE(sse2,  AC_HELP_STRING([--disable-sse2], [disable SSE2 instructions]), , enable_sse2=yes)
if test "$enable_sse2" != no; then
  dnl decide if system actually allows SSE2 instructions (MSG_CHECKING and MSG_RESULT go together to print to user)
  AC_MSG_CHECKING(for SSE2 support)
  AC_LINK_IFELSE([ AC_LANG_PROGRAM([[ 
    #include <immintrin.h>
    __m128 testtunc(float *a, float *b) {return _mm_add_ps(_mm_loadu_ps(a), _mm_loadu_ps(b));}
  ]])], [has_sse2=yes], [has_sse2=no])
  AC_MSG_RESULT($has_sse2) 
  AS_IF([test "x$has_sse2" != xno], [
       dnl gcc already defines __SSE2__ but we create an extra HAVE_SSE2 macro
       AC_DEFINE([HAVE_SSE2], [1], [Define if SSE2 support is available at build])
       AX_APPEND_COMPILE_FLAGS([-msse2], [AM_CFLAGS])
       ], [
       AC_MSG_NOTICE([SSE2 not supported by this system])
       AX_APPEND_COMPILE_FLAGS([-mno-sse2], [AM_CFLAGS])
       ])
else
  AC_MSG_NOTICE([SSE2 support disabled at request of user])
  AX_APPEND_COMPILE_FLAGS([-mno-sse2], [AM_CFLAGS])
fi

AC_ARG_ENABLE(avx2, AC_HELP_STRING([--enable-avx2], [enable AVX2 instructions (makes executable less portable)]),enable_avx2=yes, enable_avx2=no)
if test "$enable_avx2" != no; then
  AC_MSG_CHECKING(for AVX2 support)
  AC_LINK_IFELSE([AC_LANG_PROGRAM([[
    #pragma GCC target("avx2")
    #include <immintrin.h>
    ]],[[
    __m256 x = _mm256_set1_ps(7); __m256 y = _mm256_permutevar8x32_ps(x, _mm256_set1_epi32(21));
    return _mm256_movemask_ps(_mm256_cmp_ps(x, y, _CMP_NEQ_OQ));
  ]])], [has_avx2=yes], [has_avx2=no])
  AC_MSG_RESULT($enable_avx2)
  AS_IF([test "x$has_avx2" != xno], [
       AC_DEFINE([HAVE_AVX2], [1], [Define if AVX2 support is available at build])
       AX_APPEND_COMPILE_FLAGS([-mavx2], [AM_CFLAGS])
       ], [
       AC_MSG_NOTICE([AVX2 not supported by this system])
       AX_APPEND_COMPILE_FLAGS([-mno-avx2], [AM_CFLAGS])
       ])
fi

AC_MSG_CHECKING([whether to build (slower) debug code])
AC_ARG_ENABLE([debug],
    [AS_HELP_STRING([--enable-debug], [enable debugging with gdb and friends (default=no)])],
    [debugit="$enableval"], [debugit=no])
AC_MSG_RESULT([$debugit])
if test x"$debugit" = x"yes"; then
    AC_DEFINE([CURUPIXA_DEBUG],[],[Debug Mode, with assert()-like error checking])
    AX_APPEND_COMPILE_FLAGS([-g -pg -W -Wall -Werror -Wno-uninitialized -O0], [AM_CFLAGS])
else
    AC_DEFINE([NDEBUG],[],[not in debug mode, disabling all assert() macros])
    AX_APPEND_COMPILE_FLAGS([-funroll-loops -fomit-frame-pointer -finline-functions -O4], [AM_CFLAGS])
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
Configuration parameters for curupixa-lib:

 Source code location:   ${srcdir}
 Compiler     (library): ${CC}
 Compiler flags: ${AM_CFLAGS} 
 Linker  flags:  ${AM_LDFLAGS}
 Install path:   ${prefix}

----"
