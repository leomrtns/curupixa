dnl Autoconf init, autoheader output location, config and prefix directories
AC_INIT([curupixa], [0.1], [leomrtns@gmail.com])
AC_CONFIG_AUX_DIR([.])

dnl see AX_CHECK_X86_FEATURES https://www.gnu.org/software/autoconf-archive/ax_check_x86_features.html 

AC_CONFIG_HEADER([lib/config.h])
AC_PREFIX_DEFAULT(${HOME})
dnl libtoolize recomendation
AC_CONFIG_MACRO_DIR([m4])

dnl Override default O2
CFLAGS=${CFLAGS-""}
CXXFLAGS=${CXXFLAGS-""}

dnl automake initialization (completely unrelated to compiler arguments)
AM_INIT_AUTOMAKE([-Wall -Werror])
AC_USE_SYSTEM_EXTENSIONS
AM_MAINTAINER_MODE
dnl must be called after AC_USE_SYSTEM_EXTENSIONS (new restriction in autoconf2.68) 
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])

dnl Basic compiler and related tools
AC_LANG_C
AC_PROG_CC
AC_PROG_INSTALL
dnl directs the libc header files to provide the standard GNU system interfaces including all GNU extensions
AC_GNU_SOURCE
dnl Libtool (for library): in lib/Makefile.am, we include an "abstract" libfoo.la
LT_INIT
AC_PROG_LIBTOOL

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
dnl will add lines #define SIZEOF_INT to config.h etc.
AC_CHECK_SIZEOF(int)
AC_CHECK_SIZEOF(double)
AC_CHECK_SIZEOF(long long int)

AC_MSG_RESULT([                ===    configuration options specific to curupixa])

AC_MSG_CHECKING([whether to build (slower) debug code])
AC_ARG_ENABLE([debug],
    [AS_HELP_STRING([--enable-debug], [enable debugging with gdb and friends (default=no)])],
    [debugit="$enableval"], [debugit=no])
AC_MSG_RESULT([$debugit])
if test x"$debugit" = x"yes"; then
    AC_DEFINE([CURUPIXA_DEBUG],[],[Debug Mode, with assert()-like error checking])
    AM_CFLAGS="${AM_CFLAGS} -g -pg -W -Wall -Werror -Wno-uninitialized -O1"
else
    AC_DEFINE([NDEBUG],[],[not in debug mode, disabling all assert() macros])
    AM_CFLAGS="${AM_CFLAGS} -funroll-loops -fomit-frame-pointer -finline-functions -O4"
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
