# lzma and zlib are set via pkg-config, and bz2 is set via AC_CHECK (thus @AM_@)
AM_CPPFLAGS = $(GTKDEPS_CFLAGS)  @OPENMP_CPPFLAGS@ @ZLIB_CFLAGS@ @LZMA_CFLAGS@
AM_CFLAGS = @AM_CFLAGS@  @OPENMP_CFLAGS@ @ZLIB_CFLAGS@ @LZMA_CFLAGS@

SUBDIRS = global # run Makefile.am on subdirs

LOCALLIBS  = global/libcrpxglobal.la # convenience (internal) libraries

common_headers = index_arrangement.h quasi_random.h quasi_random_constants.h

common_src     = index_arrangement.c quasi_random.c 

otherincludedir = $(includedir)/curupixa
otherinclude_HEADERS = curupixa.h $(common_headers) # if headers are here (=global) should not be on SOURCES (=local)

lib_LTLIBRARIES = libcurupixa.la  # libtools 
libcurupixa_la_SOURCES = $(common_headers)
#libcurupixa_la_CPPFLAGS = $(AM_CPPFLAGS)  -I$(top_srcdir)/lib/global # not needed since libtools does it
libcurupixa_la_CPPFLAGS = $(AM_CPPFLAGS)
libcurupixa_la_CFLAGS = $(AM_CFLAGS) 
libcurupixa_la_LIBADD = $(GTKDEPS_LIBS) $(LOCALLIBS) @ZLIB_LIBS@  @LZMA_LIBS@ $(AM_LDFLAGS) 

noinst_LTLIBRARIES = libcurupixa_static.la # convenience lib (not installed globally)
libcurupixa_static_la_SOURCES = $(common_src) $(common_headers) # .c and .h listed here s.t. any change triggers recompilation 
libcurupixa_static_la_CPPFLAGS = $(libcurupixa_la_CPPFLAGS)
libcurupixa_static_la_CFLAGS = $(libcurupixa_la_CFLAGS)
libcurupixa_static_la_LIBADD = $(libcurupixa_la_LIBADD) # LOCALLIBS must be here as well

