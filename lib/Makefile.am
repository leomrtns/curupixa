# lzma and zlib are set via pkg-config, and bz2 is set via AC_CHECK (thus @AM_@)
AM_CPPFLAGS = $(GTKDEPS_CFLAGS)  @OPENMP_CPPFLAGS@ @ZLIB_CFLAGS@ @LZMA_CFLAGS@
AM_CFLAGS = @AM_CFLAGS@  @OPENMP_CFLAGS@ @ZLIB_CFLAGS@ @LZMA_CFLAGS@

common_headers = curupixa.h lowlevel.h maths_and_bits.h random_constants.h hash_functions.h random_number_generators.h

common_src     = curupixa.c lowlevel.c maths_and_bits.c random_constants.c hash_functions.c random_number_generators.c

otherincludedir = $(includedir)/curupixa
otherinclude_HEADERS = config.h $(common_headers) # if headers are here (=global) should not be on SOURCES (=local)

lib_LTLIBRARIES = libcurupixa.la  # libtools 
libcurupixa_la_SOURCES = $(common_headers) config.h
libcurupixa_la_CPPFLAGS = $(AM_CPPFLAGS)
libcurupixa_la_CFLAGS = $(AM_CFLAGS) 
libcurupixa_la_LIBADD = $(GTKDEPS_LIBS) libcurupixa_static.la @ZLIB_LIBS@  @LZMA_LIBS@ $(AM_LDFLAGS) 

noinst_LTLIBRARIES = libcurupixa_static.la # convenience lib (not installed globally)
libcurupixa_static_la_SOURCES = $(common_src) $(common_headers) config.h # listed here s.t. any change triggers recompilation 
libcurupixa_static_la_CPPFLAGS = $(libcurupixa_la_CPPFLAGS)
libcurupixa_static_la_CFLAGS = $(libcurupixa_la_CFLAGS)

#libcurupixa_la_LIBADD  = libedlib.la
#libedlib_la_SOURCES = edlib.cpp edlib.h
#libedlib_la_CXXFLAGS =  -std=c++03 $(CXXFLAGS) #c++98 c++03 c++11
