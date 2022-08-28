# lzma and zlib are set via pkg-config, and bz2 is set via AC_CHECK (thus @AM_@)
AM_CPPFLAGS = $(GTKDEPS_CFLAGS)  @OPENMP_CPPFLAGS@ @ZLIB_CFLAGS@ @LZMA_CFLAGS@
AM_CFLAGS = @AM_CFLAGS@  @OPENMP_CFLAGS@ @ZLIB_CFLAGS@ @LZMA_CFLAGS@

common_headers = global_variable.h lowlevel.h maths_and_bits.h internal_random_constants.h hash_functions_generators.h hash_functions.h \
								 random_number.h random_number_generators.h

common_src     = global_variable.c lowlevel.c maths_and_bits.c internal_random_constants.c hash_functions_generators.c hash_functions.c \
								 random_number.c random_number_generators.c

otherincludedir = $(includedir)/curupixa/global
otherinclude_HEADERS = config.h $(common_headers) # if headers are here (=global) should not be on SOURCES (=local but not convenience)

noinst_LTLIBRARIES = libcrpxglobal.la # convenience lib (not installed globally)
libcrpxglobal_la_SOURCES = $(common_src) $(common_headers) config.h # .c and .h listed here s.t. any change triggers recompilation 
libcrpxglobal_la_CPPFLAGS = $(AM_CPPFLAGS)
libcrpxglobal_la_CFLAGS = $(AM_CFLAGS) 
libcrpxglobal_la_LIBADD = $(GTKDEPS_LIBS) @ZLIB_LIBS@  @LZMA_LIBS@ $(AM_LDFLAGS) 
