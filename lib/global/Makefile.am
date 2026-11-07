AM_CPPFLAGS = -I$(top_builddir)/lib/global
AM_CFLAGS = @AM_CFLAGS@ @OPENMP_CFLAGS@
common_headers = context_types.h global_variable.h lowlevel.h maths_and_bits.h \
  hash_functions_generators.h hash_functions.h random_number.h random_number_generators.h
common_src = global_variable.c lowlevel.c maths_and_bits.c internal_random_constants.c \
  hash_functions_generators.c hash_functions.c random_number.c random_number_generators.c

otherincludedir = $(includedir)/curupixa/global
otherinclude_HEADERS = config.h $(common_headers)
noinst_LTLIBRARIES = libcrpxglobal.la
libcrpxglobal_la_SOURCES = $(common_src) $(common_headers)
noinst_HEADERS = runtime_internal.h internal_random_constants.h
