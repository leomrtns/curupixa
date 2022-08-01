AM_CPPFLAGS = $(GTKDEPS_CFLAGS)  @OPENMP_CPPFLAGS@ -I$(srcdir)/../lib -DTEST_FILE_DIR='"$(srcdir)/files/"' 
AM_CFLAGS = @AM_CFLAGS@  @OPENMP_CFLAGS@  @CHECK_CFLAGS@ @ZLIB_CFLAGS@  @LZMA_CFLAGS@
LDADD =  ../lib/libcurupixa_static.la $(GTKDEPS_LIBS) $(AM_LDFLAGS) @CHECK_LIBS@ @ZLIB_LIBS@  @LZMA_LIBS@

EXTRA_DIST = files # directory with fasta etc files (accessed with #define TEST_FILE_DIR above)

# list of programs to be compiled only with 'make check' (like noinst_PROGRAMS)
check_PROGRAMS = check_instructions check_hashfunctions dieharder_rng dieharder_hash
# list of test programs (duplicate of above, since we want all to be compiled only with 'make check'):
TESTS = $(check_PROGRAMS)

