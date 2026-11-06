SUBDIRS = global
AM_CPPFLAGS = -I$(top_builddir)/lib/global
AM_CFLAGS = @AM_CFLAGS@ @OPENMP_CFLAGS@

common_headers = index_arrangement.h quasi_random.h status.h tree.h bipartition.h char_vector.h \
  distance_matrix.h hashtable.h quickselect_quantile.h empirical_frequency.h prob_distribution.h \
  hll.h kmerhash.h rolling_hash.h suffix_tree.h topology_distance.h topology_randomise.h file_compression.h \
  alignment.h nexus_common.h gff3_format.h parsimony.h reconciliation.h minimiser_simplex.h svd.h \
  distance_generator.h clustering_goptics.h upgma.h newick_space.h lowlevel_algorithms.h splitset_distances.h
common_src = index_arrangement.c quasi_random.c status.c tree.c bipartition.c char_vector.c \
  distance_matrix.c hashtable.c quickselect_quantile.c empirical_frequency.c prob_distribution.c \
  hll.c kmerhash.c rolling_hash.c suffix_tree.c topology_distance.c topology_randomise.c file_compression.c \
  alignment.c nexus_common.c gff3_format.c parsimony.c reconciliation.c minimiser_simplex.c svd.c \
  distance_generator.c clustering_goptics.c upgma.c newick_space.c lowlevel_algorithms.c splitset_distances.c

otherincludedir = $(includedir)/curupixa
otherinclude_HEADERS = curupixa.h $(common_headers)

lib_LTLIBRARIES = libcurupixa.la
libcurupixa_la_SOURCES = $(common_src) $(common_headers) quasi_random_constants.h
libcurupixa_la_LIBADD = global/libcrpxglobal.la
libcurupixa_la_LDFLAGS = -version-info 0:0:0 @OPENMP_CFLAGS@

# Downstream Autotools projects can retain their existing convenience-library link.
noinst_LTLIBRARIES = libcurupixa_static.la
libcurupixa_static_la_SOURCES = $(libcurupixa_la_SOURCES)
libcurupixa_static_la_LIBADD = $(libcurupixa_la_LIBADD)
