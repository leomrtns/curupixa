# while lzma and zlib are set via pkg-config, bz2 is set via AC_CHECK (thus @AM_@)
AM_CPPFLAGS = $(GTKDEPS_CFLAGS)  @OPENMP_CPPFLAGS@ @ZLIB_CFLAGS@ @LZMA_CFLAGS@
AM_CFLAGS = @AM_CFLAGS@  @OPENMP_CFLAGS@ @ZLIB_CFLAGS@ @LZMA_CFLAGS@

# ideally I would hide reconciliation.h etc used only inside genetree.c 
common_headers = hashtable.h lowlevel.h random_number_gen.h constant_random_lists.h random_number.h nexus_common.h \
                 bipartition.h prob_distribution.h empirical_frequency.h biomcmc.h argtable3.h \
                 distance_matrix.h alignment.h topology_common.h parsimony.h genetree.h \
                 reconciliation.h splitset_distances.h read_newick_trees.h char_vector.h \
                 upgma.h topology_randomise.h newick_space.h topology_space.h topology_distance.h \
                 kmerhash.h hashfunctions.h distance_generator.h clustering_goptics.h \
                 quickselect_quantile.h fortune_cookies.h suffix_tree.h phylogeny.h likelihood.h \
								 gff3_format.h file_compression.h 
                 
common_src     = hashtable.c lowlevel.c random_number_gen.c constant_random_lists.c random_number.c nexus_common.c \
                 bipartition.c prob_distribution.c empirical_frequency.c argtable3.c \
                 distance_matrix.c alignment.c topology_common.c parsimony.c genetree.c \
                 reconciliation.c splitset_distances.c read_newick_trees.c char_vector.c \
                 upgma.c topology_randomise.c newick_space.c topology_space.c topology_distance.c \
                 kmerhash.c hashfunctions.c distance_generator.c clustering_goptics.c \
                 quickselect_quantile.c fortune_cookies.c suffix_tree.c phylogeny.c likelihood.c \
								 gff3_format.c file_compression.c

otherincludedir = $(includedir)/biomcmc
otherinclude_HEADERS = config.h $(common_headers) # if headers are here (=global) should not be on SOURCES (=local)

lib_LTLIBRARIES = libbiomcmc.la  # libedlib.la # libtools 
libbiomcmc_la_SOURCES = $(common_headers) config.h
libbiomcmc_la_CPPFLAGS = $(AM_CPPFLAGS)
libbiomcmc_la_CFLAGS = $(AM_CFLAGS) 
libbiomcmc_la_LIBADD = $(GTKDEPS_LIBS) libbiomcmc_static.la @ZLIB_LIBS@  @LZMA_LIBS@ $(AM_LDFLAGS) 

noinst_LTLIBRARIES = libbiomcmc_static.la # convenience lib (not installed globally)
libbiomcmc_static_la_SOURCES = $(common_src) $(common_headers) config.h # listed here s.t. any change triggers recompilation 
libbiomcmc_static_la_CPPFLAGS = $(libbiomcmc_la_CPPFLAGS)
libbiomcmc_static_la_CFLAGS = $(libbiomcmc_la_CFLAGS)

#libbiomcmc_la_LIBADD   = libedlib.la
#libedlib_la_SOURCES = edlib.cpp edlib.h
#libedlib_la_CXXFLAGS =  -std=c++03 $(CXXFLAGS) #c++98 c++03 c++11
## TODO: mast hll rolling_hash svd 
