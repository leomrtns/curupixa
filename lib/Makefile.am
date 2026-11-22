SUBDIRS = global
AM_CPPFLAGS = -I$(top_builddir)/lib/global
AM_CFLAGS = @AM_CFLAGS@ @OPENMP_CFLAGS@

phylogenetics_headers = index_arrangement.h quasi_random.h status.h tree.h bipartition.h \
  char_vector.h distance_matrix.h hashtable.h quickselect_quantile.h empirical_frequency.h \
  prob_distribution.h hll.h kmerhash.h rolling_hash.h suffix_tree.h \
  topology_distance.h topology_randomise.h file_compression.h alignment.h nexus_common.h \
  gff3_format.h parsimony.h minimiser_simplex.h svd.h distance_generator.h \
  clustering_goptics.h upgma.h newick_space.h lowlevel_algorithms.h splitset_distances.h \
  spr_distance.h ancestry_index.h character_data.h character_analysis.h substitution_model.h \
  likelihood.h continuous_traits.h reconstruction.h quartet_distance.h tanglegram.h \
  dtl.h phylogenetics.h sampling.h quasi_sequence.h tree_comparison.h triplet_distance.h \
  gaussian.h phylogenetic_regression.h trait_association.h stochastic_mapping.h tip_correspondence.h
phylogenetics_sources = index_arrangement.c quasi_random.c status.c tree.c bipartition.c \
  char_vector.c distance_matrix.c hashtable.c quickselect_quantile.c empirical_frequency.c \
  prob_distribution.c hll.c kmerhash.c rolling_hash.c suffix_tree.c \
  topology_distance.c topology_randomise.c file_compression.c alignment.c nexus_common.c \
  gff3_format.c parsimony.c minimiser_simplex.c svd.c distance_generator.c \
  clustering_goptics.c upgma.c newick_space.c lowlevel_algorithms.c splitset_distances.c \
  spr_distance.c ancestry_index.c character_data.c character_analysis.c substitution_model.c \
  likelihood.c likelihood_kernel_internal.c continuous_traits.c reconstruction.c quartet_distance.c \
  tanglegram.c dtl.c sampling.c quasi_sequence.c tree_comparison.c triplet_distance.c \
  gaussian.c phylogenetic_regression.c trait_association.c stochastic_mapping.c tip_correspondence.c
gene_species_headers = reconciliation.h reconciliation_roots.h taxon_map.h gene_species_distances.h paralog_quartets.h gene_species.h
gene_species_sources = reconciliation.c reconciliation_roots.c taxon_map.c gene_species_distances.c paralog_quartets.c
inference_headers = inference.h amalgamation.h tree_search.h tree_prior.h
inference_sources = inference.c amalgamation.c tree_search.c tree_prior.c
noinst_HEADERS = parsimony_internal.h substitution_model_internal.h likelihood_kernel_internal.h \
  topology_counts_internal.h gaussian_internal.h

otherincludedir = $(includedir)/curupixa
otherinclude_HEADERS = curupixa.h $(phylogenetics_headers) $(gene_species_headers) $(inference_headers)

lib_LTLIBRARIES = libcurupixa_phylogenetics.la
libcurupixa_phylogenetics_la_SOURCES = $(phylogenetics_sources) $(phylogenetics_headers) quasi_random_constants.h
libcurupixa_phylogenetics_la_LIBADD = global/libcrpxglobal.la
libcurupixa_phylogenetics_la_LDFLAGS = -version-info 0:0:0 @OPENMP_CFLAGS@

if BUILD_GENE_SPECIES
lib_LTLIBRARIES += libcurupixa_gene_species.la
libcurupixa_gene_species_la_SOURCES = $(gene_species_sources) $(gene_species_headers)
libcurupixa_gene_species_la_LIBADD = libcurupixa_phylogenetics.la
libcurupixa_gene_species_la_LDFLAGS = -version-info 0:0:0 @OPENMP_CFLAGS@
endif
if BUILD_INFERENCE
lib_LTLIBRARIES += libcurupixa_inference.la
libcurupixa_inference_la_SOURCES = $(inference_sources) $(inference_headers)
libcurupixa_inference_la_LIBADD = libcurupixa_gene_species.la
libcurupixa_inference_la_LDFLAGS = -version-info 0:0:0 @OPENMP_CFLAGS@
endif
