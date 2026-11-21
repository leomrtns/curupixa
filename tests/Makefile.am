AM_CPPFLAGS = -I$(top_srcdir)/lib -I$(top_builddir)/lib/global
AM_CFLAGS = @AM_CFLAGS@ @OPENMP_CFLAGS@ -UNDEBUG
LDADD = ../lib/libcurupixa_phylogenetics.la @OPENMP_CFLAGS@
check_PROGRAMS = test_core test_tree test_runtime test_statistics test_probability test_hashing \
  test_suffix_index test_topology test_alignment test_compression test_gff3 \
  test_minimiser test_svd test_clustering test_character_analysis \
  test_substitution_model test_likelihood test_likelihood_prepared test_continuous_traits test_reconstruction test_tree_extensions test_tree_allocation
check_PROGRAMS += test_quartet test_tanglegram test_dtl
if BUILD_GENE_SPECIES
check_PROGRAMS += test_phylogenetics test_algorithms test_spr test_reconciliation_roots test_paralog_quartets
# These suites exercise both ordinary and explicitly mapped operations.
test_paralog_quartets_LDADD = ../lib/libcurupixa_gene_species.la $(LDADD)
test_phylogenetics_LDADD = ../lib/libcurupixa_gene_species.la $(LDADD)
test_algorithms_LDADD = ../lib/libcurupixa_gene_species.la $(LDADD)
test_spr_LDADD = ../lib/libcurupixa_gene_species.la $(LDADD)
test_reconciliation_roots_LDADD = ../lib/libcurupixa_gene_species.la $(LDADD)
endif
if BUILD_INFERENCE
check_PROGRAMS += test_amalgamation test_inference test_tree_search
test_tree_search_LDADD = ../lib/libcurupixa_inference.la $(LDADD)
test_amalgamation_LDADD = ../lib/libcurupixa_inference.la $(LDADD)
test_inference_LDADD = ../lib/libcurupixa_inference.la $(LDADD)
endif
TESTS = $(check_PROGRAMS)
