AM_CPPFLAGS = -I$(top_srcdir)/lib -I$(top_builddir)/lib/global
AM_CFLAGS = @AM_CFLAGS@ @OPENMP_CFLAGS@ -UNDEBUG
LDADD = ../lib/libcurupixa_static.la @OPENMP_CFLAGS@
check_PROGRAMS = test_core test_tree test_runtime test_statistics test_probability test_hashing \
  test_suffix_index test_topology test_alignment test_compression test_gff3 test_phylogenetics \
  test_minimiser test_svd test_clustering test_algorithms test_character_analysis \
  test_substitution_model test_likelihood test_likelihood_prepared test_continuous_traits test_reconstruction test_tree_extensions test_tree_allocation test_spr test_reconciliation_roots
TESTS = $(check_PROGRAMS)
