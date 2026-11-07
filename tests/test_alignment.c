/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifdef NDEBUG
#undef NDEBUG
#endif
#include "alignment.h"
#include "file_compression.h"
#include "nexus_common.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define OK(expression) assert((expression) == CRPX_OK)

static void test_patterns(void)
{
  crpx_alignment *a;
  OK(crpx_alignment_from_text(">taxon one\nacgt acgt NN\n>taxon two\nacgt ACGA RY\n", true, &a));
  assert(crpx_alignment_taxa(a) == 2 && crpx_alignment_sites(a) == 10);
  assert(crpx_alignment_patterns(a) == 7);
  const unsigned char *patterns = crpx_alignment_pattern_data(a);
  const size_t *mapping = crpx_alignment_site_patterns(a), *weights = crpx_alignment_pattern_weights(a);
  size_t sum = 0, index;
  for (size_t p = 0; p < crpx_alignment_patterns(a); ++p) sum += weights[p];
  assert(sum == 10);
  for (size_t site = 0; site < 10; ++site) for (size_t taxon = 0; taxon < 2; ++taxon)
    assert(patterns[mapping[site]*2+taxon] == (unsigned char)crpx_alignment_sequence(a, taxon)[site]);
  OK(crpx_alignment_find_taxon(a, "taxon two", &index)); assert(index == 1);
  double distance;
  size_t used;
  OK(crpx_alignment_pairwise_distance(a, 0, 1, CRPX_DISTANCE_P, &distance, &used));
  assert(distance == 0.125 && used == 8);
  OK(crpx_alignment_pairwise_distance(a, 0, 1, CRPX_DISTANCE_JC69, &distance, &used));
  assert(fabs(distance + 0.75*log(1-4.0/24.0)) < 1e-14);
  double states[28];
  OK(crpx_alignment_tip_states(a, 1, states, 28));
  size_t ambiguous = mapping[8];
  assert(states[4*ambiguous] == 1 && states[4*ambiguous+1] == 0 && states[4*ambiguous+2] == 1);
  crpx_alignment_free(a);
  OK(crpx_alignment_from_text(">x\nAAA\n>y\nCCC\n", true, &a));
  OK(crpx_alignment_pairwise_distance(a, 0, 1, CRPX_DISTANCE_K2P, &distance, &used)); assert(isinf(distance));
  crpx_alignment_free(a);
  assert(crpx_alignment_from_text(">x\nAAA\n>x\nAAA\n", false, &a) == CRPX_ERROR_INVALID_ARGUMENT && !a);
  OK(crpx_alignment_from_text(">x\nA\n>y\nAC\n", false, &a));
  assert(!crpx_alignment_is_aligned(a)); crpx_alignment_free(a);
  assert(crpx_alignment_from_text(">x\nA\n>y\nAC\n", true, &a) == CRPX_ERROR_INVALID_ARGUMENT && !a);
}

static void test_nexus(void)
{
  crpx_alignment *a;
  const char *input = "#NEXUS\nBEGIN DATA;\nDIMENSIONS NTAX=2 NCHAR=8;\n"
    "FORMAT DATATYPE=DNA GAP=- MISSING=? INTERLEAVE=YES;\nMATRIX\n"
    "'taxon [one]' AC[comment[nested]]GT\n'second ''taxon' AGGT\n\n"
    "'taxon [one]' ACGT\n'second ''taxon' .C..\n;\nEND;";
  OK(crpx_alignment_from_text(input, true, &a));
  assert(!strcmp(crpx_alignment_name(a, 0), "taxon [one]"));
  assert(!strcmp(crpx_alignment_name(a, 1), "second 'taxon"));
  assert(!strcmp(crpx_alignment_sequence(a, 0), "ACGTACGT"));
  assert(!strcmp(crpx_alignment_sequence(a, 1), "AGGTACGT"));
  crpx_alignment_free(a);
  assert(crpx_alignment_from_text("#NEXUS\nBEGIN DATA; DIMENSIONS NTAX=1 NCHAR=2; MATRIX\nx AC\n", true, &a)
         == CRPX_ERROR_PARSE && !a);
  assert(crpx_alignment_from_text("#NEXUS\nBEGIN DATA; DIMENSIONS NTAX=1 NCHAR=2; FORMAT TRANSPOSE; MATRIX\nx AC\n;",
                                  true, &a) == CRPX_ERROR_UNSUPPORTED && !a);
  OK(crpx_alignment_from_text("# leading FASTA comment\n; another comment\n>A\nAC\n>B\nAG\n", true, &a));
  assert(!strcmp(crpx_alignment_sequence(a, 1), "AG")); crpx_alignment_free(a);
  OK(crpx_alignment_from_text("#NEXUS\nBEGIN TREES; TREE MATRIX = (a,b); END;\n"
    "BEGIN DATA; DIMENSIONS NTAX=2 NCHAR=8; FORMAT DATATYPE=DNA GAP='-' MISSING='?'; MATRIX\n"
    "a AC GT\nAC GT\nb AG GT AG GT\n; END;", true, &a));
  assert(!strcmp(crpx_alignment_sequence(a, 0), "ACGTACGT"));
  assert(!strcmp(crpx_alignment_sequence(a, 1), "AGGTAGGT")); crpx_alignment_free(a);
  const char *unsupported[] = {
    "#NEXUS BEGIN DATA; DIMENSIONS NTAX=1 NCHAR=2; FORMAT DATATYPE=PROTEIN; MATRIX\nx AC\n; END;",
    "#NEXUS BEGIN DATA; DIMENSIONS NTAX=1 NCHAR=2; FORMAT EQUATE=R; MATRIX\nx AC\n; END;",
    "#NEXUS BEGIN DATA; DIMENSIONS NTAX=1 NCHAR=2; FORMAT RESPECTCASE; MATRIX\nx AC\n; END;",
    "#NEXUS BEGIN DATA; DIMENSIONS NTAX=1 NCHAR=2; MATRIX\nx AC\n; MATRIX\nx GT\n; END;"
  };
  for (size_t i = 0; i < sizeof(unsupported)/sizeof(*unsupported); i++)
    assert(crpx_alignment_from_text(unsupported[i], true, &a) == CRPX_ERROR_UNSUPPORTED && !a);
  const char *malformed[] = {
    "#NEXUS BEGIN DATA; DIMENSIONS NTAX=1oops NCHAR=2; MATRIX\nx AC\n; END;",
    "#NEXUS BEGIN DATA; DIMENSIONS NTAX=1 NCHAR=2; MATRIX\nx AC\n;",
    "#NEXUS BEGIN DATA; DIMENSIONS NTAX=1 NCHAR=2; MATRIX\n'x'AC\n; END;"
  };
  for (size_t i = 0; i < sizeof(malformed)/sizeof(*malformed); i++)
    assert(crpx_alignment_from_text(malformed[i], true, &a) == CRPX_ERROR_PARSE && !a);
  char *copy = NULL;
  assert(crpx_nexus_strip_comments("unbalanced [", &copy, NULL) == CRPX_ERROR_PARSE && !copy);
}

/* CTest runs inside its project build directory. No temporary files leave it. */
static void test_io(void)
{
  const char *paths[] = {"alignment-test.fasta", "alignment-test.fasta.gz", "alignment-test.fasta.bz2"};
  crpx_alignment *a, *b;
  OK(crpx_alignment_from_text(">A\nACGTN\n>B\nACGT-\n", true, &a));
  const crpx_file_format formats[] = {CRPX_FILE_RAW, CRPX_FILE_GZIP, CRPX_FILE_BZIP2};
  for (size_t i = 0; i < 3; ++i) {
    if (!crpx_file_format_available(formats[i])) continue;
    OK(crpx_alignment_write_fasta(a, paths[i], 2));
    OK(crpx_alignment_read(paths[i], true, &b));
    assert(!strcmp(crpx_alignment_sequence(b, 1), "ACGT-"));
    assert(crpx_alignment_sites(b) == 5);
    crpx_alignment_free(b);
    assert(!remove(paths[i]));
  }
  crpx_alignment_free(a);
}

int main(void)
{
  test_patterns(); test_nexus(); test_io();
  puts("alignments, patterns and sequence I/O: passed");
  return 0;
}
