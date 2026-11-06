/* Copyright (C) 2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#include "gff3_format.h"
#include <stdio.h>
#include <locale.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define CHECK(expression) do { if (!(expression)) { \
  fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, #expression); return 1; \
} } while (0)

static int check_document(const crpx_gff3 *document)
{
  CHECK(crpx_gff3_feature_count(document) == 4 && crpx_gff3_region_count(document) == 1);
  const crpx_gff3_region *region = crpx_gff3_region_at(document, 0);
  CHECK(region && !strcmp(region->seqid, "ctg 1") && region->start == 0 && region->end == 100);
  const crpx_gff3_feature *feature = crpx_gff3_feature_at(document, 0);
  CHECK(feature && !strcmp(feature->seqid, "ctg 1") && !strcmp(feature->source, "source with spaces"));
  CHECK(feature->start == 0 && feature->end == 100 && feature->has_score && feature->score == 1.25);
  CHECK(feature->phase == -1 && feature->strand == '+');
  size_t attribute;
  CHECK(crpx_gff3_find_attribute(document, 0, "Note", &attribute) == CRPX_OK);
  CHECK(crpx_gff3_attribute_value_count(document, 0, attribute) == 2);
  CHECK(!strcmp(crpx_gff3_attribute_value(document, 0, attribute, 0), "one,two"));
  CHECK(!strcmp(crpx_gff3_attribute_value(document, 0, attribute, 1), "three;four=five&six"));
  CHECK(crpx_gff3_find_attribute(document, 0, "Alias", &attribute) == CRPX_OK);
  CHECK(!strcmp(crpx_gff3_attribute_value(document, 0, attribute, 0), "a+b"));
  CHECK(crpx_gff3_find_attribute(document, 3, "Parent", &attribute) == CRPX_OK);
  CHECK(crpx_gff3_attribute_value_count(document, 3, attribute) == 2);
  size_t indices[4] = {99, 99, 99, 99}, count;
  CHECK(crpx_gff3_query(document, "ctg 1", 29, 30, NULL, indices, 0, &count) == CRPX_ERROR_OUT_OF_RANGE);
  CHECK(count == 1 && indices[0] == 99);
  CHECK(crpx_gff3_query(document, "ctg 1", 29, 30, NULL, indices, 4, &count) == CRPX_OK);
  CHECK(count == 1 && indices[0] == 0);
  CHECK(crpx_gff3_query(document, NULL, 20, 21, NULL, indices, 4, &count) == CRPX_OK && count == 2);
  CHECK(crpx_gff3_query(document, NULL, 0, 100, "CDS", indices, 4, &count) == CRPX_OK && count == 1 && indices[0] == 3);
  CHECK(crpx_gff3_query(document, "missing", 0, 1, NULL, NULL, 0, &count) == CRPX_OK && count == 0);
  const char *sequence; size_t length;
  CHECK(crpx_gff3_sequence(document, "ctg 1", &sequence, &length) == CRPX_OK && length == 100);
  CHECK(!strncmp(sequence, "ACGT", 4));
  return 0;
}

static int test_score_locale(void)
{
  const char *current = setlocale(LC_NUMERIC, NULL);
  char *saved = malloc(strlen(current) + 1);
  CHECK(saved);
  strcpy(saved, current);
  /* Only this isolated test process changes its locale; library code never does. */
  if (!setlocale(LC_NUMERIC, "en_DK.UTF-8")) setlocale(LC_NUMERIC, "en_DK.utf8");
  const char *text = "##gff-version 3\nseq\t.\tgene\t1\t2\t1.25\t+\t.\t.\n";
  crpx_gff3 *document = NULL, *copy = NULL;
  CHECK(crpx_gff3_from_text(text, &document) == CRPX_OK);
  CHECK(crpx_gff3_feature_at(document, 0)->score == 1.25);
  char *written = NULL;
  CHECK(crpx_gff3_to_text(document, &written) == CRPX_OK);
  CHECK(strstr(written, "\t1.25\t"));
  CHECK(crpx_gff3_from_text(written, &copy) == CRPX_OK);
  CHECK(crpx_gff3_feature_at(copy, 0)->score == 1.25);
  free(written); crpx_gff3_free(copy); crpx_gff3_free(document);
  const char *scores[] = {"-1.25e+2", "+.125E1", "1e-300", "5e-324", "1.7976931348623157e308", "-0"};
  const double expected[] = {-125, 1.25, 1e-300, 0x0.0000000000001p-1022, DBL_MAX, -0.0};
  for (size_t i = 0; i < sizeof(scores)/sizeof(*scores); i++) {
    char record[200];
    snprintf(record, sizeof(record), "##gff-version 3\nseq\t.\tgene\t1\t2\t%s\t+\t.\t.\n", scores[i]);
    CHECK(crpx_gff3_from_text(record, &document) == CRPX_OK);
    CHECK(crpx_gff3_feature_at(document, 0)->score == expected[i]);
    CHECK(crpx_gff3_to_text(document, &written) == CRPX_OK);
    CHECK(crpx_gff3_from_text(written, &copy) == CRPX_OK);
    CHECK(crpx_gff3_feature_at(copy, 0)->score == expected[i]);
    free(written); crpx_gff3_free(copy); crpx_gff3_free(document);
  }
  const char *invalid[] = {"0x1.4p0", "1,25", "1.2.3", "1e", "1e999999", "1e-999999", "NaN", "inf"};
  for (size_t i = 0; i < sizeof(invalid)/sizeof(*invalid); i++) {
    char record[200];
    snprintf(record, sizeof(record), "##gff-version 3\nseq\t.\tgene\t1\t2\t%s\t+\t.\t.\n", invalid[i]);
    CHECK(crpx_gff3_from_text(record, &document) == CRPX_ERROR_PARSE && !document);
  }
  CHECK(setlocale(LC_NUMERIC, saved));
  free(saved);
  return 0;
}

int main(void)
{
  CHECK(test_score_locale() == 0);
  unsigned char sentinel;
  crpx_gff3 *invalid_output = (crpx_gff3 *)&sentinel;
  CHECK(crpx_gff3_from_text(NULL, &invalid_output) == CRPX_ERROR_INVALID_ARGUMENT && !invalid_output);
  invalid_output = (crpx_gff3 *)&sentinel;
  CHECK(crpx_gff3_from_text("invalid", &invalid_output) == CRPX_ERROR_PARSE && !invalid_output);
  const char text[] =
    "##gff-version 3\r\n"
    "# annotation comment\n"
    "##custom-directive retained\n"
    "##sequence-region ctg%201 1 100\n"
    "ctg%201\tsource with spaces\tgene\t1\t100\t1.25\t+\t.\tID=g1;Note=one%2Ctwo,three%3Bfour%3Dfive%26six;Alias=a+b\n"
    "ctg%201\t.\texon\t20\t25\t.\t+\t.\tID=exon1;Parent=g1\n"
    "ctg%201\t.\tgene\t40\t60\t.\t?\t.\tID=g2\n"
    "###\n"
    "ctg%201\t.\tCDS\t2\t9\t.\t-\t0\tParent=g1,g2\n"
    "##FASTA\n"
    ">ctg%201 description retained\n"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgt\n";
  crpx_gff3 *document = NULL;
  CHECK(crpx_gff3_from_text(text, &document) == CRPX_OK);
  CHECK(check_document(document) == 0);
  char *written;
  CHECK(crpx_gff3_to_text(document, &written) == CRPX_OK);
  CHECK(strstr(written, "# annotation comment\n##custom-directive retained\n"));
  CHECK(strstr(written, "Note=one%2Ctwo,three%3Bfour%3Dfive%26six"));
  CHECK(strstr(written, "###\nctg%201"));
  crpx_gff3 *roundtrip;
  CHECK(crpx_gff3_from_text(written, &roundtrip) == CRPX_OK && check_document(roundtrip) == 0);
  free(written); crpx_gff3_free(roundtrip);
  const char *path = "curupixa-test-gff3.gff3";
  CHECK(crpx_gff3_write(document, path) == CRPX_OK);
  CHECK(crpx_gff3_read(path, &roundtrip) == CRPX_OK && check_document(roundtrip) == 0);
  remove(path); crpx_gff3_free(roundtrip);
  const char *names[] = {"ctg%201 first", "ctg%201 second"}, *sequences[] = {"ACGT", "ACGT"};
  crpx_alignment *duplicate;
  CHECK(crpx_alignment_new(names, sequences, 2, false, &duplicate) == CRPX_OK);
  CHECK(crpx_gff3_set_sequences(document, duplicate) == CRPX_ERROR_PARSE);
  CHECK(check_document(document) == 0);
  crpx_alignment_free(duplicate);
  crpx_gff3_free(document);
  const char *invalid[] = {
    "##gff-version 2\n",
    "seq\t.\tgene\t1\t2\t.\t+\t.\t.\n",
    "##gff-version 3\nseq\t.\tgene\t0\t2\t.\t+\t.\t.\n",
    "##gff-version 3\nseq\t.\tgene\t3\t2\t.\t+\t.\t.\n",
    "##gff-version 3\nseq\t.\tgene\t1\t2\tNaN\t+\t.\t.\n",
    "##gff-version 3\nseq\t.\tgene\t1\t2\t.\t+\t3\t.\n",
    "##gff-version 3\nseq\t.\tgene\t1\t2\t.\tx\t.\t.\n",
    "##gff-version 3\nseq\t.\tCDS\t1\t2\t.\t+\t.\t.\n",
    "##gff-version 3\nseq\t.\tgene\t1\t2\t.\t+\t.\tID=x;ID=y\n",
    "##gff-version 3\nseq\t.\tgene\t1\t2\t.\t+\t.\tID=%GG\n",
    "##gff-version 3\nseq\t.\tgene\t1\t2\t.\t+\t.\tID=%0\n",
    "##gff-version 3\nseq\t.\tgene\t1\t2\t.\t+\t.\tID=x;;Name=y\n",
    "##gff-version 3\n##sequence-region seq 2 1\n"
  };
  for (size_t i = 0; i < sizeof(invalid) / sizeof(invalid[0]); i++) {
    document = NULL;
    CHECK(crpx_gff3_from_text(invalid[i], &document) != CRPX_OK && !document);
  }
  CHECK(crpx_gff3_from_text("##gff-version 3\nseq\t.\tgene\t1\t2\t.\t+\t.\tID=%00\n", &document) == CRPX_ERROR_UNSUPPORTED);
  return 0;
}
