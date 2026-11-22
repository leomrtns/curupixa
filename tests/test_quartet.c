/* SPDX-License-Identifier: MIT */
#include "quartet_distance.h"
#include <assert.h>
#include <stdio.h>
#include <string.h>
#define OK(call) assert ((call) == CRPX_OK)

static void
map_tips (const crpx_tree *tree, const char *names, size_t taxa, crpx_node_id *tips)
{
  for (size_t i = 0; i < taxa; i++) {
    tips[i] = CRPX_NODE_NONE;
    for (size_t j = 0; j < crpx_tree_node_count (tree); j++) {
      const char *label = crpx_tree_label (tree, (crpx_node_id) j);
      if (label && label[0] == names[i] && !label[1]) tips[i] = (crpx_node_id) j;
    }
    assert (tips[i] != CRPX_NODE_NONE);
  }
}

int
main (void)
{
  crpx_tree *first, *second, *star, *unary;
  crpx_node_id a[5], b[5], c[5], d[5];
  crpx_quartet_counts counts;
  OK (crpx_tree_from_newick ("((A:1,B:2),(C:3,D:4));", &first, NULL));
  OK (crpx_tree_from_newick ("((A,C),(B,D));", &second, NULL));
  OK (crpx_tree_from_newick ("(D,C,B,A);", &star, NULL));
  OK (crpx_tree_from_newick ("(((A)),B,(C,D));", &unary, NULL));
  map_tips (first, "ABCD", 4, a);
  map_tips (second, "ABCD", 4, b);
  map_tips (star, "ABCD", 4, c);
  map_tips (unary, "ABCD", 4, d);
  OK (crpx_tree_quartet_distance (first, a, second, b, 4, &counts));
  assert (counts.total == 1 && counts.resolved_conflict == 1 && crpx_quartet_disagreements (&counts) == 1);
  OK (crpx_tree_quartet_distance (first, a, unary, d, 4, &counts));
  assert (counts.total == 1 && counts.resolved_agreement == 1);
  OK (crpx_tree_quartet_distance (first, a, star, c, 4, &counts));
  assert (counts.first_resolved_only == 1 && crpx_quartet_disagreements (&counts) == 1);
  OK (crpx_tree_quartet_distance (star, c, second, b, 4, &counts));
  assert (counts.second_resolved_only == 1);
  OK (crpx_tree_quartet_distance (star, c, star, c, 4, &counts));
  assert (counts.both_unresolved == 1 && crpx_quartet_disagreements (&counts) == 0);
  OK (crpx_tree_quartet_distance (first, a, second, b, 3, &counts));
  assert (!counts.total);
  crpx_quartet_index *index = NULL;
  a[1] = a[0];
  assert (crpx_quartet_index_new (first, a, 4, &index) == CRPX_ERROR_INVALID_ARGUMENT && !index);
  crpx_tree_free (first); crpx_tree_free (second); crpx_tree_free (star); crpx_tree_free (unary);
  OK (crpx_tree_from_newick ("(((A,B),C),(D,E));", &first, NULL));
  map_tips (first, "ABCDE", 5, a);
  OK (crpx_quartet_index_new (first, a, 5, &index));
  crpx_tree_free (first);
  OK (crpx_quartet_compare (index, index, &counts));
  assert (counts.total == 5 && counts.resolved_agreement == 5);
  crpx_quartet_index_free (index);
  puts ("quartet tests passed");
  return 0;
}
