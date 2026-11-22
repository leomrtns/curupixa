/* SPDX-License-Identifier: MIT */
#include "paralog_quartets.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define OK(call) assert ((call) == CRPX_OK)

static crpx_node_id
leaf (const crpx_tree *tree, const char *name)
{
  for (size_t i = 0; i < crpx_tree_node_count (tree); i++) {
    const char *label = crpx_tree_label (tree, (crpx_node_id) i);
    if (label && !strcmp (label, name)) return (crpx_node_id) i;
  }
  return CRPX_NODE_NONE;
}

static void
check_family (const char *newick, bool two_loci, bool pectinate)
{
  crpx_tree *gene, *species, *other, *star;
  OK (crpx_tree_from_newick (newick, &gene, NULL));
  OK (crpx_tree_from_newick ("((A,B),(C,D));", &species, NULL));
  OK (crpx_tree_from_newick ("((A,C),(B,D));", &other, NULL));
  OK (crpx_tree_from_newick ("(A,B,C,D);", &star, NULL));
  size_t nodes = crpx_tree_node_count (gene);
  size_t *map = calloc (nodes, sizeof (*map)); uint8_t *tags = calloc (nodes, sizeof (*tags));
  assert (map && tags);
  for (size_t i = 0; i < nodes; i++) if (crpx_tree_first_child (gene, (crpx_node_id) i) == CRPX_NODE_NONE)
    map[i] = (size_t) (crpx_tree_label (gene, (crpx_node_id) i)[0] - 'A');
  crpx_node_id duplicate = two_loci ? crpx_tree_root (gene) : crpx_tree_parent (gene, leaf (gene, "A1"));
  tags[duplicate] = 1;
  crpx_paralog_quartets *quartets;
  OK (crpx_paralog_quartets_new (gene, map, 4, tags, &quartets));
  assert (quartets->count == (two_loci ? 2 : 1));
  if (pectinate) {
    crpx_node_id abc = crpx_tree_parent (gene, leaf (gene, "C"));
    assert (quartets->classes[0].anchor == abc && abc != crpx_tree_root (gene));
  }
  crpx_node_id tips[4];
  const char *names[] = {"A", "B", "C", "D"};
  crpx_quartet_index *index;
  for (size_t i = 0; i < 4; i++) tips[i] = leaf (species, names[i]);
  OK (crpx_quartet_index_new (species, tips, 4, &index));
  crpx_paralog_quartet_score score;
  OK (crpx_paralog_quartets_score (quartets, index, &score));
  assert (score.agreement == quartets->count && score.total == quartets->count);
  crpx_quartet_index_free (index);
  for (size_t i = 0; i < 4; i++) tips[i] = leaf (other, names[i]);
  OK (crpx_quartet_index_new (other, tips, 4, &index));
  OK (crpx_paralog_quartets_score (quartets, index, &score));
  assert (score.conflict == quartets->count && !score.agreement);
  crpx_quartet_index_free (index);
  for (size_t i = 0; i < 4; i++) tips[i] = leaf (star, names[i]);
  OK (crpx_quartet_index_new (star, tips, 4, &index));
  OK (crpx_paralog_quartets_score (quartets, index, &score));
  assert (score.unresolved_species == quartets->count);
  crpx_quartet_index_free (index); crpx_paralog_quartets_free (quartets);
  tags[duplicate] = 0;
  assert (crpx_paralog_quartets_new (gene, map, 4, tags, &quartets) == CRPX_ERROR_INVALID_ARGUMENT && !quartets);
  free (map); free (tags);
  crpx_tree_free (gene); crpx_tree_free (species); crpx_tree_free (other); crpx_tree_free (star);
}

int main (void)
{
  check_family ("(((A1,A2),B),(C,D));", false, false);
  check_family ("((((A1,A2),B),C),D);", false, true);
  check_family ("(((A1,B1),(C1,D1)),((A2,B2),(C2,D2)));", true, false);
  puts ("paralog quartet tests passed");
  return 0;
}
