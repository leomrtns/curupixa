/* SPDX-License-Identifier: GPL-3.0-or-later */
/*! \file union_find_internal.h
    Path compression and union by component size, shared by Kruskal and density hierarchies. */
#ifndef CRPX_UNION_FIND_INTERNAL_H
#define CRPX_UNION_FIND_INTERNAL_H
#include <stddef.h>

static inline size_t
crpx_component_root (size_t *parent, size_t item)
{
  size_t root = item;
  while (parent[root] != root) root = parent[root];
  while (parent[item] != item) { size_t next = parent[item]; parent[item] = root; item = next; }
  return root;
}

static inline size_t
crpx_component_join (size_t *parent, size_t *sizes, size_t first, size_t second)
{ // Arguments are distinct roots; return the root of the joined component.
  if (sizes[first] < sizes[second] || (sizes[first] == sizes[second] && first > second)) {
    size_t swap = first; first = second; second = swap;
  }
  parent[second] = first;
  sizes[first] += sizes[second];
  return first;
}
#endif
