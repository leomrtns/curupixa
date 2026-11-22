/* Copyright (C) 2026 Leonardo de Oliveira Martins
   SPDX-License-Identifier: GPL-3.0-or-later */
#include "rooting.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

struct crpx_rooting_workspace {
  size_t capacity, count;
  crpx_node_id *order, *parent, *original, *first, *next;
  long double *length, *distance, *depth, *sum, *square, *a, *b, *value;
  size_t *tips;
};
static crpx_status reserve (crpx_rooting_workspace *w, size_t n);
static void traverse (crpx_rooting_workspace *w, crpx_node_id source);
static crpx_status locate (const crpx_tree *tree, crpx_rooting_workspace *workspace,
                           crpx_root_location *out, unsigned method);
static long double clamp (long double x, long double length);

crpx_status
crpx_rooting_workspace_new (crpx_rooting_workspace **out)
{
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = calloc (1, sizeof (**out));
  return *out ? CRPX_OK : CRPX_ERROR_ALLOCATION;
}

void
crpx_rooting_workspace_free (crpx_rooting_workspace *w)
{
  if (!w) return;
  free (w->order); free (w->parent); free (w->original); free (w->first); free (w->next);
  free (w->length); free (w->distance); free (w->depth); free (w->sum); free (w->square);
  free (w->a); free (w->b); free (w->value); free (w->tips); free (w);
}

crpx_status
crpx_tree_root_midpoint (const crpx_tree *tree, crpx_rooting_workspace *workspace, crpx_root_location *out)
{ /* Diameter centre, equivalently the minimum maximum tip distance; Mai et al. (2017). */
  return locate (tree, workspace, out, 0);
}

crpx_status
crpx_tree_root_minvar (const crpx_tree *tree, crpx_rooting_workspace *workspace, crpx_root_location *out)
{ /* Reroot distance moments and minimise their quadratic on each edge; Mai et al. (2017). */
  return locate (tree, workspace, out, 1);
}

crpx_status
crpx_tree_root_mad (const crpx_tree *tree, crpx_rooting_workspace *workspace, crpx_root_location *out)
{ /* Accumulate edge quadratics for the A1-B1 objective of Tria et al. (2017). */
  return locate (tree, workspace, out, 2);
}

static crpx_status
locate (const crpx_tree *tree, crpx_rooting_workspace *workspace, crpx_root_location *out, unsigned method)
{
  if (!tree || !out) return CRPX_ERROR_INVALID_ARGUMENT;
  crpx_status status = crpx_tree_validate_policy (tree, CRPX_TREE_REQUIRE_COMPLETE_LENGTHS |
                                                      CRPX_TREE_REQUIRE_NONNEGATIVE_LENGTHS);
  if (status != CRPX_OK) return status;
  size_t n = crpx_tree_node_count (tree), ntips = crpx_tree_leaf_count (tree);
  crpx_node_id root = crpx_tree_root (tree), child = crpx_tree_first_child (tree, root);
  if (ntips < 2 || child == CRPX_NODE_NONE) return CRPX_ERROR_INVALID_ARGUMENT;
  if (crpx_tree_next_sibling (tree, child) == CRPX_NODE_NONE) return CRPX_ERROR_UNSUPPORTED;
  crpx_rooting_workspace *w = workspace;
  if (!w && (status = crpx_rooting_workspace_new (&w)) != CRPX_OK) return status;
  status = reserve (w, n);
  if (status != CRPX_OK) goto done;
  w->count = n;
  long double scale = 0;
  for (crpx_node_id v = 0; v < n; v++) {
    double length = 0;
    w->original[v] = crpx_tree_parent (tree, v);
    w->first[v] = crpx_tree_first_child (tree, v);
    w->next[v] = crpx_tree_next_sibling (tree, v);
    if (v != root) {
      status = crpx_tree_branch_length (tree, v, &length);
      if (status != CRPX_OK) goto done;
    }
    w->length[v] = length;
    if (length > scale) scale = length;
  }
  if (!scale) scale = 1;
  for (size_t v = 0; v < n; v++) w->length[v] /= scale;
  traverse (w, root);
  memcpy (w->depth, w->distance, n * sizeof (*w->depth));
  memset (w->a, 0, n * sizeof (*w->a));
  memset (w->b, 0, n * sizeof (*w->b));
  long double base = 0, denominator = 1;
  crpx_node_id endpoint = child;
  if (!method) {
    for (size_t pass = 0; pass < 2; pass++) {
      endpoint = CRPX_NODE_NONE;
      for (crpx_node_id v = 0; v < n; v++) if (w->first[v] == CRPX_NODE_NONE &&
          (endpoint == CRPX_NODE_NONE || w->distance[v] > w->distance[endpoint])) endpoint = v;
      if (!pass) traverse (w, endpoint);
    }
    /* For every edge, the farthest distances from the two diameter endpoints determine its eccentricity. */
    memcpy (w->sum, w->distance, n * sizeof (*w->sum));
    base = w->distance[endpoint] / 2;
    traverse (w, endpoint);
    memcpy (w->square, w->distance, n * sizeof (*w->square));
    traverse (w, root);
  } else if (method == 1) {
    memset (w->sum, 0, n * sizeof (*w->sum));
    memset (w->square, 0, n * sizeof (*w->square));
    memset (w->tips, 0, n * sizeof (*w->tips));
    for (size_t k = n; k-- > 0;) {
      crpx_node_id v = w->order[k], p = w->original[v];
      if (w->first[v] == CRPX_NODE_NONE) w->tips[v] = 1;
      if (p != CRPX_NODE_NONE) {
        long double l = w->length[v];
        w->tips[p] += w->tips[v];
        w->sum[p] += w->sum[v] + l * w->tips[v];
        w->square[p] += w->square[v] + 2 * l * w->sum[v] + l * l * w->tips[v];
      }
    }
    base = w->square[root] / ntips - (w->sum[root] / ntips) * (w->sum[root] / ntips);
    for (size_t k = 1; k < n; k++) {
      crpx_node_id v = w->order[k], p = w->original[v];
      long double l = w->length[v], inside = w->sum[v] + l * w->tips[v];
      long double slope = 1 - 2 * (long double)w->tips[v] / ntips;
      w->a[v] = 1 - slope * slope;
      w->b[v] = 2 * (w->sum[p] - 2 * inside) / ntips - 2 * w->sum[p] / ntips * slope;
      w->sum[v] = w->sum[p] + l * (ntips - 2 * (long double)w->tips[v]);
    }
  } else {
    denominator = 0;
    for (crpx_node_id tip = 0; tip < n; tip++) if (w->first[tip] == CRPX_NODE_NONE) {
      traverse (w, tip);
      memset (w->sum, 0, n * sizeof (*w->sum));
      memset (w->square, 0, n * sizeof (*w->square));
      for (crpx_node_id v = 0; v < n; v++) if (w->first[v] == CRPX_NODE_NONE && w->distance[v] > 0) {
        long double inverse = 1 / w->distance[v];
        w->sum[v] = inverse;
        w->square[v] = inverse * inverse;
        if (v > tip) {
          long double deviation = (w->depth[tip] - w->depth[v]) * inverse;
          base += deviation * deviation;
          denominator++;
        }
      }
      for (size_t k = n; k-- > 1;) {
        crpx_node_id v = w->order[k], p = w->parent[v];
        crpx_node_id edge = w->original[v] == p ? v : p;
        long double a = 4 * w->square[v];
        long double b = 4 * (2 * w->distance[p] * w->square[v] - w->sum[v]);
        w->a[edge] += a / 2;
        w->b[edge] += (edge == v ? b : -b - 2 * a * w->length[edge]) / 2;
        w->sum[p] += w->sum[v];
        w->square[p] += w->square[v];
      }
    }
    if (!denominator) { status = CRPX_ERROR_INVALID_ARGUMENT; goto done; }
    traverse (w, root);
  }
  w->value[root] = base;
  long double best = LDBL_MAX, fraction = 0;
  crpx_node_id selected = CRPX_NODE_NONE;
  for (size_t k = 1; k < n; k++) {
    crpx_node_id v = w->order[k];
    crpx_node_id p = w->original[v];
    long double l = w->length[v], x, score;
    if (!method) {
      /* At a diameter centre both endpoint distances equal diameter/2. */
      long double direction = w->sum[v] > w->sum[p] ? 1 : -1;
      x = clamp ((base - w->sum[p]) * direction, l);
      long double d1 = w->sum[p] + direction * x;
      long double direction2 = w->square[v] > w->square[p] ? 1 : -1;
      long double d2 = w->square[p] + direction2 * x;
      score = fmaxl (d1, d2);
    } else {
      long double a = w->a[v], b = w->b[v];
      x = a > 0 ? clamp (-b / (2 * a), l) : (b < 0 ? l : 0);
      score = w->value[p] + x * (b + a * x);
      w->value[v] = w->value[p] + l * (b + a * l);
    }
    long double f = l > 0 ? x / l : 0;
    if (!isfinite (score) || !isfinite (f)) { status = CRPX_ERROR_OVERFLOW; goto done; }
    if (score < best || (score == best && v < selected)) { best = score; selected = v; fraction = f; }
  }
  best = fmaxl (0, best);
  if (method == 2) best = sqrtl (best / denominator);
  else if (method == 1) best *= scale * scale;
  else best *= scale;
  if (!isfinite (best) || best > DBL_MAX) { status = CRPX_ERROR_OVERFLOW; goto done; }
  *out = (crpx_root_location){selected, (double)fraction, (double)best};
done:
  if (!workspace) crpx_rooting_workspace_free (w);
  return status;
}

static long double
clamp (long double x, long double length)
{
  return fminl (length, fmaxl (0, x));
}

static void
traverse (crpx_rooting_workspace *w, crpx_node_id source)
{ /* Iterative undirected traversal, using the stored orientation only to enumerate neighbours. */
  size_t count = 1;
  w->order[0] = source;
  w->parent[source] = CRPX_NODE_NONE;
  w->distance[source] = 0;
  for (size_t k = 0; k < count; k++) {
    crpx_node_id v = w->order[k], neighbour = w->original[v];
    if (neighbour != CRPX_NODE_NONE && neighbour != w->parent[v]) {
      w->parent[neighbour] = v;
      w->distance[neighbour] = w->distance[v] + w->length[v];
      w->order[count++] = neighbour;
    }
    for (neighbour = w->first[v]; neighbour != CRPX_NODE_NONE; neighbour = w->next[neighbour])
      if (neighbour != w->parent[v]) {
        w->parent[neighbour] = v;
        w->distance[neighbour] = w->distance[v] + w->length[neighbour];
        w->order[count++] = neighbour;
      }
  }
}

static crpx_status
reserve (crpx_rooting_workspace *w, size_t n)
{
  if (n <= w->capacity) return CRPX_OK;
  if (n > SIZE_MAX / sizeof (long double)) return CRPX_ERROR_OVERFLOW;
  /* Independent reallocations leave every successful allocation reusable after a later failure. */
#define GROW(field) do { void *memory = realloc (w->field, n * sizeof (*w->field)); \
  if (!memory) { return CRPX_ERROR_ALLOCATION; } w->field = memory; } while (0)
  GROW (order); GROW (parent); GROW (original); GROW (first); GROW (next);
  GROW (length); GROW (distance); GROW (depth); GROW (sum); GROW (square);
  GROW (a); GROW (b); GROW (value); GROW (tips);
#undef GROW
  w->capacity = n;
  return CRPX_OK;
}
