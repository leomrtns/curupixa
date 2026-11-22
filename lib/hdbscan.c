/* Copyright (C) 2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: MIT */
/*! \file hdbscan.c
    HDBSCAN*: core distances, a mutual-reachability MST, hierarchy condensation and excess-of-mass/leaf selection.
    Written from Campello, Moulavi, Zimek and Sander (2015), doi:10.1145/2733381. Equal-weight merges are collapsed
    together, avoiding artificial zero-lifetime clusters when observations or distances coincide. */
#include "hdbscan.h"
#include "clustering_internal.h"
#include "quickselect_quantile.h"
#include "union_find_internal.h"
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

struct density_distances { crpx_distance_generator_t generator; size_t component; double *core; };
struct density_node { size_t first, second, size; double weight; };
struct density_task { size_t node, cluster; };
struct density_work {
  size_t size, clusters;
  size_t *parents, *point_cluster, *children, *stack, *ancestor;
  long double *birth, *stability, *maximum, *point_lambda, *score;
  bool *choose, *has_children;
};
static crpx_status mutual_reachability (void *data, size_t first, size_t second, double *weight);
static crpx_status density_tree (const crpx_spanning_forest *mst, struct density_node **out);
static void density_condense (const struct density_node *nodes, size_t minimum, struct density_work *work,
                              struct density_task *tasks, crpx_hdbscan *result);
static void density_edge (struct density_work *work, crpx_hdbscan *result, size_t parent,
                          size_t child, size_t size, long double lambda);
static void density_select (struct density_work *work, crpx_hdbscan_selection selection, bool allow_root,
                            crpx_hdbscan *result);
static int density_compare (const void *first, const void *second);

void
crpx_hdbscan_free (crpx_hdbscan *result)
{
  if (!result) return;
  free (result->labels); free (result->selected_clusters); free (result->probabilities);
  free (result->stabilities); free (result->hierarchy); free (result);
}

crpx_status
crpx_hdbscan_run (crpx_distance_generator_t generator, size_t component, size_t min_points,
                  size_t min_cluster_size, crpx_hdbscan_selection selection, bool allow_single_cluster, crpx_hdbscan **out)
{ // Core distances feed the generic dense Prim kernel; no quadratic distance matrix is constructed here.
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  size_t size = crpx_distance_generator_samples (generator);
  if (!generator || !min_points || min_points > size || min_cluster_size < 2 || min_cluster_size > size ||
      component >= crpx_distance_generator_components (generator) || selection > CRPX_HDBSCAN_LEAF || selection < 0)
    return CRPX_ERROR_INVALID_ARGUMENT;
  if (size > SIZE_MAX / sizeof (double)) return CRPX_ERROR_OVERFLOW;
  double *core = malloc (size * sizeof (*core)), *row = malloc (size * sizeof (*row));
  if (!core || !row) { free (core); free (row); return CRPX_ERROR_ALLOCATION; }
  crpx_status status = CRPX_OK;
  for (size_t i = 0; i < size; i++) {
    if (min_points == 1) { core[i] = 0; continue; }
    for (size_t j = 0; j < size; j++) {
      status = crpx_cluster_distance (generator, i, j, component, row + j);
      if (status != CRPX_OK) goto done;
    }
    status = crpx_select_double (row, size, min_points - 1, core + i);
    if (status != CRPX_OK) goto done;
  }
  struct density_distances data = {generator, component, core};
  crpx_spanning_forest *mst = NULL;
  status = crpx_prim_dense (size, mutual_reachability, &data, &mst);
  if (status == CRPX_OK) status = crpx_hdbscan_from_mst (mst, min_cluster_size, selection, allow_single_cluster, out);
  crpx_spanning_forest_free (mst);
done:
  free (core); free (row);
  return status;
}

crpx_status
crpx_hdbscan_from_mst (const crpx_spanning_forest *mst, size_t min_cluster_size,
                       crpx_hdbscan_selection selection, bool allow_single_cluster, crpx_hdbscan **out)
{ // Condense a density hierarchy independently of how its MST was obtained.
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (!mst || mst->vertices < 2 || min_cluster_size < 2 || min_cluster_size > mst->vertices ||
      selection > CRPX_HDBSCAN_LEAF || selection < 0) return CRPX_ERROR_INVALID_ARGUMENT;
  size_t size = mst->vertices;
  if (size > (SIZE_MAX - 1) / 2 || 2 * size > SIZE_MAX / sizeof (crpx_hdbscan_edge) ||
      2 * size > SIZE_MAX / sizeof (struct density_node) || size > SIZE_MAX / sizeof (struct density_task) ||
      size > SIZE_MAX / sizeof (long double)) return CRPX_ERROR_OVERFLOW;
  struct density_node *nodes = NULL;
  crpx_status status = density_tree (mst, &nodes);
  if (status != CRPX_OK) return status;
  crpx_hdbscan *result = calloc (1, sizeof (*result));
  struct density_work work = {0};
  work.size = size; work.clusters = 1;
  work.parents = malloc (size * sizeof (*work.parents));
  work.point_cluster = malloc (size * sizeof (*work.point_cluster));
  work.children = malloc (size * sizeof (*work.children)); work.stack = malloc (size * sizeof (*work.stack));
  work.ancestor = malloc (size * sizeof (*work.ancestor));
  work.birth = calloc (size, sizeof (*work.birth)); work.stability = calloc (size, sizeof (*work.stability));
  work.maximum = calloc (size, sizeof (*work.maximum)); work.point_lambda = malloc (size * sizeof (*work.point_lambda));
  work.score = calloc (size, sizeof (*work.score)); work.choose = calloc (size, sizeof (*work.choose));
  work.has_children = calloc (size, sizeof (*work.has_children));
  struct density_task *tasks = malloc (size * sizeof (*tasks));
  if (!result || !work.parents || !work.point_cluster || !work.children || !work.stack || !work.ancestor ||
      !work.birth || !work.stability || !work.maximum || !work.point_lambda || !work.score || !work.choose || !work.has_children || !tasks) {
    status = CRPX_ERROR_ALLOCATION; goto done;
  }
  result->samples = size;
  result->labels = malloc (size * sizeof (*result->labels));
  result->selected_clusters = malloc (size * sizeof (*result->selected_clusters));
  result->probabilities = calloc (size, sizeof (*result->probabilities));
  result->stabilities = malloc (size * sizeof (*result->stabilities));
  result->hierarchy = malloc ((2 * size - 2) * sizeof (*result->hierarchy));
  if (!result->labels || !result->selected_clusters || !result->probabilities ||
      !result->stabilities || !result->hierarchy) {
    status = CRPX_ERROR_ALLOCATION; goto done;
  }
  work.parents[0] = SIZE_MAX;
  density_condense (nodes, min_cluster_size, &work, tasks, result);
  density_select (&work, selection, allow_single_cluster, result);
done:
  free (nodes); free (tasks); free (work.parents); free (work.point_cluster); free (work.children); free (work.stack);
  free (work.ancestor); free (work.birth); free (work.stability); free (work.maximum); free (work.point_lambda);
  free (work.score); free (work.choose); free (work.has_children);
  if (status != CRPX_OK) crpx_hdbscan_free (result);
  else *out = result;
  return status;
}

static crpx_status
mutual_reachability (void *data, size_t first, size_t second, double *weight)
{
  struct density_distances *source = data;
  crpx_status status = crpx_cluster_distance (source->generator, first, second, source->component, weight);
  if (status == CRPX_OK) *weight = fmax (*weight, fmax (source->core[first], source->core[second]));
  return status;
}

static crpx_status
density_tree (const crpx_spanning_forest *mst, struct density_node **out)
{ // Kruskal's merge order supplies the binary hierarchy; condensation flattens merges at equal weights.
  *out = NULL;
  size_t size = mst->vertices;
  if (mst->components != 1 || mst->edge_count != size - 1 || !mst->edges) return CRPX_ERROR_TOPOLOGY;
  crpx_graph_edge *edges = malloc ((size - 1) * sizeof (*edges));
  size_t *parent = malloc (size * sizeof (*parent)), *sizes = malloc (size * sizeof (*sizes));
  size_t *component_node = malloc (size * sizeof (*component_node));
  struct density_node *nodes = calloc (2 * size - 1, sizeof (*nodes));
  crpx_status status = CRPX_OK;
  if (!edges || !parent || !sizes || !component_node || !nodes) { status = CRPX_ERROR_ALLOCATION; goto done; }
  for (size_t i = 0; i < size - 1; i++) {
    edges[i] = mst->edges[i];
    if (edges[i].first >= size || edges[i].second >= size || !isfinite (edges[i].weight) || edges[i].weight < 0) {
      status = CRPX_ERROR_INVALID_ARGUMENT; goto done;
    }
    if (edges[i].first > edges[i].second) {
      size_t swap = edges[i].first; edges[i].first = edges[i].second; edges[i].second = swap;
    }
  }
  qsort (edges, size - 1, sizeof (*edges), density_compare);
  for (size_t i = 0; i < size; i++) { parent[i] = component_node[i] = i; sizes[i] = nodes[i].size = 1; }
  for (size_t i = 0; i < size - 1; i++) {
    size_t first = crpx_component_root (parent, edges[i].first), second = crpx_component_root (parent, edges[i].second);
    if (first == second) { status = CRPX_ERROR_TOPOLOGY; goto done; }
    size_t node = size + i;
    nodes[node] = (struct density_node){component_node[first], component_node[second],
                                        sizes[first] + sizes[second], edges[i].weight};
    size_t root = crpx_component_join (parent, sizes, first, second);
    component_node[root] = node;
  }
  *out = nodes; nodes = NULL;
done:
  free (edges); free (parent); free (sizes); free (component_node); free (nodes);
  return status;
}

static void
density_edge (struct density_work *work, crpx_hdbscan *result, size_t parent,
               size_t child, size_t size, long double lambda)
{
  result->hierarchy[result->edge_count++] = (crpx_hdbscan_edge){work->size + parent, child, size, (double)lambda};
  long double lifetime = lambda == work->birth[parent] ? 0 : lambda - work->birth[parent];
  work->stability[parent] += lifetime * size;
  if (child >= work->size) work->has_children[parent] = true;
  if (lambda > work->maximum[parent]) work->maximum[parent] = lambda;
}

static void
density_condense (const struct density_node *nodes, size_t minimum, struct density_work *work,
                    struct density_task *tasks, crpx_hdbscan *result)
{
  size_t pending = 1, size = work->size;
  tasks[0] = (struct density_task){2 * size - 2, 0};
  while (pending) {
    struct density_task task = tasks[--pending];
    double weight = nodes[task.node].weight;
    long double lambda = weight == 0 ? INFINITY : 1.0L / weight;
    size_t count = 0, stack_size = 1, large = 0;
    work->stack[0] = task.node;
    while (stack_size) {
      size_t node = work->stack[--stack_size];
      if (node >= size && nodes[node].weight == weight) {
        work->stack[stack_size++] = nodes[node].first; work->stack[stack_size++] = nodes[node].second;
      } else { work->children[count++] = node; if (nodes[node].size >= minimum) large++; }
    }
    for (size_t i = 0; i < count; i++) {
      size_t child = work->children[i];
      if (nodes[child].size >= minimum) {
        size_t cluster = task.cluster;
        if (large > 1) {
          cluster = work->clusters++;
          work->parents[cluster] = task.cluster; work->birth[cluster] = lambda;
          density_edge (work, result, task.cluster, size + cluster, nodes[child].size, lambda);
        }
        tasks[pending++] = (struct density_task){child, cluster};
      } else {
        stack_size = 1; work->stack[0] = child;
        while (stack_size) {
          size_t node = work->stack[--stack_size];
          if (node >= size) { work->stack[stack_size++] = nodes[node].first; work->stack[stack_size++] = nodes[node].second; }
          else {
            work->point_cluster[node] = task.cluster; work->point_lambda[node] = lambda;
            density_edge (work, result, task.cluster, node, 1, lambda);
          }
        }
      }
    }
  }
}

static void
density_select (struct density_work *work, crpx_hdbscan_selection selection, bool allow_root, crpx_hdbscan *result)
{ // Bottom-up dynamic programming chooses stability; top-down propagation enforces non-overlapping membership.
  size_t count = work->clusters, size = work->size;
  for (size_t position = count; position-- > 0;) {
    bool has_children = work->has_children[position];
    bool eligible = position != 0 || (allow_root && selection == CRPX_HDBSCAN_EOM);
    bool choose = eligible && (selection == CRPX_HDBSCAN_LEAF ? !has_children :
                               work->stability[position] > 0 && work->stability[position] >= work->score[position]);
    work->choose[position] = choose;
    if (choose) work->score[position] = work->stability[position];
    if (position) {
      size_t parent = work->parents[position];
      work->score[parent] += work->score[position];
    }
  }
  for (size_t i = 0; i < count; i++) {
    size_t inherited = i ? work->ancestor[work->parents[i]] : SIZE_MAX;
    work->ancestor[i] = inherited != SIZE_MAX ? inherited : (work->choose[i] ? i : SIZE_MAX);
    work->children[i] = SIZE_MAX;
  }
  for (size_t i = 0; i < size; i++) {
    size_t cluster = work->ancestor[work->point_cluster[i]];
    result->labels[i] = SIZE_MAX;
    if (cluster == SIZE_MAX || (!cluster && work->point_lambda[i] < work->maximum[0])) continue;
    if (work->children[cluster] == SIZE_MAX) {
      size_t label = result->clusters++;
      work->children[cluster] = label; result->selected_clusters[label] = size + cluster;
      result->stabilities[label] = (double)work->stability[cluster];
    }
    result->labels[i] = work->children[cluster];
    long double maximum = work->maximum[cluster], value = work->point_lambda[i];
    result->probabilities[i] = isinf (value) || maximum == 0 ? 1 : (double)(fminl (value, maximum) / maximum);
  }
}

static int
density_compare (const void *first, const void *second)
{
  const crpx_graph_edge *a = first, *b = second;
  if (a->weight != b->weight) return a->weight < b->weight ? -1 : 1;
  if (a->first != b->first) return a->first < b->first ? -1 : 1;
  return (a->second > b->second) - (a->second < b->second);
}
