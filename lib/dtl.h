/* Copyright (C) 2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef CURUPIXA_DTL_H
#define CURUPIXA_DTL_H
#include "tree.h"
#ifdef __cplusplus
extern "C" {
#endif

typedef struct crpx_dtl_species crpx_dtl_species;
typedef struct crpx_dtl_workspace crpx_dtl_workspace;
typedef struct { double duplication, transfer, loss; } crpx_dtl_costs;
typedef enum { CRPX_DTL_MINIMUM_COST, CRPX_DTL_LOG_PARTITION } crpx_dtl_mode;
typedef enum {
  CRPX_DTL_TIP, CRPX_DTL_SPECIATION, CRPX_DTL_DUPLICATION, CRPX_DTL_TRANSFER
} crpx_dtl_event;
typedef struct {
  crpx_dtl_event event;
  crpx_node_id species, recipient;
} crpx_dtl_mapping;
typedef struct { uint64_t duplications, transfers, losses; } crpx_dtl_counts;

/* Prepare an immutable rooted binary species snapshot. NULL ages selects an
 * undated model: transfer recipients must be incomparable species branches;
 * global temporal consistency is not imposed. Otherwise ages[node] are finite
 * times before the present, leaves have age zero and parents are strictly older.
 * origin_age must be strictly greater than the root age. Edges are sliced at
 * every distinct species age; transfers join segments in the same open interval.
 * A finite series of gene events within one interval admits times in gene order.
 * The origin is fixed at the top of the species root stem. No inputs are retained.
 * Undated origin_age is ignored. Branch lengths and labels are ignored. */
crpx_status crpx_dtl_species_new (const crpx_tree *species, const double *ages, double origin_age,
                                 crpx_dtl_species **result);
void crpx_dtl_species_free (crpx_dtl_species *species);
size_t crpx_dtl_species_states (const crpx_dtl_species *species);
size_t crpx_dtl_species_storage_size (const crpx_dtl_species *species);
crpx_status crpx_dtl_workspace_new (crpx_dtl_workspace **result);
void crpx_dtl_workspace_free (crpx_dtl_workspace *workspace);
void crpx_dtl_workspace_clear (crpx_dtl_workspace *workspace);
size_t crpx_dtl_workspace_storage_size (const crpx_dtl_workspace *workspace);

/* Finite observed-tree DTL model, with explicit tip mappings (gene-node indexed,
 * leaf entries only). Gene and species trees must be rooted and binary. All
 * costs are finite and nonnegative. Duplication/transfer consume a gene branching;
 * speciation either splits the gene children between species children, or passes
 * one gene lineage to one child and pays one loss for the other child.
 * Unobserved duplication, transfer-loss chains and extinct subtrees are excluded.
 * MINIMUM_COST returns the least total event cost; LOG_PARTITION returns
 * log(sum_history exp(-total_cost)). This is a weighted finite-history partition,
 * not a continuous-time birth-death likelihood. Dated histories distinguish
 * event segments, not continuous event times; changing the slicing changes the
 * partition. Gene children and gene copies retain distinct identities.
 * Mutable capacity is retained for reuse; one workspace per concurrent call.
 * Prepared species must remain alive until the next evaluation/clear/free if
 * traceback is required. The gene tree and map are not retained. Outputs change
 * only on success; failure invalidates the preceding traceback.
 * Memory O(gene_nodes * species_states); work additionally scans contemporaneous
 * recipient states for each gene branching and donor state. */
crpx_status crpx_dtl_evaluate (const crpx_tree *gene, const crpx_dtl_species *species,
                              const crpx_node_id *tip_map, const crpx_dtl_costs *costs,
                              crpx_dtl_mode mode, crpx_dtl_workspace *workspace, double *value);
/* One minimum-cost history, with first encountered tie. Mapping has one entry
 * per gene node; species identifies an incoming species branch (root stem for
 * species root). recipient is NONE except at a transfer. Loss counts include
 * losses between the fixed origin and the gene root. MINIMUM_COST mode required.
 * Output capacity must cover gene_node_count; no allocation is needed. */
crpx_status crpx_dtl_traceback (crpx_dtl_workspace *workspace, crpx_dtl_mapping *mapping,
                               size_t capacity, crpx_dtl_counts *counts);
#ifdef __cplusplus
}
#endif
#endif
