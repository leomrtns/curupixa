# Cladistic distances and reconciliation
These modules adapt the reusable split-pruning and reconciliation algorithms in guenomu and the newer genefam-dist. Their data belong to analysis objects; the 12-byte connectivity node has no additional fields. Include only the public headers required by the application:
```c
#include <curupixa/spr_distance.h>
#include <curupixa/reconciliation_roots.h>
#include <curupixa/taxon_map.h>
```
The [complete example](https://github.com/leomrtns/curupixa/blob/main/examples/gene_species_distances.c) shows ownership, explicit mapping, distances and root optimisation. The headers also work through `add_subdirectory()` using the unprefixed spelling, such as `#include <spr_distance.h>`.
## Choosing a distance
| Entry point | Result and assumptions |
| --- | --- |
| `crpx_tree_rf_distance()` in `topology_distance.h` | Exact rooted or unrooted RF for two trees with unique matching tip labels. |
| `crpx_splitset_compare()` in `splitset_distances.h` | Exact unrooted split RF and/or assignment h-distance on explicitly matched taxon bits. Positive smaller-side penalties for unmatched splits. Request RF alone to avoid assignment work. |
| `crpx_gene_species_split_distances()` | Same RF/positive-penalty assignment, with species leaves expanded into the gene-copy taxon space. |
| `crpx_tree_spr_distance()` in `spr_distance.h` | Unrooted greedy dSPR, exact RF, historical h-distance, common-cherry-reduced h-distance and a separate crossing-prune diagnostic. Unique matching labels. |
| `crpx_splitset_spr_distance()` | Same heuristic on precomputed splitsets sharing an explicit taxon bit order. Useful for repeated comparisons. |
| `crpx_gene_species_spr_distance()` | Same heuristic with an explicit gene-leaf to species-leaf map and expanded gene-copy taxon space. |

All distances ignore branch lengths, root stems and internal labels. Split extraction ignores trivial edges, deduplicates identical splits and handles polytomies without resolving them. One expanded species tip represents a multifurcation among its gene copies; absent species contribute no bits. Within-species gene relationships can consequently remain unresolved in the expanded species tree.
## Greedy dSPR and hdist conventions
The heuristic removes common splits, compresses common cherries, finds a minimum-cost matching between remaining splits, chooses a small matched disagreement, prunes it, and repeats. `spr` counts these pruning rounds. `spr_extra` counts rounds in which the selected prune straddles both sides of a remaining gene split. The caller decides whether and how to use that diagnostic.

This greedy estimate has no guarantee of equalling, bounding or defining a metric for exact SPR distance. Equal-cost assignment ties, taxon bit order and comparison direction can affect it. The tree adapter uses the first tree's leaf node order for bits; the explicit splitset interface gives the caller control of that order. It stops when either disagreeing split list is empty. An unresolved tree can therefore have positive RF and zero dSPR against a resolved tree.

The `hdist` and `reduced_hdist` fields preserve the historical **zero dummy-edge penalty** used within genefam-dist's dSPR algorithm. In contrast, `crpx_splitset_compare()` charges the smaller-side cardinality for an unmatched split. These agree when both trees have the same number of unique splits, including ordinary strictly binary trees with the same taxa. Use the positive-penalty interface when unequal resolution should incur a cost. Neither interface applies an implicit normalisation.

Curupixa corrects the reference's skipped one-pair assignment: discordant quartets have RF 2, h-distance 2 and dSPR 1. Tests against 160 seeded random pairs of binary trees with 4–12 taxa matched every RF result and every non-degenerate h-distance; the reference returned zero in 18 one-pair cases. SPR counts matched in 134 pairs. The remaining counts and crossing-prune diagnostics can differ because the new assignment solver and split processing choose different ties. These comparisons establish a semantic successor, with the changes stated here.
## Memory and parallel use
Construct one `crpx_spr_workspace` per concurrent worker. Construction allocates only its small control object. Subsequent calls reuse capacity; `crpx_spr_workspace_storage_size()` reports owned bytes excluding allocator overhead, and `clear()` releases retained capacity. No borrowed tree or splitset pointers survive a call. Multiple threads may read shared input trees/splitsets with separate workspaces; mutation and destruction require caller synchronisation.

For $a,b$ input splits and $T$ taxon bits, scratch costs $O((a+b)\lceil T/64\rceil + \max(a,b)^2)$ storage. It stores scalar assignment costs and recomputes bitsets only for selected matches, replacing the reference's all-pairs disagreement bitsets. One square assignment costs $O(\max(a,b)^3)$ time; pruning may require several assignments. This remains an optional expensive calculation. Tree adapters additionally extract descendant bitsets in $O(V\lceil T/64\rceil)$ storage. Precompute splitsets when appropriate, and use RF alone when assignment/SPR is unnecessary.
## Explicit taxon mapping
`crpx_gene_species_map_labels()` matches whole labels or one zero-based delimited field exactly and case-sensitively. For `Homo_sapiens|copy7`, delimiter `'|'` and field `0` select `Homo_sapiens`. Delimiter zero requires field zero and matches the whole label. Duplicate species tip names, empty fields and unmatched names return errors. On failure the caller's output map remains unchanged. On success internal gene entries contain `CRPX_NODE_NONE`; leaf entries identify species node IDs.

Applications can also construct this map themselves, for example from an annotation table. Each map must have one entry per gene node. Curupixa deliberately replaces the old case-insensitive longest-substring name guessing with explicit fields or caller-supplied identities. Species with overlapping names consequently cannot acquire silent accidental matches.
## Reconciliation at a supplied root
`reconciliation.h` provides `crpx_reconcile_lca()` for rooted, strictly binary gene and species trees. Its result owns the LCA map and per-species lineage counts, plus duplication, loss and extra-lineage totals. It uses the full species tree, within the gene family's mapped MRCA. It remains useful when the gene root is known or the downstream application needs the mapped nodes and lineage counts.
## Reconciliation over every gene root
`ancestry_index.h` provides an immutable connectivity snapshot with binary-lifting LCA queries. Building costs $O(S\log S)$ time/storage; each LCA costs $O(\log S)$. This replaces the reference's quadratic MRCA cache. The snapshot owns all its arrays, holds no tree pointer, and can be shared read-only among workers. Rebuild it after changing species connectivity; node IDs in maps must refer to the same snapshot. Labels and branch lengths do not affect it.

`crpx_reconcile_all_roots()` scores every edge of an unrooted binary gene tree using two directional message passes. Supply a strictly binary rooted gene tree as its representation: the function suppresses its artificial root and bridges the two children. It considers each original edge of that unrooted topology once, including terminal edges. Endpoints retain original gene node IDs; the bridge has the original root's children as endpoints. It returns scores rather than changing or allocating a rerooted tree. A single-tip gene has one degenerate zero-score candidate with both endpoints `CRPX_NODE_NONE`.

Choose the loss convention explicitly:
- `CRPX_RECONCILE_ALL_SPECIES` counts losses along the full species tree below the family's MRCA.
- `CRPX_RECONCILE_REPRESENTED_SPECIES` suppresses absent species and the resulting unary paths, reproducing the historical family-specific convention.

Extra lineages always use the represented-species convention and the legacy identity
$$
C = L_{\mathrm{represented}} - 2D + 2(n_{\mathrm{gene\ tips}}-n_{\mathrm{represented\ species}}).
$$
The candidate array reports all three scores together. `crpx_root_reconciliations_minima()` returns three **indices**, one for each objective. These indices can identify different rootings; combining the three minimum values does not describe a single reconciled tree. Choose an actual candidate when applying a weighted or joint objective. Polytomous reconciliation and dated transfers require a separately specified biological model and currently return `UNSUPPORTED` for these binary routines.

After building the reusable species index, all-root scoring uses $O(S+G\log S)$ work and $O(S+G)$ temporary storage. Result arrays own their scores and retain neither inputs nor the index. Independent calls can share the index and input trees while using ordinary separate results; no shared mutable reconciliation state exists.
## Other reusable material and boundaries
| Reference material | Curupixa position |
| --- | --- |
| guenomu/genefam splitsets, RF, assignment, dSPR, MRCA and reconciliation | Covered by the modules above, with explicit mapping and separate analysis ownership. |
| guenomu branch-integrated likelihood, patterns and phylogeny buffers | Redesigned as arbitrary-state character/model/likelihood objects; see [character analysis](character-analysis.md) and the [reference likelihood audit](likelihood-reference-audit.md). |
| Tree parsing, collections, distance construction, RNG, distributions, histograms and assignment | Already represented by independent curupixa modules; see the [migration inventory](migration.md). |
| Tree-sample consensus/support, credible-set summaries and indexed collection deduplication | Useful further generic work. Current collections retain integer sample weights and support exact RF compaction; these richer summaries remain gaps. |
| guenomu chain proposals, checkpoint layout, control parser and posterior hierarchy | Application-specific. A future reusable sampler or proposal layer should define its own state and callback contracts rather than attaching chain fields to trees. |

The local reference source files are `guenomu/lib/topology_splitset.c`, `guenomu/lib/topology_mrca.c` and their newer `genefam-dist/lib/` counterparts. Original authorship and GPL attribution remain in the adapted files. The implementations operate on the curupixa backbone and retain no dependency on either reference software.
