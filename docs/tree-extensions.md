# Tree roots, construction and operation-specific validation
The tree remains one compact connectivity object: each node has three 32-bit indices (parent, first child, next sibling), or 12 bytes on the tested ABI. Labels, lengths, traversal/sister arrays and analysis caches remain optional separate allocations. There are no exposed node pointers or per-node reference counters.
## Implementation roots and biological roots
Every nonempty curupixa tree has one ordinary vertex selected as its implementation root, with no parent. The first-child/next-sibling orientation permits any degree, including a bifurcation, trifurcation or unary vertex. Parsing preserves the root and vertices represented by Newick; it adds no fake tip, artificial vertex or three-struct representation per biological node.

Biological rootedness belongs to the analysis and the source interpretation. A binary rooted tree commonly has $2T-1$ vertices for $T$ tips, including its bifurcating root. The same unrooted binary topology can use an existing degree-three vertex as its implementation root, requiring $2T-2$ vertices. It can instead include an artificial binary root inserted along an edge. A persistent boolean would not change these vertex/edge counts or tell each algorithm which root convention is appropriate. RF already takes a rooted/unrooted argument; likelihood has an explicit artificial-binary-root convention. The new rooting function likewise requires an explicit structural policy.
## Changing size and retaining identities
Trees own contiguous arrays with separate count and reserved capacity. `crpx_tree_reserve()` allocates replacement arrays only when capacity must grow, copies existing entries, and commits them together after all allocations succeed. Appending a node uses spare capacity when available and otherwise grows geometrically. Applications expecting fixed-size rearrangements can reserve once. Existing NNI and SPR rearrangements reuse vertices; they do not increase node count.

Node IDs remain stable even when arrays move. Downstream code retains IDs and the tree owner, rather than addresses of private nodes. A tree revision changes after semantic mutations; optional workspaces must be refreshed, and immutable ancestry snapshots must be rebuilt for new connectivity. Stable IDs identify vertices, not a guarantee that a vertex keeps the same ancestors, descendants or biological interpretation after a rearrangement. There is no need to increment a reference counter for each node. The caller keeps borrowed input objects alive and synchronises mutation/destruction across threads.
## Rooting inside an edge
`crpx_tree_reroot_edge(tree, child, fraction, policy, &root)` places a binary root inside the edge from `child`'s parent to `child`. The fraction is measured from the parent. The policies are:
- `CRPX_TREE_REUSE_BINARY_ROOT`: suppress the old bifurcating root, join its two incident lengths, and reuse that vertex ID as the new root. This keeps node count fixed and leaves a strictly binary tree binary. The root label remains attached to its reused ID and therefore moves with it.
- `CRPX_TREE_KEEP_OLD_ROOT`: append one new root, keeping the old root vertex and its label at their former location. The old root can become unary. This also handles an implementation root of degree three or higher; a degree-three root becomes bifurcating after edge reorientation.

Neither policy deletes or renumbers existing vertices. Leaves retain their IDs and names. Internal labels remain vertex-attached; the downstream application must reassess support values or other annotations whose biological meaning depends on the original topology. Root stem metadata follows the new root and is excluded from pairwise tip distances. Unknown lengths propagate through sums/splits as unknown; signed lengths remain supported by the general representation, and unrepresentable sums return `OVERFLOW`.

The old root must have at least two children. Unary root stem chains return `UNSUPPORTED`, allowing the application to decide whether to remove, annotate or retain those vertices. Reuse also returns `UNSUPPORTED` when the old root has more than two children; select `KEEP_OLD_ROOT` to represent the additional root explicitly. There is no silent dropping of labels or unary chains.

Fractions zero and one still place a binary root on the selected edge, adding a zero-length segment when the original length is known. This preserves a tip even when rooting exactly at its edge endpoint. To root at an existing internal vertex without inserting an extra zero segment, use `crpx_tree_reroot()` instead. A zero-length edge permits any fraction but supplies no unique geometric placement. No automatic midpoint or diameter choice is made by this primitive.

The function prepares replacement connectivity and optional length/label-pointer arrays in $O(\mathrm{capacity})$ scratch, preserving the original tree, capacity, revision and output argument if any check/allocation fails. It copies label pointers without duplicating their strings and commits only after reorientation succeeds. The old tree object remains the owner. Reuse keeps the node count fixed; this transactional implementation trades temporary capacity-proportional storage for a simple failure contract. Numerical analyses that do not root trees incur none of this scratch cost.
## Bulk appending without permanent per-node tails
Ordinary `crpx_tree_add_node()` remains convenient for bounded-degree construction; appending to a parent scans its siblings. Newick parsing already tracks construction tails and remains linear. For public bulk construction, use an optional `crpx_tree_builder`:
```c
crpx_tree_builder *builder;
crpx_tree *tree;
crpx_node_id root, tip;
/* Check all returned statuses in application code. */
crpx_tree_new(&tree);
crpx_tree_builder_new(&builder);
crpx_tree_reserve(tree, 4001);
crpx_tree_builder_append(builder, tree, CRPX_NODE_NONE, &root);
for (size_t i = 0; i < 4000; ++i)
  crpx_tree_builder_append(builder, tree, root, &tip);
crpx_tree_builder_free(builder);
crpx_tree_free(tree);
```
The first bind or external topology mutation requires an $O(V)$ validation/tail refresh. Appends through that builder then use the cached tail in constant sibling work, with amortised allocation growth. Child order is preserved. Geometric tail capacity belongs to the builder, never to every tree; `storage_size()` reports it. Keep the bound tree alive until resetting, rebinding or freeing the builder. Reset drops its borrowed identity while retaining tail capacity. A private topology generation (one 64-bit value per tree) separates connectivity changes from ordinary semantic revisions: interleaving label/length annotation with every append keeps the cached tail path. External reparenting/rerooting/rearrangement triggers a safe refresh. This adds no per-node storage.

Polytomies do not require a second tree class: the existing three-index representation also stores a binary node in 12 bytes, and traversing children across the whole tree remains $O(V)$. Optional sister lookup costs one index per node only when requested; a unique sister exists for binary parents. Binary-specific reconciliation/proposal/PIC algorithms may still enforce bifurcation. Automatically resolving polytomies with zero edges would change split-based distances and parsimony assumptions, and zero edges can make continuous-trait models singular. Curupixa therefore preserves the represented topology and lets each algorithm specify its biological requirements.
## Reusable policy validation
`crpx_tree_validate()` checks connectivity alone. The general tree can legitimately retain missing names, duplicate names and finite signed lengths. `crpx_tree_validate_policy()` adds explicitly selected requirements:
| Flag | Requirement |
| --- | --- |
| `CRPX_TREE_REQUIRE_NAMED_TIPS` | Nonempty tip names. |
| `CRPX_TREE_REQUIRE_UNIQUE_TIPS` | Nonempty, exactly unique tip names; case-sensitive. |
| `CRPX_TREE_REQUIRE_COMPLETE_LENGTHS` | Every inspected edge has a length. |
| `CRPX_TREE_REQUIRE_NONNEGATIVE_LENGTHS` | Present lengths are nonnegative. |
| `CRPX_TREE_REQUIRE_POSITIVE_LENGTHS` | Present lengths are strictly positive. |
| `CRPX_TREE_INCLUDE_ROOT_STEM` | Apply the selected length checks to the root stem too. |
| `CRPX_TREE_REQUIRE_BINARY` | Every internal vertex has exactly two children; a one-tip tree is allowed. |

Sign checks and completeness are independent. For example, combine complete/nonnegative lengths for discrete likelihood, or complete/positive lengths for a Brownian calculation that requires positive edges. The root stem is excluded by default. Missing required lengths return `NOT_FOUND`; naming/sign violations return `INVALID_ARGUMENT`; unsupported degrees under the binary policy return `UNSUPPORTED`. Structural errors retain the ordinary structural-validation status.

Likelihood and data-mapping entry points continue to validate their own modelling requirements. An application can use the shared policy validator to report unsuitable input earlier, without maintaining its own tree walk. General storage and operation-specific validation serve different purposes: accepting signed lengths for distance-based estimation does not make a probabilistic transition valid at a negative evolutionary time.
## Checks
`test_tree_extensions` roots on every edge at endpoint/interior fractions for several binary and polytomous topologies, and checks all taxon-pair distances, original labels/IDs, stem transfer, revision invalidation, missing/signed/zero lengths and unsupported unary roots. Reserved 1,000/2,000/4,000-tip star builds check ordered labelled/length-annotated tail appends and external-topology refresh. Allocation counts prove no validator/tail refresh occurs between annotated appends; only occasional geometric tail growth allocates. `test_tree_allocation` intercepts allocators only in its local implementation compilation and injects every allocation failure in the rooting transaction, checking unchanged Newick, storage, revision and output; failed builder growth can be retried.
