# Library organisation
Curupixa builds three libraries from one repository. Each implementation is compiled once. Higher layers link the lower layer; the lower layer never calls a higher one. The CMake target `curupixa::curupixa` and pkg-config package `curupixa` select all enabled components without producing a fourth compiled library.
| Component | CMake target | Library | Scope |
| --- | --- | --- | --- |
| General phylogenetics | `curupixa::phylogenetics` | `libcurupixa_phylogenetics` | Trees, numerical primitives, sequence data, likelihood, PCM, ASR, parsimony and shared comparison/reconciliation kernels. |
| Gene/species methods | `curupixa::gene_species` | `libcurupixa_gene_species` | Explicit copy-to-species correspondence, LCA reconciliation and paralogy-aware comparison interfaces. |
| Inference | `curupixa::inference` | `libcurupixa_inference` | Conditional-clade probability models, sampling, probability evaluation and caller-directed optimisation. |
Inference links gene/species methods, which link general phylogenetics. A general phylogenetics client therefore has no gene/species or inference link dependency. Public headers remain in `include/curupixa`; include an individual header or the component umbrella `phylogenetics.h`, `gene_species.h` or `inference.h`. `curupixa.h` includes all three umbrellas. Headers declare interfaces; including them allocates no analysis state.
## Trees and biological interpretation
A tree records connectivity, optional names and optional lengths. Gene copies are ordinary distinct tips, so likelihood, parsimony, ancestral-state reconstruction and phylogenetic comparative methods can be applied to gene trees directly. An analysis requiring species membership receives that correspondence separately. No function attempts to infer whether a tree represents paralogs by inspecting its labels.
The storage root provides an orientation for traversal. Unrooted RF, split comparisons and quartet comparisons ignore that orientation. Rooted reconciliation and conditional-clade models assign meaning to the supplied root. Tanglegrams retain the supplied rooted subtrees when changing child order. General tree storage accepts multifurcations, unary vertices and absent lengths; each algorithm documents the narrower conditions required by its mathematics. DTL and the initial conditional-clade model require binary inputs. Likelihood requires branch parameters, whereas the topology comparisons ignore lengths.
## Shared distances without duplicate implementations
`splitset_distances.c` and `spr_distance.c` belong to general phylogenetics. `gene_species_distances.c` prepares comparable splitsets and calls those implementations. `crpx_gene_species_splitsets()` creates one pair of owned splitsets with a common gene-copy bit order. Repeated species leaf IDs expand a species into its gene copies; absent species add no bits. Callers can retain the splitsets and use `crpx_splitset_compare()` or `crpx_splitset_spr_distance()` repeatedly without remapping the original trees. Rebuild the affected splitset after topology or correspondence changes.
The ordinary entry points never prepare a species map. The mapped entry points add only the preparation needed by their input semantics; they do not copy the distance algorithm. A small function call between libraries does not repeat preparation. Dynamic libraries are linked directly, without a plugin dispatcher or runtime search for implementations.
RF and quartet measures with optimised correspondence between repeated labels are separate mathematical problems. Expanding a specified map does not claim to solve that correspondence problem. Likewise, the paralogy-aware quartet objective has its own counting convention and cannot be replaced by ordinary quartet distance on relabelled tips.
## Existing functions by component
| Functions or headers | Component | Shared dependencies |
| --- | --- | --- |
| `tree.h`, `ancestry_index.h`, `topology_randomise.h`, `newick_space.h` | General phylogenetics | Status codes; topology workspaces and numerical primitives where needed. |
| `topology_distance.h`, `splitset_distances.h`, `spr_distance.h`, `quartet_distance.h`, `triplet_distance.h`, `tree_comparison.h`, `tanglegram.h` | General phylogenetics | Tree connectivity, split/bit operations, assignment or ancestry primitives. |
| `alignment.h`, `character_data.h`, `substitution_model.h`, `likelihood.h`, `character_analysis.h` | General phylogenetics | Data storage, mappings from observations to tips and caller-owned numerical workspaces. |
| `continuous_traits.h`, `reconstruction.h`, `parsimony.h` | General phylogenetics | Trees and the existing trait or discrete-state calculations. |
| `gaussian.h`, `phylogenetic_regression.h`, `trait_association.h` | General phylogenetics | Covariance factors and Gaussian messages; tree preparation is separate from regression. |
| `sampling.h`, `quasi_sequence.h`, `stochastic_mapping.h` | General phylogenetics | Explicit random callbacks, tree builders and existing substitution generators. |
| `clock.h`, `particle_weights.h` | General phylogenetics | Tree/date conversion, existing random callbacks and allocation-free numerical weights. |
| `molecular_clock.h` | Inference | Tree traversals and the existing lognormal probability functions; no particle population. |
| `tree_prior.h` | Inference | General tree constructors, draws and comparison indices. |
| `dtl.h` | General phylogenetics | Explicit correspondence and prepared host-tree states; no inference dependency. |
| `reconciliation.h`, `reconciliation_roots.h` | Gene/species methods | Trees and ancestry indices from general phylogenetics. |
| `taxon_map.h`, `gene_species_distances.h`, `paralog_quartets.h` | Gene/species methods | Exact label matching or explicit maps; shared split and dSPR functions. |
| `amalgamation.h`, `tree_search.h`, `inference.h` | Inference | General trees, caller-owned state and the lower-level reconciliation kernel. |
The [comparative methods guide](comparative-methods.md) describes the shared scoring and reconstruction primitives, regression and association. The [sampling guide](tree-sampling.md) describes prior conditioning and the lower-level constructors. The [comparison and inference guide](comparison-and-inference.md) describes the new function contracts and explicitly records the reference implementations and remaining algorithms.
Remaining numerical, hashing, random-stream, compression, distance-matrix and clustering modules belong to general phylogenetics. Numerical minimisation in `minimiser_simplex.h` remains a general primitive. Species-tree search is a higher-level client of an explicit objective.
## Repeated evaluations and ownership
Prepare immutable data separately from mutable worker storage. Different workers may share a prepared ancestry index, quartet index or conditional-clade model; they need independent mutable workspaces. No calculation in this branch starts worker threads. OpenMP remains an optional compilation setting; concurrency is controlled by the downstream program.
`crpx_reconcile_all_roots_into()` retains scratch and score capacities. Its result is a workspace-owned view, invalidated by the next evaluation, clear or destruction. `crpx_reconciliation_root_mapping()` extracts a mapping for a selected scored root without rerooting the input. `crpx_reconcile_all_roots()` remains a one-off convenience function using the same implementation and returning an owned result. Structural validation still allocates temporary arrays, so retained workspace capacity is not a blanket guarantee of zero allocations.
The prepared likelihood interface already retains tip mapping and traversal. It recomputes numerical kernels and all downward partials on evaluation. Incremental partial updates and proposal rollback are not implemented here. Tree-search callbacks can retain their own workspaces, but must refresh prepared topology snapshots when the trial topology changes.
## Building components
```sh
cmake -S . -B build/shared -DBUILD_SHARED_LIBS=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build/shared --parallel
ctest --test-dir build/shared --output-on-failure
```
For general phylogenetics alone, add `-DCURUPIXA_BUILD_GENE_SPECIES=OFF -DCURUPIXA_BUILD_INFERENCE=OFF`. Inference requires gene/species methods in this initial build hierarchy. CMake propagates link dependencies automatically:
```cmake
find_package(curupixa CONFIG REQUIRED)
target_link_libraries(my_program PRIVATE curupixa::phylogenetics)
```
Use `curupixa::gene_species`, `curupixa::inference` or `curupixa::curupixa` when those interfaces are needed. The pkg-config packages are `curupixa-phylogenetics`, `curupixa-gene-species`, `curupixa-inference` and the umbrella `curupixa`. There is no longer a monolithic `-lcurupixa`; raw linker users must select component libraries in dependency order or obtain the flags from pkg-config. No compatibility constraint from the experimental previous interface is imposed on the split.

General graph, assignment and clustering interfaces also belong to `curupixa_phylogenetics`: see [their kernels and dependencies](graph-clustering.md). The [sequence-hashing guide](sequence-hashing.md) covers the additional split-rotation iterator and its shared integer mixer.
