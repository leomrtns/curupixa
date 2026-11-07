# Discrete ancestral reconstruction
Include `reconstruction.h` for full-data marginal likelihood reconstructions and optimal parsimony reconstructions. These functions move the upward/outside contributions and ambiguity-set calculations from branchtraits into curupixa. Score-only callers continue using `likelihood.h` or `character_analysis.h`; their calls allocate no reconstruction arrays or run upward/traceback passes.
## Ownership and basic use
Construct a `crpx_reconstruction_workspace` for each concurrent worker. Construction allocates only a small control object. Each call returns an independently owned opaque result, with its own destructor. The result remains valid after another call, after destroying the workspace, or after destroying the input tree/data/model. Result node indices refer to the input tree at the time of calculation; labels are not copied into results. Rebuild a result after changing topology, observations or model parameters.

Workspaces retain capacities but no borrowed tree/data/model identities between calls. `crpx_reconstruction_workspace_clear()` releases all retained buffers and nested workspaces without invalidating owned results; the next call allocates lazily again. All immutable inputs may be shared read-only by workers using separate mutable workspaces. Calls set their output pointer to `NULL` on failure, and leave previously returned results intact. As with other curupixa functions, the caller synchronises mutations and destruction of shared inputs.
```c
crpx_reconstruction_workspace *worker = NULL;
crpx_likelihood_reconstruction *result = NULL;
crpx_status status = crpx_reconstruction_workspace_new(&worker);
if (status == CRPX_OK)
  status = crpx_likelihood_reconstruct(tree, data, model, &options, worker, &result);
if (status == CRPX_OK) {
  const double *probability = crpx_likelihood_reconstruction_probabilities(result, node, pattern);
  /* probability[state] uses data/model state ordering. */
}
crpx_likelihood_reconstruction_free(result);
crpx_reconstruction_workspace_free(worker);
```
The [character-analysis guide](character-analysis.md) describes constructing immutable data/models and choosing branch options. Sites refer to compact patterns through `crpx_character_data_site_patterns()`. Reconstruction results use pattern indices, avoiding repeated ancestral vectors for identical columns.
## Marginal likelihood ASR
`crpx_likelihood_reconstruct()` returns $P(X_v=s\mid\text{all tip observations})$ for every node, state and pattern. It combines downward subtree likelihoods with an upward message containing the root prior and observations outside that subtree. Downward partials alone provide a different conditional quantity at non-root nodes. The stationary root prior, finite-state model, ambiguity weights, fixed versus sitewise-exponential edges and artificial binary-root collapse have exactly the same interpretation as `crpx_likelihood_evaluate()`. Root stems remain excluded. The function ignores the retention flag because reconstruction requires downward partials.

All currently provided models are reversible. The upward transition uses the transpose operator through detailed balance, preserving linear work per state for Poisson/F81. HKY uses its four-state analytic log kernel. Pruning and reconstruction share private transition-kernel code: adding a backend or improving its numerical implementation requires no corresponding downstream implementation of an upward transition. The public contract specifies probabilities and conditioning rather than exposing traversal order or numerical buffer layout.

Messages remain in log space, and sibling exclusion counts impossible contributions explicitly rather than subtracting negative infinity. Zero-length edges, unary nodes and arbitrary polytomies are supported. An impossible pattern has no conditional distribution, so reconstruction returns `INVALID_ARGUMENT`, including for a zero-weight pattern. Scoring alone continues to allow impossible data with log likelihood `-INFINITY`. Marginal vectors do not identify a single joint maximum-posterior history; choosing each node's most probable state independently need not maximise a joint probability. Joint MAP histories, stochastic histories and transition-event counts are separate future APIs.

The result score equals the likelihood scorer's weighted log likelihood. Pattern weights affect that score; they do not change a pattern's conditional node probabilities. Full output and retained likelihood partials require $O(PVK)$ memory, with $O(VK)$ extra worker scratch for $P$ patterns, $V$ nodes and $K$ states. Poisson/F81 passes use $O(PVK)$ time; no dense arbitrary-state transition matrices are allocated. A single-state model and alphabets above 64 states are supported.
## Parsimony reconstruction
`crpx_parsimony_reconstruct()` computes the hard-polytomy Sankoff objective. Select `CRPX_PARSIMONY_UNORDERED` for unit change costs, `CRPX_PARSIMONY_ORDERED` for $|i-j|$, or `CRPX_PARSIMONY_CUSTOM` with an explicit row-major directed cost matrix. The first two use linear min-plus transforms, so ordered reconstruction needs no dense $K\times K$ matrix. Custom finite nonnegative costs require zero diagonal, allow asymmetry, and make rooting relevant. Branch lengths have no role in this objective. Positive observation entries identify allowed states; their magnitudes do not change parsimony costs.

The result provides two different forms of reconstruction:
- `assignment()` gives one jointly coherent optimal assignment, choosing the smallest state on exact local ties. Its selected edge costs sum to the reported pattern score.
- `optimal_states()` gives, for each node, every state that occurs in some globally optimal assignment for that pattern. This uses downward plus outside costs. Choosing each node independently from these sets need not produce a jointly optimal assignment.

Each pattern minimises its own unweighted objective, even when its weight is zero; weights affect the total score only. Globally admissible here refers to that full-tree pattern objective, rather than the local subtree objective. Integer unit/ordered scores are exact up to the checked $2^{53}$ weighted-total limit. With fractional custom costs, accumulation and sibling subtraction can introduce floating-point round-off. The explicit relative tolerance controls only the ambiguity-set comparison:
$$
|c_s-c_{\min}|\leq\epsilon\max(1,|c_s|,|c_{\min}|).
$$
Choose zero for exact representable integer costs, or a small tolerance appropriate to the number and scale of fractional operations. A larger tolerance can deliberately admit near-optimal states; it does not alter the optimum or coherent traceback. Output storage costs $O(PVK)$ for admissibility bytes and $O(PV)$ for the selected states, plus $O(VK)$ worker scratch. Ordered/unordered time costs $O(PVK)$; arbitrary directed costs require $O(PEK^2)$ for $E$ edges.
## Validation and downstream stability
The reconstruction tests enumerate complete state assignments on tiny trees independently of pruning/messages. They compare every node probability for Poisson, F81 and HKY under fixed and sitewise-exponential branches, including collapsed roots, unary nodes, polytomies, unknown/weighted observations and zero edges. They separately enumerate every optimum under unordered, ordered and directed parsimony, checking all admissible states and the cost of the coherent assignment. Additional cases cover subnormal edges, impossible data, independent result lifetime, single-state models and 70-state characters.

Downstream applications can replace their private upward/parsimony passes with these functions while retaining application-specific fitting, output formatting and model choices. Curupixa can change its numerical implementation while preserving this mathematical/API contract; regression comparisons also check agreement with the public score-only interfaces. Scientific bug fixes may intentionally change erroneous outputs and should be documented as such.
