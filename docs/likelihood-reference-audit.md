# Likelihood reference audit
The local biomcmc-recomb and guenomu implementations both calculate an HKY transition averaged over an exponential branch-length distribution, then use pruning and weighted site-pattern log likelihoods. Guenomu adds scaling and equal-weight site categories. Curupixa uses these mathematical ideas with explicit immutable models, separate worker-owned workspaces and recoverable errors.
## Sources and implemented assumptions
| Local reference | Relevant code |
| --- | --- |
| `../biomcmc-recomb/src/likelihood.c:46–80` | HKY normalisation, analytic decay rates and exponential-averaged transitions. |
| `../biomcmc-recomb/src/likelihood.c:94–134` | Pruning and weighted site likelihoods using the same averaged matrix on each edge. |
| `../biomcmc-recomb/src/update_chain_params.c:156–179` | Proposed segment mean branch length updates the transition matrix; individual edge lengths are absent from this calculation. |
| `../guenomu/lib/phylogeny.c:285–304` | Empirical frequencies, discrete gamma category means and initial averaged matrices. |
| `../guenomu/lib/phylogeny.c:332–378` | Closed-form HKY projectors, decay rates and exponential averaging. |
| `../guenomu/lib/likelihood.c:64–134` | Per-pattern/category scaling, root-edge evaluation and log-space category averaging. |
| `../guenomu/src/main.c:168` | The importance-sampling path selects the constant likelihood implementation. |
| `../guenomu/src/test/main_likelihood.c:31–47` | The standalone test invokes the actual likelihood engine. |

The model in [Martins, Leal and Kishino (2008), DOI 10.1371/journal.pone.0002651](https://pmc.ncbi.nlm.nih.gov/articles/PMC2440540/) explicitly assumes independent branch lengths across edges and alignment positions, conditional on the segment's mean. This supports the sitewise interpretation of the executable reference code. Guenomu's category mixture makes the mean vary between site categories; its matrices still integrate exponential branches within each category.
## Transition model and units
Use A,C,G,T state order, stationary frequencies $\pi$, $\pi_R=\pi_A+\pi_G$ and $\pi_Y=\pi_C+\pi_T$. For HKY, the unnormalised off-diagonal generator entries are $q_{ij}=\pi_j\kappa$ for A↔G or C↔T, and $q_{ij}=\pi_j$ otherwise. Divide by
$$
D=2\left[\kappa(\pi_A\pi_G+\pi_C\pi_T)+\pi_R\pi_Y\right]
$$
and set each diagonal to the negative row sum. This gives stationary mean substitution rate one, so fixed lengths and exponential means use expected substitutions per site. The positive analytic decay rates are $0$, $1/D$, $(\kappa\pi_R+\pi_Y)/D$ and $(\kappa\pi_Y+\pi_R)/D$.
Writing the associated projectors as $\Pi_k$, the fixed and averaged transitions are
$$
P(t)=\sum_k\Pi_k e^{-a_k t},\qquad
\bar P(\mu)=\int_0^\infty P(t)\mu^{-1}e^{-t/\mu}\,dt
=\sum_k\frac{\Pi_k}{1+a_k\mu}=(I-\mu Q)^{-1}.
$$
Here $\mu$ is an exponential **mean**, rather than its reciprocal rate. The zero-mean boundary is the identity. The reference documentation has an erroneous positive sign in the exponential density; the executable denominator uses the proper negative-density integral.
For arbitrary-state F81, $Q_{ij}=\pi_j/(1-\sum_s\pi_s^2)$ for $i\ne j$. Uniform frequencies give the equal-rates/Poisson model; with $K>1$ states its decay rate is $K/(K-1)$. The binary equal-rates case therefore has $P_{00}(t)=(1+e^{-2t})/2$. Curupixa handles the one-state constant process separately.
## Sitewise and shared-edge marginal likelihoods
Curupixa's `CRPX_BRANCH_SITEWISE_EXPONENTIAL` matches the reference's independent-edge, independent-site integration. Pattern compression remains valid: each pattern's averaged likelihood can be raised to its integer weight. Each supplied edge may have a distinct exponential mean.
An alignment-wide shared-edge integral defines a different model:
$$
\prod_h\int L_h(\mathbf t)p(\mathbf t)\,d\mathbf t
\quad\hbox{versus}\quad
\int\prod_h L_h(\mathbf t)p(\mathbf t)\,d\mathbf t.
$$
The first expression is the implemented sitewise model. The second uses one random length per edge shared by all positions and requires a separate integration method; averaged edge matrices alone do not compute it.
For a binary process with a stationary root and two observed zero states joined by one edge, one site's joint likelihood is $f(T)=(1+e^{-2T})/4$. With two repeated sites and $T\sim\mathrm{Exp}(\text{mean }\mu)$,
$$
\mathbb E[f(T)^2]=\frac{1+2/(1+2\mu)+1/(1+4\mu)}{16},\qquad
\mathbb E[f(T)]^2=\left[\frac{1+1/(1+2\mu)}{4}\right]^2.
$$
At $\mu=0.5$, these equal $7/48$ and $9/64$ respectively. The distinction concerns the model, rather than a numerical approximation.
## Root-edge convention
The references treat an artificial binary root as superfluous. Guenomu directly connects its two root children with one transition matrix. A fixed process satisfies $P(t_1)P(t_2)=P(t_1+t_2)$, so splitting an edge preserves its path likelihood. Independently averaged exponential edges generally satisfy $\bar P(\mu_1)\bar P(\mu_2)\ne\bar P(\mu_1+\mu_2)$.
The default curupixa convention attaches one process to each supplied non-root edge. The explicit binary-root collapse option instead uses one specified effective root-edge parameter and an identity process on the other half. This preserves the reference's single root-edge interpretation without silently inferring an exponential prior from two supplied half-edges. A root at a degree-three internal node avoids artificial edge splitting naturally. The root stem never contributes an additional transition.
## Implementation hazards avoided
- Biomcmc-recomb stores loop indices and partial-vector pointers at file scope (`likelihood.c:29–32`), so concurrent evaluations share mutable scratch state. Its pruning lacks scaling and can underflow on large trees.
- Guenomu allocates `lnmax` arrays with `malloc` (`phylogeny.c:200–203`), then initialises tip scales only for category indices greater than zero (`phylogeny.c:49–55`). Category-zero tip scales are read uninitialised by likelihood evaluation.
- Guenomu exposes a nominal state-count argument but requires four states (`phylogeny.h:60`). Generic state handling requires new loops and storage contracts.
- Zero purine or pyrimidine mass causes division by zero in the legacy analytic projectors. Curupixa requires positive finite frequency weights and checks representability after normalisation.
- Global real/dummy likelihood function pointers and cyclic proposal-vector ownership belong to the original applications. Curupixa keeps scoring mode, data ownership and cache validity in explicit objects.

Both reference sources credit Leonardo de Oliveira Martins and carry GPL-2.0-or-later notices (2006 and 2009 respectively). Derived curupixa sources retain that attribution and use the project's GPL-3.0-or-later designation.
