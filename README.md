<img src="recipe/curupixa.png" height="220" align="right" alt="Curupixá">
# Curupixá
Curupixa provides low-level C objects and algorithms for phylogenetic software. It is a modern reimplementation of [biomcmc-lib](https://github.com/quadram-institute-bioscience/biomcmc-lib), with explicit state and recoverable errors for use from C, C++, language bindings, MPI processes and worker threads. The [migration inventory](docs/migration.md) records implemented capabilities, API changes and remaining gaps against each biomcmc-lib header.

**Author:** Leonardo de Oliveira Martins, Quadram Institute Bioscience, Norwich, UK.

**Etymology:** Curupixá is a Brazilian tree, *Micropholis venulosa*.
## Developer build
CMake is the primary build system, using version 3.16 or later and a C17 compiler. The current platform baseline is a POSIX system with GCC or Clang. The core needs the C runtime and maths library. It downloads no dependencies and does not require OpenMP, compression libraries, Python or libcheck. Both CMake and Autotools pass local build, test and installation checks, including XZ; CMake fits recursive developer clones and submodule embedding more directly. The [build-system comparison](docs/build-systems.md) records the decision and the role retained for Autotools.

Run these commands from the checkout:
```sh
cmake -S . -B build/cmake -DCMAKE_BUILD_TYPE=Release
cmake --build build/cmake --parallel
(cd build/cmake && ctest --output-on-failure)
```
Generated files stay inside the build directory. Separate build directories allow different compilers, build types and feature selections without modifying the source checkout. Use a CMake toolchain file for cross-compilation; configuration does not probe the build machine's instruction set.

The default library is static and uses position-independent code, so a downstream Python extension can link it. To build a shared library for dynamic linking or `ctypes`/CFFI:
```sh
cmake -S . -B build/shared -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=ON
cmake --build build/shared --parallel
```
The C API uses `extern "C"` guards for C++ callers. New object APIs expose explicit ownership and status returns; a binding should translate these into its own exceptions or result values. There is no Python binding in this repository yet. The [validation report](docs/validation.md) records the tested build variants, installed consumers, Python smoke check and platform/tooling limits.
## Embed as a submodule
A downstream project can obtain the repository and its submodules with SSH:
```sh
git clone --recurse-submodules git@github.com:OWNER/DOWNSTREAM.git
```
To add curupixa to a project you maintain:
```sh
git submodule add git@github.com:leomrtns/curupixa.git external/curupixa
```
A downstream CMake project links the same target used by the standalone build:
```cmake
add_subdirectory(external/curupixa)
add_executable(my_phylogeny main.c)
target_link_libraries(my_phylogeny PRIVATE curupixa::curupixa)
```
Tests and installation rules default to off when curupixa is a subdirectory. No separate installation, bootstrap or network request is needed during the build. The parent controls `BUILD_SHARED_LIBS`, compiler flags and toolchain selection.

Include the umbrella header with `#include <curupixa.h>`. Installed packages also support `#include <curupixa/curupixa.h>`.
## Install into a development environment
Installation is optional. For conda, container images or another managed prefix, set the destination explicitly:
```sh
cmake -S . -B build/package -DCMAKE_BUILD_TYPE=Release \
  -DCURUPIXA_BUILD_TESTS=OFF -DCMAKE_INSTALL_PREFIX="$CONDA_PREFIX"
cmake --build build/package --parallel
cmake --install build/package
```
For another prefix, replace `$CONDA_PREFIX` with the intended directory. Packaging systems can stage installation with `DESTDIR`; no interactive steps or system-wide commands are required. CMake honours the normal compiler environment and `CMAKE_INSTALL_LIBDIR`, including distribution choices such as `lib64`.

An installed CMake consumer uses:
```cmake
find_package(curupixa CONFIG REQUIRED)
target_link_libraries(my_phylogeny PRIVATE curupixa::curupixa)
```
Point `CMAKE_PREFIX_PATH` at the prefix when it is outside the toolchain's search path. Other build systems can use `pkg-config --cflags --libs curupixa`; add `--static` when linking the static archive so that private dependencies are included. A package installs headers below `include/curupixa`, the library in the configured library directory, and CMake and pkg-config metadata alongside it.

A conda recipe or a Galaxy dependency package can use the configure, build and install steps above with its supplied prefix. MPI callers link their own MPI implementation; curupixa does not initialise MPI or add an MPI dependency.
## Build options
| CMake option | Default | Purpose |
| --- | --- | --- |
| `BUILD_SHARED_LIBS` | `OFF` | Build a shared library instead of a static archive. |
| `CURUPIXA_BUILD_TESTS` | `ON` standalone; `OFF` embedded | Build dependency-free regression tests. |
| `CURUPIXA_INSTALL` | `ON` standalone; `OFF` embedded | Generate installation and package metadata rules. |
| `CURUPIXA_ENABLE_OPENMP` | `OFF` | Enable optional OpenMP compilation. |
| `CURUPIXA_ENABLE_SSE42` | `OFF` | Require SSE4.2 on target CPUs. |
| `CURUPIXA_ENABLE_AVX2` | `OFF` | Require AVX2 and FMA on target CPUs. |
| `CURUPIXA_ENABLE_RDRAND` | `OFF` | Require Intel RDRAND on target CPUs. |
| `CURUPIXA_ENABLE_ZLIB` | `OFF` | Read and write gzip file streams. |
| `CURUPIXA_ENABLE_BZIP2` | `OFF` | Read and write bzip2 file streams. |
| `CURUPIXA_ENABLE_LZMA` | `OFF` | Read and write XZ file streams. |

Leave hardware options disabled for general conda packages, heterogeneous clusters and portable containers. They add target compiler flags to the entire library and require compatible deployment CPUs. They do not provide runtime dispatch. Context RNG draws are serialised across host threads with or without OpenMP. Use a separate explicit random stream per worker for reproducible parallel calculations; sequence assignment from a shared context depends on scheduling.
## Existing Autotools consumers
Autotools support remains available. A developer checkout needs Autoconf, Automake and Libtool for bootstrapping; the CMake path avoids these tools. Run from the checkout:
```sh
./autogen.sh
mkdir -p build/autotools
cd build/autotools
../../configure --prefix="$CONDA_PREFIX"
make
make check
make install
```
The installed `libcurupixa` now includes all implemented modules. Existing downstream Autotools projects can continue linking `lib/libcurupixa_static.la`, the uninstalled convenience library. Regenerate the configuration after changing `configure.ac` or `Makefile.am`; old generated `configure` and makefiles do not contain the updated build definitions.

OpenMP, SSE4.2, AVX2 and RDRAND are opt-in through `--enable-openmp`, `--enable-sse`, `--enable-avx2` and `--enable-rdrnd`. Compression uses `--enable-zlib`, `--enable-bzip2` and `--enable-lzma`. Compiler optimisation and instrumentation remain under the caller's control through `CFLAGS` and `LDFLAGS`. The modern regression tests do not require libcheck; older exploratory `check_*` and Dieharder programs remain in `tests/` for reference.
## Migration choices
This update changes the API and ABI. Opaque objects with explicit ownership replace exposed reference-counted structures; packed symmetric distances replace overloaded matrix triangles. The SVD backend uses dense Jacobi workspaces and materialises sparse input. Mapped split assignment gives unmatched edges a positive penalty; the dSPR result separately exposes the historical zero-penalty h-distance convention. Corrected BioNJ and probability calculations can intentionally change results; the [tree design](docs/tree-design.md) records BioNJ and SPR branch-length policies. Greedy dSPR and all-root duplication/loss/extra-lineage reconciliation now have separate analysis objects. Exact SPR, specialised move proposals and branch fitting remain gaps. The new likelihood/parsimony interface is described below; see the [migration inventory](docs/migration.md) and [function-level map](docs/api-migration.csv).
## Cladistic distances and reconciliation
The `spr_distance`, `ancestry_index`, `reconciliation_roots` and `taxon_map` headers adapt reusable algorithms from guenomu and genefam-dist. They provide greedy dSPR, RF and hdist, explicit gene-copy/species mapping, a compact reusable LCA index and duplication/loss/extra-lineage scores over every gene-tree edge rooting. All scratch and results stay outside the tree; workers can share immutable species indices and input trees. The [distance/reconciliation guide](docs/cladistic-distances.md) and [example](examples/gene_species_distances.c) specify assumptions, complexity and ownership.

Important conventions are explicit: dSPR is a heuristic, equal-cost choices can change it, and its hdist retains the reference's zero dummy penalty. The general splitset assignment API uses a positive unmatched penalty. Reconciliation supports full or represented-species losses, and returns separate optimal-root indices for each objective.
## Header structure and memory
Include individual headers such as `<curupixa/tree.h>` or `<curupixa/spr_distance.h>` as needed. Include guards support repeated inclusion, and the public dependency graph is acyclic. A separate `context_types.h` keeps hashing/RNG declarations from importing the full runtime. The umbrella remains available for convenience.

Adding algorithms increases available library code; per-tree and numerical buffers are allocated only when their objects/functions are used. Static linking selects the required object modules. Shared libraries map code/data pages, with resident memory depending on the pages used; loading does not construct analysis workspaces. The [footprint guide](docs/library-footprint.md) gives measured sizes, link granularity and the calculations that require substantial working memory.
## API development
Keep mutable state in an object or caller-owned workspace. An MPI rank should own its contexts and workspaces; worker threads should use the ownership rules documented in each header. New functions report failure to the caller and do not terminate the process. Topology, labels, branch lengths and algorithm caches can have different lifetimes; the tree interface keeps optional data outside its compact connectivity backbone.

See [runtime and hashing ownership](docs/runtime-and-hashing.md) for stream reproducibility, synchronisation and hashing interfaces. Statistical functions expose typed rankings, checked histograms and quantiles; [probability functions](lib/prob_distribution.h) cover normal, lognormal, gamma, Poisson and truncated discrete exponential distributions with caller-owned random streams.

## Discrete-state likelihood and parsimony
The new `character_data`, `substitution_model`, `likelihood` and `character_analysis` objects provide fixed-length likelihood, analytically sitewise-marginalised likelihood and parsimony through one data/workspace interface. Binary and arbitrary finite state counts are supported, including alphabets above 64 states. Poisson/equal-rates and general-frequency F81 work for any state count; HKY uses A,C,G,T. All models have closed-form transitions and require no numerical eigensystem. GTR is left for a separate optional backend so these simple paths retain their linear state-count cost.

Each fixed edge can use its own stored length or an explicit override. In the branchless mode, each edge can use a supplied exponential prior mean. Both reference programs integrate branch lengths independently across edges **and sites**; this differs from a whole-alignment integral with one random length shared across all positions. An explicit binary-root option reproduces their single root-edge prior without silently assigning two priors to that edge. The [reference audit](docs/likelihood-reference-audit.md) derives these assumptions and records corrections to the original numerical handling.

Worker-owned workspaces expose per-pattern log likelihoods and downward partials, with log-space pruning and an optional streaming mode. The same immutable data and analysis wrapper support equal-cost or general directed-cost Sankoff parsimony. The [developer guide](docs/character-analysis.md) and [binary example](examples/discrete_characters.c) show all three scoring modes and ownership. The [likelihood design](docs/likelihood-design.md) records implementation costs and remaining extensions. No libPLL dependency is needed.
## Build references
The build follows CMake's [library target conventions](https://cmake.org/cmake/help/latest/command/add_library.html), [installation and export support](https://cmake.org/cmake/help/latest/command/install.html) and [OpenMP imported targets](https://cmake.org/cmake/help/latest/module/FindOpenMP.html).
## Licence
SPDX-License-Identifier: [GPL-3.0-or-later](LICENSE).

Copyright (C) 2022–today Leonardo de Oliveira Martins. See [references](docs/references.md) for incorporated implementations and attribution.

The [logo owl](https://www.flickr.com/photos/biodivlibrary/51267563412/in/album-72157719506234015/) derives from *Fauna japonica*, volume 4, made available by the Biodiversity Heritage Library in the public domain.
