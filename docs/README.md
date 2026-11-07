# Curupixá
Curupixa is a low-level C library for phylogenetic software, with explicit object ownership, recoverable errors and caller-owned state for MPI and threaded applications. It is being developed as the successor to biomcmc-lib. The [migration inventory](migration.md) records implemented capabilities, deliberate changes and remaining algorithms.
## Developer build
Use CMake 3.16 or later and a C17 GCC/Clang toolchain on a POSIX system. From the checkout root:
```sh
cmake -S . -B build/cmake -DCMAKE_BUILD_TYPE=Release
cmake --build build/cmake --parallel
(cd build/cmake && ctest --output-on-failure)
```
The default static library uses position-independent code and needs no optional dependencies. Add `-DBUILD_SHARED_LIBS=ON` for a shared library. Compression codecs and CPU-specific instructions are optional and disabled by default.
## Downstream use
Clone downstream projects with `git clone --recurse-submodules git@github.com:OWNER/DOWNSTREAM.git`. A parent CMake project can embed curupixa directly:
```cmake
add_subdirectory(external/curupixa)
target_link_libraries(my_phylogeny PRIVATE curupixa::curupixa)
```
Embedding needs no separate installation or bootstrap. Alternatively, configure with the package manager's `CMAKE_INSTALL_PREFIX`, build, then run `cmake --install` for conda, Galaxy dependency packages or container images. Installed CMake consumers use `find_package(curupixa CONFIG REQUIRED)`; other build systems can use pkg-config.
## API guides
- [Source and installed files](general_structure.md) — Repository layout, public headers and installation paths.
- [Migration coverage](migration.md) — Module coverage and API changes from biomcmc-lib.
- [Tree design](tree-design.md) — Connectivity, ownership, workspaces and rearrangements.
- [Downstream trait API improvements](downstream-trait-api.md) — Prepared analyses and the changes prompted by branchtraits.
- [Tree rooting, builders and validation](tree-extensions.md) — Stable-ID rooting, bulk construction and operation-specific checks.
- [Discrete ancestral reconstruction](reconstruction.md) — Full marginal ASR and optimal parsimony assignments and state sets.
- [Continuous traits](continuous-traits.md) — Brownian fitting, ancestral reconstruction and independent contrasts.
- [Cladistic distances and reconciliation](cladistic-distances.md) — RF, hdist, approximate SPR and gene/species reconciliation.
- [Headers and memory footprint](library-footprint.md) — Header dependencies, linking and analysis storage costs.
- [Runtime and hashing](runtime-and-hashing.md) — Contexts, random streams, sequence indexes and threading.
- [Likelihood and parsimony](character-analysis.md) — Discrete-state data, models and scoring interfaces.
- [Likelihood design](likelihood-design.md) — Pruning costs, partial buffers and backend boundaries.
- [Likelihood reference audit](likelihood-reference-audit.md) — Branch-integration assumptions in the reference implementations.
- [Local validation](validation.md) — Numerical tests, package checks and tested platforms.
- [Build-system decision](build-systems.md) — CMake and Autotools installation choices.
- [Source references and attribution](references.md) — Scientific references and implementation provenance.
## Licence
Curupixa is licensed under GPL-3.0-or-later.
