# Local validation
Validated on 6–7 November 2026 with GCC 13.3.0, CMake 3.28.3 and Python 3.12.8 on Linux. All generated builds, test fixtures and temporary consumers stayed inside this checkout's `build/` directory.
## Build and regression matrix
| Configuration | Result |
| --- | --- |
| Static Release, optional dependencies disabled (6 November) | 19/19 CTest suites passed; current source also builds without optional dependencies. |
| Static Release, zlib 1.3, bzip2 1.0.8 and liblzma 5.4.5 enabled | 19/19 suites passed. |
| Shared Release, all three codecs enabled, strict warnings (7 November) | 21/21 suites passed. |
| Debug, optional codecs disabled, `-Wall -Wextra -Wpedantic -Werror`, AddressSanitizer and UndefinedBehaviourSanitizer (7 November) | 21/21 suites passed. Earlier all-codec sanitised run passed the then-current 19 suites. |
| OpenMP 4.5, Release, all codecs enabled | 19/19 suites passed. |
| Autotools, all codecs, OpenMP and strict warnings (7 November) | 21/21 `make check` suites passed; bootstrap, install and XZ-enabled `make distcheck` passed. |
| CMake build from the Autotools release tarball, XZ enabled (7 November) | 21/21 suites passed. |

The suites cover containers, trees, runtime ownership/threading, statistics, probability functions, hashing, suffix indexing, topology operations/construction/collections, alignment patterns and distances, compression, GFF3, parsimony/reconciliation, minimisation, SVD, clustering, assignment/edit/split distances, analytic substitution models, generic character data and phylogenetic likelihoods. Scientific regressions include probability tails and quantile support, BioNJ branch distances, OPTICS border points, weighted parsimony overflow, mapped multi-copy splits and scaled matrix reconstruction. Parser tests cover malformed input, quotes/comments, deep trees, large polytomies, locale-independent decimal scores and compressed-stream integrity.
Cladistic regressions compare RF and assignment costs with independent implementations, test one-pair assignments, polytomies, repeated copies, word boundaries and workspace reuse. A legacy genefam-dist driver and curupixa agreed on every RF and non-degenerate h-distance for 160 seeded binary-tree pairs; documented one-pair corrections and heuristic tie differences appear in the [distance guide](cladistic-distances.md). All-root reconciliation tests rebuild every edge-rooted gene tree independently and compare duplication/full-loss scores with the existing rooted reconciler, and represented-loss/extra-lineage scores with independently pruned species trees. Cases include unequal species depths, absent species, repeated gene copies, one-tip inputs and unsupported polytomies. Ancestry snapshots match parent-walk LCA for all node pairs in the test tree and remain usable after destroying that tree.

New model tests compare analytic kernels against independently constructed generators, uniformisation exponentials and pivoted resolvents. Likelihood tests enumerate complete small-tree state assignments independently of pruning, check fixed and sitewise exponential modes, artificial-root semantics, partial vectors, and deep-tree log-space stability. Generic character/parsimony tests cover more than 64 states, custom costs, ambiguous observations, overflow and concurrent workers sharing immutable inputs.
XZ now passes the same codec tests as gzip/bzip2: concatenated members, empty streams, mixed line/block reads, truncation, trailing garbage and payload corruption. XZ-specific tests also cover valid four-byte padding between/after streams. Cross-API read/write checks with Python's gzip, bz2 and lzma modules passed for a 1,000,003-byte input.
The sanitised run used `ASAN_OPTIONS=detect_leaks=0`: LeakSanitizer cannot operate under the execution environment's ptrace wrapper. Address/undefined-behaviour checks ran; this run does not establish leak-check coverage.
## Consumer and package checks
- Fresh installed static/shared libraries linked and ran from separate C and C++ consumer projects. The new gene/species example compiled and ran as C17 and C++17 against both installed variants.
- All 49 installed headers compiled independently and repeatedly as C17 and C++17 with strict warnings; combined forward/reverse/shuffled orders also passed (104 compilation checks). Both source and installed header dependency graphs were acyclic. Private runtime, constant and substitution-model headers were absent from the installed API.
- The public discrete-character example compiled and ran against freshly installed static/shared CMake packages as both C17 and C++17, including fixed/sitewise likelihood and the same-data parsimony interface.
- A Python `ctypes` smoke test evaluated binary likelihood and parsimony and destroyed all objects. The fixed example returned log likelihood `-5.3920926373562` and parsimony score `3`.
- Installed OpenMP-built static/shared packages worked from C++-only CMake projects with OpenMP discovery disabled. An installed all-codec shared consumer also worked with codec package discovery disabled; codec development headers remain private build requirements.
- A parent project built and ran through `add_subdirectory()` and `curupixa::curupixa`, with curupixa tests and installation disabled by default.
- CMake and pkg-config consumers worked with an overridden installation prefix. pkg-config also worked after moving a prefix with a nested relative library directory, and with a separately configured absolute library directory.
- The default shared object measured 345,624 bytes on disk and 120 KiB resident immediately after loading in the measured process. Static tree-only clients had identical text/data/bss sizes with narrow and umbrella includes, and no likelihood/reconciliation/SPR/SVD/GFF3 symbols; see [memory measurements](library-footprint.md).
- `git diff --check` passed. New required CMake templates are visible to Git; generated build products remain ignored.
## Reproduce the default checks
From the checkout root:
```sh
cmake -S . -B build/cmake -DCMAKE_BUILD_TYPE=Release
cmake --build build/cmake --parallel
(cd build/cmake && ctest --output-on-failure)
```
The parenthesised CTest invocation also works with the declared CMake 3.16 minimum. Separate build directories preserve other configurations. Optional feature switches and installation instructions appear in the developer README.
## Limits of this validation
Autoconf 2.71, Automake 1.16.5, Libtool 2.4.7 and liblzma 5.4.5 were available for this update. Autotools bootstrapping and distribution checks ran in an isolated source snapshot inside `build/`, so generated scripts did not modify the working source checkout. The strict Autotools check also exposed and corrected a maths-library probe: the real `sqrt` prototype and a volatile function pointer now avoid builtin declaration warnings and optimiser elimination when the caller enables `-Werror`. CMake remains the primary route for recursive developer clones; see the [build-system decision](build-systems.md). Optional codecs remain disabled by default. No system packages were installed.
This run did not test an MPI application across ranks, macOS, Windows, another compiler architecture or a cross-compiled runtime. Threaded tests exercised arbitrary POSIX worker threads and explicit streams. The library itself introduces no MPI runtime dependency. Algorithmic limits and remaining biomcmc-lib functionality are recorded separately in the [migration inventory](migration.md).
