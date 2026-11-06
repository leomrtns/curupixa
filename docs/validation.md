# Local validation
Validated on 6 November 2026 with GCC 13.3.0, CMake 3.28.3 and Python 3.12.8 on Linux. All generated builds, test fixtures and temporary consumers stayed inside this checkout's `build/` directory.
## Build and regression matrix
| Configuration | Result |
| --- | --- |
| Static Release, optional dependencies disabled | 16/16 CTest suites passed. |
| Shared Release, zlib 1.3 and bzip2 1.0.8 enabled | 16/16 suites passed. |
| Debug, `-Wall -Wextra -Werror`, AddressSanitizer and UndefinedBehaviourSanitizer | 16/16 suites passed. |
| OpenMP 4.5 enabled | 16/16 suites passed. |

The suites cover containers, trees, runtime ownership/threading, statistics, probability functions, hashing, suffix indexing, topology operations/construction/collections, alignment patterns and distances, compression, GFF3, parsimony/reconciliation, minimisation, SVD, clustering and assignment/edit/split distances. Scientific regressions include probability tails and quantile support, BioNJ branch distances, OPTICS border points, weighted parsimony overflow, mapped multi-copy splits and scaled matrix reconstruction. Parser tests cover malformed input, quotes/comments, deep trees, large polytomies, locale-independent decimal scores and compressed-stream integrity.
The sanitised run used `ASAN_OPTIONS=detect_leaks=0`: LeakSanitizer cannot operate under the execution environment's ptrace wrapper. Address/undefined-behaviour checks ran; this run does not establish leak-check coverage.
## Consumer and package checks
- Fresh installed static/shared libraries linked and ran from separate C and C++ consumer projects.
- All 40 installed headers compiled independently as C17 and C++17 with strict warnings. Private runtime/constant headers were absent from the installed API.
- A Python `ctypes` smoke test created a tree, added and labelled tips, wrote `(A,B);`, released the returned string/tree and queried a status string.
- A parent project built and ran through `add_subdirectory()` and `curupixa::curupixa`, with curupixa tests and installation disabled by default.
- CMake and pkg-config consumers worked with an overridden installation prefix. pkg-config also worked after moving a prefix with a nested relative library directory, and with a separately configured absolute library directory.
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
Autoconf/Automake/Libtool bootstrapping tools were unavailable, so the updated Autotools source definitions were neither regenerated nor validated. CMake is the validated build route. The optional XZ implementation was uncompiled because liblzma development headers were unavailable; XZ remains disabled by default. No system packages were installed.
This run did not test an MPI application across ranks, macOS, Windows, another compiler architecture or a cross-compiled runtime. Threaded tests exercised arbitrary POSIX worker threads and explicit streams. The library itself introduces no MPI runtime dependency. Algorithmic limits and remaining biomcmc-lib functionality are recorded separately in the [migration inventory](migration.md).
