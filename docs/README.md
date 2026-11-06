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
- [Source and installed files](general_structure.md)
- [Module coverage and intentional changes](migration.md)
- [Compact tree backbone and analysis workspaces](tree-design.md)
- [Runtime, random streams, hashing and annotation I/O](runtime-and-hashing.md)
- [Future likelihood backend](likelihood-design.md)
- [Local validation and platform limits](validation.md)
- [Source references and attribution](references.md)
Likelihood calculation remains a future backend. New code should keep numerical buffers and caches in separate worker-owned workspaces, using stable tree node IDs and immutable alignment patterns.
## Licence
Copyright (C) 2022–today Leonardo de Oliveira Martins. Curupixa uses GPL-3.0-or-later; incorporated source files retain their attribution. See the repository licence and references for details.
