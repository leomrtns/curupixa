# Headers, linking and memory
Curupixa keeps algorithm buffers in explicitly constructed objects and separates algorithms into compiled source modules. Adding another header/function expands the API and library file; it does not add per-node fields or allocate its analysis buffers when a program loads the library.
## Public dependency structure
Each public header has an include guard and declares its own required standard and curupixa types. Guards make repeated inclusion safe. The include graph is also acyclic: algorithm headers depend on primitive objects, while primitive objects do not include algorithm headers. The umbrella `curupixa.h` is a convenient top-level aggregate.

For example, `reconciliation_roots.h` includes `ancestry_index.h`, which includes `tree.h`, which includes standard types and `status.h`. `spr_distance.h` includes `splitset_distances.h`, which includes `tree.h` and `bipartition.h`. None of these paths imports likelihood, SVD, file compression or the global runtime. `context_types.h` provides the context forward declaration for RNG/hash declarations. Raw RNG, hash-generator and bit-operation headers now use only their required standard headers; their source files include implementation dependencies explicitly.

The installed API has 51 headers, including the generated configuration header. All compile independently, twice in succession, as C17 and C++17 with strict warnings. Combined tests use forward, reverse and seeded shuffled orders. The optional maintainer check also rejects missing local includes and cycles:
```sh
python3 tests/check_headers.py build/package-prefix/include/curupixa
```
Python is needed only to run this check. Configuration, library compilation and normal regression suites have no Python dependency. Individual headers improve compile time and make dependencies clearer; including the umbrella does not itself select more compiled library code.
## Static and shared linking
The default static archive contains one object member per source module. The linker selects members needed to resolve the application's symbols. A tree-only program does not link likelihood, reconciliation, SPR, SVD or GFF3 objects merely because it includes the umbrella. Selection has **source-module granularity**: requesting one function may bring other functions from that same object. This is the normal [GNU linker archive behaviour](https://www.sourceware.org/binutils/docs/ld/Options.html).

Shared builds place the available code and data in one shared object. Linux maps the file into virtual memory; resident memory depends on which pages have been accessed. File size, mapped virtual size, resident size and private dirty pages measure different things. Dependencies such as libm and optionally enabled compression libraries have their own mappings. Curupixa has no library-load constructor that creates trees, random contexts or analysis workspaces. The dynamic loader still uses memory for mappings, relocations and its own bookkeeping.

Static builds already avoid unused algorithm modules. A downstream toolchain can additionally enable function/data sections and linker garbage collection when reducing within-module code matters. Such flags are platform-specific and remain under downstream control. Splitting every algorithm into separate shared libraries would add packaging/dependency complexity without reducing the dominant numerical workspace allocations.
## Trait-API update measurement
The branchtraits-driven additions retain the 12-byte node layout and add one 64-bit topology generation per tree. The default Release shared library now measures 381,760 bytes on disk (about 373 KiB), mapping 336 KiB virtual, 120 KiB resident and 8 KiB private dirty immediately after loading in the measured process. These figures exclude optional codec libraries and other process dependencies.

Tree-only static clients have identical text/data/bss figures of 22,572/712/8 bytes whether they include `tree.h` or the umbrella. Their executable files differ by eight bytes of file metadata (35,320/35,328 bytes). Neither includes likelihood, ASR, Brownian, PIC or reconciliation symbols. New tree primitives increase the linked tree object, while the optional analysis modules remain absent. The earlier measurements below record the pre-trait-API baseline.

The allocation-instrumented 100-rate test records zero allocation/mapping/validation calls with an already warmed prepared binding. The ordinary repeated one-off interface records 1,000/100/100 respectively under the tested static build. These count library call sites rather than arbitrary loader/libm internals, and measure setup work rather than guaranteeing a machine-independent timing speedup. The prepared path pays its snapshot/preparation cost once; one-off callers do not allocate a snapshot.
## Measured earlier footprint
On Linux x86-64, GCC 13.3.0, Release `-O3 -DNDEBUG`, default hardware options and optional codecs disabled, on 7 November 2026:
| Measurement | Value |
| --- | --- |
| Unstripped shared-library file | 345,624 bytes (about 338 KiB) |
| `.text` section | 199,921 bytes |
| Read-only `.rodata` | 17,360 bytes |
| `.data` / `.bss` | 8 / 8 bytes, excluding relocation/dynamic-loader sections |
| Library mappings immediately after `ctypes.CDLL(..., RTLD_NOW)` | 308 KiB virtual; 120 KiB resident; 8 KiB private dirty |
| Static tree-constructor/destructor client, individual header | 30,752-byte executable; text/data/bss 17,540/704/8 bytes |
| Same static client, umbrella header | Identical executable size and text/data/bss; unused algorithm symbols absent |

The mapping figures cover this shared object alone, excluding Python, libm and loader bookkeeping. They are one local observation, not a fixed requirement: compiler, platform, instrumentation, enabled dependencies, loader policy and page use change them. Sanitised/debug builds are larger. No analysis objects had been constructed for the load measurement. `size -A` and `/proc/self/smaps` supplied the section and mapping figures.
## Where working memory grows
| Object or calculation | Storage driver |
| --- | --- |
| Tree connectivity | 12 bytes per allocated node on the tested ABI; optional labels and lengths have separate allocations. |
| Tree traversal/sister/depth/update workspace | Only arrays requested by flags; capacity reused until freed. |
| Ancestry snapshot | $O(S\log S)$ indices, immutable and shareable across gene-family workers. |
| All-root reconciliation | $O(S+G)$ scratch/result storage beyond the ancestry snapshot; released at the end of the call except owned result scores. |
| Split extraction | $O(V\lceil T/64\rceil)$ descendant bitsets. |
| RF-only split comparison | Hash matching of split words, without the assignment matrix. |
| Assignment / greedy dSPR | Quadratic scalar cost matrix; dSPR additionally retains linear lists of split bitsets. No all-pairs disagreement bitsets. |
| Likelihood and parsimony | Character data, model and per-worker partial buffers; see [character analysis](character-analysis.md) for streaming and stored-partial options. |
| Dense SVD, packed distances, NJ/BioNJ and collections | Explicit dense/packed matrices or owned trees, constructed only on request. |

For MPI, each rank owns ordinary local objects. Within a process, workers can share immutable data/models/ancestry snapshots, and use separate mutable workspaces. Reused capacities stay allocated intentionally; `crpx_spr_workspace_clear()` releases its large buffers, and destructors release all object-owned arrays. Storage-size getters include retained capacity while excluding allocator overhead. Callers control how many simultaneous analyses exist.
