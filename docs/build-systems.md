# Build-system decision
CMake is the primary build system for curupixa. The developer workflow starts from a recursive Git clone and commonly embeds curupixa in another program. CMake builds that checkout directly and provides the same `curupixa::curupixa` target for embedded and installed use. It needs no generated source files, bootstrap step or network access during configuration. Keep Autotools as a tested alternative for existing downstream Autotools projects and conventional source-release tarballs.
## Comparison for this library
| Requirement | CMake | Autotools |
| --- | --- | --- |
| Fresh Git checkout | CMake and C17 compiler; no bootstrap. | Autoconf, Automake and Libtool generate the build scripts, then configure/make. |
| Embed a submodule | `add_subdirectory()` and an imported-style target carrying usage requirements. | Recursive configuration or downstream Makefile integration with the convenience library. |
| Installed developer use | Exported CMake config/target plus pkg-config. | pkg-config and normal installed library/headers. |
| Static/shared and Python extension linkage | Static PIC by default; shared switch; dependencies propagated by targets. | Libtool supplies static/shared and convenience libraries; use `--with-pic` when an archive will enter a shared extension. |
| Conda, Galaxy dependency package or container prefix | Non-interactive configure/build/install with an explicit prefix; supports staging. | Non-interactive configure/make/install with an explicit prefix; supports staging. |
| Source-release archive | Compiler and CMake required by the documented workflow. | Distributed generated `configure` supports compiler/make without installing bootstrap tools. |
| Maintenance | Source/header/test lists and optional dependency switches. | Equivalent lists/switches plus generated-script bootstrap and distribution checks. |

Both systems successfully build and test the same implementation. The recommendation follows curupixa's developer/submodule use case, rather than a numerical performance difference: both compile the same C sources. New downstream projects should use CMake. Existing Autotools integrations do not need to migrate solely for this update.
OpenMP is an optional private compilation feature. The current C implementation owns its worker state and calls no OpenMP runtime functions, so installed packages do not force OpenMP discovery, consumer compiler flags or an extra C compiler in a C++-only consuming project. Applications that use OpenMP select it in their own build. Introducing OpenMP runtime calls inside the library would require updating the exported link requirements.
## Optional codecs
The dependency-free default uses raw file streams. Enable individual codecs with `CURUPIXA_ENABLE_ZLIB`, `CURUPIXA_ENABLE_BZIP2` and `CURUPIXA_ENABLE_LZMA` in CMake, or `--enable-zlib`, `--enable-bzip2` and `--enable-lzma` in Autotools. The build checks for the selected development headers and libraries and fails clearly if they are unavailable. It downloads or installs nothing.
Static installed consumers inherit codec link requirements. Shared CMake consumers require the codec runtime libraries but need no codec development headers: the codec types remain private to the implementation. pkg-config places codec requirements in `Libs.private`, used with `--static`.
## Verification
Autoconf 2.71, Automake 1.16.5, Libtool 2.4.7 and liblzma 5.4.5 are available in the current environment. Bootstrapping, out-of-source configuration, compilation, tests, installation and `make distcheck` all succeeded with gzip, bzip2 and XZ enabled. Bootstrapping ran in an isolated source snapshot below `build/`, preserving the working checkout's generated files.
The resulting source-release tarball also configured, compiled and passed all 19 suites through CMake with XZ enabled, checking that the archive includes the CMake templates, private headers and new character-analysis sources.
CMake all-codec Release and strict AddressSanitizer/UndefinedBehaviourSanitizer builds passed the regression suites. Codec tests cover round trips, concatenated members, empty streams, mixed line/block reads, truncation, trailing garbage and payload corruption. Independent read/write interoperability checks with Python's gzip, bz2 and lzma interfaces passed for a 1,000,003-byte input. See the [validation report](validation.md) for the full implementation test matrix and remaining platform limits.
## Reproduce an Autotools release check
From a disposable source snapshot inside the project's `build/` directory, bootstrap with `autoreconf -fi`, then configure into a separate objects directory:
```sh
../source/configure --prefix="$PWD/install" --enable-zlib --enable-bzip2 --enable-lzma
make -j4
make -j4 check
make install
make -j4 distcheck DISTCHECK_CONFIGURE_FLAGS='--enable-zlib --enable-bzip2 --enable-lzma'
```
Use a prefix inside the project or the package manager's staging area during local checks. Ordinary downstream builds need neither a system installation nor administrator privileges.
## Primary references
CMake documents [subdirectory builds](https://cmake.org/cmake/help/latest/command/add_subdirectory.html) and [exported/relocatable package targets](https://cmake.org/cmake/help/latest/guide/importing-exporting/index.html). GNU Automake documents the [distribution check workflow](https://www.gnu.org/software/automake/manual/html_node/Checking-the-Distribution.html).
