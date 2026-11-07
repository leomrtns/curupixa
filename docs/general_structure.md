# Source and installed files
## Repository layout
`lib/` contains public object interfaces and their implementations. `lib/global/` contains context, logging, random-number and hash primitives. Private synchronisation and constants stay outside the installed API.

`tests/test_*.c` contains the regression suites. `examples/` demonstrates public API use, and `docs/` contains these guides.
## Build files
Keep generated files below `build/`. CMake places configuration headers and package metadata inside the selected build directory, so different build configurations can coexist.

See [Build-system decision](build-systems.md) for the CMake and Autotools workflows.
## Installed files
An installation places public headers in `${prefix}/include/curupixa` and the library in `${prefix}/${CMAKE_INSTALL_LIBDIR}`. The library directory is usually `lib`, but packaging toolchains can select another location.

Static builds produce `libcurupixa.a`. Shared builds use the platform's normal shared-library names. CMake package configuration is installed below the library directory's `cmake/curupixa` folder; pkg-config metadata is installed as `pkgconfig/curupixa.pc`.
## Using the API
CMake consumers link `curupixa::curupixa`, either through `add_subdirectory()` or `find_package(curupixa CONFIG REQUIRED)`. This target supplies the include paths and link requirements.

Include `<curupixa.h>` when using the build target. Installed consumers can also use `<curupixa/curupixa.h>`. Individual object headers work independently, allowing each source file to include the interfaces it needs.

Objects expose read-only views with documented lifetimes and return recoverable statuses. See [Migration coverage](migration.md) for the relationship between biomcmc-lib and the redesigned interfaces.
