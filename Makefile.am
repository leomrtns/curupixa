ACLOCAL_AMFLAGS = -I m4
SUBDIRS = lib tests
pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = curupixa.pc curupixa-phylogenetics.pc
if BUILD_GENE_SPECIES
pkgconfig_DATA += curupixa-gene-species.pc
endif
if BUILD_INFERENCE
pkgconfig_DATA += curupixa-inference.pc
endif
EXTRA_DIST = CMakeLists.txt cmake/config.h.in cmake/curupixaConfig.cmake.in cmake/InstallPkgConfig.cmake.in \
  curupixa.pc.in curupixa-phylogenetics.pc.in curupixa-gene-species.pc.in curupixa-inference.pc.in \
  docs/261107.mit-licence-audit.md docs/261107.licence-provenance-findings.csv \
  docs/random-hashing.md docs/probability-numerics.md docs/molecular-clocks.md docs/comparative-methods.md docs/tree-sampling.md \
  docs/runtime-and-hashing.md docs/migration.md docs/api-migration.csv docs/likelihood-design.md \
  docs/tree-design.md docs/README.md docs/general_structure.md docs/_sidebar.md docs/references.md docs/validation.md \
  docs/language-bindings.md docs/library-architecture.md docs/comparison-and-inference.md docs/build-systems.md docs/character-analysis.md docs/likelihood-reference-audit.md \
  docs/cladistic-distances.md docs/library-footprint.md docs/continuous-traits.md docs/reconstruction.md docs/tree-extensions.md docs/downstream-trait-api.md \
  tests/check_probability_reference.py examples/python_api.py examples/discrete_characters.c examples/gene_species_distances.c examples/trait_analysis.c tests/check_headers.py

dist_pkgdata_DATA = LICENSES/Apache-2.0.txt \
  LICENSES/Blackman-Vigna.txt \
  LICENSES/Cephes-BSD-3-Clause.txt \
  LICENSES/FastHash-MIT.txt \
  LICENSES/HashProspector-Unlicense.txt \
  LICENSES/MersenneTwister-BSD-3-Clause.txt \
  LICENSES/PCG-MIT.txt \
  LICENSES/RNG-hash-NOTICE.txt \
  LICENSES/SciPy-BSD-3-Clause.txt \
  LICENSES/SipHash-CC0.txt \
  LICENSES/SuperFastHash-BSD-3-Clause.txt
