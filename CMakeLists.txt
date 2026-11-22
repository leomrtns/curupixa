cmake_minimum_required(VERSION 3.16)
project(curupixa VERSION 0.1.0 LANGUAGES C)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
include(CheckCCompilerFlag)

if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
  set(CURUPIXA_TOP_LEVEL ON)
else()
  set(CURUPIXA_TOP_LEVEL OFF)
endif()

option(BUILD_SHARED_LIBS "Build shared curupixa libraries" OFF)
option(CURUPIXA_BUILD_GENE_SPECIES "Build gene/species-tree algorithms" ON)
option(CURUPIXA_BUILD_INFERENCE "Build inference algorithms" ${CURUPIXA_BUILD_GENE_SPECIES})
if(CURUPIXA_BUILD_INFERENCE AND NOT CURUPIXA_BUILD_GENE_SPECIES)
  message(FATAL_ERROR "Inference requires CURUPIXA_BUILD_GENE_SPECIES")
endif()
option(CURUPIXA_BUILD_TESTS "Build dependency-free regression tests" ${CURUPIXA_TOP_LEVEL})
option(CURUPIXA_ENABLE_OPENMP "Enable optional OpenMP compilation" OFF)
option(CURUPIXA_ENABLE_SSE42 "Require SSE4.2 on the target CPU" OFF)
option(CURUPIXA_ENABLE_AVX2 "Require AVX2 and FMA on the target CPU" OFF)
option(CURUPIXA_ENABLE_RDRAND "Require RDRAND on the target CPU" OFF)
option(CURUPIXA_ENABLE_ZLIB "Enable gzip file streams" OFF)
option(CURUPIXA_ENABLE_BZIP2 "Enable bzip2 file streams" OFF)
option(CURUPIXA_ENABLE_LZMA "Enable XZ file streams" OFF)
option(CURUPIXA_INSTALL "Generate installation rules" ${CURUPIXA_TOP_LEVEL})

set(curupixa_phylogenetics_sources
  lib/index_arrangement.c
  lib/quasi_random.c
  lib/status.c
  lib/tree.c
  lib/quartet_distance.c
  lib/tanglegram.c
  lib/dtl.c
  lib/bipartition.c
  lib/char_vector.c
  lib/distance_matrix.c
  lib/hashtable.c
  lib/quickselect_quantile.c
  lib/empirical_frequency.c
  lib/file_compression.c
  lib/alignment.c
  lib/nexus_common.c
  lib/gff3_format.c
  lib/parsimony.c
  lib/minimiser_simplex.c
  lib/svd.c
  lib/distance_generator.c
  lib/clustering_goptics.c
  lib/upgma.c
  lib/newick_space.c
  lib/lowlevel_algorithms.c
  lib/spr_distance.c
  lib/ancestry_index.c
  lib/splitset_distances.c
  lib/character_data.c
  lib/character_analysis.c
  lib/substitution_model.c
  lib/likelihood.c
  lib/likelihood_kernel_internal.c
  lib/continuous_traits.c
  lib/reconstruction.c
  lib/probability_sampling.c
  lib/probability_extra.c
  lib/probability_special_internal.c
  lib/prob_distribution.c
  lib/hll.c
  lib/kmerhash.c
  lib/rolling_hash.c
  lib/suffix_tree.c
  lib/topology_distance.c
  lib/topology_randomise.c
  lib/global/global_variable.c
  lib/global/lowlevel.c
  lib/global/maths_and_bits.c
  lib/global/internal_random_constants.c
  lib/global/hash_bytes.c
  lib/global/hash_metro.c
  lib/global/hash_functions_generators.c
  lib/global/hash_functions.c
  lib/global/random_number.c
  lib/global/random_number_generators.c)
set(curupixa_gene_species_sources
  lib/reconciliation.c lib/reconciliation_roots.c lib/taxon_map.c lib/gene_species_distances.c lib/paralog_quartets.c)
set(curupixa_inference_sources lib/inference.c lib/amalgamation.c lib/tree_search.c)
set(curupixa_headers
  lib/curupixa.h lib/phylogenetics.h lib/gene_species.h lib/inference.h lib/amalgamation.h
  lib/tree_search.h lib/paralog_quartets.h lib/gene_species_distances.h lib/quartet_distance.h lib/tanglegram.h lib/dtl.h lib/status.h lib/tree.h lib/bipartition.h lib/char_vector.h
  lib/distance_matrix.h lib/hashtable.h lib/index_arrangement.h lib/quasi_random.h
  lib/quickselect_quantile.h lib/empirical_frequency.h lib/prob_distribution.h lib/file_compression.h
  lib/hll.h lib/kmerhash.h lib/rolling_hash.h lib/suffix_tree.h
  lib/topology_distance.h lib/topology_randomise.h lib/alignment.h lib/nexus_common.h
  lib/gff3_format.h lib/parsimony.h lib/reconciliation.h lib/minimiser_simplex.h lib/svd.h
  lib/distance_generator.h lib/clustering_goptics.h lib/upgma.h lib/newick_space.h lib/lowlevel_algorithms.h
  lib/splitset_distances.h lib/spr_distance.h lib/ancestry_index.h lib/reconciliation_roots.h lib/taxon_map.h
  lib/character_data.h lib/character_analysis.h
  lib/substitution_model.h lib/likelihood.h lib/continuous_traits.h lib/reconstruction.h)
set(curupixa_global_headers
  lib/global/context_types.h lib/global/global_variable.h lib/global/lowlevel.h lib/global/maths_and_bits.h
  lib/global/hash_functions_generators.h
  lib/global/hash_functions.h lib/global/random_number.h lib/global/random_number_generators.h)

list(APPEND curupixa_phylogenetics_sources
  lib/clock.c lib/particle_weights.c lib/tip_correspondence.c lib/stochastic_mapping.c lib/sampling.c lib/quasi_sequence.c
  lib/tree_comparison.c lib/triplet_distance.c lib/gaussian.c lib/phylogenetic_regression.c lib/trait_association.c)
list(APPEND curupixa_inference_sources lib/tree_prior.c lib/molecular_clock.c)
list(APPEND curupixa_headers
  lib/clock.h lib/molecular_clock.h lib/particle_weights.h lib/tip_correspondence.h lib/stochastic_mapping.h
  lib/sampling.h lib/quasi_sequence.h lib/tree_comparison.h
  lib/triplet_distance.h lib/gaussian.h lib/phylogenetic_regression.h lib/trait_association.h lib/tree_prior.h)

# Compile each implementation once. Higher layers link the lower layer directly.
set(curupixa_libraries curupixa_phylogenetics)
add_library(curupixa_phylogenetics ${curupixa_phylogenetics_sources})
add_library(curupixa::phylogenetics ALIAS curupixa_phylogenetics)
set_target_properties(curupixa_phylogenetics PROPERTIES EXPORT_NAME phylogenetics)
if(CURUPIXA_BUILD_GENE_SPECIES)
  add_library(curupixa_gene_species ${curupixa_gene_species_sources})
  add_library(curupixa::gene_species ALIAS curupixa_gene_species)
  set_target_properties(curupixa_gene_species PROPERTIES EXPORT_NAME gene_species)
  target_link_libraries(curupixa_gene_species PUBLIC curupixa_phylogenetics)
  list(APPEND curupixa_libraries curupixa_gene_species)
endif()
if(CURUPIXA_BUILD_INFERENCE)
  add_library(curupixa_inference ${curupixa_inference_sources})
  add_library(curupixa::inference ALIAS curupixa_inference)
  set_target_properties(curupixa_inference PROPERTIES EXPORT_NAME inference)
  target_link_libraries(curupixa_inference PUBLIC curupixa_gene_species)
  list(APPEND curupixa_libraries curupixa_inference)
endif()
add_library(curupixa INTERFACE)
add_library(curupixa::curupixa ALIAS curupixa)
list(GET curupixa_libraries -1 curupixa_top_library)
target_link_libraries(curupixa INTERFACE ${curupixa_top_library})
foreach(library IN LISTS curupixa_libraries)
  target_compile_features(${library} PUBLIC c_std_17)
  set_target_properties(${library} PROPERTIES C_EXTENSIONS OFF POSITION_INDEPENDENT_CODE ON
    VERSION ${PROJECT_VERSION} SOVERSION 0)
  target_compile_definitions(${library} PRIVATE _DEFAULT_SOURCE)
  target_include_directories(${library} PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/lib>
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include/curupixa/global>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/curupixa>)
endforeach()

if(UNIX)
  target_link_libraries(curupixa_phylogenetics PUBLIC m)
endif()
set(CURUPIXA_PC_LIBS_PRIVATE "")
set(CURUPIXA_PC_CFLAGS "")
if(UNIX)
  string(APPEND CURUPIXA_PC_LIBS_PRIVATE " -lm")
endif()
if(CURUPIXA_ENABLE_OPENMP)
  find_package(OpenMP REQUIRED COMPONENTS C)
  # Current implementation carries/synchronises its own state and calls no
  # OpenMP runtime API. Keep this compiler feature private to this build; an
  # installed C++ consumer should not need to enable the C compiler or OpenMP.
  target_link_libraries(curupixa_phylogenetics PRIVATE $<BUILD_INTERFACE:OpenMP::OpenMP_C>)
endif()
if(CURUPIXA_ENABLE_ZLIB)
  find_package(ZLIB REQUIRED)
  set(HAVE_ZLIB 1)
  target_link_libraries(curupixa_phylogenetics PRIVATE ZLIB::ZLIB)
  string(APPEND CURUPIXA_PC_LIBS_PRIVATE " -lz")
endif()
if(CURUPIXA_ENABLE_BZIP2)
  find_package(BZip2 REQUIRED)
  set(HAVE_BZIP2 1)
  target_link_libraries(curupixa_phylogenetics PRIVATE BZip2::BZip2)
  string(APPEND CURUPIXA_PC_LIBS_PRIVATE " -lbz2")
endif()
if(CURUPIXA_ENABLE_LZMA)
  find_package(LibLZMA REQUIRED)
  set(HAVE_LZMA 1)
  target_link_libraries(curupixa_phylogenetics PRIVATE LibLZMA::LibLZMA)
  string(APPEND CURUPIXA_PC_LIBS_PRIVATE " -llzma")
endif()

foreach(feature SSE42 AVX2 RDRAND)
  if(CURUPIXA_ENABLE_${feature})
    if(feature STREQUAL "SSE42")
      set(feature_flags -msse4.2)
      set(HAVE_SSE 1)
    elseif(feature STREQUAL "AVX2")
      set(feature_flags -mavx2 -mfma)
      set(HAVE_AVX 1)
    else()
      set(feature_flags -mrdrnd)
      set(HAVE_RDRND 1)
    endif()
    foreach(flag IN LISTS feature_flags)
      string(MAKE_C_IDENTIFIER "${flag}" flag_identifier)
      check_c_compiler_flag("${flag}" "CURUPIXA_HAS_${flag_identifier}")
      if(NOT CURUPIXA_HAS_${flag_identifier})
        message(FATAL_ERROR "The selected compiler does not support ${flag}")
      endif()
      target_compile_options(curupixa_phylogenetics PRIVATE "${flag}")
    endforeach()
  endif()
endforeach()

configure_file(cmake/config.h.in include/curupixa/global/config.h)

if(CURUPIXA_BUILD_TESTS)
  enable_testing()
  find_package(Threads REQUIRED)
  foreach(test_name core tree runtime statistics probability hashing rng_hash suffix_index topology alignment
      compression gff3 phylogenetics minimiser svd clustering algorithms
      character_analysis substitution_model likelihood likelihood_prepared continuous_traits reconstruction tree_extensions tree_allocation
      spr reconciliation_roots quartet tanglegram dtl amalgamation inference tree_search paralog_quartets comparative tree_priors asr_sampling clocks)
    if(NOT CURUPIXA_BUILD_GENE_SPECIES AND test_name MATCHES "^(algorithms|phylogenetics|spr|reconciliation_roots|paralog_quartets)$")
      continue()
    endif()
    if(NOT CURUPIXA_BUILD_INFERENCE AND test_name MATCHES "^(amalgamation|inference|tree_search|tree_priors|clocks)$")
      continue()
    endif()
    add_executable(test_${test_name} tests/test_${test_name}.c)
    set(test_library curupixa_phylogenetics)
    if(test_name MATCHES "^(algorithms|phylogenetics|spr|reconciliation_roots|paralog_quartets)$")
      set(test_library curupixa_gene_species)
    elseif(test_name MATCHES "^(amalgamation|inference|tree_search|tree_priors|clocks)$")
      set(test_library curupixa_inference)
    endif()
    target_link_libraries(test_${test_name} PRIVATE ${test_library} Threads::Threads)
    if(MSVC)
      target_compile_options(test_${test_name} PRIVATE /UNDEBUG)
    else()
      target_compile_options(test_${test_name} PRIVATE -UNDEBUG)
    endif()
    if(test_name STREQUAL "likelihood_prepared" AND NOT BUILD_SHARED_LIBS AND CMAKE_SYSTEM_NAME STREQUAL "Linux"
        AND CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
      target_compile_definitions(test_${test_name} PRIVATE CRPX_TEST_ALLOCATION_WRAPPERS)
      target_link_options(test_${test_name} PRIVATE
        -Wl,--wrap=malloc -Wl,--wrap=calloc -Wl,--wrap=realloc
        -Wl,--wrap=crpx_character_data_map_tips -Wl,--wrap=crpx_tree_validate)
    endif()
    add_test(NAME ${test_name} COMMAND test_${test_name})
  endforeach()
endif()

# Keep component dependencies discoverable after moving a shared installation.
if(APPLE)
  set_target_properties(${curupixa_libraries} PROPERTIES INSTALL_RPATH "@loader_path")
elseif(UNIX)
  set_target_properties(${curupixa_libraries} PROPERTIES INSTALL_RPATH "$ORIGIN")
endif()

if(CURUPIXA_INSTALL)
  install(FILES
    LICENSES/Apache-2.0.txt
    LICENSES/Blackman-Vigna.txt
    LICENSES/Cephes-BSD-3-Clause.txt
    LICENSES/FastHash-MIT.txt
    LICENSES/HashProspector-Unlicense.txt
    LICENSES/MersenneTwister-BSD-3-Clause.txt
    LICENSES/PCG-MIT.txt
    LICENSES/RNG-hash-NOTICE.txt
    LICENSES/SciPy-BSD-3-Clause.txt
    LICENSES/SipHash-CC0.txt
    LICENSES/SuperFastHash-BSD-3-Clause.txt
    DESTINATION ${CMAKE_INSTALL_DATADIR}/licenses/curupixa)
  install(TARGETS curupixa ${curupixa_libraries} EXPORT curupixaTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
  install(FILES ${curupixa_headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/curupixa)
  install(FILES ${curupixa_global_headers}
    "${PROJECT_BINARY_DIR}/include/curupixa/global/config.h"
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/curupixa/global)
  install(EXPORT curupixaTargets NAMESPACE curupixa::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/curupixa)
  configure_package_config_file(cmake/curupixaConfig.cmake.in
    "${PROJECT_BINARY_DIR}/curupixaConfig.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/curupixa)
  write_basic_package_version_file("${PROJECT_BINARY_DIR}/curupixaConfigVersion.cmake"
    VERSION ${PROJECT_VERSION} COMPATIBILITY SameMajorVersion)
  install(FILES "${PROJECT_BINARY_DIR}/curupixaConfig.cmake"
    "${PROJECT_BINARY_DIR}/curupixaConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/curupixa)
  if(IS_ABSOLUTE "${CMAKE_INSTALL_LIBDIR}")
    set(prefix "${CMAKE_INSTALL_PREFIX}")
    set(libdir "${CMAKE_INSTALL_LIBDIR}")
  else()
    file(RELATIVE_PATH curupixa_pc_prefix_relative
      "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/pkgconfig" "${CMAKE_INSTALL_PREFIX}")
    set(prefix "\${pcfiledir}/${curupixa_pc_prefix_relative}")
    set(libdir "\${prefix}/${CMAKE_INSTALL_LIBDIR}")
  endif()
  set(exec_prefix "\${prefix}")
  if(IS_ABSOLUTE "${CMAKE_INSTALL_INCLUDEDIR}")
    set(includedir "${CMAKE_INSTALL_INCLUDEDIR}")
  else()
    set(includedir "\${prefix}/${CMAKE_INSTALL_INCLUDEDIR}")
  endif()
  set(PACKAGE_VERSION "${PROJECT_VERSION}")
  set(CURUPIXA_PC_TOP curupixa-phylogenetics)
  set(curupixa_pc_files curupixa-phylogenetics.pc)
  configure_file(curupixa-phylogenetics.pc.in curupixa-phylogenetics.pc @ONLY)
  if(CURUPIXA_BUILD_GENE_SPECIES)
    configure_file(curupixa-gene-species.pc.in curupixa-gene-species.pc @ONLY)
    list(APPEND curupixa_pc_files curupixa-gene-species.pc)
    set(CURUPIXA_PC_TOP curupixa-gene-species)
  endif()
  if(CURUPIXA_BUILD_INFERENCE)
    configure_file(curupixa-inference.pc.in curupixa-inference.pc @ONLY)
    list(APPEND curupixa_pc_files curupixa-inference.pc)
    set(CURUPIXA_PC_TOP curupixa-inference)
  endif()
  configure_file(curupixa.pc.in curupixa.pc @ONLY)
  list(APPEND curupixa_pc_files curupixa.pc)
  if(IS_ABSOLUTE "${CMAKE_INSTALL_LIBDIR}")
    # An absolute libdir stays fixed when --prefix changes. Record the actual
    # installation prefix for relative headers while preserving that libdir.
    configure_file(cmake/InstallPkgConfig.cmake.in InstallPkgConfig.cmake @ONLY)
    install(SCRIPT "${PROJECT_BINARY_DIR}/InstallPkgConfig.cmake")
  else()
    foreach(pc_file IN LISTS curupixa_pc_files)
      install(FILES "${PROJECT_BINARY_DIR}/${pc_file}" DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)
    endforeach()
  endif()
endif()
