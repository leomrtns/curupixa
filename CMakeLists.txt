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

option(BUILD_SHARED_LIBS "Build a shared curupixa library" OFF)
option(CURUPIXA_BUILD_TESTS "Build dependency-free regression tests" ${CURUPIXA_TOP_LEVEL})
option(CURUPIXA_ENABLE_OPENMP "Enable optional OpenMP compilation" OFF)
option(CURUPIXA_ENABLE_SSE42 "Require SSE4.2 on the target CPU" OFF)
option(CURUPIXA_ENABLE_AVX2 "Require AVX2 and FMA on the target CPU" OFF)
option(CURUPIXA_ENABLE_RDRAND "Require RDRAND on the target CPU" OFF)
option(CURUPIXA_ENABLE_ZLIB "Enable gzip file streams" OFF)
option(CURUPIXA_ENABLE_BZIP2 "Enable bzip2 file streams" OFF)
option(CURUPIXA_ENABLE_LZMA "Enable XZ file streams" OFF)
option(CURUPIXA_INSTALL "Generate installation rules" ${CURUPIXA_TOP_LEVEL})

set(curupixa_sources
  lib/index_arrangement.c
  lib/quasi_random.c
  lib/status.c
  lib/tree.c
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
  lib/reconciliation.c
  lib/minimiser_simplex.c
  lib/svd.c
  lib/distance_generator.c
  lib/clustering_goptics.c
  lib/upgma.c
  lib/newick_space.c
  lib/lowlevel_algorithms.c
  lib/spr_distance.c
  lib/ancestry_index.c
  lib/reconciliation_roots.c
  lib/taxon_map.c
  lib/splitset_distances.c
  lib/character_data.c
  lib/character_analysis.c
  lib/substitution_model.c
  lib/likelihood.c
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
  lib/global/hash_functions_generators.c
  lib/global/hash_functions.c
  lib/global/random_number.c
  lib/global/random_number_generators.c)
set(curupixa_headers
  lib/curupixa.h lib/status.h lib/tree.h lib/bipartition.h lib/char_vector.h
  lib/distance_matrix.h lib/hashtable.h lib/index_arrangement.h lib/quasi_random.h
  lib/quickselect_quantile.h lib/empirical_frequency.h lib/prob_distribution.h lib/file_compression.h
  lib/hll.h lib/kmerhash.h lib/rolling_hash.h lib/suffix_tree.h
  lib/topology_distance.h lib/topology_randomise.h lib/alignment.h lib/nexus_common.h
  lib/gff3_format.h lib/parsimony.h lib/reconciliation.h lib/minimiser_simplex.h lib/svd.h
  lib/distance_generator.h lib/clustering_goptics.h lib/upgma.h lib/newick_space.h lib/lowlevel_algorithms.h
  lib/splitset_distances.h lib/spr_distance.h lib/ancestry_index.h lib/reconciliation_roots.h lib/taxon_map.h
  lib/character_data.h lib/character_analysis.h
  lib/substitution_model.h lib/likelihood.h)
set(curupixa_global_headers
  lib/global/context_types.h lib/global/global_variable.h lib/global/lowlevel.h lib/global/maths_and_bits.h
  lib/global/hash_functions_generators.h
  lib/global/hash_functions.h lib/global/random_number.h lib/global/random_number_generators.h)

add_library(curupixa ${curupixa_sources})
add_library(curupixa::curupixa ALIAS curupixa)
target_compile_features(curupixa PUBLIC c_std_17)
set_target_properties(curupixa PROPERTIES
  C_EXTENSIONS OFF
  POSITION_INDEPENDENT_CODE ON
  VERSION ${PROJECT_VERSION}
  SOVERSION 0)
target_compile_definitions(curupixa PRIVATE _DEFAULT_SOURCE)
target_include_directories(curupixa PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/lib>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include/curupixa/global>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/curupixa>)

if(UNIX)
  target_link_libraries(curupixa PUBLIC m)
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
  target_link_libraries(curupixa PRIVATE $<BUILD_INTERFACE:OpenMP::OpenMP_C>)
endif()
if(CURUPIXA_ENABLE_ZLIB)
  find_package(ZLIB REQUIRED)
  set(HAVE_ZLIB 1)
  target_link_libraries(curupixa PRIVATE ZLIB::ZLIB)
  string(APPEND CURUPIXA_PC_LIBS_PRIVATE " -lz")
endif()
if(CURUPIXA_ENABLE_BZIP2)
  find_package(BZip2 REQUIRED)
  set(HAVE_BZIP2 1)
  target_link_libraries(curupixa PRIVATE BZip2::BZip2)
  string(APPEND CURUPIXA_PC_LIBS_PRIVATE " -lbz2")
endif()
if(CURUPIXA_ENABLE_LZMA)
  find_package(LibLZMA REQUIRED)
  set(HAVE_LZMA 1)
  target_link_libraries(curupixa PRIVATE LibLZMA::LibLZMA)
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
      target_compile_options(curupixa PRIVATE "${flag}")
    endforeach()
  endif()
endforeach()

configure_file(cmake/config.h.in include/curupixa/global/config.h)

if(CURUPIXA_BUILD_TESTS)
  enable_testing()
  find_package(Threads REQUIRED)
  foreach(test_name core tree runtime statistics probability hashing suffix_index topology alignment
      compression gff3 phylogenetics minimiser svd clustering algorithms
      character_analysis substitution_model likelihood spr reconciliation_roots)
    add_executable(test_${test_name} tests/test_${test_name}.c)
    target_link_libraries(test_${test_name} PRIVATE curupixa::curupixa Threads::Threads)
    if(MSVC)
      target_compile_options(test_${test_name} PRIVATE /UNDEBUG)
    else()
      target_compile_options(test_${test_name} PRIVATE -UNDEBUG)
    endif()
    add_test(NAME ${test_name} COMMAND test_${test_name})
  endforeach()
endif()

if(CURUPIXA_INSTALL)
  install(TARGETS curupixa EXPORT curupixaTargets
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
  configure_file(curupixa.pc.in curupixa.pc @ONLY)
  if(IS_ABSOLUTE "${CMAKE_INSTALL_LIBDIR}")
    # An absolute libdir stays fixed when --prefix changes. Record the actual
    # installation prefix for relative headers while preserving that libdir.
    configure_file(cmake/InstallPkgConfig.cmake.in InstallPkgConfig.cmake @ONLY)
    install(SCRIPT "${PROJECT_BINARY_DIR}/InstallPkgConfig.cmake")
  else()
    install(FILES "${PROJECT_BINARY_DIR}/curupixa.pc" DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)
  endif()
endif()
