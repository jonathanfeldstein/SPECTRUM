cmake_minimum_required(VERSION 3.16)

project(SPECTRUM LANGUAGES CXX)

option(ENABLE_DEBUG "Enable profiling/debug code" OFF)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_BUILD_TYPE Release)

add_compile_options(
    -O3
    -march=native
)

# ------------------------------------------------------------------
# Homebrew packages (Apple Silicon)
# ------------------------------------------------------------------

if(APPLE)
    list(APPEND CMAKE_PREFIX_PATH
        "/opt/homebrew/opt/libomp"
        "/opt/homebrew/opt/jemalloc"
    )
endif()

# ------------------------------------------------------------------
# Dependencies
# ------------------------------------------------------------------

find_package(OpenMP REQUIRED)
find_package(Boost REQUIRED)

find_path(
    JEMALLOC_INCLUDE_DIR
    NAMES jemalloc/jemalloc.h
)

find_library(
    JEMALLOC_LIBRARY
    NAMES jemalloc
)

if(NOT JEMALLOC_INCLUDE_DIR OR NOT JEMALLOC_LIBRARY)
    message(FATAL_ERROR "Could not find jemalloc.")
endif()

message(STATUS "Found OpenMP")
message(STATUS "Found Boost")
message(STATUS "Found jemalloc:")
message(STATUS "    include: ${JEMALLOC_INCLUDE_DIR}")
message(STATUS "    library: ${JEMALLOC_LIBRARY}")

# ------------------------------------------------------------------
# Project paths
# ------------------------------------------------------------------

set(EIGEN_PATH "${CMAKE_SOURCE_DIR}/dependencies/eigen-3.4.0")
set(SPECTRA_PATH "${CMAKE_SOURCE_DIR}/dependencies/spectra-1.0.1/include")

# ------------------------------------------------------------------
# Executable
# ------------------------------------------------------------------

add_executable(
    SPECTRUM

    src/main.cpp

    src/Utils/vector_utils.h
    src/Utils/map_utils.cpp
    src/Utils/map_utils.h
    src/Utils/set_utils.h
    src/Utils/profiling_utils.h
    src/Utils/read_files_utils.cpp
    src/Utils/read_files_utils.h
    src/Utils/exceptions.cpp
    src/Utils/exceptions.h

    src/Fragments/Fragments.cpp
    src/Fragments/Fragments.h

    src/CandidateRules/CandidateRules.cpp
    src/CandidateRules/CandidateRules.h

    src/HyperGraphs/Relation.h
    src/HyperGraphs/HyperGraph.cpp
    src/HyperGraphs/HyperGraph.h
)

target_include_directories(
    SPECTRUM
    PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${EIGEN_PATH}
        ${SPECTRA_PATH}
        ${Boost_INCLUDE_DIRS}
        ${JEMALLOC_INCLUDE_DIR}
)

target_link_libraries(
    SPECTRUM
    PRIVATE
        OpenMP::OpenMP_CXX
        ${JEMALLOC_LIBRARY}
)

if(ENABLE_DEBUG)
    target_compile_definitions(
        SPECTRUM
        PRIVATE
        ENABLE_DEBUG
    )
endif()
# ------------------------------------------------------------------
# Random-walk variant of SPECTRUM (ablation): fragments found by random walks
# instead of pattern mining; everything else identical
# ------------------------------------------------------------------

add_executable(
    SPECTRUM_RW

    src/main.cpp

    src/Utils/vector_utils.h
    src/Utils/map_utils.cpp
    src/Utils/map_utils.h
    src/Utils/set_utils.h
    src/Utils/profiling_utils.h
    src/Utils/read_files_utils.cpp
    src/Utils/read_files_utils.h
    src/Utils/exceptions.cpp
    src/Utils/exceptions.h

    src/Fragments/Fragments_rw.cpp
    src/Fragments/Fragments_rw.h

    src/CandidateRules/CandidateRules.cpp
    src/CandidateRules/CandidateRules.h

    src/HyperGraphs/Relation.h
    src/HyperGraphs/HyperGraph.cpp
    src/HyperGraphs/HyperGraph.h
)

target_include_directories(
    SPECTRUM_RW
    PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${EIGEN_PATH}
        ${SPECTRA_PATH}
        ${Boost_INCLUDE_DIRS}
        ${JEMALLOC_INCLUDE_DIR}
)

target_link_libraries(
    SPECTRUM_RW
    PRIVATE
        OpenMP::OpenMP_CXX
        ${JEMALLOC_LIBRARY}
)

if(ENABLE_DEBUG)
    target_compile_definitions(
        SPECTRUM_RW
        PRIVATE
        ENABLE_DEBUG
    )
endif()

target_compile_definitions(
    SPECTRUM_RW
    PRIVATE
        USE_RANDOM_WALKS
)

# ------------------------------------------------------------------
# Evaluator (NCRL protocol)
# ------------------------------------------------------------------

add_executable(
    spectrum_eval

    evaluator/main.cpp
)

target_link_libraries(
    spectrum_eval
    PRIVATE
        OpenMP::OpenMP_CXX
)

# Scores must be computed with the same floating-point operations as NCRL (numpy/scipy),
# otherwise exact ties between candidates break differently: no fused multiply-add.
target_compile_options(
    spectrum_eval
    PRIVATE
        -ffp-contract=off
)
