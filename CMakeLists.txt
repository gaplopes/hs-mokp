# ─────────────────────────────────────────────────────────────────────────────
# HS-MOCO — Hypervolume Scalarization for Multiobjective Combinatorial
#            Optimization (header-only C++17 library)
# ─────────────────────────────────────────────────────────────────────────────
cmake_minimum_required(VERSION 3.16)

project(hs_moco
  VERSION   0.1.0
  LANGUAGES CXX
  DESCRIPTION "Header-only C++17 library for hypervolume scalarization in MOCO"
)

# ── Options ──────────────────────────────────────────────────────────────────

# Build examples only when this is the top-level project (not a subdirectory).
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  set(_HS_MOCO_IS_TOPLEVEL ON)
else()
  set(_HS_MOCO_IS_TOPLEVEL OFF)
endif()

option(HS_MOCO_BUILD_EXAMPLES "Build example executables" ${_HS_MOCO_IS_TOPLEVEL})

# ── Find Gurobi ─────────────────────────────────────────────────────────────

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules")
find_package(GUROBI REQUIRED)

message(STATUS "[hs_moco] Gurobi ${GUROBI_VERSION} found at ${GUROBI_INCLUDE_DIRS}")

# ── Find mooutils ───────────────────────────────────────────────────────────

find_package(mooutils QUIET)
if(NOT mooutils_FOUND)
  include(FetchContent)
  FetchContent_Declare(
    mooutils
    GIT_REPOSITORY https://github.com/adbjesus/mooutils.git
    GIT_TAG        main
  )
  FetchContent_MakeAvailable(mooutils)
  if(NOT TARGET mooutils::mooutils)
    add_library(mooutils::mooutils ALIAS mooutils)
  endif()
  message(STATUS "[hs_moco] mooutils fetched from GitHub")
else()
  message(STATUS "[hs_moco] mooutils found (installed)")
endif()

# ── Header-only INTERFACE library ────────────────────────────────────────────

add_library(hs_moco INTERFACE)
add_library(hs_moco::hs_moco ALIAS hs_moco)

target_compile_features(hs_moco INTERFACE cxx_std_20)

target_include_directories(hs_moco
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_link_libraries(hs_moco
  INTERFACE
    GUROBI::GurobiCXX
    mooutils::mooutils
)

# ── Examples ─────────────────────────────────────────────────────────────────

if(HS_MOCO_BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()

# ── Install rules ────────────────────────────────────────────────────────────

include(GNUInstallDirs)

install(
  TARGETS hs_moco
  EXPORT  hs_moco-targets
)

install(
  DIRECTORY   include/hs_moco
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
)

install(
  EXPORT      hs_moco-targets
  FILE        hs_moco-targets.cmake
  NAMESPACE   hs_moco::
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/hs_moco"
)

# ── Package config for find_package(hs_moco) ────────────────────────────────

include(CMakePackageConfigHelpers)

configure_package_config_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/hs_moco-config.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/hs_moco-config.cmake"
  INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/hs_moco"
)

write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/hs_moco-config-version.cmake"
  VERSION       ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)

install(
  FILES
    "${CMAKE_CURRENT_BINARY_DIR}/hs_moco-config.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/hs_moco-config-version.cmake"
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules/FindGUROBI.cmake"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/hs_moco"
)
