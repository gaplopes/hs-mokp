# FindGUROBI.cmake
# ─────────────────────────────────────────────────────────────────────────────
# Locates the Gurobi Optimizer installation.
#
# Usage:
#   find_package(GUROBI REQUIRED)
#
# This module reads hints from:
#   - GUROBI_HOME  environment variable  (preferred)
#   - GUROBI_DIR   CMake variable
#
# Imported targets:
#   GUROBI::GurobiC     – C   runtime library  (libgurobiXX)
#   GUROBI::GurobiCXX   – C++ wrapper library  (libgurobi_c++)
#
# Result variables:
#   GUROBI_FOUND            – TRUE if both libraries and headers are found.
#   GUROBI_INCLUDE_DIRS     – Path to gurobi_c.h / gurobi_c++.h.
#   GUROBI_C_LIBRARY        – Full path to the C library.
#   GUROBI_CXX_LIBRARY      – Full path to the C++ library.
#   GUROBI_VERSION           – Version string (e.g. "11.0").
# ─────────────────────────────────────────────────────────────────────────────

# ── 1. Determine search root ─────────────────────────────────────────────────

if(NOT GUROBI_HOME AND DEFINED ENV{GUROBI_HOME})
  set(GUROBI_HOME "$ENV{GUROBI_HOME}")
endif()

set(_GUROBI_SEARCH_PATHS
  ${GUROBI_HOME}
  ${GUROBI_DIR}
  /opt/gurobi
  /opt/gurobi*/linux64
  /Library/gurobi*/macos_universal2
  "C:/gurobi*"
)

# ── 2. Find headers ──────────────────────────────────────────────────────────

find_path(GUROBI_INCLUDE_DIRS
  NAMES gurobi_c.h
  HINTS ${_GUROBI_SEARCH_PATHS}
  PATH_SUFFIXES include
)

# ── 3. Auto-detect Gurobi version from the header ────────────────────────────

set(GUROBI_VERSION "")
set(_GUROBI_VERSION_MAJOR "")
set(_GUROBI_VERSION_MINOR "")

if(GUROBI_INCLUDE_DIRS AND EXISTS "${GUROBI_INCLUDE_DIRS}/gurobi_c.h")
  file(STRINGS "${GUROBI_INCLUDE_DIRS}/gurobi_c.h" _GRB_VERSION_LINES
       REGEX "#define GRB_VERSION_(MAJOR|MINOR|TECHNICAL)")
  foreach(_line ${_GRB_VERSION_LINES})
    if(_line MATCHES "#define GRB_VERSION_MAJOR[ \t]+([0-9]+)")
      set(_GUROBI_VERSION_MAJOR "${CMAKE_MATCH_1}")
    elseif(_line MATCHES "#define GRB_VERSION_MINOR[ \t]+([0-9]+)")
      set(_GUROBI_VERSION_MINOR "${CMAKE_MATCH_1}")
    elseif(_line MATCHES "#define GRB_VERSION_TECHNICAL[ \t]+([0-9]+)")
      set(_GUROBI_VERSION_PATCH "${CMAKE_MATCH_1}")
    endif()
  endforeach()
  if(_GUROBI_VERSION_MAJOR AND DEFINED _GUROBI_VERSION_MINOR)
    set(GUROBI_VERSION "${_GUROBI_VERSION_MAJOR}.${_GUROBI_VERSION_MINOR}.${_GUROBI_VERSION_PATCH}")
    # Library name suffix: gurobi120 for v12.0, gurobi110 for v11.0, etc.
    set(_GUROBI_LIB_SUFFIX "${_GUROBI_VERSION_MAJOR}${_GUROBI_VERSION_MINOR}")
  endif()
endif()

# ── 4. Find C runtime library (libgurobiXX) ──────────────────────────────────

find_library(GUROBI_C_LIBRARY
  NAMES
    gurobi${_GUROBI_LIB_SUFFIX}
    gurobi                          # fallback
  HINTS ${_GUROBI_SEARCH_PATHS}
  PATH_SUFFIXES lib
)

# ── 5. Find C++ wrapper library (libgurobi_c++) ──────────────────────────────

find_library(GUROBI_CXX_LIBRARY
  NAMES gurobi_c++
  HINTS ${_GUROBI_SEARCH_PATHS}
  PATH_SUFFIXES lib
)

# ── 6. Handle standard find_package arguments ────────────────────────────────

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GUROBI
  REQUIRED_VARS
    GUROBI_INCLUDE_DIRS
    GUROBI_C_LIBRARY
    GUROBI_CXX_LIBRARY
  VERSION_VAR GUROBI_VERSION
)

# ── 7. Create imported targets ───────────────────────────────────────────────

if(GUROBI_FOUND AND NOT TARGET GUROBI::GurobiC)
  add_library(GUROBI::GurobiC UNKNOWN IMPORTED)
  set_target_properties(GUROBI::GurobiC PROPERTIES
    IMPORTED_LOCATION             "${GUROBI_C_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${GUROBI_INCLUDE_DIRS}"
  )
endif()

if(GUROBI_FOUND AND NOT TARGET GUROBI::GurobiCXX)
  add_library(GUROBI::GurobiCXX UNKNOWN IMPORTED)
  set_target_properties(GUROBI::GurobiCXX PROPERTIES
    IMPORTED_LOCATION             "${GUROBI_CXX_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${GUROBI_INCLUDE_DIRS}"
    INTERFACE_LINK_LIBRARIES      GUROBI::GurobiC
  )
endif()

mark_as_advanced(
  GUROBI_INCLUDE_DIRS
  GUROBI_C_LIBRARY
  GUROBI_CXX_LIBRARY
)
