# ─────────────────────────────────────────────────────────────────────────────
# HS-MOKP examples
# ─────────────────────────────────────────────────────────────────────────────

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

# ── Fetch MOKP benchmark instances ───────────────────────────────────────────

include(FetchContent)
FetchContent_Declare(
  mobkp_instances
  GIT_REPOSITORY https://github.com/gaplopes/mobkp-instances.git
  GIT_TAG        main
  GIT_SHALLOW    TRUE
  SOURCE_SUBDIR  _none_  # Don't add_subdirectory; we only want the data files
)
FetchContent_MakeAvailable(mobkp_instances)

set(MOBKP_INSTANCES_DIR "${mobkp_instances_SOURCE_DIR}/instances")
message(STATUS "[hs_moco] MOKP instances at ${MOBKP_INSTANCES_DIR}")

# ── Example: MOKP ────────────────────────────────────────────────────────────

add_executable(example_mokp mokp.cpp)
target_link_libraries(example_mokp PRIVATE hs_moco::hs_moco mooutils::mooutils)
target_compile_definitions(example_mokp PRIVATE
  MOBKP_INSTANCES_DIR="${MOBKP_INSTANCES_DIR}"
)
