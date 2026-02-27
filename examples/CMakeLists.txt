# ─────────────────────────────────────────────────────────────────────────────
# HS-MOCO examples
# ─────────────────────────────────────────────────────────────────────────────

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
target_link_libraries(example_mokp PRIVATE hs_moco::hs_moco)
target_compile_definitions(example_mokp PRIVATE
  MOBKP_INSTANCES_DIR="${MOBKP_INSTANCES_DIR}"
)
