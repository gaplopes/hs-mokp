#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# test_mokp.sh — Run the MOKP strategy-comparison example across instances.
#
# Usage:
#   ./examples/test_mokp.sh                  # defaults
#   ./examples/test_mokp.sh --max-k 3        # override max reference-point set size
#   ./examples/test_mokp.sh --build-dir build # custom build directory
#   ./examples/test_mokp.sh --quick           # small subset for a quick smoke test
#
# The script:
#   1. Builds the project (if needed).
#   2. Discovers MOKP benchmark instances from the FetchContent cache.
#   3. Runs example_mokp on each instance, collecting pass/fail results.
#      Instances are grouped by dimension (e.g. all 2D first, then 3D, etc.).
#   4. Prints a summary table at the end.
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail

# ── Defaults ─────────────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
BUILD_DIR="${PROJECT_DIR}/build"
MAX_K=3
QUICK=false
TYPES=("random")           # instance families to test
DIMS=("2D" "3D")           # objective dimensions to test
MAX_ITEMS=120               # skip instances with n > this (for speed)
# MAX_INSTANCES_PER_GROUP=1  # max instances per (type, dim, n) group
TIMEOUT=3600                # per-instance timeout in seconds

# ── Parse arguments ──────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
  case "$1" in
    --max-k)       MAX_K="$2"; shift 2 ;;
    --build-dir)   BUILD_DIR="$2"; shift 2 ;;
    --max-items)   MAX_ITEMS="$2"; shift 2 ;;
    --dims)        IFS=',' read -ra DIMS <<< "$2"; shift 2 ;;
    --types)       IFS=',' read -ra TYPES <<< "$2"; shift 2 ;;
    --timeout)     TIMEOUT="$2"; shift 2 ;;
    --quick)       QUICK=true; MAX_K=2; MAX_ITEMS=25; MAX_INSTANCES_PER_GROUP=1; shift ;;
    --all)         MAX_ITEMS=9999; MAX_INSTANCES_PER_GROUP=999; DIMS=("2D" "3D" "4D"); shift ;;
    -h|--help)
      echo "Usage: $0 [OPTIONS]"
      echo ""
      echo "Options:"
      echo "  --max-k N          Max reference-point set size (default: 3)"
      echo "  --build-dir DIR    Build directory (default: build)"
      echo "  --max-items N      Skip instances with n > N (default: 50)"
      echo "  --dims D1,D2,...   Dimensions to test (default: 2D,3D)"
      echo "  --types T1,T2,... Instance types to test (default: random)"
      echo "  --timeout SEC      Per-instance timeout (default: 300)"
      echo "  --quick            Quick smoke test (k<=2, n<=25, 1 inst/group)"
      echo "  --all              Run all instances (no limits)"
      echo "  -h, --help         Show this help"
      exit 0
      ;;
    *) echo "Unknown option: $1"; exit 1 ;;
  esac
done

if [[ ! "$MAX_K" =~ ^[0-9]+$ || ! "$MAX_K" =~ [1-9] ]]; then
  echo "Error: max-k must be a positive integer."
  exit 1
fi

# ── 1. Build ─────────────────────────────────────────────────────────────────
EXECUTABLE="${BUILD_DIR}/examples/example_mokp"

if [[ ! -x "$EXECUTABLE" ]]; then
  echo "Building project..."
  cmake -S "$PROJECT_DIR" -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Release 2>&1 | tail -5
  cmake --build "$BUILD_DIR" -j"$(nproc)" 2>&1 | tail -5
  echo ""
fi

if [[ ! -x "$EXECUTABLE" ]]; then
  echo "Error: executable not found at ${EXECUTABLE}"
  exit 1
fi

# ── 2. Discover instances (grouped by dimension) ─────────────────────────────
INSTANCES_DIR="${BUILD_DIR}/_deps/mobkp_instances-src/instances"
if [[ ! -d "$INSTANCES_DIR" ]]; then
  echo "Error: instances not found at ${INSTANCES_DIR}"
  echo "Run cmake configure first to fetch benchmark instances."
  exit 1
fi

# Build the instance list dimension-by-dimension so that all instances of one
# dim are consecutive before moving to the next dim.
declare -a INSTANCE_LIST=()
declare -a INSTANCE_DIM=()   # parallel array: dimension label per instance

for dim in "${DIMS[@]}"; do
  for type in "${TYPES[@]}"; do
    dir="${INSTANCES_DIR}/${type}/${dim}"
    [[ -d "$dir" ]] || continue

    count=0
    while IFS= read -r f; do
      [[ -f "$f" ]] || continue
      basename_f="$(basename "$f")"
      n="${basename_f%%_*}"
      if [[ "$n" -le "$MAX_ITEMS" ]]; then
        INSTANCE_LIST+=("$f")
        INSTANCE_DIM+=("$dim")
        count=$((count + 1))
        # [[ $count -ge $MAX_INSTANCES_PER_GROUP ]] && break
      fi
    done < <(ls "$dir"/*.in 2>/dev/null | sort -t/ -k1 -V)
  done
done

NUM_INSTANCES=${#INSTANCE_LIST[@]}
if [[ $NUM_INSTANCES -eq 0 ]]; then
  echo "No instances found matching criteria."
  exit 1
fi

echo "HS-MOCO Strategy Test"
echo "──────────────────────────────────────────────────────────────────────"
echo "  Instances  : ${NUM_INSTANCES}"
echo "  Max k      : ${MAX_K}"
echo "  Max items  : ${MAX_ITEMS}"
echo "  Dimensions : ${DIMS[*]}"
echo "  Timeout    : ${TIMEOUT}s"
echo "──────────────────────────────────────────────────────────────────────"
echo ""

# ── 3. Run tests ─────────────────────────────────────────────────────────────
TOTAL_PASS=0
TOTAL_FAIL=0
TOTAL_TIMEOUT=0
TOTAL_ERROR=0

declare -a RESULTS=()

# Set up log file
LOG_DIR="${BUILD_DIR}/test_results"
mkdir -p "$LOG_DIR"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOG_FILE="${LOG_DIR}/mokp_test_${TIMESTAMP}.log"

# Write header to log file
{
  echo "HS-MOCO Strategy Test Results"
  echo "Date: $(date)"
  echo "──────────────────────────────────────────────────────────────────────"
  echo "  Instances  : ${NUM_INSTANCES}"
  echo "  Max k      : ${MAX_K}"
  echo "  Max items  : ${MAX_ITEMS}"
  echo "  Dimensions : ${DIMS[*]}"
  echo "  Timeout    : ${TIMEOUT}s"
  echo "──────────────────────────────────────────────────────────────────────"
  echo ""
} > "$LOG_FILE"

echo "Log file: ${LOG_FILE}"
echo ""

CURRENT_DIM=""

for ((i=0; i<NUM_INSTANCES; i++)); do
  instance="${INSTANCE_LIST[$i]}"
  dim="${INSTANCE_DIM[$i]}"
  rel_path="${instance#"$INSTANCES_DIR"/}"

  # Print a dimension header when switching to a new dimension group
  if [[ "$dim" != "$CURRENT_DIM" ]]; then
    CURRENT_DIM="$dim"
    echo ""
    echo "── Dimension: ${dim} ──────────────────────────────────────────────────"
    {
      echo ""
      echo "── Dimension: ${dim} ──────────────────────────────────────────────────"
    } >> "$LOG_FILE"
  fi

  printf "[%d/%d] %-40s " "$((i+1))" "$NUM_INSTANCES" "$rel_path"

  # Run with timeout, capture output and exit code
  set +e
  output=$(timeout "$TIMEOUT" "$EXECUTABLE" "$instance" "$MAX_K" 2>&1)
  exit_code=$?
  set -e

  if [[ $exit_code -eq 124 ]]; then
    # Timeout
    printf "TIMEOUT\n"
    TOTAL_TIMEOUT=$((TOTAL_TIMEOUT + 1))
    RESULTS+=("TIMEOUT  $rel_path")
    {
      echo "[TIMEOUT] $rel_path"
      echo "--- Full output ---"
      echo "$output"
      echo "--- End output ---"
      echo ""
    } >> "$LOG_FILE"
  elif [[ $exit_code -ne 0 ]]; then
    # Error (GRBException, crash, etc.)
    error_line=$(echo "$output" | grep -iE 'GRBException|Error|error' | head -1 || true)
    printf "ERROR (exit %d) %s\n" "$exit_code" "$error_line"
    TOTAL_ERROR=$((TOTAL_ERROR + 1))
    RESULTS+=("ERROR    $rel_path  $error_line")
    {
      echo "[ERROR] $rel_path (exit $exit_code) $error_line"
      echo "--- Full output ---"
      echo "$output"
      echo "--- End output ---"
      echo ""
    } >> "$LOG_FILE"
  else
    # Parse summary line: "Summary: X/Y passed[, Z FAILED]"
    summary=$(printf '%s\n' "$output" | grep '^Summary:' || true)
    summary=${summary#Summary: }
    passed=0
    total=0
    failed=0
    if [[ "$summary" =~ ^([0-9]+)/([0-9]+)\ passed(,\ ([0-9]+)\ FAILED)?$ ]]; then
      passed=${BASH_REMATCH[1]}
      total=${BASH_REMATCH[2]}
      failed=${BASH_REMATCH[4]:-0}
    fi

    if [[ "$total" -gt 0 && "$passed" -eq "$total" && "$failed" -eq 0 ]]; then
      printf "PASS (%s/%s)\n" "$passed" "$total"
      TOTAL_PASS=$((TOTAL_PASS + 1))
      RESULTS+=("PASS     $rel_path  ($passed/$total)")
      {
        echo "[PASS] $rel_path ($passed/$total)"
        echo "--- Full output ---"
        echo "$output"
        echo "--- End output ---"
        echo ""
      } >> "$LOG_FILE"
    else
      printf "FAIL (%s)\n" "${summary:-missing summary}"
      TOTAL_FAIL=$((TOTAL_FAIL + 1))
      RESULTS+=("FAIL     $rel_path  (${summary:-missing summary})")
      {
        echo "[FAIL] $rel_path (${summary:-missing summary})"
        echo "--- Full output ---"
        echo "$output"
        echo "--- End output ---"
        echo ""
      } >> "$LOG_FILE"
    fi
  fi
done

# ── 4. Summary ───────────────────────────────────────────────────────────────
echo ""
echo "══════════════════════════════════════════════════════════════════════"
echo "Results Summary"
echo "──────────────────────────────────────────────────────────────────────"

# Also write summary to log
{
  echo ""
  echo "══════════════════════════════════════════════════════════════════════"
  echo "Results Summary"
  echo "──────────────────────────────────────────────────────────────────────"
} >> "$LOG_FILE"

for r in "${RESULTS[@]}"; do
  echo "  $r"
  echo "  $r" >> "$LOG_FILE"
done

SUMMARY_LINE=$(printf "  Pass: %d  Fail: %d  Timeout: %d  Error: %d  Total: %d" \
  "$TOTAL_PASS" "$TOTAL_FAIL" "$TOTAL_TIMEOUT" "$TOTAL_ERROR" "$NUM_INSTANCES")

echo "──────────────────────────────────────────────────────────────────────"
echo "$SUMMARY_LINE"
echo "══════════════════════════════════════════════════════════════════════"
echo ""
echo "Full results saved to: ${LOG_FILE}"

{
  echo "──────────────────────────────────────────────────────────────────────"
  echo "$SUMMARY_LINE"
  echo "══════════════════════════════════════════════════════════════════════"
} >> "$LOG_FILE"

# Exit with failure if anything went wrong
[[ $TOTAL_FAIL -eq 0 && $TOTAL_TIMEOUT -eq 0 && $TOTAL_ERROR -eq 0 ]] && exit 0 || exit 1
