#ifndef HS_MOCO_CORE_TYPES_HPP
#define HS_MOCO_CORE_TYPES_HPP

/**
 * @file core/types.hpp
 * @brief Core vocabulary types for the MOKP hypervolume scalarization (HS) library.
 *
 * Defines solver options, result containers, and status codes shared across
 * all problem implementations and optimization strategies.
 */

#include <cstdint>
#include <string>
#include <vector>

namespace hs_moco {

// ─────────────────────────── Solve Status ────────────────────────────────────

/**
 * @brief Terminal status of a single solver invocation.
 */
enum class SolveStatus : uint8_t {
  Optimal,           ///< Proven global optimum found.
  TimeLimitReached,  ///< Time budget exhausted; incumbent may be available.
  InfeasibleOrUnbounded,        ///< No feasible solution exists for the given reference points.
  NodeLimitReached,  ///< B&B node limit hit; bound is returned.
  Error,             ///< Solver or model error.
};

/**
 * @brief Human-readable conversion of SolveStatus.
 */
[[nodiscard]] inline std::string to_string(SolveStatus s) {
  switch (s) {
    case SolveStatus::Optimal:                return "Optimal";
    case SolveStatus::TimeLimitReached:       return "TimeLimitReached";
    case SolveStatus::InfeasibleOrUnbounded:  return "InfeasibleOrUnbounded";
    case SolveStatus::NodeLimitReached:       return "NodeLimitReached";
    case SolveStatus::Error:                  return "Error";
  }
  return "Unknown";
}

// ─────────────────────────── Solver Options ──────────────────────────────────

/**
 * @brief Tunable parameters forwarded to the underlying MIP solver.
 *
 * All fields have sensible defaults; users only override what they need.
 */
struct SolverOptions {
  double time_limit     = 3600.0;   ///< Wall-clock time limit in seconds.
  double mip_gap        = 1e-6;    ///< Relative MIP gap tolerance.
  int    threads        = 4;       ///< Maximum solver threads.
  int    numeric_focus  = 3;       ///< 0..3 — higher = more careful numerics.
  bool   console_log    = false;   ///< Print solver log to stdout.
  double feasibility_tol = 1e-9;   ///< Primal feasibility tolerance.

  /// @brief Strategy-specific flag: enable Gurobi NonConvex=2 (required by MIQPStrategy).
  bool   allow_nonconvex = false;

  /// @brief Strategy-specific flag: enable Gurobi FuncNonlinear=1 (required by LogStrategy).
  bool   func_nonlinear  = false;
};

// ─────────────────────────── Solve Result ────────────────────────────────────

/**
 * @brief Aggregated result returned from every MOKPSolver::solve() call.
 */
struct SolveResult {
  SolveStatus            status;          ///< Terminal solver status.
  int64_t                objective_value; ///< Best HS value found (or surrogate).
  int64_t                node_count;      ///< B&B nodes explored.
  std::vector<int64_t>   solution_point;  ///< Objective-space image f(x*).
  std::vector<double>    decision_values; ///< Decision variable values x* (optional).

  SolveResult()
      : status(SolveStatus::Error),
        objective_value(0),
        node_count(0) {}
};

}  // namespace hs_moco

#endif  // HS_MOCO_CORE_TYPES_HPP
