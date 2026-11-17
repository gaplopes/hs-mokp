#ifndef HS_MOCO_CORE_TYPES_HPP
#define HS_MOCO_CORE_TYPES_HPP

/**
 * @file core/types.hpp
 * @brief Core vocabulary types for the MOKP hypervolume scalarization (HS) library.
 *
 * Defines solver options, result containers, and status codes shared across
 * all problem implementations and optimization strategies.
 */

#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <optional>
#include <stdexcept>
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

  /// @brief Explicit Gurobi overrides, applied during solver construction.
  std::map<std::string, std::string> gurobi_params;
};

// ─────────────────────────── Solve Result ────────────────────────────────────

/**
 * @brief Aggregated result returned from every MOKPSolver::solve() call.
 */
struct SolveResult {
  SolveStatus            status;          ///< Terminal solver status.
  double                 objective_value; ///< Approximate decoded HS value (or bound/surrogate).
  int64_t                node_count;      ///< B&B nodes explored.
  std::vector<int64_t>   solution_point;  ///< Objective-space image f(x*).
  std::vector<double>    decision_values; ///< Decision variable values x* (optional).
  std::optional<int64_t> recomputed_hs;   ///< HS of solution_point, only when requested.

  SolveResult()
      : status(SolveStatus::Error),
        objective_value(0),
        node_count(0) {}

  /**
   * @brief Round the stored objective to a bounded signed binary integer.
   *
   * Halfway cases round away from zero; the stored objective is unchanged.
   * @throws std::range_error if the value is non-finite or rounds out of range.
   */
  template<class Integer>
  [[nodiscard]] Integer objective_as() const {
    using Limits = std::numeric_limits<Integer>;
    static_assert(Limits::is_specialized && Limits::is_integer &&
                      Limits::is_signed && Limits::is_bounded && Limits::radix == 2,
                  "objective_as requires a bounded signed binary integer type");
    static_assert(Limits::lowest() < -Limits::max() &&
                      Limits::digits < std::numeric_limits<double>::max_exponent,
                  "objective_as requires a two's-complement range within double's exponent range");
    if (!std::isfinite(objective_value))
      throw std::range_error("objective value is not finite");
    const double rounded = std::round(objective_value);
    // An integer maximum can round up in double; use an exclusive power-of-two bound.
    const double bound = std::ldexp(1.0, Limits::digits);
    if (rounded < -bound || rounded >= bound)
      throw std::range_error("objective value is outside the target integer range");
    return static_cast<Integer>(rounded);
  }
};

}  // namespace hs_moco

#endif  // HS_MOCO_CORE_TYPES_HPP
