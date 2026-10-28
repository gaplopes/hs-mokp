#ifndef HS_MOCO_CORE_SOLVER_HPP
#define HS_MOCO_CORE_SOLVER_HPP

/**
 * @file core/solver.hpp
 * @brief MOKPSolver — the central orchestrator for MOKP hypervolume scalarization.
 *
 * Usage synopsis
 * ──────────────
 * @code
 *   auto problem  = std::make_unique<hs_moco::MOKProblem>(items, W);
 *   auto strategy = std::make_unique<hs_moco::MIQPStrategy>();
 *   hs_moco::SolverOptions opts;
 *   opts.time_limit = 300.0;
 *
 *   hs_moco::MOKPSolver solver(std::move(problem), std::move(strategy), opts);
 *
 *   // Can be called repeatedly with different reference points:
 *   auto r1 = solver.solve({{0, 0, 0}});
 *   auto r2 = solver.solve({{1, 2, 3}, {4, 5, 6}});
 * @endcode
 *
 * Lifecycle
 * ─────────
 * 1. **Construction** (expensive, done once):
 *    - Creates GRBEnv + GRBModel.
 *    - Calls Problem::setup_variables() → decision vars x_j.
 *    - Calls Problem::setup_constraints() → feasibility constraints.
 *    - Calls Problem::build_objective_exprs() → f_i(x).
 *    - Binds each f_i(x) to a continuous GRBVar (obj_vars).
 *
 * 2. **solve(refs)** (repeatable):
 *    - Clears previous temporary state, retaining reusable strategy structure.
 *    - Adds domination constraints: f_i(x) > max_j(refs_j^i).
 *    - Calls Strategy::build_objective() to set up the HS formulation.
 *    - Optimizes.
 *    - Extracts SolveResult.
 *    - Retains state until cleanup before the next call.
 */

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <vector>

#include "gurobi_c++.h"
#include "hs_moco/problems/mokp.hpp"
#include "hs_moco/core/strategy.hpp"
#include "hs_moco/core/types.hpp"

namespace hs_moco {

class MOKPSolver {
 public:
  /**
   * @brief Construct a solver for the given problem and strategy.
   *
   * Performs all one-time model setup (variables, constraints, f_i(x) binding).
   * Throws GRBException on Gurobi errors.
   *
   * @param problem   Concrete MOKP problem (ownership transferred).
   * @param strategy  Concrete HS strategy (ownership transferred).
   * @param opts      Solver parameters.
   * @param var_type  GRB_BINARY (default) or GRB_CONTINUOUS for LP relaxation.
   */
  MOKPSolver(std::unique_ptr<MOKProblem> problem,
             std::unique_ptr<HSStrategy> strategy,
             SolverOptions opts = {},
             char var_type = GRB_BINARY)
      : problem_(std::move(problem)),
        strategy_(std::move(strategy)),
        opts_(opts),
        M_(0),
        N_(0) {
    if (!problem_) throw std::invalid_argument("MOKPSolver: problem must not be null");
    if (!strategy_) throw std::invalid_argument("MOKPSolver: strategy must not be null");

    M_ = problem_->num_objectives();
    N_ = problem_->num_variables();

    // Let the strategy override options it requires
    strategy_->configure_options(opts_);

    // ── 1. Environment ────────────────────────────────────────────────────
    env_ = std::make_unique<GRBEnv>(true);
    env_->set(GRB_IntParam_OutputFlag, opts_.console_log ? 1 : 0);
    env_->set(GRB_IntParam_Threads, opts_.threads);
    env_->set(GRB_DoubleParam_TimeLimit, opts_.time_limit);
    env_->set(GRB_DoubleParam_MIPGap, opts_.mip_gap);
    env_->set(GRB_IntParam_NumericFocus, opts_.numeric_focus);
    env_->set(GRB_DoubleParam_FeasibilityTol, opts_.feasibility_tol);

    for (const auto& [name, value] : opts_.gurobi_params)
      env_->set(name, value);
    env_->start();

    // ── 2. Model ──────────────────────────────────────────────────────────
    model_ = std::make_unique<GRBModel>(*env_);

    // ── 3. Decision variables ─────────────────────────────────────────────
    x_vars_ = problem_->setup_variables(*model_, var_type);
    if (static_cast<int>(x_vars_.size()) != N_)
      throw std::logic_error("setup_variables() returned wrong number of variables");

    // ── 4. Feasibility constraints ────────────────────────────────────────
    problem_->setup_constraints(*model_, x_vars_);

    // ── 5. Objective expressions f_i(x) → bound to continuous vars ───────
    auto f_exprs = problem_->build_objective_exprs(*model_, x_vars_);
    if (static_cast<int>(f_exprs.size()) != M_)
      throw std::logic_error("build_objective_exprs() returned wrong number of expressions");

    obj_vars_.resize(M_);
    for (int i = 0; i < M_; ++i) {
      obj_vars_[i] = model_->addVar(-GRB_INFINITY, GRB_INFINITY, 0.0,
                                    GRB_CONTINUOUS, "f_" + std::to_string(i));
      model_->addConstr(obj_vars_[i] == f_exprs[i], "f_def_" + std::to_string(i));
    }

    // ── 6. Cache coefficients and upper bounds for strategies ──────────────
    coeff_matrix_ = problem_->coefficient_matrix();
    objective_upper_bounds_ = problem_->objective_upper_bounds();

    model_->update();
  }

  /// Non-copyable (Gurobi handles are not copyable).
  MOKPSolver(const MOKPSolver&) = delete;
  MOKPSolver& operator=(const MOKPSolver&) = delete;

  /// Movable.
  MOKPSolver(MOKPSolver&&) noexcept = default;
  MOKPSolver& operator=(MOKPSolver&&) noexcept = default;

  // ─────────────────────────── Public API ─────────────────────────────────

  /**
   * @brief Solve the HS sub-problem for the given reference points.
   *
   * Can be called repeatedly — auxiliary state is cleaned up automatically.
   *
   * @param refs  Non-empty vector of reference points (each M-dimensional).
   *              Counts at or above std::numeric_limits<int>::digits return Error.
   * @return SolveResult with status, objective value, and solution point.
   */
  [[nodiscard]] SolveResult solve(const std::vector<std::vector<int64_t>>& refs) {
    SolveResult result;
    result.solution_point.assign(M_, 0);

    if (!validate_refs(refs)) {
      result.status = SolveStatus::Error;
      return result;
    }

    try {
      // ── 1. Clean up any previous auxiliary state ────────────────────────
      strategy_->clear_auxiliary(*model_, aux_);
      domination_aux_.clear(*model_);

      // ── 2. Domination constraints: f_i(x) > max_j(refs_j^i) ───────────
      add_domination_constraints(refs);

      // ── 3. Build the HS objective via the strategy ─────────────────────
      StrategyContext ctx{M_, N_, obj_vars_, x_vars_, coeff_matrix_, objective_upper_bounds_};
      strategy_->build_objective(*model_, ctx, refs, aux_);

      model_->update();

      // ── 4. Optimize ────────────────────────────────────────────────────
      model_->optimize();

      // ── 5. Extract result ──────────────────────────────────────────────
      int status = model_->get(GRB_IntAttr_Status);
      result = extract_result(status);

    } catch (GRBException& e) {
      std::cerr << "[MOKPSolver] Gurobi error " << e.getErrorCode()
                << ": " << e.getMessage() << std::endl;
      result.status = SolveStatus::Error;
    } catch (std::exception& e) {
      std::cerr << "[MOKPSolver] Error: " << e.what() << std::endl;
      result.status = SolveStatus::Error;
    }

    return result;
  }

  // ─────────────────────── Accessors ──────────────────────────────────────

  /** @brief Number of objectives. */
  [[nodiscard]] int num_objectives() const noexcept { return M_; }

  /** @brief Number of decision variables. */
  [[nodiscard]] int num_variables() const noexcept { return N_; }

  /** @brief Strategy name. */
  [[nodiscard]] std::string strategy_name() const { return strategy_->name(); }

  /** @brief Read-only access to the underlying Gurobi model (for diagnostics). */
  [[nodiscard]] const GRBModel& model() const { return *model_; }

  /** @brief Print strategy/model statistics for the current model state. */
  void print_statistics(std::ostream& os = std::cout) const {
    strategy_->print_statistics(*model_, os);
  }

 private:
  // ── Owned state ─────────────────────────────────────────────────────────
  std::unique_ptr<const MOKProblem> problem_;
  std::unique_ptr<HSStrategy> strategy_;
  SolverOptions opts_;

  std::unique_ptr<GRBEnv> env_;
  std::unique_ptr<GRBModel> model_;

  int M_;  ///< Number of objectives.
  int N_;  ///< Number of decision variables.

  std::vector<GRBVar> x_vars_;    ///< Decision variables (immutable after construction).
  std::vector<GRBVar> obj_vars_;  ///< f_i(x) continuous variables (immutable).

  std::vector<std::vector<int64_t>> coeff_matrix_;  ///< Cached coefficient matrix.
  std::vector<double> objective_upper_bounds_;     ///< Cached fractional-knapsack bounds.

  AuxiliaryState aux_;             ///< Strategy-owned auxiliary state.
  AuxiliaryState domination_aux_;  ///< Domination constraint auxiliary state.

  // ── Internal helpers ────────────────────────────────────────────────────

  [[nodiscard]] bool validate_refs(const std::vector<std::vector<int64_t>>& refs) const noexcept {
    if (refs.empty() || refs.size() >= std::numeric_limits<int>::digits) {
      return false;
    }
    for (const auto& r : refs) {
      if (static_cast<int>(r.size()) != M_) {
        return false;
      }
    }
    return true;
  }

  /**
   * @brief Add constraints f_i(x) >= max_j(refs_j^i) + 1 for strict domination.
   *
   * These ensure the solution strictly dominates all reference points.
   */
  void add_domination_constraints(const std::vector<std::vector<int64_t>>& refs) {
    constexpr int64_t MIN = std::numeric_limits<int64_t>::min();
    std::vector<int64_t> max_ref(M_, MIN);
    for (const auto& r : refs) {
      for (int i = 0; i < M_; ++i) {
        max_ref[i] = std::max(max_ref[i], r[i]);
      }
    }
    for (int i = 0; i < M_; ++i) {
      domination_aux_.addConstr(*model_,
                                obj_vars_[i] >= max_ref[i] + 1,
                                "domination_" + std::to_string(i));
    }
  }

  /**
   * @brief Map Gurobi status to SolveResult and extract solution.
   */
  SolveResult extract_result(int grb_status) {
    SolveResult result;
    result.solution_point.assign(M_, 0);
    result.node_count = static_cast<int64_t>(model_->get(GRB_DoubleAttr_NodeCount));

    switch (grb_status) {
      case GRB_OPTIMAL:
        result.status = SolveStatus::Optimal;
        result.objective_value = strategy_->decode_objective(
            model_->get(GRB_DoubleAttr_ObjVal));
        for (int i = 0; i < M_; ++i) {
          result.solution_point[i] =
              static_cast<int64_t>(std::llround(obj_vars_[i].get(GRB_DoubleAttr_X)));
        }
        // Optionally extract decision values
        result.decision_values.resize(N_);
        for (int j = 0; j < N_; ++j) {
          result.decision_values[j] = x_vars_[j].get(GRB_DoubleAttr_X);
        }
        break;

      case GRB_TIME_LIMIT:
        result.status = SolveStatus::TimeLimitReached;
        // Try to get incumbent if available
        if (model_->get(GRB_IntAttr_SolCount) > 0) {
          result.objective_value = strategy_->decode_objective(
              model_->get(GRB_DoubleAttr_ObjVal));
          for (int i = 0; i < M_; ++i) {
            result.solution_point[i] =
                static_cast<int64_t>(std::llround(obj_vars_[i].get(GRB_DoubleAttr_X)));
          }
        }
        break;

      case GRB_INF_OR_UNBD:
      case GRB_INFEASIBLE:
      case GRB_UNBOUNDED:
        result.status = SolveStatus::InfeasibleOrUnbounded;
        break;

      case GRB_NODE_LIMIT:
        result.status = SolveStatus::NodeLimitReached;
        if (model_->get(GRB_IntAttr_SolCount) > 0) {
          result.objective_value = strategy_->decode_objective(
              model_->get(GRB_DoubleAttr_ObjBound));
        }
        break;

      default:
        result.status = SolveStatus::Error;
        break;
    }

    return result;
  }
};

}  // namespace hs_moco

#endif  // HS_MOCO_CORE_SOLVER_HPP
