#ifndef HS_MOCO_CORE_STRATEGY_HPP
#define HS_MOCO_CORE_STRATEGY_HPP

/**
 * @file core/strategy.hpp
 * @brief Abstract HS-strategy interface and auxiliary-state tracker.
 *
 * A strategy encodes *how* the hypervolume scalarization (HS) is formulated
 * in terms of the Gurobi model.  Three concrete strategies are provided:
 *
 *  - MIQPStrategy  : direct polynomial product via nonlinear constraints.
 *  - LogStrategy   : log-sum-exp reformulation.
 *  - MILPStrategy  : McCormick-envelope linearization.
 *
 * ──── Auxiliary-state lifecycle ────
 * Every call to build_objective() may add auxiliary variables, linear /
 * quadratic / general constraints to the model.  These are recorded in an
 * AuxiliaryState and passed to clear_auxiliary() before the next solve().
 * Strategies may retain reference-independent structure; the immutable base
 * model (decision vars + feasibility constraints) remains intact.
 */

#include <cmath>
#include <cstdint>
#include <iostream>
#include <ostream>
#include <string>
#include <vector>

#include "gurobi_c++.h"

namespace hs_moco {

// ─────────────────────── Auxiliary State Tracker ─────────────────────────────

/**
 * @brief RAII-ish tracker for temporary Gurobi objects created by a strategy.
 *
 * Records all auxiliary GRBVars, GRBConstrs, GRBQConstrs, and GRBGenConstrs
 * added during build_objective().  Calling clear() removes them from the model
 * in the correct order (constraints first, then variables).
 */
struct AuxiliaryState {
  std::vector<GRBVar>       vars;
  std::vector<GRBConstr>    constrs;
  std::vector<GRBQConstr>   qconstrs;
  std::vector<GRBGenConstr> gen_constrs;

  /** @brief True if nothing has been recorded. */
  [[nodiscard]] bool empty() const noexcept {
    return vars.empty() && constrs.empty() && qconstrs.empty() && gen_constrs.empty();
  }

  /**
   * @brief Remove every tracked object from @p model, then reset.
   *
   * Removal order: general constraints → quadratic constraints →
   * linear constraints → variables.  This avoids dangling references.
   */
  void clear(GRBModel& model) {
    for (auto& c : gen_constrs) model.remove(c);
    for (auto& c : qconstrs)    model.remove(c);
    for (auto& c : constrs)     model.remove(c);
    for (auto& v : vars)        model.remove(v);
    gen_constrs.clear();
    qconstrs.clear();
    constrs.clear();
    vars.clear();
    model.update();
  }

  // ── Convenience recorders ───────────────────────────────────────────────

  GRBVar addVar(GRBModel& model, double lb, double ub, double obj,
                char vtype, const std::string& name = "") {
    GRBVar v = model.addVar(lb, ub, obj, vtype, name);
    vars.push_back(v);
    return v;
  }

  GRBConstr addConstr(GRBModel& model, const GRBTempConstr& tc,
                      const std::string& name = "") {
    GRBConstr c = model.addConstr(tc, name);
    constrs.push_back(c);
    return c;
  }

  GRBQConstr addQConstr(GRBModel& model, const GRBTempConstr& tc,
                        const std::string& name = "") {
    GRBQConstr c = model.addQConstr(tc, name);
    qconstrs.push_back(c);
    return c;
  }

  void addGenConstrNL(GRBModel& model, GRBVar resvar, int nnodes,
                      const int* opcode, const double* data, const int* parent,
                      const std::string& name = "") {
    GRBGenConstr gc = model.addGenConstrNL(resvar, nnodes, opcode, data, parent, name);
    gen_constrs.push_back(gc);
  }

  void addGenConstrLog(GRBModel& model, GRBVar xvar, GRBVar yvar,
                       const std::string& name = "") {
    GRBGenConstr gc = model.addGenConstrLog(xvar, yvar, name);
    gen_constrs.push_back(gc);
  }

  void addGenConstrExp(GRBModel& model, GRBVar xvar, GRBVar yvar,
                       const std::string& name = "") {
    GRBGenConstr gc = model.addGenConstrExp(xvar, yvar, name);
    gen_constrs.push_back(gc);
  }

  void addGenConstrMax(GRBModel& model, GRBVar resvar, GRBVar* vars,
                       int nvars, double constant,
                       const std::string& name = "") {
    GRBGenConstr gc = model.addGenConstrMax(resvar, vars, nvars, constant, name);
    gen_constrs.push_back(gc);
  }
};

// ──────────────────────── Strategy Context ───────────────────────────────────

/**
 * @brief Read-only context passed to a strategy's build_objective().
 *
 * Bundles all the information a strategy needs without coupling it to the
 * solver's internals.
 */
struct StrategyContext {
  int                          M;           ///< Number of objectives.
  int                          N;           ///< Number of decision variables.
  const std::vector<GRBVar>&   obj_vars;    ///< f_i(x) as Gurobi continuous variables.
  const std::vector<GRBVar>&   x_vars;      ///< Decision variables x_j.

  /// Per-variable, per-objective coefficients (may be empty).
  const std::vector<std::vector<int64_t>>& coeff_matrix;

  /// Optional upper bounds for each objective (used for normalization).
  std::vector<double> objective_upper_bounds{}; ///< Optional bounds for objective normalization.
};

// ──────────────────────── Abstract Strategy ──────────────────────────────────

/**
 * @brief Interface that every HS-formulation strategy implements.
 */
class HSStrategy {
 public:
  virtual ~HSStrategy() = default;

  /**
   * @brief Human-readable name (e.g. "MIQP", "Log", "MILP").
   */
  [[nodiscard]] virtual std::string name() const = 0;

  /**
   * @brief Apply strategy-specific parameter overrides to @p opts.
   *
   * Called once during solver construction, before creating the environment
   * and model.
   */
  virtual void configure_options(struct SolverOptions& opts) const {}

  /** @brief Clear temporary state before the next solve; reusable structure may remain. */
  virtual void clear_auxiliary(GRBModel& model, AuxiliaryState& aux) {
    aux.clear(model);
  }

  /**
   * @brief Build the HS objective for the given reference points.
   *
   * Adds auxiliary variables / constraints to @p model (tracked via @p aux)
   * and sets the model objective.
   *
   * @param model  Gurobi model (base constraints already present).
   * @param ctx    Read-only strategy context.
   * @param refs   Reference points (each is an M-dimensional coordinate vector).
   * @param aux    Auxiliary state — record every added var / constr here.
   */
  virtual void build_objective(GRBModel& model,
                               const StrategyContext& ctx,
                               const std::vector<std::vector<int64_t>>& refs,
                               AuxiliaryState& aux) = 0;

  /**
   * @brief Print strategy/model statistics.
   *
   * Intended for diagnostics after build_objective() and before/after optimize().
   */
  virtual void print_statistics(const GRBModel& model,
                                std::ostream& os = std::cout) const {
    os << "[" << name() << "] Model statistics\n";
    print_model_statistics(model, os);
  }

  /**
   * @brief Post-process the raw Gurobi objective value.
   *
   * Some strategies (e.g. Log) store log(HS); this hook converts back.
   * Default: identity (returns @p raw_value as-is).
   */
  [[nodiscard]] virtual int64_t decode_objective(double raw_value) const noexcept {
    return static_cast<int64_t>(std::llround(raw_value));
  }

 protected:
  static void print_model_statistics(const GRBModel& model, std::ostream& os) {
    os << "  NumVars      : " << model.get(GRB_IntAttr_NumVars) << "\n";
    os << "  NumConstrs   : " << model.get(GRB_IntAttr_NumConstrs) << "\n";
    os << "  NumQConstrs  : " << model.get(GRB_IntAttr_NumQConstrs) << "\n";
    os << "  NumGenConstrs: " << model.get(GRB_IntAttr_NumGenConstrs) << "\n";
    int64_t total_constrs = model.get(GRB_IntAttr_NumConstrs) +
                           model.get(GRB_IntAttr_NumQConstrs) +
                           model.get(GRB_IntAttr_NumGenConstrs);
    os << "  TotalConstrs : " << total_constrs << "\n";
  }
};

}  // namespace hs_moco

#endif  // HS_MOCO_CORE_STRATEGY_HPP
