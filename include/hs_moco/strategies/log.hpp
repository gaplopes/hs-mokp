#ifndef HS_MOCO_LOG_STRATEGY_HPP
#define HS_MOCO_LOG_STRATEGY_HPP

/**
 * @file strategies/log.hpp
 * @brief LogStrategy — log-sum-exp reformulation of the HS objective.
 *
 * ═══════════════════════════════════════════════════════════════════════════
 *  1. MOTIVATION
 * ═══════════════════════════════════════════════════════════════════════════
 *
 * The HS product Π_d (f_d(x) − r_d) can span many orders of magnitude,
 * causing numerical difficulties for solvers working in the original
 * (product) space.  By taking the logarithm we convert products into sums:
 *
 *   log Π_d (f_d − r_d)  =  Σ_d log(f_d − r_d)
 *
 * This dramatically compresses the objective range and improves solver
 * numerics, at the cost of introducing non-linear general constraints
 * (log, exp, max) that Gurobi handles via piecewise-linear approximation
 * (FuncNonlinear=1).
 *
 * ═══════════════════════════════════════════════════════════════════════════
 *  2. SINGLE REFERENCE POINT (k = 1)
 * ═══════════════════════════════════════════════════════════════════════════
 *
 * For a single reference point r, the formulation is straightforward:
 *
 *   maximize  Σ_{d=0}^{M-1} log(f_d(x) − r_d)
 *
 * Auxiliary variables:
 *   z_d = f_d(x) − r_d    (continuous, lb = DOMINATION_EPS = 0.5)
 *   y_d = log(z_d)         (via Gurobi GenConstrLog)
 *
 * The z lower bound of 0.5 implicitly enforces strict domination:
 * since f_d(x) is integer-valued, z_d ≥ 0.5 implies f_d(x) ≥ r_d + 1.
 *
 * Model size: 2M variables, M linear constraints, M general constraints.
 *
 * ═══════════════════════════════════════════════════════════════════════════
 *  3. MULTIPLE REFERENCE POINTS (k ≥ 2) — INCLUSION–EXCLUSION IN LOG SPACE
 * ═══════════════════════════════════════════════════════════════════════════
 *
 * The inclusion–exclusion HS in the original space is:
 *
 *   HS = Σ_{S⊆R, S≠∅} (−1)^{|S|+1} Π_d (f_d − max_S r_d)
 *      = (sum of positive terms) − (sum of negative terms)
 *
 * With reduce_corners=true, first aggregate subset signs by component-wise
 * maximum reference corner c. Otherwise, keep every nonempty subset with
 * weight +1 or -1.
 * In log space, each retained term becomes S_c = log|w_c| + Σ_d log(f_d − c_d).
 * The full HS is:
 *
 *   HS = Σ_{w_c > 0} exp(S_c) − Σ_{w_c < 0} exp(S_c)
 *
 * We want to maximize log(HS), which becomes:
 *
 *   log(HS) = log( Σ_pos exp(S_i) − Σ_neg exp(S_j) )
 *
 * This is computed via the numerically stable identity:
 *
 *   log(A − B) = log(A) + log(1 − B/A)
 *              = log(A) + log(1 − exp(log(B) − log(A)))
 *
 * where A = Σ_pos exp(S_i) and B = Σ_neg exp(S_j), with A > B > 0.
 *
 * ═══════════════════════════════════════════════════════════════════════════
 *  4. LOG-SUM-EXP (LSE) TRICK FOR NUMERICAL STABILITY
 * ═══════════════════════════════════════════════════════════════════════════
 *
 * Computing log(Σ exp(s_i)) directly risks overflow/underflow.  The
 * standard max-shift trick avoids this:
 *
 *   LSE(s_0, …, s_n) = max(s) + log( Σ exp(s_i − max(s)) )
 *
 * Since s_i − max(s) ≤ 0, all exp arguments are non-positive, so
 * exp values lie in (0, 1] — no overflow possible.
 *
 * Implementation via Gurobi general constraints:
 *   1. max_var = max(s_0, …, s_n)        (GenConstrMax)
 *   2. diff_i = s_i − max_var             (linear constraint)
 *   3. e_i   = exp(diff_i)  ∈ (0, 1]     (GenConstrExp)
 *   4. sumExp = Σ e_i        ∈ [1, n+1]  (linear constraint)
 *   5. logSum = log(sumExp)               (GenConstrLog)
 *   6. LSE   = max_var + logSum           (linear constraint)
 *
 * Model cost per LSE(n terms): 2n+4 variables, n+2 linear constraints,
 *                               n+2 general constraints.
 * A singleton LSE returns its input directly, adding no model objects.
 *
 * ═══════════════════════════════════════════════════════════════════════════
 *  5. FINAL ASSEMBLY (k ≥ 2)
 * ═══════════════════════════════════════════════════════════════════════════
 *
 * Given Lpos = LSE(positive S terms) and Lneg = LSE(negative S terms):
 *
 *   1. Enforce Lpos ≥ Lneg + ε            (feasibility guard)
 *   2. B_minus_A = Lneg − Lpos  ≤ −ε      (continuous var)
 *   3. exp_bma   = exp(B_minus_A)  ∈ (0,1) (GenConstrExp)
 *   4. one_minus = 1 − exp_bma    ∈ (0,1] (linear)
 *   5. log_om    = log(one_minus)          (GenConstrLog)
 *   6. Maximize Lpos + log_om directly     (linear objective)
 *
 * No resultant variable or defining equality is needed.
 * decode_objective() recovers HS by exponentiating the raw objective.
 *
 * ═══════════════════════════════════════════════════════════════════════════
 *  6. MODEL SIZE SUMMARY
 * ═══════════════════════════════════════════════════════════════════════════
 *
 *   k = 1:  2M vars, M linear constrs, M general constrs.
 *
 *   With reduce_corners=true:
 *   C = 1 retained corner: same formulation as k = 1.
 *   C > 1: G shared coordinate-gap logs (each 2 vars, 1 lin, 1 gen), G <= M k
 *           + C sum-of-logs bindings (each 1 var, 1 lin)
 *           + up to 2 LSE blocks (singleton lists add no objects)
 *           + 4 final-assembly vars, 3 linear constrs, 2 general constrs
 *
 *   Without corner reduction, C = 2^k - 1 submitted subset terms; the same
 *   shared gap/log construction and LSE assembly apply.
 *
 *   General constraints scale as O(M k + C).
 *
 *   Adaptive FeasibilityTol: 1e-9 for k=1 (steep log gradients),
 *                            1e-6 for k>1 (noise accumulation).
 */

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include "gurobi_c++.h"

#include "hs_moco/core/references.hpp"
#include "hs_moco/core/strategy.hpp"
#include "hs_moco/core/types.hpp"

namespace hs_moco {

/**
 * @class LogStrategy
 * @brief Reformulates the HS objective using logarithm transformation and log-sum-exp.
 */
class LogStrategy final : public HSStrategy {
 public:
  explicit LogStrategy(bool reduce_corners = false) : reduce_corners_(reduce_corners) {}

  [[nodiscard]] std::string name() const override { return "Log"; }

  void print_statistics(const GRBModel& model,
                        std::ostream& os = std::cout) const override {
    os << "[Log] Strategy statistics\n";
    print_model_statistics(model, os);
  }

  /**
   * @brief The raw objective is log(HS); decode by exponentiating.
   */
  [[nodiscard]] int64_t decode_objective(double raw_value) const noexcept override {
    return static_cast<int64_t>(std::llround(std::exp(raw_value)));
  }

  /**
   * @brief Build the log-space HS objective for the given reference points.
   *
   * Requires FuncNonlinear=1 for log/exp general constraints.
   *
   * Dispatches to:
   *  - build_single_ref() for k=1: direct Σ log(f_d − r_d).
   *  - build_multi_ref()  for k≥2: inclusion–exclusion via LSE.
   */
  void build_objective(GRBModel& model,
                       const StrategyContext& ctx,
                       const std::vector<std::vector<int64_t>>& refs,
                       AuxiliaryState& aux) override {
    model.set(GRB_IntParam_FuncNonlinear, 1);
    const int M = ctx.M;
    const int k = static_cast<int>(refs.size());

    // Adaptive FeasibilityTol:
    // k=1: High precision (1e-9) for steep log gradients.
    // k>1: Standard precision (1e-6) to prevent noise accumulation.
    if (k == 1) {
      model.set(GRB_DoubleParam_FeasibilityTol, 1e-9);
    } else {
      model.set(GRB_DoubleParam_FeasibilityTol, 1e-6);
    }

    if (k == 1) {
      build_single_ref(model, ctx, refs[0], aux);
    } else {
      build_multi_ref(model, ctx, refs, aux);
    }
  }

 private:
  bool reduce_corners_ = false;

  static constexpr double EPS            = 1e-12;
  static constexpr double DOMINATION_EPS = 0.5;

  // ── k=1: Direct sum of logs ────────────────────────────────────────────

  void build_single_ref(GRBModel& model,
                        const StrategyContext& ctx,
                        const std::vector<int64_t>& ref,
                        AuxiliaryState& aux) const {
    const int M = ctx.M;
    GRBLinExpr S = 0;

    for (int i = 0; i < M; ++i) {
      // z_i = f_i(x) − r_i   (must be > 0 for log)
      GRBVar z = aux.addVar(model, DOMINATION_EPS, GRB_INFINITY, 0.0,
                            GRB_CONTINUOUS, "z_" + std::to_string(i));
      aux.addConstr(model, z == ctx.obj_vars[i] - ref[i],
                    "zdef_" + std::to_string(i));

      // y_i = log(z_i)
      GRBVar y = aux.addVar(model, -GRB_INFINITY, GRB_INFINITY, 0.0,
                            GRB_CONTINUOUS, "y_" + std::to_string(i));
      aux.addGenConstrLog(model, z, y, "log_" + std::to_string(i));

      S += y;
    }

    model.setObjective(S, GRB_MAXIMIZE);
  }

  // ── k≥2: Inclusion-exclusion via LSE ───────────────────────────────────

  void build_multi_ref(GRBModel& model,
                       const StrategyContext& ctx,
                       const std::vector<std::vector<int64_t>>& refs,
                       AuxiliaryState& aux) const {
    assert(refs.size() >= 2);

    std::vector<GRBVar> positives;  // positive aggregate weights
    std::vector<GRBVar> negatives;  // negative aggregate weights
    // Share gaps and logarithms within this solve; auxiliary handles expire at cleanup.
    std::vector<std::map<int64_t, GRBVar>> log_vars(ctx.M);
    auto add_term = [&](const std::vector<int64_t>& corner, int64_t weight, int index) {
      GRBVar S_mask = build_sum_of_logs(model, ctx.obj_vars, corner, weight, log_vars,
                                        "S_" + std::to_string(index), aux);
      if (weight > 0)
        positives.push_back(S_mask);
      else
        negatives.push_back(S_mask);
    };

    if (reduce_corners_) {
      const auto corners = detail::inclusion_exclusion_corners(refs);
      if (corners.size() == 1) {
        build_single_ref(model, ctx, corners.begin()->first, aux);
        return;
      }
      // Build S = log(|weight|) + Σ_i log(f_i − corner_i) for each retained corner.
      int index = 0;
      for (const auto& [corner, weight] : corners)
        add_term(corner, weight, index++);
    } else {
      const int k = static_cast<int>(refs.size());
      const int total_masks = (1 << k);
      for (int mask = 1; mask < total_masks; ++mask) {
        // Compute component-wise max reference for this subset.
        std::vector<int64_t> max_ref(ctx.M, std::numeric_limits<int64_t>::min());
        int bits = 0;
        for (int j = 0; j < k; ++j) {
          if (mask & (1 << j)) {
            ++bits;
            for (int i = 0; i < ctx.M; ++i)
              max_ref[i] = std::max(max_ref[i], refs[j][i]);
          }
        }
        int sign = (bits % 2 == 1) ? +1 : -1;
        add_term(max_ref, sign, mask);
      }
    }

    // Compute: result = log( Σ_pos exp(S_i) − Σ_neg exp(S_j) )
    //        = Lpos + log(1 − exp(Lneg − Lpos))

    // Compute Lpos = log(Σ_pos exp(S_i))
    GRBVar Lpos = log_sum_exp(model, positives, "pos_lse", aux);
    // Compute Lneg = log(Σ_neg exp(S_j))
    GRBVar Lneg = log_sum_exp(model, negatives, "neg_lse", aux);

    // Enforce Lpos > Lneg for numerical stability
    aux.addConstr(model, Lpos >= Lneg + EPS, "ensure_Lpos_gt_Lneg");

    // B_minus_A = Lneg − Lpos  (≤ −EPS)
    GRBVar B_minus_A = aux.addVar(model, -GRB_INFINITY, 0.0, 0.0,
                                  GRB_CONTINUOUS, "B_minus_A");
    aux.addConstr(model, B_minus_A == Lneg - Lpos, "def_B_minus_A");

    // exp(B_minus_A) ∈ (0, 1)
    GRBVar exp_bma = aux.addVar(model, 0.0, 1.0, 0.0, GRB_CONTINUOUS, "exp_bma");
    aux.addGenConstrExp(model, B_minus_A, exp_bma, "exp_bma_constr");

    // one_minus = 1 − exp(B_minus_A) ∈ (0, 1]
    GRBVar one_minus = aux.addVar(model, 0.0, 1.0, 0.0, GRB_CONTINUOUS, "one_minus");
    aux.addConstr(model, one_minus == 1.0 - exp_bma, "def_one_minus");

    // log_one_minus = log(one_minus)
    GRBVar log_om = aux.addVar(model, -GRB_INFINITY, GRB_INFINITY, 0.0,
                               GRB_CONTINUOUS, "log_one_minus");
    aux.addGenConstrLog(model, one_minus, log_om, "log_one_minus_constr");

    // result = Lpos + log(1 − exp(Lneg − Lpos))
    model.setObjective(GRBLinExpr(Lpos) + log_om, GRB_MAXIMIZE);
  }

  // ── Helper: S = log|weight| + Σ_i log(f_i − ref_i) ─────────────────────
  //
  // Builds the log magnitude of one aggregated inclusion–exclusion term.
  //
  // For each objective d:
  //   z_d = f_d(x) − ref_d   (lb = DOMINATION_EPS → strict domination)
  //   y_d = log(z_d)         (Gurobi GenConstrLog)
  //   S   = log|weight| + Σ_d y_d
  // Gap/log pairs are reused whenever the objective and reference coordinate match.

  static GRBVar build_sum_of_logs(GRBModel& model,
                                  const std::vector<GRBVar>& obj_vars,
                                  const std::vector<int64_t>& ref,
                                  int64_t weight,
                                  std::vector<std::map<int64_t, GRBVar>>& log_vars,
                                  const std::string& tag,
                                  AuxiliaryState& aux) {
    const int M = static_cast<int>(obj_vars.size());
    GRBLinExpr sum_y = std::log(std::abs(static_cast<double>(weight)));

    for (int i = 0; i < M; ++i) {
      auto found = log_vars[i].find(ref[i]);
      if (found == log_vars[i].end()) {
        GRBVar z = aux.addVar(model, DOMINATION_EPS, GRB_INFINITY, 0.0,
                              GRB_CONTINUOUS, tag + "_z_" + std::to_string(i));
        aux.addConstr(model, z == obj_vars[i] - ref[i],
                      tag + "_zdef_" + std::to_string(i));

        GRBVar y = aux.addVar(model, -GRB_INFINITY, GRB_INFINITY, 0.0,
                              GRB_CONTINUOUS, tag + "_y_" + std::to_string(i));
        aux.addGenConstrLog(model, z, y, tag + "_log_" + std::to_string(i));
        found = log_vars[i].emplace(ref[i], y).first;
      }

      sum_y += found->second;
    }

    GRBVar S = aux.addVar(model, -GRB_INFINITY, GRB_INFINITY, 0.0,
                          GRB_CONTINUOUS, tag + "_S");
    aux.addConstr(model, S == sum_y, tag + "_Sdef");
    return S;
  }

  // ── Helper: LSE = log(Σ exp(v_i)) via max-shift trick ─────────────────
  //
  // Computes LSE(v_0, …, v_{n-1}) = log(Σ_i exp(v_i)) in a numerically
  // stable way by shifting all arguments by the maximum:
  //
  //   LSE = max(v) + log( Σ exp(v_i − max(v)) )
  //
  // Since v_i − max(v) ≤ 0, all exp values are in (0, 1], avoiding overflow.
  // At least one diff is 0, so sumExp ≥ 1 and logSum ≥ 0.

  static GRBVar log_sum_exp(GRBModel& model,
                            const std::vector<GRBVar>& vars,
                            const std::string& tag,
                            AuxiliaryState& aux) {
    if (vars.empty())
      throw std::runtime_error("log_sum_exp: empty var list");
    if (vars.size() == 1) return vars.front();

    // maxVar = max(vars)
    GRBVar max_var = aux.addVar(model, -GRB_INFINITY, GRB_INFINITY, 0.0,
                                GRB_CONTINUOUS, tag + "_max");
    // Gurobi's addGenConstrMax expects a non-const pointer to GRBVar array
    std::vector<GRBVar> vars_copy(vars);
    aux.addGenConstrMax(model, max_var, vars_copy.data(),
                        static_cast<int>(vars_copy.size()), -GRB_INFINITY,
                        tag + "_maxConstr");

    // For each var: diff = var − maxVar; e = exp(diff) ∈ (0, 1]
    std::vector<GRBVar> exp_vars;
    exp_vars.reserve(vars.size());
    for (size_t i = 0; i < vars.size(); ++i) {
      GRBVar diff = aux.addVar(model, -GRB_INFINITY, GRB_INFINITY, 0.0,
                               GRB_CONTINUOUS,
                               tag + "_diff_" + std::to_string(i));
      aux.addConstr(model, diff == vars[i] - max_var,
                    tag + "_diff_def_" + std::to_string(i));

      GRBVar e = aux.addVar(model, 0.0, 1.0, 0.0, GRB_CONTINUOUS,
                            tag + "_exp_" + std::to_string(i));
      aux.addGenConstrExp(model, diff, e,
                          tag + "_expConstr_" + std::to_string(i));
      exp_vars.push_back(e);
    }

    // sumExp ∈ [1, |vars|]  (at least one exp(0)=1)
    GRBVar sum_exp = aux.addVar(model, 1.0, static_cast<double>(vars.size()),
                                0.0, GRB_CONTINUOUS, tag + "_sumExp");
    GRBLinExpr tmp = 0.0;
    for (auto& ev : exp_vars) tmp += ev;
    aux.addConstr(model, sum_exp == tmp, tag + "_sumExp_def");

    // logSum = log(sumExp)
    GRBVar log_sum = aux.addVar(model, -GRB_INFINITY, GRB_INFINITY, 0.0,
                                GRB_CONTINUOUS, tag + "_logSum");
    aux.addGenConstrLog(model, sum_exp, log_sum, tag + "_logSumConstr");

    // result = maxVar + logSum
    GRBVar result = aux.addVar(model, -GRB_INFINITY, GRB_INFINITY, 0.0,
                               GRB_CONTINUOUS, tag + "_LSE");
    aux.addConstr(model, result == max_var + log_sum, tag + "_LSE_def");
    return result;
  }
};

}  // namespace hs_moco

#endif  // HS_MOCO_LOG_STRATEGY_HPP
