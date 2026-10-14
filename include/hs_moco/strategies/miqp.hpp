#ifndef HS_MOCO_MIQP_STRATEGY_HPP
#define HS_MOCO_MIQP_STRATEGY_HPP

/**
 * @file strategies/miqp.hpp
 * @brief MIQPStrategy — direct polynomial-product HS formulation.
 *
 * ═══════════════════════════════════════════════════════════════════════════
 *  1. MATHEMATICAL FORMULATION
 * ═══════════════════════════════════════════════════════════════════════════
 *
 * The Hypervolume Scalarization (HS) for reference points R is:
 *
 *   HS(x, R) = Σ_{S⊆R, S≠∅} (−1)^{|S|+1} Π_{d=0}^{M-1} (f_d(x) − max_{r∈S} r_d)
 *
 * This strategy models each product Π_d (f_d(x) − r_d) directly as a
 * nonlinear equality using addGenConstrNL (Gurobi 12 or newer).
 * Gurobi handles the product through spatial branch-and-bound.
 * The balanced bilinear implementation is retained below and uses NonConvex=2.
 * PrecisionOptimal normalization divides each factor by
 * norm_d = ub_d / (abs_tol · Π_i ub_i)^(1/M), using Dantzig upper bounds.
 * abs_tol = max(1e-8, 10 · FeasibilityTol),
 * so the modeled objective is abs_tol · HS;
 * decode_objective() restores the original units. Nonpositive bounds leave
 * PrecisionOptimal formulation unscaled.
 * NormalizationMethod::None is the default and leaves factors unscaled.
 * Select a normalization method from strategies/norm.hpp in the constructor.
 *
 * ═══════════════════════════════════════════════════════════════════════════
 *  2. NONLINEAR PRODUCT AND RETAINED BINARY-TREE REDUCTION
 * ═══════════════════════════════════════════════════════════════════════════
 *
 * By default, M >= 2 terms are bound to continuous variables and passed as
 * children of one GRB_OPCODE_MULTIPLY node in a nonlinear expression tree.
 * One resultant variable P and one general constraint define P = Π_d t_d.
 * Gurobi chooses the internal product formulation. Zero or one term uses
 * the retained helper's constant or linear binding, respectively.
 *
 * The retained balanced binary-tree reduction computes the product of
 * M linear terms t_0, t_1, …, t_{M-1} and minimizes the longest chain of
 * dependent quadratic constraints:
 *
 *   Level 0:  p_{01} = t_0 · t_1,   p_{23} = t_2 · t_3,   …
 *   Level 1:  p_{0123} = p_{01} · p_{23},  …
 *   …
 *   Final:    P = single product variable
 *
 * Each multiplication introduces:
 *   - 1 auxiliary continuous variable for the product result
 *   - 1 quadratic equality constraint (GRBQConstr)
 *
 * Tree depth is ⌈log₂(M)⌉, giving O(M) aux variables and QConstrs total.
 * If M is odd at any level, the unpaired term is carried to the next level.
 *
 *   Example (M = 3, terms a, b, c):
 *     Level 0: ab = a · b        (1 var, 1 QConstr)
 *     Level 1: P  = ab · c       (1 var, 1 QConstr)
 *     Total: 2 aux vars, 2 QConstrs
 *
 *   Example (M = 4, terms a, b, c, d):
 *     Level 0: ab = a · b, cd = c · d    (2 vars, 2 QConstrs)
 *     Level 1: P  = ab · cd              (1 var, 1 QConstr)
 *     Total: 3 aux vars, 3 QConstrs
 *
 * Additionally, each normalized term t_d = (f_d(x) − r_d) / norm_d is bound to an
 * auxiliary variable via a linear equality, adding M vars and M linear
 * constraints per subset.
 *
 * ═══════════════════════════════════════════════════════════════════════════
 *  3. INCLUSION–EXCLUSION ASSEMBLY
 * ═══════════════════════════════════════════════════════════════════════════
 *
 * For K reference points, there are 2^K − 1 non-empty subsets S.
 * Each subset contributes one product term P_S with sign (−1)^{|S|+1}.
 * The final objective is a linear combination of all P_S:
 *
 *   obj = Σ_S sign_S · P_S
 *
 * This is set as a GRBLinExpr over the product variables (which are
 * continuous), so the objective itself is linear — all non-linearity
 * resides in the nonlinear constraints that define the products.
 *
 * ═══════════════════════════════════════════════════════════════════════════
 *  4. MODEL SIZE (per-subset S, M objectives)
 * ═══════════════════════════════════════════════════════════════════════════
 *
 *   Default nonlinear formulation (M >= 2), per subset:
 *     Auxiliary variables:     M (term bindings) + 1 (product result)
 *     Linear constraints:      M (term definitions)
 *     General constraints:     1 (nonlinear product)
 *     Quadratic constraints:   0
 *
 *   Retained balanced binary-tree formulation, per subset:
 *     Auxiliary variables:     M (term bindings) + M − 1 (product levels) = 2M − 1
 *     Linear constraints:      M (term definitions)
 *     Quadratic constraints:   M − 1 (pairwise products)
 *
 *   Balanced binary-tree total (K reference points, 2^K − 1 subsets):
 *     Variables:      (2^K − 1) · (2M − 1)
 *     Linear constrs: (2^K − 1) · M
 *     Quad constrs:   (2^K − 1) · (M − 1)
 */

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

#include "gurobi_c++.h"

#include "hs_moco/core/strategy.hpp"
#include "hs_moco/core/types.hpp"
#include "hs_moco/strategies/norm.hpp"

namespace hs_moco {

/**
 * @class MIQPStrategy
 * @brief Models the HS product directly via nonlinear general constraints.
 */
class MIQPStrategy final : public HSStrategy {
 public:
  MIQPStrategy() = default;
  explicit MIQPStrategy(NormalizationMethod method) : normalization_method_(method) {}

  [[nodiscard]] std::string name() const override { return "MIQP"; }

  void print_statistics(const GRBModel& model,
                        std::ostream& os = std::cout) const override {
    os << "[MIQP] Strategy statistics\n";
    print_model_statistics(model, os);
  }

  /**
   * @brief Requires NonConvex=2 for bilinear equality constraints.
   */
  void configure_options(SolverOptions& opts) const override {
    opts.allow_nonconvex = true;
    opts.feasibility_tol = 1e-9;
  }

  [[nodiscard]] int64_t decode_objective(double raw_value) const noexcept override {
    return static_cast<int64_t>(std::llround(raw_value * objective_multiplier_));
  }

  /**
   * @brief Build the inclusion–exclusion polynomial product objective.
   *
   * For each non-empty subset S ⊆ {0..K-1}:
   *  1. Compute the component-wise max reference: r_S^d = max_{j∈S} refs[j][d].
   *  2. Build M normalized linear terms: t_d = (f_d(x) − r_S^d) / norm_d.
   *  3. Define P_S = Π_d t_d with one nonlinear product constraint.
   *  4. Weight by inclusion–exclusion sign: (−1)^{|S|+1}.
   *
   * The final objective is the linear combination: obj = Σ_S sign_S · P_S.
   */
  void build_objective(GRBModel& model,
                       const StrategyContext& ctx,
                       const std::vector<std::vector<int64_t>>& refs,
                       AuxiliaryState& aux) override {
    const int M = ctx.M;
    const int k = static_cast<int>(refs.size());
    std::vector<double> norm(M, 1.0);
    if (normalization_method_ != NormalizationMethod::None) {
      const double abs_tol = std::max(1e-8, 10.0 * model.get(GRB_DoubleParam_FeasibilityTol));
      norm = normalization_factors(normalization_method_, ctx.objective_upper_bounds, refs, abs_tol);
    }
    double log_multiplier = 0;
    for (double value : norm) log_multiplier += std::log(value);
    objective_multiplier_ = std::exp(log_multiplier);

    // Inclusion–exclusion: Σ over non-empty subsets S ⊆ {0..k-1}
    //   sign = (−1)^{|S|+1}, product = Π_i (f_i − max_S r_i)
    GRBLinExpr obj_expr = 0;
    const int total_masks = (1 << k);

    for (int mask = 1; mask < total_masks; ++mask) {
      // ── Compute component-wise max reference for this subset ────────
      std::vector<int64_t> max_ref(M, std::numeric_limits<int64_t>::min());
      int bits = 0;
      for (int j = 0; j < k; ++j) {
        if (mask & (1 << j)) {
          ++bits;
          for (int i = 0; i < M; ++i)
            max_ref[i] = std::max(max_ref[i], refs[j][i]);
        }
      }

      // ── Build product terms: (f_i − max_ref_i) for each objective ───
      std::vector<GRBLinExpr> terms;
      terms.reserve(M);
      for (int i = 0; i < M; ++i) {
        terms.emplace_back((ctx.obj_vars[i] - max_ref[i]) / norm[i]);
      }

      // ── Reduce to a single product variable ────────────────────────
      GRBVar P = product_of_terms_nl(model, terms,
                                  "PS_" + std::to_string(mask), aux);

      // ── Inclusion–exclusion sign ───────────────────────────────────
      int sign = (bits % 2 == 1) ? +1 : -1;
      obj_expr += sign * GRBLinExpr(P);
    }

    model.setObjective(obj_expr, GRB_MAXIMIZE);
  }

 private:
  NormalizationMethod normalization_method_ = NormalizationMethod::None;
  double objective_multiplier_ = 1.0;

  /**
   * @brief Tie a linear expression to a new continuous variable.
   */
  static GRBVar var_from_expr(GRBModel& model, const GRBLinExpr& e,
                              const std::string& name, AuxiliaryState& aux) {
    GRBVar v = aux.addVar(model, -GRB_INFINITY, GRB_INFINITY, 0.0,
                          GRB_CONTINUOUS, name);
    aux.addConstr(model, v == e, name + "_def");
    return v;
  }

  /**
   * @brief Define a product with one n-ary nonlinear multiplication constraint.
   */
  static GRBVar product_of_terms_nl(GRBModel& model,
                                    const std::vector<GRBLinExpr>& terms,
                                    const std::string& base_name,
                                    AuxiliaryState& aux) {
    if (terms.size() < 2) return product_of_terms(model, terms, base_name, aux);

    std::vector<GRBVar> factors;
    factors.reserve(terms.size());
    for (size_t i = 0; i < terms.size(); ++i) {
      factors.push_back(
          var_from_expr(model, terms[i], base_name + "_t" + std::to_string(i), aux));
    }
    GRBVar P = aux.addVar(model, -GRB_INFINITY, GRB_INFINITY, 0.0,
                          GRB_CONTINUOUS, base_name);
    // Nonlinear variable nodes require indices assigned by a model update.
    model.update();

    const int nnodes = static_cast<int>(terms.size()) + 1;
    std::vector<int> opcode(nnodes, GRB_OPCODE_VARIABLE);
    std::vector<double> data(nnodes, -1.0);
    std::vector<int> parent(nnodes, 0);
    opcode[0] = GRB_OPCODE_MULTIPLY;
    parent[0] = -1;
    for (size_t i = 0; i < factors.size(); ++i) data[i + 1] = factors[i].index();
    aux.addGenConstrNL(model, P, nnodes, opcode.data(), data.data(), parent.data(),
                       base_name + "_nl");
    return P;
  }

  /**
   * @brief Balanced binary-tree reduction of M linear terms into one product variable.
   *
   * For 2 terms (a,b): P = a*b  (1 aux var, 1 QConstr).
   * For 3 terms (a,b,c): ab = a*b, P = ab*c.
   * For 4 terms: ab = a*b, cd = c*d, P = ab*cd.
   * General: O(M) aux vars and QConstrs.
   */
  static GRBVar product_of_terms(GRBModel& model,
                                 const std::vector<GRBLinExpr>& terms,
                                 const std::string& base_name,
                                 AuxiliaryState& aux) {
    if (terms.empty()) {
      // Degenerate: constant 1
      return aux.addVar(model, 1.0, 1.0, 0.0, GRB_CONTINUOUS, base_name + "_one");
    }

    // Wrap each term in a variable
    std::vector<GRBVar> partials;
    partials.reserve(terms.size());
    for (size_t i = 0; i < terms.size(); ++i) {
      partials.push_back(
          var_from_expr(model, terms[i], base_name + "_t" + std::to_string(i), aux));
    }

    // Pairwise reduction
    int level = 0;
    while (partials.size() > 1) {
      std::vector<GRBVar> next;
      next.reserve((partials.size() + 1) / 2);
      for (size_t i = 0; i < partials.size(); i += 2) {
        if (i + 1 < partials.size()) {
          std::string vn = base_name + "_L" + std::to_string(level) +
                           "_p" + std::to_string(i / 2);
          std::string qn = base_name + "_L" + std::to_string(level) +
                           "_q" + std::to_string(i / 2);
          GRBVar pij = aux.addVar(model, -GRB_INFINITY, GRB_INFINITY, 0.0,
                                  GRB_CONTINUOUS, vn);
          aux.addQConstr(model, pij == partials[i] * partials[i + 1], qn);
          next.push_back(pij);
        } else {
          next.push_back(partials[i]);  // carry odd element
        }
      }
      partials.swap(next);
      ++level;
    }

    return partials.front();
  }
};

}  // namespace hs_moco

#endif  // HS_MOCO_MIQP_STRATEGY_HPP
