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
 * sequence of bilinear (quadratic) equalities. Because Gurobi's default
 * QP solver only handles convex objectives, the NonConvex=2 flag is
 * required to enable the spatial branch-and-bound for non-convex MIQCP.
 *
 * ═══════════════════════════════════════════════════════════════════════════
 *  2. BALANCED BINARY-TREE PRODUCT REDUCTION
 * ═══════════════════════════════════════════════════════════════════════════
 *
 * Given M linear terms t_0, t_1, …, t_{M-1}, their product is computed by
 * a balanced binary-tree reduction that minimizes the longest chain of
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
 * Additionally, each linear term t_d = f_d(x) − r_d is bound to an
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
 * resides in the quadratic constraints that define the products.
 *
 * ═══════════════════════════════════════════════════════════════════════════
 *  4. MODEL SIZE (per-subset S, M objectives)
 * ═══════════════════════════════════════════════════════════════════════════
 *
 *   Per subset:
 *     Auxiliary variables:     M (term bindings) + M − 1 (product levels) = 2M − 1
 *     Linear constraints:      M (term definitions)
 *     Quadratic constraints:   M − 1 (pairwise products)
 *
 *   Total (K reference points, 2^K − 1 subsets):
 *     Variables:      (2^K − 1) · (2M − 1)
 *     Linear constrs: (2^K − 1) · M
 *     Quad constrs:   (2^K − 1) · (M − 1)
 */

#include <algorithm>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

#include "gurobi_c++.h"

#include "hs_moco/core/strategy.hpp"
#include "hs_moco/core/types.hpp"

namespace hs_moco {

/**
 * @class MIQPStrategy
 * @brief Models the HS product directly via bilinear quadratic constraints.
 */
class MIQPStrategy final : public HSStrategy {
 public:
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

  /**
   * @brief Build the inclusion–exclusion polynomial product objective.
   *
   * For each non-empty subset S ⊆ {0..K-1}:
   *  1. Compute the component-wise max reference: r_S^d = max_{j∈S} refs[j][d].
   *  2. Build M linear terms: t_d = f_d(x) − r_S^d.
   *  3. Reduce to a single product P_S = Π_d t_d via binary-tree reduction.
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
        terms.emplace_back(ctx.obj_vars[i] - max_ref[i]);
      }

      // ── Reduce to a single product variable ────────────────────────
      GRBVar P = product_of_terms(model, terms,
                                  "PS_" + std::to_string(mask), aux);

      // ── Inclusion–exclusion sign ───────────────────────────────────
      int sign = (bits % 2 == 1) ? +1 : -1;
      obj_expr += sign * GRBLinExpr(P);
    }

    model.setObjective(obj_expr, GRB_MAXIMIZE);
  }

 private:
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
