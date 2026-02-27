#ifndef HS_MOCO_MILP_STRATEGY_HPP
#define HS_MOCO_MILP_STRATEGY_HPP

/**
 * @file strategies/milp.hpp
 * @brief MILPStrategy — McCormick-envelope linearization of the HS product.
 *
 * ═══════════════════════════════════════════════════════════════════════════
 *  1. MATHEMATICAL BACKGROUND
 * ═══════════════════════════════════════════════════════════════════════════
 *
 * The Hypervolume Scalarization (HS) for a single reference point r is:
 *
 *   HS(x, r) = Π_{d=0}^{M-1} (f_d(x) − r_d)
 *
 * where f_d(x) = Σ_j C[d][j] · x_j  is the d-th objective and C[d][j] is
 * the coefficient of decision variable x_j in objective d.
 *
 * For multiple reference points (inclusion–exclusion over subsets S ⊆ R):
 *
 *   HS(x, R) = Σ_{S⊆R, S≠∅} (−1)^{|S|+1} Π_d (f_d(x) − max_{r∈S} r_d)
 *
 * ═══════════════════════════════════════════════════════════════════════════
 *  2. POLYNOMIAL EXPANSION
 * ═══════════════════════════════════════════════════════════════════════════
 *
 * Each product Π_d (f_d(x) − r_d) can be expanded into a multilinear
 * polynomial in the individual decision variables x_j. For each subset, we
 * substitute f_d(x) = Σ_j C[d][j]·x_j and expand the product of M linear
 * expressions. This yields terms of degree 0 through M:
 *
 *   Π_d (f_d(x) − r_d) = Σ_{m=0}^{M}  (−1)^{M−m}  Σ_{t ∈ Tuples(N,m)}
 *                           coeff(t, r) · x_{t_0} · x_{t_1} · … · x_{t_{m-1}}
 *
 * The coefficient coeff(t, r) for a tuple t = (j_0, …, j_{m-1}) and
 * reference r is computed by enumerating which m of the M dimensions take a
 * variable factor vs. a reference-point factor:
 *
 *   coeff(t, r) = Σ_{D ⊂ {0..M-1}, |D|=m}
 *                   ( Π_{d ∉ D} r_d )                     ← reference part
 *                 · ( Σ_{π ∈ Perm(t)} Π_{p} C[D_p][π_p] ) ← variable part
 *
 * The outer sum iterates over all (M choose m) ways to select which
 * dimensions are "filled" by variables.  The inner sum iterates over all
 * distinct permutations of the tuple, assigning each variable to one of
 * the selected dimensions.
 *
 * ═══════════════════════════════════════════════════════════════════════════
 *  3. BINARY-VARIABLE SIMPLIFICATION
 * ═══════════════════════════════════════════════════════════════════════════
 *
 * Since x_j ∈ {0,1}, we have x_j² = x_j.  Therefore any multilinear
 * monomial with repeated indices collapses:
 *
 *   x_0 · x_0 · x_1  =  x_0 · x_1
 *
 * This means we only need linearization variables for *strictly increasing*
 * index tuples {i < j < k < …}, not combinations with replacement.
 *
 * ═══════════════════════════════════════════════════════════════════════════
 *  4. McCORMICK ENVELOPE LINEARIZATION
 * ═══════════════════════════════════════════════════════════════════════════
 *
 * Every product of binary variables is replaced by a continuous auxiliary
 * variable and McCormick envelope constraints. These are exact (tight) for
 * binary variables, turning a multilinear program into a pure MILP.
 *
 *   Degree 2:  y_{ij} ≈ x_i · x_j
 *     y_{ij} ≤ x_i
 *     y_{ij} ≤ x_j
 *     y_{ij} ≥ x_i + x_j − 1
 *     (y_{ij} ≥ 0 is implicit in variable bounds)
 *
 *   Degree 3:  w_{ijk} ≈ y_{ij} · x_k  (and symmetrically via y_{ik}, y_{jk})
 *     For each decomposition (t \ {p}, {p}):
 *       w_{ijk} ≤ y_{..}    (the sub-product variable)
 *       w_{ijk} ≤ x_p
 *       w_{ijk} ≥ y_{..} + x_p − 1
 *
 *     Three decompositions per triple (p = i, j, k) give 9 constraints,
 *     tightening the LP relaxation from all directions.
 *
 *   Higher degrees generalize analogously, decomposing a degree-m variable
 *   into a degree-(m−1) sub-product times a single x, for each of the m
 *   possible single variables to split off.
 *
 * ═══════════════════════════════════════════════════════════════════════════
 *  5. IMPLEMENTATION OUTLINE
 * ═══════════════════════════════════════════════════════════════════════════
 *
 *   Step 1 — Precompute inclusion–exclusion subsets: for each non-empty
 *            S ⊆ {0..K-1}, store (max_S r_d, |S|).
 *
 *   Step 2 — Generate all combinations-with-replacement of N items taken
 *            m at a time (for objective terms), and all strict combinations
 *            (for linearization variables).
 *
 *   Step 3 — calc_coeff(t, r): the universal coefficient generator described
 *            above, using bitmask enumeration of dimension subsets and
 *            permutation iteration.
 *
 *   Step 4 — Create linearization variables: y_{ij} for degree 2,
 *            w_{ijk} for degree 3, etc. — only strict tuples.
 *
 *   Step 5 — Build objective: for each degree m = 1..M, iterate over all
 *            tuples (with replacement), compute the combined inclusion–
 *            exclusion coefficient, collapse repeated indices (x²=x), and
 *            add to the GRBLinExpr. Also add the constant term (m = 0).
 *
 *   Step 6 — Add McCormick envelope constraints for all linearization vars.
 *
 * ═══════════════════════════════════════════════════════════════════════════
 *  6. MODEL SIZE (for M objectives, N variables, K reference points)
 * ═══════════════════════════════════════════════════════════════════════════
 *
 *   Linearization variables:  Σ_{m=2}^{M} C(N, m)
 *     M=2:  C(N,2) = N(N−1)/2
 *     M=3:  C(N,2) + C(N,3) = N(N−1)/2 + N(N−1)(N−2)/6
 *
 *   McCormick constraints:
 *     Degree 2:  3 · C(N,2)
 *     Degree 3:  9 · C(N,3)    (3 decompositions × 3 constraints each)
 *     Degree m:  3m · C(N,m)
 *
 *   Pure MILP — no quadratic or general constraints.
 *
 *   Key requirement: Needs per-variable, per-objective coefficients from
 *                    the problem (via coefficient_matrix()).
 */

#include <algorithm>
#include <cstdint>
#include <functional>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <string>
#include <vector>

#include "gurobi_c++.h"

#include "hs_moco/core/strategy.hpp"
#include "hs_moco/core/types.hpp"

namespace hs_moco {

/**
 * @class MILPStrategy
 * @brief Linearizes the multilinear HS polynomial via McCormick envelopes.
 */
class MILPStrategy final : public HSStrategy {
 public:
  [[nodiscard]] std::string name() const override { return "MILP"; }

  void print_statistics(const GRBModel& model,
                        std::ostream& os = std::cout) const override {
    os << "[MILP] Strategy statistics\n";
    print_model_statistics(model, os);
  }

  /**
   * @brief No special model-level flags needed (pure MILP).
   */
  void configure_options(SolverOptions& /*opts*/) const override {
    // Nothing to override — standard MILP settings suffice.
  }

  /**
   * @brief Build the McCormick-linearized polynomial HS objective.
   *
   * Performs six steps (see file-level documentation for the math):
   *  1. Precompute inclusion–exclusion subsets of reference points.
   *  2. Generate index tuples (with-replacement for coefficients,
   *     strict for linearization variables).
   *  3. Compute polynomial coefficients via calc_coeff().
   *  4. Create linearization variables for each strict tuple.
   *  5. Assemble the GRBLinExpr objective (collapsing x²=x).
   *  6. Add McCormick envelope constraints.
   *
   * @pre ctx.coeff_matrix must be non-empty (M × N matrix).
   * @throws std::runtime_error if coefficient_matrix is missing.
   */
  void build_objective(GRBModel& model,
                       const StrategyContext& ctx,
                       const std::vector<std::vector<int64_t>>& refs,
                       AuxiliaryState& aux) override {
    const int M = ctx.M;
    const int N = ctx.N;

    if (ctx.coeff_matrix.empty()) {
      throw std::runtime_error(
          "MILPStrategy requires coefficient_matrix() from the problem.");
    }

    // ── 1. Precompute inclusion-exclusion subsets ─────────────────────────
    const int K = static_cast<int>(refs.size());
    struct SubsetMax {
      std::vector<int64_t> r;
      int subsetSize;
    };
    std::vector<SubsetMax> subsets;
    for (int mask = 1; mask < (1 << K); ++mask) {
      std::vector<int64_t> r(M, std::numeric_limits<int64_t>::min());
      int bits = 0;
      for (int m = 0; m < K; ++m) {
        if (mask & (1 << m)) {
          ++bits;
          for (int d = 0; d < M; ++d) r[d] = std::max(r[d], refs[m][d]);
        }
      }
      subsets.push_back({std::move(r), bits});
    }

    // ── 2. Combination generators ────────────────────────────────────────
    auto get_combinations_with_replacement = [](int n, int m) {
      std::vector<std::vector<int>> res;
      std::vector<int> current;
      std::function<void(int, int)> bt = [&](int start, int depth) {
        if (depth == m) { res.push_back(current); return; }
        for (int i = start; i < n; ++i) {
          current.push_back(i);
          bt(i, depth + 1);
          current.pop_back();
        }
      };
      bt(0, 0);
      return res;
    };

    auto get_strict_combinations = [](int n, int m) {
      std::vector<std::vector<int>> res;
      std::vector<int> current;
      std::function<void(int, int)> bt = [&](int start, int depth) {
        if (depth == m) { res.push_back(current); return; }
        for (int i = start; i < n; ++i) {
          current.push_back(i);
          bt(i + 1, depth + 1);
          current.pop_back();
        }
      };
      bt(0, 0);
      return res;
    };

    // ── 3. Universal polynomial coefficient generator ────────────────────
    //
    // calc_coeff(t, r) computes the coefficient of the monomial x_{t_0}·…·x_{t_{m-1}}
    // arising from expanding Π_d (f_d(x) − r_d).
    //
    // Algorithm:
    //   For each way to choose m dimensions out of M (bitmask over dims):
    //     - The chosen dimensions take variable factors from the tuple.
    //     - The remaining (M − m) dimensions contribute r_d factors.
    //     - Sum over all permutations of the tuple across chosen dimensions.
    //
    // Example (M=3, tuple = {i,j}, r = {r0,r1,r2}):
    //   Bitmask 011 → dims {0,1} for vars, dim {2} for r:
    //     r2 · Σ_{perm(i,j)} C[0][π0]·C[1][π1]
    //     = r2 · (C[0][i]·C[1][j] + C[0][j]·C[1][i])
    //   Bitmask 101 → dims {0,2}, dim {1}: r1 · (...)
    //   Bitmask 110 → dims {1,2}, dim {0}: r0 · (...)
    //
    const auto& C = ctx.coeff_matrix;  // C[d][j] = coefficient of x_j in objective d
    auto calc_coeff = [&](const std::vector<int>& var_idx,
                          const std::vector<int64_t>& r) -> long long {
      long long total = 0;
      int num_vars = static_cast<int>(var_idx.size());

      for (int mask = 0; mask < (1 << M); ++mask) {
        int bits = 0;
        for (int d = 0; d < M; ++d)
          if (mask & (1 << d)) bits++;
        if (bits != num_vars) continue;

        std::vector<int> var_dims, r_dims;
        for (int d = 0; d < M; ++d) {
          if (mask & (1 << d)) var_dims.push_back(d);
          else                 r_dims.push_back(d);
        }

        long long r_prod = 1;
        for (int d : r_dims) r_prod *= r[d];

        std::vector<int> perm = var_idx;
        std::sort(perm.begin(), perm.end());
        do {
          long long ip = 1;
          for (size_t p = 0; p < perm.size(); ++p)
            ip *= C[var_dims[p]][perm[p]];
          total += ip * r_prod;
        } while (std::next_permutation(perm.begin(), perm.end()));
      }
      return total;
    };

    // ── 4. Variable generation (strict index tuples only) ────────────────
    //
    // Since x_j ∈ {0,1} ⇒ x_j² = x_j, any monomial with repeated indices
    // collapses to a lower-degree monomial over the *unique* indices.
    // Therefore we only create linearization variables for strictly
    // increasing tuples {i < j}, {i < j < k}, etc.
    //
    // The map key is the sorted tuple; degree-1 maps to the original x_j.
    //
    std::map<std::vector<int>, GRBVar> vars;

    // Degree-1 variables are the decision variables themselves
    for (int i = 0; i < N; ++i) {
      vars[{i}] = ctx.x_vars[i];
    }

    // Higher-degree linearization variables: y_{ij}, w_{ijk}, ...
    for (int m = 2; m <= M; ++m) {
      auto strict_tuples = get_strict_combinations(N, m);
      for (const auto& t : strict_tuples) {
        vars[t] = aux.addVar(model, 0.0, 1.0, 0.0, GRB_CONTINUOUS,
                             "lin_" + tuple_name(t));
      }
    }

    // ── 5. Objective construction ────────────────────────────────────────
    //
    // Assemble:  obj = Σ_{m=0}^{M} (−1)^{M−m} · Σ_t coeff_sum(t) · vars[collapse(t)]
    //
    // For the highest-degree terms (m = M), the reference part is empty,
    // so coeff is independent of the reference points. For lower degrees
    // (m < M), inclusion–exclusion over subsets modulates the coefficients.
    //
    // The "collapse" step maps tuples with repeats to their unique-index
    // equivalent (e.g., {0,0,1} → {0,1}), exploiting x² = x for binaries.
    //
    GRBLinExpr obj = 0;
    std::vector<int64_t> zero_r(M, 0);

    for (int m = 1; m <= M; ++m) {
      int sign_M = ((M - m) % 2 == 0) ? +1 : -1;  // (−1)^{M−m}
      auto all_tuples = get_combinations_with_replacement(N, m);

      for (const auto& t : all_tuples) {
        long long coeff_sum = 0;
        if (m == M) {
          coeff_sum = calc_coeff(t, zero_r);
        } else {
          for (const auto& sub : subsets) {
            int sign_sub = (sub.subsetSize % 2 == 1) ? +1 : -1;
            coeff_sum += sign_sub * calc_coeff(t, sub.r);
          }
        }
        coeff_sum *= sign_M;

        if (coeff_sum != 0) {
          // Collapse repeated vars: x_0·x_0·x_1 → y_{0,1} (since x_i²=x_i for binary)
          std::set<int> unique(t.begin(), t.end());
          std::vector<int> collapsed(unique.begin(), unique.end());
          obj += coeff_sum * vars[collapsed];
        }
      }
    }

    // Constant term (m = 0)
    long long const_sum = 0;
    int sign_M_0 = (M % 2 == 0) ? +1 : -1;
    for (const auto& sub : subsets) {
      int sign_sub = (sub.subsetSize % 2 == 1) ? +1 : -1;
      const_sum += sign_sub * calc_coeff({}, sub.r);
    }
    obj += sign_M_0 * const_sum;

    model.setObjective(obj, GRB_MAXIMIZE);

    // ── 6. McCormick envelopes ───────────────────────────────────────────
    //
    // For each linearization variable (representing a product of binaries),
    // add McCormick envelope constraints that make the LP relaxation tight.
    //
    // General rule for z ≈ a · b where a, b ∈ [0,1]:
    //   z ≤ a,  z ≤ b,  z ≥ a + b − 1,  z ≥ 0
    //
    // These are exact (not a relaxation) when a, b are binary.
    //
    // For degree-m variables, we decompose into (degree-(m-1)) × (degree-1)
    // and apply McCormick to each decomposition, giving 3·m constraints
    // per variable. Multiple decompositions tighten the LP relaxation.
    //
    for (int m = 2; m <= M; ++m) {
      auto strict_tuples = get_strict_combinations(N, m);
      for (const auto& t : strict_tuples) {
        if (m == 2) {
          // For degree-2, the two possible splits {i}×{j} and {j}×{i}
          // produce identical McCormick constraints (product is symmetric).
          // Emit each constraint exactly once:
          //   lin_{ij} <= x_i,  lin_{ij} <= x_j,
          //   lin_{ij} >= x_i + x_j - 1
          std::string tag = "mc_" + tuple_name(t);
          aux.addConstr(model, vars[t] <= vars[{t[0]}], tag + "_ub1");
          aux.addConstr(model, vars[t] <= vars[{t[1]}], tag + "_ub2");
          aux.addConstr(model, vars[t] >= vars[{t[0]}] + vars[{t[1]}] - 1,
                        tag + "_lb");
        } else {
          // For degree >= 3, each split through a different sub-tuple
          // produces distinct constraints, so emit all of them.
          for (size_t p = 0; p < t.size(); ++p) {
            std::vector<int> t_minus_p = t;
            t_minus_p.erase(t_minus_p.begin() + static_cast<long>(p));
            int x_idx = t[p];

            std::string tag = "mc_" + tuple_name(t) + "_" + std::to_string(p);
            aux.addConstr(model, vars[t] <= vars[t_minus_p], tag + "_ub1");
            aux.addConstr(model, vars[t] <= vars[{x_idx}],   tag + "_ub2");
            aux.addConstr(model, vars[t] >= vars[t_minus_p] + vars[{x_idx}] - 1,
                          tag + "_lb");
          }
        }
      }
    }
  }

 private:
  static std::string tuple_name(const std::vector<int>& t) {
    std::string s;
    for (size_t i = 0; i < t.size(); ++i) {
      if (i > 0) s += "_";
      s += std::to_string(t[i]);
    }
    return s;
  }
};

}  // namespace hs_moco

#endif  // HS_MOCO_MILP_STRATEGY_HPP
