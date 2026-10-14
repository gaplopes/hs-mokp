#ifndef HS_MOCO_MOKP_PROBLEM_HPP
#define HS_MOCO_MOKP_PROBLEM_HPP

/**
 * @file mokp.hpp
 * @brief Multiobjective Knapsack Problem (MOKP) definition and setup.
 *
 * The MOKP is defined by:
 *   - N items, each with a weight w_j and M objective values v_{ij}.
 *   - A knapsack capacity W.
 *   - Decision variables x_j ∈ {0,1} (select / skip item j).
 *   - Objectives f_i(x) = Σ_j v_{ij} · x_j  (maximize all i simultaneously).
 *   - Constraint: Σ_j w_j · x_j ≤ W.
 */

#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "gurobi_c++.h"

namespace hs_moco {

/**
 * @brief Item representation for the MOKP.
 *
 * Lightweight value type — no virtual dispatch, cheap to copy.
 */
struct MOKPItem {
  int32_t              idx;     ///< Item index.
  int64_t              weight;  ///< Item weight.
  std::vector<int64_t> values;  ///< Per-objective values [v_0, …, v_{M-1}].

  MOKPItem(int32_t idx, int64_t weight, std::vector<int64_t> values)
      : idx(idx), weight(weight), values(std::move(values)) {}
};

/**
 * @brief Multiobjective Knapsack Problem.
 *
 * Implements the problem interface for the standard 0-1 MOKP.
 *
 * @code
 *   std::vector<hs_moco::MOKPItem> items = { ... };
 *   auto problem = std::make_unique<hs_moco::MOKProblem>(items, W);
 * @endcode
 */
class MOKProblem final {
 public:
  /**
   * @brief Construct from items and capacity.
   *
   * @param items  Vector of MOKPItem (moved into the problem).
   * @param W      Knapsack capacity.
   * @throws std::invalid_argument if items is empty or objectives are inconsistent.
   */
  MOKProblem(std::vector<MOKPItem> items, int64_t W)
      : items_(std::move(items)),
        W_(W),
        N_(static_cast<int>(items_.size())),
        M_(items_.empty() ? 0 : static_cast<int>(items_[0].values.size())) {
    if (items_.empty())
      throw std::invalid_argument("MOKProblem: items must not be empty");
    // Validate consistent objective dimensionality
    for (const auto& it : items_) {
      if (static_cast<int>(it.values.size()) != M_)
        throw std::invalid_argument(
            "MOKProblem: all items must have the same number of objective values");
    }
  }

  // ── Problem interface ──────────────────────────────────────────────

  [[nodiscard]] int num_objectives() const noexcept { return M_; }
  [[nodiscard]] int num_variables()  const noexcept { return N_; }

  /**
   * @brief Add N binary decision variables x_0 … x_{N-1}.
   */
  std::vector<GRBVar> setup_variables(GRBModel& model,
                                      char var_type) const {
    std::vector<GRBVar> x(N_);
    for (int j = 0; j < N_; ++j) {
      x[j] = model.addVar(0.0, 1.0, 0.0, var_type,
                           "x_" + std::to_string(j));
    }
    return x;
  }

  /**
   * @brief Add the knapsack capacity constraint: Σ w_j·x_j ≤ W.
   */
  void setup_constraints(GRBModel& model,
                         const std::vector<GRBVar>& x) const {
    GRBLinExpr cap = 0;
    for (int j = 0; j < N_; ++j)
      cap += items_[j].weight * x[j];
    model.addConstr(cap <= W_, "knapsack_capacity");
  }

  /**
   * @brief Build M linear objective expressions f_i(x) = Σ_j v_{ij}·x_j.
   */
  std::vector<GRBLinExpr> build_objective_exprs(
      GRBModel& /*model*/,
      const std::vector<GRBVar>& x) const {
    std::vector<GRBLinExpr> f(M_);
    for (int i = 0; i < M_; ++i) {
      f[i] = 0;
      for (int j = 0; j < N_; ++j)
        f[i] += items_[j].values[i] * x[j];
    }
    return f;
  }

  /**
   * @brief Return the per-objective, per-variable coefficient matrix.
   *
   * C[i][j] = value of item j in objective i.
   * Required by MILPStrategy.
   */
  [[nodiscard]] std::vector<std::vector<int64_t>> coefficient_matrix() const {
    std::vector<std::vector<int64_t>> C(M_, std::vector<int64_t>(N_));
    for (int i = 0; i < M_; ++i)
      for (int j = 0; j < N_; ++j)
        C[i][j] = items_[j].values[i];
    return C;
  }

  /** @brief Dantzig bounds for nonnegative weights/capacity; otherwise positive-profit sums. */
  [[nodiscard]] std::vector<double> objective_upper_bounds() const {
    std::vector<double> bounds(M_, 0.0);
    const bool nonstandard_weights = W_ < 0 || std::any_of(
        items_.begin(), items_.end(), [](const MOKPItem& item) { return item.weight < 0; });
    for (int d = 0; d < M_; ++d) {
      long double bound = 0;
      std::vector<int> order;
      for (int j = 0; j < N_; ++j) {
        if (items_[j].values[d] <= 0) continue;
        // The fractional greedy bound assumes nonnegative weights and capacity.
        if (nonstandard_weights || items_[j].weight == 0)
          bound += items_[j].values[d];
        else
          order.push_back(j);
      }
      std::sort(order.begin(), order.end(), [&](int a, int b) {
        return static_cast<long double>(items_[a].values[d]) / items_[a].weight >
               static_cast<long double>(items_[b].values[d]) / items_[b].weight;
      });
      long double remaining = W_;
      for (int j : order) {
        if (remaining <= 0) break;
        const long double fraction = std::min(1.0L, remaining / items_[j].weight);
        bound += fraction * items_[j].values[d];
        remaining -= fraction * items_[j].weight;
        if (fraction < 1) break;
      }
      bounds[d] = static_cast<double>(bound);
    }
    return bounds;
  }

  // ── MOKP-specific accessors ────────────────────────────────────────────

  [[nodiscard]] int64_t capacity() const noexcept { return W_; }
  [[nodiscard]] const std::vector<MOKPItem>& items() const noexcept { return items_; }

 private:
  const std::vector<MOKPItem> items_;
  const int64_t W_;
  const int N_;
  const int M_;
};

}  // namespace hs_moco

#endif  // HS_MOCO_MOKP_PROBLEM_HPP
