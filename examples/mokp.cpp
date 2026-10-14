/**
 * @file mokp.cpp
 * @brief Compares every HS-MOCO strategy on an MOKP instance.
 *
 * For reference-point set sizes k = 1 … K.
 * For each k, it generates 30 random nondominated reference-point sets (each of size k):
 *   1. Computes the expected solution via exhaustive enumeration over the
 *      known Pareto front (simulate_solution — the validation oracle).
 *   2. Runs each strategy (MIQP, Log, MILP) on the same reference points.
 *   3. Measures wall-clock time and verifies the solver output matches.
 *
 * Usage:
 *   ./example_mokp                              # default instance, K=5
 *   ./example_mokp /path/to/instance.in         # custom instance
 *   ./example_mokp /path/to/instance.in 3       # custom instance, K=3
 */

#include <algorithm>
#include <charconv>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <functional>
#include <hs_moco/hs_moco.hpp>
#include <iomanip>
#include <iostream>
#include <memory>
#include <mooutils/indicators.hpp>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

// ── Instance data loaded from an .in file ───────────────────────────────────

struct MOKPInstance {
  std::vector<hs_moco::MOKPItem> items;
  int64_t capacity;
  int num_objs;
  std::vector<std::vector<int64_t>> pareto_front;

  [[nodiscard]] std::unique_ptr<hs_moco::MOKProblem> to_problem() const {
    auto items_copy = items;
    return std::make_unique<hs_moco::MOKProblem>(std::move(items_copy), capacity);
  }
};

/// Load a MOKP instance from an .in file (mobkp-instances format).
///
/// Format:
///   n m
///   W
///   w_1 v^1_1 … v^m_1
///   …
///   w_n v^1_n … v^m_n
///   nd
///   p^1_1 … p^m_1
///   …
///   p^1_nd … p^m_nd
[[nodiscard]] MOKPInstance load_mokp_instance(const std::string& filepath) {
  std::ifstream fin(filepath);
  if (!fin.is_open())
    throw std::runtime_error("load_mokp_instance: cannot open file: " + filepath);

  MOKPInstance inst;

  int n = 0;
  fin >> n >> inst.num_objs;
  if (!fin || n <= 0 || inst.num_objs <= 0)
    throw std::invalid_argument("load_mokp_instance: invalid header (n, m)");

  fin >> inst.capacity;
  if (!fin)
    throw std::invalid_argument("load_mokp_instance: invalid capacity W");

  inst.items.reserve(n);
  for (int j = 0; j < n; ++j) {
    int64_t weight = 0;
    fin >> weight;
    std::vector<int64_t> values(inst.num_objs);
    for (int i = 0; i < inst.num_objs; ++i)
      fin >> values[i];
    if (!fin)
      throw std::invalid_argument(
          "load_mokp_instance: failed to read item " + std::to_string(j));
    inst.items.emplace_back(static_cast<int32_t>(j), weight, std::move(values));
  }

  int nd = 0;
  if (fin >> nd && nd > 0) {
    inst.pareto_front.reserve(nd);
    for (int k = 0; k < nd; ++k) {
      std::vector<int64_t> point(inst.num_objs);
      for (int i = 0; i < inst.num_objs; ++i)
        fin >> point[i];
      if (!fin)
        throw std::invalid_argument(
            "load_mokp_instance: failed to read Pareto point " + std::to_string(k));
      inst.pareto_front.push_back(std::move(point));
    }
  }

  return inst;
}

// ── HypervolumeIndicator wrapper around mooutils ────────────────────────────

template <typename T, typename Solution>
class HypervolumeIndicator {
 public:
  HypervolumeIndicator(Solution ref_point, bool is_maximization)
      : ref_point_(std::move(ref_point)),
        is_maximization_(is_maximization),
        hv_space_(initialize_hv_space(ref_point_, is_maximization_)) {}

  T insert(const Solution& sol) {
    Solution aux = is_maximization_ ? sol : negative_vector(sol);
    return hv_space_.insert(aux);
  }

  T contribution(const Solution& sol) const {
    Solution aux = is_maximization_ ? sol : negative_vector(sol);
    return hv_space_.contribution(aux);
  }

  T value() const { return hv_space_.value(); }

  /// Insert every solution; return total HV if all contribute, else 0.
  T set_hypervolume(const std::vector<Solution>& solutions) {
    T hvc = 0;
    for (const auto& sol : solutions) {
      T c = insert(sol);
      if (c == 0) return 0;
      hvc += c;
    }
    return hvc;
  }

 private:
  static mooutils::incremental_hv<T, Solution> initialize_hv_space(
      const Solution& ref, bool is_max) {
    return is_max ? mooutils::incremental_hv<T, Solution>(ref)
                  : mooutils::incremental_hv<T, Solution>(negative_vector(ref));
  }

  static Solution negative_vector(const Solution& v) {
    Solution res(v.size());
    std::transform(v.begin(), v.end(), res.begin(), std::negate<T>());
    return res;
  }

  Solution ref_point_;
  bool is_maximization_;
  mooutils::incremental_hv<T, Solution> hv_space_;
};

// ── Expected-solution simulation ────────────────────────────────────────────

/// Given the known Pareto front and a set of reference points, compute the
/// expected best solution (and its HS value) by exhaustive enumeration.
std::pair<std::vector<int64_t>, int64_t> simulate_solution(
    const std::vector<std::vector<int64_t>>& pareto_front,
    int num_objs,
    const std::vector<std::vector<int64_t>>& ref_points) {
  auto strongly_dominates = [](const std::vector<int64_t>& v1,
                               const std::vector<int64_t>& v2) {
    for (size_t i = 0; i < v1.size(); ++i)
      if (v1[i] <= v2[i]) return false;
    return true;
  };

  // Union (component-wise max) of all reference points.
  std::vector<int64_t> union_point(static_cast<size_t>(num_objs), 0);
  for (const auto& rp : ref_points)
    for (size_t j = 0; j < union_point.size(); ++j)
      union_point[j] = std::max(union_point[j], rp[j]);

  std::vector<int64_t> origin(static_cast<size_t>(num_objs), 0);
  std::pair<std::vector<int64_t>, int64_t> best = {origin, 0};

  for (const auto& sol : pareto_front) {
    if (!strongly_dominates(sol, union_point)) continue;
    HypervolumeIndicator<int64_t, std::vector<int64_t>> hvi(sol, false);
    int64_t hv = hvi.set_hypervolume(ref_points);
    if (hv > best.second) {
      best.first = sol;
      best.second = hv;
    }
  }
  return best;
}

// ── Random nondominated reference-point generation ──────────────────────────

/// Compute the component-wise maximum across the Pareto front.
std::vector<int64_t> compute_max_coords(
    const std::vector<std::vector<int64_t>>& pareto_front, int num_objs) {
  std::vector<int64_t> max_coords(static_cast<size_t>(num_objs), 0);
  for (const auto& pt : pareto_front)
    for (size_t j = 0; j < max_coords.size(); ++j)
      max_coords[j] = std::max(max_coords[j], pt[j]);
  return max_coords;
}

/// Generate `k` random nondominated points in [0, max_coords] inclusive.
/// Points are filtered so that the returned set is mutually nondominated.
/// A fixed seed is used for reproducibility.
std::vector<std::vector<int64_t>> generate_nondominated_refs(
    const std::vector<int64_t>& max_coords, int k, unsigned seed = 42) {
  if (k <= 0 || max_coords.empty() ||
      std::any_of(max_coords.begin(), max_coords.end(),
                  [](int64_t coordinate) { return coordinate < 0; }))
    throw std::invalid_argument(
        "generate_nondominated_refs: k must be positive and coordinate bounds nonempty and nonnegative");

  const auto m = max_coords.size();
  std::mt19937_64 rng(seed);

  std::vector<std::uniform_int_distribution<int64_t>> dists;
  dists.reserve(m);
  for (size_t j = 0; j < m; ++j)
    dists.emplace_back(int64_t{0}, max_coords[j]);

  auto weakly_dominates = [](const std::vector<int64_t>& a,
                             const std::vector<int64_t>& b) {
    for (size_t i = 0; i < a.size(); ++i)
      if (a[i] < b[i]) return false;
    return true;
  };

  std::vector<std::vector<int64_t>> refs;
  // shortcut: sampling may reject feasible large sets; increase the budget if larger sets are needed.
  constexpr int max_attempts = 100000;
  for (int attempt = 0; attempt < max_attempts && static_cast<int>(refs.size()) < k;
       ++attempt) {
    std::vector<int64_t> pt(m);
    for (size_t j = 0; j < m; ++j)
      pt[j] = dists[j](rng);

    bool dominated = false;
    for (const auto& r : refs) {
      if (weakly_dominates(r, pt)) {
        dominated = true;
        break;
      }
    }
    if (dominated) continue;

    refs.erase(std::remove_if(refs.begin(), refs.end(),
                              [&](const std::vector<int64_t>& r) {
                                return weakly_dominates(pt, r);
                              }),
               refs.end());
    refs.push_back(std::move(pt));
  }
  if (static_cast<int>(refs.size()) < k)
    throw std::runtime_error(
        "generate_nondominated_refs: could not generate " + std::to_string(k) +
        " nondominated points after " + std::to_string(max_attempts) + " attempts");
  return refs;
}

// ── Helpers ─────────────────────────────────────────────────────────────────

void print_vec(const std::vector<int64_t>& v) {
  std::cout << "(";
  for (size_t i = 0; i < v.size(); ++i)
    std::cout << v[i] << (i + 1 < v.size() ? ", " : "");
  std::cout << ")";
}

/// Names of the three strategies we compare.
const std::vector<std::string> STRATEGY_NAMES = {"miqp", "log", "milp"};

/// Build a fresh strategy by name.
std::unique_ptr<hs_moco::HSStrategy> make_strategy(
    const std::string& name, hs_moco::NormalizationMethod normalization_method) {
  if (name == "log") return std::make_unique<hs_moco::LogStrategy>();
  if (name == "milp") return std::make_unique<hs_moco::MILPStrategy>();
  return std::make_unique<hs_moco::MIQPStrategy>(normalization_method);
}

/// Build a fresh solver for the given instance + strategy.
hs_moco::MOKPSolver make_solver(const MOKPInstance& inst,
                                const std::string& strategy_name,
                                const hs_moco::SolverOptions& opts,
                                hs_moco::NormalizationMethod normalization_method) {
  return hs_moco::MOKPSolver(
      inst.to_problem(), make_strategy(strategy_name, normalization_method), opts);
}

}  // namespace

int main(int argc, char* argv[]) {
  try {
    // ────────────────────────────────────────────────────────────────────────
    // 1. Parse arguments
    // ────────────────────────────────────────────────────────────────────────
    std::string instance_path;
    if (argc >= 2) {
      instance_path = argv[1];
    } else {
#ifdef MOBKP_INSTANCES_DIR
      instance_path = std::string(MOBKP_INSTANCES_DIR) + "/random/3D/50_1.in";
#else
      std::cerr << "Usage: " << argv[0] << " <instance.in> [max-k]\n";
      return 1;
#endif
    }

    if (!std::filesystem::exists(instance_path)) {
      std::cerr << "Error: instance file not found: " << instance_path << "\n";
      return 1;
    }

    int max_k = 5;
    if (argc >= 3) {
      const std::string argument = argv[2];
      const auto [end, error] = std::from_chars(
          argument.data(), argument.data() + argument.size(), max_k);
      if (error != std::errc{} || end != argument.data() + argument.size() || max_k <= 0)
        throw std::invalid_argument("max-k must be a positive integer");
    }

    // ────────────────────────────────────────────────────────────────────────
    // 2. Load instance
    // ────────────────────────────────────────────────────────────────────────
    std::cout << "Loading instance: " << instance_path << "\n";
    auto instance = load_mokp_instance(instance_path);

    std::cout << "  Items     : " << instance.items.size() << "\n";
    std::cout << "  Objectives: " << instance.num_objs << "\n";
    std::cout << "  Capacity  : " << instance.capacity << "\n";
    if (!instance.pareto_front.empty())
      std::cout << "  Nondominated set  : " << instance.pareto_front.size() << " points\n";
    std::cout << "\n";

    if (instance.pareto_front.empty()) {
      std::cerr << "Error: Pareto front required for validation.\n";
      return 1;
    }

    // ────────────────────────────────────────────────────────────────────────
    // 3. Solver options (shared across all strategies)
    // ────────────────────────────────────────────────────────────────────────
    hs_moco::SolverOptions opts;
    opts.time_limit = 600.0;
    opts.threads = 4;

    // Keep one method assignment active; None disables MIQP normalization.
    auto normalization_method = hs_moco::NormalizationMethod::None;

    auto max_coords = compute_max_coords(instance.pareto_front, instance.num_objs);
    std::cout << "Max. coords : ";
    print_vec(max_coords);
    std::cout << "\n";
    std::cout << "Strategies : ";
    for (const auto& s : STRATEGY_NAMES) std::cout << s << " ";
    std::cout << "\nMax. k      : " << max_k << "\n";
    std::cout << std::string(72, '=') << "\n\n";

    int total_tests = 0, total_pass = 0, total_fail = 0;

    // ────────────────────────────────────────────────────────────────────────
    // 4. For each k = 1 … max_k, generate refs and compare all strategies
    // ────────────────────────────────────────────────────────────────────────
    const int num_runs = 1;

    for (int k = 1; k <= max_k; ++k) {
      std::cout << "==== k=" << k << " (" << num_runs << " runs) ";
      std::cout << std::string(50, '=') << "\n\n";

      if (k > 1 && instance.num_objs == 2) {
        std::cerr << "Warning: skipping k=" << k
                  << " since instance has only " << instance.num_objs
                  << " objective(s).\n\n";
        continue;
      }

      for (int run = 0; run < num_runs; ++run) {
        // Generate k nondominated reference points (different seed per run).
        unsigned seed = static_cast<unsigned>(k * 10000 + run * 100 + 7);
        auto refs = generate_nondominated_refs(max_coords, k, seed);
        const int actual_k = static_cast<int>(refs.size());

        std::cout << "---- run=" << (run + 1) << "/" << num_runs << " ---- \n\n";
        std::cout << "  Refs: ";
        for (const auto& r : refs) {
          print_vec(r);
          std::cout << " ";
        }
        std::cout << "\n";

        // 4a. Compute expected solution.
        auto t0_sim = std::chrono::steady_clock::now();
        auto [exp_point, exp_hs] = simulate_solution(
            instance.pareto_front, instance.num_objs, refs);
        auto t1_sim = std::chrono::steady_clock::now();
        double sim_ms = std::chrono::duration<double, std::milli>(t1_sim - t0_sim).count();

        std::cout << "  Expected HS : " << exp_hs << "\n";
        std::cout << "  Expected pt : ";
        print_vec(exp_point);
        std::cout << "  (" << std::fixed << std::setprecision(2) << sim_ms << " ms)\n\n";

        // 4b. Run each strategy.
        for (const auto& sname : STRATEGY_NAMES) {
          auto solver = make_solver(instance, sname, opts, normalization_method);
          ++total_tests;

          auto t0 = std::chrono::steady_clock::now();
          auto result = solver.solve(refs);
          auto t1 = std::chrono::steady_clock::now();
          double solve_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

          std::cout << "  [" << std::setw(4) << sname << "]  "
                    << "Status=" << std::setw(16) << std::left
                    << hs_moco::to_string(result.status) << std::right
                    << "  HV=" << std::setw(10) << result.objective_value
                    << "  Nodes=" << std::setw(8) << result.node_count
                    << "  Time=" << std::fixed << std::setprecision(1)
                    << std::setw(9) << solve_ms << " ms"
                    << "  f(x*)=";
          print_vec(result.solution_point);

          // solver.print_statistics(std::cout);

          // Validate.
          if (result.status == hs_moco::SolveStatus::Optimal ||
              result.status == hs_moco::SolveStatus::InfeasibleOrUnbounded) {
            if (result.solution_point == exp_point) {
              std::cout << "  ✔ PASS";
              ++total_pass;
            } else {
              std::cout << "  ✘ FAIL (expected ";
              print_vec(exp_point);
              std::cout << ")";
              ++total_fail;
            }
          } else if (result.status == hs_moco::SolveStatus::TimeLimitReached) {
            std::cout << "  ✘ FAIL (" << hs_moco::to_string(result.status) << ")";
            ++total_fail;
          } else {
            std::cout << "  ✘ FAIL (unexpected status: "
                      << hs_moco::to_string(result.status) << ")";
            ++total_fail;
          }
          std::cout << "\n";
          std::cout.flush();
        }
        std::cout << "\n";
        std::cout.flush();
      }
    }

    // ────────────────────────────────────────────────────────────────────────
    // 5. Summary
    // ────────────────────────────────────────────────────────────────────────
    std::cout << std::string(72, '=') << "\n";
    std::cout << "Summary: " << total_pass << "/" << total_tests << " passed";
    if (total_fail > 0)
      std::cout << ", " << total_fail << " FAILED";
    std::cout << "\n";

    return total_tests > 0 && total_pass == total_tests ? 0 : 1;
  } catch (GRBException& e) {
    std::cerr << "GRBException: " << e.getMessage()
              << " (error code " << e.getErrorCode() << ")\n";
    return 1;
  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << "\n";
    return 1;
  }
}
