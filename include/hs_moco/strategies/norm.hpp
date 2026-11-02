#ifndef HS_MOCO_STRATEGIES_NORM_HPP
#define HS_MOCO_STRATEGIES_NORM_HPP

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace hs_moco {

enum class NormalizationMethod {
  None,
  PrecisionOptimal,  ///< Legacy name for a repository-specific objective-scaling heuristic.
  UpperBound,
  ReferenceGap,
};

namespace detail {

inline bool valid_norm_bounds(const std::vector<double>& bounds) {
  return !bounds.empty() && std::all_of(bounds.begin(), bounds.end(),
      [](double value) { return value > 0 && std::isfinite(value); });
}

// Equalize coordinate bounds at exp(log_volume / M); invalid inputs leave units unchanged.
inline std::vector<double> norm_with_log_volume(const std::vector<double>& bounds,
                                               double log_volume) {
  if (!valid_norm_bounds(bounds) || !std::isfinite(log_volume))
    return std::vector<double>(bounds.size(), 1.0);
  const double scale = std::exp(log_volume / bounds.size());
  if (!(scale > 0) || !std::isfinite(scale))
    return std::vector<double>(bounds.size(), 1.0);
  auto result = bounds;
  for (double& value : result) value /= scale;
  if (!valid_norm_bounds(result)) return std::vector<double>(bounds.size(), 1.0);
  return result;
}

}  // namespace detail

/// Divide by U_i, giving coordinate upper bounds of one.
inline std::vector<double> upper_bound_norm(const std::vector<double>& upper_bound) {
  return detail::norm_with_log_volume(upper_bound, 0.0);
}

/// Repository-specific heuristic; `objective_scale` is a positive objective scale
/// factor, not a claim of optimal scaling or an HS error/gap guarantee.
/// The returned divisors make the modeled objective `objective_scale * HS`.
inline std::vector<double> precision_optimal_norm(const std::vector<double>& upper_bound,
                                                 double objective_scale = 1e-8) {
  if (!detail::valid_norm_bounds(upper_bound) || !(objective_scale > 0) ||
      !std::isfinite(objective_scale))
    return std::vector<double>(upper_bound.size(), 1.0);
  double log_volume = std::log(objective_scale);
  for (double value : upper_bound) log_volume += std::log(value);
  return detail::norm_with_log_volume(upper_bound, log_volume);
}

/// Divide by (U_i - min_r r_i) / target_volume^(1/M), bounding each subset product.
inline std::vector<double> reference_gap_norm(
    const std::vector<double>& upper_bound,
    const std::vector<std::vector<int64_t>>& refs,
    double target_volume = 1e4) {
  if (refs.empty() || !(target_volume > 0) || !std::isfinite(target_volume) ||
      std::any_of(refs.begin(), refs.end(), [&](const auto& ref) {
        return ref.size() != upper_bound.size();
      }))
    return std::vector<double>(upper_bound.size(), 1.0);
  auto gaps = upper_bound;
  for (size_t i = 0; i < gaps.size(); ++i) {
    int64_t minimum = refs.front()[i];
    for (const auto& ref : refs) minimum = std::min(minimum, ref[i]);
    gaps[i] -= static_cast<double>(minimum);
  }
  return detail::norm_with_log_volume(gaps, std::log(target_volume));
}

/// Return positive divisors for objective gaps; references must use the same units as bounds.
inline std::vector<double> normalization_factors(
    NormalizationMethod method,
    const std::vector<double>& upper_bound,
    const std::vector<std::vector<int64_t>>& refs,
    double objective_scale = 1e-8) {
  switch (method) {
    case NormalizationMethod::None:
      return std::vector<double>(upper_bound.size(), 1.0);
    case NormalizationMethod::PrecisionOptimal:
      return precision_optimal_norm(upper_bound, objective_scale);
    case NormalizationMethod::UpperBound:
      return upper_bound_norm(upper_bound);
    case NormalizationMethod::ReferenceGap:
      return reference_gap_norm(upper_bound, refs);
  }
  return std::vector<double>(upper_bound.size(), 1.0);
}

/// A normalization method and its validated objective scale (used by PrecisionOptimal).
class Normalization {
 public:
  Normalization(NormalizationMethod method = NormalizationMethod::None,
                double objective_scale = 1e-8)
      : method_(method), objective_scale_(objective_scale) {
    if (!(objective_scale_ > 0) || !std::isfinite(objective_scale_))
      throw std::invalid_argument("objective_scale must be positive and finite");
  }

  [[nodiscard]] NormalizationMethod method() const noexcept { return method_; }
  [[nodiscard]] double objective_scale() const noexcept { return objective_scale_; }

  [[nodiscard]] std::vector<double> factors(
      const std::vector<double>& upper_bound,
      const std::vector<std::vector<int64_t>>& refs) const {
    return normalization_factors(method_, upper_bound, refs, objective_scale_);
  }

 private:
  NormalizationMethod method_;
  double objective_scale_;
};

}  // namespace hs_moco

#endif  // HS_MOCO_STRATEGIES_NORM_HPP
