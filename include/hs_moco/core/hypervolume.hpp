#ifndef HS_MOCO_CORE_HYPERVOLUME_HPP
#define HS_MOCO_CORE_HYPERVOLUME_HPP

#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

#include "hs_moco/core/references.hpp"

namespace hs_moco::detail {

inline int64_t checked_hs_add(int64_t a, int64_t b) {
  constexpr int64_t lo = std::numeric_limits<int64_t>::min();
  constexpr int64_t hi = std::numeric_limits<int64_t>::max();
  if ((b > 0 && a > hi - b) || (b < 0 && a < lo - b))
    throw std::overflow_error("exact HS addition overflow");
  return a + b;
}

inline int64_t checked_hs_subtract(int64_t a, int64_t b) {
  constexpr int64_t lo = std::numeric_limits<int64_t>::min();
  constexpr int64_t hi = std::numeric_limits<int64_t>::max();
  if ((b < 0 && a > hi + b) || (b > 0 && a < lo + b))
    throw std::overflow_error("exact HS subtraction overflow");
  return a - b;
}

inline int64_t checked_hs_multiply(int64_t a, int64_t b) {
  constexpr int64_t lo = std::numeric_limits<int64_t>::min();
  constexpr int64_t hi = std::numeric_limits<int64_t>::max();
  if ((a > 0 && ((b > 0 && a > hi / b) || (b < 0 && b < lo / a))) ||
      (a < 0 && ((b > 0 && a < lo / b) || (b < 0 && a < hi / b))))
    throw std::overflow_error("exact HS multiplication overflow");
  return a * b;
}

/**
 * @brief Compute exact hypervolume for a point strictly above every reference.
 *
 * Requires a nonempty point and reference set of matching dimensions, with
 * strict domination in every coordinate. The reference count is capped below
 * int64's signed value-bit count; gaps, products, and weighted sums are checked.
 * A representable final hypervolume can still fail if an intermediate term
 * overflows before cancellation.
 *
 * shortcut: intermediate overflow may reject a representable result; use wider exact arithmetic if that must be supported.
 */
inline int64_t hypervolume(const std::vector<int64_t>& point,
                           const std::vector<std::vector<int64_t>>& refs) {
  if (point.empty() || refs.empty() ||
      refs.size() >= static_cast<std::size_t>(std::numeric_limits<int64_t>::digits))
    throw std::invalid_argument("invalid exact HS input");
  for (const auto& ref : refs) {
    if (ref.size() != point.size())
      throw std::invalid_argument("reference objective dimension mismatch");
    for (std::size_t d = 0; d < point.size(); ++d)
      if (point[d] <= ref[d])
        throw std::invalid_argument("point does not strictly dominate references");
  }

  int64_t volume = 0;
  for (const auto& [corner, weight] : inclusion_exclusion_corners(refs)) {
    int64_t product = 1;
    for (std::size_t d = 0; d < point.size(); ++d)
      product = checked_hs_multiply(product, checked_hs_subtract(point[d], corner[d]));
    volume = checked_hs_add(volume, checked_hs_multiply(weight, product));
  }
  return volume;
}

}  // namespace hs_moco::detail

#endif  // HS_MOCO_CORE_HYPERVOLUME_HPP
