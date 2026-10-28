#ifndef HS_MOCO_CORE_REFERENCES_HPP
#define HS_MOCO_CORE_REFERENCES_HPP

#include <algorithm>
#include <cstdint>
#include <map>
#include <utility>
#include <vector>

namespace hs_moco::detail {

/**
 * @brief Aggregate inclusion-exclusion weights by component-wise maximum.
 *
 * Each key is a box's lower corner; its value is the sum of subset signs.
 * All boxes share the upper corner f(x), so an intersection's lower corner
 * is the component-wise maximum of its inputs.
 * Adding a reference includes its box and subtracts its intersection with
 * the existing union. Identical corners cancel before the next iteration,
 * so no reference-subset bitmask enumeration is needed.
 * References must have matching dimensions and satisfy the solver's count limit.
 */
inline std::map<std::vector<int64_t>, int64_t> inclusion_exclusion_corners(
    const std::vector<std::vector<int64_t>>& refs) {
  std::map<std::vector<int64_t>, int64_t> corners;
  for (const auto& ref : refs) {
    auto next = corners;
    ++next[ref];
    for (const auto& [corner, weight] : corners) {
      auto maximum = corner;
      for (size_t d = 0; d < maximum.size(); ++d) maximum[d] = std::max(maximum[d], ref[d]);
      next[maximum] -= weight;
    }
    std::erase_if(next, [](const auto& entry) { return entry.second == 0; });
    corners = std::move(next);
  }
  return corners;
}

}  // namespace hs_moco::detail

#endif  // HS_MOCO_CORE_REFERENCES_HPP
