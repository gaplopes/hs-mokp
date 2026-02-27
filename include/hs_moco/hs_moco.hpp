#ifndef HS_MOCO_HPP
#define HS_MOCO_HPP

/**
 * @file hs_moco.hpp
 * @brief Umbrella header for the HS-MOCO hypervolume scalarization library.
 *
 * @code
 *   #include <hs_moco/hs_moco.hpp>
 *   - Gives you: hs_moco::MOKPSolver, hs_moco::MOKProblem,
 *   -            hs_moco::MIQPStrategy, hs_moco::LogStrategy,
 *   -            hs_moco::MILPStrategy, hs_moco::SolverOptions,
 *   -            hs_moco::SolveResult, hs_moco::SolveStatus
 * @endcode
 */

// Core
#include "hs_moco/core/types.hpp"
#include "hs_moco/problems/mokp.hpp"
#include "hs_moco/core/strategy.hpp"
#include "hs_moco/core/solver.hpp"

// Strategies
#include "hs_moco/strategies/miqp.hpp"
#include "hs_moco/strategies/log.hpp"
#include "hs_moco/strategies/milp.hpp"

// Concrete problems
#include "hs_moco/problems/mokp.hpp"

#endif  // HS_MOCO_HPP
