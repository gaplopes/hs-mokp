# HS-MOKP

> Header-only C++20 library for **Hypervolume Scalarization** in
> the Multiobjective Knapsack Problem (MOKP).

## Background

In multiobjective optimization we seek solutions that are Pareto-optimal — no
objective can be improved without degrading another.  The **hypervolume
indicator** measures the volume of objective space that is dominated by a set of
solutions and bounded by a reference point.

**Hypervolume Scalarization (HS)** turns this set-level indicator into a
scalar objective for a *single* solution.  Given a set of reference points
$R = \{r^1, \ldots, r^m\}$ in objective space, the HS of a feasible
solution $x$ is the hypervolume of $f(x)$ w.r.t. the reference points, which can be
computed via inclusion–exclusion over all non-empty subsets $S \subseteq R$:

$$
\text{HS}(x, R) = \sum_{\emptyset \neq S \subseteq R} (-1)^{|S|+1}
\prod_{i=1}^{m} \bigl(f_i(x) - \max_{r \in S}\, r_i\bigr)
$$

HS-MOKP encodes this objective as a Mixed-Integer Program (MIP)
specifically parameterized for the Multiobjective Knapsack Problem (MOKP),
and provides three formulation strategies, with optional corner reduction in MIQP and Log,
trading off model size and numerical behaviour directly tailored for Gurobi formulation patterns.
See the [Strategies](#strategies) section below for details.

> [!NOTE]
> **Solver backend.** This library requires **[Gurobi Optimizer](https://www.gurobi.com/)** 
> to function, directly relying on its C++ API capabilities to interpret linear 
> constraints (MILP), nonlinear product constraints (`addGenConstrNL`), and non-linear function
> features (`FuncNonLinear`) like logs and exponential operations utilized in the formulations.

> [!NOTE]
> **Maximization only.** The HS formulation as implemented assumes a
> **maximization** problem (all objectives are to be maximized). For
> minimization problems, this can be easily adapted by negating the
> objective values, i.e., replacing $f_d(x)$ with $-f_d(x)$ and adjusting
> the reference points accordingly.

## Strategies

MIQP and Log with corner reduction enabled, and MILP, aggregate inclusion–exclusion
weights by the component-wise maximum reference point. Identical corners combine,
and zero weights disappear. MIQP and Log default to subset enumeration:

$$
\text{HS}(x,R) = \sum_c w_c \prod_d (f_d(x)-c_d).
$$

This preprocessing builds intersections incrementally, without enumerating all
$2^K-1$ reference subsets. Duplicate or contained boxes add no retained terms.
The number of distinct corners can still grow rapidly with many objectives.

Every reference $r$ defines a box $[r,f(x)]$, since the solver requires $f(x)$
to strictly dominate all references. Boxes share an upper corner, so their
intersection has lower corner $\max(r^1,r^2)$, computed coordinate by coordinate.
Different subsets with the same lower corner define the same product for every
feasible $x$; adding their signs therefore preserves the entire objective.

For a three-objective example, take the mutually nondominated references
$A=(1,1,0)$, $B=(1,0,1)$, and $C=(0,1,1)$:

| Subsets | Maximum corner | Combined weight |
|---|---|---:|
| $A$ | $(1,1,0)$ | $+1$ |
| $B$ | $(1,0,1)$ | $+1$ |
| $C$ | $(0,1,1)$ | $+1$ |
| $AB$, $AC$, $BC$, $ABC$ | $(1,1,1)$ | $-1-1-1+1=-2$ |

Seven subset terms reduce to four nonzero corner terms. At $f(x)=(2,2,2)$,
each singleton box has volume $2$ and every intersection has volume $1$:
both expressions yield $2+2+2-2\cdot1=4$.

`inclusion_exclusion_corners` builds this map incrementally. When adding a new
box, it retains the old map, adds the box with weight $+1$, then subtracts each
weighted intersection with the old map. This is the identity
$\operatorname{vol}(U\cup B)=\operatorname{vol}(U)+\operatorname{vol}(B)-\operatorname{vol}(U\cap B)$.
Weights are signed integer sums, not necessarily just $+1$ or $-1$.

### 1. MIQP — Original Nonlinear Product

`MIQPStrategy` in `strategies/miqp.hpp` preserves the original formulation.
It enumerates all $2^K-1$ nonempty reference subsets and builds a product for
each one. For $M\geq2$, each subset adds $M+1$ variables, $M$ linear constraints,
and one nonlinear general constraint. Its balanced binary-tree helper is also
retained unchanged.

Pass `true` as the second constructor argument to aggregate identical corners
before building products:

```cpp
auto aggregated_strategy = std::make_unique<hs_moco::MIQPStrategy>(
    hs_moco::NormalizationMethod::None, true);
```

This opt-in mode uses the same unshifted factors and product helper as the
original formulation. Each retained corner contributes its full signed integer
weight. With $C$ nonzero corners and $M\geq2$, it adds $(M+1)C$ variables,
$MC$ linear equalities, and $C$ nonlinear constraints. The default remains
`false`, preserving subset enumeration. The MOKP example compares both modes
on the same references and reports term counts, model sizes, and solve times.

Both MIQP modes default to `NormalizationMethod::None`, leaving factors unscaled.
Per-objective Dantzig upper bounds are computed once during solver construction
from the fractional knapsack relaxation. Selecting another method uses these
bounds for normalization. With `PrecisionOptimal` and `abs_tol = 1e-8`,
the modeled objective is `1e-8 * HS`; returned results restore the original HS units.
The normalization tolerance is at least ten times the model's effective
`FeasibilityTol` (normally `1e-9`).
`PrecisionOptimal` leaves the formulation unscaled for nonpositive bounds.

Upper bounds alone do not enable normalization. All methods return HS in its original units.

`strategies/norm.hpp` provides `precision_optimal_norm`, `upper_bound_norm`,
and `reference_gap_norm`. Both MIQP modes accept the same normalization method:

```cpp
auto strategy = std::make_unique<hs_moco::MIQPStrategy>(
    hs_moco::NormalizationMethod::ReferenceGap);
auto aggregated_strategy = std::make_unique<hs_moco::MIQPStrategy>(
    hs_moco::NormalizationMethod::ReferenceGap, true);
```

`UpperBound` divides by each upper bound, and `ReferenceGap` scales the
reference-gap bounding box to volume `1e4`. These choices can produce different
numerical errors despite preserving the exact mathematical maximizer.

### 2.Logarithm transformation and Log-Sum-Exp (LSE) reformulation

Takes the logarithm of the HS product, converting products into sums:

$$
\log \prod_d (f_d - r_d) = \sum_d \log(f_d - r_d)
$$

For $k = 1$ reference point, this is the complete formulation. By default,
`LogStrategy` enumerates every nonempty reference subset, as MIQP does. Pass
`true` to enable the same corner-reduction approach:

```cpp
auto strategy = std::make_unique<hs_moco::LogStrategy>();
auto aggregated_strategy = std::make_unique<hs_moco::LogStrategy>(true);
```

The default `reduce_corners` value is `false`. Both modes use the same log/LSE
construction and shared gap/log pairs. With reduction enabled, multiple references
that reduce to one retained corner use the single-reference formulation.
The MOKP example compares both Log modes alongside both MIQP modes.

For multiple retained corners, each sum of logs includes $\log|w_c|$ to account
for the aggregated weight. Inclusion–exclusion requires computing
$\log(A - B)$ where $A$ and $B$ are sums of exponentials. This uses the
numerically stable identity:

$$
\log(A - B) = \underbrace{\log A}_{L_\text{pos}} + \log\!\bigl(1 - \exp(\underbrace{\log B}_{L_\text{neg}} - L_\text{pos})\bigr)
$$

Each $\log A$ and $\log B$ is computed via the **max-shift LSE trick**:

$$
\text{LSE}(s_0, \ldots, s_n) = \max(s) + \log\!\Bigl(\sum_i \exp(s_i - \max(s))\Bigr)
$$

Gap variables and their logarithms are shared by objective index and reference
coordinate within each solve. There are at most $MK$ distinct gap/log pairs for
$M$ objectives and $K$ references, regardless of the retained corner count.
A singleton LSE returns its input directly. The final objective maximizes
$L_\text{pos}+\log(1-\exp(L_\text{neg}-L_\text{pos}))$ without another resultant
variable or defining equality.

### 3. MILP — McCormick Linearization

Expands the product $\prod_d (f_d(x) - r_d)$ into a **multilinear
polynomial** over individual decision variables $x_j$, then linearizes every
product of binaries via McCormick envelopes:

$$
y_{ij} \approx x_i \cdot x_j \quad \Longleftrightarrow \quad
y_{ij} \leq x_i, \;
y_{ij} \leq x_j, \;
y_{ij} \geq x_i + x_j - 1
$$

Since $x_j \in \{0,1\} \Rightarrow x_j^2 = x_j$, degenerate monomials
collapse (e.g., $x_0 \cdot x_0 \cdot x_1 = x_0 \cdot x_1$), so only
$\binom{N}{m}$ **strict** combinations are needed at each degree $m$.

Higher-degree products use $m+1$ constraints: $z \leq x_j$ for each factor and
$z \geq \sum_j x_j - m + 1$. This remains exact for binary decisions, but its
continuous relaxation can be weaker than the recursive formulation. Product
variables and constraints are retained across solves with different references.

> [!TIP]
> Each strategy header file
> ([`miqp.hpp`](include/hs_moco/strategies/miqp.hpp),
> [`log.hpp`](include/hs_moco/strategies/log.hpp),
> [`milp.hpp`](include/hs_moco/strategies/milp.hpp))
> contains detailed documentation on the mathematical formulation,
> model size analysis, and implementation steps. Refer to the file-level
> docstrings for a comprehensive deep-dive beyond the overview above.

## Features

- **Built exclusively on top of the Gurobi C++ API.**
- **Three MIP formulation strategies for MOKP instances:**
  1. MIQP (nonlinear product, with optional corner reduction)
  2. Logarithm transformation and Log-Sum-Exp (LSE) reformulation (with optional corner reduction)
  3. MILP (McCormick linearization)
- **Build once, solve many times:** separate model construction from solving, allowing efficient re-use across multiple reference-point sets
- **Pure header-only:** no compiled library, just `#include <hs_moco/hs_moco.hpp>`

## Requirements

### Core library (`hs_moco`)

| Dependency | Version | Notes |
|------------|---------|-------|
| C++ compiler | C++20 support | GCC ≥ 10, Clang ≥ 13, Microsoft Visual C++ ≥ 19.29 |
| CMake | ≥ 3.21 | |
| [Gurobi Optimizer](https://www.gurobi.com/) | ≥ 12 | Set `GUROBI_HOME` to your installation directory |

Set `GUROBI_HOME` to your Gurobi installation directory (e.g. `/opt/gurobi1203/linux64`).

### Example programs only

**Note that:** Example programs (and their dependencies) are configured and built **only when `hs-moco` is the top-level project**. If used as a subproject (e.g., via `add_subdirectory` or `FetchContent`), examples are skipped by default.

| Dependency | Version | Notes |
|------------|---------|-------|
| [`mooutils`](https://github.com/adbjesus/mooutils) | latest | Needed only for examples; auto-fetched via `FetchContent` if not found |

## Building

```bash
# Configure (uses CMake presets)
cmake --preset release

# Build
cmake --build --preset release

# Or manually:
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build .
```

The example executable `example_mokp` will be in `build/release/examples/`.

### Running the MOKP Example

The MOKP example uses MOKP instances from
[`gaplopes/mobkp-instances`](https://github.com/gaplopes/mobkp-instances),
which are downloaded automatically via CMake `FetchContent` at configure time.
No manual setup is required. You can also provide any compatible `.in` instance
file.

```bash
# Use the bundled default instance (max_k = 5)
./build/release/examples/example_mokp

# Use a custom instance file
./build/release/examples/example_mokp /path/to/instance.in

# Use a custom instance file and a custom maximum number of reference points
./build/release/examples/example_mokp /path/to/instance.in 3
```

The example compares MIQP, MIQP-Aggregated, Log, and Log-Aggregated for
reference-point set sizes \(k = 1, \ldots, K\). For each \(k\), it generates 30
(default) random non-dominated reference-point sets (each of size \(k\)),
measures wall-clock time, and validates each result against exhaustive
enumeration over the known Pareto front.

#### Batch Testing

The script `examples/test_mokp.sh` automates running `example_mokp` across all
bundled benchmark instances and collects pass/fail results:

```bash
# Quick sanity check (small instances only, max_k = 3)
./examples/test_mokp.sh --quick

# Full run over all instances (max_k = 5, the default)
./examples/test_mokp.sh --all

# Custom options
./examples/test_mokp.sh --max-k 4 --dims 3 --timeout 120
```

Run `./examples/test_mokp.sh --help` for the complete list of options.

The results are saved in `examples/test_results/` (inside the build directory) with timestamps, and a summary table is printed to the console.

## Usage

```cpp
#include <hs_moco/hs_moco.hpp>

// 1. Define items
std::vector<hs_moco::MOKPItem> items = {
    {0, 2, {8, 3, 5}},
    {1, 3, {4, 7, 6}},
    // ...
};

// 2. Create problem + strategy + solver
auto problem  = std::make_unique<hs_moco::MOKProblem>(std::move(items), /*W=*/10);
auto strategy = std::make_unique<hs_moco::MIQPStrategy>();

hs_moco::SolverOptions opts;
opts.time_limit = 60.0;

hs_moco::MOKPSolver solver(std::move(problem), std::move(strategy), opts);

// 3. Solve repeatedly with different reference points
auto result = solver.solve({{0, 0, 0}});

if (result.status == hs_moco::SolveStatus::Optimal) {
    std::cout << "HS = " << result.objective_value << "\n";
}
```

Set additional Gurobi parameters through `SolverOptions::gurobi_params` before
constructing the solver. Both names and values are strings:

```cpp
hs_moco::SolverOptions opts;
opts.gurobi_params = {
    {"Threads", "8"},
    {"MIPFocus", "1"},
    {"TimeLimit", "60"},
    {"MIPGap", "1e-4"},
    {"LogFile", "gurobi.log"},
};
```

Entries override existing option fields and strategy defaults during construction.
When building its objective, MIQP sets `NonConvex=2` and Log sets `FuncNonlinear=1`.
These strategy requirements override conflicting `gurobi_params` entries.
On every solve, the Log strategy sets the model's `FeasibilityTol` to `1e-9` for
`k=1` and `1e-6` for `k>1`, regardless of `feasibility_tol` or `gurobi_params`.

Parameters apply before environment startup and model creation, so startup
parameters such as `ThreadLimit` are supported. Gurobi validates names, values,
and applicability; invalid entries throw `GRBException` during construction.
Parameter compatibility with the chosen strategy remains Gurobi's responsibility.

The existing `SolveStatus` mapping is unchanged. Additional stopping criteria can
produce termination statuses currently mapped to `SolveStatus::Error`.


## Consuming in Your CMake Project

Gurobi 12 or newer must be findable at configure time. When examples are enabled,
`mooutils` is automatically fetched via `FetchContent` if it is not already installed.

### Via `add_subdirectory`

```cmake
add_subdirectory(path/to/hs-moco)
target_link_libraries(your_target PRIVATE hs_moco::hs_moco)
```

### Via `FetchContent`

```cmake
include(FetchContent)
FetchContent_Declare(
  hs_moco
  GIT_REPOSITORY https://github.com/gaplopes/hs-moco.git
  GIT_TAG        main
)
FetchContent_MakeAvailable(hs_moco)
target_link_libraries(your_target PRIVATE hs_moco::hs_moco)
```

### Via `find_package` (after install)

```bash
cmake --install build/release --prefix /usr/local
```

```cmake
find_package(hs_moco REQUIRED)
target_link_libraries(your_target PRIVATE hs_moco::hs_moco)
```

## Project Structure

```
hs-moco/
├── CMakeLists.txt
├── CMakePresets.json
├── cmake/
│   ├── hs_moco-config.cmake.in
│   └── modules/
│       └── FindGUROBI.cmake
├── include/
│   └── hs_moco/
│       ├── hs_moco.hpp              ← umbrella header
│       ├── core/
│       │   ├── types.hpp            ← SolveStatus, SolverOptions, SolveResult
│       │   ├── strategy.hpp         ← HSStrategy, AuxiliaryState, StrategyContext
│       │   └── solver.hpp           ← MOKPSolver orchestrator
│       ├── strategies/
│       │   ├── miqp.hpp             ← nonlinear product (addGenConstrNL)
│       │   ├── log.hpp              ← log-sum-exp (FuncNonlinear=1)
│       │   └── milp.hpp             ← McCormick linearization
│       └── problems/
│           └── mokp.hpp             ← Multiobjective Knapsack Problem
├── examples/
│   ├── CMakeLists.txt
│   ├── mokp.cpp                 ← MOKP strategy-comparison example
│   └── test_mokp.sh             ← batch test script across instances
└── README.md
```

## Contributing

Contributions are welcome!  If you have bug fixes, new problem implementations,
additional strategies, or other improvements, please fork the repository and
open a pull request.  For major changes, consider opening an issue first to
discuss the approach.

## Citation

```bibtex
@software{lopes2026hsmokp,
  author       = {Lopes, Gon{\c{c}}alo},
  title        = {{HS-MOKP}: Hypervolume Scalarization for the Multiobjective
                  Knapsack Problem},
  year         = {2026},
  url          = {https://github.com/gaplopes/hs-mokp}
}
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

For any questions or issues, please open an issue on the repository or contact the author at galopes@dei.uc.pt or via GitHub or LinkedIn (see profile for contact information).
