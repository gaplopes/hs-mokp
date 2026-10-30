#include <cstdint>
#include <exception>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "hs_moco/hs_moco.hpp"

namespace py = pybind11;

PYBIND11_MODULE(hs_moco, m) {
  using namespace hs_moco;

  m.doc() = "Hypervolume scalarization solvers for the multiobjective knapsack problem";

  py::register_local_exception_translator([](std::exception_ptr error) {
    try {
      if (error) std::rethrow_exception(error);
    } catch (const GRBException& e) {
      const auto message = "Gurobi error " + std::to_string(e.getErrorCode()) +
                           ": " + e.getMessage();
      PyErr_SetString(PyExc_RuntimeError, message.c_str());
    }
  });

  py::enum_<SolveStatus>(m, "SolveStatus")
      .value("Optimal", SolveStatus::Optimal)
      .value("TimeLimitReached", SolveStatus::TimeLimitReached)
      .value("InfeasibleOrUnbounded", SolveStatus::InfeasibleOrUnbounded)
      .value("NodeLimitReached", SolveStatus::NodeLimitReached)
      .value("Error", SolveStatus::Error);

  py::enum_<NormalizationMethod>(m, "NormalizationMethod")
      .value("None_", NormalizationMethod::None)
      .value("PrecisionOptimal", NormalizationMethod::PrecisionOptimal)
      .value("UpperBound", NormalizationMethod::UpperBound)
      .value("ReferenceGap", NormalizationMethod::ReferenceGap);

  py::class_<MOKPItem>(m, "MOKPItem")
      .def(py::init<int32_t, int64_t, std::vector<int64_t>>(),
           py::arg("idx"), py::arg("weight"), py::arg("values"))
      .def_readwrite("idx", &MOKPItem::idx)
      .def_readwrite("weight", &MOKPItem::weight)
      .def_readwrite("values", &MOKPItem::values);

  py::class_<SolverOptions>(m, "SolverOptions")
      .def(py::init<>())
      .def_readwrite("time_limit", &SolverOptions::time_limit)
      .def_readwrite("mip_gap", &SolverOptions::mip_gap)
      .def_readwrite("threads", &SolverOptions::threads)
      .def_readwrite("numeric_focus", &SolverOptions::numeric_focus)
      .def_readwrite("console_log", &SolverOptions::console_log)
      .def_readwrite("feasibility_tol", &SolverOptions::feasibility_tol)
      .def_readwrite("gurobi_params", &SolverOptions::gurobi_params);

  py::class_<SolveResult>(m, "SolveResult")
      .def_readonly("status", &SolveResult::status)
      .def_readonly("objective_value", &SolveResult::objective_value)
      .def_readonly("node_count", &SolveResult::node_count)
      .def_readonly("solution_point", &SolveResult::solution_point)
      .def_readonly("decision_values", &SolveResult::decision_values);

  py::class_<MOKPSolver>(m, "MOKPSolver")
      .def(py::init([](std::vector<MOKPItem> items, int64_t capacity,
                       const std::string& strategy, SolverOptions options,
                       bool reduce_corners, NormalizationMethod normalization) {
        if (!items.empty() && items.front().values.empty())
          throw std::invalid_argument("items must have at least one objective");

        std::unique_ptr<HSStrategy> selected_strategy;
        if (strategy == "miqp") {
          selected_strategy = std::make_unique<MIQPStrategy>(normalization, reduce_corners);
        } else {
          if (normalization != NormalizationMethod::None)
            throw std::invalid_argument("normalization is only supported by 'miqp'");
          if (strategy == "log") {
            selected_strategy = std::make_unique<LogStrategy>(reduce_corners);
          } else if (strategy == "milp") {
            if (reduce_corners)
              throw std::invalid_argument("MILP always reduces corners; omit reduce_corners");
            selected_strategy = std::make_unique<MILPStrategy>();
          } else {
            throw std::invalid_argument("strategy must be 'miqp', 'log', or 'milp'");
          }
        }

        auto problem = std::make_unique<MOKProblem>(std::move(items), capacity);
        return std::make_unique<MOKPSolver>(std::move(problem),
                                            std::move(selected_strategy), options);
      }), py::arg("items"), py::arg("capacity"), py::kw_only(),
          py::arg("strategy") = "miqp", py::arg("options") = SolverOptions{},
          py::arg("reduce_corners") = false,
          py::arg("normalization") = NormalizationMethod::None,
          "Build a reusable solver. Items and options are copied at construction.")
      .def("solve", &MOKPSolver::solve, py::arg("refs"),
           py::call_guard<py::gil_scoped_release>(),
           "Solve for integer reference points. Do not call concurrently on the same solver.")
      .def_property_readonly("num_objectives", &MOKPSolver::num_objectives)
      .def_property_readonly("num_variables", &MOKPSolver::num_variables)
      .def_property_readonly("strategy_name", &MOKPSolver::strategy_name);
}
