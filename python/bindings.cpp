#include <bottleneck/core.hpp>
#include <bottleneck/wasserstein.hpp>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <cstddef>
#include <memory>
#include <span>
#include <vector>

namespace py = pybind11;

namespace {

struct NativePreparedDiagram {
  explicit NativePreparedDiagram(const py::array_t<double, py::array::c_style>& values)
      : point_count(static_cast<std::size_t>(values.shape(0))),
        prepared(make_diagram(values)) {}

  static bottleneck::Diagram make_diagram(
      const py::array_t<double, py::array::c_style>& values) {
    const auto view = values.unchecked<2>();
    bottleneck::Diagram diagram;
    diagram.reserve(static_cast<std::size_t>(view.shape(0)));
    for (py::ssize_t row = 0; row < view.shape(0); ++row) {
      diagram.push_back({view(row, 0), view(row, 1)});
    }
    return diagram;
  }

  [[nodiscard]] std::size_t finite_point_count() const noexcept {
    return prepared.finite_points().size();
  }

  std::size_t point_count;
  bottleneck::PreparedDiagram prepared;
};

std::vector<bottleneck::PreparedDiagram> copy_prepared(
    const std::vector<std::shared_ptr<NativePreparedDiagram>>& diagrams) {
  std::vector<bottleneck::PreparedDiagram> result;
  result.reserve(diagrams.size());
  for (const auto& diagram : diagrams) {
    result.push_back(diagram->prepared);
  }
  return result;
}

void require_output(const py::array_t<double, py::array::c_style>& output,
                    std::size_t size) {
  if (output.ndim() != 1 || static_cast<std::size_t>(output.shape(0)) != size) {
    throw py::value_error("output shape must match the number of diagrams");
  }
  if (!output.writeable()) {
    throw py::value_error("output array must be writeable");
  }
}

}  // namespace

PYBIND11_MODULE(_core, module) {
  module.doc() = "Private native extension for topp";

  py::class_<NativePreparedDiagram, std::shared_ptr<NativePreparedDiagram>>(
      module, "PreparedDiagram")
      .def_property_readonly("n_points",
                             [](const NativePreparedDiagram& value) {
                               return value.point_count;
                             })
      .def_property_readonly("n_finite_points",
                             &NativePreparedDiagram::finite_point_count)
      .def("__len__", [](const NativePreparedDiagram& value) {
        return value.point_count;
      });

  module.def("prepare_diagram",
             [](const py::array_t<double, py::array::c_style>& values) {
               return std::make_shared<NativePreparedDiagram>(values);
             });

  module.def("bottleneck_distance",
             [](const NativePreparedDiagram& first,
                const NativePreparedDiagram& second) {
               py::gil_scoped_release release;
               return bottleneck::bottleneck_distance(first.prepared,
                                                      second.prepared);
             });

  module.def("wasserstein_distance",
             [](const NativePreparedDiagram& first,
                const NativePreparedDiagram& second, int metric) {
               bottleneck::WassersteinConfig config;
               config.metric = metric == 0 ? bottleneck::WassersteinMetric::w1_linf
                                           : bottleneck::WassersteinMetric::w2_l2;
               py::gil_scoped_release release;
               return bottleneck::wasserstein_distance(first.prepared,
                                                       second.prepared, config);
             });

  module.def("bottleneck_within",
             [](const NativePreparedDiagram& first,
                const NativePreparedDiagram& second, double threshold) {
               py::gil_scoped_release release;
               return bottleneck::bottleneck_within(first.prepared, second.prepared,
                                                    threshold);
             });

  module.def("bottleneck_distances",
             [](const NativePreparedDiagram& query,
                const std::vector<std::shared_ptr<NativePreparedDiagram>>& diagrams,
                py::array_t<double, py::array::c_style> output) {
               require_output(output, diagrams.size());
               auto prepared = copy_prepared(diagrams);
               auto buffer = output.mutable_unchecked<1>();
               {
                 py::gil_scoped_release release;
                 bottleneck::bottleneck_distances(
                     query.prepared, prepared,
                     std::span<double>(buffer.mutable_data(0), diagrams.size()));
               }
             });

  module.def("wasserstein_distances",
             [](const NativePreparedDiagram& query,
                const std::vector<std::shared_ptr<NativePreparedDiagram>>& diagrams,
                int metric, py::array_t<double, py::array::c_style> output) {
               require_output(output, diagrams.size());
               auto prepared = copy_prepared(diagrams);
               auto buffer = output.mutable_unchecked<1>();
               bottleneck::WassersteinConfig config;
               config.metric = metric == 0 ? bottleneck::WassersteinMetric::w1_linf
                                           : bottleneck::WassersteinMetric::w2_l2;
               bottleneck::WassersteinWorkspace workspace;
               {
                 py::gil_scoped_release release;
                 bottleneck::wasserstein_distances(
                     query.prepared, prepared,
                     std::span<double>(buffer.mutable_data(0), diagrams.size()),
                     workspace, config);
               }
             });
}
