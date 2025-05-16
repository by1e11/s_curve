#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "s_curve.h"

namespace py = pybind11;

PYBIND11_MODULE(s_curve, m) {
    py::enum_<Derivative>(m, "Derivative")
        .value("POSITION", Derivative::Position)
        .value("VELOCITY", Derivative::Velocity)
        .value("ACCELERATION", Derivative::Acceleration)
        .value("JERK", Derivative::Jerk);

    py::class_<SCurveConstraints>(m, "SCurveConstraints")
        .def(py::init<>())
        .def_readwrite("max_jerk", &SCurveConstraints::max_jerk)
        .def_readwrite("max_acceleration", &SCurveConstraints::max_acceleration)
        .def_readwrite("max_velocity", &SCurveConstraints::max_velocity);

    py::class_<SCurveTimeIntervals>(m, "SCurveTimeIntervals")
        .def(py::init<>())
        .def_readwrite("t_j1", &SCurveTimeIntervals::t_j1)
        .def_readwrite("t_j2", &SCurveTimeIntervals::t_j2)
        .def_readwrite("t_a", &SCurveTimeIntervals::t_a)
        .def_readwrite("t_v", &SCurveTimeIntervals::t_v)
        .def_readwrite("t_d", &SCurveTimeIntervals::t_d)
        .def("total_duration", &SCurveTimeIntervals::total_duration)
        .def("is_max_acceleration_not_reached", &SCurveTimeIntervals::is_max_acceleration_not_reached);

    py::class_<SCurveStartConditions>(m, "SCurveStartConditions")
        .def(py::init<>())
        .def_readwrite("q0", &SCurveStartConditions::q0)
        .def_readwrite("q1", &SCurveStartConditions::q1)
        .def_readwrite("v0", &SCurveStartConditions::v0)
        .def_readwrite("v1", &SCurveStartConditions::v1)
        .def("h", &SCurveStartConditions::h)
        .def("dir", &SCurveStartConditions::dir);

    py::class_<SCurveInput>(m, "SCurveInput")
        .def(py::init<>())
        .def_readwrite("constraints", &SCurveInput::constraints)
        .def_readwrite("start_conditions", &SCurveInput::start_conditions)
        .def("calc_intervals", &SCurveInput::calc_intervals)
        .def("is_trajectory_feasible", &SCurveInput::is_trajectory_feasible);

    py::class_<SCurveParameters>(m, "SCurveParameters")
        .def(py::init<>())
        .def_readwrite("time_intervals", &SCurveParameters::time_intervals)
        .def_readwrite("j_max", &SCurveParameters::j_max)
        .def_readwrite("j_min", &SCurveParameters::j_min)
        .def_readwrite("a_lim_a", &SCurveParameters::a_lim_a)
        .def_readwrite("a_lim_d", &SCurveParameters::a_lim_d)
        .def_readwrite("v_lim", &SCurveParameters::v_lim)
        .def_readwrite("conditions", &SCurveParameters::conditions)
        .def_static("create", &SCurveParameters::create);

    m.def("s_curve_generator", &s_curve_generator);
}