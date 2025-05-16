#include "s_curve.h"

#include <cmath>
#include <functional>

double SCurveTimeIntervals::total_duration() const {
    return t_a + t_d + t_v;
}

bool SCurveTimeIntervals::is_max_acceleration_not_reached() const {
    return t_a < 2. * t_j1 || t_d < 2. * t_j2;
}

double SCurveStartConditions::h() const {
    return std::abs(q1 - q0);
}

double SCurveStartConditions::dir() const {
    return q1 < q0 ? -1.0 : 1.0;
}

SCurveParameters SCurveParameters::create(const SCurveTimeIntervals& times, const SCurveInput& p) {
    double a_lim_a = p.constraints.max_jerk * times.t_j1;
    double a_lim_d = -p.constraints.max_jerk * times.t_j2;
    double v_lim = p.start_conditions.dir() * p.start_conditions.v0 + 
                   (times.t_a - times.t_j1) * a_lim_a;
    
    return SCurveParameters{
        times,
        p.constraints.max_jerk,
        -p.constraints.max_jerk,
        a_lim_a,
        a_lim_d,
        v_lim,
        p.start_conditions
    };
}

SCurveTimeIntervals SCurveInput::calc_intervals() const {
    return calc_times_case_1();
}

bool SCurveInput::is_trajectory_feasible() const {
    double t_j_star = std::min(
        std::sqrt(
            std::abs(start_conditions.v1 - start_conditions.v0) / 
            constraints.max_jerk
        ),
        constraints.max_acceleration / constraints.max_jerk
    );
    
    if (std::abs(t_j_star - constraints.max_acceleration / constraints.max_jerk) < 1e-10) {
        return start_conditions.h() > 
            0.5 * (start_conditions.v1 + start_conditions.v0) * 
            (t_j_star + std::abs(start_conditions.v1 - start_conditions.v0) / 
            constraints.max_acceleration);
    }
    
    if (t_j_star < constraints.max_acceleration / constraints.max_jerk) {
        return start_conditions.h() > 
            t_j_star * (start_conditions.v0 + start_conditions.v1);
    }
    
    return false;
}

bool SCurveInput::is_a_max_not_reached() const {
    return (constraints.max_velocity - start_conditions.v0) * constraints.max_jerk 
            < std::pow(constraints.max_acceleration, 2);
}

bool SCurveInput::is_a_min_not_reached() const {
    return (constraints.max_velocity - start_conditions.v1) * constraints.max_jerk 
            < std::pow(constraints.max_acceleration, 2);
}

SCurveTimeIntervals SCurveInput::calc_times_case_1() const {
    SCurveTimeIntervals times;
    auto new_input = *this;
    double dir = start_conditions.dir();
    
    if (is_a_max_not_reached()) {
        times.t_j1 = std::sqrt(
            (constraints.max_velocity - start_conditions.v0) / constraints.max_jerk
        );
        times.t_a = 2. * times.t_j1;
    } else {
        times.t_j1 = constraints.max_acceleration / constraints.max_jerk;
        times.t_a = times.t_j1 + 
                    (constraints.max_velocity - dir * start_conditions.v0) / constraints.max_acceleration;
    }

    if (is_a_min_not_reached()) {
        times.t_j2 = std::sqrt(
            (constraints.max_velocity - start_conditions.v1) / constraints.max_jerk
        );
        times.t_d = 2. * times.t_j2;
    } else {
        times.t_j2 = constraints.max_acceleration / constraints.max_jerk;
        times.t_d = times.t_j2 + 
                    (constraints.max_velocity - dir * start_conditions.v1) / constraints.max_acceleration;
    }

    times.t_v = start_conditions.h() / constraints.max_velocity
                - times.t_a / 2. * (1. + dir * start_conditions.v0 / constraints.max_velocity)
                - times.t_d / 2. * (1. + dir * start_conditions.v1 / constraints.max_velocity);
    
    if (times.t_v <= 0.) {
        return calc_times_case_2();
    }
    
    if (times.is_max_acceleration_not_reached()) {
        new_input.constraints.max_acceleration *= 0.5;
        if (new_input.constraints.max_acceleration > 0.01) {
            return new_input.calc_times_case_2();
        }
        new_input.constraints.max_acceleration = 0.;
    }
    
    handle_negative_acceleration_time(&times, &new_input);
    return times;
}

SCurveTimeIntervals SCurveInput::calc_times_case_2() const {
    SCurveTimeIntervals times;
    auto new_input = *this;
    int recursion_depth = 0;
    
    times = get_times_case_2();
    if (times.is_max_acceleration_not_reached()) {
        new_input.constraints.max_acceleration *= 0.5;
        if (new_input.constraints.max_acceleration > 0.01) {
            return new_input.calc_times_case_2();
        }
        new_input.constraints.max_acceleration = 0.;
    }
    handle_negative_acceleration_time(&times, &new_input);
    
    if (recursion_depth != 1) {
        new_input.constraints.max_acceleration *= 2.;
    }
    return calc_times_case_2_precise();
}

SCurveTimeIntervals SCurveInput::get_times_case_2() const {
    double t_j1 = constraints.max_acceleration / constraints.max_jerk;
    double t_j2 = constraints.max_acceleration / constraints.max_jerk;
    
    double delta = std::pow(constraints.max_acceleration, 4) / std::pow(constraints.max_jerk, 2)
                    + 2. * (std::pow(start_conditions.v0, 2) + std::pow(start_conditions.v1, 2))
                    + constraints.max_acceleration * (4. * start_conditions.h()
                    - 2. * constraints.max_acceleration / constraints.max_jerk 
                    * (start_conditions.v0 + start_conditions.v1));
    
    double t_a = (std::pow(constraints.max_acceleration, 2) / constraints.max_jerk
                - 2. * start_conditions.v0 + std::sqrt(delta))
                / (2. * constraints.max_acceleration);
    
    double t_d = (std::pow(constraints.max_acceleration, 2) / constraints.max_jerk
                - 2. * start_conditions.v1 + std::sqrt(delta))
                / (2. * constraints.max_acceleration);
    
    return SCurveTimeIntervals{t_j1, t_j2, t_a, 0., t_d};
}

SCurveTimeIntervals SCurveInput::calc_times_case_2_precise() const {
    SCurveTimeIntervals times = get_times_case_2();
    auto new_input = *this;
    
    if (times.is_max_acceleration_not_reached()) {
        new_input.constraints.max_acceleration *= 0.99;
        if (new_input.constraints.max_acceleration > 0.01) {
            return new_input.calc_times_case_2_precise();
        }
        new_input.constraints.max_acceleration = 0.;
    }
    handle_negative_acceleration_time(&times, &new_input);
    return times;
}

void SCurveInput::handle_negative_acceleration_time(SCurveTimeIntervals* times, const SCurveInput* new_input) const {
    if (times->t_a < 0.) {
        times->t_j1 = 0.;
        times->t_a = 0.;
        times->t_d = 2. * start_conditions.h() 
                    / (start_conditions.v0 + start_conditions.v1);
        times->t_j2 = (new_input->constraints.max_jerk * start_conditions.h()
                        - std::sqrt(
                            new_input->constraints.max_jerk 
                            * (new_input->constraints.max_jerk * std::pow(start_conditions.h(), 2)
                            + std::pow(start_conditions.v0 + start_conditions.v1, 2)
                            * (start_conditions.v1 - start_conditions.v0))
                        ))
                        / (new_input->constraints.max_jerk 
                        * (start_conditions.v1 + start_conditions.v0));
    }
    
    if (times->t_d < 0.) {
        times->t_j2 = 0.;
        times->t_d = 0.;
        times->t_a = 2. * start_conditions.h() 
                    / (start_conditions.v0 + start_conditions.v1);
        times->t_j2 = (new_input->constraints.max_jerk * start_conditions.h()
                        - std::sqrt(
                            new_input->constraints.max_jerk 
                            * (new_input->constraints.max_jerk * std::pow(start_conditions.h(), 2)
                            - std::pow(start_conditions.v0 + start_conditions.v1, 2)
                            * (start_conditions.v1 - start_conditions.v0))
                        ))
                        / (new_input->constraints.max_jerk 
                        * (start_conditions.v1 + start_conditions.v0));
    }
}

double eval_position(const SCurveParameters& p, double t) {
    const auto& times = p.time_intervals;
    if (t < 0.) {
        return p.conditions.q0;
    }
    double dir = p.conditions.dir();
    if (t <= times.t_j1) {
        return p.conditions.q0 + p.conditions.v0 * t + dir * p.j_max * std::pow(t, 3) / 6.;
    } else if (t <= times.t_a - times.t_j1) {
        return p.conditions.q0 + p.conditions.v0 * t + dir * p.j_max * std::pow(times.t_j1, 2) * (t - times.t_j1 / 3.) / 2.;
    } else if (t <= times.t_a) {
        return p.conditions.q0 + (p.conditions.v0 + dir * p.a_lim_a * (times.t_a - times.t_j1)) * t 
               - dir * p.a_lim_a * (times.t_a - t) * (times.t_a - t) / 2. 
               + dir * p.j_min * std::pow(times.t_a - t, 3) / 6.;
    } else if (t <= times.t_a + times.t_v) {
        return p.conditions.q1 - dir * (p.v_lim + dir * p.conditions.v1) * times.t_d / 2.
               + dir * p.v_lim * (t - times.t_a - times.t_v);
    } else if (t <= times.total_duration() - times.t_d + times.t_j2) {
        return p.conditions.q1 - dir * (p.v_lim + dir * p.conditions.v1) * times.t_d / 2.
               + dir * p.v_lim * (t - times.total_duration() + times.t_d)
               - dir * p.j_max * std::pow(t - times.total_duration() + times.t_d, 3) / 6.;
    } else if (t <= times.total_duration() - times.t_j2) {
        return p.conditions.q1 - dir * (p.v_lim + dir * p.conditions.v1) * times.t_d / 2.
               + dir * p.v_lim * (t - times.total_duration() + times.t_d)
               + dir * p.a_lim_d / 6. * (3. * std::pow(t - times.total_duration() + times.t_d, 2)
               - 3. * times.t_j2 * (t - times.total_duration() + times.t_d) + std::pow(times.t_j2, 2));
    } else if (t <= times.total_duration()) {
        return p.conditions.q1 - p.conditions.v1 * (times.total_duration() - t)
               - dir * p.j_max * std::pow(times.total_duration() - t, 3) / 6.;
    } else {
        return p.conditions.q1;
    }
}

double eval_velocity(const SCurveParameters& p, double t) {
    const auto& times = p.time_intervals;
    if (t < 0.) {
        return p.conditions.v0;
    }
    double dir = p.conditions.dir();
    if (t <= times.t_j1) {
        return p.conditions.v0 + dir * p.j_max * std::pow(t, 2) / 2.;
    } else if (t <= times.t_a - times.t_j1) {
        return p.conditions.v0 + dir * p.a_lim_a * (t - times.t_j1 / 2.);
    } else if (t <= times.t_a) {
        return dir * p.v_lim + dir * p.j_min * std::pow(times.t_a - t, 2) / 2.;
    } else if (t <= times.t_a + times.t_v) {
        return dir * p.v_lim;
    } else if (t <= times.total_duration() - times.t_d + times.t_j2) {
        return dir * p.v_lim - dir * p.j_max * std::pow(t - times.total_duration() + times.t_d, 2) / 2.;
    } else if (t <= times.total_duration() - times.t_j2) {
        return dir * p.v_lim + dir * p.a_lim_d * (t - times.total_duration() + times.t_d - times.t_j2 / 2.);
    } else if (t <= times.total_duration()) {
        return p.conditions.v1 + dir * p.j_max * std::pow(times.total_duration() - t, 2) / 2.;
    } else {
        return p.conditions.v1;
    }
}

double eval_acceleration(const SCurveParameters& p, double t) {
    const auto& times = p.time_intervals;
    double dir = p.conditions.dir();
    if (t < 0.) {
        return 0.;
    } else if (t <= times.t_j1) {
        return dir * p.j_max * t;
    } else if (t <= times.t_a - times.t_j1) {
        return dir * p.a_lim_a;
    } else if (t <= times.t_a) {
        return dir * (-p.j_min) * (times.t_a - t);
    } else if (t <= times.t_a + times.t_v) {
        return 0.;
    } else if (t <= times.total_duration() - times.t_d + times.t_j2) {
        return dir * (-p.j_max) * (t - times.total_duration() + times.t_d);
    } else if (t <= times.total_duration() - times.t_j2) {
        return dir * p.a_lim_d;
    } else if (t <= times.total_duration()) {
        return dir * (-p.j_max) * (times.total_duration() - t);
    } else {
        return 0.;
    }
}

double eval_jerk(const SCurveParameters& p, double t) {
    const auto& times = p.time_intervals;
    double dir = p.conditions.dir();
    if (t < 0.) {
        return 0.;
    } else if (t <= times.t_j1) {
        return dir * p.j_max;
    } else if (t <= times.t_a - times.t_j1) {
        return 0.;
    } else if (t <= times.t_a) {
        return dir * p.j_min;
    } else if (t <= times.t_a + times.t_v) {
        return 0.;
    } else if (t <= times.total_duration() - times.t_d + times.t_j2) {
        return dir * p.j_min;
    } else if (t <= times.total_duration() - times.t_j2) {
        return 0.;
    } else {
        return dir * p.j_max;
    }
}

std::pair<SCurveParameters, std::function<double(double)>> s_curve_generator(
    const SCurveInput& input_parameters,
    Derivative derivative) {
    
    auto times = input_parameters.calc_intervals();
    auto params = SCurveParameters::create(times, input_parameters);

    switch (derivative) {
        case Derivative::Position:
            return {params, [&params](double t) { return eval_position(params, t); }};
        case Derivative::Velocity:
            return {params, [&params](double t) { return eval_velocity(params, t); }};
        case Derivative::Acceleration:
            return {params, [&params](double t) { return eval_acceleration(params, t); }};
        case Derivative::Jerk:
            return {params, [&params](double t) { return eval_jerk(params, t); }};
    }
}