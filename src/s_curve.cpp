#include "s_curve.h"
#include <cmath>
#include <functional>

double SCurveTimeIntervals::total_duration() const {
    return t_a + t_d + t_v;
}

bool SCurveTimeIntervals::is_max_acceleration_not_reached() const {
    return t_a < 2. * t_j1 || t_d < 2. * t_j2;
}

SCurveTimeIntervals SCurveInput::calc_intervals() const {
    return calc_times_case_1();
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

bool SCurveInput::is_a_max_not_reached() const {
    return (constraints.max_velocity - start_conditions.v0) * constraints.max_jerk < 
           std::pow(constraints.max_acceleration, 2);
}

bool SCurveInput::is_a_min_not_reached() const {
    return (constraints.max_velocity - start_conditions.v1) * constraints.max_jerk < 
           std::pow(constraints.max_acceleration, 2);
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
        return start_conditions.h() > 0.5 * (start_conditions.v1 + start_conditions.v0) * (
            t_j_star + std::abs(start_conditions.v1 - start_conditions.v0) / constraints.max_acceleration
        );
    }
    
    if (t_j_star < constraints.max_acceleration / constraints.max_jerk) {
        return start_conditions.h() > t_j_star * (start_conditions.v0 + start_conditions.v1);
    }
    
    return false;
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
    
    double t_j1 = constraints.max_acceleration / constraints.max_jerk;
    double t_j2 = constraints.max_acceleration / constraints.max_jerk;
    
    double delta = std::pow(constraints.max_acceleration, 4) / std::pow(constraints.max_jerk, 2) 
        + 2. * (std::pow(start_conditions.v0, 2) + std::pow(start_conditions.v1, 2))
        + constraints.max_acceleration * (4. * start_conditions.h() 
            - 2. * constraints.max_acceleration / constraints.max_jerk 
            * (start_conditions.v0 + start_conditions.v1));
    
    times.t_a = (std::pow(constraints.max_acceleration, 2) / constraints.max_jerk 
        - 2. * start_conditions.v0 + std::sqrt(delta)) 
        / (2. * constraints.max_acceleration);
    
    times.t_d = (std::pow(constraints.max_acceleration, 2) / constraints.max_jerk 
        - 2. * start_conditions.v1 + std::sqrt(delta)) 
        / (2. * constraints.max_acceleration);
    
    times.t_j1 = t_j1;
    times.t_j2 = t_j2;
    times.t_v = 0.;
    
    return times;
}

void SCurveInput::handle_negative_acceleration_time(
    SCurveTimeIntervals* times, 
    const SCurveInput* new_input) const {
    
    if (times->t_a < 0.) {
        times->t_j1 = 0.;
        times->t_a = 0.;
        times->t_d = 2. * start_conditions.h() / 
            (start_conditions.v0 + start_conditions.v1);
        times->t_j2 = (new_input->constraints.max_jerk * start_conditions.h() 
            - std::sqrt(new_input->constraints.max_jerk * 
                (new_input->constraints.max_jerk * std::pow(start_conditions.h(), 2) 
                + std::pow(start_conditions.v0 + start_conditions.v1, 2) 
                * (start_conditions.v1 - start_conditions.v0))))
            / (new_input->constraints.max_jerk * (start_conditions.v1 + start_conditions.v0));
    }
    
    if (times->t_d < 0.) {
        times->t_j2 = 0.;
        times->t_d = 0.;
        times->t_a = 2. * start_conditions.h() / 
            (start_conditions.v0 + start_conditions.v1);
        times->t_j1 = (new_input->constraints.max_jerk * start_conditions.h() 
            - std::sqrt(new_input->constraints.max_jerk * 
                (new_input->constraints.max_jerk * std::pow(start_conditions.h(), 2) 
                - std::pow(start_conditions.v0 + start_conditions.v1, 2) 
                * (start_conditions.v1 - start_conditions.v0))))
            / (new_input->constraints.max_jerk * (start_conditions.v1 + start_conditions.v0));
    }
}

static double eval_position(const SCurveParameters& params, double t) {
    const auto& times = params.time_intervals;
    double dir = params.conditions.dir();
    double q0 = params.conditions.q0;
    double v0 = params.conditions.v0;
    double j_max = params.j_max;
    double a_lim_a = params.a_lim_a;
    double v_lim = params.v_lim;
    double j_min = params.j_min;
    double a_lim_d = params.a_lim_d;
    double q1 = params.conditions.q1;

    if (t <= times.t_j1) {
        // jerk 上升阶段
        return q0 + dir * (v0 * t + j_max * t * t * t / 6.);
    } else if (t <= times.t_a - times.t_j1) {
        double t1 = t - times.t_j1;
        // 计算 jerk 上升阶段结束时的位置
        double q_at_tj1 = q0 + dir * (v0 * times.t_j1 + j_max * times.t_j1 * times.t_j1 * times.t_j1 / 6.);
        // 计算 jerk 上升阶段结束时的速度
        double v_at_tj1 = v0 + 0.5 * j_max * times.t_j1 * times.t_j1;
        return q_at_tj1 + dir * (v_at_tj1 * t1 + a_lim_a * t1 * t1 / 2.);
    } else if (t <= times.t_a) {
        double t1 = t - (times.t_a - times.t_j1);
        // 计算加速度恒定阶段结束时的位置
        double t_constant_end = times.t_a - times.t_j1;
        double t_constant = t_constant_end - times.t_j1;
        double q_at_constant_end = q0 + dir * (v0 * times.t_j1 + j_max * times.t_j1 * times.t_j1 * times.t_j1 / 6.)
                                  + dir * ((v0 + 0.5 * j_max * times.t_j1 * times.t_j1) * t_constant + a_lim_a * t_constant * t_constant / 2.);
        // 计算加速度恒定阶段结束时的速度
        double v_at_constant_end = v0 + 0.5 * j_max * times.t_j1 * times.t_j1 + a_lim_a * t_constant;
        return q_at_constant_end + dir * (v_at_constant_end * t1 + a_lim_a * t1 * t1 / 2. - j_max * t1 * t1 * t1 / 6.);
    } else if (t <= times.t_a + times.t_v) {
        double t1 = t - times.t_a;
        // 计算加速阶段结束时的位置
        double q_at_acc_end = q0 + dir * (v0 * times.t_j1 + j_max * times.t_j1 * times.t_j1 * times.t_j1 / 6.)
                             + dir * ((v0 + 0.5 * j_max * times.t_j1 * times.t_j1) * (times.t_a - 2 * times.t_j1) + a_lim_a * (times.t_a - 2 * times.t_j1) * (times.t_a - 2 * times.t_j1) / 2.)
                             + dir * ((v0 + 0.5 * j_max * times.t_j1 * times.t_j1 + a_lim_a * (times.t_a - 2 * times.t_j1)) * times.t_j1 + a_lim_a * times.t_j1 * times.t_j1 / 2. - j_max * times.t_j1 * times.t_j1 * times.t_j1 / 6.);
        return q_at_acc_end + dir * (v_lim * t1);
    } else if (t <= times.t_a + times.t_v + times.t_j2) {
        double t1 = t - (times.t_a + times.t_v);
        // 计算匀速阶段结束时的位置
        double q_at_vel_end = q0 + dir * (v0 * times.t_j1 + j_max * times.t_j1 * times.t_j1 * times.t_j1 / 6.)
                             + dir * ((v0 + 0.5 * j_max * times.t_j1 * times.t_j1) * (times.t_a - 2 * times.t_j1) + a_lim_a * (times.t_a - 2 * times.t_j1) * (times.t_a - 2 * times.t_j1) / 2.)
                             + dir * ((v0 + 0.5 * j_max * times.t_j1 * times.t_j1 + a_lim_a * (times.t_a - 2 * times.t_j1)) * times.t_j1 + a_lim_a * times.t_j1 * times.t_j1 / 2. - j_max * times.t_j1 * times.t_j1 * times.t_j1 / 6.)
                             + dir * (v_lim * times.t_v);
        // 计算匀速阶段结束时的速度
        double v_at_vel_end = v_lim;
        return q_at_vel_end + dir * (v_at_vel_end * t1 + j_min * t1 * t1 * t1 / 6.);
    } else if (t <= times.total_duration() - times.t_j2) {
        double t1 = t - (times.t_a + times.t_v + times.t_j2);
        // 计算减速 jerk 上升阶段结束时的位置
        double q_at_jerk_end = q0 + dir * (v0 * times.t_j1 + j_max * times.t_j1 * times.t_j1 * times.t_j1 / 6.)
                             + dir * ((v0 + 0.5 * j_max * times.t_j1 * times.t_j1) * (times.t_a - 2 * times.t_j1) + a_lim_a * (times.t_a - 2 * times.t_j1) * (times.t_a - 2 * times.t_j1) / 2.)
                             + dir * ((v0 + 0.5 * j_max * times.t_j1 * times.t_j1 + a_lim_a * (times.t_a - 2 * times.t_j1)) * times.t_j1 + a_lim_a * times.t_j1 * times.t_j1 / 2. - j_max * times.t_j1 * times.t_j1 * times.t_j1 / 6.)
                             + dir * (v_lim * times.t_v)
                             + dir * (v_lim * times.t_j2 + j_min * times.t_j2 * times.t_j2 * times.t_j2 / 6.);
        // 计算减速 jerk 上升阶段结束时的速度
        double v_at_jerk_end = v_lim + 0.5 * j_min * times.t_j2 * times.t_j2;
        return q_at_jerk_end + dir * (v_at_jerk_end * t1 + a_lim_d * t1 * t1 / 2.);
    } else {
        double t1 = t - (times.total_duration() - times.t_j2);
        // 计算减速加速度恒定阶段结束时的位置
        double t_constant_end = times.total_duration() - times.t_j2;
        double t_constant = t_constant_end - (times.t_a + times.t_v + times.t_j2);
        double q_at_constant_end = q0 + dir * (v0 * times.t_j1 + j_max * times.t_j1 * times.t_j1 * times.t_j1 / 6.)
                                  + dir * ((v0 + 0.5 * j_max * times.t_j1 * times.t_j1) * (times.t_a - 2 * times.t_j1) + a_lim_a * (times.t_a - 2 * times.t_j1) * (times.t_a - 2 * times.t_j1) / 2.)
                                  + dir * ((v0 + 0.5 * j_max * times.t_j1 * times.t_j1 + a_lim_a * (times.t_a - 2 * times.t_j1)) * times.t_j1 + a_lim_a * times.t_j1 * times.t_j1 / 2. - j_max * times.t_j1 * times.t_j1 * times.t_j1 / 6.)
                                  + dir * (v_lim * times.t_v)
                                  + dir * (v_lim * times.t_j2 + j_min * times.t_j2 * times.t_j2 * times.t_j2 / 6.)
                                  + dir * ((v_lim + 0.5 * j_min * times.t_j2 * times.t_j2) * t_constant + a_lim_d * t_constant * t_constant / 2.);
        // 计算减速加速度恒定阶段结束时的速度
        double v_at_constant_end = v_lim + 0.5 * j_min * times.t_j2 * times.t_j2 + a_lim_d * t_constant;
        return q_at_constant_end + dir * (v_at_constant_end * t1 + a_lim_d * t1 * t1 / 2. + (-j_min) * t1 * t1 * t1 / 6.);
    }
}

static double eval_velocity(const SCurveParameters& params, double t) {
    const auto& times = params.time_intervals;
    double dir = params.conditions.dir();
    double v0 = params.conditions.v0;
    double j_max = params.j_max;
    double a_lim_a = params.a_lim_a;
    double v_lim = params.v_lim;  // 匀速段速度
    double j_min = params.j_min;   // 减速阶段jerk值（通常为负）
    double a_lim_d = params.a_lim_d; // 减速阶段最大加速度（通常为负）
    double total_duration = times.total_duration();

    if (t <= times.t_j1) {
        // jerk 上升阶段：v = v0 + 0.5 * j_max * t²
        return dir * (v0 + 0.5 * j_max * t * t);
    } else if (t <= times.t_a - times.t_j1) {
        double t1 = t - times.t_j1;
        // 加速度恒定阶段：v = v(t_j1) + a_lim_a * t1
        double v_at_tj1 = v0 + 0.5 * j_max * times.t_j1 * times.t_j1;
        return dir * (v_at_tj1 + a_lim_a * t1);
    } else if (t <= times.t_a) {
        double t1 = t - (times.t_a - times.t_j1);
        // jerk 下降阶段：v = v(t_a - t_j1) + a_lim_a * t1 - 0.5 * j_max * t1²
        double t_constant_end = times.t_a - times.t_j1;
        double v_at_constant_end = v0 + 0.5 * j_max * times.t_j1 * times.t_j1 + a_lim_a * (t_constant_end - times.t_j1);
        return dir * (v_at_constant_end + a_lim_a * t1 - 0.5 * j_max * t1 * t1);
    } 
    else if (t <= times.t_a + times.t_v) {
        return dir * v_lim;  // 匀速段速度
    } 
    else if (t <= times.t_a + times.t_v + times.t_j2) {
        // 减速阶段jerk上升段（加速度从0线性增加到a_lim_d）
        double t1 = t - (times.t_a + times.t_v);
        // 初始速度为匀速段的v_lim
        return dir * (v_lim + 0.5 * j_min * t1 * t1);
    } 
    else if (t <= total_duration - times.t_j2) {
        // 减速阶段加速度恒定段（加速度保持a_lim_d）
        double t1 = t - (times.t_a + times.t_v + times.t_j2);
        // 初始速度为jerk上升段结束时的速度
        double v_at_jerk_end = v_lim + 0.5 * j_min * times.t_j2 * times.t_j2;
        return dir * (v_at_jerk_end + a_lim_d * t1);
    } 
    else {
        // 减速阶段jerk下降段（加速度从a_lim_d线性减小到0）
        double t1 = t - (times.total_duration() - times.t_j2);
        // 初始速度为恒定加速度段结束时的速度
        double t_constant_end = times.total_duration() - times.t_j2;
        double t_constant = t_constant_end - (times.t_a + times.t_v + times.t_j2);
        double v_at_constant_end = v_lim + 0.5 * params.j_min * times.t_j2 * times.t_j2 
                                  + params.a_lim_d * t_constant;
        // 速度 = v_at_constant_end + a_lim_d * t1 + 0.5 * jerk * t1²（jerk此时为正）
        return dir * (v_at_constant_end + params.a_lim_d * t1 + 0.5 * (-params.j_min) * t1 * t1);  // 修正jerk符号
    }
}

static double eval_acceleration(const SCurveParameters& params, double t) {
    const auto& times = params.time_intervals;
    double dir = params.conditions.dir();
    
    if (t <= times.t_j1) {
        return dir * params.j_max * t;
    } else if (t <= times.t_a - times.t_j1) {
        return dir * params.a_lim_a;
    } else if (t <= times.t_a) {
        double t1 = t - (times.t_a - times.t_j1);
        return dir * (params.a_lim_a - params.j_max * t1);
    } else if (t <= times.t_a + times.t_v) {
        return 0.;
    } else if (t <= times.t_a + times.t_v + times.t_j2) {
        double t1 = t - (times.t_a + times.t_v);
        return dir * params.j_min * t1;
    } else if (t <= times.total_duration() - times.t_j2) {
        return dir * params.a_lim_d;
    } else {
        double t1 = t - (times.total_duration() - times.t_j2);
        // 加速度 = a_lim_d + jerk * t1（jerk此时为正，故用+）
        return dir * (params.a_lim_d + (-params.j_min) * t1);  // 修正jerk符号
    }
}

static double eval_jerk(const SCurveParameters& params, double t) {
    const auto& times = params.time_intervals;
    
    if (t <= times.t_j1) {
        return params.j_max;
    } else if (t <= times.t_a - times.t_j1) {
        return 0.;
    } else if (t <= times.t_a) {
        return -params.j_max;
    } else if (t <= times.t_a + times.t_v) {
        return 0.;
    } else if (t <= times.t_a + times.t_v + times.t_j2) {
        return params.j_min;  // 减速初始段jerk（负）
    } else if (t <= times.total_duration() - times.t_j2) {
        return 0.;
    } else {
        return -params.j_min;  // 最后一段jerk（正，与初始段相反）
    }
}


std::pair<SCurveParameters, std::function<double(double)>> s_curve_generator(
    const SCurveInput& input_parameters,
    Derivative derivative) {
    
    auto times = input_parameters.calc_intervals();
    auto params = SCurveParameters::create(times, input_parameters);
    
    switch (derivative) {
        case Derivative::Position:
            return {params, [params](double t) { return eval_position(params, t); }};
        case Derivative::Velocity:
            return {params, [params](double t) { return eval_velocity(params, t); }};
        case Derivative::Acceleration:
            return {params, [params](double t) { return eval_acceleration(params, t); }};
        case Derivative::Jerk:
            return {params, [params](double t) { return eval_jerk(params, t); }};
        default:
            throw std::invalid_argument("Invalid derivative type");
    }
}

