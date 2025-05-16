#ifndef S_CURVE_H
#define S_CURVE_H

#include <functional>
#include <utility>

class SCurveConstraints {
public:
    double max_jerk;
    double max_acceleration;
    double max_velocity;
};

enum class Derivative {
    Position = 0,
    Velocity = 1,
    Acceleration = 2,
    Jerk = 3
};

class SCurveTimeIntervals {
public:
    double t_j1;
    double t_j2;
    double t_a;
    double t_v;
    double t_d;

    double total_duration() const;
    bool is_max_acceleration_not_reached() const;
};

class SCurveStartConditions {
public:
    double q0;
    double q1;
    double v0;
    double v1;

    double h() const;
    double dir() const;
};

class SCurveInput {
public:
    SCurveConstraints constraints;
    SCurveStartConditions start_conditions;

    SCurveTimeIntervals calc_intervals() const;
    bool is_trajectory_feasible() const;

private:
    bool is_a_max_not_reached() const;
    bool is_a_min_not_reached() const;
    SCurveTimeIntervals calc_times_case_1() const;
    SCurveTimeIntervals calc_times_case_2() const;
    SCurveTimeIntervals get_times_case_2() const;
    SCurveTimeIntervals calc_times_case_2_precise() const;
    void handle_negative_acceleration_time(SCurveTimeIntervals* times, const SCurveInput* new_input) const;
};

class SCurveParameters {
public:
    SCurveTimeIntervals time_intervals;
    double j_max;
    double j_min;
    double a_lim_a;
    double a_lim_d;
    double v_lim;
    SCurveStartConditions conditions;

    static SCurveParameters create(const SCurveTimeIntervals& times, const SCurveInput& p);
};

std::pair<SCurveParameters, std::function<double(double)>> s_curve_generator(
    const SCurveInput& input_parameters,
    Derivative derivative);

#endif // S_CURVE_H