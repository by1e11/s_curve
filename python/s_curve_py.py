import math
from enum import Enum
from dataclasses import dataclass
from typing import Callable, Tuple

class Derivative(Enum):
    POSITION = 0
    VELOCITY = 1
    ACCELERATION = 2
    JERK = 3

@dataclass
class SCurveConstraints:
    max_jerk: float
    max_acceleration: float 
    max_velocity: float

@dataclass
class SCurveTimeIntervals:
    t_j1: float
    t_j2: float
    t_a: float
    t_v: float
    t_d: float

    def total_duration(self) -> float:
        return self.t_a + self.t_d + self.t_v
    
    def is_max_acceleration_not_reached(self) -> bool:
        return self.t_a < 2. * self.t_j1 or self.t_d < 2. * self.t_j2

@dataclass
class SCurveStartConditions:
    q0: float
    q1: float
    v0: float
    v1: float

    def h(self) -> float:
        return abs(self.q1 - self.q0)
    
    def dir(self) -> float:
        return -1.0 if self.q1 < self.q0 else 1.0

@dataclass
class SCurveParameters:
    time_intervals: SCurveTimeIntervals
    j_max: float
    j_min: float
    a_lim_a: float
    a_lim_d: float
    v_lim: float
    conditions: SCurveStartConditions

    @staticmethod
    def new(times: SCurveTimeIntervals, p: 'SCurveInput') -> 'SCurveParameters':
        a_lim_a = p.constraints.max_jerk * times.t_j1
        a_lim_d = -p.constraints.max_jerk * times.t_j2
        v_lim = p.start_conditions.dir() * p.start_conditions.v0 + (times.t_a - times.t_j1) * a_lim_a
        
        return SCurveParameters(
            time_intervals=times,
            j_max=p.constraints.max_jerk,
            j_min=-p.constraints.max_jerk,
            a_lim_a=a_lim_a,
            a_lim_d=a_lim_d,
            v_lim=v_lim,
            conditions=p.start_conditions
        )

@dataclass
class SCurveInput:
    constraints: SCurveConstraints
    start_conditions: SCurveStartConditions

    def calc_intervals(self) -> SCurveTimeIntervals:
        return self.calc_times_case_1()
    
    def is_trajectory_feasible(self) -> bool:
        t_j_star = min(
            math.sqrt(
                abs(self.start_conditions.v1 - self.start_conditions.v0) / 
                self.constraints.max_jerk
            ),
            self.constraints.max_acceleration / self.constraints.max_jerk
        )
        
        if abs(t_j_star - self.constraints.max_acceleration / self.constraints.max_jerk) < 1e-10:
            return self.start_conditions.h() > 0.5 * (self.start_conditions.v1 + self.start_conditions.v0) * (
                t_j_star + abs(self.start_conditions.v1 - self.start_conditions.v0) / self.constraints.max_acceleration
            )
        
        if t_j_star < self.constraints.max_acceleration / self.constraints.max_jerk:
            return self.start_conditions.h() > t_j_star * (self.start_conditions.v0 + self.start_conditions.v1)
        
        return False
    
    def is_a_max_not_reached(self) -> bool:
        return (self.constraints.max_velocity - self.start_conditions.v0) * self.constraints.max_jerk < math.pow(
            self.constraints.max_acceleration, 2
        )
    
    def is_a_min_not_reached(self) -> bool:
        return (self.constraints.max_velocity - self.start_conditions.v1) * self.constraints.max_jerk < math.pow(
            self.constraints.max_acceleration, 2
        )
    
    def calc_times_case_1(self) -> SCurveTimeIntervals:
        times = SCurveTimeIntervals(0, 0, 0, 0, 0)
        new_input = self
        dir = self.start_conditions.dir()
        
        if self.is_a_max_not_reached():
            times.t_j1 = math.sqrt(
                (self.constraints.max_velocity - self.start_conditions.v0) / self.constraints.max_jerk
            )
            times.t_a = 2. * times.t_j1
        else:
            times.t_j1 = self.constraints.max_acceleration / self.constraints.max_jerk
            times.t_a = times.t_j1 + (
                self.constraints.max_velocity - dir * self.start_conditions.v0
            ) / self.constraints.max_acceleration

        if self.is_a_min_not_reached():
            times.t_j2 = math.sqrt(
                (self.constraints.max_velocity - self.start_conditions.v1) / self.constraints.max_jerk
            )
            times.t_d = 2. * times.t_j2
        else:
            times.t_j2 = self.constraints.max_acceleration / self.constraints.max_jerk
            times.t_d = times.t_j2 + (
                self.constraints.max_velocity - dir * self.start_conditions.v1
            ) / self.constraints.max_acceleration

        times.t_v = self.start_conditions.h() / self.constraints.max_velocity - times.t_a / 2. * (
            1. + dir * self.start_conditions.v0 / self.constraints.max_velocity
        ) - times.t_d / 2. * (
            1. + dir * self.start_conditions.v1 / self.constraints.max_velocity
        )
        
        if times.t_v <= 0.:
            return self.calc_times_case_2(0)
        
        if times.is_max_acceleration_not_reached():
            new_input.constraints.max_acceleration *= 0.5
            if new_input.constraints.max_acceleration > 0.01:
                return new_input.calc_times_case_2(0)
            new_input.constraints.max_acceleration = 0.
        
        self.handle_negative_acceleration_time(times, new_input)
        return times
    
    def calc_times_case_2(self, recursion_depth: int) -> SCurveTimeIntervals:
        recursion_depth += 1
        times = self.get_times_case_2()
        new_input = self
        
        if times.is_max_acceleration_not_reached():
            new_input.constraints.max_acceleration *= 0.5
            if new_input.constraints.max_acceleration > 0.01:
                return new_input.calc_times_case_2(recursion_depth)
            new_input.constraints.max_acceleration = 0.
        
        self.handle_negative_acceleration_time(times, new_input)
        
        if recursion_depth != 1:
            new_input.constraints.max_acceleration *= 2.
        
        return new_input.calc_times_case_2_precise(recursion_depth)
    
    def get_times_case_2(self) -> SCurveTimeIntervals:
        t_j1 = self.constraints.max_acceleration / self.constraints.max_jerk
        t_j2 = self.constraints.max_acceleration / self.constraints.max_jerk
        
        delta = math.pow(self.constraints.max_acceleration, 4) / math.pow(self.constraints.max_jerk, 2) + 2. * (
            math.pow(self.start_conditions.v0, 2) + math.pow(self.start_conditions.v1, 2)
        ) + self.constraints.max_acceleration * (
            4. * self.start_conditions.h() - 2. * self.constraints.max_acceleration / self.constraints.max_jerk * (
                self.start_conditions.v0 + self.start_conditions.v1
            )
        )
        
        t_a = (
            math.pow(self.constraints.max_acceleration, 2) / self.constraints.max_jerk - 2. * self.start_conditions.v0 + math.sqrt(delta)
        ) / (2. * self.constraints.max_acceleration)
        
        t_d = (
            math.pow(self.constraints.max_acceleration, 2) / self.constraints.max_jerk - 2. * self.start_conditions.v1 + math.sqrt(delta)
        ) / (2. * self.constraints.max_acceleration)
        
        return SCurveTimeIntervals(t_j1, t_j2, t_a, 0., t_d)
    
    def calc_times_case_2_precise(self, recursion_depth: int) -> SCurveTimeIntervals:
        recursion_depth += 1
        times = self.get_times_case_2()
        new_input = self
        
        if times.is_max_acceleration_not_reached():
            new_input.constraints.max_acceleration *= 0.99
            if new_input.constraints.max_acceleration > 0.01:
                return new_input.calc_times_case_2_precise(recursion_depth)
            new_input.constraints.max_acceleration = 0.
        
        self.handle_negative_acceleration_time(times, new_input)
        return times
    
    def handle_negative_acceleration_time(self, times: SCurveTimeIntervals, new_input: 'SCurveInput'):
        if times.t_a < 0.:
            times.t_j1 = 0.
            times.t_a = 0.
            times.t_d = 2. * self.start_conditions.h() / (
                self.start_conditions.v0 + self.start_conditions.v1
            )
            times.t_j2 = (
                new_input.constraints.max_jerk * self.start_conditions.h() - math.sqrt(
                    new_input.constraints.max_jerk * (
                        new_input.constraints.max_jerk * math.pow(self.start_conditions.h(), 2) + math.pow(
                            self.start_conditions.v0 + self.start_conditions.v1, 2
                        ) * (self.start_conditions.v1 - self.start_conditions.v0)
                    )
                )
            ) / (
                new_input.constraints.max_jerk * (self.start_conditions.v1 + self.start_conditions.v0)
            )
        
        if times.t_d < 0.:
            times.t_j2 = 0.
            times.t_d = 0.
            times.t_a = 2. * self.start_conditions.h() / (
                self.start_conditions.v0 + self.start_conditions.v1
            )
            times.t_j2 = (
                new_input.constraints.max_jerk * self.start_conditions.h() - math.sqrt(
                    new_input.constraints.max_jerk * (
                        new_input.constraints.max_jerk * math.pow(self.start_conditions.h(), 2) - math.pow(
                            self.start_conditions.v0 + self.start_conditions.v1, 2
                        ) * (self.start_conditions.v1 - self.start_conditions.v0)
                    )
                )
            ) / (
                new_input.constraints.max_jerk * (self.start_conditions.v1 + self.start_conditions.v0)
            )

def eval_position(p: SCurveParameters, t: float) -> float:
    times = p.time_intervals
    if t < 0.:
        return p.conditions.q0
    
    dir = p.conditions.dir()
    
    if t <= times.t_j1:
        return p.conditions.q0 + p.conditions.v0 * t + dir * p.j_max * math.pow(t, 3) / 6.
    elif t <= times.t_a - times.t_j1:
        return p.conditions.q0 + p.conditions.v0 * t + dir * p.a_lim_a / 6. * (
            3. * math.pow(t, 2) - 3. * times.t_j1 * t + math.pow(times.t_j1, 2)
        )
    elif t <= times.t_a:
        return p.conditions.q0 + dir * (p.v_lim + dir * p.conditions.v0) * times.t_a / 2. - dir * p.v_lim * (
            times.t_a - t
        ) - dir * p.j_min * math.pow(times.t_a - t, 3) / 6.
    elif t <= times.t_a + times.t_v:
        return p.conditions.q0 + dir * (p.v_lim + dir * p.conditions.v0) * times.t_a / 2. + dir * p.v_lim * (
            t - times.t_a
        )
    elif t <= times.total_duration() - times.t_d + times.t_j2:
        return p.conditions.q1 - dir * (p.v_lim + dir * p.conditions.v1) * times.t_d / 2. + dir * p.v_lim * (
            t - times.total_duration() + times.t_d
        ) - dir * p.j_max * math.pow(t - times.total_duration() + times.t_d, 3) / 6.
    elif t <= times.total_duration() - times.t_j2:
        return p.conditions.q1 - dir * (p.v_lim + dir * p.conditions.v1) * times.t_d / 2. + dir * p.v_lim * (
            t - times.total_duration() + times.t_d
        ) + dir * p.a_lim_d / 6. * (
            3. * math.pow(t - times.total_duration() + times.t_d, 2) - 3. * times.t_j2 * (
                t - times.total_duration() + times.t_d
            ) + math.pow(times.t_j2, 2)
        )
    elif t <= times.total_duration():
        return p.conditions.q1 - p.conditions.v1 * (times.total_duration() - t) - dir * p.j_max * math.pow(
            times.total_duration() - t, 3
        ) / 6.
    else:
        return p.conditions.q1

def eval_velocity(p: SCurveParameters, t: float) -> float:
    times = p.time_intervals
    if t < 0.:
        return p.conditions.v0
    
    dir = p.conditions.dir()
    
    if t <= times.t_j1:
        return p.conditions.v0 + dir * p.j_max * math.pow(t, 2) / 2.
    elif t <= times.t_a - times.t_j1:
        return p.conditions.v0 + dir * p.a_lim_a * (t - times.t_j1 / 2.)
    elif t <= times.t_a:
        return dir * p.v_lim + dir * p.j_min * math.pow(times.t_a - t, 2) / 2.
    elif t <= times.t_a + times.t_v:
        return dir * p.v_lim
    elif t <= times.total_duration() - times.t_d + times.t_j2:
        return dir * p.v_lim - dir * p.j_max * math.pow(t - times.total_duration() + times.t_d, 2) / 2.
    elif t <= times.total_duration() - times.t_j2:
        return dir * p.v_lim + dir * p.a_lim_d * (t - times.total_duration() + times.t_d - times.t_j2 / 2.)
    elif t <= times.total_duration():
        return p.conditions.v1 + dir * p.j_max * math.pow(times.total_duration() - t, 2) / 2.
    else:
        return p.conditions.v1

def eval_acceleration(p: SCurveParameters, t: float) -> float:
    times = p.time_intervals
    dir = p.conditions.dir()
    
    if t < 0.:
        return 0.
    elif t <= times.t_j1:
        return dir * p.j_max * t
    elif t <= times.t_a - times.t_j1:
        return dir * p.a_lim_a
    elif t <= times.t_a:
        return dir * (-p.j_min) * (times.t_a - t)
    elif t <= times.t_a + times.t_v:
        return 0.
    elif t <= times.total_duration() - times.t_d + times.t_j2:
        return dir * (-p.j_max) * (t - times.total_duration() + times.t_d)
    elif t <= times.total_duration() - times.t_j2:
        return dir * p.a_lim_d
    elif t <= times.total_duration():
        return dir * (-p.j_max) * (times.total_duration() - t)
    else:
        return 0.

def eval_jerk(p: SCurveParameters, t: float) -> float:
    times = p.time_intervals
    dir = p.conditions.dir()
    
    if t < times.t_j1:
        return dir * p.j_max
    elif t <= times.t_a - times.t_j1:
        return 0.
    elif t <= times.t_a:
        return dir * p.j_min
    elif t <= times.t_a + times.t_v:
        return 0.
    elif t <= times.total_duration() - times.t_d + times.t_j2:
        return dir * p.j_min
    elif t <= times.total_duration() - times.t_j2:
        return 0.
    else:
        return dir * p.j_max

def s_curve_generator(
    input_parameters: SCurveInput,
    derivative: Derivative
) -> Tuple[SCurveParameters, Callable[[float], float]]:
    times = input_parameters.calc_intervals()
    params = SCurveParameters.new(times, input_parameters)
    
    if derivative == Derivative.POSITION:
        return params, lambda t: eval_position(params, t)
    elif derivative == Derivative.VELOCITY:
        return params, lambda t: eval_velocity(params, t)
    elif derivative == Derivative.ACCELERATION:
        return params, lambda t: eval_acceleration(params, t)
    elif derivative == Derivative.JERK:
        return params, lambda t: eval_jerk(params, t)
    else:
        raise ValueError("Invalid derivative type")