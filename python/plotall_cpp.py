import matplotlib.pyplot as plt
import numpy as np
import s_curve

def plot_s_curve():
    constraints = s_curve.SCurveConstraints()
    constraints.max_jerk = 3.0
    constraints.max_acceleration = 2.0
    constraints.max_velocity = 3.0

    start_conditions = s_curve.SCurveStartConditions()
    start_conditions.q0 = 0.0
    start_conditions.q1 = 10.0
    start_conditions.v0 = 0.0
    start_conditions.v1 = 0.0

    input_params = s_curve.SCurveInput()
    input_params.constraints = constraints
    input_params.start_conditions = start_conditions

    # 计算时间间隔
    times = input_params.calc_intervals()
    print(f"计算结果: t_j1={times.t_j1}, t_a={times.t_a}, t_v={times.t_v}, t_d={times.t_d}")
    if times.t_v <= 0:
        print("警告: 无法达到最大速度阶段，请增加距离或减小加速度")
    total_time = times.total_duration()
    
    # 生成时间序列
    t = np.linspace(0, total_time, 1000)
    
    # 获取所有导数曲线
    params_pos, func_pos = s_curve.s_curve_generator(input_params, s_curve.Derivative.POSITION)
    params_vel, func_vel = s_curve.s_curve_generator(input_params, s_curve.Derivative.VELOCITY)
    params_acc, func_acc = s_curve.s_curve_generator(input_params, s_curve.Derivative.ACCELERATION)
    params_jerk, func_jerk = s_curve.s_curve_generator(input_params, s_curve.Derivative.JERK)

    # 计算各参数曲线
    position = [func_pos(ti) for ti in t]
    velocity = [func_vel(ti) for ti in t]
    acceleration = [func_acc(ti) for ti in t]
    jerk = [func_jerk(ti) for ti in t]
    
    # 创建图形
    plt.figure(figsize=(12, 8))
    
    # 位置曲线
    plt.subplot(4, 1, 1)
    plt.plot(t, position, 'b-', linewidth=2)
    plt.ylabel('Position (m)')
    plt.grid(True)
    
    # 速度曲线
    plt.subplot(4, 1, 2)
    plt.plot(t, velocity, 'r-', linewidth=2)
    plt.ylabel('Velocity (m/s)')
    plt.grid(True)
    
    # 加速度曲线
    plt.subplot(4, 1, 3)
    plt.plot(t, acceleration, 'g-', linewidth=2)
    plt.ylabel('Acceleration (m/s²)')
    plt.grid(True)
    
    # 加加速度曲线
    plt.subplot(4, 1, 4)
    plt.plot(t, jerk, 'm-', linewidth=2)
    plt.ylabel('Jerk (m/s³)')
    plt.xlabel('Time (s)')
    plt.grid(True)
    
    plt.suptitle('S-Curve Motion Profile', fontsize=14)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    plot_s_curve()