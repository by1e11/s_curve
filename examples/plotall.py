import matplotlib.pyplot as plt
import numpy as np
from s_curve import *

def main():
    # 设置S曲线参数
    constraints = SCurveConstraints(
        max_jerk=3.0,
        max_acceleration=2.0,
        max_velocity=3.0
    )
    
    start_conditions = SCurveStartConditions(
        q0=0.0,
        q1=10.0,
        v0=0.0,
        v1=0.0
    )
    
    input_params = SCurveInput(constraints, start_conditions)
    
    # 计算时间间隔
    times = input_params.calc_intervals()
    total_time = times.total_duration()
    
    # 生成时间序列
    t = np.linspace(0, total_time, 1000)
    
    # 计算各参数曲线
    params = SCurveParameters.new(times, input_params)
    
    # 计算位置、速度、加速度和加加速度
    position = [eval_position(params, ti) for ti in t]
    velocity = [eval_velocity(params, ti) for ti in t]
    acceleration = [eval_acceleration(params, ti) for ti in t]
    jerk = [eval_jerk(params, ti) for ti in t]
    
    # 绘制图形
    plt.figure(figsize=(10, 8))
    
    plt.subplot(4, 1, 1)
    plt.plot(t, position, label='Position')
    plt.ylabel('Position (m)')
    plt.legend()
    
    plt.subplot(4, 1, 2)
    plt.plot(t, velocity, label='Velocity')
    plt.ylabel('Velocity (m/s)')
    plt.legend()
    
    plt.subplot(4, 1, 3)
    plt.plot(t, acceleration, label='Acceleration')
    plt.ylabel('Acceleration (m/s²)')
    plt.legend()
    
    plt.subplot(4, 1, 4)
    plt.plot(t, jerk, label='Jerk')
    plt.ylabel('Jerk (m/s³)')
    plt.xlabel('Time (s)')
    plt.legend()
    
    plt.suptitle('S-Curve Velocity Motion Profile')
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()