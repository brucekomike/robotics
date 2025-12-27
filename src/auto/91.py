import numpy as np
import matplotlib.pyplot as plt
import control # 导入 control 库

# 定义系统参数
omega_n = 6
xi_values = [0.2, 0.4, 0.6, 0.8, 1.0, 2.0]
time_end_step = 5 # 阶跃响应仿真结束时间
time_end_impulse = 3 # 脉冲响应仿真结束时间

# --- 单位阶跃响应曲线 ---
plt.figure(figsize=(12, 8))

for xi in xi_values:
    # 构建传递函数 G(s) = omega_n^2 / (s^2 + 2*xi*omega_n*s + omega_n^2)
    numerator = [omega_n**2]
    denominator = [1, 2 * xi * omega_n, omega_n**2]
    system = control.tf(numerator, denominator) # 使用 control.tf 创建传递函数

    # 绘制单位阶跃响应
    # control.step_response 返回时间和输出
    t_step, y_step = control.step_response(system, T=np.linspace(0, time_end_step, 500))
    plt.plot(t_step, y_step, label=f'$\\xi = {xi}$')

plt.title('Unit Step Response for Different Damping Ratios ($\omega_n = 6$)')
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
plt.grid(True)
plt.legend()
plt.xlim(0, time_end_step)
plt.ylim(0, 1.6) # 调整y轴范围，以便更好地观察不同曲线
plt.show()

# --- 单位脉冲响应曲线 ---
plt.figure(figsize=(12, 8))

for xi in xi_values:
    # 构建传递函数 G(s) = omega_n^2 / (s^2 + 2*xi*omega_n*s + omega_n^2)
    numerator = [omega_n**2]
    denominator = [1, 2 * xi * omega_n, omega_n**2]
    system = control.tf(numerator, denominator) # 使用 control.tf 创建传递函数

    # 绘制单位脉冲响应
    # control.impulse_response 返回时间和输出
    t_impulse, y_impulse = control.impulse_response(system, T=np.linspace(0, time_end_impulse, 500))
    plt.plot(t_impulse, y_impulse, label=f'$\\xi = {xi}$')

plt.title('Unit Impulse Response for Different Damping Ratios ($\omega_n = 6$)')
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
plt.grid(True)
plt.legend()
plt.xlim(0, time_end_impulse)
plt.ylim(-3, 4.7) # 调整y轴范围
plt.show()