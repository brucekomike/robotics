import matplotlib.pyplot as plt
import numpy as np
from scipy import signal

# --- 参数设置 ---

# 校正前的开环传递函数
# G(s) = k / (s * (s + 2))
k_pre_correction = 40  # 根据 Kv = 20 计算得出
numerator_pre = [k_pre_correction]
denominator_pre = [1, 2, 0]  # s^2 + 2s + 0 (这里是 s*(s+2) = s^2 + 2s)

# 创建校正前的传递函数系统
system_pre_correction = signal.TransferFunction(numerator_pre, denominator_pre)

# 频率范围
# 为了计算裕量，通常需要包含一些穿越频率的范围。
# 这里的频率范围可以根据实际情况调整，但需要足够宽以覆盖主导的频率响应。
# np.logspace(-1, 2, 500) 表示从 10^-1 (0.1) 到 10^2 (100) 共 500 个对数间隔的点。
omega = np.logspace(-1, 2, 500)

# --- 绘制 Bode 图 ---

# 创建图形和子图
plt.figure(figsize=(10, 8))

# 幅频特性
plt.subplot(2, 1, 1)
w_pre, mag_linear_pre, phase_rad_pre = signal.bode(system_pre_correction, omega)
mag_db_pre = 20 * np.log10(np.clip(mag_linear_pre, a_min=1e-10, a_max=None)) # 避免 log10(0)
plt.semilogx(w_pre, mag_db_pre, label=f'G(s) = {k_pre_correction}/(s(s+2)) (校正前)')
plt.title('Bode Plot - Magnitude Response (Before Correction)')
plt.xlabel('Frequency (rad/s)')
plt.ylabel('Magnitude (dB)')
plt.grid(True, which="both", ls="-")
plt.legend()

# 相频特性
plt.subplot(2, 1, 2)
phase_deg_pre = np.rad2deg(phase_rad_pre)
plt.semilogx(w_pre, phase_deg_pre, label=f'G(s) = {k_pre_correction}/(s(s+2)) (校正前)')
plt.title('Bode Plot - Phase Response (Before Correction)')
plt.xlabel('Frequency (rad/s)')
plt.ylabel('Phase (degrees)')
plt.grid(True, which="both", ls="-")
plt.legend()

plt.tight_layout() # 调整子图布局，防止重叠
plt.show()

# --- 计算并输出裕度 ---

# 1. 幅值裕量 (Magnitude Margin, Mg)
# 幅值裕量是增益为 0 dB 时的频率 (穿越频率) 处的相位。
# 或者，找到相位为 -180 度的频率，此时的幅值（dB）就是幅值裕量。
# 在 scipy.signal.bode 中，mag_db 和 phase_deg 是对应频率 w_pre 下的值。

# 找到相位从 -180度 附近的点（要注意相位可能不止一次经过 -180度，这里取第一次）
# 理想情况下，我们寻找相位等于 -180 度时的幅度，但实际可能找不到精确的 -180 度。
# 我们可以找到相位小于 -180 度的点，然后插值或取其近似。
# 更准确的方法是找到相位穿越 -180度 的频率，以及对应的幅度。

# 寻找相位为 -180 度时的频率点 (w_pm)
# 假设我们关心的是在某个频率点 mag_db = 0 时的相位，这对应相位裕量。
# 这里我们反过来，寻找 phase_deg = -180 时的幅度（dB）。
# 寻找相位达到 -180 度 (或接近 -180 度) 的索引
# 注意：对于本例 G(s) = 40 / (s(s+2))，其相位在低频时为 -90 度，在高频时趋向于 -180 度。
# 实际上，由于存在 s 因子，它在 s=0 时是 -90 度，在 s 趋于无穷时是 -180 度。
# 这里我们寻找相位接近 -180 度的频率点。
# 假设我们寻找相位为 -180 度的点（更准确说是接近-180度）
# 注意：对于纯粹的 G(s) = k/(s(s+2))，相位不会精确达到-180度，而是渐近。
# 我们可以查找相位小于-180度的第一个点，并考虑其附近的幅值。
# 或者，更常见的是，寻找增益穿越频率 (mag_db = 0) 时的相位，这才是相位裕量。
# 这里我们先计算相位裕量。

# 寻找增益穿越频率 (mag_db = 0 dB)
gain_crossover_indices = np.where(np.abs(mag_db_pre) < 0.5)[0] # 寻找幅度接近 0 dB 的点
if len(gain_crossover_indices) > 0:
    gain_crossover_freq = w_pre[gain_crossover_indices[0]]
    phase_at_gain_crossover = phase_deg_pre[gain_crossover_indices[0]]
    phase_margin = 180 + phase_at_gain_crossover # 相位裕量 Ph = 180 + arg(G(j*w_gc))
    print(f"增益穿越频率 (w_gc): {gain_crossover_freq:.2f} rad/s")
    print(f"相位裕量 (Phase Margin, PM): {phase_margin:.2f} degrees")
else:
    print("未找到增益穿越频率 (0 dB)。")
    phase_margin = None

# 2. 相位裕量 (Phase Margin, PM)
# 相位裕量是当开环增益 G(jω) 的幅值等于 1 (0 dB) 时的相位值，加上 180 度。
# 也就是说，找到 mag_db = 0 时的频率点 (增益穿越频率 w_gc)，然后查看该频率点的相位。

# 找到相位穿越 -180度 的频率点 (w_pm)
# 找到相位接近 -180 度的索引
# 注意：在这个系统中，相位会从 -90度 渐近到 -180度。
# 我们需要找到相位为 -180 度时的幅度（dB）。
# 寻找相位为 -180 度的频率点
phase_crossover_indices = np.where(np.abs(phase_deg_pre + 180) < 0.5)[0] # 寻找相位接近 -180 度的点

if len(phase_crossover_indices) > 0:
    phase_crossover_freq = w_pre[phase_crossover_indices[0]]
    mag_at_phase_crossover = mag_db_pre[phase_crossover_indices[0]]
    # 幅值裕量是在相位达到-180度时的幅度（dB）。
    # 注意：如果相位会多次穿越 -180度，需要考虑第一个或最相关的点。
    # 对于 G(s) = k/(s(s+2))，相位在 s=0 时是-90度，s->inf时是-180度。
    # 所以，相位只会渐近到-180度，不会精确等于-180度。
    # 通常，我们取相位为-180度时的幅度，或者在相位-180度附近的幅值。
    # 如果找不到精确-180度，可能需要插值。
    # 在此简单取接近 -180 度的幅度
    magnitude_margin_db = mag_at_phase_crossover
    print(f"相位穿越频率 (w_pc, 接近 -180 度): {phase_crossover_freq:.2f} rad/s")
    print(f"幅值裕量 (Magnitude Margin, Mg): {magnitude_margin_db:.2f} dB")
else:
    # 寻找相位在 -180 度附近的幅度
    # 如果精确的 -180 度点不存在，可以尝试寻找相位接近 -180 度的幅度
    # 例如，找到相位小于 -179 度的点
    indices_near_neg180 = np.where(phase_deg_pre < -179)[0]
    if len(indices_near_neg180) > 0:
        # 取第一个这样的点的幅度
        magnitude_margin_db = mag_db_pre[indices_near_neg180[0]]
        print(f"近似找到相位接近 -180 度的频率点，其幅度为: {magnitude_margin_db:.2f} dB (幅值裕量估计)")
    else:
        print("未找到相位穿越 -180 度或接近 -180 度的频率点。")
        magnitude_margin_db = None


# --- 输出结果 ---
print("\n--- 校正前系统性能指标 ---")
print(f"静态速度误差系数 Kv: {k_pre_correction/2:.2f} (1/s)") # Kv = k/2
if phase_margin is not None:
    print(f"相位裕量 PM: {phase_margin:.2f} degrees")
if magnitude_margin_db is not None:
    print(f"幅值裕量 Mg: {magnitude_margin_db:.2f} dB")