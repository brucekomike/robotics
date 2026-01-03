import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from credit import get_credit

try:
    plt.rcParams['font.family'] = 'Hiragino Sans GB'
    plt.rcParams['axes.unicode_minus'] = False
except Exception as e:
    print(f"Could not set PingFang SC directly: {e}")
    print("You may need to use matplotlib.font_manager to register the font.")

def plot_bode_and_nyquist(omega_n, damping_ratios):
    """
    绘制不同阻尼比下二阶系统的 Bode 图和 Nyquist 图。

    Args:
        omega_n (float): 无阻尼自然振荡频率。
        damping_ratios (list): 阻尼比的列表。
    """

    # 频率范围
    omega = np.logspace(-1, 2, 500) # 从 0.1 rad/s 到 100 rad/s

    # 绘制 Bode 图
    plt.figure(figsize=(12, 8))

    # 幅频特性
    plt.subplot(2, 1, 1)
    for xi in damping_ratios:
        numerator = [omega_n**2]
        denominator = [1, 2 * xi * omega_n, omega_n**2]
        system = signal.TransferFunction(numerator, denominator)
        w, mag_linear, phase_rad = signal.bode(system, omega)
        mag_db = 20 * np.log10(mag_linear)
        # phase_deg = np.rad2deg(phase_rad) # 在幅频图中不需要相位


        plt.semilogx(w, mag_db, label=f'ξ = {xi}')

    plt.title(f'Bode Plot - Magnitude Response for $\\omega_n$ = {omega_n} - {get_credit()}')
    plt.xlabel('Frequency (rad/s)')
    plt.ylabel('Magnitude (dB)')
    plt.grid(True, which="both", ls="-")
    plt.legend()

    # 相频特性
    plt.subplot(2, 1, 2)
    for xi in damping_ratios:
        numerator = [omega_n**2]
        denominator = [1, 2 * xi * omega_n, omega_n**2]
        system = signal.TransferFunction(numerator, denominator)
        w, mag_linear, phase_rad = signal.bode(system, omega)
        phase_deg = np.rad2deg(phase_rad)

        plt.semilogx(w, phase_deg, label=f'ξ = {xi}')

    plt.title(f'Bode Plot - Phase Response for $\\omega_n$ = {omega_n} - {get_credit()}')
    plt.xlabel('Frequency (rad/s)')
    plt.ylabel('Phase (degrees)')
    plt.grid(True, which="both", ls="-")
    plt.legend()

    plt.tight_layout()
    plt.show()

    # 绘制 Nyquist 图
    plt.figure(figsize=(8, 8))
    for xi in damping_ratios:
        numerator = [omega_n**2]
        denominator = [1, 2 * xi * omega_n, omega_n**2]
        system = signal.TransferFunction(numerator, denominator)

        freqs_nyquist, H_jw = signal.freqresp(system, omega)

        plt.plot(np.real(H_jw), np.imag(H_jw), label=f'ξ = {xi}')

    plt.title(f'Nyquist Plot for $\\omega_n$ = {omega_n} - {get_credit()}')
    plt.xlabel('Real Part')
    plt.ylabel('Imaginary Part')
    plt.grid(True)
    plt.axhline(0, color='black',linewidth=0.5)
    plt.axvline(0, color='black',linewidth=0.5)
    plt.legend()
    plt.axis('equal')
    plt.show()

# 设置参数
omega_n = 3
damping_ratios = [0.2, 0.4, 0.6, 0.8, 1.0]

# 绘制图形
plot_bode_and_nyquist(omega_n, damping_ratios)