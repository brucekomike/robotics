import matplotlib.pyplot as plt
import numpy as np

# 定义参数 (示例值)
a_s_val = 0.2  # meters
b_s_val = 0.5  # meters
l_4_val = 0.1  # meters

# A_i0 坐标 (在 P-xyz 系)
a10 = np.array([0, -a_s_val, l_4_val])
a20 = np.array([a_s_val, 0, l_4_val])
a30 = np.array([-a_s_val, 0, l_4_val])
A_coords = np.array([a10, a20, a30])
# --- Matplotlib 中文显示设置 ---
plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'STHeiti']
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['axes.unicode_minus'] = False
# B_i 坐标 (在 B-xyz 系)
b1 = np.array([0, -b_s_val, 0])
b2 = np.array([b_s_val, 0, 0])
b3 = np.array([-b_s_val, 0, 0])
B_coords = np.array([b1, b2, b3])

# 绘制 B_i 点 (基座连接点)
fig1, ax1 = plt.subplots()
ax1.scatter(B_coords[:, 0], B_coords[:, 1], c='blue', marker='o', label='基座连接点 $B_i$')
for i, txt in enumerate(['$B_1$', '$B_2$', '$B_3$']):
    ax1.annotate(txt, (B_coords[i, 0], B_coords[i, 1]), textcoords="offset points", xytext=(0,5), ha='center')
ax1.set_xlabel('$x_B$ (m)')
ax1.set_ylabel('$y_B$ (m)')
ax1.set_title('基座连接点 $B_i$ 在 $B-xyz$ 平面投影 (俯视图)')
ax1.axis('equal')
ax1.grid(True)
ax1.legend()
plt.show() # 在实际执行时取消注释

# 绘制 A_i0 点 (动平台连接点，在 P-xyz 系的 x_p-y_p 平面投影)
fig2, ax2 = plt.subplots()
ax2.scatter(A_coords[:, 0], A_coords[:, 1], c='red', marker='s', label='动平台连接点 $A_i$ (在 $P$ 系)')
for i, txt in enumerate(['$A_1$', '$A_2$', '$A_3$']):
    ax2.annotate(txt, (A_coords[i, 0], A_coords[i, 1]), textcoords="offset points", xytext=(0,5), ha='center')
ax2.set_xlabel('$x_P$ (m)')
ax2.set_ylabel('$y_P$ (m)')
ax2.set_title('动平台连接点 $A_i$ 在 $P-x_p y_p z_p$ 系中 ($x_p-y_p$ 平面投影)')
ax2.text(0,0, f'$P$ (原点)\n$z_p$ 坐标均为 $l_4={l_4_val}$m', ha='center', va='center', color='gray', fontsize=9)
ax2.axis('equal')
ax2.grid(True)
ax2.legend()
plt.show() # 在实际执行时取消注释
