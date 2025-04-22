# %% 导入库 (Import Libraries)
import numpy as np
import sympy as sp
import math
import matplotlib.pyplot as plt
# mpl_toolkits 在新版本 matplotlib 中不再需要显式导入 Axes3D
# from mpl_toolkits.mplot3d import Axes3D # 用于 3D 绘图 (For 3D plotting)
import time # 用于动画暂停 (For animation pause)

# 使用 SymPy 的美观打印 (Use SymPy's pretty printing)
sp.init_printing(use_unicode=True)

# %% 机器人参数 (Robot Parameters)
def get_robot_params():
    """返回包含机器人尺度参数的字典"""
    return {
        'a11': 135.0,   # 移动平台上连接点到平台中心的距离
        'a22': 135.0,   # 同上
        'a33': 135.0,   # 同上
        'b11': 570.0,   # 固定平台上连接点到平台中心的距离
        'b22': 320.0,   # 同上
        'b33': 320.0,   # 同上
        'e': 345.0      # 中心柱固定长度
    }

# %% 位置逆解函数 (Inverse Kinematics Function)
def inverse_kinematics(xp, yp, zp, params):
    """
    计算机器人逆运动学。
    输入:
        xp, yp, zp: 末端执行器目标位置 (float)
        params: 包含机器人尺寸的字典 (dict)
    输出:
        一个包含逆解结果的字典，如果无法求解则返回 None:
        {
            'q1', 'q2', 'q3': 驱动杆计算长度 (float) - 注意其定义可能不标准
            'q4': 中心杆可变长度 (float)
            'theta1', 'theta2': 平台姿态角 (rad) (float)
            'rp': 末端位置向量 (numpy.ndarray, shape (3,1))
            'R4': 移动平台旋转矩阵 (numpy.ndarray, shape (3,3))
            'w4': 中心杆方向向量 (numpy.ndarray, shape (3,1))
            'error': 错误信息 (str or None)
        }
    """
    a11 = params['a11']
    a22 = params['a22']
    a33 = params['a33']
    b11 = params['b11']
    b22 = params['b22']
    b33 = params['b33']
    e = params['e']

    rp = np.array([[xp], [yp], [zp]]) # 列向量
    norm_rp = np.linalg.norm(rp)

    if norm_rp < e - 1e-6: # 添加容差以避免浮点问题
        return {'error': f"警告：末端位置 (范数 {norm_rp:.2f}) 在中心柱固定长度 (e={e}) 范围内，无法达到。"}
        # 返回 None 或包含错误信息的字典

    q4 = norm_rp - e
    if q4 < 0: q4 = 0 # q4 代表可变长度，不能为负

    # 避免除以零
    denominator = q4 + e
    if abs(denominator) < 1e-9:
        return {'error': "错误：q4 + e 接近零，无法计算角度。"}

    # 计算 theta2
    sin_theta2 = xp / denominator
    sin_theta2 = np.clip(sin_theta2, -1.0, 1.0) # 限制范围
    theta2 = np.arcsin(sin_theta2)

    # 计算 theta1
    cos_theta2 = np.cos(theta2)
    if abs(cos_theta2) < 1e-9:
        # 当 theta2 接近 +/- pi/2 (平台几乎垂直于 Z 轴)
        # 此时偏航角定义可能不稳定，通常发生在奇异点附近
        # 如果 yp 和 zp 同时接近 0，则 theta1 未定义，可设为 0
        if abs(yp) < 1e-9 and abs(zp) < 1e-9:
            theta1 = 0.0
            # print("警告：cos(theta2) 接近零，且 yp, zp 接近零，theta1 设为 0。")
        else:
            # 使用 atan2 更稳健地处理四象限
            theta1 = np.arctan2(-yp, zp)
            # print("警告：cos(theta2) 接近零，使用 atan2(-yp, zp) 计算 theta1。")
    else:
        arg1 = (-yp / denominator) / cos_theta2
        arg2 = (zp / denominator) / cos_theta2
        theta1 = np.arctan2(arg1, arg2)

    # 计算旋转矩阵 R4
    c1, s1 = np.cos(theta1), np.sin(theta1)
    c2, s2 = np.cos(theta2), np.sin(theta2)
    R4 = np.array([
        [c2,   0,  s2],
        [s1*s2, c1, -s1*c2],
        [-c1*s2, s1, c1*c2]
    ])

    # 固定平台连接点 (在固定坐标系)
    b1 = np.array([[0], [-b11], [0]])
    b2 = np.array([[b22], [0], [0]])
    b3 = np.array([[-b33], [0], [0]])

    # 移动平台连接点 (在移动平台自身坐标系)
    a10 = np.array([[0], [-a11], [0]])
    a20 = np.array([[a22], [0], [0]])
    a30 = np.array([[-a33], [0], [0]])

    # 移动平台连接点 (在固定坐标系中的姿态，相对平台中心 rp)
    a1 = R4 @ a10
    a2 = R4 @ a20
    a3 = R4 @ a30

    # 中心杆方向向量 w4 (固定坐标系)
    w4 = np.array([[s2], [-s1*c2], [c1*c2]])
    # 验证 w4 是否等于 rp / norm(rp)
    # if norm_rp > 1e-6:
    #     w4_check = rp / norm_rp
    #     if not np.allclose(w4, w4_check):
    #          print(f"警告: w4 计算不一致! w4={w4.flatten()}, rp/norm(rp)={w4_check.flatten()}")

    # 计算驱动杆长度 q1, q2, q3 (根据原始 MATLAB 公式)
    # 再次注意：这个公式可能不是标准的物理连杆长度
    # 标准长度 L_i = || (rp + a_i) - b_i ||
    q1 = np.linalg.norm(rp - b1 + a1 - e*w4)
    q2 = np.linalg.norm(rp - b2 + a2 - e*w4)
    q3 = np.linalg.norm(rp - b3 + a3 - e*w4)

    return {
        'q1': q1, 'q2': q2, 'q3': q3,
        'q4': q4, 'theta1': theta1, 'theta2': theta2,
        'rp': rp, 'R4': R4, 'w4': w4,
        'error': None
    }

# %% 位置正解 - 符号设置与数值函数生成 (Forward Kinematics - Symbolic Setup)
def define_symbolic_fk(params, q1_in, q2_in, q3_in):
    """
    定义正解的符号方程、雅可比矩阵，并返回数值计算函数。
    输入:
        params: 机器人参数字典
        q1_in, q2_in, q3_in: 输入的驱动杆计算长度
    输出:
        (Func_sym, J_sym, num_func, num_jac, X_sym)
        Func_sym: 符号方程组 (SymPy Matrix)
        J_sym: 符号雅可比矩阵 (SymPy Matrix)
        num_func: Func_sym 的数值计算函数 (callable)
        num_jac: J_sym 的数值计算函数 (callable)
        X_sym: 符号变量向量 [q4, theta1, theta2] (SymPy Matrix)
    """
    a11 = params['a11']
    a22 = params['a22']
    a33 = params['a33']
    b11 = params['b11']
    b22 = params['b22']
    b33 = params['b33']
    # e = params['e'] # e 不直接出现在 F, G, H 方程中

    # 定义符号变量
    q4_sym, theta1_sym, theta2_sym = sp.symbols('q4_sym theta1_sym theta2_sym')
    X_sym = sp.Matrix([q4_sym, theta1_sym, theta2_sym])

    # 定义约束方程 F, G, H (来自连杆闭环)
    # 这些方程将输入的 q1_in, q2_in, q3_in 与未知的 q4, theta1, theta2 联系起来
    F_sym = a11**2 + b11**2 + q4_sym**2 - q1_in**2 \
            - 2*a11*b11*sp.cos(theta1_sym) \
            - 2*b11*q4_sym*sp.sin(theta1_sym)*sp.cos(theta2_sym)

    G_sym = a22**2 + b22**2 + q4_sym**2 - q2_in**2 \
            - 2*b22*q4_sym*sp.sin(theta2_sym) \
            - 2*a22*b22*sp.cos(theta2_sym)

    H_sym = a33**2 + b33**2 + q4_sym**2 - q3_in**2 \
            + 2*b33*q4_sym*sp.sin(theta2_sym) \
            - 2*a33*b33*sp.cos(theta2_sym)

    Func_sym = sp.Matrix([F_sym, G_sym, H_sym])

    # 计算雅可比矩阵 J = d(Func)/d(X)
    J_sym = Func_sym.jacobian(X_sym)

    # 转换为快速数值函数
    num_func = sp.lambdify((q4_sym, theta1_sym, theta2_sym), Func_sym, 'numpy')
    num_jac = sp.lambdify((q4_sym, theta1_sym, theta2_sym), J_sym, 'numpy')

    return Func_sym, J_sym, num_func, num_jac, X_sym

# %% 位置正解 - 牛顿迭代求解器 (Forward Kinematics - Newton Solver)
def forward_kinematics(q1_in, q2_in, q3_in, params, initial_guess, max_iter=100, tolerance=1e-6):
    """
    使用牛顿法迭代求解正运动学。
    输入:
        q1_in, q2_in, q3_in: 驱动杆计算长度 (float)
        params: 机器人参数字典 (dict)
        initial_guess: [q4_guess, theta1_guess, theta2_guess] (list or tuple)
        max_iter: 最大迭代次数 (int)
        tolerance: 收敛阈值 (float)
    输出:
        一个包含正解结果的字典，如果失败则包含错误信息:
        {
            'xp_fk', 'yp_fk', 'zp_fk': 计算得到的末端位置 (float)
            'q4_k', 'theta1_k', 'theta2_k': 计算得到的姿态变量 (float)
            'iterations': 实际迭代次数 (int)
            'final_error_norm': 最终方程误差范数 (float)
            'converged': 是否收敛 (bool)
            'error': 错误信息 (str or None)
        }
    """
    e = params['e']

    # 获取符号定义和数值函数
    Func_sym, J_sym, num_func, num_jac, X_sym = define_symbolic_fk(params, q1_in, q2_in, q3_in)

    # 初始猜测值
    q4_k, theta1_k, theta2_k = initial_guess
    X_k = np.array(initial_guess, dtype=float) # 使用 numpy 数组

    iterations = 0
    converged = False
    final_error_norm = float('inf')

    # print("\n--- 开始牛顿迭代求解正解 ---")
    # print(f"初始猜测: q4={q4_k:.4f}, theta1={theta1_k:.4f}, theta2={theta2_k:.4f}")
    # print(f"目标 q1={q1_in:.4f}, q2={q2_in:.4f}, q3={q3_in:.4f}")

    for i in range(max_iter):
        iterations = i
        # 计算当前点的函数值 F(X_k) (需要是列向量)
        F_k_mat = num_func(X_k[0], X_k[1], X_k[2])
        F_k = F_k_mat.flatten() # 转换为 1D 数组

        # 检查收敛性
        current_error_norm = np.linalg.norm(F_k)
        # print(f"Iter {i}: q4={X_k[0]:.6f}, t1={X_k[1]:.6f}, t2={X_k[2]:.6f}, ErrorNorm={current_error_norm:.6e}")
        if current_error_norm < tolerance:
            converged = True
            final_error_norm = current_error_norm
            # print(f"\n在 {i} 次迭代后收敛。")
            break

        # 计算当前点的雅可比矩阵 J(X_k)
        J_k = num_jac(X_k[0], X_k[1], X_k[2])

        # 求解线性方程组 J(X_k) * delta_X = -F(X_k)
        try:
            delta_X = np.linalg.solve(J_k, -F_k)
        except np.linalg.LinAlgError:
            # print(f"\n迭代 {i} 时雅可比矩阵奇异，无法求解。检查初始猜测值或模型。")
            return {
                'xp_fk': None, 'yp_fk': None, 'zp_fk': None,
                'q4_k': X_k[0], 'theta1_k': X_k[1], 'theta2_k': X_k[2],
                'iterations': i, 'final_error_norm': current_error_norm,
                'converged': False, 'error': "雅可比矩阵奇异"
            }

        # 更新变量 X_{k+1} = X_k + delta_X
        X_k += delta_X

    else: # while 循环正常结束 (未 break)，意味着达到最大迭代次数
        final_error_norm = np.linalg.norm(num_func(X_k[0], X_k[1], X_k[2]).flatten())
        # print(f"\n达到最大迭代次数 {max_iter}，可能未收敛。")
        # print(f"当前误差范数: {final_error_norm:.6e}")
        return {
            'xp_fk': None, 'yp_fk': None, 'zp_fk': None,
            'q4_k': X_k[0], 'theta1_k': X_k[1], 'theta2_k': X_k[2],
            'iterations': max_iter, 'final_error_norm': final_error_norm,
            'converged': False, 'error': "达到最大迭代次数"
        }

    # 迭代成功，计算最终的末端执行器位置
    q4_final, theta1_final, theta2_final = X_k
    xp_fk = (q4_final + e) * np.sin(theta2_final)
    yp_fk = (q4_final + e) * (-np.sin(theta1_final) * np.cos(theta2_final))
    zp_fk = (q4_final + e) * np.cos(theta1_final) * np.cos(theta2_final)

    return {
        'xp_fk': xp_fk, 'yp_fk': yp_fk, 'zp_fk': zp_fk,
        'q4_k': q4_final, 'theta1_k': theta1_final, 'theta2_k': theta2_final,
        'iterations': iterations, 'final_error_norm': final_error_norm,
        'converged': True, 'error': None
    }

# %% 可视化函数 (Visualization Function)
def plot_robot(ax, rp, R4, params, title="Robot Configuration"):
    """
    在给定的 3D 坐标轴上绘制机器人结构。
    输入:
        ax: Matplotlib 3D 坐标轴对象 (Axes3D)
        rp: 末端位置向量 (numpy.ndarray, shape (3,1))
        R4: 移动平台旋转矩阵 (numpy.ndarray, shape (3,3))
        params: 机器人参数字典 (dict)
        title: 图表标题 (str)
    """
    a11 = params['a11']
    a22 = params['a22']
    a33 = params['a33']
    b11 = params['b11']
    b22 = params['b22']
    b33 = params['b33']
    e = params['e']

    # 清除之前的绘图
    ax.cla()

    # --- 计算关键点坐标 ---
    # 原点 O
    O = np.array([0, 0, 0])
    # 固定平台连接点 B1, B2, B3 (在固定坐标系)
    B1 = np.array([0, -b11, 0])
    B2 = np.array([b22, 0, 0])
    B3 = np.array([-b33, 0, 0])
    # 移动平台连接点 A10, A20, A30 (在移动平台自身坐标系)
    A10 = np.array([0, -a11, 0])
    A20 = np.array([a22, 0, 0])
    A30 = np.array([-a33, 0, 0])
    # 移动平台连接点 A1, A2, A3 (相对平台中心 rp 的向量，在固定坐标系)
    A1_vec = R4 @ A10.reshape(3, 1)
    A2_vec = R4 @ A20.reshape(3, 1)
    A3_vec = R4 @ A30.reshape(3, 1)
    # 移动平台连接点 A1', A2', A3' 的绝对坐标 (在固定坐标系)
    P = rp.flatten() # 末端平台中心点
    A1_prime = P + A1_vec.flatten()
    A2_prime = P + A2_vec.flatten()
    A3_prime = P + A3_vec.flatten()
    # 中心柱顶端 E (固定坐标系)
    # 需要 w4，可以从 IK 结果获取，或从 FK 结果 (theta1, theta2) 重新计算
    s2 = np.sin(np.arcsin(np.clip(P[0] / np.linalg.norm(P), -1, 1))) # 近似 theta2
    # 更准确地，如果知道 theta1, theta2
    ik_result_temp = inverse_kinematics(P[0], P[1], P[2], params) # 临时计算获取角度
    if ik_result_temp and ik_result_temp['error'] is None:
        theta1_vis = ik_result_temp['theta1']
        theta2_vis = ik_result_temp['theta2']
        c1_vis, s1_vis = np.cos(theta1_vis), np.sin(theta1_vis)
        c2_vis, s2_vis = np.cos(theta2_vis), np.sin(theta2_vis)
        w4_vis = np.array([s2_vis, -s1_vis*c2_vis, c1_vis*c2_vis])
        E = e * w4_vis
    else: # 如果逆解失败或 P 在原点，用 Z 轴近似
         w4_vis = np.array([0, 0, 1]) if np.linalg.norm(P) < 1e-6 else P / np.linalg.norm(P)
         E = e * w4_vis


    # --- 绘制结构 ---
    # 固定平台 (灰色虚线)
    fixed_platform_pts = np.array([B1, B2, B3, B1])
    ax.plot(fixed_platform_pts[:, 0], fixed_platform_pts[:, 1], fixed_platform_pts[:, 2], 'k--', label='Fixed Platform', color='gray')
    ax.scatter(fixed_platform_pts[:-1, 0], fixed_platform_pts[:-1, 1], fixed_platform_pts[:-1, 2], c='black', marker='s', s=50) # 连接点 B_i

    # 移动平台 (蓝色实线)
    moving_platform_pts = np.array([A1_prime, A2_prime, A3_prime, A1_prime])
    ax.plot(moving_platform_pts[:, 0], moving_platform_pts[:, 1], moving_platform_pts[:, 2], 'b-', label='Moving Platform')
    ax.scatter(moving_platform_pts[:-1, 0], moving_platform_pts[:-1, 1], moving_platform_pts[:-1, 2], c='blue', marker='o', s=50) # 连接点 A'_i

    # 驱动杆 (红色实线) - 绘制物理连杆 B_i 到 A'_i
    ax.plot([B1[0], A1_prime[0]], [B1[1], A1_prime[1]], [B1[2], A1_prime[2]], 'r-', label='Driving Rods')
    ax.plot([B2[0], A2_prime[0]], [B2[1], A2_prime[1]], [B2[2], A2_prime[2]], 'r-')
    ax.plot([B3[0], A3_prime[0]], [B3[1], A3_prime[1]], [B3[2], A3_prime[2]], 'r-')

    # 中心柱 (绿色实线)
    ax.plot([O[0], E[0]], [O[1], E[1]], [O[2], E[2]], 'g-', linewidth=3, label='Central Column (Fixed Part)')
    ax.plot([E[0], P[0]], [E[1], P[1]], [E[2], P[2]], 'g--', linewidth=2, label='Central Column (Variable Part)') # 可变部分

    # 标出末端点 P
    ax.scatter(P[0], P[1], P[2], c='magenta', marker='*', s=100, label='End Effector (P)')
    # 标出原点 O
    ax.scatter(O[0], O[1], O[2], c='black', marker='x', s=50, label='Origin (O)')

    # --- 设置绘图属性 ---
    ax.set_xlabel('X (mm)')
    ax.set_ylabel('Y (mm)')
    ax.set_zlabel('Z (mm)')
    ax.set_title(title)

    # 设置大致的坐标轴范围，可以根据需要调整
    max_range = max(b11, b22, b33, np.linalg.norm(rp) + max(a11,a22,a33)) + 50
    ax.set_xlim(-max_range, max_range)
    ax.set_ylim(-max_range, max_range)
    ax.set_zlim(0, max(e+500, np.linalg.norm(rp)+100)) # Z 轴通常从 0 开始

    # 尝试设置等比例轴 (在 3D 中可能效果有限)
    try:
        ax.set_aspect('equal', adjustable='box')
    except NotImplementedError:
        print("警告: 3D 等比例轴可能不受支持。使用自动缩放。")
        # ax.axis('auto') # 默认行为
        # 手动设置比例接近
        x_limits = ax.get_xlim()
        y_limits = ax.get_ylim()
        z_limits = ax.get_zlim()
        x_range = abs(x_limits[1] - x_limits[0])
        y_range = abs(y_limits[1] - y_limits[0])
        z_range = abs(z_limits[1] - z_limits[0])
        plot_radius = 0.5 * max([x_range, y_range, z_range])
        mid_x = 0.5 * (x_limits[0] + x_limits[1])
        mid_y = 0.5 * (y_limits[0] + y_limits[1])
        mid_z = 0.5 * (z_limits[0] + z_limits[1])
        ax.set_xlim(mid_x - plot_radius, mid_x + plot_radius)
        ax.set_ylim(mid_y - plot_radius, mid_y + plot_radius)
        ax.set_zlim(mid_z - plot_radius, mid_z + plot_radius)


    ax.legend()
    ax.grid(True)

# %% --- 主程序与测试用例 ---

if __name__ == "__main__":

    # 获取机器人参数
    robot_params = get_robot_params()
    print("--- 机器人参数 ---")
    for key, value in robot_params.items():
        print(f"{key}: {value}")

    # === 测试 1: 单点逆解与正解验证 ===
    print("\n\n=== 测试 1: 单点逆解与正解 ===")
    # 目标末端位置 (与原始代码相同)
    xp_target = 10.0
    yp_target = -20.0
    zp_target = 1234.0
    print(f"\n--- 目标位置 (IK Input) ---")
    print(f"xp = {xp_target:.4f}")
    print(f"yp = {yp_target:.4f}")
    print(f"zp = {zp_target:.4f}")

    # --- 运行逆解 ---
    ik_result = inverse_kinematics(xp_target, yp_target, zp_target, robot_params)

    if ik_result['error']:
        print(f"\n逆解错误: {ik_result['error']}")
    else:
        print("\n--- 逆解结果 (IK Output) ---")
        print(f"计算得到的 q1: {ik_result['q1']:.4f} (注意: 此 q1 定义可能非标准)")
        print(f"计算得到的 q2: {ik_result['q2']:.4f} (注意: 此 q2 定义可能非标准)")
        print(f"计算得到的 q3: {ik_result['q3']:.4f} (注意: 此 q3 定义可能非标准)")
        print(f"计算得到的 q4 (中心杆可变长度): {ik_result['q4']:.4f}")
        print(f"计算得到的 theta1 (rad): {ik_result['theta1']:.4f}")
        print(f"计算得到的 theta2 (rad): {ik_result['theta2']:.4f}")

        # --- 使用逆解结果运行正解 ---
        q1_ik = ik_result['q1']
        q2_ik = ik_result['q2']
        q3_ik = ik_result['q3']

        # 正解需要初始猜测值，可以使用逆解结果加扰动，或一个通用值
        # 使用逆解结果作为猜测值通常收敛最快
        fk_initial_guess = [ik_result['q4'], ik_result['theta1'], ik_result['theta2']]
        # 或者使用一个固定的猜测值
        # fk_initial_guess = [1000.0, 0.1, 0.1]

        print("\n--- 运行正解 (FK) ---")
        print(f"输入 q1={q1_ik:.4f}, q2={q2_ik:.4f}, q3={q3_ik:.4f}")
        print(f"初始猜测 [q4, t1, t2]: [{fk_initial_guess[0]:.4f}, {fk_initial_guess[1]:.4f}, {fk_initial_guess[2]:.4f}]")

        # 显示符号方程 (只显示一次)
        if 'fk_sym_defined' not in globals(): # 检查变量是否已定义
             Func_sym, J_sym, _, _, X_sym = define_symbolic_fk(robot_params, q1_ik, q2_ik, q3_ik)
             print("\n正解约束方程 F(q4, t1, t2) = 0:")
             sp.pprint(Func_sym)
             print("\n正解雅可比矩阵 J = dF/dX:")
             sp.pprint(J_sym)
             fk_sym_defined = True # 标记已定义

        fk_result = forward_kinematics(q1_ik, q2_ik, q3_ik, robot_params, fk_initial_guess)

        if fk_result['error']:
            print(f"\n正解错误: {fk_result['error']}")
            print(f"迭代次数: {fk_result['iterations']}")
            print(f"最终误差范数: {fk_result['final_error_norm']:.6e}")
        else:
            print("\n--- 正解结果 (FK Output) ---")
            print(f"迭代次数: {fk_result['iterations']}")
            print(f"是否收敛: {fk_result['converged']}")
            print(f"最终误差范数: {fk_result['final_error_norm']:.6e}")
            print(f"计算得到的 q4: {fk_result['q4_k']:.4f}")
            print(f"计算得到的 theta1 (rad): {fk_result['theta1_k']:.4f}")
            print(f"计算得到的 theta2 (rad): {fk_result['theta2_k']:.4f}")
            print(f"计算得到的 xp: {fk_result['xp_fk']:.4f}")
            print(f"计算得到的 yp: {fk_result['yp_fk']:.4f}")
            print(f"计算得到的 zp: {fk_result['zp_fk']:.4f}")

            # --- 对比 ---
            error_fk = np.sqrt((fk_result['xp_fk'] - xp_target)**2 +
                               (fk_result['yp_fk'] - yp_target)**2 +
                               (fk_result['zp_fk'] - zp_target)**2)
            print("\n--- 对比: 目标位置 vs 正解位置 ---")
            print(f"目标 (IK Input): xp={xp_target:.4f}, yp={yp_target:.4f}, zp={zp_target:.4f}")
            print(f"正解 (FK Output): xp={fk_result['xp_fk']:.4f}, yp={fk_result['yp_fk']:.4f}, zp={fk_result['zp_fk']:.4f}")
            print(f"位置误差 (欧氏距离): {error_fk:.6f}")

        # --- 可视化单点结果 ---
        print("\n--- 绘制单点姿态 ---")
        fig_single = plt.figure(figsize=(8, 8))
        ax_single = fig_single.add_subplot(111, projection='3d')
        plot_robot(ax_single, ik_result['rp'], ik_result['R4'], robot_params,
                   title=f"Robot Pose at ({xp_target:.1f}, {yp_target:.1f}, {zp_target:.1f})")
        plt.show(block=False) # 显示图形，非阻塞

    # === 测试 2: 轨迹示例 - 垂直升降 ===
    print("\n\n=== 测试 2: 轨迹示例 - 垂直升降 ===")
    fig_traj = plt.figure(figsize=(8, 8))
    ax_traj = fig_traj.add_subplot(111, projection='3d')
    plt.ion() # 开启交互模式，用于动画

    # 轨迹参数
    x_traj = 50.0
    y_traj = 50.0
    z_start = robot_params['e'] + 100 # 起始高度，略高于中心柱
    z_end = z_start + 600           # 结束高度
    num_steps = 50                  # 轨迹点数量
    pause_time = 0.05               # 每帧暂停时间

    print(f"轨迹: 从 ({x_traj}, {y_traj}, {z_start:.1f}) 到 ({x_traj}, {y_traj}, {z_end:.1f})")

    trajectory_points = [] # 存储轨迹点用于绘制

    try:
        for i in range(num_steps + 1):
            # 计算当前目标位置
            z_current = z_start + (z_end - z_start) * (i / num_steps)
            xp_traj = x_traj
            yp_traj = y_traj
            zp_traj = z_current

            # 运行逆解
            ik_traj_result = inverse_kinematics(xp_traj, yp_traj, zp_traj, robot_params)

            if ik_traj_result['error']:
                print(f"轨迹点 {i} 逆解错误: {ik_traj_result['error']} @ Z={zp_traj:.1f}")
                time.sleep(1) # 暂停一下让用户看到错误
                continue # 跳过这个点

            # 存储轨迹点
            trajectory_points.append(ik_traj_result['rp'].flatten())

            # 更新绘图
            plot_robot(ax_traj, ik_traj_result['rp'], ik_traj_result['R4'], robot_params,
                       title=f"Trajectory Step {i}/{num_steps} - Z={zp_traj:.1f}")

            # 绘制已经过的轨迹线
            if len(trajectory_points) > 1:
                traj_array = np.array(trajectory_points)
                ax_traj.plot(traj_array[:, 0], traj_array[:, 1], traj_array[:, 2], 'm:', label='Trajectory Path' if i==1 else "") # 只加一次标签

            plt.draw()
            plt.pause(pause_time)

    except KeyboardInterrupt:
        print("\n轨迹动画被用户中断。")
    finally:
        plt.ioff() # 关闭交互模式
        print("\n轨迹动画结束。关闭绘图窗口以继续。")
        plt.show() # 保持最后一个窗口打开，直到用户关闭


    # === 测试 3: (可选) 轨迹示例 - 水平圆周 ===
    # print("\n\n=== 测试 3: 轨迹示例 - 水平圆周 ===")
    # ... (类似测试 2 的结构，但改变 xp_traj 和 yp_traj 的计算方式)
    # radius = 200
    # z_circle = robot_params['e'] + 400
    # num_steps_circle = 60
    # fig_circle = plt.figure(figsize=(8, 8))
    # ax_circle = fig_circle.add_subplot(111, projection='3d')
    # plt.ion()
    # trajectory_circle = []
    # try:
    #     for i in range(num_steps_circle + 1):
    #         angle = 2 * np.pi * i / num_steps_circle
    #         xp_traj = radius * np.cos(angle)
    #         yp_traj = radius * np.sin(angle)
    #         zp_traj = z_circle
    #         # ... (运行IK, 绘图, 暂停) ...
    # except KeyboardInterrupt:
    #     print("\n圆周轨迹动画被用户中断。")
    # finally:
    #     plt.ioff()
    #     print("\n圆周轨迹动画结束。")
    #     plt.show()