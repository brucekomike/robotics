# %% Import Libraries (Existing + New)
import numpy as np
import sympy as sp
import math
import matplotlib.pyplot as plt
# mpl_toolkits is no longer needed for Axes3D in recent matplotlib versions
# from mpl_toolkits.mplot3d import Axes3D # For 3D plotting
import time # For animation pause

# Use SymPy's pretty printing
sp.init_printing(use_unicode=True)

# --- Paste your existing functions here ---
# get_robot_params()
# inverse_kinematics(xp, yp, zp, params)
# define_symbolic_fk(params, q1_in, q2_in, q3_in)
# forward_kinematics(q1_in, q2_in, q3_in, params, initial_guess, ...)
# plot_robot(ax, rp, R4, params, title="Robot Configuration")
# --- End of pasted functions ---

# %% 机器人参数 (Robot Parameters) - Copied from your code
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

# %% 位置逆解函数 (Inverse Kinematics Function) - Copied from your code
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

    # Check if target is reachable (outside the fixed central column)
    # Add a small tolerance for floating point comparisons
    min_dist = e
    if norm_rp < min_dist - 1e-9:
        return {'error': f"警告：末端位置 (范数 {norm_rp:.2f}) 在中心柱固定长度 (e={e}) 范围内，无法达到。"}

    q4 = norm_rp - e
    # Ensure q4 (variable length) is not negative
    if q4 < 0:
        # This case should ideally be caught by the norm_rp check above,
        # but as a safeguard:
        q4 = 0
        # Recalculate norm_rp if we force q4=0, meaning rp is exactly on the sphere of radius e
        norm_rp = e
        if norm_rp < 1e-9: # Avoid division by zero if rp is at the origin
             return {'error': "错误：末端位置在原点，姿态未定义。"}
        rp = (rp / np.linalg.norm(rp)) * e # Project onto the sphere surface


    # Avoid division by zero when calculating angles
    denominator = q4 + e # This is norm_rp
    if abs(denominator) < 1e-9:
        # This should only happen if rp is the zero vector, handled above.
        return {'error': "错误：q4 + e 接近零 (rp 在原点)，无法计算角度。"}

    # Calculate theta2
    # Clip argument to [-1, 1] to handle potential floating point inaccuracies
    sin_theta2_arg = xp / denominator
    sin_theta2_arg = np.clip(sin_theta2_arg, -1.0, 1.0)
    theta2 = np.arcsin(sin_theta2_arg)

    # Calculate theta1
    cos_theta2 = np.cos(theta2)
    if abs(cos_theta2) < 1e-9:
        # When theta2 is close to +/- pi/2, the platform is near horizontal.
        # Use atan2 for robustness, considering the signs of -yp and zp.
        theta1 = np.arctan2(-yp, zp)
    else:
        # Calculate arguments for atan2 carefully
        arg1 = (-yp / denominator) / cos_theta2
        arg2 = (zp / denominator) / cos_theta2
        # Clip arguments for safety before atan2, although atan2 handles quadrants correctly
        # arg1 = np.clip(arg1, -1.0, 1.0) # Sin component
        # arg2 = np.clip(arg2, -1.0, 1.0) # Cos component
        theta1 = np.arctan2(arg1, arg2)


    # Calculate rotation matrix R4
    c1, s1 = np.cos(theta1), np.sin(theta1)
    c2, s2 = np.cos(theta2), np.sin(theta2)
    R4 = np.array([
        [c2,   0,  s2],
        [s1*s2, c1, -s1*c2],
        [-c1*s2, s1, c1*c2]
    ])

    # Fixed platform connection points (in fixed frame)
    b1 = np.array([[0], [-b11], [0]])
    b2 = np.array([[b22], [0], [0]])
    b3 = np.array([[-b33], [0], [0]])

    # Moving platform connection points (in moving frame {4})
    a10 = np.array([[0], [-a11], [0]])
    a20 = np.array([[a22], [0], [0]])
    a30 = np.array([[-a33], [0], [0]])

    # Moving platform connection points relative vector (in fixed frame)
    a1 = R4 @ a10
    a2 = R4 @ a20
    a3 = R4 @ a30

    # Central column direction vector w4 (in fixed frame)
    # w4 should be rp / norm(rp)
    if norm_rp > 1e-9:
        w4 = rp / norm_rp
    else:
        # Handle the case where rp is at the origin (though ideally prevented earlier)
        # Default to pointing up Z-axis, though this case is problematic
        w4 = np.array([[0], [0], [1.0]])


    # Calculate driving rod lengths q1, q2, q3 (using the formula from the original code)
    # Note: This definition might differ from the physical length || (rp + a_i) - b_i ||
    # The formula seems to represent the length of the vector from B_i to A_i'
    # where A_i' is the connection point on the moving platform projected down
    # along the central column's variable part.
    # L_i = (rp + a_i) - b_i = (e*w4 + q4*w4 + a_i) - b_i
    # The formula used is || q4*w4 + a_i - b_i || = || rp - e*w4 + a_i - b_i ||
    q1_vec = rp - b1 + a1 - e*w4
    q2_vec = rp - b2 + a2 - e*w4
    q3_vec = rp - b3 + a3 - e*w4
    q1 = np.linalg.norm(q1_vec)
    q2 = np.linalg.norm(q2_vec)
    q3 = np.linalg.norm(q3_vec)

    return {
        'q1': q1, 'q2': q2, 'q3': q3,
        'q4': q4, 'theta1': theta1, 'theta2': theta2,
        'rp': rp, 'R4': R4, 'w4': w4,
        'error': None
    }


# %% 位置正解 - 符号设置与数值函数生成 (Forward Kinematics - Symbolic Setup) - Copied from your code
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
    # Based on the inverse kinematics formula q_i = || q4*w4 + a_i - b_i ||
    # where w4 = [s2; -s1*c2; c1*c2]
    # a1 = R4*a10 = R4*[0; -a11; 0] = [-a11*0; a11*c1; a11*s1] ? No, R4 is complex.
    # Let's use the formula provided directly in the MATLAB symbolic part if possible,
    # assuming it corresponds to the || q4*w4 + a_i - b_i ||^2 = q_i^2 relation.

    # R4 = [[c2, 0, s2], [s1*s2, c1, -s1*c2], [-c1*s2, s1, c1*c2]]
    # a10 = [0; -a11; 0] -> a1 = R4*a10 = [0; -a11*c1; -a11*s1]
    # a20 = [a22; 0; 0] -> a2 = R4*a20 = [a22*c2; a22*s1*s2; -a22*c1*s2]
    # a30 = [-a33; 0; 0] -> a3 = R4*a30 = [-a33*c2; -a33*s1*s2; a33*c1*s2]
    # w4 = [s2; -s1*c2; c1*c2]
    # b1 = [0; -b11; 0]
    # b2 = [b22; 0; 0]
    # b3 = [-b33; 0; 0]

    # Vector L1 = q4*w4 + a1 - b1
    # L1 = [q4*s2 + 0 - 0;
    #      q4*(-s1*c2) - a11*c1 - (-b11);
    #      q4*(c1*c2) - a11*s1 - 0]
    # L1 = [q4*s2; -q4*s1*c2 - a11*c1 + b11; q4*c1*c2 - a11*s1]
    # ||L1||^2 = (q4*s2)^2 + (-q4*s1*c2 - a11*c1 + b11)^2 + (q4*c1*c2 - a11*s1)^2 = q1^2

    # Vector L2 = q4*w4 + a2 - b2
    # L2 = [q4*s2 + a22*c2 - b22;
    #      q4*(-s1*c2) + a22*s1*s2 - 0;
    #      q4*(c1*c2) - a22*c1*s2 - 0]
    # L2 = [q4*s2 + a22*c2 - b22; -q4*s1*c2 + a22*s1*s2; q4*c1*c2 - a22*c1*s2]
    # ||L2||^2 = (q4*s2 + a22*c2 - b22)^2 + (-q4*s1*c2 + a22*s1*s2)^2 + (q4*c1*c2 - a22*c1*s2)^2 = q2^2

    # Vector L3 = q4*w4 + a3 - b3
    # L3 = [q4*s2 - a33*c2 - (-b33);
    #      q4*(-s1*c2) - a33*s1*s2 - 0;
    #      q4*(c1*c2) + a33*c1*s2 - 0]
    # L3 = [q4*s2 - a33*c2 + b33; -q4*s1*c2 - a33*s1*s2; q4*c1*c2 + a33*c1*s2]
    # ||L3||^2 = (q4*s2 - a33*c2 + b33)^2 + (-q4*s1*c2 - a33*s1*s2)^2 + (q4*c1*c2 + a33*c1*s2)^2 = q3^2

    # Let's check the provided F_sym, G_sym, H_sym if they match the expansion of ||Li||^2 - qi^2 = 0
    # F_sym = a11**2 + b11**2 + q4_sym**2 - q1_in**2 \
    #         - 2*a11*b11*sp.cos(theta1_sym) \
    #         - 2*b11*q4_sym*sp.sin(theta1_sym)*sp.cos(theta2_sym)
    # G_sym = a22**2 + b22**2 + q4_sym**2 - q2_in**2 \
    #         - 2*b22*q4_sym*sp.sin(theta2_sym) \
    #         - 2*a22*b22*sp.cos(theta2_sym)
    # H_sym = a33**2 + b33**2 + q4_sym**2 - q3_in**2 \
    #         + 2*b33*q4_sym*sp.sin(theta2_sym) \
    #         - 2*a33*b33*sp.cos(theta2_sym)

    # These symbolic equations seem simpler than the full expansion of ||Li||^2.
    # They might be derived from a different geometric constraint or simplification.
    # IMPORTANT: We will use the F, G, H equations provided in the original code,
    # assuming they are correct for the specific robot mechanism intended.
    c1_sym, s1_sym = sp.cos(theta1_sym), sp.sin(theta1_sym)
    c2_sym, s2_sym = sp.cos(theta2_sym), sp.sin(theta2_sym)

    F_sym = a11**2 + b11**2 + q4_sym**2 - q1_in**2 \
            - 2*a11*b11*c1_sym \
            - 2*b11*q4_sym*s1_sym*c2_sym # Check signs: w4_y = -s1*c2. Term is -2*b1_y*(q4*w4_y)? No.

    # Let's re-verify the terms in F_sym based on ||L1||^2
    # ||L1||^2 = (q4*s2)**2 + (-q4*s1*c2 - a11*c1 + b11)**2 + (q4*c1*c2 - a11*s1)**2
    # = q4^2*s2^2
    # + (q4*s1*c2 + a11*c1 - b11)^2 + (q4*c1*c2 - a11*s1)^2
    # + q4^2*s1^2*c2^2 + a11^2*c1^2 + b11^2 + 2*q4*s1*c2*a11*c1 - 2*q4*s1*c2*b11 - 2*a11*c1*b11
    # + q4^2*c1^2*c2^2 - 2*q4*c1*c2*a11*s1 + a11^2*s1^2
    # Group terms:
    # q4^2 terms: q4^2*s2^2 + q4^2*s1^2*c2^2 + q4^2*c1^2*c2^2 = q4^2*s2^2 + q4^2*c2^2*(s1^2+c1^2) = q4^2*s2^2 + q4^2*c2^2 = q4^2*(s2^2+c2^2) = q4^2
    # a11^2 terms: a11^2*c1^2 + a11^2*s1^2 = a11^2*(c1^2+s1^2) = a11^2
    # b11^2 terms: b11^2
    # q4*a11 terms: 2*q4*s1*c2*a11*c1 - 2*q4*c1*c2*a11*s1 = 0
    # q4*b11 terms: -2*q4*s1*c2*b11
    # a11*b11 terms: -2*a11*c1*b11
    # So, ||L1||^2 = q4^2 + a11^2 + b11^2 - 2*q4*s1*c2*b11 - 2*a11*c1*b11
    # Equation F = ||L1||^2 - q1^2 = 0 becomes:
    # q4^2 + a11^2 + b11^2 - 2*b11*q4*s1*c2 - 2*a11*b11*c1 - q1^2 = 0
    # Comparing with the provided F_sym:
    # F_sym = a11**2 + b11**2 + q4_sym**2 - q1_in**2 - 2*a11*b11*c1_sym - 2*b11*q4_sym*s1_sym*c2_sym
    # They MATCH! Good.

    # Now check G_sym based on ||L2||^2
    # L2 = [q4*s2 + a22*c2 - b22; -q4*s1*c2 + a22*s1*s2; q4*c1*c2 - a22*c1*s2]
    # ||L2||^2 = (q4*s2 + a22*c2 - b22)^2 + (-q4*s1*c2 + a22*s1*s2)^2 + (q4*c1*c2 - a22*c1*s2)^2
    # = (q4*s2 + a22*c2 - b22)^2 + s1^2*(-q4*c2 + a22*s2)^2 + c1^2*(q4*c2 - a22*s2)^2
    # = (q4*s2 + a22*c2 - b22)^2 + (s1^2+c1^2)*(q4*c2 - a22*s2)^2
    # = (q4*s2 + a22*c2 - b22)^2 + (q4*c2 - a22*s2)^2
    # Expand:
    # (q4*s2)^2 + (a22*c2)^2 + b22^2 + 2*q4*s2*a22*c2 - 2*q4*s2*b22 - 2*a22*c2*b22
    # + (q4*c2)^2 - 2*q4*c2*a22*s2 + (a22*s2)^2
    # Group terms:
    # q4^2 terms: q4^2*s2^2 + q4^2*c2^2 = q4^2
    # a22^2 terms: a22^2*c2^2 + a22^2*s2^2 = a22^2
    # b22^2 terms: b22^2
    # q4*a22 terms: 2*q4*s2*a22*c2 - 2*q4*c2*a22*s2 = 0
    # q4*b22 terms: -2*q4*s2*b22
    # a22*b22 terms: -2*a22*c2*b22
    # So, ||L2||^2 = q4^2 + a22^2 + b22^2 - 2*q4*s2*b22 - 2*a22*c2*b22
    # Equation G = ||L2||^2 - q2^2 = 0 becomes:
    # q4^2 + a22^2 + b22^2 - 2*b22*q4*s2 - 2*a22*b22*c2 - q2^2 = 0
    # Comparing with provided G_sym:
    # G_sym = a22**2 + b22**2 + q4_sym**2 - q2_in**2 - 2*b22*q4_sym*s2_sym - 2*a22*b22*c2_sym
    # They MATCH! Good.

    # Now check H_sym based on ||L3||^2
    # L3 = [q4*s2 - a33*c2 + b33; -q4*s1*c2 - a33*s1*s2; q4*c1*c2 + a33*c1*s2]
    # ||L3||^2 = (q4*s2 - a33*c2 + b33)^2 + (-q4*s1*c2 - a33*s1*s2)^2 + (q4*c1*c2 + a33*c1*s2)^2
    # = (q4*s2 - a33*c2 + b33)^2 + s1^2*(-q4*c2 - a33*s2)^2 + c1^2*(q4*c2 + a33*s2)^2
    # = (q4*s2 - a33*c2 + b33)^2 + (s1^2+c1^2)*(q4*c2 + a33*s2)^2
    # = (q4*s2 - a33*c2 + b33)^2 + (q4*c2 + a33*s2)^2
    # Expand:
    # (q4*s2)^2 + (a33*c2)^2 + b33^2 - 2*q4*s2*a33*c2 + 2*q4*s2*b33 - 2*a33*c2*b33
    # + (q4*c2)^2 + 2*q4*c2*a33*s2 + (a33*s2)^2
    # Group terms:
    # q4^2 terms: q4^2*s2^2 + q4^2*c2^2 = q4^2
    # a33^2 terms: a33^2*c2^2 + a33^2*s2^2 = a33^2
    # b33^2 terms: b33^2
    # q4*a33 terms: -2*q4*s2*a33*c2 + 2*q4*c2*a33*s2 = 0
    # q4*b33 terms: 2*q4*s2*b33
    # a33*b33 terms: -2*a33*c2*b33
    # So, ||L3||^2 = q4^2 + a33^2 + b33^2 + 2*q4*s2*b33 - 2*a33*c2*b33
    # Equation H = ||L3||^2 - q3^2 = 0 becomes:
    # q4^2 + a33^2 + b33^2 + 2*b33*q4*s2 - 2*a33*b33*c2 - q3^2 = 0
    # Comparing with provided H_sym:
    # H_sym = a33**2 + b33**2 + q4_sym**2 - q3_in**2 + 2*b33*q4_sym*s2_sym - 2*a33*b33*c2_sym
    # They MATCH! Excellent.

    # The symbolic equations provided are correct based on the IK definition q_i = || rp - e*w4 + a_i - b_i ||

    F_sym = a11**2 + b11**2 + q4_sym**2 - q1_in**2 \
            - 2*a11*b11*c1_sym \
            - 2*b11*q4_sym*s1_sym*c2_sym

    G_sym = a22**2 + b22**2 + q4_sym**2 - q2_in**2 \
            - 2*b22*q4_sym*s2_sym \
            - 2*a22*b22*c2_sym

    H_sym = a33**2 + b33**2 + q4_sym**2 - q3_in**2 \
            + 2*b33*q4_sym*s2_sym \
            - 2*a33*b33*c2_sym

    Func_sym = sp.Matrix([F_sym, G_sym, H_sym])

    # 计算雅可比矩阵 J = d(Func)/d(X)
    J_sym = Func_sym.jacobian(X_sym)

    # 转换为快速数值函数
    num_func = sp.lambdify((q4_sym, theta1_sym, theta2_sym), Func_sym, 'numpy')
    num_jac = sp.lambdify((q4_sym, theta1_sym, theta2_sym), J_sym, 'numpy')

    return Func_sym, J_sym, num_func, num_jac, X_sym

# %% 位置正解 - 牛顿迭代求解器 (Forward Kinematics - Newton Solver) - Copied from your code
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
    # Note: Defining symbolic functions every call can be slow.
    # Consider defining them once outside if performance is critical.
    Func_sym, J_sym, num_func, num_jac, X_sym = define_symbolic_fk(params, q1_in, q2_in, q3_in)

    # 初始猜测值
    # Ensure initial_guess is a mutable numpy array for in-place updates
    X_k = np.array(initial_guess, dtype=float)

    iterations = 0
    converged = False
    final_error_norm = float('inf')
    error_message = None

    for i in range(max_iter):
        iterations = i
        q4_k, theta1_k, theta2_k = X_k[0], X_k[1], X_k[2]

        # Ensure angles stay within reasonable bounds if needed, e.g., -pi to pi for theta1
        # X_k[1] = (X_k[1] + np.pi) % (2 * np.pi) - np.pi
        # However, Newton's method might jump outside, rely on trig functions in F, J

        # 计算当前点的函数值 F(X_k) (需要是列向量 for numpy solve, but lambdify gives matrix)
        try:
            F_k_mat = num_func(q4_k, theta1_k, theta2_k)
            # Ensure F_k is a flat array (vector)
            F_k = np.array(F_k_mat).flatten()
            if F_k.shape != (3,):
                 raise ValueError(f"Function output shape mismatch: {F_k.shape}")

        except (TypeError, ValueError, NameError) as e_eval:
            error_message = f"Error evaluating F at iteration {i}: {e_eval}. Check inputs/lambdify."
            # print(f"\n{error_message}")
            # print(f"X_k = {X_k}")
            break # Exit loop on evaluation error

        # 检查收敛性
        current_error_norm = np.linalg.norm(F_k)
        if current_error_norm < tolerance:
            converged = True
            final_error_norm = current_error_norm
            break

        # 计算当前点的雅可比矩阵 J(X_k)
        try:
            J_k = num_jac(q4_k, theta1_k, theta2_k)
             # Ensure J_k is a 2D numpy array
            J_k = np.array(J_k).astype(float)
            if J_k.shape != (3, 3):
                 raise ValueError(f"Jacobian output shape mismatch: {J_k.shape}")

        except (TypeError, ValueError, NameError) as e_eval:
            error_message = f"Error evaluating J at iteration {i}: {e_eval}. Check inputs/lambdify."
            # print(f"\n{error_message}")
            # print(f"X_k = {X_k}")
            break # Exit loop on evaluation error


        # 求解线性方程组 J(X_k) * delta_X = -F(X_k)
        try:
            # Use pseudo-inverse for potentially singular Jacobian (more robust but slower)
            # delta_X = np.linalg.pinv(J_k) @ (-F_k)
            # Or use standard solve and catch singularity
            delta_X = np.linalg.solve(J_k, -F_k)

            # Add damping factor if steps are too large (optional)
            # step_norm = np.linalg.norm(delta_X)
            # max_step = 1.0 # Adjust as needed
            # if step_norm > max_step:
            #     delta_X = delta_X * (max_step / step_norm)

        except np.linalg.LinAlgError:
            error_message = f"Jacobian matrix singular at iteration {i}. Cannot solve."
            # print(f"\n{error_message}")
            # print(f"J_k = {J_k}")
            final_error_norm = current_error_norm # Record error before failing
            converged = False
            break # Exit loop if singular

        # 更新变量 X_{k+1} = X_k + delta_X
        X_k += delta_X

        # Check for NaN or Inf values
        if not np.all(np.isfinite(X_k)):
            error_message = f"Non-finite values encountered at iteration {i}. Divergence?"
            # print(f"\n{error_message}")
            converged = False
            break


    # After loop finishes (break or max_iter)
    if converged:
        error_message = None
        iterations += 1 # Increment iterations to count the last successful one
    elif error_message is None: # Reached max_iter without converging or explicit error
        error_message = "Reached maximum iterations without converging."
        final_error_norm = np.linalg.norm(num_func(X_k[0], X_k[1], X_k[2]).flatten())

    # Calculate final end effector position if converged
    xp_fk, yp_fk, zp_fk = None, None, None
    if converged:
        q4_final, theta1_final, theta2_final = X_k
        # Use the forward kinematics equations: rp = (e + q4) * w4
        norm_rp_final = e + q4_final
        c1_f, s1_f = np.cos(theta1_final), np.sin(theta1_final)
        c2_f, s2_f = np.cos(theta2_final), np.sin(theta2_final)
        w4_final = np.array([s2_f, -s1_f*c2_f, c1_f*c2_f])
        rp_final = norm_rp_final * w4_final
        xp_fk, yp_fk, zp_fk = rp_final[0], rp_final[1], rp_final[2]

    return {
        'xp_fk': xp_fk, 'yp_fk': yp_fk, 'zp_fk': zp_fk,
        'q4_k': X_k[0], 'theta1_k': X_k[1], 'theta2_k': X_k[2],
        'iterations': iterations, 'final_error_norm': final_error_norm,
        'converged': converged, 'error': error_message
    }

# %% 可视化函数 (Visualization Function) - Copied from your code
def plot_robot(ax, rp, R4, params, title="Robot Configuration"):
    """
    在给定的 3D 坐标轴上绘制机器人结构。
    输入:
        ax: Matplotlib 3D 坐标轴对象 (Axes3D)
        rp: 末端位置向量 (numpy.ndarray, shape (3,1) or (3,))
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

    # Ensure rp is a flat array for easier indexing
    P = np.array(rp).flatten()
    if P.shape != (3,):
        print(f"Warning: Unexpected rp shape in plot_robot: {rp.shape}")
        # Attempt to reshape if it's (3,1)
        if P.shape == (3, 1):
            P = P.flatten()
        else:
            # Cannot plot if shape is wrong
            ax.set_title(f"{title}\n(Error: Invalid rp input)")
            return

    # --- 计算关键点坐标 ---
    # 原点 O
    O = np.array([0, 0, 0])
    # 固定平台连接点 B1, B2, B3 (在固定坐标系)
    B1 = np.array([0, -b11, 0])
    B2 = np.array([b22, 0, 0])
    B3 = np.array([-b33, 0, 0])
    # 移动平台连接点 A10, A20, A30 (在移动平台自身坐标系 {4})
    A10 = np.array([0, -a11, 0])
    A20 = np.array([a22, 0, 0])
    A30 = np.array([-a33, 0, 0])
    # 移动平台连接点 A1, A2, A3 (相对平台中心 rp 的向量，在固定坐标系)
    # Ensure R4 is a 3x3 matrix
    if R4.shape != (3, 3):
        print(f"Warning: Unexpected R4 shape in plot_robot: {R4.shape}")
        ax.set_title(f"{title}\n(Error: Invalid R4 input)")
        return

    A1_vec = R4 @ A10.reshape(3, 1)
    A2_vec = R4 @ A20.reshape(3, 1)
    A3_vec = R4 @ A30.reshape(3, 1)
    # 移动平台连接点 A1', A2', A3' 的绝对坐标 (在固定坐标系)
    A1_prime = P + A1_vec.flatten()
    A2_prime = P + A2_vec.flatten()
    A3_prime = P + A3_vec.flatten()

    # 中心柱顶端 E (固定坐标系)
    norm_P = np.linalg.norm(P)
    if norm_P > 1e-9:
        w4_vis = P / norm_P
    else: # Handle case where P is at the origin
        w4_vis = np.array([0, 0, 1.0]) # Default to Z-axis
    E_pt = e * w4_vis # Coordinates of the top of the fixed part

    # --- 清除之前的绘图 ---
    ax.cla()

    # --- 绘制结构 ---
    # 固定平台 (灰色虚线)
    fixed_platform_pts = np.array([B1, B2, B3, B1])
    ax.plot(fixed_platform_pts[:, 0], fixed_platform_pts[:, 1], fixed_platform_pts[:, 2], 'k--', lw=1, color='gray', label='_Fixed Platform Base') # Underscore hides from legend unless specified
    ax.scatter(fixed_platform_pts[:-1, 0], fixed_platform_pts[:-1, 1], fixed_platform_pts[:-1, 2], c='black', marker='s', s=40, label='Base Joints (B1-B3)')

    # 移动平台 (蓝色实线)
    moving_platform_pts = np.array([A1_prime, A2_prime, A3_prime, A1_prime])
    ax.plot(moving_platform_pts[:, 0], moving_platform_pts[:, 1], moving_platform_pts[:, 2], 'b-', lw=2, label='Moving Platform')
    ax.scatter(moving_platform_pts[:-1, 0], moving_platform_pts[:-1, 1], moving_platform_pts[:-1, 2], c='blue', marker='o', s=40, label='Platform Joints (A1\'-A3\')')

    # 驱动杆 (红色实线) - 绘制物理连杆 B_i 到 A'_i
    # These are the physical links whose lengths might NOT be q1, q2, q3 from IK
    ax.plot([B1[0], A1_prime[0]], [B1[1], A1_prime[1]], [B1[2], A1_prime[2]], 'r-', lw=1.5, label='Driving Rods (Physical)')
    ax.plot([B2[0], A2_prime[0]], [B2[1], A2_prime[1]], [B2[2], A2_prime[2]], 'r-', lw=1.5)
    ax.plot([B3[0], A3_prime[0]], [B3[1], A3_prime[1]], [B3[2], A3_prime[2]], 'r-', lw=1.5)

    # 中心柱 (绿色)
    ax.plot([O[0], E_pt[0]], [O[1], E_pt[1]], [O[2], E_pt[2]], 'g-', linewidth=4, label='Central Column (Fixed Part)') # Fixed part
    ax.plot([E_pt[0], P[0]], [E_pt[1], P[1]], [E_pt[2], P[2]], 'g--', linewidth=2, label='Central Column (Variable Part)') # Variable part (length q4)

    # 标出末端点 P
    ax.scatter(P[0], P[1], P[2], c='magenta', marker='*', s=150, label='End Effector (P)', depthshade=False, zorder=10)
    # 标出原点 O
    ax.scatter(O[0], O[1], O[2], c='black', marker='x', s=50, label='Origin (O)')

    # --- 设置绘图属性 ---
    ax.set_xlabel('X (mm)')
    ax.set_ylabel('Y (mm)')
    ax.set_zlabel('Z (mm)')
    ax.set_title(title)

    # 设置大致的坐标轴范围，可以根据需要调整
    # Find max extent dynamically
    all_points = np.vstack([O, B1, B2, B3, P, A1_prime, A2_prime, A3_prime, E_pt])
    max_vals = np.max(all_points, axis=0)
    min_vals = np.min(all_points, axis=0)
    center = (max_vals + min_vals) / 2
    max_range = np.max(max_vals - min_vals)
    if max_range < 100: max_range = max(b11, b22, b33, e) # Ensure minimum range if points are close
    plot_radius = max_range * 0.6 # Add some padding

    ax.set_xlim(center[0] - plot_radius, center[0] + plot_radius)
    ax.set_ylim(center[1] - plot_radius, center[1] + plot_radius)
    # Ensure Z starts from 0 or slightly below
    z_min = min(0, min_vals[2]) - 50
    z_max = max(e + 500, max_vals[2]) + 50 # Ensure Z range is sufficient
    ax.set_zlim(z_min, z_max)

    # Attempt to set aspect ratio to 'equal' for a more realistic view
    try:
        # Use 'equal' for axes box, not data limits necessarily
        ax.set_aspect('equal', adjustable='box')
        # print("Applied 'equal' aspect ratio.")
    except NotImplementedError:
        # Fallback: Manually set limits to be cubic around the center
        # This was the logic previously used, keep it as fallback
        x_limits = ax.get_xlim()
        y_limits = ax.get_ylim()
        z_limits = ax.get_zlim()
        x_range = abs(x_limits[1] - x_limits[0])
        y_range = abs(y_limits[1] - y_limits[0])
        z_range = abs(z_limits[1] - z_limits[0])
        plot_radius_manual = 0.5 * max([x_range, y_range, z_range])
        mid_x = 0.5 * (x_limits[0] + x_limits[1])
        mid_y = 0.5 * (y_limits[0] + y_limits[1])
        mid_z = 0.5 * (z_limits[0] + z_limits[1])
        ax.set_xlim(mid_x - plot_radius_manual, mid_x + plot_radius_manual)
        ax.set_ylim(mid_y - plot_radius_manual, mid_y + plot_radius_manual)
        ax.set_zlim(mid_z - plot_radius_manual, mid_z + plot_radius_manual)
        print("Warning: 3D equal aspect ratio not fully supported. Attempted manual scaling.")


    ax.legend(fontsize='small')
    ax.grid(True)


# %% S-Curve Velocity Profile Generation (Replicating MATLAB's Svelocity_motion)
def s_velocity_motion(vmax, TA, sall, ts):
    """
    Generates an S-curve (trapezoidal acceleration) velocity profile.
    This function aims to replicate the specific MATLAB implementation provided,
    including its calculation of intermediate times and jerk.
    Note: The calculation of T1, T2, T4, T5 and je seems simplified/non-standard
          compared to typical 7-segment S-curve planning, but the piecewise
          formulas *are* used. We replicate the MATLAB code's behavior.

    Inputs:
        vmax: Maximum velocity (mm/s)
        TA: Total time for acceleration phase (s) - Used to derive amax and je
        sall: Total distance to travel (mm)
        ts: Time step (s)

    Outputs:
        ac: numpy array of acceleration values (mm/s^2)
        vc: numpy array of velocity values (mm/s)
        sc: numpy array of displacement values (mm)
        time_vec: numpy array of time points (s)
        params: Dictionary containing intermediate calculation parameters (T1-T5, je, etc.)
    """
    if vmax <= 0 or TA <= 0 or sall <= 0 or ts <= 0:
        raise ValueError("vmax, TA, sall, and ts must be positive.")

    # --- Calculations as per MATLAB code ---
    amax = vmax / TA  # Max acceleration derived from vmax and TA
    pi_val = np.pi # Use numpy's pi

    # These T calculations seem specific to the MATLAB code's logic
    T2 = vmax / amax # This just results in T2 = TA
    T1 = T2          # T1 = TA
    T4 = T2          # T4 = TA (Deceleration time)
    T5 = T1          # T5 = TA (Deceleration jerk time?)

    # Jerk calculation based on the derived amax and T1 (=TA)
    # This implies a relationship je = amax / T1 = (vmax/TA) / TA = vmax / TA^2
    je = amax / T1

    # Calculate constant velocity time T3
    # Check if the distance is sufficient for reaching full velocity
    s_accel_decel = vmax * TA # Distance covered during accel and decel phases in a simple trapezoid
                              # The S-curve formula used in MATLAB is different:
    s_ad_matlab = 2 * (je * T1**3 / 6 + je * T1**2 / 2 * T2 + amax * T2**2 / 2) # Approx?
    # Let's use the formula from MATLAB for T3 calculation:
    # T3 = (sall - 2 * je * T1^3) / vmax # This formula seems highly suspect.
    # Let's use the standard trapezoidal calculation for T3, consistent with total time calc:
    # T3 = (sall / vmax) - TA # Time at constant velocity if accel/decel phases take TA each.
    # Let's re-examine the MATLAB T3 calc: T3=(sall-2*je*T1^3)/vmax;
    # This looks related to s1 = (1/6)je*T1^3, maybe sall - 2*s1? Unlikely.
    # Let's use the T3 derived from total time calculation in MATLAB:
    # t5 = T1+T2+T3+T4+T5 => t5 = TA+TA+T3+TA+TA = 4*TA + T3
    # Also, standard total time = TA + T3 + TA = 2*TA + T3
    # If T3 = (sall/vmax) - 2*TA (standard trapezoid if TA is accel *and* decel time)
    # Then total time t5 = 2*TA + (sall/vmax - 2*TA) = sall/vmax. This is wrong.
    #
    # Let's assume TA is the time for *each* phase (accel phase, decel phase).
    # Total time = TA (accel) + T3 (const vel) + TA (decel)
    # Distance during accel = Distance during decel = s_a
    # sall = 2*s_a + vmax * T3
    # T3 = (sall - 2*s_a) / vmax
    # For a simple trapezoid, s_a = 0.5 * vmax * TA. T3 = (sall - vmax*TA)/vmax = sall/vmax - TA.
    # Total time = TA + (sall/vmax - TA) + TA = sall/vmax + TA.
    #
    # Let's look at the MATLAB time boundaries:
    # t1=T1=TA, t2=T1+T2=2*TA, t3=t2+T3, t4=t3+T4=t3+TA, t5=t4+T5=t4+TA=t3+2*TA
    # If we use T3 = sall/vmax - TA (from simple trapezoid total time)
    # t3 = 2*TA + sall/vmax - TA = TA + sall/vmax
    # t5 = t3 + 2*TA = TA + sall/vmax + 2*TA = 3*TA + sall/vmax. Still seems wrong.
    #
    # What if the MATLAB code intended TA to be the *jerk* time T_j?
    # Then amax = je * TA. vmax = amax * TA_const + ...
    #
    # **Reverting to exact MATLAB calculation replication:**
    # Assume the MATLAB calculation for T3 is intended, despite looking odd.
    # T3 = (sall - 2 * je * T1**3) / vmax # Use this formula directly.
    # Check if T3 is negative, meaning vmax is not reached (triangular profile needed)
    if T3 < 0:
        # This profile generation assumes vmax is reached.
        # A more robust implementation would handle the triangular case.
        print(f"Warning: Calculated T3 ({T3:.3f}) is negative. vmax may not be reached for given sall, vmax, TA.")
        print("Proceeding with T3=0 (assuming minimal constant velocity time). Results might be inaccurate.")
        # Adjust parameters for a triangular profile (approximate here)
        # Need to recalculate TA or vmax based on sall.
        # Or, simply set T3 = 0 and let the profile run.
        T3 = 0
        # Recalculate total time based on this T3.

    # Time boundaries based on MATLAB logic (T1=T2=T4=T5=TA)
    t1 = T1
    t2 = T1 + T2
    t3 = T1 + T2 + T3
    t4 = T1 + T2 + T3 + T4
    t5 = T1 + T2 + T3 + T4 + T5 # Total time
    # t1 = TA
    # t2 = 2*TA
    # t3 = 2*TA + T3
    # t4 = 3*TA + T3
    # t5 = 4*TA + T3

    # Store parameters
    params = {'vmax': vmax, 'TA': TA, 'sall': sall, 'ts': ts,
              'amax': amax, 'je': je, 'T1': T1, 'T2': T2, 'T3': T3, 'T4': T4, 'T5': T5,
              't1': t1, 't2': t2, 't3': t3, 't4': t4, 't5': t5}

    # --- Generate Time Vector ---
    # nt = int(np.fix((t5 * 1000) / (ts * 1000))) # MATLAB's nt calculation
    # A more direct way:
    nt = int(np.ceil(t5 / ts)) + 1 # Number of points including t=0 and t=t5
    time_vec = np.linspace(0, t5, nt)
    # Adjust actual time step slightly if needed to land exactly on t5
    actual_ts = t5 / (nt - 1) if nt > 1 else 0
    params['actual_ts'] = actual_ts
    params['nt'] = nt


    # --- Initialize Output Arrays ---
    ac = np.zeros(nt)
    vc = np.zeros(nt)
    sc = np.zeros(nt)

    # --- Loop through time and apply piecewise formulas ---
    # Use the formulas exactly as in the MATLAB snippet images 5 & 6
    for i, t in enumerate(time_vec):
        # Note: Using <= for boundaries to match typical closed intervals [t_start, t_end]
        # The MATLAB code uses combinations of <= and <. We'll try to match it.
        # MATLAB: if(t>=0 && t<=t1)
        if 0 <= t <= t1:
            ac[i] = je * t
            vc[i] = 0.5 * je * t**2
            sc[i] = (1/6) * je * t**3
        # MATLAB: else if (t>t1 && t<=t2)
        elif t1 < t <= t2:
            # Calculate values at t1 first
            vc_t1 = 0.5 * je * T1**2
            sc_t1 = (1/6) * je * T1**3
            # Apply formulas for t > t1
            ac[i] = amax # Constant acceleration phase?
            vc[i] = vc_t1 + amax * (t - T1)
            # sc[i] = sc_t1 + vc_t1 * (t - T1) + 0.5 * amax * (t - T1)**2 # Standard formula
            # MATLAB formula: sc(ntt)=(1/6)*je*T1^3+(1/2)*je*T1^2*(t-t1)+(1/2)*amax*(t-t1)^2;
            sc[i] = sc_t1 + 0.5 * je * T1**2 * (t - T1) + 0.5 * amax * (t - T1)**2
        # MATLAB: else if (t>t2 && t<=t3)
        elif t2 < t <= t3:
             # Calculate values at t2 first (end of previous segment)
            vc_t1 = 0.5 * je * T1**2
            sc_t1 = (1/6) * je * T1**3
            vc_t2 = vc_t1 + amax * (T2) # T2 is duration t2-t1
            # sc_t2 = sc_t1 + 0.5 * je * T1**2 * (T2) + 0.5 * amax * (T2)**2 # Using MATLAB sc formula
            # Let's recalculate sc_t2 based on vc integral:
            # sc_t2 = sc_t1 + integral(vc_t1 + amax*(tau-T1)) dtau from T1 to T2
            # sc_t2 = sc_t1 + vc_t1*(T2-T1) + 0.5*amax*(T2-T1)^2. Here T2-T1 = T2 (duration)
            sc_t2 = sc_t1 + vc_t1*T2 + 0.5*amax*T2**2

            # Apply formulas for t > t2
            ac[i] = 0 # Constant velocity phase
            vc[i] = vmax # Should be vmax if calculations are consistent
            # Verify vc_t2 == vmax
            # vc_t2 = 0.5*je*T1^2 + amax*T2 = 0.5*(amax/T1)*T1^2 + amax*T2
            #       = 0.5*amax*T1 + amax*T2
            # Since T1=T2=TA, vc_t2 = 0.5*amax*TA + amax*TA = 1.5*amax*TA = 1.5*vmax. This is wrong!
            # This confirms the initial T1,T2 calculation is inconsistent with the piecewise formulas
            # assuming they represent a standard profile reaching vmax.
            #
            # **Decision:** We MUST stick to replicating the MATLAB formulas exactly, even if inconsistent.
            vc[i] = vc_t2 # Use the calculated velocity at t2 as constant velocity? No, MATLAB uses vmax.
            vc[i] = vmax # Use vmax as specified in the MATLAB formula for this segment.

            # MATLAB formula: sc(ntt)=je*T1^3/6+vmax*T2+vmax*(t-t2); This looks very wrong.
            # It uses vmax*T2, which assumes constant velocity vmax during the *previous* segment.
            # Let's try the formula from Image 6:
            # sc(ntt)=je*T1^3/6 + je*T1^2/2*T2 + amax*T2^2/2 + vmax*(t-t2)
            # This is sc_t2 (calculated using MATLAB's formula) + vmax*(t-t2)
            sc_t2_matlab = (1/6)*je*T1**3 + 0.5*je*T1**2*T2 + 0.5*amax*T2**2
            sc[i] = sc_t2_matlab + vmax * (t - t2)

        # MATLAB: else if (t>t3 && t<=t4)
        elif t3 < t <= t4:
            # Calculate values at t3 first
            sc_t2_matlab = (1/6)*je*T1**3 + 0.5*je*T1**2*T2 + 0.5*amax*T2**2
            sc_t3 = sc_t2_matlab + vmax * (T3) # T3 is duration t3-t2

            # Apply formulas for t > t3
            ac[i] = -je * (t - t3) # Start deceleration jerk
            vc[i] = vmax - 0.5 * je * (t - t3)**2
            # MATLAB formula: sc(ntt)=je*T1^3/6+vmax*T2+vmax*T3+vmax*(t-t3)-(1/6)*je*(t-t3)^3; Again, seems wrong.
            # Let's use formula from Image 6:
            # sc(ntt) = sc_t3 + vmax*(t-t3) - (1/6)*je*(t-t3)^3
            sc[i] = sc_t3 + vmax * (t - t3) - (1/6) * je * (t - t3)**3

        # MATLAB: else if (t>t4 && t<=t5)
        elif t4 < t <= t5:
            # Calculate values at t4 first
            sc_t2_matlab = (1/6)*je*T1**3 + 0.5*je*T1**2*T2 + 0.5*amax*T2**2
            sc_t3 = sc_t2_matlab + vmax * (T3)
            sc_t4 = sc_t3 + vmax * (T4) - (1/6) * je * (T4)**3 # T4 is duration t4-t3
            vc_t4 = vmax - 0.5 * je * (T4)**2 # T4 is duration t4-t3

            # Apply formulas for t > t4
            ac[i] = -amax # Constant deceleration?
            vc[i] = vc_t4 - amax * (t - t4)
            # MATLAB formula: sc(ntt)=je*T1^3/6+vmax*T2+vmax*T3+vmax*T4-(1/6)*je*T4^3+vc_t4*(t-t4)-(1/2)*amax*(t-t4)^2;
            # This is sc_t4 + vc_t4*(t-t4) - 0.5*amax*(t-t4)^2 (Standard integral)
            sc[i] = sc_t4 + vc_t4 * (t - t4) - 0.5 * amax * (t - t4)**2
        # Note: The MATLAB code snippet doesn't explicitly show the 6th and 7th segments
        # (final jerk phase to bring velocity to zero). It stops at constant deceleration.
        # This means the final velocity might not be zero.
        # We will replicate this behavior. If a full stop is needed, the profile needs extension.

    # Final checks:
    # Ensure final displacement is close to sall
    if nt > 0 and abs(sc[-1] - sall) > 1e-3: # Tolerance
         print(f"Warning: Final displacement sc[-1] ({sc[-1]:.4f}) differs significantly from target sall ({sall:.4f}).")
         print(f"Profile parameters: {params}")
    # Ensure final velocity is close to zero (or whatever the profile implies)
    # Since the profile seems incomplete (missing final jerk phase), vc[-1] likely won't be 0.
    if nt > 0:
        print(f"Note: Final calculated velocity vc[-1] = {vc[-1]:.4f} (Profile might be incomplete based on MATLAB code)")


    return ac, vc, sc, time_vec, params

# %% --- Main Trajectory Planning and Execution ---

if __name__ == "__main__":

    # --- Get Robot Parameters ---
    robot_params = get_robot_params()
    print("--- Robot Parameters ---")
    for key, value in robot_params.items():
        print(f"{key}: {value}")

    # --- Trajectory Parameters (from MATLAB snippet) ---
    print("\n--- Trajectory Parameters ---")
    r = 200.0  # Radius (mm)
    vmax = 1000.0 # Max velocity (mm/s) (60000/60)
    TA = 0.2   # Acceleration phase time (s) (200 * 0.001)
    ts = 0.01  # Time step (s) (10 * 0.001)
    pi = np.pi # Use numpy pi

    # Circular path parameters
    sall = 2 * pi * r # Total path length (circumference)
    z_level = 1200.0  # Constant Z height for the circle
    y_offset = -190.0 # Center of circle offset in Y

    print(f"Radius (r): {r:.1f} mm")
    print(f"Max Velocity (vmax): {vmax:.1f} mm/s")
    print(f"Acceleration Time (TA): {TA:.3f} s")
    print(f"Time Step (ts): {ts:.3f} s")
    print(f"Total Path Length (sall): {sall:.3f} mm")
    print(f"Z Level: {z_level:.1f} mm")
    print(f"Y Offset: {y_offset:.1f} mm")
    print(f"Circle Center: (0, {y_offset:.1f}, {z_level:.1f})")

    # --- Generate Motion Profile ---
    print("\n--- Generating S-Curve Motion Profile ---")
    try:
        ac, vc, sc, time_vec, profile_params = s_velocity_motion(vmax, TA, sall, ts)
        nt = len(time_vec)
        print(f"Generated {nt} points.")
        print(f"Calculated total time (t5): {profile_params['t5']:.3f} s")
        print(f"Actual time step: {profile_params['actual_ts']:.5f} s")

        # Plot the motion profile
        fig_prof, axs_prof = plt.subplots(3, 1, sharex=True, figsize=(10, 8))
        axs_prof[0].plot(time_vec, sc, label='Displacement (s)')
        axs_prof[0].set_ylabel('s (mm)')
        axs_prof[0].grid(True)
        axs_prof[0].legend()
        axs_prof[0].set_title('S-Curve Motion Profile (Replicating MATLAB)')

        axs_prof[1].plot(time_vec, vc, label='Velocity (v)')
        axs_prof[1].set_ylabel('v (mm/s)')
        axs_prof[1].grid(True)
        axs_prof[1].legend()

        axs_prof[2].plot(time_vec, ac, label='Acceleration (a)')
        axs_prof[2].set_ylabel('a (mm/s^2)')
        axs_prof[2].set_xlabel('Time (s)')
        axs_prof[2].grid(True)
        axs_prof[2].legend()

        fig_prof.tight_layout()
        plt.show(block=False)

    except ValueError as e:
        print(f"Error generating motion profile: {e}")
        exit()
    except Exception as e:
        print(f"An unexpected error occurred during profile generation: {e}")
        exit()


    # --- Calculate Trajectory Points and Joint Values ---
    print("\n--- Calculating Inverse Kinematics for Trajectory ---")

    # Initialize arrays to store results
    xp = np.zeros(nt)
    yp = np.zeros(nt)
    zp = np.zeros(nt)
    vxp = np.zeros(nt)
    vyp = np.zeros(nt)
    vzp = np.zeros(nt)
    axp = np.zeros(nt)
    ayp = np.zeros(nt)
    azp = np.zeros(nt)

    q1 = np.zeros(nt)
    q2 = np.zeros(nt)
    q3 = np.zeros(nt)
    q4 = np.zeros(nt)
    theta1 = np.zeros(nt)
    theta2 = np.zeros(nt)

    # Store R4 and rp for animation/plotting
    R4_traj = [np.eye(3)] * nt
    rp_traj = [np.zeros((3,1))] * nt

    ik_errors = 0
    start_time_ik = time.time()

    for ntt in range(nt):
        # Current displacement along the path
        s_current = sc[ntt]
        # Current angle on the circle (starts from -pi/2 direction based on MATLAB)
        angle = s_current / r - (pi / 2) # Adjust start angle if needed

        # Calculate target Cartesian position (matches MATLAB formulas)
        xp[ntt] = r * np.cos(angle) # MATLAB uses sin(sc/r), which corresponds to cos(angle) here
        yp[ntt] = r * np.sin(angle) + y_offset # MATLAB uses -r*cos(sc/r)-190 -> r*sin(angle)+y_offset
        zp[ntt] = z_level

        # Calculate target Cartesian velocity (using chain rule: d(pos)/dt = d(pos)/d(angle) * d(angle)/dt)
        # d(angle)/dt = (d(angle)/ds) * (ds/dt) = (1/r) * vc[ntt]
        v_current = vc[ntt]
        vxp[ntt] = -r * np.sin(angle) * (v_current / r) # = -v_current * sin(angle)
        vyp[ntt] =  r * np.cos(angle) * (v_current / r) # =  v_current * cos(angle)
        vzp[ntt] = 0
        # Check against MATLAB formulas:
        # vxp = vc * cos(sc/r) -> vc * cos(angle+pi/2) = -vc * sin(angle). Matches.
        # vyp = vc * sin(sc/r) -> vc * sin(angle+pi/2) = vc * cos(angle). Matches.

        # Calculate target Cartesian acceleration (tangential + centripetal)
        # Tangential acceleration vector: a_t = ac[ntt] * unit_tangent
        # Centripetal acceleration vector: a_c = (vc[ntt]^2 / r) * unit_normal_inward
        a_current = ac[ntt]
        unit_tangent = np.array([-np.sin(angle), np.cos(angle), 0])
        unit_normal_inward = np.array([-np.cos(angle), -np.sin(angle), 0]) # Points towards circle center (relative)

        accel_tangential = a_current * unit_tangent
        accel_centripetal = (v_current**2 / r) * unit_normal_inward if r > 1e-6 else np.zeros(3)

        accel_total = accel_tangential + accel_centripetal
        axp[ntt], ayp[ntt], azp[ntt] = accel_total[0], accel_total[1], accel_total[2]
        # Check against MATLAB formulas:
        # axp = ac*cos(sc/r) - (vc^2/r)*sin(sc/r) -> ac*(-sin(angle)) - (vc^2/r)*cos(angle). Matches tangential_x + centripetal_x.
        # ayp = ac*sin(sc/r) + (vc^2/r)*cos(sc/r) -> ac*(cos(angle)) + (vc^2/r)*(-sin(angle)). Matches tangential_y + centripetal_y.

        # Perform Inverse Kinematics
        ik_result = inverse_kinematics(xp[ntt], yp[ntt], zp[ntt], robot_params)

        if ik_result['error']:
            print(f"IK Error at step {ntt}/{nt} (t={time_vec[ntt]:.3f}s): {ik_result['error']}")
            ik_errors += 1
            # Handle error: Use previous joint values? Stop? For now, store NaN or zeros.
            q1[ntt], q2[ntt], q3[ntt], q4[ntt] = np.nan, np.nan, np.nan, np.nan
            theta1[ntt], theta2[ntt] = np.nan, np.nan
            # Keep previous R4 and rp for visualization continuity
            if ntt > 0:
                R4_traj[ntt] = R4_traj[ntt-1]
                rp_traj[ntt] = rp_traj[ntt-1]
            else:
                R4_traj[ntt] = np.eye(3)
                rp_traj[ntt] = np.array([[0],[0],[robot_params['e']]]) # Default start pose
        else:
            q1[ntt] = ik_result['q1']
            q2[ntt] = ik_result['q2']
            q3[ntt] = ik_result['q3']
            q4[ntt] = ik_result['q4']
            theta1[ntt] = ik_result['theta1']
            theta2[ntt] = ik_result['theta2']
            R4_traj[ntt] = ik_result['R4']
            rp_traj[ntt] = ik_result['rp']

        # Optional: Print progress
        if (ntt + 1) % 100 == 0 or ntt == nt - 1:
             print(f"  Processed point {ntt+1}/{nt}")

    end_time_ik = time.time()
    print(f"IK calculation finished in {end_time_ik - start_time_ik:.2f} seconds.")
    if ik_errors > 0:
        print(f"Warning: Encountered {ik_errors} IK errors.")

    # --- Calculate Joint Velocities (using Finite Differences) ---
    print("\n--- Calculating Joint Velocities (Finite Differences) ---")
    vq1 = np.zeros(nt)
    vq2 = np.zeros(nt)
    vq3 = np.zeros(nt)

    # Forward difference for the first point, backward for the last, central for others
    if nt > 1:
        # Central difference for interior points
        vq1[1:-1] = (q1[2:] - q1[:-2]) / (2 * profile_params['actual_ts'])
        vq2[1:-1] = (q2[2:] - q2[:-2]) / (2 * profile_params['actual_ts'])
        vq3[1:-1] = (q3[2:] - q3[:-2]) / (2 * profile_params['actual_ts'])
        # Forward difference for the first point
        vq1[0] = (q1[1] - q1[0]) / profile_params['actual_ts']
        vq2[0] = (q2[1] - q2[0]) / profile_params['actual_ts']
        vq3[0] = (q3[1] - q3[0]) / profile_params['actual_ts']
        # Backward difference for the last point
        vq1[-1] = (q1[-1] - q1[-2]) / profile_params['actual_ts']
        vq2[-1] = (q2[-1] - q2[-2]) / profile_params['actual_ts']
        vq3[-1] = (q3[-1] - q3[-2]) / profile_params['actual_ts']
    elif nt == 1:
        # Only one point, velocity is zero
        pass # Already initialized to zero

    # Handle NaNs from IK errors if any
    vq1[np.isnan(q1)] = np.nan
    vq2[np.isnan(q2)] = np.nan
    vq3[np.isnan(q3)] = np.nan
    print("Joint velocities calculated.")


    # --- Plotting Results ---
    print("\n--- Plotting Trajectory Results ---")

    # 1. Cartesian Path Plot
    fig_path = plt.figure(figsize=(8, 8))
    ax_path = fig_path.add_subplot(111, projection='3d')
    ax_path.plot(xp, yp, zp, 'b-', label='Desired Path')
    ax_path.scatter(xp[0], yp[0], zp[0], c='g', marker='o', s=100, label='Start Point')
    ax_path.scatter(xp[-1], yp[-1], zp[-1], c='r', marker='x', s=100, label='End Point')
    ax_path.set_xlabel('X (mm)')
    ax_path.set_ylabel('Y (mm)')
    ax_path.set_zlabel('Z (mm)')
    ax_path.set_title('Desired Cartesian Trajectory (Circular Path)')
    # Make axes equal for better visualization of the circle
    try:
        ax_path.set_aspect('equal', adjustable='box')
    except NotImplementedError:
         # Fallback if equal aspect ratio fails
        max_range = np.max([xp.max()-xp.min(), yp.max()-yp.min(), zp.max()-zp.min()])
        mid_x = (xp.max()+xp.min())/2
        mid_y = (yp.max()+yp.min())/2
        mid_z = (zp.max()+zp.min())/2
        ax_path.set_xlim(mid_x - max_range/2, mid_x + max_range/2)
        ax_path.set_ylim(mid_y - max_range/2, mid_y + max_range/2)
        ax_path.set_zlim(mid_z - max_range/2, mid_z + max_range/2)

    ax_path.legend()
    ax_path.grid(True)
    plt.show(block=False)


    # 2. Joint Space Plots
    fig_joint, axs_joint = plt.subplots(2, 1, sharex=True, figsize=(10, 8))

    # Joint Positions
    axs_joint[0].plot(time_vec, q1, label='q1')
    axs_joint[0].plot(time_vec, q2, label='q2')
    axs_joint[0].plot(time_vec, q3, label='q3')
    # axs_joint[0].plot(time_vec, q4, label='q4 (variable length)') # Optional
    axs_joint[0].set_ylabel('Joint Position (mm or rad)')
    axs_joint[0].set_title('Joint Space Trajectory')
    axs_joint[0].grid(True)
    axs_joint[0].legend()

    # Joint Velocities
    axs_joint[1].plot(time_vec, vq1, label='vq1')
    axs_joint[1].plot(time_vec, vq2, label='vq2')
    axs_joint[1].plot(time_vec, vq3, label='vq3')
    axs_joint[1].set_ylabel('Joint Velocity (mm/s or rad/s)')
    axs_joint[1].set_xlabel('Time (s)')
    axs_joint[1].grid(True)
    axs_joint[1].legend()

    fig_joint.tight_layout()
    plt.show(block=False)


    # --- Animation (Optional) ---
    print("\n--- Preparing Animation (Close other plots to start) ---")
    input("Press Enter to start animation...")

    fig_anim = plt.figure(figsize=(9, 9))
    ax_anim = fig_anim.add_subplot(111, projection='3d')
    plt.ion() # Turn on interactive mode

    try:
        # Determine plot limits based on the whole trajectory
        all_rp = np.array([r.flatten() for r in rp_traj if r is not None and not np.isnan(r).any()])
        if len(all_rp) > 0:
            min_coords = np.min(all_rp, axis=0) - 200 # Add padding
            max_coords = np.max(all_rp, axis=0) + 200 # Add padding
        else: # Default limits if no valid points
            min_coords = np.array([-500, -500, robot_params['e'] - 100])
            max_coords = np.array([500, 500, robot_params['e'] + 1000])

        # Fixed plot limits for animation stability
        ax_anim.set_xlim(min(min_coords[0], -robot_params['b22']-100), max(max_coords[0], robot_params['b22']+100))
        ax_anim.set_ylim(min(min_coords[1], -robot_params['b11']-100), max(max_coords[1], 100)) # Y often negative
        ax_anim.set_zlim(0, max(max_coords[2], robot_params['e'] + 500))


        for i in range(nt):
            if np.isnan(q1[i]): # Skip frames with IK errors
                 print(f"Skipping frame {i} due to IK error.")
                 continue

            rp_current = rp_traj[i]
            R4_current = R4_traj[i]

            plot_robot(ax_anim, rp_current, R4_current, robot_params,
                       title=f"Robot Trajectory Animation (t={time_vec[i]:.2f}s)")

            # Draw trajectory path up to current point
            ax_anim.plot(xp[:i+1], yp[:i+1], zp[:i+1], 'm:', lw=1, label='_Trajectory Path') # Underscore hides label repetition

            plt.draw()
            # Adjust pause time for desired animation speed
            # Faster than real-time: pause_time = 0.01
            # Closer to real-time: pause_time = max(0.001, ts - elapsed_draw_time)
            plt.pause(0.01) # Small pause to allow rendering

            # Stop animation if figure is closed
            if not plt.fignum_exists(fig_anim.number):
                print("Animation window closed by user.")
                break

    except KeyboardInterrupt:
        print("\nAnimation interrupted by user.")
    finally:
        plt.ioff() # Turn off interactive mode
        print("\nAnimation finished. Close the final plot window to exit.")
        plt.show() # Keep the last frame visible