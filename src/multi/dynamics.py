import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root
from matplotlib.animation import FuncAnimation # Import FuncAnimation

# --- 参数定义 ---
L1 = 0.3  # 杆1 (AB) 长度 (m)
L2 = 0.6  # 杆2 (BC) 长度 (m)
L3 = 0.4  # 杆3 (CD) 长度 (m)
L4 = 0.5  # 机架 (AD) 长度 (m)

# 假设质心位于杆的中点
# COM (Center of Mass) positions relative to the start of the link
r_C1_local = L1 / 2
r_C2_local = L2 / 2
r_C3_local = L3 / 2

# --- 输入运动 ---
omega1 = np.pi / 2  # 杆1的角速度 (rad/s) - 恒定
# theta1_initial = np.pi / 4 # Original - Likely causes singularity at t=0
theta1_initial = np.pi / 4 + 0.01 # Slightly offset initial angle to avoid starting exactly on singularity
T_total = 4.0       # 总仿真时间 (s)
dt = 0.01           # 时间步长 (s)
time_steps = int(T_total / dt)
t = np.linspace(0, T_total, time_steps)

# --- 存储结果的数组 ---
theta1_hist = np.zeros(time_steps)
theta2_hist = np.zeros(time_steps)
theta3_hist = np.zeros(time_steps)

omega1_hist = np.zeros(time_steps)
omega2_hist = np.zeros(time_steps)
omega3_hist = np.zeros(time_steps)

alpha1_hist = np.zeros(time_steps)
alpha2_hist = np.zeros(time_steps)
alpha3_hist = np.zeros(time_steps)

# 质心运动学 (COM Kinematics)
# L1 COM (C1)
r_C1_hist = np.zeros((time_steps, 2))
v_C1_hist = np.zeros((time_steps, 2))
a_C1_hist = np.zeros((time_steps, 2))
# L2 COM (C2)
r_C2_hist = np.zeros((time_steps, 2))
v_C2_hist = np.zeros((time_steps, 2))
a_C2_hist = np.zeros((time_steps, 2))
# L3 COM (C3)
r_C3_hist = np.zeros((time_steps, 2))
v_C3_hist = np.zeros((time_steps, 2))
a_C3_hist = np.zeros((time_steps, 2))

# --- 约束方程函数 ---
def constraint_equations(vars, theta1_val):
    """ 位置约束方程 Phi(q) = 0 """
    theta2, theta3 = vars
    eq1 = L1 * np.cos(theta1_val) + L2 * np.cos(theta2) - L3 * np.cos(theta3) - L4
    eq2 = L1 * np.sin(theta1_val) + L2 * np.sin(theta2) - L3 * np.sin(theta3)
    return [eq1, eq2]

# --- 雅可比矩阵函数 ---
def jacobian(theta1, theta2, theta3):
    """ 约束雅可比矩阵 Phi_q """
    # Note: The Jacobian should be for the constraint equations C(q) = 0
    # C1 = L1 cos(th1) + L2 cos(th2) - L3 cos(th3) - L4 = 0
    # C2 = L1 sin(th1) + L2 sin(th2) - L3 sin(th3)      = 0
    # Phi_q = dC/dq where q = [th1, th2, th3]
    J = np.array([
        [-L1 * np.sin(theta1), -L2 * np.sin(theta2), L3 * np.sin(theta3)],
        [ L1 * np.cos(theta1),  L2 * np.cos(theta2), -L3 * np.cos(theta3)]
    ])
    return J

# --- 二次速度项 gamma ---
# gamma = - (d(Phi_q)/dt * dot_q) = - (d(Phi_q)/dq * dot_q) * dot_q
# This is equivalent to - (Phi_qq * dot_q) * dot_q
# Or more simply, differentiating the velocity constraint Phi_q * dot_q = 0 w.r.t time:
# d(Phi_q)/dt * dot_q + Phi_q * ddot_q = 0
# So Phi_q * ddot_q = - d(Phi_q)/dt * dot_q = gamma
def gamma_vector(theta1, theta2, theta3, omega1, omega2, omega3):
    """ 计算 gamma = -d(Phi_q)/dt * dq/dt """
    # d(Phi_q)/dt * dot_q = [d/dt(-L1s1) d/dt(-L2s2) d/dt(L3s3)] * [w1]
    #                      [d/dt( L1c1) d/dt( L2c2) d/dt(-L3c3)] * [w2]
    #                                                             [w3]
    # d/dt(-L1s1) = -L1c1*w1; d/dt(-L2s2) = -L2c2*w2; d/dt(L3s3) = L3c3*w3
    # d/dt( L1c1) = -L1s1*w1; d/dt( L2c2) = -L2s2*w2; d/dt(-L3c3) = L3s3*w3
    # Row1 = (-L1c1*w1)*w1 + (-L2c2*w2)*w2 + (L3c3*w3)*w3
    # Row2 = (-L1s1*w1)*w1 + (-L2s2*w2)*w2 + (L3s3*w3)*w3
    # gamma = - result
    gamma = np.array([
        L1 * np.cos(theta1) * omega1**2 + L2 * np.cos(theta2) * omega2**2 - L3 * np.cos(theta3) * omega3**2,
        L1 * np.sin(theta1) * omega1**2 + L2 * np.sin(theta2) * omega2**2 - L3 * np.sin(theta3) * omega3**2
    ])
    return gamma

# --- 初始猜测 ---
# Provide a reasonable starting guess, e.g., based on visual inspection or expected range
initial_guess = [np.pi/2, np.pi/3] # Initial guess [theta2, theta3]
# Try to find a better initial guess for the very first step using the solver
print("Solving for initial position...")
sol_init = root(constraint_equations, initial_guess, args=(theta1_initial,), method='lm', tol=1e-9)
if sol_init.success:
    initial_guess = sol_init.x
    print(f"Initial position found: theta2={np.rad2deg(initial_guess[0]):.2f} deg, theta3={np.rad2deg(initial_guess[1]):.2f} deg")
else:
    print(f"Warning: Initial position solver failed. Message: {sol_init.message}")
    print("Using default initial guess, simulation might fail or be inaccurate at the start.")
    # Check if linkage can physically reach this theta1
    dist_ad_sq = L4**2 # Should be L4^2, but we need distance AC
    # Point A=(0,0), D=(L4,0)
    # Point B=(L1*c1, L1*s1)
    # Point C needs to be found. Distance AC^2 = (Cx-Ax)^2 + (Cy-Ay)^2
    # Distance DC^2 = (Cx-Dx)^2 + (Cy-Dy)^2 = L3^2
    # Distance BC^2 = (Cx-Bx)^2 + (Cy-By)^2 = L2^2
    # Check distance AD vs possible range L2+L3 and |L2-L3|
    dist_ac_sq = (L4 - L1*np.cos(theta1_initial))**2 + (-L1*np.sin(theta1_initial))**2
    print(f"  Checking assembly for theta1={np.rad2deg(theta1_initial):.2f} deg:")
    print(f"  Distance AC^2 = {dist_ac_sq:.4f}")
    print(f"  (L2+L3)^2 = {(L2+L3)**2:.4f}")
    print(f"  (L2-L3)^2 = {(L2-L3)**2:.4f}")
    if dist_ac_sq > (L2+L3)**2 + 1e-6 or dist_ac_sq < (L2-L3)**2 - 1e-6: # Add tolerance
         print("  Reason: Linkage cannot assemble (Grashof condition violation or limit reached).")
         # Exit if initial position cannot be found
         exit()


# --- 数值仿真循环 ---
print("Starting simulation...")
simulation_successful = True
for i in range(time_steps):
    current_t = t[i]
    # 1. 计算当前输入运动
    theta1 = theta1_initial + omega1 * current_t
    dot_theta1 = omega1 # 恒定角速度
    ddot_theta1 = 0     # 恒定角速度，角加速度为0

    # Ensure theta1 is within [-pi, pi] or [0, 2pi] if needed, though not strictly necessary for cos/sin
    theta1 = np.arctan2(np.sin(theta1), np.cos(theta1))
    theta1_hist[i] = theta1
    omega1_hist[i] = dot_theta1
    alpha1_hist[i] = ddot_theta1

    # 2. 求解位置约束 (theta2, theta3)
    # Use the previous step's result as the initial guess for the current step
    if i > 0:
        # Check if previous step was valid
        if not np.isnan(theta2_hist[i-1]) and not np.isnan(theta3_hist[i-1]):
            initial_guess = [theta2_hist[i-1], theta3_hist[i-1]]
        # else: keep the last known good guess or the default initial_guess
        # (Using the initial_guess from the start might be safer if NaNs appear)

    sol = root(constraint_equations, initial_guess, args=(theta1,), method='lm', tol=1e-9) # Levenberg-Marquardt
    if not sol.success:
        print(f"FATAL: Position solver failed at step {i}, t={current_t:.3f}. Message: {sol.message}")
        print(f"  theta1={np.rad2deg(theta1):.2f} deg, guess=[{np.rad2deg(initial_guess[0]):.2f}, {np.rad2deg(initial_guess[1]):.2f}] deg")
        # Check if linkage can physically reach this theta1
        dist_ac_sq = (L4 - L1*np.cos(theta1))**2 + (-L1*np.sin(theta1))**2
        if dist_ac_sq > (L2+L3)**2 + 1e-6 or dist_ac_sq < (L2-L3)**2 - 1e-6: # Add tolerance
             print("  Reason: Linkage cannot assemble (limit reached).")
        else:
             print("  Reason: Solver failed to converge (check tolerances, initial guess, or proximity to singularity).")

        simulation_successful = False
        # Fill remaining steps with NaN
        theta2_hist[i:] = np.nan
        theta3_hist[i:] = np.nan
        omega2_hist[i:] = np.nan
        omega3_hist[i:] = np.nan
        alpha2_hist[i:] = np.nan
        alpha3_hist[i:] = np.nan
        r_C1_hist[i:,:] = np.nan
        v_C1_hist[i:,:] = np.nan
        a_C1_hist[i:,:] = np.nan
        r_C2_hist[i:,:] = np.nan
        v_C2_hist[i:,:] = np.nan
        a_C2_hist[i:,:] = np.nan
        r_C3_hist[i:,:] = np.nan
        v_C3_hist[i:,:] = np.nan
        a_C3_hist[i:,:] = np.nan
        break # Stop simulation

    theta2, theta3 = sol.x
    # Normalize angles to prevent excessive wrapping
    theta2_hist[i] = np.arctan2(np.sin(theta2), np.cos(theta2))
    theta3_hist[i] = np.arctan2(np.sin(theta3), np.cos(theta3))
    theta2 = theta2_hist[i] # use normalized value for consistency
    theta3 = theta3_hist[i]

    # 3. 求解速度约束 (omega2, omega3)
    # Phi_q * dot_q = 0  =>  Phi_q[:, 1:] * [dot_theta2; dot_theta3] = - Phi_q[:, 0] * dot_theta1
    J = jacobian(theta1, theta2, theta3)
    A_vel = J[:, 1:3] # Columns corresponding to theta2, theta3
    b_vel = -J[:, 0] * dot_theta1 # Column corresponding to theta1

    dot_theta2 = np.nan # Initialize as NaN
    dot_theta3 = np.nan

    # --- Singularity Check ---
    # Determinant = L2*L3*sin(theta2 - theta3)
    det_A_vel = np.linalg.det(A_vel)
    # Use condition number as a more robust check for near-singularity
    # cond_A_vel = np.linalg.cond(A_vel)

    # Define a threshold for singularity (can be tuned)
    singularity_threshold = 1e-6

    if abs(det_A_vel) < singularity_threshold:
         print(f"Warning: Near singularity detected (det={det_A_vel:.2e}) at step {i}, t={current_t:.3f}")
         print(f"  theta1={np.rad2deg(theta1):.2f}, theta2={np.rad2deg(theta2):.2f}, theta3={np.rad2deg(theta3):.2f}")
         print(f"  theta2 - theta3 = {np.rad2deg(theta2 - theta3):.2f} deg")
         # Keep velocities as NaN, cannot reliably solve
    else:
        try:
            dot_theta23 = np.linalg.solve(A_vel, b_vel)
            dot_theta2 = dot_theta23[0]
            dot_theta3 = dot_theta23[1]
        except np.linalg.LinAlgError as e:
            # This might still happen if matrix is ill-conditioned even if det isn't exactly zero
            print(f"Warning: Velocity solver failed ({e}) at step {i}, t={current_t:.3f}")
            # Keep velocities as NaN

    omega2_hist[i] = dot_theta2
    omega3_hist[i] = dot_theta3

    # 4. 求解加速度约束 (alpha2, alpha3)
    # Phi_q * ddot_q = gamma => Phi_q[:, 1:] * [ddot_theta2; ddot_theta3] = gamma - Phi_q[:, 0] * ddot_theta1
    ddot_theta2 = np.nan # Initialize as NaN
    ddot_theta3 = np.nan

    # Cannot calculate acceleration if velocity is unknown or at singularity
    if np.isnan(dot_theta2) or np.isnan(dot_theta3) or abs(det_A_vel) < singularity_threshold:
        # Keep accelerations as NaN
        pass # Already initialized to NaN
    else:
        gamma = gamma_vector(theta1, theta2, theta3, dot_theta1, dot_theta2, dot_theta3)
        A_acc = A_vel # Matrix is the same as for velocity
        b_acc = gamma - J[:, 0] * ddot_theta1

        try:
            # We already checked the determinant/singularity for velocity
            ddot_theta23 = np.linalg.solve(A_acc, b_acc)
            ddot_theta2 = ddot_theta23[0]
            ddot_theta3 = ddot_theta23[1]
        except np.linalg.LinAlgError as e:
            # Should ideally not happen if velocity solve succeeded, but check anyway
            print(f"Warning: Acceleration solver failed ({e}) at step {i}, t={current_t:.3f}")
            # Keep accelerations as NaN

    alpha2_hist[i] = ddot_theta2
    alpha3_hist[i] = ddot_theta3

    # 5. 计算质心运动学 (COM Kinematics)
    # Check if angles/velocities/accelerations are valid before calculating
    # Note: We allow calculation even if acc is NaN, but vel/pos must be valid
    if np.isnan(theta1) or np.isnan(theta2) or np.isnan(theta3) or \
       np.isnan(dot_theta1) or np.isnan(dot_theta2) or np.isnan(dot_theta3):
         r_C1_hist[i, :] = [np.nan, np.nan]
         v_C1_hist[i, :] = [np.nan, np.nan]
         a_C1_hist[i, :] = [np.nan, np.nan]
         r_C2_hist[i, :] = [np.nan, np.nan]
         v_C2_hist[i, :] = [np.nan, np.nan]
         a_C2_hist[i, :] = [np.nan, np.nan]
         r_C3_hist[i, :] = [np.nan, np.nan]
         v_C3_hist[i, :] = [np.nan, np.nan]
         a_C3_hist[i, :] = [np.nan, np.nan]
    else:
        # --- 杆 1 质心 (C1) ---
        r_C1_hist[i, :] = [r_C1_local * np.cos(theta1), r_C1_local * np.sin(theta1)]
        v_C1_hist[i, :] = [-r_C1_local * np.sin(theta1) * dot_theta1, r_C1_local * np.cos(theta1) * dot_theta1]
        # Check if acceleration is valid before using it
        if not np.isnan(ddot_theta1):
             a_C1_hist[i, :] = [-r_C1_local * (np.cos(theta1) * dot_theta1**2 + np.sin(theta1) * ddot_theta1),
                                r_C1_local * (-np.sin(theta1) * dot_theta1**2 + np.cos(theta1) * ddot_theta1)]
        else:
             a_C1_hist[i, :] = [np.nan, np.nan]


        # --- 杆 2 质心 (C2) ---
        r_B = np.array([L1 * np.cos(theta1), L1 * np.sin(theta1)])
        r_C2_rel_B = np.array([r_C2_local * np.cos(theta2), r_C2_local * np.sin(theta2)])
        r_C2_hist[i, :] = r_B + r_C2_rel_B

        v_B = np.array([-L1 * np.sin(theta1) * dot_theta1, L1 * np.cos(theta1) * dot_theta1])
        v_C2_rel_B = np.array([-r_C2_local * np.sin(theta2) * dot_theta2, r_C2_local * np.cos(theta2) * dot_theta2])
        v_C2_hist[i, :] = v_B + v_C2_rel_B

        # Check if accelerations are valid before using them
        if not np.isnan(ddot_theta1) and not np.isnan(ddot_theta2):
            a_B = np.array([-L1 * (np.cos(theta1) * dot_theta1**2 + np.sin(theta1) * ddot_theta1),
                            L1 * (-np.sin(theta1) * dot_theta1**2 + np.cos(theta1) * ddot_theta1)])
            a_C2_rel_B = np.array([-r_C2_local * (np.cos(theta2) * dot_theta2**2 + np.sin(theta2) * ddot_theta2),
                                   r_C2_local * (-np.sin(theta2) * dot_theta2**2 + np.cos(theta2) * ddot_theta2)])
            a_C2_hist[i, :] = a_B + a_C2_rel_B
        else:
            a_C2_hist[i, :] = [np.nan, np.nan]


        # --- 杆 3 质心 (C3) ---
        # Point D is at (L4, 0)
        # r_C3 is position of COM of Link 3 relative to Origin A
        # r_C3 = r_D + r_C3_rel_D
        # r_C3_rel_D is vector from D to C3 (COM of link 3)
        # Vector DC is L3 long at angle theta3. Vector DC3 is r_C3_local long at angle theta3.
        # Position of C = (L4 + L3*cos(theta3), L3*sin(theta3))
        # Position of C3 = (L4 + r_C3_local*cos(theta3), r_C3_local*sin(theta3)) -> This seems incorrect.
        # C3 is r_C3_local along the link CD *from C towards D*.
        # Let's define r_C3 relative to D, pointing towards C.
        # Vector DC = [L3*cos(theta3), L3*sin(theta3)]
        # Vector C3 relative to D = (L3 - r_C3_local) * unit_vector_DC ? No, COM is usually from start of link.
        # Let's assume r_C3_local is distance from D along link DC.
        r_D = np.array([L4, 0])
        r_C3_rel_D = np.array([r_C3_local * np.cos(theta3), r_C3_local * np.sin(theta3)])
        # This assumes theta3 is angle of vector DC relative to horizontal *originating at D*.
        # Let's verify this from constraint: L1c1+L2c2 = L4+L3c3, L1s1+L2s2=L3s3. Yes, theta3 is angle of DC.
        r_C3_hist[i, :] = r_D + r_C3_rel_D # Position of C3 relative to global origin A

        # Velocity of D is zero.
        v_C3_rel_D = np.array([-r_C3_local * np.sin(theta3) * dot_theta3, r_C3_local * np.cos(theta3) * dot_theta3])
        v_C3_hist[i, :] = v_C3_rel_D # v_D = 0

        # Acceleration of D is zero. Check if acceleration is valid before using it
        if not np.isnan(ddot_theta3):
             a_C3_rel_D = np.array([-r_C3_local * (np.cos(theta3) * dot_theta3**2 + np.sin(theta3) * ddot_theta3),
                                    r_C3_local * (-np.sin(theta3) * dot_theta3**2 + np.cos(theta3) * ddot_theta3)])
             a_C3_hist[i, :] = a_C3_rel_D # a_D = 0
        else:
             a_C3_hist[i, :] = [np.nan, np.nan]


if simulation_successful:
    print("Simulation finished successfully.")
else:
    print("Simulation failed or encountered singularities. Plots may be incomplete or show NaN values.")

# --- 绘图 ---
# (Plotting code remains the same as provided, it should handle NaNs gracefully)
print("Generating plots...")

# 图2: 角度、角速度、角加速度 (Theta1, Theta2, Theta3)
fig1, axs1 = plt.subplots(3, 1, figsize=(10, 12), sharex=True)
fig1.suptitle('Link Angles, Angular Velocities, and Angular Accelerations')

# 角度
axs1[0].plot(t, np.rad2deg(theta1_hist), label=r'$\theta_1$')
axs1[0].plot(t, np.rad2deg(theta2_hist), label=r'$\theta_2$')
axs1[0].plot(t, np.rad2deg(theta3_hist), label=r'$\theta_3$')
axs1[0].set_ylabel('Angle (deg)')
axs1[0].legend()
axs1[0].grid(True)

# 角速度
axs1[1].plot(t, omega1_hist, label=r'$\omega_1$')
axs1[1].plot(t, omega2_hist, label=r'$\omega_2$')
axs1[1].plot(t, omega3_hist, label=r'$\omega_3$')
axs1[1].set_ylabel('Angular Velocity (rad/s)')
axs1[1].legend()
axs1[1].grid(True)

# 角加速度
axs1[2].plot(t, alpha1_hist, label=r'$\alpha_1$')
axs1[2].plot(t, alpha2_hist, label=r'$\alpha_2$')
axs1[2].plot(t, alpha3_hist, label=r'$\alpha_3$')
axs1[2].set_ylabel('Angular Acceleration (rad/s$^2$)')
axs1[2].set_xlabel('Time (s)')
axs1[2].legend()
axs1[2].grid(True)

plt.tight_layout(rect=[0, 0.03, 1, 0.96]) # Adjust layout


# 图3 & 图4: 各杆质心运动学
link_names = ['L1 (AB)', 'L2 (BC)', 'L3 (CD)']
com_r_data = [r_C1_hist, r_C2_hist, r_C3_hist]
com_v_data = [v_C1_hist, v_C2_hist, v_C3_hist]
com_a_data = [a_C1_hist, a_C2_hist, a_C3_hist]

for i_link in range(3): # Loop through links 1, 2, 3 (use different index name)
    fig, axs = plt.subplots(3, 2, figsize=(12, 10), sharex='col') # Share x-axis per column
    fig.suptitle(f'Center of Mass Kinematics for Link {i_link+1} ({link_names[i_link]})')

    r_com = com_r_data[i_link]
    v_com = com_v_data[i_link]
    a_com = com_a_data[i_link]

    # 位移 (Position)
    axs[0, 0].plot(t, r_com[:, 0], label='x')
    axs[0, 0].set_ylabel('Position x (m)')
    axs[0, 0].grid(True)
    axs[0, 0].legend()
    axs[0, 1].plot(t, r_com[:, 1], label='y')
    axs[0, 1].set_ylabel('Position y (m)')
    axs[0, 1].grid(True)
    axs[0, 1].legend()

    # 速度 (Velocity)
    axs[1, 0].plot(t, v_com[:, 0], label='vx')
    axs[1, 0].set_ylabel('Velocity x (m/s)')
    axs[1, 0].grid(True)
    axs[1, 0].legend()
    axs[1, 1].plot(t, v_com[:, 1], label='vy')
    axs[1, 1].set_ylabel('Velocity y (m/s)')
    axs[1, 1].grid(True)
    axs[1, 1].legend()

    # 加速度 (Acceleration)
    axs[2, 0].plot(t, a_com[:, 0], label='ax')
    axs[2, 0].set_ylabel('Acceleration x (m/s$^2$)')
    axs[2, 0].set_xlabel('Time (s)')
    axs[2, 0].grid(True)
    axs[2, 0].legend()
    axs[2, 1].plot(t, a_com[:, 1], label='ay')
    axs[2, 1].set_ylabel('Acceleration y (m/s$^2$)')
    axs[2, 1].set_xlabel('Time (s)')
    axs[2, 1].grid(True)
    axs[2, 1].legend()

    plt.tight_layout(rect=[0, 0.03, 1, 0.96]) # Adjust layout

# --- Visualization ---
# Only run animation if simulation was successful (or mostly successful)
# Check if there are any valid (non-NaN) results to animate
if np.any(~np.isnan(theta2_hist)):
    print("Setting up animation...")

    fig_anim, ax_anim = plt.subplots(figsize=(8, 7))
    ax_anim.set_aspect('equal')
    # Determine axis limits based on max reach, add some padding
    # Max X extent could be L4+L3 or L1+L2 depending on configuration
    # Max Y extent could be L1+L2 or L3
    max_r = L1 + L2 + L3 + L4 # Overestimate
    ax_anim.set_xlim([min(-L1, L4-L3) - 0.2, max(L1+L2, L4+L3) + 0.2])
    ax_anim.set_ylim([-max(L1, L3)*1.1 - 0.2, max(L1+L2, L3)*1.1 + 0.2]) # Adjusted Y limits


    ax_anim.set_xlabel("X (m)")
    ax_anim.set_ylabel("Y (m)")
    ax_anim.set_title("Four-Bar Linkage Animation")
    ax_anim.grid(True)

    # Plot the fixed link AD
    ax_anim.plot([0, L4], [0, 0], 'k-', lw=4, label='Link 4 (AD)', zorder=1) # Ground link

    # Initialize plot elements that will be updated
    link1, = ax_anim.plot([], [], 'r-', lw=3, label='Link 1 (AB)', zorder=2)
    link2, = ax_anim.plot([], [], 'g-', lw=3, label='Link 2 (BC)', zorder=2)
    link3, = ax_anim.plot([], [], 'b-', lw=3, label='Link 3 (CD)', zorder=2)
    joint_A, = ax_anim.plot(0, 0, 'ko', ms=6, zorder=3)
    joint_B, = ax_anim.plot([], [], 'ko', ms=6, zorder=3)
    joint_C, = ax_anim.plot([], [], 'ko', ms=6, zorder=3)
    joint_D, = ax_anim.plot(L4, 0, 'ko', ms=6, zorder=3)
    # Optional: COM markers
    com1, = ax_anim.plot([], [], 'rx', ms=5, mew=2, label='COM 1', zorder=4)
    com2, = ax_anim.plot([], [], 'gx', ms=5, mew=2, label='COM 2', zorder=4)
    com3, = ax_anim.plot([], [], 'bx', ms=5, mew=2, label='COM 3', zorder=4)

    time_template = 'Time = %.2fs'
    time_text = ax_anim.text(0.05, 0.95, '', transform=ax_anim.transAxes, fontsize=12,
                             bbox=dict(boxstyle='round,pad=0.3', fc='wheat', alpha=0.5))

    ax_anim.legend(loc='upper right')

    # Store artists that need updating
    artists = [link1, link2, link3, joint_B, joint_C, com1, com2, com3, time_text]

    def animate(i_frame):
        """Update the plot for frame i_frame."""
        # Check for NaN values from failed simulation steps
        if np.isnan(theta1_hist[i_frame]) or np.isnan(theta2_hist[i_frame]) or np.isnan(theta3_hist[i_frame]):
            # If NaN, make the moving parts invisible for this frame
            link1.set_data([], [])
            link2.set_data([], [])
            link3.set_data([], [])
            joint_B.set_data([], [])
            joint_C.set_data([], [])
            com1.set_data([], [])
            com2.set_data([], [])
            com3.set_data([], [])
            time_text.set_text(time_template % t[i_frame] + " (Solver Failed)")
            return artists # Return all artists for blitting

        th1 = theta1_hist[i_frame]
        th2 = theta2_hist[i_frame]
        th3 = theta3_hist[i_frame]

        # Calculate joint positions
        Ax, Ay = 0, 0
        Bx = L1 * np.cos(th1)
        By = L1 * np.sin(th1)
        Dx, Dy = L4, 0
        # Calculate C from D using theta3 for consistency with how theta3 is defined
        Cx = Dx + L3 * np.cos(th3)
        Cy = Dy + L3 * np.sin(th3)
        # Note: Due to solver tolerance, C calculated from B might differ slightly.
        # Cx_alt = Bx + L2 * np.cos(th2)
        # Cy_alt = By + L2 * np.sin(th2)

        # Update link lines
        link1.set_data([Ax, Bx], [Ay, By])
        link2.set_data([Bx, Cx], [By, Cy])
        link3.set_data([Cx, Dx], [Cy, Dy])

        # Update joint markers
        joint_B.set_data([Bx], [By]) # Pass as lists or arrays
        joint_C.set_data([Cx], [Cy])

        # Update COM markers (using the pre-calculated history)
        if not np.isnan(r_C1_hist[i_frame, 0]):
            com1.set_data([r_C1_hist[i_frame, 0]], [r_C1_hist[i_frame, 1]])
        else:
            com1.set_data([],[])
        if not np.isnan(r_C2_hist[i_frame, 0]):
            com2.set_data([r_C2_hist[i_frame, 0]], [r_C2_hist[i_frame, 1]])
        else:
            com2.set_data([],[])
        if not np.isnan(r_C3_hist[i_frame, 0]):
            com3.set_data([r_C3_hist[i_frame, 0]], [r_C3_hist[i_frame, 1]])
        else:
            com3.set_data([],[])

        time_text.set_text(time_template % t[i_frame])

        # Return list of artists modified
        return artists

    # Create animation
    # Use a subset of frames if simulation is long/dense for smoother animation playback
    # Aim for roughly 100-200 frames in the animation.
    frame_skip = max(1, time_steps // 200)
    frames_to_animate = range(0, time_steps, frame_skip)

    # interval is delay between frames in milliseconds
    # dt is in seconds, so interval = dt * frame_skip * 1000
    ani = FuncAnimation(fig_anim, animate, frames=frames_to_animate,
                        interval=max(1, dt*frame_skip*1000), blit=True, repeat=False) # Ensure interval > 0

    # To save the animation (optional, requires ffmpeg or other writer)
    # print("Saving animation (may take a while)...")
    # try:
    #     ani.save('four_bar_linkage_kinematics.mp4', writer='ffmpeg', fps=30)
    #     print("Animation saved to four_bar_linkage_kinematics.mp4")
    # except Exception as e:
    #     print(f"Could not save animation: {e}")
    #     print("Showing animation directly instead.")

else:
    print("Skipping animation due to simulation failure or no valid data.")


# --- Display Plots and Animation ---
print("Showing plots (close figure windows to exit)...")
plt.show()

print("Script finished.")