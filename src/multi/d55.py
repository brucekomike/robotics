# %% Imports
import sympy as sp
import numpy as np
import matplotlib.pyplot as plt # Plotting needed
import matplotlib.animation as animation # Animation needed
from scipy.optimize import fsolve
from sympy.physics.mechanics import ReferenceFrame, Point
import time

print("Starting FIVE-bar linkage setup - VARIANT 2 (Two Drivers)...")
start_time = time.time()

# --------------------------------------------------------------------------
# 1. Symbolic Setup - FIVE BAR (Two Drivers)
# --------------------------------------------------------------------------
print("Setting up symbolic variables (a1, a2 are independent)...")

# Time and input
t = sp.Symbol('t')
# --- Independent Coordinates (a1, a2) ---
a1, a2 = sp.symbols('a1 a2')
da1, da2 = sp.symbols('da1 da2')
dda1, dda2 = sp.symbols('dda1 dda2')

# --- Dependent Coordinates ---
# Link 1
r11, r12 = sp.symbols('r11 r12')
dr11, dr12 = sp.symbols('dr11 dr12')
ddr11, ddr12 = sp.symbols('ddr11 ddr12')
# Link 2
r21, r22 = sp.symbols('r21 r22')
dr21, dr22 = sp.symbols('dr21 dr22')
ddr21, ddr22 = sp.symbols('ddr21 ddr22')
# Link 3
r31, r32, a3 = sp.symbols('r31 r32 a3')
dr31, dr32, da3 = sp.symbols('dr31 dr32 da3')
ddr31, ddr32, dda3 = sp.symbols('ddr31 ddr32 dda3')
# Link 4
r41, r42, a4 = sp.symbols('r41 r42 a4')
dr41, dr42, da4 = sp.symbols('dr41 dr42 da4')
ddr41, ddr42, dda4 = sp.symbols('ddr41 ddr42 dda4')

# Parameters
L1, L2, L3, L4, L5 = sp.symbols('L1 L2 L3 L4 L5') # Lengths

# --- Define generalized coordinates (q) and their derivatives (dq, ddq) ---
# Independent coordinates: q_i = [a1, a2] (2 DOFs specified)
# Dependent coordinates: q_d = [r11, r12, r21, r22, r31, r32, a3, r41, r42, a4] (10 vars)
q_i = sp.Matrix([a1, a2])
q_d = sp.Matrix([r11, r12, r21, r22, r31, r32, a3, r41, r42, a4])
q = sp.Matrix.vstack(q_i, q_d) # Full vector (12 vars)

dq_i = sp.Matrix([da1, da2])
dq_d = sp.Matrix([dr11, dr12, dr21, dr22, dr31, dr32, da3, dr41, dr42, da4])
dq = sp.Matrix.vstack(dq_i, dq_d) # Full vector (12 vars)

ddq_i = sp.Matrix([dda1, dda2])
ddq_d = sp.Matrix([ddr11, ddr12, ddr21, ddr22, ddr31, ddr32, dda3, ddr41, ddr42, dda4])
ddq = sp.Matrix.vstack(ddq_i, ddq_d) # Full vector (12 vars)

# Rotation Matrices
A1 = sp.Matrix([[sp.cos(a1), -sp.sin(a1)], [sp.sin(a1), sp.cos(a1)]])
A2 = sp.Matrix([[sp.cos(a2), -sp.sin(a2)], [sp.sin(a2), sp.cos(a2)]])
A3 = sp.Matrix([[sp.cos(a3), -sp.sin(a3)], [sp.sin(a3), sp.cos(a3)]])
A4 = sp.Matrix([[sp.cos(a4), -sp.sin(a4)], [sp.sin(a4), sp.cos(a4)]])

# Position Vectors (Centers of Mass)
R1 = sp.Matrix([r11, r12]); R2 = sp.Matrix([r21, r22])
R3 = sp.Matrix([r31, r32]); R4 = sp.Matrix([r41, r42])

# Local vectors from CoM to joints
u_1a = sp.Matrix([-L1/2, 0]); u_1b = sp.Matrix([ L1/2, 0])
u_2b = sp.Matrix([-L2/2, 0]); u_2c = sp.Matrix([ L2/2, 0])
u_3c = sp.Matrix([-L3/2, 0]); u_3d = sp.Matrix([ L3/2, 0])
u_4c = sp.Matrix([-L4/2, 0]); u_4d = sp.Matrix([ L4/2, 0])

# Fixed points
P0 = sp.Matrix([0, 0])
P3_p = sp.Matrix([L5, 0]) # Ground link length is L5

# --- Constraint Equations C(q) = 0 --- (Still 10 equations)
r1a = R1 + A1 * u_1a; fun1 = r1a - P0
r1b = R1 + A1 * u_1b; r2b = R2 + A2 * u_2b; fun2 = r2b - r1b
r2c = R2 + A2 * u_2c; r3c = R3 + A3 * u_3c; fun3 = r3c - r2c
r3d = R3 + A3 * u_3d; r4c = R4 + A4 * u_4c; fun4 = r4c - r3d
r4d = R4 + A4 * u_4d; fun5 = r4d - P3_p

func = sp.Matrix([fun1, fun2, fun3, fun4, fun5]) # Shape (10, 1)
print(f"Symbolic constraints defined ({func.shape[0]} equations).")
# sp.pprint(func) # Skip printing constraints
print("-" * 30)

# --------------------------------------------------------------------------
# 2. Symbolic Derivations for Kinematics
# --------------------------------------------------------------------------
print("Deriving symbolic Jacobians...")

# Jacobian of constraints w.r.t. all coordinates q (10x12)
Cq = func.jacobian(q)

# Partition Jacobian
Cq_d = Cq[:, 2:] # Jacobian w.r.t. q_d (10x10) -> SQUARE!
Cq_i = Cq[:, :2] # Jacobian w.r.t. q_i (a1, a2) (10x2)

print("Jacobians Cq_d and Cq_i derived.")
print(f"\nJacobian Cq_d (w.r.t. dependent coordinates): Shape {Cq_d.shape}")
print(f"\nJacobian Cq_i (w.r.t. independent coordinates): Shape {Cq_i.shape}")
print("-" * 30)

# --- Gamma Term Calculation ---
print("Calculating symbolic gamma term (may take significant time)...")
gamma = sp.zeros(func.shape[0], 1)
try:
    Cq_dq = Cq * dq
    gamma_term_jacobian = Cq_dq.jacobian(q)
    gamma = -gamma_term_jacobian * dq
    gamma = sp.simplify(gamma)
    print("Symbolic gamma term calculated.")
    print(f"\nGamma Term (gamma): Shape {gamma.shape}")
    print("-" * 30)
except Exception as e:
    print(f"WARNING: Symbolic calculation of gamma failed: {e}")
    gamma = sp.zeros(func.shape[0],1)
    print("\nGamma Term (gamma): Calculation Failed.")
    print("-" * 30)

# --------------------------------------------------------------------------
# 3. Lambdify Symbolic Expressions
# --------------------------------------------------------------------------
print("Lambdifying symbolic functions (constraints, Jacobians, gamma)...")

params_kin = [L1, L2, L3, L4, L5]
state_q = list(q) # 12 vars
state_dq = list(dq) # 12 vars

# Lambdify constraint function (12 state vars + 5 params)
func_num = sp.lambdify(state_q + params_kin, func, modules=['numpy'])
print("Lambdified: func_num")

# Lambdify Jacobians (12 state vars + 5 params)
Cq_d_num = sp.lambdify(state_q + params_kin, Cq_d, modules=['numpy'])
Cq_i_num = sp.lambdify(state_q + params_kin, Cq_i, modules=['numpy'])
print("Lambdified: Cq_d_num (10x10), Cq_i_num (10x2)")

# Lambdify gamma term (12 state + 12 velocity + 5 params)
gamma_num = sp.lambdify(state_q + state_dq + params_kin, gamma, modules=['numpy'])
print("Lambdified: gamma_num")

setup_time = time.time() - start_time
print(f"Symbolic setup and lambdification complete in {setup_time:.2f} seconds.")
print("-" * 20)

# --------------------------------------------------------------------------
# 4. Numerical Simulation Parameters
# --------------------------------------------------------------------------
print("Setting up simulation parameters...")
t_start = 0.0; t_end = 10.0; dt = 0.02
t_m = np.arange(t_start, t_end + dt, dt)
nn = len(t_m)
L1_val, L2_val, L3_val, L4_val, L5_val = 0.7, 1.4, 1.3, 1.2, 1.5 # Example lengths
param_vals_kin_num = [L1_val, L2_val, L3_val, L4_val, L5_val]
w1 = 100.0 # Freq for a1
#w2 = 5.0
w2 = 50.0  # Freq for a2

print(f"Time steps: {nn}")
print("-" * 20)

# --------------------------------------------------------------------------
# 5. Define Input Motion (a1 and a2)
# --------------------------------------------------------------------------
print("Defining input motion functions (a1 and a2)...")
# a1(t)
def input_a1(t, w): return np.exp(-t/5) * np.cos(t * w) # Slower decay
def input_da1(t, w):
    exp_t = np.exp(-t/5); return -0.2*exp_t*np.cos(t*w) - w*exp_t*np.sin(t*w)
def input_dda1(t, w):
    exp_t = np.exp(-t/5); cos_tw = np.cos(t*w); sin_tw = np.sin(t*w); w_sq = w**2
    return (0.04 - w_sq)*exp_t*cos_tw + 0.4*w*exp_t*sin_tw

# a2(t) - Example: Constant velocity + small oscillation
amp2 = np.pi / 6
vel2 = np.pi / 4
def input_a2(t, w): return vel2 * t + amp2 * np.sin(t * w)
def input_da2(t, w): return vel2 + amp2 * w * np.cos(t * w)
def input_dda2(t, w): return -amp2 * w**2 * np.sin(t * w)

print("Input motion functions defined.")
print("-" * 20)

# --------------------------------------------------------------------------
# 6. Run Numerical Simulation Loop
# --------------------------------------------------------------------------
print("Running numerical simulation...")

# Initialize storage arrays
results = {}
var_names_q = ['a1', 'a2', 'r11', 'r12', 'r21', 'r22', 'r31', 'r32', 'a3', 'r41', 'r42', 'a4']
var_names_dq = ['da1', 'da2', 'dr11', 'dr12', 'dr21', 'dr22', 'dr31', 'dr32', 'da3', 'dr41', 'dr42', 'da4']
var_names_ddq = ['dda1', 'dda2', 'ddr11', 'ddr12', 'ddr21', 'ddr22', 'ddr31', 'ddr32', 'dda3', 'ddr41', 'ddr42', 'dda4']
all_var_names = var_names_q + var_names_dq + var_names_ddq

for name in all_var_names:
    results[name + '_m'] = np.full(nn, np.nan)
results['fsolve_info'] = [{}] * nn

# --- Initial Guess (q_d - 10 vars) ---
a1_0 = input_a1(t_start, w1)
a2_0 = input_a2(t_start, w2)
# Rough guess (similar to before, but using a2_0)
r11_g, r12_g = L1_val/2*np.cos(a1_0), L1_val/2*np.sin(a1_0)
r21_g, r22_g = r11_g+L1_val/2*np.cos(a1_0)+L2_val/2*np.cos(a2_0), r12_g+L1_val/2*np.sin(a1_0)+L2_val/2*np.sin(a2_0)
a3_g = a1_0 + a2_0 # Arbitrary combination
r31_g, r32_g = r21_g+L2_val/2*np.cos(a2_0)+L3_val/2*np.cos(a3_g), r22_g+L2_val/2*np.sin(a2_0)+L3_val/2*np.sin(a3_g)
a4_g = a2_0 # Arbitrary
r41_g, r42_g = L5_val-L4_val/2*np.cos(a4_g), 0.0-L4_val/2*np.sin(a4_g)
q_d_guess = np.array([r11_g, r12_g, r21_g, r22_g, r31_g, r32_g, a3_g, r41_g, r42_g, a4_g])

print(f"Attempting to find initial configuration for a1(0)={a1_0:.4f}, a2(0)={a2_0:.4f}")

# Function for fsolve (10 variables, 10 equations)
def position_residuals_2dof(qd_unknown, a1_val, a2_val, params_k):
    # qd_unknown has 10 elements
    q_full = np.concatenate(([a1_val, a2_val], qd_unknown)) # Construct full q (12 elements)
    res = func_num(*q_full, *params_k) # Returns 10 residuals
    return res.flatten()

# Solve for initial q_d
sol_init, info_init, ier_init, msg_init = fsolve(position_residuals_2dof, q_d_guess, args=(a1_0, a2_0, param_vals_kin_num), full_output=True)

if ier_init == 1:
    q_d_guess = sol_init
    print(f"Found initial q_d: {q_d_guess}")
else:
    print(f"ERROR: Position solver failed for initial condition. Msg: {msg_init}")
    exit()

# --- Simulation Loop ---
sim_start_time = time.time()
for k in range(nn):
    t_k = t_m[k]
    # Get independent inputs
    a1_k = input_a1(t_k, w1)
    da1_k = input_da1(t_k, w1)
    dda1_k = input_dda1(t_k, w1)
    a2_k = input_a2(t_k, w2)
    da2_k = input_da2(t_k, w2)
    dda2_k = input_dda2(t_k, w2)

    # Store independent variables
    results['a1_m'][k] = a1_k; results['da1_m'][k] = da1_k; results['dda1_m'][k] = dda1_k
    results['a2_m'][k] = a2_k; results['da2_m'][k] = da2_k; results['dda2_m'][k] = dda2_k
    dq_i_k = np.array([[da1_k], [da2_k]]) # 2x1 vector
    ddq_i_k = np.array([[dda1_k], [dda2_k]]) # 2x1 vector

    # --- Position Analysis ---
    sol, info, ier, msg = fsolve(position_residuals_2dof, q_d_guess, args=(a1_k, a2_k, param_vals_kin_num), full_output=True)

    if ier == 1:
        q_d_k = sol
        q_d_guess = q_d_k
        results['r11_m'][k], results['r12_m'][k], results['r21_m'][k], results['r22_m'][k], \
        results['r31_m'][k], results['r32_m'][k], results['a3_m'][k], \
        results['r41_m'][k], results['r42_m'][k], results['a4_m'][k] = q_d_k
    else:
        print(f"WARNING: Position solver failed at step k={k}, t={t_k:.2f}. Msg: {msg}")
        break

    results['fsolve_info'][k] = info
    q_k = np.concatenate(([a1_k, a2_k], q_d_k)) # Full q (12 elements)

    # --- Velocity Analysis ---
    try:
        Cq_d_k = Cq_d_num(*q_k, *param_vals_kin_num) # 10x10
        Cq_i_k = Cq_i_num(*q_k, *param_vals_kin_num) # 10x2

        cond_num = np.linalg.cond(Cq_d_k)
        if cond_num > 1e9:
             print(f"WARNING: Cq_d ill-conditioned at k={k}, t={t_k:.2f} (cond={cond_num:.2e})")

        # Solve Cq_d * dq_d = -Cq_i * dq_i
        rhs_vel = -Cq_i_k @ dq_i_k # 10x2 @ 2x1 -> 10x1
        dq_d_k = np.linalg.solve(Cq_d_k, rhs_vel).flatten() # Solve for 10 velocities

        results['dr11_m'][k], results['dr12_m'][k], results['dr21_m'][k], results['dr22_m'][k], \
        results['dr31_m'][k], results['dr32_m'][k], results['da3_m'][k], \
        results['dr41_m'][k], results['dr42_m'][k], results['da4_m'][k] = dq_d_k

    except np.linalg.LinAlgError:
        print(f"ERROR: Velocity analysis failed at k={k}, t={t_k:.2f} (Singular Cq_d)")
        break
    except Exception as e:
        print(f"ERROR: Velocity analysis failed at k={k}, t={t_k:.2f}. Error: {e}")
        break

    dq_k = np.concatenate((dq_i_k.flatten(), dq_d_k)) # Full dq (12 elements)

    # --- Acceleration Analysis ---
    try:
        gamma_k = gamma_num(*q_k, *dq_k, *param_vals_kin_num) # 10x1

        # Solve Cq_d * ddq_d = gamma - Cq_i * ddq_i
        rhs_accel = gamma_k - Cq_i_k @ ddq_i_k # 10x1 - (10x2 @ 2x1) -> 10x1
        ddq_d_k = np.linalg.solve(Cq_d_k, rhs_accel).flatten() # Solve for 10 accelerations

        results['ddr11_m'][k], results['ddr12_m'][k], results['ddr21_m'][k], results['ddr22_m'][k], \
        results['ddr31_m'][k], results['ddr32_m'][k], results['dda3_m'][k], \
        results['ddr41_m'][k], results['ddr42_m'][k], results['dda4_m'][k] = ddq_d_k

    except np.linalg.LinAlgError:
        print(f"ERROR: Acceleration analysis failed at k={k}, t={t_k:.2f} (Singular Cq_d)")
        break
    except Exception as e:
        print(f"ERROR: Acceleration analysis failed at k={k}, t={t_k:.2f}. Error: {e}")
        if 'gamma_num' in str(e): print("   -> Likely symbolic gamma issue.")
        break

    # Progress update
    if (k + 1) % (max(1, nn // 20)) == 0:
        print(f"Simulation progress: {100 * (k + 1) / nn:.0f}% (t={t_k:.2f}s)")

sim_time = time.time() - sim_start_time
print(f"Numerical simulation loop finished in {sim_time:.2f} seconds.")
print("-" * 20)

# Find last valid index (check a dependent variable like a4)
last_valid_index = np.where(~np.isnan(results['a4_m']))[0]
if len(last_valid_index) > 0:
    nn_valid = last_valid_index[-1] + 1
    if nn_valid < nn:
        print(f"Simulation stopped early. Only {nn_valid} / {nn} steps completed.")
        t_m = t_m[:nn_valid]
        for name in results:
             if name != 'fsolve_info': results[name] = results[name][:nn_valid]
             else: results[name] = results[name][:nn_valid]
        nn = nn_valid
    else:
        print("All simulation steps completed successfully.")
else:
    print("ERROR: No steps completed successfully.")
    exit()

# --------------------------------------------------------------------------
# 7. Animation Setup - FIVE BAR (Two Drivers)
# --------------------------------------------------------------------------
print("Setting up animation figure...")

fig_anim, ax_anim = plt.subplots()
ax_anim.set_aspect('equal')
max_reach = L1_val + L2_val + L3_val + L4_val
ax_anim.set_xlim(-0.5 * (L1_val + L2_val), L5_val + 0.5 * L4_val + 0.5)
ax_anim.set_ylim(-0.6 * max_reach, 0.6 * max_reach) # Adjusted ylim
ax_anim.set_xlabel("X position")
ax_anim.set_ylabel("Y position")
ax_anim.set_title("Five-Bar Linkage Animation (Inputs: a1, a2)")
ax_anim.grid(True)

P0_np = np.array([0, 0])
P3_ground_np = np.array([L5_val, 0])

line1, = ax_anim.plot([], [], 'o-', lw=2, color='blue', label='L1')
line2, = ax_anim.plot([], [], 'o-', lw=2, color='red', label='L2')
line3, = ax_anim.plot([], [], 'o-', lw=2, color='green', label='L3')
line4, = ax_anim.plot([], [], 'o-', lw=2, color='purple', label='L4')
line5, = ax_anim.plot([P0_np[0], P3_ground_np[0]], [P0_np[1], P3_ground_np[1]], 's-', lw=1, color='black', label=f'L5 (Gnd)')

com1, = ax_anim.plot([], [], 'bo', ms=4, label='CoM1'); com2, = ax_anim.plot([], [], 'ro', ms=4, label='CoM2')
com3, = ax_anim.plot([], [], 'go', ms=4, label='CoM3'); com4, = ax_anim.plot([], [], 'mo', ms=4, label='CoM4')

time_template = 'time = %.2fs'
time_text = ax_anim.text(0.05, 0.9, '', transform=ax_anim.transAxes)
ax_anim.legend(loc='upper right')

def get_joint_positions_5bar_2dof(k):
    """Calculates joint positions for frame k (uses a2 from results)."""
    r11_k=results['r11_m'][k]; r12_k=results['r12_m'][k]; a1_k=results['a1_m'][k]
    r21_k=results['r21_m'][k]; r22_k=results['r22_m'][k]; a2_k=results['a2_m'][k] # Use stored a2
    r31_k=results['r31_m'][k]; r32_k=results['r32_m'][k]; a3_k=results['a3_m'][k]
    r41_k=results['r41_m'][k]; r42_k=results['r42_m'][k]; a4_k=results['a4_m'][k]

    A1_k = np.array([[np.cos(a1_k),-np.sin(a1_k)],[np.sin(a1_k),np.cos(a1_k)]])
    A2_k = np.array([[np.cos(a2_k),-np.sin(a2_k)],[np.sin(a2_k),np.cos(a2_k)]]) # Use a2_k
    A3_k = np.array([[np.cos(a3_k),-np.sin(a3_k)],[np.sin(a3_k),np.cos(a3_k)]])
    A4_k = np.array([[np.cos(a4_k),-np.sin(a4_k)],[np.sin(a4_k),np.cos(a4_k)]])

    u_1b_local = np.array([0.5*L1_val,0]); u_2c_local = np.array([0.5*L2_val,0])
    u_3d_local = np.array([0.5*L3_val,0])

    R1_k = np.array([r11_k,r12_k]); R2_k = np.array([r21_k,r22_k])
    R3_k = np.array([r31_k,r32_k]); R4_k = np.array([r41_k,r42_k])

    P1_k = R1_k + A1_k @ u_1b_local
    P2_k = R2_k + A2_k @ u_2c_local
    P3j_k = R3_k + A3_k @ u_3d_local

    CoM1_k=R1_k; CoM2_k=R2_k; CoM3_k=R3_k; CoM4_k=R4_k
    return P1_k, P2_k, P3j_k, CoM1_k, CoM2_k, CoM3_k, CoM4_k

def init_anim():
    line1.set_data([], []); line2.set_data([], []); line3.set_data([], [])
    line4.set_data([], []); com1.set_data([], []); com2.set_data([], [])
    com3.set_data([], []); com4.set_data([], []); time_text.set_text('')
    return line1, line2, line3, line4, com1, com2, com3, com4, time_text

def update_anim(frame):
    if frame >= nn: return line1, line2, line3, line4, com1, com2, com3, com4, time_text
    P1_k, P2_k, P3j_k, CoM1_k, CoM2_k, CoM3_k, CoM4_k = get_joint_positions_5bar_2dof(frame)
    line1.set_data([P0_np[0], P1_k[0]], [P0_np[1], P1_k[1]])
    line2.set_data([P1_k[0], P2_k[0]], [P1_k[1], P2_k[1]])
    line3.set_data([P2_k[0], P3j_k[0]], [P2_k[1], P3j_k[1]])
    line4.set_data([P3j_k[0], P3_ground_np[0]], [P3j_k[1], P3_ground_np[1]])
    com1.set_data([CoM1_k[0]], [CoM1_k[1]]); com2.set_data([CoM2_k[0]], [CoM2_k[1]])
    com3.set_data([CoM3_k[0]], [CoM3_k[1]]); com4.set_data([CoM4_k[0]], [CoM4_k[1]])
    time_text.set_text(time_template % t_m[frame])
    return line1, line2, line3, line4, com1, com2, com3, com4, time_text

# Create animation object
ani = animation.FuncAnimation(fig_anim, update_anim, frames=nn,
                              init_func=init_anim, blit=True, interval=dt*1000)

print("Animation object created.")
print("-" * 20)

# --------------------------------------------------------------------------
# 8. Display Static Plots and Start Animation (Revised Pause)
# --------------------------------------------------------------------------
print("Generating static plots...")

# Angles a1, a2 vs Time
plt.figure()
plt.plot(t_m, results['a1_m'] * 180 / np.pi, label='a1')
plt.plot(t_m, results['a2_m'] * 180 / np.pi, label='a2')
plt.xlabel('Time (s)'); plt.ylabel('Angles (degrees)')
plt.title('Input Angles vs Time (Five-Bar, 2 Drivers)')
plt.legend(); plt.grid(True)

# Angles a3, a4 vs Time
plt.figure()
plt.plot(t_m, results['a3_m'] * 180 / np.pi, label='a3')
plt.plot(t_m, results['a4_m'] * 180 / np.pi, label='a4')
plt.xlabel('Time (s)'); plt.ylabel('Angles (degrees)')
plt.title('Dependent Angles vs Time (Five-Bar, 2 Drivers)')
plt.legend(); plt.grid(True)

# Link 3 CoM Trajectory
plt.figure()
plt.plot(results['r31_m'], results['r32_m'])
plt.xlabel('r31 (m)'); plt.ylabel('r32 (m)')
plt.title('Link 3 CoM Trajectory (Five-Bar, 2 Drivers)')
plt.axis('equal'); plt.grid(True)

# Link 4 CoM Trajectory
plt.figure()
plt.plot(results['r41_m'], results['r42_m'])
plt.xlabel('r41 (m)'); plt.ylabel('r42 (m)')
plt.title('Link 4 CoM Trajectory (Five-Bar, 2 Drivers)')
plt.axis('equal'); plt.grid(True)

# Input Kinematics a1 vs Time (Subplots)
fig_kin1, axs_kin1 = plt.subplots(3, 1, sharex=True); fig_kin1.suptitle('Input a1 Kinematics')
axs_kin1[0].plot(t_m, results['dda1_m'] * 180 / np.pi); axs_kin1[0].set_ylabel('dda1 (deg/s^2)'); axs_kin1[0].grid(True)
axs_kin1[1].plot(t_m, results['da1_m'] * 180 / np.pi); axs_kin1[1].set_ylabel('da1 (deg/s)'); axs_kin1[1].grid(True)
axs_kin1[2].plot(t_m, results['a1_m'] * 180 / np.pi); axs_kin1[2].set_ylabel('a1 (deg)'); axs_kin1[2].grid(True)
axs_kin1[2].set_xlabel('Time (s)')

# Input Kinematics a2 vs Time (Subplots)
fig_kin2, axs_kin2 = plt.subplots(3, 1, sharex=True); fig_kin2.suptitle('Input a2 Kinematics')
axs_kin2[0].plot(t_m, results['dda2_m'] * 180 / np.pi); axs_kin2[0].set_ylabel('dda2 (deg/s^2)'); axs_kin2[0].grid(True)
axs_kin2[1].plot(t_m, results['da2_m'] * 180 / np.pi); axs_kin2[1].set_ylabel('da2 (deg/s)'); axs_kin2[1].grid(True)
axs_kin2[2].plot(t_m, results['a2_m'] * 180 / np.pi); axs_kin2[2].set_ylabel('a2 (deg)'); axs_kin2[2].grid(True)
axs_kin2[2].set_xlabel('Time (s)')


# --- Show static plots FIRST (blocking) ---
print("Displaying static plots... Close all plot windows to continue.")
plt.show() # This blocks until windows are closed

# --- Wait for key press before showing animation ---
print("-" * 20)
input("Static plots closed. Press Enter to start the animation...")
print("Starting animation...")

# --- Show animation figure (blocking) ---
# Need to show the specific animation figure now
plt.figure(fig_anim.number) # Bring animation figure to front
plt.show() # Start event loop for animation

print("Script finished.")