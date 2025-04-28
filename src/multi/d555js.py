# %% Imports
import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.optimize import fsolve
from sympy.physics.mechanics import ReferenceFrame, Point
import time
import sys # To exit cleanly on error

# --------------------------------------------------------------------------
# <<< CONFIGURATION SECTION >>>
# --------------------------------------------------------------------------

# --- Choose Driven Joints ---
# Select EXACTLY TWO angles from ['a1', 'a2', 'a3', 'a4'] to be the drivers
DRIVEN_JOINTS = ['a1', 'a2'] # Example 1: Original configuration
#DRIVEN_JOINTS = ['a1', 'a3'] # Example 2: Drive link 1 and link 3
#DRIVEN_JOINTS = ['a2', 'a4'] # Example 3: Drive link 2 and link 4
# DRIVEN_JOINTS = ['a1', 'a4'] # Example 4: Drive link 1 and link 4

# --- Link Lengths ---
L1_val = 1.3
L2_val = 0.6
L3_val = 1.8
L4_val = 1.8
L5_val = 2.9 # Ground link length

# --- Simulation Time ---
t_start = 0.0
t_end = 10.0 # Duration in seconds
dt = 0.001  # Time step in seconds

# --- Input Motion Parameters (Define for ALL potential drivers) ---
# 'w':   Angular frequency (rad/s) - Controls "spinning speed"
# 'amp': Amplitude of oscillation (rad)
# 'vel': Constant angular velocity component (rad/s)
# 'decay': Exponential decay time constant (s) (0 or large value means no decay)
input_params = {
    'a1': {'w': 20.0,  'amp': 1.0,       'decay': 5.0, 'vel': 0.0},       # Slower spin, large swing with decay
    'a2': {'w': 300.0, 'amp': np.pi/10,  'decay': 0.0, 'vel': np.pi/0.1},   # Faster spin, small oscillation, constant drift
    'a3': {'w': 8.0,  'amp': np.pi/6,   'decay': 8.0, 'vel': -np.pi/7},  # Medium spin, medium swing with slow decay, drift back
    'a4': {'w': 3.0,  'amp': np.pi/4,   'decay': 0.0, 'vel': np.pi/8}    # Slow spin, large swing, slow drift
}

# --------------------------------------------------------------------------
# <<< END CONFIGURATION SECTION >>>
# --------------------------------------------------------------------------

# --- Validate Configuration ---
if len(DRIVEN_JOINTS) != 2 or not all(j in ['a1', 'a2', 'a3', 'a4'] for j in DRIVEN_JOINTS):
    print("ERROR: Please specify exactly two distinct driven joints from ['a1', 'a2', 'a3', 'a4'] in DRIVEN_JOINTS list.")
    sys.exit()
if not all(j in input_params for j in DRIVEN_JOINTS):
     print(f"ERROR: Input parameters missing for one or both driven joints: {DRIVEN_JOINTS}")
     sys.exit()
print(f"Selected driven joints: {DRIVEN_JOINTS}")
print(f"Link Lengths: L1={L1_val}, L2={L2_val}, L3={L3_val}, L4={L4_val}, L5={L5_val}")
print(f"Simulation Time: {t_start}s to {t_end}s, dt={dt}s")
print("-" * 30)

# --------------------------------------------------------------------------
# 1. Symbolic Setup - FIVE BAR (Flexible Drivers)
# --------------------------------------------------------------------------
print("Setting up symbolic variables (all angles initially symbolic)...")
script_start_time = time.time() # Use a different variable name

# Time
t = sp.Symbol('t')
# --- Define ALL potential generalized coordinates ---
a1, a2, a3, a4 = sp.symbols('a1 a2 a3 a4')
r11, r12 = sp.symbols('r11 r12'); r21, r22 = sp.symbols('r21 r22')
r31, r32 = sp.symbols('r31 r32'); r41, r42 = sp.symbols('r41 r42')

# --- Full coordinate vector q_all (12 variables) ---
q_all_sym = sp.Matrix([a1, a2, a3, a4, r11, r12, r21, r22, r31, r32, r41, r42])
coord_names = [str(s) for s in q_all_sym]

# --- Define ALL corresponding velocities and accelerations ---
dq_all_sym = sp.Matrix([sp.Symbol('d' + s) for s in coord_names])
ddq_all_sym = sp.Matrix([sp.Symbol('dd' + s) for s in coord_names])

# Parameters
L1, L2, L3, L4, L5 = sp.symbols('L1 L2 L3 L4 L5') # Symbolic Lengths

# --- Helper function to get indices based on DRIVEN_JOINTS ---
def get_coord_indices(all_coords_list, driven_list):
    driven_indices = sorted([all_coords_list.index(j) for j in driven_list])
    dependent_indices = sorted([i for i, name in enumerate(all_coords_list) if name not in driven_list])
    return driven_indices, dependent_indices

driven_idx, dependent_idx = get_coord_indices(coord_names, DRIVEN_JOINTS)
print(f"Driven coordinate indices: {driven_idx}")
print(f"Dependent coordinate indices: {dependent_idx}")

# Rotation Matrices
A1 = sp.Matrix([[sp.cos(a1), -sp.sin(a1)], [sp.sin(a1), sp.cos(a1)]])
A2 = sp.Matrix([[sp.cos(a2), -sp.sin(a2)], [sp.sin(a2), sp.cos(a2)]])
A3 = sp.Matrix([[sp.cos(a3), -sp.sin(a3)], [sp.sin(a3), sp.cos(a3)]])
A4 = sp.Matrix([[sp.cos(a4), -sp.sin(a4)], [sp.sin(a4), sp.cos(a4)]])

# Position Vectors
R1 = sp.Matrix([r11, r12]); R2 = sp.Matrix([r21, r22])
R3 = sp.Matrix([r31, r32]); R4 = sp.Matrix([r41, r42])

# Local vectors
u_1a = sp.Matrix([-L1/2, 0]); u_1b = sp.Matrix([ L1/2, 0])
u_2b = sp.Matrix([-L2/2, 0]); u_2c = sp.Matrix([ L2/2, 0])
u_3c = sp.Matrix([-L3/2, 0]); u_3d = sp.Matrix([ L3/2, 0])
u_4c = sp.Matrix([-L4/2, 0]); u_4d = sp.Matrix([ L4/2, 0])

# Fixed points
P0 = sp.Matrix([0, 0]); P3_p = sp.Matrix([L5, 0])

# --- Constraint Equations C(q) = 0 ---
r1a = R1 + A1 * u_1a; fun1 = r1a - P0
r1b = R1 + A1 * u_1b; r2b = R2 + A2 * u_2b; fun2 = r2b - r1b
r2c = R2 + A2 * u_2c; r3c = R3 + A3 * u_3c; fun3 = r3c - r2c
r3d = R3 + A3 * u_3d; r4c = R4 + A4 * u_4c; fun4 = r4c - r3d
r4d = R4 + A4 * u_4d; fun5 = r4d - P3_p
func = sp.Matrix([fun1, fun2, fun3, fun4, fun5])
print(f"Symbolic constraints defined ({func.shape[0]} equations).")
print("-" * 30)

# --------------------------------------------------------------------------
# 2. Symbolic Derivations for Kinematics (using q_all)
# --------------------------------------------------------------------------
print("Deriving symbolic FULL Jacobian Cq...")
Cq_full = func.jacobian(q_all_sym)
print(f"Full Jacobian Cq_full derived. Shape: {Cq_full.shape}")
print("-" * 30)

print("Calculating symbolic FULL gamma term (may take significant time)...")
gamma_full = sp.zeros(func.shape[0], 1)
try:
    # Optional: Use CSE for potentially faster Jacobian calculation if needed
    # Cq_dq_full_rep, Cq_dq_full_red = sp.cse(Cq_full * dq_all_sym)
    # Cq_dq_full_opt = sp.Matrix([r[1] for r in Cq_dq_full_red])
    # gamma_term_jacobian_full = Cq_dq_full_opt.jacobian(q_all_sym)
    Cq_dq_full = Cq_full * dq_all_sym
    gamma_term_jacobian_full = Cq_dq_full.jacobian(q_all_sym)
    gamma_full = -gamma_term_jacobian_full * dq_all_sym
    gamma_full = sp.simplify(gamma_full)
    print("Symbolic FULL gamma term calculated.")
    print(f"\nGamma Term (gamma_full): Shape {gamma_full.shape}")
    print("-" * 30)
except Exception as e:
    print(f"WARNING: Symbolic calculation of gamma failed: {e}")
    gamma_full = sp.zeros(func.shape[0],1)
    print("\nGamma Term (gamma_full): Calculation Failed.")
    print("-" * 30)

# --------------------------------------------------------------------------
# 3. Lambdify Symbolic Expressions (Full versions)
# --------------------------------------------------------------------------
print("Lambdifying symbolic functions (constraints, FULL Jacobian, FULL gamma)...")

params_kin_sym = [L1, L2, L3, L4, L5] # Symbolic params for lambdify
state_q_all = list(q_all_sym)
state_dq_all = list(dq_all_sym)

func_num = sp.lambdify(state_q_all + params_kin_sym, func, modules=['numpy'])
print("Lambdified: func_num")
Cq_full_num = sp.lambdify(state_q_all + params_kin_sym, Cq_full, modules=['numpy'])
print("Lambdified: Cq_full_num (10x12)")
gamma_full_num = sp.lambdify(state_q_all + state_dq_all + params_kin_sym, gamma_full, modules=['numpy'])
print("Lambdified: gamma_full_num (10x1)")

setup_time = time.time() - script_start_time
print(f"Symbolic setup and lambdification complete in {setup_time:.2f} seconds.")
print("-" * 20)

# --------------------------------------------------------------------------
# 4. Numerical Simulation Parameters (Using Config Vars)
# --------------------------------------------------------------------------
print("Setting up simulation parameters from configuration...")
t_m = np.arange(t_start, t_end + dt, dt)
nn = len(t_m)
param_vals_kin_num = [L1_val, L2_val, L3_val, L4_val, L5_val] # Numerical values

print(f"Time steps: {nn}")
print("-" * 20)

# --------------------------------------------------------------------------
# 5. Define Input Motion (Flexible for any angle)
# --------------------------------------------------------------------------
print(f"Defining input motion functions for driven joints: {DRIVEN_JOINTS}...")

def get_input_motion(t, name, params_dict):
    """Calculates angle, velocity, acceleration for a given joint name."""
    # Fetch parameters for the specific joint 'name'
    params = params_dict[name]
    w = params['w']
    amp = params['amp']
    decay = params['decay']
    vel = params['vel']

    if decay > 1e-6: # Exponential decay + Sine + Constant Velocity
        exp_t = np.exp(-t / decay)
        cos_tw = np.cos(t * w)
        sin_tw = np.sin(t * w)
        decay_inv = 1.0 / decay
        angle = vel * t + amp * exp_t * cos_tw
        ang_vel = vel + amp * exp_t * (-decay_inv * cos_tw - w * sin_tw)
        ang_acc = amp * exp_t * ((decay_inv**2 - w**2) * cos_tw + 2 * w * decay_inv * sin_tw)
    else: # Constant Velocity + Sine (no decay)
        angle = vel * t + amp * np.sin(t * w)
        ang_vel = vel + amp * w * np.cos(t * w)
        ang_acc = -amp * w**2 * np.sin(t * w)
    return angle, ang_vel, ang_acc

# Get the parameters for the chosen driven joints from the config dict
driver1_name = DRIVEN_JOINTS[0]
driver2_name = DRIVEN_JOINTS[1]
# No need to store params separately, just use input_params[driver_name]

print("Input motion functions defined.")
print("-" * 20)

# --------------------------------------------------------------------------
# 6. Run Numerical Simulation Loop
# --------------------------------------------------------------------------
print("Running numerical simulation...")

# Initialize storage arrays for ALL coordinates
results = {}
for name in coord_names:
    results[name + '_m'] = np.full(nn, np.nan)
    results['d' + name + '_m'] = np.full(nn, np.nan)
    results['dd' + name + '_m'] = np.full(nn, np.nan)
results['fsolve_info'] = [{}] * nn

# --- Initial Guess (q_d - 10 vars) ---
# Get initial values for *all* potential drivers to build the guess
q_all_guess_t0 = np.zeros(12)
for name in ['a1', 'a2', 'a3', 'a4']:
    angle_t0, _, _ = get_input_motion(t_start, name, input_params)
    q_all_guess_t0[coord_names.index(name)] = angle_t0
a1_t0 = q_all_guess_t0[coord_names.index('a1')] # Extract for clarity
a2_t0 = q_all_guess_t0[coord_names.index('a2')]
a3_t0 = q_all_guess_t0[coord_names.index('a3')]
a4_t0 = q_all_guess_t0[coord_names.index('a4')]

# Rough guess for r variables (same as before)
q_all_guess_t0[coord_names.index('r11')] = L1_val/2*np.cos(a1_t0)
q_all_guess_t0[coord_names.index('r12')] = L1_val/2*np.sin(a1_t0)
q_all_guess_t0[coord_names.index('r21')] = q_all_guess_t0[coord_names.index('r11')] + L1_val/2*np.cos(a1_t0)+L2_val/2*np.cos(a2_t0)
q_all_guess_t0[coord_names.index('r22')] = q_all_guess_t0[coord_names.index('r12')] + L1_val/2*np.sin(a1_t0)+L2_val/2*np.sin(a2_t0)
q_all_guess_t0[coord_names.index('r31')] = q_all_guess_t0[coord_names.index('r21')] + L2_val/2*np.cos(a2_t0)+L3_val/2*np.cos(a3_t0)
q_all_guess_t0[coord_names.index('r32')] = q_all_guess_t0[coord_names.index('r22')] + L2_val/2*np.sin(a2_t0)+L3_val/2*np.sin(a3_t0)
q_all_guess_t0[coord_names.index('r41')] = L5_val-L4_val/2*np.cos(a4_t0)
q_all_guess_t0[coord_names.index('r42')] = 0.0-L4_val/2*np.sin(a4_t0)

# Extract the guess for the *dependent* variables
q_d_guess = q_all_guess_t0[dependent_idx]
# Get the known values for the *independent* variables at t=0
q_i_t0 = q_all_guess_t0[driven_idx]

print(f"Attempting to find initial configuration for {driver1_name}(0)={q_i_t0[0]:.4f}, {driver2_name}(0)={q_i_t0[1]:.4f}")

# Function for fsolve (10 variables, 10 equations)
def position_residuals_flexible(qd_unknown, qi_known, params_k, driven_indices, dependent_indices, num_all_coords):
    q_full = np.zeros(num_all_coords)
    q_full[driven_indices] = qi_known
    q_full[dependent_indices] = qd_unknown
    res = func_num(*q_full, *params_k)
    return res.flatten()

# Solve for initial q_d
sol_init, info_init, ier_init, msg_init = fsolve(
    position_residuals_flexible, q_d_guess,
    args=(q_i_t0, param_vals_kin_num, driven_idx, dependent_idx, len(q_all_sym)),
    full_output=True, xtol=1e-8 # Adjust tolerance if needed
)

if ier_init == 1:
    q_d_guess = sol_init
    print(f"Found initial dependent q_d: {q_d_guess}")
else:
    print(f"ERROR: Position solver failed for initial condition. Msg: {msg_init}")
    print("Check linkage dimensions, initial guess, or driver combination/motion.")
    sys.exit()

# Store the complete initial state q_all_k
q_all_k = np.zeros(len(q_all_sym))
q_all_k[driven_idx] = q_i_t0
q_all_k[dependent_idx] = q_d_guess

# --- Simulation Loop ---
sim_start_time = time.time()
# Initialize dq_all_k for gamma calculation
dq_all_k = np.zeros(len(q_all_sym))
_, dq_i_t0_val1, _ = get_input_motion(t_start, driver1_name, input_params)
_, dq_i_t0_val2, _ = get_input_motion(t_start, driver2_name, input_params)
dq_all_k[driven_idx] = [dq_i_t0_val1, dq_i_t0_val2]
# Solve for initial dq_d (optional but good practice)
try:
    Cq_full_k0 = Cq_full_num(*q_all_k, *param_vals_kin_num)
    Cq_d_k0 = Cq_full_k0[:, dependent_idx]
    Cq_i_k0 = Cq_full_k0[:, driven_idx]
    rhs_vel0 = -Cq_i_k0 @ dq_all_k[driven_idx].reshape(-1,1)
    dq_d_k0 = np.linalg.solve(Cq_d_k0, rhs_vel0).flatten()
    dq_all_k[dependent_idx] = dq_d_k0
    print("Calculated initial velocities for dependent coordinates.")
except np.linalg.LinAlgError:
    print("Warning: Could not solve for initial dependent velocities (singular Cq_d at t=0?). Using zeros.")


for k in range(nn):
    t_k = t_m[k]
    # Get independent inputs q_i, dq_i, ddq_i for this step
    a_i_k = np.zeros(2); da_i_k = np.zeros(2); dda_i_k = np.zeros(2)
    a_i_k[0], da_i_k[0], dda_i_k[0] = get_input_motion(t_k, driver1_name, input_params)
    a_i_k[1], da_i_k[1], dda_i_k[1] = get_input_motion(t_k, driver2_name, input_params)
    dq_i_k_vec = da_i_k.reshape(-1, 1); ddq_i_k_vec = dda_i_k.reshape(-1, 1)

    # Store independent variables directly
    results[driver1_name+'_m'][k]=a_i_k[0]; results['d'+driver1_name+'_m'][k]=da_i_k[0]; results['dd'+driver1_name+'_m'][k]=dda_i_k[0]
    results[driver2_name+'_m'][k]=a_i_k[1]; results['d'+driver2_name+'_m'][k]=da_i_k[1]; results['dd'+driver2_name+'_m'][k]=dda_i_k[1]

    # --- Position Analysis ---
    q_d_guess = q_all_k[dependent_idx] # Use previous step's solution
    sol, info, ier, msg = fsolve(
        position_residuals_flexible, q_d_guess,
        args=(a_i_k, param_vals_kin_num, driven_idx, dependent_idx, len(q_all_sym)),
        full_output=True, xtol=1e-8
    )

    if ier == 1:
        q_d_k = sol
        q_all_k[driven_idx] = a_i_k
        q_all_k[dependent_idx] = q_d_k
        for idx, name in enumerate(coord_names): results[name + '_m'][k] = q_all_k[idx]
    else:
        print(f"WARNING: Position solver failed at step k={k}, t={t_k:.2f}. Msg: {msg}")
        break

    results['fsolve_info'][k] = info

    # --- Velocity Analysis ---
    try:
        Cq_full_k = Cq_full_num(*q_all_k, *param_vals_kin_num)
        Cq_d_k = Cq_full_k[:, dependent_idx]
        Cq_i_k = Cq_full_k[:, driven_idx]
        cond_num = np.linalg.cond(Cq_d_k)
        if cond_num > 1e10: print(f"WARNING: Cq_d ill-conditioned at k={k}, t={t_k:.2f} (cond={cond_num:.2e})")
        rhs_vel = -Cq_i_k @ dq_i_k_vec
        dq_d_k = np.linalg.solve(Cq_d_k, rhs_vel).flatten()
        dq_all_k[driven_idx] = da_i_k
        dq_all_k[dependent_idx] = dq_d_k
        for idx, name in enumerate(coord_names): results['d' + name + '_m'][k] = dq_all_k[idx]
    except np.linalg.LinAlgError: print(f"ERROR: Velocity analysis failed at k={k}, t={t_k:.2f} (Singular Cq_d)"); break
    except Exception as e: print(f"ERROR: Velocity analysis failed at k={k}, t={t_k:.2f}. Error: {e}"); break

    # --- Acceleration Analysis ---
    try:
        gamma_full_k = gamma_full_num(*q_all_k, *dq_all_k, *param_vals_kin_num)
        rhs_accel = gamma_full_k - Cq_i_k @ ddq_i_k_vec
        ddq_d_k = np.linalg.solve(Cq_d_k, rhs_accel).flatten()
        ddq_all_k = np.zeros(len(q_all_sym))
        ddq_all_k[driven_idx] = dda_i_k
        ddq_all_k[dependent_idx] = ddq_d_k
        for idx, name in enumerate(coord_names): results['dd' + name + '_m'][k] = ddq_all_k[idx]
    except np.linalg.LinAlgError: print(f"ERROR: Acceleration analysis failed at k={k}, t={t_k:.2f} (Singular Cq_d)"); break
    except Exception as e:
        print(f"ERROR: Acceleration analysis failed at k={k}, t={t_k:.2f}. Error: {e}")
        if 'gamma_full_num' in str(e): print("   -> Likely symbolic gamma issue.")
        break

    # Progress update
    if (k + 1) % (max(1, nn // 20)) == 0: print(f"Simulation progress: {100 * (k + 1) / nn:.0f}% (t={t_k:.2f}s)")

sim_time = time.time() - sim_start_time
print(f"Numerical simulation loop finished in {sim_time:.2f} seconds.")
print("-" * 20)

# --- Post-processing and Trimming Results ---
last_valid_idx = -1
for i in range(nn - 1, -1, -1):
    dep_var_name = coord_names[dependent_idx[0]]
    if not np.isnan(results[dep_var_name + '_m'][i]): last_valid_idx = i; break
if last_valid_idx != -1:
    nn_valid = last_valid_idx + 1
    if nn_valid < nn:
        print(f"Simulation stopped early. Only {nn_valid} / {nn} steps completed.")
        t_m = t_m[:nn_valid]
        for res_key in results:
             if res_key != 'fsolve_info': results[res_key] = results[res_key][:nn_valid]
             else: results[res_key] = results[res_key][:nn_valid]
        nn = nn_valid
    else: print("All simulation steps completed successfully.")
else: print("ERROR: No steps completed successfully."); sys.exit()

# --------------------------------------------------------------------------
# 7. Animation Setup - FIVE BAR (Flexible Drivers)
# --------------------------------------------------------------------------
print("Setting up animation figure...")
# (Setup is identical to previous version)
fig_anim, ax_anim = plt.subplots()
ax_anim.set_aspect('equal')
max_reach = L1_val + L2_val + L3_val + L4_val
ax_anim.set_xlim(-0.5 * (L1_val + L2_val) - 0.2, L5_val + 0.5 * L4_val + 0.5)
ax_anim.set_ylim(-0.6 * max_reach, 0.6 * max_reach)
ax_anim.set_xlabel("X position"); ax_anim.set_ylabel("Y position")
ax_anim.set_title(f"Five-Bar Linkage Animation (Driven: {', '.join(DRIVEN_JOINTS)})")
ax_anim.grid(True)
P0_np = np.array([0, 0]); P3_ground_np = np.array([L5_val, 0])
line1, = ax_anim.plot([], [], 'o-', lw=2, color='blue', label='L1')
line2, = ax_anim.plot([], [], 'o-', lw=2, color='red', label='L2')
line3, = ax_anim.plot([], [], 'o-', lw=2, color='green', label='L3')
line4, = ax_anim.plot([], [], 'o-', lw=2, color='purple', label='L4')
line5, = ax_anim.plot([P0_np[0], P3_ground_np[0]], [P0_np[1], P3_ground_np[1]], 's-', lw=1, color='black', label=f'L5 (Gnd)')
com1, = ax_anim.plot([], [], 'bo', ms=4, label='CoM1'); com2, = ax_anim.plot([], [], 'ro', ms=4, label='CoM2')
com3, = ax_anim.plot([], [], 'go', ms=4, label='CoM3'); com4, = ax_anim.plot([], [], 'mo', ms=4, label='CoM4')
time_template = 'time = %.2fs'; time_text = ax_anim.text(0.05, 0.9, '', transform=ax_anim.transAxes)
ax_anim.legend(loc='upper right', fontsize='small')

def get_joint_positions_5bar_flexible(k):
    """Calculates joint positions for frame k (uses results dict)."""
    a1_k=results['a1_m'][k]; a2_k=results['a2_m'][k]; a3_k=results['a3_m'][k]; a4_k=results['a4_m'][k]
    r11_k=results['r11_m'][k]; r12_k=results['r12_m'][k]; r21_k=results['r21_m'][k]; r22_k=results['r22_m'][k]
    r31_k=results['r31_m'][k]; r32_k=results['r32_m'][k]; r41_k=results['r41_m'][k]; r42_k=results['r42_m'][k]
    A1_k=np.array([[np.cos(a1_k),-np.sin(a1_k)],[np.sin(a1_k),np.cos(a1_k)]])
    A2_k=np.array([[np.cos(a2_k),-np.sin(a2_k)],[np.sin(a2_k),np.cos(a2_k)]])
    A3_k=np.array([[np.cos(a3_k),-np.sin(a3_k)],[np.sin(a3_k),np.cos(a3_k)]])
    A4_k=np.array([[np.cos(a4_k),-np.sin(a4_k)],[np.sin(a4_k),np.cos(a4_k)]])
    u_1b_local=np.array([0.5*L1_val,0]); u_2c_local=np.array([0.5*L2_val,0]); u_3d_local=np.array([0.5*L3_val,0])
    R1_k=np.array([r11_k,r12_k]); R2_k=np.array([r21_k,r22_k]); R3_k=np.array([r31_k,r32_k]); R4_k=np.array([r41_k,r42_k])
    P1_k=R1_k+A1_k@u_1b_local; P2_k=R2_k+A2_k@u_2c_local; P3j_k=R3_k+A3_k@u_3d_local
    CoM1_k=R1_k; CoM2_k=R2_k; CoM3_k=R3_k; CoM4_k=R4_k
    return P1_k, P2_k, P3j_k, CoM1_k, CoM2_k, CoM3_k, CoM4_k

def init_anim():
    line1.set_data([],[]); line2.set_data([],[]); line3.set_data([],[]); line4.set_data([],[])
    com1.set_data([],[]); com2.set_data([],[]); com3.set_data([],[]); com4.set_data([],[])
    time_text.set_text('')
    return line1, line2, line3, line4, com1, com2, com3, com4, time_text

def update_anim(frame):
    if frame >= nn: return line1, line2, line3, line4, com1, com2, com3, com4, time_text
    P1_k, P2_k, P3j_k, CoM1_k, CoM2_k, CoM3_k, CoM4_k = get_joint_positions_5bar_flexible(frame)
    line1.set_data([P0_np[0], P1_k[0]], [P0_np[1], P1_k[1]])
    line2.set_data([P1_k[0], P2_k[0]], [P1_k[1], P2_k[1]])
    line3.set_data([P2_k[0], P3j_k[0]], [P2_k[1], P3j_k[1]])
    line4.set_data([P3j_k[0], P3_ground_np[0]], [P3j_k[1], P3_ground_np[1]])
    com1.set_data([CoM1_k[0]], [CoM1_k[1]]); com2.set_data([CoM2_k[0]], [CoM2_k[1]])
    com3.set_data([CoM3_k[0]], [CoM3_k[1]]); com4.set_data([CoM4_k[0]], [CoM4_k[1]])
    time_text.set_text(time_template % t_m[frame])
    return line1, line2, line3, line4, com1, com2, com3, com4, time_text

# Create animation object
ani = animation.FuncAnimation(fig_anim, update_anim, frames=nn, init_func=init_anim, blit=True, interval=dt*1000)
print("Animation object created.")
print("-" * 20)

# --------------------------------------------------------------------------
# 8. Display Static Plots and Start Animation (Revised Pause)
# --------------------------------------------------------------------------
print("Generating static plots...")
# (Plotting setup is identical to previous version)
all_angles = ['a1', 'a2', 'a3', 'a4']
dependent_angles = sorted([a for a in all_angles if a not in DRIVEN_JOINTS])
plt.figure(); plt.plot(t_m, results[driver1_name+'_m']*180/np.pi, label=driver1_name); plt.plot(t_m, results[driver2_name+'_m']*180/np.pi, label=driver2_name); plt.xlabel('Time (s)'); plt.ylabel('Angles (degrees)'); plt.title('Input Angles vs Time'); plt.legend(); plt.grid(True)
if dependent_angles: plt.figure(); plt.plot(t_m, results[dependent_angles[0]+'_m']*180/np.pi, label=dependent_angles[0]); plt.plot(t_m, results[dependent_angles[1]+'_m']*180/np.pi, label=dependent_angles[1]); plt.xlabel('Time (s)'); plt.ylabel('Angles (degrees)'); plt.title('Dependent Angles vs Time'); plt.legend(); plt.grid(True)
plt.figure(); plt.plot(results['r31_m'], results['r32_m']); plt.xlabel('r31 (m)'); plt.ylabel('r32 (m)'); plt.title('Link 3 CoM Trajectory'); plt.axis('equal'); plt.grid(True)
plt.figure(); plt.plot(results['r41_m'], results['r42_m']); plt.xlabel('r41 (m)'); plt.ylabel('r42 (m)'); plt.title('Link 4 CoM Trajectory'); plt.axis('equal'); plt.grid(True)
fig_kin1, axs_kin1 = plt.subplots(3, 1, sharex=True); fig_kin1.suptitle(f'Input {driver1_name} Kinematics'); axs_kin1[0].plot(t_m, results['dd'+driver1_name+'_m']*180/np.pi); axs_kin1[0].set_ylabel('Accel (deg/s^2)'); axs_kin1[0].grid(True); axs_kin1[1].plot(t_m, results['d'+driver1_name+'_m']*180/np.pi); axs_kin1[1].set_ylabel('Vel (deg/s)'); axs_kin1[1].grid(True); axs_kin1[2].plot(t_m, results[driver1_name+'_m']*180/np.pi); axs_kin1[2].set_ylabel('Angle (deg)'); axs_kin1[2].grid(True); axs_kin1[2].set_xlabel('Time (s)')
fig_kin2, axs_kin2 = plt.subplots(3, 1, sharex=True); fig_kin2.suptitle(f'Input {driver2_name} Kinematics'); axs_kin2[0].plot(t_m, results['dd'+driver2_name+'_m']*180/np.pi); axs_kin2[0].set_ylabel('Accel (deg/s^2)'); axs_kin2[0].grid(True); axs_kin2[1].plot(t_m, results['d'+driver2_name+'_m']*180/np.pi); axs_kin2[1].set_ylabel('Vel (deg/s)'); axs_kin2[1].grid(True); axs_kin2[2].plot(t_m, results[driver2_name+'_m']*180/np.pi); axs_kin2[2].set_ylabel('Angle (deg)'); axs_kin2[2].grid(True); axs_kin2[2].set_xlabel('Time (s)')

# --- Show static plots FIRST (blocking) ---
print("Displaying static plots... Close all plot windows to continue.")
plt.show()

# --- Wait for key press before showing animation ---
print("-" * 20)
input("Static plots closed. Press Enter to start the animation...")
print("Starting animation...")

# --- Show animation figure (blocking) ---
plt.figure(fig_anim.number)
plt.show()

print("Script finished.")