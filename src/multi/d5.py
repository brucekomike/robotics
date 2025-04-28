# %% Imports
import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.optimize import fsolve
# Corrected import: removed unused physics imports
from sympy.physics.mechanics import ReferenceFrame, Point # Keep Point for origin
import time # To time the setup

print("Starting general FIVE-bar linkage simulation setup...")
start_time = time.time()

# --------------------------------------------------------------------------
# 1. Symbolic Setup (General Case using Sympy) - FIVE BAR
# --------------------------------------------------------------------------
print("Setting up symbolic variables and equations for 5-bar linkage...")

# Time and input
t = sp.Symbol('t')
# --- Independent Coordinate ---
a1 = sp.Symbol('a1')
da1 = sp.Symbol('da1')
dda1 = sp.Symbol('dda1')

# --- Dependent Coordinates ---
# Link 1
r11, r12 = sp.symbols('r11 r12')
dr11, dr12 = sp.symbols('dr11 dr12')
ddr11, ddr12 = sp.symbols('ddr11 ddr12')
# Link 2 (r21, r22 are dependent, a2 = 0 is constraint)
r21, r22 = sp.symbols('r21 r22')
dr21, dr22 = sp.symbols('dr21 dr22')
ddr21, ddr22 = sp.symbols('ddr21 ddr22')
a2 = 0 # Apply constraint directly
da2 = 0
dda2 = 0
# Link 3
r31, r32, a3 = sp.symbols('r31 r32 a3')
dr31, dr32, da3 = sp.symbols('dr31 dr32 da3')
ddr31, ddr32, dda3 = sp.symbols('ddr31 ddr32 dda3')
# Link 4 (NEW)
r41, r42, a4 = sp.symbols('r41 r42 a4')
dr41, dr42, da4 = sp.symbols('dr41 dr42 da4')
ddr41, ddr42, dda4 = sp.symbols('ddr41 ddr42 dda4')


# Parameters
L1, L2, L3, L4, L5, g = sp.symbols('L1 L2 L3 L4 L5 g') # Added L5
# Masses/Inertia (not used for kinematics, but keep symbols for potential extension)
m1, m2, m3, m4 = sp.symbols('m1 m2 m3 m4')
I1zz, I2zz, I3zz, I4zz = sp.symbols('I1zz I2zz I3zz I4zz')

# --- Define generalized coordinates (q) and their derivatives (dq, ddq) ---
# Independent coordinate: q_i = [a1]
# Dependent coordinates: q_d = [r11, r12, r21, r22, r31, r32, a3, r41, r42, a4] (10 vars)
q_i = sp.Matrix([a1])
q_d = sp.Matrix([r11, r12, r21, r22, r31, r32, a3, r41, r42, a4])
q = sp.Matrix.vstack(q_i, q_d) # Full vector (11 vars)

dq_i = sp.Matrix([da1])
dq_d = sp.Matrix([dr11, dr12, dr21, dr22, dr31, dr32, da3, dr41, dr42, da4])
dq = sp.Matrix.vstack(dq_i, dq_d) # Full vector (11 vars)

ddq_i = sp.Matrix([dda1])
ddq_d = sp.Matrix([ddr11, ddr12, ddr21, ddr22, ddr31, ddr32, dda3, ddr41, ddr42, dda4])
ddq = sp.Matrix.vstack(ddq_i, ddq_d) # Full vector (11 vars)

# Rotation Matrices
A1 = sp.Matrix([[sp.cos(a1), -sp.sin(a1)], [sp.sin(a1), sp.cos(a1)]])
A2 = sp.eye(2) # Since a2 = 0
A3 = sp.Matrix([[sp.cos(a3), -sp.sin(a3)], [sp.sin(a3), sp.cos(a3)]])
A4 = sp.Matrix([[sp.cos(a4), -sp.sin(a4)], [sp.sin(a4), sp.cos(a4)]]) # NEW

# Position Vectors (Centers of Mass)
R1 = sp.Matrix([r11, r12])
R2 = sp.Matrix([r21, r22])
R3 = sp.Matrix([r31, r32])
R4 = sp.Matrix([r41, r42]) # NEW

# Local vectors from CoM to joints
u_1a = sp.Matrix([-L1/2, 0])
u_1b = sp.Matrix([ L1/2, 0])
u_2b = sp.Matrix([-L2/2, 0])
u_2c = sp.Matrix([ L2/2, 0])
u_3c = sp.Matrix([-L3/2, 0])
u_3d = sp.Matrix([ L3/2, 0])
u_4c = sp.Matrix([-L4/2, 0]) # NEW (Link 4 CoM to joint with Link 3)
u_4d = sp.Matrix([ L4/2, 0]) # NEW (Link 4 CoM to fixed point P3)
# L5 is the ground link length

# Fixed points
P0 = sp.Matrix([0, 0])
P3_p = sp.Matrix([L5, 0]) # Ground link length is L5

# --- Constraint Equations C(q) = 0 --- (10 equations for 10 dependent vars)
# fun1: Link 1 'a' end to P0 (2 eqns)
r1a = R1 + A1 * u_1a
fun1 = r1a - P0

# fun2: Link 1 'b' end to Link 2 'b' end (2 eqns)
r1b = R1 + A1 * u_1b
r2b = R2 + A2 * u_2b # A2 is identity
fun2 = r2b - r1b

# fun3: Link 2 'c' end to Link 3 'c' end (2 eqns)
r2c = R2 + A2 * u_2c # A2 is identity
r3c = R3 + A3 * u_3c
fun3 = r3c - r2c

# fun4: Link 3 'd' end to Link 4 'c' end (NEW, 2 eqns)
r3d = R3 + A3 * u_3d
r4c = R4 + A4 * u_4c
fun4 = r4c - r3d

# fun5: Link 4 'd' end to Ground P3_p (NEW, 2 eqns)
r4d = R4 + A4 * u_4d
fun5 = r4d - P3_p

# Full constraint vector (10 equations)
func = sp.Matrix([fun1, fun2, fun3, fun4, fun5])
print(f"Symbolic constraints defined ({func.shape[0]} equations for {q_d.shape[0]} dependent vars).")
print("\nConstraint Equations (func):")
sp.pprint(func)
print("-" * 30)

# --------------------------------------------------------------------------
# 2. Symbolic Derivations for Kinematics (Jacobians, Gamma term)
# --------------------------------------------------------------------------
print("Deriving symbolic Jacobians and acceleration term...")

# Jacobian of constraints w.r.t. all coordinates q (now 10x11)
Cq = func.jacobian(q)

# Partition Jacobian into dependent and independent parts
Cq_d = Cq[:, 1:] # Jacobian w.r.t. q_d (10x10)
Cq_i = Cq[:, 0]  # Jacobian w.r.t. q_i (a1) (10x1)

print("Jacobians Cq_d and Cq_i derived.")
print("\nJacobian Cq_d (w.r.t. dependent coordinates):")
sp.pprint(Cq_d) # This will be large
print("\nJacobian Cq_i (w.r.t. independent coordinate a1):")
sp.pprint(Cq_i)
print("-" * 30)


# Derive the acceleration term gamma = - (d(Cq)/dt * dq)
# gamma = - (jacobian(Cq * dq, q) * dq)
print("Calculating symbolic gamma term (may take significant time)...")
gamma = sp.zeros(func.shape[0], 1) # Initialize with correct size (10x1)
try:
    Cq_dq = Cq * dq
    gamma_term_jacobian = Cq_dq.jacobian(q)
    gamma = -gamma_term_jacobian * dq
    gamma = sp.simplify(gamma) # Simplify if possible (might also be slow)
    print("Symbolic gamma term calculated.")
    print("\nGamma Term (gamma):")
    sp.pprint(gamma) # This will be very large
    print("-" * 30)
except Exception as e:
    print(f"WARNING: Symbolic calculation of gamma failed: {e}")
    print("Acceleration calculation will likely fail later.")
    gamma = sp.zeros(func.shape[0],1) # Placeholder
    print("\nGamma Term (gamma): Calculation Failed, using placeholder:")
    sp.pprint(gamma)
    print("-" * 30)

# --------------------------------------------------------------------------
# 3. Inverse Dynamics Setup REMOVED
#    (Calculating Tw1 is ill-defined for this constrained 2-DOF system
#     driven by only one input 'a1')
# --------------------------------------------------------------------------
print("Inverse dynamics (KE, PE, Power, Torque) calculations are REMOVED.")
print("Focusing on kinematics for the constrained 5-bar linkage.")
print("-" * 30)


# --------------------------------------------------------------------------
# 4. Lambdify Symbolic Expressions (Kinematics Only)
# --------------------------------------------------------------------------
print("Lambdifying symbolic functions (constraints, Jacobians, gamma)...")

# Arguments needed for numerical functions
# Parameters now include L5
params_kin = [L1, L2, L3, L4, L5] # Only lengths needed for kinematics
state_q = list(q)
state_dq = list(dq)
# state_ddq = list(ddq) # Not needed for lambdify args without dynamics

# Lambdify constraint function
# Args: a1, r11, r12, r21, r22, r31, r32, a3, r41, r42, a4, L1, L2, L3, L4, L5
func_num = sp.lambdify(state_q + params_kin, func, modules=['numpy'])
print("Lambdified: func_num")

# Lambdify Jacobians
# Args: a1, r11, r12, r21, r22, r31, r32, a3, r41, r42, a4, L1, L2, L3, L4, L5
Cq_d_num = sp.lambdify(state_q + params_kin, Cq_d, modules=['numpy'])
Cq_i_num = sp.lambdify(state_q + params_kin, Cq_i, modules=['numpy'])
print("Lambdified: Cq_d_num, Cq_i_num")

# Lambdify gamma term
# Args: a1, r11..a4, da1, dr11..da4, L1..L5
gamma_num = sp.lambdify(state_q + state_dq + params_kin, gamma, modules=['numpy'])
print("Lambdified: gamma_num")

setup_time = time.time() - start_time
print(f"Symbolic setup and lambdification complete in {setup_time:.2f} seconds.")
print("-" * 20)


# --------------------------------------------------------------------------
# 5. Numerical Simulation Parameters
# --------------------------------------------------------------------------
print("Setting up simulation parameters...")
t_start = 0.0
t_end = 10.0 # Duration
dt = 0.02 # Time step
t_m = np.arange(t_start, t_end + dt, dt)
nn = len(t_m)

# Parameter Values
L1_val, L2_val, L3_val, L4_val, L5_val = 3.0, 2.0, 3.0, 4.0, 4.0 # Lengths (L5 is ground)
# Masses/Inertia (define but not used)
m1_val, m2_val, m3_val, m4_val = 1.0, 1.0, 1.0, 1.0
g_val = 10.0
I1zz_val = m1_val * L1_val**2 / 12.0
I2zz_val = m2_val * L2_val**2 / 12.0
I3zz_val = m3_val * L3_val**2 / 12.0
I4zz_val = m4_val * L4_val**2 / 12.0 # New

param_vals_kin_num = [L1_val, L2_val, L3_val, L4_val, L5_val]

w = 10.0  # Angular frequency for input motion

print(f"Time steps: {nn}")
print("-" * 20)

# --------------------------------------------------------------------------
# 6. Define Input Motion (a1 only)
# --------------------------------------------------------------------------
print("Defining input motion functions (a1)...")
def input_a1(t, w):
    return np.exp(-t) * np.cos(t * w)

def input_da1(t, w):
    exp_t = np.exp(-t)
    return -exp_t * np.cos(t * w) - w * exp_t * np.sin(t * w)

def input_dda1(t, w):
    exp_t = np.exp(-t)
    cos_tw = np.cos(t * w)
    sin_tw = np.sin(t * w)
    w_sq = w**2
    return exp_t * cos_tw - w_sq * exp_t * cos_tw + 2 * w * exp_t * sin_tw

print("Input motion functions defined.")
print("-" * 20)

# --------------------------------------------------------------------------
# 7. Run Numerical Simulation Loop (Kinematics Only)
# --------------------------------------------------------------------------
print("Running numerical simulation (kinematics only)...")

# Initialize storage arrays
results = {}
# Updated variable lists (a2 removed, a4 added)
var_names_q = ['a1', 'r11', 'r12', 'r21', 'r22', 'r31', 'r32', 'a3', 'r41', 'r42', 'a4']
var_names_dq = ['da1', 'dr11', 'dr12', 'dr21', 'dr22', 'dr31', 'dr32', 'da3', 'dr41', 'dr42', 'da4']
var_names_ddq = ['dda1', 'ddr11', 'ddr12', 'ddr21', 'ddr22', 'ddr31', 'ddr32', 'dda3', 'ddr41', 'ddr42', 'dda4']
all_var_names = var_names_q + var_names_dq + var_names_ddq

for name in all_var_names:
    results[name + '_m'] = np.full(nn, np.nan)
# results['Tw1_m'] = np.full(nn, np.nan) # REMOVED
results['fsolve_info'] = [{}] * nn # Store fsolve output dict

# --- Initial Guess (q_d) ---
# Need a guess for the 10 dependent variables
a1_0 = input_a1(t_start, w)
# Rough guess based on previous 4-bar and adding link 4 stretched out
r11_g = L1_val / 2.0 * np.cos(a1_0)
r12_g = L1_val / 2.0 * np.sin(a1_0)
a2_g = 0.0 # Fixed
r21_g = r11_g + L1_val / 2.0 * np.cos(a1_0) + L2_val / 2.0 * np.cos(a2_g)
r22_g = r12_g + L1_val / 2.0 * np.sin(a1_0) + L2_val / 2.0 * np.sin(a2_g)
a3_g = a1_0 # Guess a3 follows a1 initially
r31_g = r21_g + L2_val / 2.0 * np.cos(a2_g) + L3_val / 2.0 * np.cos(a3_g)
r32_g = r22_g + L2_val / 2.0 * np.sin(a2_g) + L3_val / 2.0 * np.sin(a3_g)
# Guess link 4 based on P3_p location (L5, 0)
a4_g = 0.0 # Guess a4 is horizontal initially
r41_g = L5_val - L4_val / 2.0 * np.cos(a4_g)
r42_g = 0.0 - L4_val / 2.0 * np.sin(a4_g)

q_d_guess = np.array([r11_g, r12_g, r21_g, r22_g, r31_g, r32_g, a3_g, r41_g, r42_g, a4_g])

print(f"Attempting to find initial configuration for a1(0) = {a1_0:.4f} using guess (10 vars): {q_d_guess}")

# Define the function to solve for q_d (needs updated signature)
def position_residuals(qd_unknown, a1_val, params_k):
    # qd_unknown has 10 elements
    # params_k has 5 elements (L1..L5)
    q_full = np.concatenate(([a1_val], qd_unknown)) # Construct full q (11 elements)
    res = func_num(*q_full, *params_k) # Pass 11 + 5 = 16 arguments
    return res.flatten() # fsolve needs a 1D array (10 residuals)

# Solve for initial q_d
sol_init, info_init, ier_init, msg_init = fsolve(position_residuals, q_d_guess, args=(a1_0, param_vals_kin_num), full_output=True)

if ier_init == 1:
    q_d_guess = sol_init # Use the found initial solution as the starting guess
    print(f"Found initial q_d: {q_d_guess}")
else:
    print(f"ERROR: Position solver (fsolve) failed for initial condition (t=0) with msg: {msg_init}")
    print("Cannot start simulation. Check initial guess or linkage parameters/constraints.")
    exit()


# --- Simulation Loop ---
sim_start_time = time.time()
for k in range(nn):
    t_k = t_m[k]
    a1_k = input_a1(t_k, w)
    da1_k = input_da1(t_k, w)
    dda1_k = input_dda1(t_k, w)

    # Store independent variables
    results['a1_m'][k] = a1_k
    results['da1_m'][k] = da1_k
    results['dda1_m'][k] = dda1_k

    # --- Position Analysis ---
    sol, info, ier, msg = fsolve(position_residuals, q_d_guess, args=(a1_k, param_vals_kin_num), full_output=True)

    if ier == 1:
        q_d_k = sol
        q_d_guess = q_d_k # Use current solution as guess for next step
        # Unpack results into dictionary (10 variables)
        results['r11_m'][k], results['r12_m'][k], results['r21_m'][k], results['r22_m'][k], \
        results['r31_m'][k], results['r32_m'][k], results['a3_m'][k], \
        results['r41_m'][k], results['r42_m'][k], results['a4_m'][k] = q_d_k
    else:
        print(f"WARNING: Position solver (fsolve) failed at step k={k}, t={t_k:.2f} with msg: {msg}")
        break # Stop simulation if position fails

    results['fsolve_info'][k] = info

    # Assemble full state vector q_k (11 elements)
    q_k = np.concatenate(([a1_k], q_d_k))

    # --- Velocity Analysis ---
    try:
        Cq_d_k = Cq_d_num(*q_k, *param_vals_kin_num) # 10x10 matrix
        Cq_i_k = Cq_i_num(*q_k, *param_vals_kin_num) # 10x1 vector

        cond_num = np.linalg.cond(Cq_d_k)
        if cond_num > 1e9: # Increased threshold slightly
             print(f"WARNING: Cq_d is ill-conditioned at step k={k}, t={t_k:.2f} (cond={cond_num:.2e}). Potential singularity.")

        # Solve Cq_d * dq_d = -Cq_i * da1
        rhs_vel = -Cq_i_k * da1_k
        dq_d_k = np.linalg.solve(Cq_d_k, rhs_vel).flatten() # Solve for 10 dependent velocities

        # Unpack results (10 variables)
        results['dr11_m'][k], results['dr12_m'][k], results['dr21_m'][k], results['dr22_m'][k], \
        results['dr31_m'][k], results['dr32_m'][k], results['da3_m'][k], \
        results['dr41_m'][k], results['dr42_m'][k], results['da4_m'][k] = dq_d_k

    except np.linalg.LinAlgError:
        print(f"ERROR: Velocity analysis failed at step k={k}, t={t_k:.2f} (Singular Cq_d matrix)")
        break
    except Exception as e:
        print(f"ERROR: Velocity analysis failed at step k={k}, t={t_k:.2f} with error: {e}")
        break

    # Assemble full velocity vector dq_k (11 elements)
    dq_k = np.concatenate(([da1_k], dq_d_k))

    # --- Acceleration Analysis ---
    try:
        # Evaluate gamma term
        gamma_k = gamma_num(*q_k, *dq_k, *param_vals_kin_num) # 10x1 vector

        # Solve Cq_d * ddq_d = gamma - Cq_i * dda1
        rhs_accel = gamma_k - Cq_i_k * dda1_k
        ddq_d_k = np.linalg.solve(Cq_d_k, rhs_accel).flatten() # Solve for 10 dependent accelerations

        # Unpack results (10 variables)
        results['ddr11_m'][k], results['ddr12_m'][k], results['ddr21_m'][k], results['ddr22_m'][k], \
        results['ddr31_m'][k], results['ddr32_m'][k], results['dda3_m'][k], \
        results['ddr41_m'][k], results['ddr42_m'][k], results['dda4_m'][k] = ddq_d_k

    except np.linalg.LinAlgError:
        print(f"ERROR: Acceleration analysis failed at step k={k}, t={t_k:.2f} (Singular Cq_d matrix)")
        break
    except Exception as e:
        print(f"ERROR: Acceleration analysis failed at step k={k}, t={t_k:.2f} with error: {e}")
        if 'gamma_num' in str(e):
             print("   -> Failure likely due to symbolic gamma calculation issue earlier.")
        break

    # Assemble full acceleration vector ddq_k (11 elements)
    # ddq_k = np.concatenate(([dda1_k], ddq_d_k)) # Not strictly needed as dynamics removed

    # --- Inverse Dynamics REMOVED ---

    # Progress update
    if (k + 1) % (max(1, nn // 20)) == 0:
        print(f"Simulation progress: {100 * (k + 1) / nn:.0f}% (t={t_k:.2f}s)")

sim_time = time.time() - sim_start_time
print(f"Numerical simulation loop finished in {sim_time:.2f} seconds.")
print("-" * 20)

# Find the last successfully computed index (check a kinematic variable like a4)
last_valid_index = np.where(~np.isnan(results['a4_m']))[0]
if len(last_valid_index) > 0:
    nn_valid = last_valid_index[-1] + 1
    if nn_valid < nn:
        print(f"Simulation stopped early. Only {nn_valid} out of {nn} steps completed successfully.")
        # Trim results to valid range
        t_m = t_m[:nn_valid]
        for name in results:
             if name != 'fsolve_info':
                 results[name] = results[name][:nn_valid]
             else:
                 results[name] = results[name][:nn_valid]
        nn = nn_valid # Update nn to valid length for plotting/animation
    else:
        print("All simulation steps completed successfully.")

else:
    print("ERROR: No steps completed successfully.")
    exit()

# --------------------------------------------------------------------------
# 8. Animation Setup (Matplotlib) - FIVE BAR
# --------------------------------------------------------------------------
print("Setting up animation figure...")

fig_anim, ax_anim = plt.subplots()
ax_anim.set_aspect('equal')
# Adjust limits based on potential motion range (rough estimate)
max_reach = L1_val + L2_val + L3_val + L4_val
ax_anim.set_xlim(-0.5 * (L1_val + L2_val), L5_val + 0.5 * L4_val + 0.5) # Adjust x limits
ax_anim.set_ylim(-0.5 * max_reach, 0.5 * max_reach)
ax_anim.set_xlabel("X position")
ax_anim.set_ylabel("Y position")
ax_anim.set_title("Five-Bar Linkage Animation (a2=0 constraint)")
ax_anim.grid(True)

# Define fixed points
P0_np = np.array([0, 0])
P3_ground_np = np.array([L5_val, 0]) # Renamed for clarity

# Initialize lines for the links
line1, = ax_anim.plot([], [], 'o-', lw=2, color='blue', label='Link 1 (L1)') # P0 to P1
line2, = ax_anim.plot([], [], 'o-', lw=2, color='red', label='Link 2 (L2)')  # P1 to P2
line3, = ax_anim.plot([], [], 'o-', lw=2, color='green', label='Link 3 (L3)') # P2 to P3 (joint)
line4, = ax_anim.plot([], [], 'o-', lw=2, color='purple', label='Link 4 (L4)') # P3 (joint) to P4 (joint) -> NEW
line5, = ax_anim.plot([P0_np[0], P3_ground_np[0]], [P0_np[1], P3_ground_np[1]], 's-', lw=1, color='black', label=f'Link 5 (Fixed L{L5})') # Ground link P0 to P3_ground

# Add CoM markers
com1, = ax_anim.plot([], [], 'bo', ms=4, label='CoM 1')
com2, = ax_anim.plot([], [], 'ro', ms=4, label='CoM 2')
com3, = ax_anim.plot([], [], 'go', ms=4, label='CoM 3')
com4, = ax_anim.plot([], [], 'mo', ms=4, label='CoM 4') # New marker (magenta)


# Add a time text
time_template = 'time = %.2fs'
time_text = ax_anim.text(0.05, 0.9, '', transform=ax_anim.transAxes)

ax_anim.legend(loc='upper right')

# Function to calculate joint positions P1, P2, P3j (joint 3), P4j (joint 4) - FIVE BAR
def get_joint_positions_5bar(k):
    """Calculates joint positions for frame k using general coords."""
    # Get CoM and angles for this frame
    r11_k = results['r11_m'][k]
    r12_k = results['r12_m'][k]
    a1_k = results['a1_m'][k]
    r21_k = results['r21_m'][k]
    r22_k = results['r22_m'][k]
    a2_k = 0.0 # Constraint
    r31_k = results['r31_m'][k]
    r32_k = results['r32_m'][k]
    a3_k = results['a3_m'][k]
    r41_k = results['r41_m'][k] # NEW
    r42_k = results['r42_m'][k] # NEW
    a4_k = results['a4_m'][k] # NEW

    # Calculate rotation matrices numerically
    A1_k = np.array([[np.cos(a1_k), -np.sin(a1_k)], [np.sin(a1_k), np.cos(a1_k)]])
    A2_k = np.identity(2) # Since a2=0
    A3_k = np.array([[np.cos(a3_k), -np.sin(a3_k)], [np.sin(a3_k), np.cos(a3_k)]])
    A4_k = np.array([[np.cos(a4_k), -np.sin(a4_k)], [np.sin(a4_k), np.cos(a4_k)]]) # NEW

    # Define local vectors to joints relative to CoM (as numpy arrays)
    u_1b_local = np.array([0.5 * L1_val, 0]) # P1 relative to R1
    u_2c_local = np.array([0.5 * L2_val, 0]) # P2 relative to R2
    u_3d_local = np.array([0.5 * L3_val, 0]) # P3j relative to R3
    u_4c_local = np.array([-0.5 * L4_val, 0]) # P3j relative to R4 (matches fun4)

    # Calculate global CoM positions
    R1_k = np.array([r11_k, r12_k])
    R2_k = np.array([r21_k, r22_k])
    R3_k = np.array([r31_k, r32_k])
    R4_k = np.array([r41_k, r42_k]) # NEW

    # Calculate global joint positions
    P1_k = R1_k + A1_k @ u_1b_local # Joint between L1 and L2
    P2_k = R2_k + A2_k @ u_2c_local # Joint between L2 and L3
    P3j_k = R3_k + A3_k @ u_3d_local # Joint between L3 and L4 (j for joint)
    # P4 can be calculated from R4 as well: P3j_k_alt = R4_k + A4_k @ u_4c_local (should be same)

    # Also return CoM positions
    CoM1_k = R1_k
    CoM2_k = R2_k
    CoM3_k = R3_k
    CoM4_k = R4_k # NEW

    return P1_k, P2_k, P3j_k, CoM1_k, CoM2_k, CoM3_k, CoM4_k

# Initialization function for the animation
def init():
    line1.set_data([], [])
    line2.set_data([], [])
    line3.set_data([], [])
    line4.set_data([], []) # Init line4
    com1.set_data([], [])
    com2.set_data([], [])
    com3.set_data([], [])
    com4.set_data([], []) # Init com4
    time_text.set_text('')
    # Return all animated artists
    return line1, line2, line3, line4, com1, com2, com3, com4, time_text

# Update function for the animation
def update(frame):
    if frame >= nn:
        return line1, line2, line3, line4, com1, com2, com3, com4, time_text

    P1_k, P2_k, P3j_k, CoM1_k, CoM2_k, CoM3_k, CoM4_k = get_joint_positions_5bar(frame)

    line1.set_data([P0_np[0], P1_k[0]], [P0_np[1], P1_k[1]])
    line2.set_data([P1_k[0], P2_k[0]], [P1_k[1], P2_k[1]])
    line3.set_data([P2_k[0], P3j_k[0]], [P2_k[1], P3j_k[1]])
    line4.set_data([P3j_k[0], P3_ground_np[0]], [P3j_k[1], P3_ground_np[1]]) # Link 4 connects P3j to P3_ground
    com1.set_data([CoM1_k[0]], [CoM1_k[1]])
    com2.set_data([CoM2_k[0]], [CoM2_k[1]])
    com3.set_data([CoM3_k[0]], [CoM3_k[1]])
    com4.set_data([CoM4_k[0]], [CoM4_k[1]]) # Update com4
    time_text.set_text(time_template % t_m[frame])
    # Return all animated artists
    return line1, line2, line3, line4, com1, com2, com3, com4, time_text

# Create animation object
ani = animation.FuncAnimation(fig_anim, update, frames=nn,
                              init_func=init, blit=True, interval=dt*1000)

print("Animation object created.")
print("-" * 20)


# --------------------------------------------------------------------------
# 9. Display Static Plots (Kinematics Only) and Start Animation
# --------------------------------------------------------------------------
print("Generating static plots...")

# --- REMOVED Torque Plots ---

# Angle a1 vs Time
plt.figure()
plt.plot(t_m, results['a1_m'] * 180 / np.pi)
plt.xlabel('Time (s)')
plt.ylabel('Angle a1 (degrees)')
plt.title('Input Angle vs Time')
plt.grid(True)

# Angles a3, a4 vs Time (a2 is fixed at 0)
plt.figure()
plt.plot(t_m, results['a3_m'] * 180 / np.pi, label='a3')
plt.plot(t_m, results['a4_m'] * 180 / np.pi, label='a4')
plt.xlabel('Time (s)')
plt.ylabel('Angles (degrees)')
plt.title('Dependent Angles vs Time (Five-Bar, a2=0)')
plt.legend()
plt.grid(True)

# Link 3 CoM Position, Velocity, Acceleration (X-component) vs Time
plt.figure()
plt.plot(t_m, results['r31_m'], label='r31 (Pos)')
plt.plot(t_m, results['dr31_m'], label='dr31 (Vel)')
plt.plot(t_m, results['ddr31_m'], label='ddr31 (Acc)')
plt.xlabel('Time (s)')
plt.ylabel('Link 3 CoM X-Kinematics')
plt.title('Link 3 CoM X-Component Kinematics (Five-Bar)')
plt.legend()
plt.grid(True)

# Link 4 CoM Position, Velocity, Acceleration (X-component) vs Time (NEW)
plt.figure()
plt.plot(t_m, results['r41_m'], label='r41 (Pos)')
plt.plot(t_m, results['dr41_m'], label='dr41 (Vel)')
plt.plot(t_m, results['ddr41_m'], label='ddr41 (Acc)')
plt.xlabel('Time (s)')
plt.ylabel('Link 4 CoM X-Kinematics')
plt.title('Link 4 CoM X-Component Kinematics (Five-Bar)')
plt.legend()
plt.grid(True)


# Link 3 CoM Trajectory
plt.figure()
plt.plot(results['r31_m'], results['r32_m'])
plt.xlabel('r31 (m)')
plt.ylabel('r32 (m)')
plt.title('Link 3 CoM Trajectory (Five-Bar)')
plt.axis('equal')
plt.grid(True)

# Link 4 CoM Trajectory (NEW)
plt.figure()
plt.plot(results['r41_m'], results['r42_m'])
plt.xlabel('r41 (m)')
plt.ylabel('r42 (m)')
plt.title('Link 4 CoM Trajectory (Five-Bar)')
plt.axis('equal')
plt.grid(True)


# Input Kinematics vs Time (Subplots) - Same as before
fig_kin, axs_kin = plt.subplots(3, 1, sharex=True)
fig_kin.suptitle('Input Kinematics vs Time')
axs_kin[0].plot(t_m, results['dda1_m'] * 180 / np.pi)
axs_kin[0].set_ylabel('dda1 (deg/s^2)')
axs_kin[0].grid(True)
axs_kin[1].plot(t_m, results['da1_m'] * 180 / np.pi)
axs_kin[1].set_ylabel('da1 (deg/s)')
axs_kin[1].grid(True)
axs_kin[2].plot(t_m, results['a1_m'] * 180 / np.pi)
axs_kin[2].set_ylabel('a1 (deg)')
axs_kin[2].set_xlabel('Time (s)')
axs_kin[2].grid(True)

# --- Show static plots non-blockingly ---
print("Displaying static plots...")
plt.show(block=False)

# --- Wait for key press before showing animation ---
print("-" * 20)
input("Press Enter to start the animation...")
print("Starting animation...")

# --- Show animation figure (blocking) ---
plt.show() # This will now show the fig_anim and run the animation

print("Script finished.")