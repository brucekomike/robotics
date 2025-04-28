# %% Imports
import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.optimize import fsolve
# Corrected import: removed angular_velocity, inertia_of_point_mass
from sympy.physics.mechanics import (ReferenceFrame, Point, RigidBody, inertia,
                                     dynamicsymbols, KanesMethod, Particle)
import time # To time the setup

print("Starting general four-bar linkage simulation setup...")
start_time = time.time()

# --------------------------------------------------------------------------
# 1. Symbolic Setup (General Case using Sympy)
# --------------------------------------------------------------------------
print("Setting up symbolic variables and equations...")

# Time and input
t = sp.Symbol('t')
# Using dynamicsymbols for time-varying quantities is standard in sympy.physics
# However, for direct substitution later, regular symbols might be easier if not using KanesMethod EOM directly for simulation
# Let's use regular symbols matching the original approach more closely for kinematics first
a1, r11, r12, r21, r22, a2, r31, r32, a3 = sp.symbols('a1 r11 r12 r21 r22 a2 r31 r32 a3')
da1, dr11, dr12, dr21, dr22, da2, dr31, dr32, da3 = sp.symbols('da1 dr11 dr12 dr21 dr22 da2 dr31 dr32 da3')
dda1, ddr11, ddr12, ddr21, ddr22, dda2, ddr31, ddr32, dda3 = sp.symbols('dda1 ddr11 ddr12 ddr21 ddr22 dda2 ddr31 ddr32 dda3')

# Parameters
L1, L2, L3, L4, g = sp.symbols('L1 L2 L3 L4 g')
m1, m2, m3 = sp.symbols('m1 m2 m3')
I1zz, I2zz, I3zz = sp.symbols('I1zz I2zz I3zz') # Moment of inertia about CoM Z-axis

# Define generalized coordinates (q) and their derivatives (dq, ddq)
# Independent coordinate: q_i = [a1]
# Dependent coordinates: q_d = [r11, r12, r21, r22, a2, r31, r32, a3]
q = sp.Matrix([a1, r11, r12, r21, r22, a2, r31, r32, a3]) # Full vector for easier jacobian indexing
q_i = sp.Matrix([a1])
q_d = sp.Matrix([r11, r12, r21, r22, a2, r31, r32, a3])

dq = sp.Matrix([da1, dr11, dr12, dr21, dr22, da2, dr31, dr32, da3])
dq_i = sp.Matrix([da1])
dq_d = sp.Matrix([dr11, dr12, dr21, dr22, da2, dr31, dr32, da3])

ddq = sp.Matrix([dda1, ddr11, ddr12, ddr21, ddr22, dda2, ddr31, ddr32, dda3])
ddq_i = sp.Matrix([dda1])
ddq_d = sp.Matrix([ddr11, ddr12, ddr21, ddr22, dda2, ddr31, ddr32, dda3])

# Rotation Matrices
A1 = sp.Matrix([[sp.cos(a1), -sp.sin(a1)], [sp.sin(a1), sp.cos(a1)]])
A2 = sp.Matrix([[sp.cos(a2), -sp.sin(a2)], [sp.sin(a2), sp.cos(a2)]])
A3 = sp.Matrix([[sp.cos(a3), -sp.sin(a3)], [sp.sin(a3), sp.cos(a3)]])

# Position Vectors (Centers of Mass)
R1 = sp.Matrix([r11, r12])
R2 = sp.Matrix([r21, r22])
R3 = sp.Matrix([r31, r32])

# Local vectors from CoM to joints
u_1a = sp.Matrix([-L1/2, 0])
u_1b = sp.Matrix([ L1/2, 0])
u_2b = sp.Matrix([-L2/2, 0])
u_2c = sp.Matrix([ L2/2, 0])
u_3c = sp.Matrix([-L3/2, 0])
u_3d = sp.Matrix([ L3/2, 0]) # Corrected sign from original MATLAB u_3d

# Fixed points
P0 = sp.Matrix([0, 0])
P3_p = sp.Matrix([L4, 0]) # Use P3_p to avoid conflict with Body 3

# Constraint Equations C(q) = 0
# fun1: Link 1 connected to P0 at its 'a' end
r1a = R1 + A1 * u_1a
fun1 = r1a - P0

# fun2: Link 1 'b' end connected to Link 2 'b' end
r1b = R1 + A1 * u_1b
r2b = R2 + A2 * u_2b
fun2 = r2b - r1b

# fun3: Link 2 'c' end connected to Link 3 'c' end
r2c = R2 + A2 * u_2c
r3c = R3 + A3 * u_3c
fun3 = r3c - r2c

# fun4: Link 3 'd' end connected to P3_p at its 'd' end
r3d = R3 + A3 * u_3d
fun4 = r3d - P3_p

# Full constraint vector (8 equations)
func = sp.Matrix([fun1, fun2, fun3, fun4])
print(f"Symbolic constraints defined ({func.shape[0]} equations).")
print("\nConstraint Equations (func):")
sp.pprint(func)
print("-" * 30)

# --------------------------------------------------------------------------
# 2. Symbolic Derivations for Kinematics (Jacobians, Gamma term)
# --------------------------------------------------------------------------
print("Deriving symbolic Jacobians and acceleration term...")

# Jacobian of constraints w.r.t. all coordinates q
Cq = func.jacobian(q)

# Partition Jacobian into dependent and independent parts
Cq_d = Cq[:, 1:] # Jacobian w.r.t. q_d (8x8)
Cq_i = Cq[:, 0]  # Jacobian w.r.t. q_i (a1) (8x1)

print("Jacobians Cq_d and Cq_i derived.")
print("\nJacobian Cq_d (w.r.t. dependent coordinates):")
sp.pprint(Cq_d)
print("\nJacobian Cq_i (w.r.t. independent coordinate a1):")
sp.pprint(Cq_i)
print("-" * 30)


# Derive the acceleration term gamma = - (d(Cq)/dt * dq)
# gamma = - (jacobian(Cq * dq, q) * dq)
print("Calculating symbolic gamma term (may take time)...")
gamma = sp.zeros(8, 1)
try:
    Cq_dq = Cq * dq
    gamma_term_jacobian = Cq_dq.jacobian(q)
    gamma = -gamma_term_jacobian * dq
    gamma = sp.simplify(gamma) # Simplify if possible (might also be slow)
    print("Symbolic gamma term calculated.")
    print("\nGamma Term (gamma):")
    sp.pprint(gamma)
    print("-" * 30)
except Exception as e:
    print(f"WARNING: Symbolic calculation of gamma failed: {e}")
    print("Acceleration calculation will likely fail later.")
    gamma = sp.zeros(8,1) # Placeholder
    print("\nGamma Term (gamma): Calculation Failed, using placeholder:")
    sp.pprint(gamma)
    print("-" * 30)


# --------------------------------------------------------------------------
# 3. Symbolic Setup for Inverse Dynamics (using sympy.physics.mechanics)
# --------------------------------------------------------------------------
print("Setting up symbolic inverse dynamics using sympy.physics.mechanics...")
# (Setup for frames, points, bodies remains the same as before)
# Define Reference Frames
N = ReferenceFrame('N') # Inertial frame
A1_f = N.orientnew('A1_f', 'Axis', [a1, N.z])
A2_f = N.orientnew('A2_f', 'Axis', [a2, N.z])
A3_f = N.orientnew('A3_f', 'Axis', [a3, N.z])

# Define Points
P0_pt = Point('P0_pt')
P0_pt.set_vel(N, 0) # Fixed origin

R1_pt = P0_pt.locatenew('R1_pt', r11*N.x + r12*N.y)
R2_pt = P0_pt.locatenew('R2_pt', r21*N.x + r22*N.y)
R3_pt = P0_pt.locatenew('R3_pt', r31*N.x + r32*N.y)

# Set velocities of CoMs based on symbolic derivatives
R1_pt.set_vel(N, dr11*N.x + dr12*N.y)
R2_pt.set_vel(N, dr21*N.x + dr22*N.y)
R3_pt.set_vel(N, dr31*N.x + dr32*N.y)

# Set angular velocities
A1_f.set_ang_vel(N, da1 * N.z)
A2_f.set_ang_vel(N, da2 * N.z)
A3_f.set_ang_vel(N, da3 * N.z)

# Define Inertia Dyadics (using Izz symbols)
I1zz_dyadic = inertia(A1_f, 0, 0, I1zz)
I2zz_dyadic = inertia(A2_f, 0, 0, I2zz)
I3zz_dyadic = inertia(A3_f, 0, 0, I3zz)


# Define Rigid Bodies
Link1 = RigidBody('Link1', R1_pt, A1_f, m1, (I1zz_dyadic, R1_pt))
Link2 = RigidBody('Link2', R2_pt, A2_f, m2, (I2zz_dyadic, R2_pt))
Link3 = RigidBody('Link3', R3_pt, A3_f, m3, (I3zz_dyadic, R3_pt))

# --- Energy / Power Method ---
# KE = 0.5 * m1 * (dr11**2 + dr12**2) + 0.5 * I1zz * da1**2 + ...
# PE = m1 * g * r12 + m2 * g * r22 + m3 * g * r32

KE = sp.Rational(1,2) * (m1*(dr11**2 + dr12**2) + I1zz*da1**2 +
                         m2*(dr21**2 + dr22**2) + I2zz*da2**2 +
                         m3*(dr31**2 + dr32**2) + I3zz*da3**2)

PE = g * (m1*r12 + m2*r22 + m3*r32)

# Calculate time derivatives using chain rule
# dKE/dt = jacobian(KE, q)*dq + jacobian(KE, dq)*ddq
# dPE/dt = jacobian(PE, q)*dq

print("Calculating symbolic dKE/dt and dPE/dt...")
vars_energy = list(q) + list(dq)
KE_grad_q = sp.Matrix([KE]).jacobian(q)
KE_grad_dq = sp.Matrix([KE]).jacobian(dq)
PE_grad_q = sp.Matrix([PE]).jacobian(q)

dKE_dt = (KE_grad_q * dq + KE_grad_dq * ddq)[0] # Scalar result
dPE_dt = (PE_grad_q * dq)[0] # Scalar result

# Tw1 * da1 = dKE_dt + dPE_dt
# P_total = dKE_dt + dPE_dt
P_total_expr = sp.simplify(dKE_dt + dPE_dt)

# Tw1 = P_total / da1 (Symbolic expression for torque)
Tw1_expr = sp.simplify(P_total_expr / da1)


print("Symbolic power and torque expressions derived.")
print("\nKinetic Energy (KE):")
sp.pprint(KE)
print("\nPotential Energy (PE):")
sp.pprint(PE)
print("\nTime Derivative of Kinetic Energy (dKE/dt):")
sp.pprint(dKE_dt)
print("\nTime Derivative of Potential Energy (dPE/dt):")
sp.pprint(dPE_dt)
print("\nTotal Power (P_total = dKE/dt + dPE/dt):")
sp.pprint(P_total_expr)
print("\nDriving Torque (Tw1 = P_total / da1):")
sp.pprint(Tw1_expr)
print("-" * 30)


# --------------------------------------------------------------------------
# 4. Lambdify Symbolic Expressions
# --------------------------------------------------------------------------
print("Lambdifying symbolic functions (constraints, Jacobians, gamma, power)...")

# Arguments needed for numerical functions
params = [L1, L2, L3, L4, m1, m2, m3, g, I1zz, I2zz, I3zz]
state_q = list(q)
state_dq = list(dq)
state_ddq = list(ddq)

# Lambdify constraint function
func_num = sp.lambdify(state_q + [L1, L2, L3, L4], func, modules=['numpy'])
print("Lambdified: func_num")

# Lambdify Jacobians
Cq_d_num = sp.lambdify(state_q + [L1, L2, L3, L4], Cq_d, modules=['numpy'])
Cq_i_num = sp.lambdify(state_q + [L1, L2, L3, L4], Cq_i, modules=['numpy'])
print("Lambdified: Cq_d_num, Cq_i_num")

# Lambdify gamma term
gamma_num = sp.lambdify(state_q + state_dq + [L1, L2, L3, L4], gamma, modules=['numpy'])
print("Lambdified: gamma_num")

# Lambdify Power term (using P_total_expr)
P_total_num = sp.lambdify(state_q + state_dq + state_ddq + params, P_total_expr, modules=['numpy'])
print("Lambdified: P_total_num")

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
#L1_val, L2_val, L3_val, L4_val = 1.3, 1.0, 1.2, 1.2
L1_val, L2_val, L3_val, L4_val = 1.3, 1.2, 1.2, 1.0

m1_val, m2_val, m3_val = 1.0, 1.0, 1.0
g_val = 10.0
# Inertias (slender rods)
I1zz_val = m1_val * L1_val**2 / 12.0
I2zz_val = m2_val * L2_val**2 / 12.0
I3zz_val = m3_val * L3_val**2 / 12.0

param_vals = [L1_val, L2_val, L3_val, L4_val, m1_val, m2_val, m3_val, g_val, I1zz_val, I2zz_val, I3zz_val]
param_vals_kin = [L1_val, L2_val, L3_val, L4_val] # For kinematics only funcs

w = 10.0  # Angular frequency for input motion

print(f"Time steps: {nn}")
print("-" * 20)

# --------------------------------------------------------------------------
# 6. Define Input Motion
# --------------------------------------------------------------------------
print("Defining input motion functions...")
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
# 7. Run Numerical Simulation Loop
# --------------------------------------------------------------------------
print("Running numerical simulation...")

# Initialize storage arrays
results = {}
var_names_q = ['a1', 'r11', 'r12', 'r21', 'r22', 'a2', 'r31', 'r32', 'a3']
var_names_dq = ['da1', 'dr11', 'dr12', 'dr21', 'dr22', 'da2', 'dr31', 'dr32', 'da3']
var_names_ddq = ['dda1', 'ddr11', 'ddr12', 'ddr21', 'ddr22', 'dda2', 'ddr31', 'ddr32', 'dda3']
all_var_names = var_names_q + var_names_dq + var_names_ddq

for name in all_var_names:
    results[name + '_m'] = np.full(nn, np.nan)
results['Tw1_m'] = np.full(nn, np.nan)
results['fsolve_info'] = [{}] * nn # Store fsolve output dict

# --- Initial Guess (q_d) ---
a1_0 = input_a1(t_start, w)
# Use a configuration near a1=0 as a rough guess for fsolve at t=0
r11_g = L1_val / 2.0 * np.cos(a1_0) # Slightly better guess using a1_0
r12_g = L1_val / 2.0 * np.sin(a1_0)
a2_g = 0.0 # Assume a2 starts near 0
r21_g = r11_g + L1_val / 2.0 * np.cos(a1_0) + L2_val / 2.0 * np.cos(a2_g)
r22_g = r12_g + L1_val / 2.0 * np.sin(a1_0) + L2_val / 2.0 * np.sin(a2_g)
a3_g = a1_0 # Assume a3 starts near a1
r31_g = L4_val + L3_val / 2.0 * np.cos(a3_g) # Guess based on fixed end
r32_g = L3_val / 2.0 * np.sin(a3_g)
q_d_guess = np.array([r11_g, r12_g, r21_g, r22_g, a2_g, r31_g, r32_g, a3_g])

print(f"Attempting to find initial configuration for a1(0) = {a1_0:.4f} using guess: {q_d_guess}")

# Define the function to solve for q_d (needed outside loop for initial step)
def position_residuals(qd_unknown, a1_val, params_k):
    q_full = np.concatenate(([a1_val], qd_unknown))
    res = func_num(*q_full, *params_k)
    return res.flatten() # fsolve needs a 1D array

# Solve for initial q_d
sol_init, info_init, ier_init, msg_init = fsolve(position_residuals, q_d_guess, args=(a1_0, param_vals_kin), full_output=True)

if ier_init == 1:
    q_d_guess = sol_init # Use the found initial solution as the starting guess
    print(f"Found initial q_d: {q_d_guess}")
else:
    print(f"ERROR: Position solver (fsolve) failed for initial condition (t=0) with msg: {msg_init}")
    print("Cannot start simulation. Check initial guess or linkage parameters.")
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
    # Solve C(q_d, a1_k) = 0 for q_d using previous step's solution as guess
    sol, info, ier, msg = fsolve(position_residuals, q_d_guess, args=(a1_k, param_vals_kin), full_output=True)

    if ier == 1:
        q_d_k = sol
        q_d_guess = q_d_k # Use current solution as guess for next step
        results['r11_m'][k], results['r12_m'][k], results['r21_m'][k], results['r22_m'][k], \
        results['a2_m'][k], results['r31_m'][k], results['r32_m'][k], results['a3_m'][k] = q_d_k
    else:
        print(f"WARNING: Position solver (fsolve) failed at step k={k}, t={t_k:.2f} with msg: {msg}")
        break # Stop simulation if position fails

    results['fsolve_info'][k] = info # Store fsolve info dictionary

    # Assemble full state vector q_k
    q_k = np.concatenate(([a1_k], q_d_k))

    # --- Velocity Analysis ---
    try:
        Cq_d_k = Cq_d_num(*q_k, *param_vals_kin)
        Cq_i_k = Cq_i_num(*q_k, *param_vals_kin)

        # Check condition number of Cq_d_k
        cond_num = np.linalg.cond(Cq_d_k)
        if cond_num > 1e8: # Threshold for ill-conditioning (adjust as needed)
             print(f"WARNING: Cq_d is ill-conditioned at step k={k}, t={t_k:.2f} (cond={cond_num:.2e}). Potential singularity.")

        # Solve Cq_d * dq_d = -Cq_i * da1
        rhs_vel = -Cq_i_k * da1_k
        dq_d_k = np.linalg.solve(Cq_d_k, rhs_vel).flatten()

        results['dr11_m'][k], results['dr12_m'][k], results['dr21_m'][k], results['dr22_m'][k], \
        results['da2_m'][k], results['dr31_m'][k], results['dr32_m'][k], results['da3_m'][k] = dq_d_k

    except np.linalg.LinAlgError:
        print(f"ERROR: Velocity analysis failed at step k={k}, t={t_k:.2f} (Singular Cq_d matrix)")
        break
    except Exception as e:
        print(f"ERROR: Velocity analysis failed at step k={k}, t={t_k:.2f} with error: {e}")
        break

    # Assemble full velocity vector dq_k
    dq_k = np.concatenate(([da1_k], dq_d_k))

    # --- Acceleration Analysis ---
    try:
        # Evaluate gamma term
        gamma_k = gamma_num(*q_k, *dq_k, *param_vals_kin)

        # Solve Cq_d * ddq_d = gamma - Cq_i * dda1
        rhs_accel = gamma_k - Cq_i_k * dda1_k
        ddq_d_k = np.linalg.solve(Cq_d_k, rhs_accel).flatten()

        results['ddr11_m'][k], results['ddr12_m'][k], results['ddr21_m'][k], results['ddr22_m'][k], \
        results['dda2_m'][k], results['ddr31_m'][k], results['ddr32_m'][k], results['dda3_m'][k] = ddq_d_k

    except np.linalg.LinAlgError:
        print(f"ERROR: Acceleration analysis failed at step k={k}, t={t_k:.2f} (Singular Cq_d matrix)")
        break
    except Exception as e:
        print(f"ERROR: Acceleration analysis failed at step k={k}, t={t_k:.2f} with error: {e}")
        if 'gamma_num' in str(e):
             print("   -> Failure likely due to symbolic gamma calculation issue earlier.")
        break

    # Assemble full acceleration vector ddq_k
    ddq_k = np.concatenate(([dda1_k], ddq_d_k))

    # --- Inverse Dynamics (Torque Calculation) ---
    try:
        # Calculate Total Power = dKE/dt + dPE/dt
        P_total_k = P_total_num(*q_k, *dq_k, *ddq_k, *param_vals)

        # Tw1 * da1 = P_total
        if abs(da1_k) > 1e-9: # Avoid division by zero
            Tw1_k = P_total_k / da1_k
        elif abs(P_total_k) < 1e-9: # If da1 is zero, power must also be zero unless infinite torque
             Tw1_k = 0.0 # Or could be indeterminate/arbitrary if system is at rest
        else:
             print(f"WARNING: da1 is near zero ({da1_k:.2e}) but P_total is not ({P_total_k:.2e}) at step {k}. Setting Tw1 to NaN.")
             Tw1_k = np.nan # Indicates potential issue or singularity

        results['Tw1_m'][k] = Tw1_k

    except Exception as e:
        print(f"ERROR: Torque calculation failed at step k={k}, t={t_k:.2f} with error: {e}")
        break

    # Progress update
    if (k + 1) % (max(1, nn // 20)) == 0: # Update progress more often
        print(f"Simulation progress: {100 * (k + 1) / nn:.0f}% (t={t_k:.2f}s)")

sim_time = time.time() - sim_start_time
print(f"Numerical simulation loop finished in {sim_time:.2f} seconds.")
print("-" * 20)

# Find the last successfully computed index
last_valid_index = np.where(~np.isnan(results['Tw1_m']))[0] # Check torque as it's the last step
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
# 8. Animation Setup (Matplotlib) - GENERALIZED
# --------------------------------------------------------------------------
print("Setting up animation figure...")

fig_anim, ax_anim = plt.subplots() # Use a different figure name
ax_anim.set_aspect('equal')
max_reach = L1_val + L2_val + L3_val
ax_anim.set_xlim(-0.5 * L1_val - 0.5*L2_val, L4_val + 1.0 * L3_val + 0.2)
ax_anim.set_ylim(-0.5 * max_reach, 0.5 * max_reach)
ax_anim.set_xlabel("X position")
ax_anim.set_ylabel("Y position")
ax_anim.set_title("Four-Bar Linkage Animation (General Numerical Solution)")
ax_anim.grid(True)

# Define fixed points
P0_np = np.array([0, 0])
P3_np = np.array([L4_val, 0])

# Initialize lines for the links
line1, = ax_anim.plot([], [], 'o-', lw=2, color='blue', label='Link 1 (L1)')
line2, = ax_anim.plot([], [], 'o-', lw=2, color='red', label='Link 2 (L2)')
line3, = ax_anim.plot([], [], 'o-', lw=2, color='green', label='Link 3 (L3)')
line4, = ax_anim.plot([P0_np[0], P3_np[0]], [P0_np[1], P3_np[1]], 's-', lw=1, color='black', label='Link 4 (Fixed L4)')
# Add CoM markers
com1, = ax_anim.plot([], [], 'bo', ms=4, label='CoM 1')
com2, = ax_anim.plot([], [], 'ro', ms=4, label='CoM 2')
com3, = ax_anim.plot([], [], 'go', ms=4, label='CoM 3')

# Add a time text
time_template = 'time = %.2fs'
time_text = ax_anim.text(0.05, 0.9, '', transform=ax_anim.transAxes)

ax_anim.legend(loc='upper right')

# Function to calculate joint positions P1 and P2 - GENERALIZED
def get_joint_positions_general(k):
    """Calculates joint positions P1 and P2 for frame k using general coords."""
    r11_k = results['r11_m'][k]
    r12_k = results['r12_m'][k]
    a1_k = results['a1_m'][k]
    r21_k = results['r21_m'][k]
    r22_k = results['r22_m'][k]
    a2_k = results['a2_m'][k]
    r31_k = results['r31_m'][k]
    r32_k = results['r32_m'][k]
    a3_k = results['a3_m'][k]

    A1_k = np.array([[np.cos(a1_k), -np.sin(a1_k)], [np.sin(a1_k), np.cos(a1_k)]])
    A2_k = np.array([[np.cos(a2_k), -np.sin(a2_k)], [np.sin(a2_k), np.cos(a2_k)]])

    u_1b_local = np.array([0.5 * L1_val, 0])
    u_2c_local = np.array([0.5 * L2_val, 0])

    R1_k = np.array([r11_k, r12_k])
    R2_k = np.array([r21_k, r22_k])

    P1_k = R1_k + A1_k @ u_1b_local
    P2_k = R2_k + A2_k @ u_2c_local

    CoM1_k = R1_k
    CoM2_k = R2_k
    CoM3_k = np.array([r31_k, r32_k])

    return P1_k, P2_k, CoM1_k, CoM2_k, CoM3_k

# Initialization function for the animation
def init():
    line1.set_data([], [])
    line2.set_data([], [])
    line3.set_data([], [])
    com1.set_data([], []) # Initialize CoM markers
    com2.set_data([], [])
    com3.set_data([], [])
    time_text.set_text('')
    return line1, line2, line3, com1, com2, com3, time_text

# Update function for the animation - CORRECTED CoM set_data
def update(frame):
    if frame >= nn:
        return line1, line2, line3, com1, com2, com3, time_text

    P1_k, P2_k, CoM1_k, CoM2_k, CoM3_k = get_joint_positions_general(frame)

    line1.set_data([P0_np[0], P1_k[0]], [P0_np[1], P1_k[1]])
    line2.set_data([P1_k[0], P2_k[0]], [P1_k[1], P2_k[1]])
    line3.set_data([P2_k[0], P3_np[0]], [P2_k[1], P3_np[1]])
    # --- FIX: Wrap CoM coordinates in lists ---
    com1.set_data([CoM1_k[0]], [CoM1_k[1]])
    com2.set_data([CoM2_k[0]], [CoM2_k[1]])
    com3.set_data([CoM3_k[0]], [CoM3_k[1]])
    # --- End Fix ---
    time_text.set_text(time_template % t_m[frame])
    return line1, line2, line3, com1, com2, com3, time_text

# Create animation object (but don't show it yet)
ani = animation.FuncAnimation(fig_anim, update, frames=nn,
                              init_func=init, blit=True, interval=dt*1000)

print("Animation object created.")
print("-" * 20)


# --------------------------------------------------------------------------
# 9. Display Static Plots and Start Animation
# --------------------------------------------------------------------------
print("Generating static plots...")

# Torque vs Time
plt.figure()
plt.plot(t_m, results['Tw1_m'])
plt.xlabel('Time (s)')
plt.ylabel('Driving Torque Tw1 (Nm)')
plt.title('Driving Torque vs Time (General Solution)')
plt.grid(True)

# Angle a1 vs Time
plt.figure()
plt.plot(t_m, results['a1_m'] * 180 / np.pi)
plt.xlabel('Time (s)')
plt.ylabel('Angle a1 (degrees)')
plt.title('Input Angle vs Time')
plt.grid(True)

# Angles a2, a3 vs Time
plt.figure()
plt.plot(t_m, results['a2_m'] * 180 / np.pi, label='a2')
plt.plot(t_m, results['a3_m'] * 180 / np.pi, label='a3')
plt.xlabel('Time (s)')
plt.ylabel('Angles (degrees)')
plt.title('Dependent Angles vs Time (General Solution)')
plt.legend()
plt.grid(True)

# Torque vs Angle
plt.figure()
plt.plot(results['a1_m'] * 180 / np.pi, results['Tw1_m'])
plt.xlabel('Angle a1 (degrees)')
plt.ylabel('Driving Torque Tw1 (Nm)')
plt.title('Driving Torque vs Input Angle (General Solution)')
plt.grid(True)

# Link 2 CoM Position, Velocity, Acceleration (X-component) vs Time
plt.figure()
plt.plot(t_m, results['r21_m'], label='r21 (Pos)')
plt.plot(t_m, results['dr21_m'], label='dr21 (Vel)')
plt.plot(t_m, results['ddr21_m'], label='ddr21 (Acc)')
plt.xlabel('Time (s)')
plt.ylabel('Link 2 CoM X-Kinematics')
plt.title('Link 2 CoM X-Component Kinematics (General Solution)')
plt.legend()
plt.grid(True)

# Link 2 CoM Trajectory
plt.figure()
plt.plot(results['r21_m'], results['r22_m'])
plt.xlabel('r21 (m)')
plt.ylabel('r22 (m)')
plt.title('Link 2 CoM Trajectory (General Solution)')
plt.axis('equal')
plt.grid(True)

# Input Kinematics vs Time (Subplots)
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