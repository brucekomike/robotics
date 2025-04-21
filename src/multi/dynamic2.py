import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# --------------------------------------------------------------------------
# 1. Symbolic Setup (Mirroring MATLAB - primarily for reference/verification)
#    NOTE: We will use the simplified analytical solutions provided later
#          in the MATLAB script for the actual calculations, as solving
#          the full system symbolically is complex and the script itself
#          switches to analytical forms.
# --------------------------------------------------------------------------
print("Setting up symbolic variables (for reference)...")
a1, r11, r12, da1, dr11, dr12, dda1, ddr11, ddr12 = sp.symbols('a1 r11 r12 da1 dr11 dr12 dda1 ddr11 ddr12')
a2, r21, r22, da2, dr21, dr22, dda2, ddr21, ddr22 = sp.symbols('a2 r21 r22 da2 dr21 dr22 dda2 ddr21 ddr22')
a3, r31, r32, da3, dr31, dr32, dda3, ddr31, ddr32 = sp.symbols('a3 r31 r32 da3 dr31 dr32 dda3 ddr31 ddr32')
L1, L2, L3, L4, t = sp.symbols('L1 L2 L3 L4 t')

# Assign link lengths (as constants in the MATLAB script)
L1_val, L2_val, L3_val, L4_val = 1, 1, 1, 1

# Define rotation matrices and vectors symbolically (example)
# A1 = sp.Matrix([[sp.cos(a1), -sp.sin(a1)], [sp.sin(a1), sp.cos(a1)]])
# R1 = sp.Matrix([r11, r12])
# ... and so on for constraint equations fun1 to fun4 ...
# func = sp.Matrix([...])
#
# Symbolic solve (can be very slow or fail for complex systems)
# resu = sp.solve(func.subs({L1: L1_val, ...}), [r11, r12, ...])
# ...
# Symbolic differentiation for velocity/acceleration
# q_i = a1
# q_d = sp.Matrix([r11, r12, r21, r22, a2, r31, r32, a3])
# Cq_d = func.jacobian(q_d)
# Cq_i = func.diff(q_i)
# dq_d = -Cq_d.inv() * Cq_i * da1
# ... etc ...
print("Symbolic setup reference complete.")
print("-" * 20)

# --------------------------------------------------------------------------
# 2. Use Provided Analytical Solutions (from MATLAB comments/assignments)
#    Define functions for these solutions using numpy for numerical evaluation.
# --------------------------------------------------------------------------
print("Defining numerical functions based on provided analytical solutions...")

def calculate_kinematics(a1_num, da1_num, dda1_num, L1=L1_val, L2=L2_val, L3=L3_val, L4=L4_val):
    """
    Calculates position, velocity, and acceleration based on the input angle a1
    and its derivatives, using the simplified analytical solutions provided
    in the MATLAB script (where L1=L2=L3=L4=1, a2=0, a3=a1).

    Args:
        a1_num (float or np.ndarray): Input angle a1 in radians.
        da1_num (float or np.ndarray): Input angular velocity da1 in rad/s.
        dda1_num (float or np.ndarray): Input angular acceleration dda1 in rad/s^2.
        L1, L2, L3, L4 (float): Link lengths (default to 1).

    Returns:
        dict: A dictionary containing all kinematic variables numerically.
              Keys match the variable names (e.g., 'r11', 'dr11', 'ddr11', etc.).
    """
    # Ensure using numpy functions for potential array inputs
    cos_a1 = np.cos(a1_num)
    sin_a1 = np.sin(a1_num)
    # tan_a1_half = np.tan(a1_num / 2) # Avoid tan(pi/2) issues if possible

    # --- Position ---
    # Based on the simplified forms in MATLAB:
    # r11 = cos(a1)/2;
    # r12 = sin(a1)/2;
    # r21 = -(tan(a1/2)^2 - 3)/(2*(tan(a1/2)^2 + 1)); -> simplifies to cos(a1)/2 + 1/2
    # r22 = (2*tan(a1/2))/(tan(a1/2)^2 + 1); -> simplifies to sin(a1)
    # a2 = 0;
    # r31 = (tan(a1/2)^2 + 3)/(2*tan(a1/2)^2 + 2); -> simplifies to 1/2 + cos(a1)/2
    # r32 = tan(a1/2)/(tan(a1/2)^2 + 1); -> simplifies to sin(a1)/2
    # a3 = a1;

    # Recalculating simplified forms using cos/sin directly (more stable)
    # From tan(x/2) = sin(x) / (1 + cos(x))
    # tan(a1/2)^2 = sin(a1)^2 / (1 + cos(a1))^2 = (1-cos(a1)^2) / (1 + cos(a1))^2 = (1-cos(a1))/(1+cos(a1))
    # tan(a1/2)^2 + 1 = (1-cos(a1))/(1+cos(a1)) + 1 = (1-cos(a1) + 1+cos(a1))/(1+cos(a1)) = 2/(1+cos(a1))
    # -(tan^2-3)/(2*(tan^2+1)) = -((1-c)/(1+c) - 3) / (2 * 2/(1+c))
    #                        = -( (1-c - 3(1+c)) / (1+c) ) * ( (1+c) / 4 )
    #                        = -( 1 - c - 3 - 3c ) / 4 = -( -2 - 4c ) / 4 = (2 + 4c) / 4 = 1/2 + c
    # (tan^2+3)/(2*(tan^2+2)) ??? Typo in matlab? Should be 2*(tan^2+1)? Assuming yes.
    # (tan^2+3)/(2*(tan^2+1)) = ( (1-c)/(1+c) + 3 ) / ( 4 / (1+c) )
    #                        = ( (1-c + 3(1+c)) / (1+c) ) * ( (1+c) / 4 )
    #                        = ( 1 - c + 3 + 3c ) / 4 = ( 4 + 2c ) / 4 = 1 + c/2 -> Check r31 again
    # Let's use the structure definitions:
    # r1a = R1 + A1*[-L1/2, 0]' = [0,0]' => R1 = -A1*[-L1/2, 0]' = A1*[L1/2, 0]'
    # R1 = [cos(a1), -sin(a1); sin(a1), cos(a1)] * [L1/2, 0]' = [L1/2*cos(a1), L1/2*sin(a1)]'
    r11 = (L1 / 2.0) * cos_a1
    r12 = (L1 / 2.0) * sin_a1
    # r3d = R3 + A3*[-L3/2, 0]' = [L4, 0]'
    # With a3=a1, L3=1: R3 = [L4, 0]' - A1*[-L3/2, 0]' = [L4, 0]' + A1*[L3/2, 0]'
    # R3 = [L4, 0]' + [L3/2*cos(a1), L3/2*sin(a1)]' = [L4 + L3/2*cos(a1), L3/2*sin(a1)]'
    r31 = L4 + (L3 / 2.0) * cos_a1
    r32 = (L3 / 2.0) * sin_a1
    a3 = a1_num
    # r2b = R2 + A2*[-L2/2, 0]' = r1b = R1 + A1*[L1/2, 0]'
    # With a2=0, A2=I: R2 = R1 + A1*[L1/2, 0]' - I*[-L2/2, 0]'
    # R2 = [L1/2*c, L1/2*s]' + [L1/2*c, L1/2*s]' - [-L2/2, 0]' = [L1*c + L2/2, L1*s]'
    r21 = L1 * cos_a1 + L2 / 2.0
    r22 = L1 * sin_a1
    a2 = 0.0 # Constant based on the simplification

    # --- Velocity ---
    # Based on the simplified forms in MATLAB:
    # dr11 = -(da1*sin(a1))/2;
    # dr12 = (da1*cos(a1))/2;
    # dr21 = -(da1*(...))/(2*sin(a2 - a3)); -> Simplified with a2=0, a3=a1
    # dr22 = -(da1*(...))/(2*sin(a2 - a3)); -> Simplified with a2=0, a3=a1
    # da2 = 0;
    # dr31 = (da1*sin(a3)*sin(a1 - a2))/(2*sin(a2 - a3)); -> Simplified
    # dr32 = -(da1*cos(a3)*sin(a1 - a2))/(2*sin(a2 - a3)); -> Simplified
    # da3 = da1;

    # Calculate simplified velocities directly by differentiating positions
    dr11 = -(L1 / 2.0) * sin_a1 * da1_num
    dr12 = (L1 / 2.0) * cos_a1 * da1_num
    dr21 = -L1 * sin_a1 * da1_num
    dr22 = L1 * cos_a1 * da1_num
    da2 = 0.0
    dr31 = -(L3 / 2.0) * sin_a1 * da1_num # Since a3=a1
    dr32 = (L3 / 2.0) * cos_a1 * da1_num  # Since a3=a1
    da3 = da1_num

    # --- Acceleration ---
    # Based on the simplified forms in MATLAB:
    # ddr11 = - (dda1*sin(a1))/2 - (da1^2*cos(a1))/2;
    # ddr12 = (dda1*cos(a1))/2 - (da1^2*sin(a1))/2;
    # ddr21 = (...)/(2*sin(a2 - a3)); -> Simplified
    # ddr22 = -(...)/(2*sin(a2 - a3)); -> Simplified
    # dda2 = 0;
    # ddr31 = (...)/(2*sin(a2 - a3)); -> Simplified
    # ddr32 = -(...)/(2*sin(a2 - a3)); -> Simplified
    # dda3 = dda1;

    # Calculate simplified accelerations directly by differentiating velocities
    da1_sq = da1_num**2
    ddr11 = -(L1 / 2.0) * (dda1_num * sin_a1 + da1_sq * cos_a1)
    ddr12 = (L1 / 2.0) * (dda1_num * cos_a1 - da1_sq * sin_a1)
    ddr21 = -L1 * (dda1_num * sin_a1 + da1_sq * cos_a1)
    ddr22 = L1 * (dda1_num * cos_a1 - da1_sq * sin_a1)
    dda2 = 0.0
    ddr31 = -(L3 / 2.0) * (dda1_num * sin_a1 + da1_sq * cos_a1) # Since a3=a1, da3=da1, dda3=dda1
    ddr32 = (L3 / 2.0) * (dda1_num * cos_a1 - da1_sq * sin_a1) # Since a3=a1, da3=da1, dda3=dda1
    dda3 = dda1_num

    return {
        'r11': r11, 'r12': r12, 'a1': a1_num,
        'r21': r21, 'r22': r22, 'a2': a2,
        'r31': r31, 'r32': r32, 'a3': a3,
        'dr11': dr11, 'dr12': dr12, 'da1': da1_num,
        'dr21': dr21, 'dr22': dr22, 'da2': da2,
        'dr31': dr31, 'dr32': dr32, 'da3': da3,
        'ddr11': ddr11, 'ddr12': ddr12, 'dda1': dda1_num,
        'ddr21': ddr21, 'ddr22': ddr22, 'dda2': dda2,
        'ddr31': ddr31, 'ddr32': ddr32, 'dda3': dda3,
    }
# here using a erion with corrected type handling
def calculate_kinematics(a1_num, da1_num, dda1_num, L1=L1_val, L2=L2_val, L3=L3_val, L4=L4_val):
    """
    Calculates position, velocity, and acceleration based on the input angle a1
    and its derivatives, using the simplified analytical solutions provided
    in the MATLAB script (where L1=L2=L3=L4=1, a2=0, a3=a1).

    Args:
        a1_num (float or np.ndarray): Input angle a1 in radians.
        da1_num (float or np.ndarray): Input angular velocity da1 in rad/s.
        dda1_num (float or np.ndarray): Input angular acceleration dda1 in rad/s^2.
        L1, L2, L3, L4 (float): Link lengths (default to 1).

    Returns:
        dict: A dictionary containing all kinematic variables numerically.
              Keys match the variable names (e.g., 'r11', 'dr11', 'ddr11', etc.).
    """
    # Ensure using numpy functions for potential array inputs
    cos_a1 = np.cos(a1_num)
    sin_a1 = np.sin(a1_num)

    # --- Position ---
    r11 = (L1 / 2.0) * cos_a1
    r12 = (L1 / 2.0) * sin_a1
    r31 = L4 + (L3 / 2.0) * cos_a1
    r32 = (L3 / 2.0) * sin_a1
    a3 = a1_num # Array if a1_num is array
    r21 = L1 * cos_a1 + L2 / 2.0
    r22 = L1 * sin_a1
    # a2 = 0.0 # OLD - Scalar
    a2 = np.zeros_like(a1_num) # NEW - Array

    # --- Velocity ---
    dr11 = -(L1 / 2.0) * sin_a1 * da1_num
    dr12 = (L1 / 2.0) * cos_a1 * da1_num
    dr21 = -L1 * sin_a1 * da1_num
    dr22 = L1 * cos_a1 * da1_num
    # da2 = 0.0 # OLD - Scalar
    da2 = np.zeros_like(a1_num) # NEW - Array
    dr31 = -(L3 / 2.0) * sin_a1 * da1_num # Since a3=a1
    dr32 = (L3 / 2.0) * cos_a1 * da1_num  # Since a3=a1
    da3 = da1_num # Array if da1_num is array

    # --- Acceleration ---
    da1_sq = da1_num**2
    ddr11 = -(L1 / 2.0) * (dda1_num * sin_a1 + da1_sq * cos_a1)
    ddr12 = (L1 / 2.0) * (dda1_num * cos_a1 - da1_sq * sin_a1)
    ddr21 = -L1 * (dda1_num * sin_a1 + da1_sq * cos_a1)
    ddr22 = L1 * (dda1_num * cos_a1 - da1_sq * sin_a1)
    # dda2 = 0.0 # OLD - Scalar
    dda2 = np.zeros_like(a1_num) # NEW - Array
    ddr31 = -(L3 / 2.0) * (dda1_num * sin_a1 + da1_sq * cos_a1) # Since a3=a1, da3=da1, dda3=dda1
    ddr32 = (L3 / 2.0) * (dda1_num * cos_a1 - da1_sq * sin_a1) # Since a3=a1, da3=da1, dda3=dda1
    dda3 = dda1_num # Array if dda1_num is array

    return {
        'r11': r11, 'r12': r12, 'a1': a1_num,
        'r21': r21, 'r22': r22, 'a2': a2, # Now returns array
        'r31': r31, 'r32': r32, 'a3': a3,
        'dr11': dr11, 'dr12': dr12, 'da1': da1_num,
        'dr21': dr21, 'dr22': dr22, 'da2': da2, # Now returns array
        'dr31': dr31, 'dr32': dr32, 'da3': da3,
        'ddr11': ddr11, 'ddr12': ddr12, 'dda1': dda1_num,
        'ddr21': ddr21, 'ddr22': ddr22, 'dda2': dda2, # Now returns array
        'ddr31': ddr31, 'ddr32': ddr32, 'dda3': dda3,
    }
print("Numerical functions defined.")
print("-" * 20)

# --------------------------------------------------------------------------
# 3. Numerical Simulation Parameters
# --------------------------------------------------------------------------
print("Setting up simulation parameters...")
t_start = 0.0
t_end = 10.0
dt = 0.01
t_m = np.arange(t_start, t_end + dt, dt)
nn = len(t_m)

w = 10.0  # Angular frequency for input motion
g = 10.0  # Gravitational acceleration
m1, m2, m3 = 1.0, 1.0, 1.0 # Masses

print(f"Time steps: {nn}")
print("-" * 20)

# --------------------------------------------------------------------------
# 4. Define Input Motion
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
# 5. Run Numerical Simulation Loop (Pre-calculation)
# --------------------------------------------------------------------------
print("Running numerical simulation to pre-calculate states...")

# Calculate input motion values
a1_vals = input_a1(t_m, w)
da1_vals = input_da1(t_m, w)
dda1_vals = input_dda1(t_m, w)

# Initialize storage arrays (using NaN for potential debugging)
kinematics_results = {}
var_names = [
    'r11', 'r12', 'a1', 'r21', 'r22', 'a2', 'r31', 'r32', 'a3',
    'dr11', 'dr12', 'da1', 'dr21', 'dr22', 'da2', 'dr31', 'dr32', 'da3',
    'ddr11', 'ddr12', 'dda1', 'ddr21', 'ddr22', 'dda2', 'ddr31', 'ddr32', 'dda3'
]
for name in var_names:
    kinematics_results[name + '_m'] = np.full(nn, np.nan)

# Calculate all kinematics using the vectorized function
kinematics_all_steps = calculate_kinematics(a1_vals, da1_vals, dda1_vals, L1_val, L2_val, L3_val, L4_val)

# Store results
for name in var_names:
     kinematics_results[name + '_m'] = kinematics_all_steps[name]

# --- Calculate Driving Torque Tw1 (Inverse Dynamics) ---
# Translating the complex Tw1_m formula from MATLAB
# CAUTION: This formula seems specific to the simplified a2=0, a3=a1 case
#          and assumes L1=L2=L3=L4=1. Double-check its derivation if needed.
# Using the second formula provided in the MATLAB script
print("Calculating driving torque...")
Tw1_m = np.full(nn, np.nan)

# Use calculated values directly
a1_n = kinematics_results['a1_m']
a2_n = kinematics_results['a2_m'] # Should be 0
a3_n = kinematics_results['a3_m'] # Should be a1
da1_n = kinematics_results['da1_m']
da2_n = kinematics_results['da2_m'] # Should be 0
da3_n = kinematics_results['da3_m'] # Should be da1
dda1_n = kinematics_results['dda1_m']
dda2_n = kinematics_results['dda2_m'] # Should be 0
dda3_n = kinematics_results['dda3_m'] # Should be dda1

cos = np.cos
sin = np.sin

# Simplified calculation since a2=0, a3=a1, da2=0, dda2=0, da3=da1, dda3=dda1
# Denominator: 12*(cos(2*a2 - 2*a3) - 1) = 12*(cos(-2*a1) - 1) = 12*(cos(2*a1) - 1)
# Avoid division by zero when cos(2*a1)=1, i.e., 2*a1 = 2*k*pi => a1 = k*pi
denominator = 12.0 * (cos(2 * a1_n) - 1.0)
# Add a small epsilon to prevent division by zero, or handle singularity if needed
denominator[np.abs(denominator) < 1e-9] = 1e-9

term1 = 3*g*m2*cos(a1_n) # cos(a1 - 2*a2) = cos(a1)
term2 = -10*dda1_n*m2
term3 = -3*dda1_n*m3
term4 = -5*dda1_n*m1
term5 = -3*g*m2*cos(a1_n - 2*a3_n) # cos(a1 - 2*a1) = cos(-a1) = cos(a1)
term6 = 3*g*m3*cos(a1_n) # cos(a1 - 2*a2) = cos(a1)
term7 = -3*g*m3*cos(a1_n - 2*a3_n) # cos(a1 - 2*a1) = cos(-a1) = cos(a1)
term8 = dda1_n*m1*cos(2*a1_n) # cos(2*a1 - 2*a2) = cos(2*a1)
term9 = 6*dda1_n*m2*cos(2*a1_n) # cos(2*a1 - 2*a2) = cos(2*a1)
term10 = 4*dda1_n*m1*cos(-2*a1_n) # cos(2*a2 - 2*a3) = cos(-2*a1) = cos(2*a1)
term11 = -2*dda1_n*m2*cos(2*a1_n - 2*a3_n) # cos(2*a1 - 2*a1) = cos(0) = 1
term12 = 3*dda1_n*m3*cos(2*a1_n) # cos(2*a1 - 2*a2) = cos(2*a1)
term13 = 6*dda1_n*m2*cos(-2*a1_n) # cos(2*a2 - 2*a3) = cos(-2*a1) = cos(2*a1)
term14 = 0 # -2*da2^2*... = 0
term15 = 0 # -10*da2^2*... = 0
term16 = 0 # -6*da2^2*... = 0
term17 = da3_n**2*m1*sin(a1_n - a3_n) # sin(a1-a1) = sin(0) = 0
term18 = 2*da3_n**2*m2*sin(a1_n - a3_n) # sin(0) = 0
term19 = 3*da3_n**2*m3*sin(a1_n - a3_n) # sin(0) = 0
term20 = 3*g*m1*cos(a1_n - 2*a2_n + 2*a3_n) # cos(a1 + 2*a1) = cos(3*a1)
term21 = 3*g*m1*cos(a1_n + 2*a2_n - 2*a3_n) # cos(a1 - 2*a1) = cos(-a1) = cos(a1)
term22 = 6*g*m2*cos(a1_n - 2*a2_n + 2*a3_n) # cos(3*a1)
term23 = 3*g*m2*cos(a1_n + 2*a2_n - 2*a3_n) # cos(a1)
term24 = 3*g*m3*cos(a1_n - 2*a2_n + 2*a3_n) # cos(3*a1)
term25 = -da1_n**2*m1*sin(2*a1_n) # sin(2*a1 - 2*a2) = sin(2*a1)
term26 = -6*da1_n**2*m2*sin(2*a1_n) # sin(2*a1 - 2*a2) = sin(2*a1)
term27 = 2*da1_n**2*m2*sin(2*a1_n - 2*a3_n) # sin(2*a1 - 2*a1) = sin(0) = 0
term28 = -3*da1_n**2*m3*sin(2*a1_n) # sin(2*a1 - 2*a2) = sin(2*a1)
term29 = 0 # 2*da2^2*... = 0
term30 = da3_n**2*m1*sin(a1_n - 2*a2_n + a3_n) # sin(a1 + a1) = sin(2*a1)
term31 = 6*da3_n**2*m2*sin(a1_n - 2*a2_n + a3_n) # sin(2*a1)
term32 = 3*da3_n**2*m3*sin(a1_n - 2*a2_n + a3_n) # sin(2*a1)
term33 = -6*g*m1*cos(a1_n)
term34 = -9*g*m2*cos(a1_n)
term35 = -3*g*m3*cos(a1_n)

# Summing terms that don't cancel due to a2=0, a3=a1 etc.
numerator = (term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8 + term9 +
             term10 + term11 + term12 + term13 + # terms 14-19 are 0
             term20 + term21 + term22 + term23 + term24 +
             term25 + term26 + # term 27 is 0
             term28 + # term 29 is 0
             term30 + term31 + term32 + term33 + term34 + term35)

Tw1_m = numerator / denominator

print("Simulation complete.")
print("-" * 20)

# --------------------------------------------------------------------------
# 6. Animation Setup (Matplotlib)
# --------------------------------------------------------------------------
print("Setting up animation...")

fig, ax = plt.subplots()
ax.set_aspect('equal')
ax.set_xlim(-1.5 * L1_val, 1.5 * (L4_val + L3_val))
ax.set_ylim(-1.5 * L1_val, 1.5 * L1_val)
ax.set_xlabel("X position")
ax.set_ylabel("Y position")
ax.set_title("Four-Bar Linkage Animation (Numerical Solution)")
ax.grid(True)

# Define fixed points
P0 = np.array([0, 0])
P3 = np.array([L4_val, 0])

# Initialize lines for the links
line1, = ax.plot([], [], 'o-', lw=2, color='blue', label='Link 1 (L1)') # P0 to P1
line2, = ax.plot([], [], 'o-', lw=2, color='red', label='Link 2 (L2)')  # P1 to P2
line3, = ax.plot([], [], 'o-', lw=2, color='green', label='Link 3 (L3)') # P2 to P3
line4, = ax.plot([P0[0], P3[0]], [P0[1], P3[1]], 's-', lw=1, color='black', label='Link 4 (Fixed L4)') # Base P0 to P3

# Add a time text
time_template = 'time = %.2fs'
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)

ax.legend(loc='upper right')

# Function to calculate joint positions P1 and P2
def get_joint_positions(k):
    """Calculates joint positions P1 and P2 for frame k."""
    # Get CoM and angles for this frame
    r11_k = kinematics_results['r11_m'][k]
    r12_k = kinematics_results['r12_m'][k]
    a1_k = kinematics_results['a1_m'][k]
    r21_k = kinematics_results['r21_m'][k]
    r22_k = kinematics_results['r22_m'][k]
    a2_k = kinematics_results['a2_m'][k] # This is 0

    # Calculate rotation matrices numerically
    A1_k = np.array([[np.cos(a1_k), -np.sin(a1_k)],
                     [np.sin(a1_k), np.cos(a1_k)]])
    # A2_k = np.array([[np.cos(a2_k), -np.sin(a2_k)],
    #                  [np.sin(a2_k), np.cos(a2_k)]]) # This is Identity

    # Define local vectors to joints relative to CoM
    u_1b_local = np.array([0.5 * L1_val, 0]) # Joint P1 relative to R1
    u_2c_local = np.array([0.5 * L2_val, 0]) # Joint P2 relative to R2

    # Calculate global joint positions
    R1_k = np.array([r11_k, r12_k])
    R2_k = np.array([r21_k, r22_k])

    P1_k = R1_k + A1_k @ u_1b_local
    # P2_k = R2_k + A2_k @ u_2c_local # Since A2=I
    P2_k = R2_k + u_2c_local

    return P1_k, P2_k

# Initialization function for the animation
def init():
    line1.set_data([], [])
    line2.set_data([], [])
    line3.set_data([], [])
    time_text.set_text('')
    return line1, line2, line3, time_text

# Update function for the animation
def update(frame):
    P1_k, P2_k = get_joint_positions(frame)

    line1.set_data([P0[0], P1_k[0]], [P0[1], P1_k[1]])
    line2.set_data([P1_k[0], P2_k[0]], [P1_k[1], P2_k[1]])
    line3.set_data([P2_k[0], P3[0]], [P2_k[1], P3[1]])
    time_text.set_text(time_template % t_m[frame])
    return line1, line2, line3, time_text

# Create animation
# Interval is dt in milliseconds
ani = animation.FuncAnimation(fig, update, frames=nn,
                              init_func=init, blit=True, interval=dt*1000)

print("Animation setup complete. Starting plots and animation...")
print("-" * 20)

# --------------------------------------------------------------------------
# 7. Display Static Plots (from MATLAB script)
# --------------------------------------------------------------------------

# Torque vs Time
plt.figure()
plt.plot(t_m, Tw1_m)
plt.xlabel('Time (s)')
plt.ylabel('Driving Torque Tw1 (Nm)')
plt.title('Driving Torque vs Time')
plt.grid(True)

# Angle vs Time
plt.figure()
plt.plot(t_m, kinematics_results['a1_m'] * 180 / np.pi)
plt.xlabel('Time (s)')
plt.ylabel('Angle a1 (degrees)')
plt.title('Input Angle vs Time')
plt.grid(True)

# Torque vs Angle
plt.figure()
plt.plot(kinematics_results['a1_m'] * 180 / np.pi, Tw1_m)
plt.xlabel('Angle a1 (degrees)')
plt.ylabel('Driving Torque Tw1 (Nm)')
plt.title('Driving Torque vs Input Angle')
plt.grid(True)

# Link 2 CoM Position, Velocity, Acceleration (X-component) vs Time
plt.figure()
plt.plot(t_m, kinematics_results['r21_m'], label='r21 (Position)')
plt.plot(t_m, kinematics_results['dr21_m'], label='dr21 (Velocity)')
plt.plot(t_m, kinematics_results['ddr21_m'], label='ddr21 (Acceleration)')
plt.xlabel('Time (s)')
plt.ylabel('Link 2 CoM X-Kinematics')
plt.title('Link 2 Center of Mass X-Component Kinematics')
plt.legend()
plt.grid(True)

# Link 2 CoM Position (X) vs Time
plt.figure()
plt.plot(t_m, kinematics_results['r21_m'])
plt.xlabel('Time (s)')
plt.ylabel('r21 (m)')
plt.title('Link 2 CoM X-Position vs Time')
plt.grid(True)

# Link 2 CoM Position (Y) vs Time
plt.figure()
plt.plot(t_m, kinematics_results['r22_m'])
plt.xlabel('Time (s)')
plt.ylabel('r22 (m)')
plt.title('Link 2 CoM Y-Position vs Time')
plt.grid(True)

# Link 2 CoM Trajectory
plt.figure()
plt.plot(kinematics_results['r21_m'], kinematics_results['r22_m'])
plt.xlabel('r21 (m)')
plt.ylabel('r22 (m)')
plt.title('Link 2 CoM Trajectory')
plt.axis('equal')
plt.grid(True)

# Input Angle, Velocity, Acceleration vs Time (Subplots)
fig_kin, axs_kin = plt.subplots(3, 1, sharex=True)
fig_kin.suptitle('Input Kinematics vs Time')
axs_kin[0].plot(t_m, kinematics_results['dda1_m'] * 180 / np.pi)
axs_kin[0].set_ylabel('dda1 (deg/s^2)')
axs_kin[0].grid(True)
axs_kin[1].plot(t_m, kinematics_results['da1_m'] * 180 / np.pi)
axs_kin[1].set_ylabel('da1 (deg/s)')
axs_kin[1].grid(True)
axs_kin[2].plot(t_m, kinematics_results['a1_m'] * 180 / np.pi)
axs_kin[2].set_ylabel('a1 (deg)')
axs_kin[2].set_xlabel('Time (s)')
axs_kin[2].grid(True)

# Phase Plot (Input Velocity vs Angle)
plt.figure()
plt.plot(kinematics_results['a1_m'] * 180 / np.pi, kinematics_results['da1_m'] * 180 / np.pi)
plt.xlabel('Angle a1 (degrees)')
plt.ylabel('Angular Velocity da1 (deg/s)')
plt.title('Input Phase Plot (Velocity vs Angle)')
plt.grid(True)

# Show the animation window and all plot windows
plt.show()

print("Script finished.")