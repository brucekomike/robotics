import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import argparse
from matplotlib import gridspec

# Constants
l1 = 29  # Length of first arm segment
l2 = 26  # Length of second arm segment
max_reach = l1 + l2

def inverse_kinematics(x, y, elbow_up=True):
    cos_theta2 = (x**2 + y**2 - l1**2 - l2**2) / (2 * l1 * l2)
    if cos_theta2 < -1 or cos_theta2 > 1:
        return None, None
    
    sin_theta2 = np.sqrt(1 - cos_theta2**2) if elbow_up else -np.sqrt(1 - cos_theta2**2)
    theta2 = np.arctan2(sin_theta2, cos_theta2)

    k1 = l1 + l2 * cos_theta2
    k2 = l2 * sin_theta2
    theta1 = np.arctan2(y, x) - np.arctan2(k2, k1)
    return theta1, theta2

def adjust_target(init_pos, theta):
    init_x, init_y = init_pos
    a = 1.0
    b = 2 * (init_x * np.cos(theta) + init_y * np.sin(theta))
    c = init_x**2 + init_y**2 - max_reach**2

    discriminant = b**2 - 4 * a * c
    if discriminant < 0:
        return init_x, init_y

    sqrt_discriminant = np.sqrt(discriminant)
    t_max = max((-b + sqrt_discriminant)/(2*a), (-b - sqrt_discriminant)/(2*a)) - 0.01
    return (init_x + t_max * np.cos(theta), 
            init_y + t_max * np.sin(theta))

def simulate_motion(init_pos, target_pos, step_ms, total_time, elbow_up=True, profile='linear'):
    if total_time <= 0 or step_ms <= 0:
        raise ValueError("Time parameters must be positive")

    num_steps = max(int(np.ceil(total_time * 1000 / step_ms)), 2)
    time = np.linspace(0, total_time, num_steps)
    dt_value = time[1] - time[0] if len(time) > 1 else total_time

    init_x, init_y = init_pos
    target_x, target_y = target_pos
    dx = target_x - init_x
    dy = target_y - init_y
    
    x_positions = []
    y_positions = []
    for t in time:
        if profile == 'linear':
            s = t / total_time
        elif profile == 'trapezoidal':
            t_ramp = total_time / 4
            if t < t_ramp:
                s = 0.5 * (16/(3*total_time**2)) * t**2
            elif t < 3*t_ramp:
                s = (4/(3*total_time)) * t - 1/6
            else:
                s = 1 - 0.5 * (16/(3*total_time**2)) * (total_time - t)**2
        elif profile == 's_curve':
            u = t / total_time
            s = 10*u**3 - 15*u**4 + 6*u**5
        else:
            raise ValueError(f"Invalid profile: {profile}")
        
        x_positions.append(init_x + s * dx)
        y_positions.append(init_y + s * dy)

    theta1_values = []
    theta2_values = []
    for x, y in zip(x_positions, y_positions):
        theta1, theta2 = inverse_kinematics(x, y, elbow_up)
        if theta1 is None or theta2 is None:
            raise ValueError(f"Unreachable position at ({x:.2f}, {y:.2f})")
        theta1_values.append(theta1)
        theta2_values.append(theta2)

    theta1 = np.array(theta1_values)
    theta2 = np.array(theta2_values)
    
    with np.errstate(divide='ignore', invalid='ignore'):
        theta1_vel = np.gradient(theta1, dt_value)
        theta2_vel = np.gradient(theta2, dt_value)
        theta1_acc = np.gradient(theta1_vel, dt_value)
        theta2_acc = np.gradient(theta2_vel, dt_value)

    return (time, x_positions, y_positions, 
            theta1, theta2, theta1_vel, theta2_vel, 
            theta1_acc, theta2_acc)

def animate_motion(data, elbow_up, step_ms, profile):
    (time, x, y, theta1, theta2, 
     theta1_vel, theta2_vel, theta1_acc, theta2_acc) = data
    
    fig = plt.figure(figsize=(18, 12))
    gs = gridspec.GridSpec(4, 2, figure=fig, width_ratios=[1.5, 1], height_ratios=[1,1,1,1])
    
    # Initialize plot elements
    ax1 = fig.add_subplot(gs[0:2, 0])
    ax1.set(xlim=(-max_reach, max_reach), ylim=(-max_reach, max_reach))
    ax1.set_title(f'Arm Motion - {profile} Profile (Elbow {"Up" if elbow_up else "Down"})')
    arm_line, = ax1.plot([], [], 'o-', lw=2)
    path_line, = ax1.plot([], [], 'g--', alpha=0.5)
    
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.set_title('X Position vs Time')
    x_line, = ax2.plot([], [], 'r-')
    
    ax3 = fig.add_subplot(gs[1, 1])
    ax3.set_title('Y Position vs Time')
    y_line, = ax3.plot([], [], 'b-')
    
    ax4 = fig.add_subplot(gs[2, 0])
    ax4.set_title('Joint Angles')
    ax4.set_ylabel('Radians')
    theta1_line, = ax4.plot([], [], 'm-', label='θ₁')
    theta2_line, = ax4.plot([], [], 'c-', label='θ₂')
    ax4.legend()
    
    ax5 = fig.add_subplot(gs[2, 1])
    ax5.set_title('Joint Velocities')
    ax5.set_ylabel('rad/s')
    vel1_line, = ax5.plot([], [], 'm-', label='ω₁')
    vel2_line, = ax5.plot([], [], 'c-', label='ω₂')
    ax5.legend()
    
    ax6 = fig.add_subplot(gs[3, :])
    ax6.set_title('Joint Accelerations')
    ax6.set_ylabel('rad/s²')
    ax6.set_xlabel('Time (s)')
    acc1_line, = ax6.plot([], [], 'm-', label='α₁')
    acc2_line, = ax6.plot([], [], 'c-', label='α₂')
    ax6.legend()
    
    plt.tight_layout()

    def init():
        """Initialize animation elements to empty"""
        arm_line.set_data([], [])
        path_line.set_data([], [])
        x_line.set_data([], [])
        y_line.set_data([], [])
        theta1_line.set_data([], [])
        theta2_line.set_data([], [])
        vel1_line.set_data([], [])
        vel2_line.set_data([], [])
        acc1_line.set_data([], [])
        acc2_line.set_data([], [])
        return (arm_line, path_line, x_line, y_line, 
                theta1_line, theta2_line, vel1_line, vel2_line, 
                acc1_line, acc2_line)

    def update(frame):
        """Update animation frame"""
        x_arm = [0, l1 * np.cos(theta1[frame]), x[frame]]
        y_arm = [0, l1 * np.sin(theta1[frame]), y[frame]]
        arm_line.set_data(x_arm, y_arm)
        path_line.set_data(x[:frame+1], y[:frame+1])
        
        x_line.set_data(time[:frame+1], x[:frame+1])
        y_line.set_data(time[:frame+1], y[:frame+1])
        
        theta1_line.set_data(time[:frame+1], theta1[:frame+1])
        theta2_line.set_data(time[:frame+1], theta2[:frame+1])
        
        vel1_line.set_data(time[:frame+1], theta1_vel[:frame+1])
        vel2_line.set_data(time[:frame+1], theta2_vel[:frame+1])
        
        acc1_line.set_data(time[:frame+1], theta1_acc[:frame+1])
        acc2_line.set_data(time[:frame+1], theta2_acc[:frame+1])
        
        for ax in [ax2, ax3, ax4, ax5, ax6]:
            ax.relim()
            ax.autoscale_view()
        
        return (arm_line, path_line, x_line, y_line, 
                theta1_line, theta2_line, vel1_line, vel2_line, 
                acc1_line, acc2_line)

    ani = FuncAnimation(fig, update, frames=len(time), 
                       init_func=init, interval=step_ms, blit=True)
    plt.show()
    
    # Generate summary plots after animation window closes
    create_summary_plots(data, elbow_up, profile)
    create_joint_angle_plot(data, elbow_up, profile)

def create_summary_plots(data, elbow_up, profile):
    """Create comprehensive summary plots of all variables"""
    time, x, y, theta1, theta2, theta1_vel, theta2_vel, theta1_acc, theta2_acc = data
    
    fig = plt.figure(figsize=(12, 16))
    gs = gridspec.GridSpec(4, 1, figure=fig)
    
    # Position plots
    ax1 = fig.add_subplot(gs[0])
    ax1.plot(time, x, 'r-', label='X Position')
    ax1.plot(time, y, 'b-', label='Y Position')
    ax1.set_ylabel('Position (units)')
    ax1.legend()
    
    # Joint angles
    ax2 = fig.add_subplot(gs[1])
    ax2.plot(time, theta1, 'm-', label='θ₁')
    ax2.plot(time, theta2, 'c-', label='θ₂')
    ax2.set_ylabel('Joint Angles (rad)')
    ax2.legend()
    
    # Velocities
    ax3 = fig.add_subplot(gs[2])
    ax3.plot(time, theta1_vel, 'm--', label='ω₁')
    ax3.plot(time, theta2_vel, 'c--', label='ω₂')
    ax3.set_ylabel('Angular Velocity (rad/s)')
    ax3.legend()
    
    # Accelerations
    ax4 = fig.add_subplot(gs[3])
    ax4.plot(time, theta1_acc, 'm-.', label='α₁')
    ax4.plot(time, theta2_acc, 'c-.', label='α₂')
    ax4.set_xlabel('Time (s)')
    ax4.set_ylabel('Angular Acceleration (rad/s²)')
    ax4.legend()
    
    fig.suptitle(f'Motion Summary - {profile} Profile (Elbow {"Up" if elbow_up else "Down"})')
    plt.tight_layout()
    plt.show()

def create_joint_angle_plot(data, elbow_up, profile):
    """Create dedicated joint angle vs time plot"""
    time, _, _, theta1, theta2, *_ = data
    
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.plot(time, theta1, 'm-', label='θ₁')
    ax.plot(time, theta2, 'c-', label='θ₂')
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Joint Angle (rad)')
    ax.set_title(f'Joint Angles vs Time - {profile} Profile (Elbow {"Up" if elbow_up else "Down"})')
    ax.legend()
    plt.tight_layout()
    plt.show()

def main():
    parser = argparse.ArgumentParser(description="Robotic Arm Motion Analysis")
    parser.add_argument('--init_x', type=float, default=0.0, help='Initial X position')
    parser.add_argument('--init_y', type=float, default=49.0, help='Initial Y position')
    parser.add_argument('--theta', type=float, required=True, help='Movement direction in degrees')
    parser.add_argument('--step', type=int, default=50, help='Time step in milliseconds (minimum 10)')
    parser.add_argument('--time', type=float, default=5.0, help='Total motion time in seconds (minimum 0.1)')
    parser.add_argument('--profile', choices=['linear', 'trapezoidal', 's_curve'], 
                       default='linear', help='Acceleration profile type')
    args = parser.parse_args()

    if args.time < 0.1 or args.step < 10:
        raise ValueError("Invalid time parameters")

    init_pos = (args.init_x, args.init_y)
    theta = np.radians(args.theta)
    target_pos = adjust_target(init_pos, theta)

    if np.hypot(*init_pos) > max_reach:
        raise ValueError("Initial position is unreachable")

    config_data = []
    for elbow_up in [True, False]:
        try:
            motion_data = simulate_motion(
                init_pos, target_pos, args.step, args.time, 
                elbow_up, args.profile
            )
            config_data.append((motion_data[3], motion_data[4], 
                              f'Elbow {"Up" if elbow_up else "Down"}'))
            animate_motion(motion_data, elbow_up, args.step, args.profile)
        except ValueError as e:
            print(f"Skipping {'Elbow Up' if elbow_up else 'Elbow Down'}: {e}")

    # Configuration space plot
    plt.figure(figsize=(10, 8))
    for theta1, theta2, label in config_data:
        plt.plot(theta1, theta2, label=label)
    plt.title('Joint Angle Correlation')
    plt.xlabel('θ₁ (rad)')
    plt.ylabel('θ₂ (rad)')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
