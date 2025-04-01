% Parameters
l1 = 0.3; % Length of the first link
l2 = 0.65; % Length of the second link
b = 0.4; % Length of the line
H = 0.6; % Vertical distance from root

% Time parameters
T = 2; % Total time for the movement
dt = 0.01; % Time step
t = 0:dt:T; % Time vector

% Cycloid parameters
r = b / (2 * pi); % Radius to match the total distance
a = 1; % Ratio for normal cycloid

% Cycloid velocity profile
v_traj = r * (1 - cos(2 * pi * t / T));

% Position trajectory using cycloid velocity profile
x_traj = x_start + cumtrapz(t, v_traj);
y_traj = -H * ones(size(x_traj));

% Preallocate arrays for joint angles
theta1 = zeros(size(t));
theta2 = zeros(size(t));
phi1 = zeros(size(t));
phi2 = zeros(size(t));

% Calculate inverse kinematics
for i = 1:length(t)
    % Position of the end-effector
    x = x_traj(i);
    y = y_traj(i);
    
    % Inverse kinematics for branch A1-B1-E
    D = (x^2 + y^2 - l1^2 - l2^2) / (2 * l1 * l2);
    theta2(i) = atan2(sqrt(1 - D^2), D);
    theta1(i) = atan2(y, x) - atan2(l2 * sin(theta2(i)), l1 + l2 * cos(theta2(i)));
    
    % Inverse kinematics for branch A2-B2-E (mirrored)
    phi2(i) = -theta2(i);
    phi1(i) = atan2(y, x) + atan2(l2 * sin(phi2(i)), l1 + l2 * cos(phi2(i)));
end

% Calculate speed and acceleration
vx = diff(x_traj) / dt;
vy = diff(y_traj) / dt;
ax = diff(vx) / dt;
ay = diff(vy) / dt;

% Plotting
figure;
subplot(4, 1, 1);
plot(t, x_traj, 'r', t, y_traj, 'b');
title('End-Effector Trajectory');
xlabel('Time (s)');
ylabel('Position (m)');
legend('x', 'y');

subplot(4, 1, 2);
plot(t(1:end-1), vx, 'r', t(1:end-1), vy, 'b');
title('End-Effector Speed');
xlabel('Time (s)');
ylabel('Speed (m/s)');
legend('v_x', 'v_y');

subplot(4, 1, 3);
plot(t(1:end-2), ax, 'r', t(1:end-2), ay, 'b');
title('End-Effector Acceleration');
xlabel('Time (s)');
ylabel('Acceleration (m/s^2)');
legend('a_x', 'a_y');

subplot(4, 1, 4);
plot(t, theta1, 'r', t, theta2, 'b');
title('Joint Angles for Branch A1-B1-E');
xlabel('Time (s)');
ylabel('Angle (rad)');
legend('\theta_1', '\theta_2');

figure;
plot(t, phi1, 'r', t, phi2, 'b');
title('Joint Angles for Branch A2-B2-E');
xlabel('Time (s)');
ylabel('Angle (rad)');
legend('\phi_1', '\phi_2');