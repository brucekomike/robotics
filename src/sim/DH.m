clear
clc
pkg load symbolic
% Define symbolic variables

syms theta1 theta2 theta3 l1 l2 l3 l4


#theta1=35*pi/180;theta2=40*pi/180;theta3=-25*pi/180;l1=240;l2=60;l3=150;l4=70

% Rotation matrix around Z-axis
Rz = @(theta) [[cos(theta), -sin(theta), 0, 0];
               [sin(theta), cos(theta), 0, 0];
               [0, 0, 1, 0];
               [0, 0, 0, 1]];

% Rotation matrix around Y-axis
Ry = @(theta) [[cos(theta), 0, sin(theta), 0];
               [0, 1, 0, 0];
               [-sin(theta), 0, cos(theta), 0];
               [0, 0, 0, 1]];

% Rotation matrix around X-axis
Rx = @(theta) [[1, 0, 0, 0];
               [0, cos(theta), -sin(theta), 0];
               [0, sin(theta), cos(theta), 0];
               [0, 0, 0, 1]];

% Translation matrix
T = @(x, y, z) [[1, 0, 0, x];
                [0, 1, 0, y];
                [0, 0, 1, z];
                [0, 0, 0, 1]];

% Calculate the total transformation matrix
T_total = Rz(theta1) * T(-l2, 0, l1) * Rx(theta2) * T(0, l3, 0) * Rx(theta3) * T(0, l4, 0)

% Extract the position vector
position = T_total(1:3, 4);

% Display the position
disp('position');
disp(position);