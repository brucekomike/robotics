%% 位置正逆解模型 (Position Forward and Inverse Kinematics Model)
% 2017.11.14
clc            % 清除命令窗口 (Clear command window)
clear all      % 清除工作区所有变量 (Clear all variables in workspace)
close all      % 关闭所有图窗 (Close all figure windows)

%% 尺度参数 (Dimensional Parameters)
a11=135;       % 移动平台上连接点到平台中心的距离 (Distance from connection point to center of moving platform)
a22=135;       % 同上，不同连接点 (Same as above, different connection point)
a33=135;       % 同上，不同连接点 (Same as above, different connection point)
b11=570;       % 固定平台上连接点到平台中心的距离 (Distance from connection point to center of fixed platform)
b22=320;       % 同上，不同连接点 (Same as above, different connection point)
b33=320;       % 同上，不同连接点 (Same as above, different connection point)
e=345;         % 中心柱长度 (Length of central column)
dv=120;        % 可能是运动范围参数 (Possibly a motion range parameter)
pi=3.1415926;  % 圆周率 (Pi value)

%% 位置逆解 (Inverse Kinematics)
% 已知末端执行器位置，求各驱动杆的长度
% (Given end-effector position, calculate the length of each actuating rod)

xp=10          % 末端执行器位置X坐标 (X-coordinate of end-effector position)
yp=-20         % 末端执行器位置Y坐标 (Y-coordinate of end-effector position)
zp=1234        % 末端执行器位置Z坐标 (Z-coordinate of end-effector position)
rp=[xp;yp;zp]; % 末端执行器位置向量 (End-effector position vector)

q4=sqrt(xp*xp+yp*yp+zp*zp)-e;  % 计算q4，为中心连接杆的长度 (Calculate q4, length of central connecting rod)
theta2=asin(xp/(q4+e));        % 计算theta2，机构俯仰角 (Calculate theta2, pitch angle of mechanism)
theta1=atan2((-yp/(q4+e))/(cos(theta2)),(zp/(q4+e))/(cos(theta2)));  % 计算theta1，机构偏航角 (Calculate theta1, yaw angle of mechanism)

% 旋转矩阵R4，描述移动平台相对于固定平台的姿态
% (Rotation matrix R4, describing the attitude of moving platform relative to fixed platform)
R4=[cos(theta2)  0 sin(theta2);
    sin(theta1)*sin(theta2)   cos(theta1)    -sin(theta1)*cos(theta2);
   -cos(theta1)*sin(theta2)  sin(theta1)  cos(theta1)*cos(theta2)];

% 固定平台上三个连接点的位置矢量（在固定坐标系中）
% (Position vectors of three connection points on fixed platform in fixed coordinate system)
a10=[0;-a11;0];
a20=[a22;0;0];
a30=[-a33;0;0];

% 移动平台上三个连接点的位置矢量（在固定坐标系中）
% (Position vectors of three connection points on moving platform in fixed coordinate system)
b1=[0;-b11;0];
b2=[b22;0;0];
b3=[-b33;0;0];

% 将移动平台上的连接点位置变换到当前姿态
% (Transform connection points on moving platform to current attitude)
a1=R4*a10;
a2=R4*a20;
a3=R4*a30;

% 计算方向向量w4，表示中心连接杆的空间方向
% (Calculate direction vector w4, representing spatial direction of central connecting rod)
w4=[sin(theta2);-sin(theta1)*cos(theta2);cos(theta1)*cos(theta2)];

% 计算三条驱动杆的长度（q1,q2,q3）
% (Calculate the lengths of three driving rods (q1, q2, q3))
q1=norm(rp-b1+a1-e*w4)  % 第一条驱动杆长度 (Length of first driving rod)
q2=norm(rp-b2+a2-e*w4)  % 第二条驱动杆长度 (Length of second driving rod)
q3=norm(rp-b3+a3-e*w4)  % 第三条驱动杆长度 (Length of third driving rod)
qi=[q1;q2;q3];          % 驱动杆长度向量 (Vector of driving rod lengths)
 
%% 位置正解 (Forward Kinematics)
%% 迭代法 (Iterative Method)
% 已知三个驱动杆长度(q1,q2,q3)，求末端执行器位置(xp,yp,zp)
% (Given three driving rod lengths (q1,q2,q3), calculate end-effector position (xp,yp,zp))

q4=1000;       % 初始猜测值 - 中心连接杆长度 (Initial guess - length of central connecting rod)
theta1=0.5;    % 初始猜测值 - 偏航角 (Initial guess - yaw angle)
theta2=0.5;    % 初始猜测值 - 俯仰角 (Initial guess - pitch angle)

% 计算约束方程F,G,H的初始值
% (Calculate initial values of constraint equations F, G, H)
F=a11*a11+b11*b11+q4*q4-q1*q1-2*a11*b11*cos(theta1)-2*b11*q4*sin(theta1)*cos(theta2);
G=a22*a22+b22*b22+q4*q4-q2*q2-2*b22*q4*sin(theta2)-2*a22*b22*cos(theta2);
H=a33*a33+b33*b33+q4*q4-q3*q3+2*b33*q4*sin(theta2)-2*a33*b33*cos(theta2);

i=0;  % 迭代计数器 (Iteration counter)

% 牛顿-拉夫森迭代求解方程组
% (Newton-Raphson iteration to solve the system of equations)
while(abs(F)>=0.0001|abs(G)>=0.0001|abs(H)>=0.0001)  % 当任一方程误差超过阈值时继续迭代 (Continue iteration when any equation error exceeds threshold)
    % 计算各约束方程对各变量的偏导数
    % (Calculate partial derivatives of each constraint equation with respect to each variable)
    dF1=2*q4-2*b11*sin(theta1)*cos(theta2);         % dF/dq4
    dF2=-2*b11*q4*cos(theta1)*cos(theta2)+2*a11*b11*sin(theta1);  % dF/dtheta1
    dF3=2*q4*b11*sin(theta1)*sin(theta2);           % dF/dtheta2
    
    dG1=2*q4-2*b22*sin(theta2);                     % dG/dq4
    dG2=0;                                          % dG/dtheta1
    dG3=-2*b22*q4*cos(theta2)+2*a22*b22*sin(theta2); % dG/dtheta2
    
    dH1=2*q4+2*b33*sin(theta2);                     % dH/dq4
    dH2=0;                                          % dH/dtheta1
    dH3=2*b33*q4*cos(theta2)+2*a33*b33*sin(theta2); % dH/dtheta2
    
    % 利用克莱默法则计算变量的增量
    % (Using Cramer's rule to calculate variable increments)
    dq4=(G*(dF2*dH3 - dF3*dH2))/(dF1*dG2*dH3 - dF1*dG3*dH2 - dF2*dG1*dH3 + dF2*dG3*dH1 + dF3*dG1*dH2 - dF3*dG2*dH1) - (F*(dG2*dH3 - dG3*dH2))/(dF1*dG2*dH3 - dF1*dG3*dH2 - dF2*dG1*dH3 + dF2*dG3*dH1 + dF3*dG1*dH2 - dF3*dG2*dH1) - (H*(dF2*dG3 - dF3*dG2))/(dF1*dG2*dH3 - dF1*dG3*dH2 - dF2*dG1*dH3 + dF2*dG3*dH1 + dF3*dG1*dH2 - dF3*dG2*dH1);
    dtheta1=(F*(dG1*dH3 - dG3*dH1))/(dF1*dG2*dH3 - dF1*dG3*dH2 - dF2*dG1*dH3 + dF2*dG3*dH1 + dF3*dG1*dH2 - dF3*dG2*dH1) - (G*(dF1*dH3 - dF3*dH1))/(dF1*dG2*dH3 - dF1*dG3*dH2 - dF2*dG1*dH3 + dF2*dG3*dH1 + dF3*dG1*dH2 - dF3*dG2*dH1) + (H*(dF1*dG3 - dF3*dG1))/(dF1*dG2*dH3 - dF1*dG3*dH2 - dF2*dG1*dH3 + dF2*dG3*dH1 + dF3*dG1*dH2 - dF3*dG2*dH1);
    dtheta2=(G*(dF1*dH2 - dF2*dH1))/(dF1*dG2*dH3 - dF1*dG3*dH2 - dF2*dG1*dH3 + dF2*dG3*dH1 + dF3*dG1*dH2 - dF3*dG2*dH1) - (F*(dG1*dH2 - dG2*dH1))/(dF1*dG2*dH3 - dF1*dG3*dH2 - dF2*dG1*dH3 + dF2*dG3*dH1 + dF3*dG1*dH2 - dF3*dG2*dH1) - (H*(dF1*dG2 - dF2*dG1))/(dF1*dG2*dH3 - dF1*dG3*dH2 - dF2*dG1*dH3 + dF2*dG3*dH1 + dF3*dG1*dH2 - dF3*dG2*dH1);
    
    % 更新变量值
    % (Update variable values)
    q4=q4+dq4;
    theta1=theta1+dtheta1;
    theta2=theta2+dtheta2;
    i=i+1;  % 迭代次数加1 (Increment iteration count)
    
    % 用更新后的变量重新计算约束方程的值
    % (Recalculate values of constraint equations with updated variables)
    F=a11*a11+b11*b11+q4*q4-q1*q1-2*a11*b11*cos(theta1)-2*b11*q4*sin(theta1)*cos(theta2);
    G=a22*a22+b22*b22+q4*q4-q2*q2-2*b22*q4*sin(theta2)-2*a22*b22*cos(theta2);
    H=a33*a33+b33*b33+q4*q4-q3*q3+2*b33*q4*sin(theta2)-2*a33*b33*cos(theta2);
end

% 计算并输出末端执行器的位置坐标
% (Calculate and output position coordinates of end-effector)
xp=(q4+e)*sin(theta2)
yp=(q4+e)*(-sin(theta1)*cos(theta2))
zp=(q4+e)*cos(theta1)*cos(theta2)