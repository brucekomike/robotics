# 曲柄滑块
机构具备多个内容
```matlab
clear all
clc
pkg load symbolic
% pkg install -forge symbolic

% 函数封装
Rx = @(r1, r2)[r1; r2]; % 位置变量
qx = @(R, a)[R; a]; % 构型变量
Ax = @(a)[cos(a) -sin(a); sin(a) cos(a)];
ux = @(L, factor)[factor * L; 0];
Ux = @(R, A, L, factor)R + A * ux(L, factor);

syms r11 r12 a1 r21 r22 a2 r31 r32 a3 L1 L2
syms dr11 dr12 da1 dr21 dr22 da2 dr31 dr32 da3
syms ddr11 ddr12 dda1 ddr21 ddr22 dda2 ddr31 ddr32 dda3

R1 = Rx(r11, r12);
A1 = Ax(a1);
fun1 = Ux(R1, A1, L1, 0.5);

R2 = Rx(r21, r22);
R3 = Rx(r31, r32);
A2 = Ax(a2);
A3 = Ax(a3);
fun2 = Ux(R1, A1, L1, 0.5) - Ux(R2, A2, L2, -0.5);

dR1 = Rx(dr11, dr12);
dR2 = Rx(dr21, dr22);
dR3 = Rx(dr31, dr32);
q1 = qx(R1, a1);
q2 = qx(R2, a2);
q3 = qx(R3, a3);
q = [q1; q2; q3];

U2b = Ux(R2, A2, L2, 0.5);
U3b = Ux(R3, A3, 0, 0);
fun3=U2b-U3b;

fun4=[r32;a3];

fun_c=[fun1;fun2;fun3;fun4]

qi=a1;
qd=[R1;R2;a2;R3;a3]

dqi=da1;
dqd=[dR1;dR2;da2;dR3;da3];

qn=[qi;qd];
dqn=[dqi;dqd];
ddqi=dda1;

yd1=solve(fun_c,qd)
/*
r11=yd1.r11
r12=yd1.r12
r21=yd1.r21
r22=yd1.r22
a2=yd1.a2
r31=yd1.r31
r32=yd1.r32
a3=yd1.a3
*/

Cqd=jacobian(fun_c,qd);
Qci=jacobian(fun_c,qi);


```