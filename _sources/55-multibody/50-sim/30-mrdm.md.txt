# 移动机器人运动学建模
```matlab
clear all
clc

pkg load symbolic
% pkg install -forge symbolic

% 函数封装
Rx= @(r1, r2)[r1;r2]; % 位置变量
qx= @(R, a) [R; a];% 构型变量
Ax= @(a) [[cos(a) -sin(a)];[sin(a) cos(a)]];
u_x = @(L) [(1/2)*L; 0];
rx = @(r1, r2, a, L) Rx(r1, r2) + Ax(a)*u_x(L);

syms r31 r32 a3 L3 r41 r42 a4 L4
syms r11 r12 a1 L1 r21 r22 a2 L2
r1c = rx(r11, r12, a1, L1);
r2c = rx(r21, r22, a2, L2);
r3c = rx(r31, r32, a3, L3);
r4c = rx(r41, r42, a4, L4);

fun4= r3c-r4c;
fun4=simplify(fun4);
fun3= r3c-r2c;
fun3 =simplify(fun3);
fun2= r3c-r2c;
fun2 =simplify(fun2);
fun1= r1c-r2c;
fun1 =simplify(fun1);
fun=[fun1;fun2;fun3;fun4];
fun=simplify(fun);
```