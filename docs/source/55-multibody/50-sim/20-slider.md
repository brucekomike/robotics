# 曲柄滑块
机构具备多个内容
## slider 函数
````{tab} slider
```matlab
% slider.m
function yd=slider(L1,L2)
%pkg load symbolic
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
qd=[R1;R2;a2;R3;a3];

dqi=da1;
dqd=[dR1;dR2;da2;dR3;da3];

qn=[qi;qd];
dqn=[dqi;dqd];
ddqi=dda1;

yd1=solve(fun_c,qd);
%{
r11=yd1.r11
r12=yd1.r12
r21=yd1.r21
r22=yd1.r22
a2=yd1.a2
r31=yd1.r31
r32=yd1.r32
a3=yd1.a3
%}

Cqd=jacobian(fun_c,qd);
Cqi=jacobian(fun_c,qi);

Qc =jacobian(jacobian(fun_c,qn)*dqn, qn)*dqn;
Qdi=-simplify(inv(Cqd)*Cqi);

yd2=Cqi*dqi;

yd3=Cqi*ddqi+inv(Cqd)*Qc;

%yd={yd1;yd2;yd3};
yd=yd1;
end
```
````
````{tab} slider2
```matlab
function yd=slider(L1,L2)
% clear
% clc
syms r11 r12 a1 r21 r22 a2 r31 r32 a3 
R1=[r11;r12];R2=[r21;r22];R3=[r31;r32];
q1=[R1;a1];q2=[R2;a2];q3=[R3;a3];
q=[q1;q2;q3];
A1=[cos(a1) -sin(a1);sin(a1) cos(a1)];
A2=[cos(a2) -sin(a2);sin(a2) cos(a2)];
A3=[cos(a3) -sin(a3);sin(a3) cos(a3)];
u1_o=[-0.5*L1;0];
u10=R1+A1*u1_o;
fun1=R1+A1*u1_o;

u1_a=[0.5*L1;0];
u2_a=[-0.5*L2;0];
u1a=R1+A1*u1_a;
u2a=R2+A2*u2_a;
fun2=u1a-u2a;
u2_b=[0.5*L2;0];
u3_b=[0;0];
u2b=R2+A2*u2_b;
u3b=R3+A3*u3_b;
fun3=u2b-u3b;
fun4=[r32;a3];
fun_c=[fun1;fun2;fun3;fun4];
qi=a1;
qd=[R1;R2;a2;R3;a3];
% dqi=da1;
% dqd=[dR1;dR2;da2;dR3;da3];
% qn=[qi;qd];
% dqn=[dqi;dqd];
% ddqi=dda1;
%% 位置正解
yd1=solve(fun_c,qd);
% r11=yd1.r11;
% r12=yd1.r12;
% r21=yd1.r21;
% r22=yd1.r22;
% a2=yd1.a2;
% r31=yd1.r31;
% r32=yd1.r32;
% a3=yd1.a3
yd=yd1;
end

```
````

## 主函数
```matlab
clear all
clc

pkg load symbolic
% pkg install -forge symbolic

syms r11 r12 a1 r21 r22 a2 r31 r32 a3 
syms dr11 dr12 da1 dr21 dr22 da2 dr31 dr32 da3
syms ddr11 ddr12 dda1 ddr21 ddr22 dda2 ddr31 ddr32 dda3
L1=0.5; L2=1.5;
yd=slider(L1,L2);
% 位移
qd1=yd{1,1};
qd1=qd1{1,1};
r11=qd1.r11
r12=qd1.r12
a1=qd1.a1
r21=qd1.r21
r22=qd1.r22
a2=qd1.a2
r31=qd1.r32
r32=qd1.r32

R1=[r11;r12];
R2=[r21;r22];
R3=[r31;r32];
%velocity
dr11=diff(r11,a1)*da1;
dr12=diff(r12,a1)*da1;
dr21=diff(r21,a1)*da1;
dr22=diff(r22,a1)*da1;
dr31=diff(r31,a1)*da1;
dr32=diff(r32,a1)*da1;
da3=diff(a3,a1)*da1;

% acceleration
function ddd=ddiff(dr11,a1,da1,dda1)
  ddd=diff(dr11,a1)*da1+diff(dr11,da1)*dda1;
end

ddr11=ddiff(r11,a1,da1,dda1);
ddr12=ddiff(r12,a1,da1,dda1);
dda2=ddiff(a2,a1,da1,dda1);
ddr21=ddiff(r21,a1,da1,dda1);
ddr22=ddiff(r22,a1,da1,dda1);
ddr31=ddiff(r31,a1,da1,dda1);
ddr32=ddiff(r32,a1,da1,dda1);
dda3=ddiff(a3,a1,da1,dda1);

```