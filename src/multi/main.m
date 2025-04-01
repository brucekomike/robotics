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
