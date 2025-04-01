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