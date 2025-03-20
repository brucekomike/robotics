# 雅可比
```matlab
clear all
clc

syms x1 x2 x3;
fun1=x1+x3^2;
fun2=x2-x1;
fun3=x1^2-x3;
fun=[fun1;fun2;fun3];
x=[x1;x2;x3];
jacobian(fun,x)
```