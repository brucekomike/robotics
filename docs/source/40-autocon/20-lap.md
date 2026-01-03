# 第二章 拉普拉斯变换

## 复变量与复变函数
复变量

$$
s = \sigma + j\omega
$$

复变函数

$$ 
G (s) = U (\sigma, \omega) + jV(\sigma,\omega)
$$

欧拉公式

$$
cos \theta = 
\frac{1}{2}(e^{j\theta}+ e^{−j\theta})
$$

$$
sin \theta =
 \frac{1}{2j}(e^{j\theta}−e^{−j\theta})
$$

## 拉普拉斯变换概念
### 原理
- 时域中的微分方程变换成复数域中的代数方程
- 傅里叶变换：
  - 区域内连续、有限个第一类间断点
  - 有限个极大极小值
  - 绝对可积
- 拉普拉斯变换：
$$
\varphi(t)\rightarrow
\varphi(t)\mu(t)e^{-\beta t}\\
\Downarrow\\
\int_{-\infty}^{+\infty} \varphi(t) u(t)
e^{-\beta t} e^{-j\omega t} \text{d}t
$$
### 拉普拉斯变化 拉氏变换
- $f(t) = \varphi(t)u(t)$ 为时间 $t$ 的函数,
并且当$t<0$时$f(t)=0$
- $s=\beta +j \omega $为复变量
- $L$ 为运算符号，放在某个时间之前，表示该时间函数用拉氏积分进行变换；
- $F(s)$为时间函数$f(t)$的拉氏变换
- $F(s)$成为像函数，$f(t)$成为原函数

$$ F(s) = L[f(t)]
= \int_{0}^{\infty }e^{-st}\text{d}t[f(t)]
= \int_{0}^{\infty }f(t)e^{-st}\text{d}t
$$

### 拉氏逆变换
拉氏逆变换：
从拉氏变换$F(s)$求时间函数$f(t)$的过程。

$$
L^{-1}[F(s)]=f(t)
= \frac{1}{2\pi j}
\int_{\beta-j\infty }^{\beta*j\infty }F(s)e^{st}\text{d}s (t\geq0)
$$
### 存在定理
1. 当$t\geq 0$时，$f(t)$在任一有限区间上是分段连续的；
2. 当$t\rightarrow \infty $时$f(t)$的增长速度不超过某一指数函数，即存在常数$M>0$,及$c\geq 0$，使得：
$$|f(t)|\leq M e^{ct},0\leq t<0$$
$$\int_{0}^{\infty }f(t)e^{-st}\text{d}t,f(t)=\varphi(t)u(t)$$
### 常用的拉氏变换
#### 1指数函数
$$
f(t) =
\begin{cases}
0 & t < 0 \\
Ae^{-\alpha t} & t \geq 0
\end{cases}
$$

$A$和$\alpha$为常数

$$L[f(t)]=L[Ae^{-\alpha t}]\\
=\int_{0}^{+\infty }Ae^{-e t}e^{-s t}\text{d}t=A\int_{0}^{+\infty }e^{-(\alpha+s)t}\text{d}t=\frac{A}{s+\alpha}$$

其中，$Re(s+\alpha)>0$,即$Re(s)>-\alpha$

指数函数在复平面内将产生一个极点。

#### 2阶跃函数

$$
f(t) =
\begin{cases}
0 & t < 0 \\
A &  t \geq 0
\end{cases}
$$

$A$为常数

$$
L[f(t)]=L[A]=\int_{0}^{+\infty }Ae^{-st}\text{d}t=\frac{A}{s}$$

> 单位阶跃函数 $A=1 \rightarrow u(t) $ 
> $$
f(t) =
\begin{cases}
0 & t < 0 \\
1 &  t \geq 0
\end{cases} $$
> $$
L[u(t)]=\frac{1}{s}$$

#### 3斜坡函数

$$
f(t) =
\begin{cases}
0 & t < 0 \\
At &  t \geq 0
\end{cases}
$$

$A$为常数

$$
L[At] = \int_{0}^{+\infty} Ate^{-st} dt \\
= At \frac{e^{-st}}{-s}\bigg|_{0}^{+\infty} - \int_{0}^{+\infty} \frac{Ae^{-st}}{-s} dt \\
= \frac{A}{s} \int_{0}^{+\infty} e^{-st} dt = \frac{A}{s^2}
$$

> 单位斜坡信号 $A=1 \rightarrow u(t)$ 
> $$
f(t) =
\begin{cases}
0 & t < 0 \\
t &  t \geq 0
\end{cases} $$
> $$
L[r(t)]=\frac{1}{s^2}$$
#### 4正弦

$$
f(t) =
\begin{cases}
0 & t < 0 \\
A\sin{\omega t} &  t \geq 0
\end{cases}
$$

$A,\omega$为常数,
> 欧拉公式：$\sin{\omega t}=\frac{1}{2j}(e^{j\omega t}-e^{-j\omega t})$
$$
L[A\sin{\omega t}]=\frac{A}{2j}\int_{0}^{+\infty }(e^{j\omega t}-e^{-j\omega t})e^{-st}\text{d}t\\
=\frac{A}{2j}\frac{1}{s-j\omega}-\frac{A}{2j}\frac{1}{s+j\omega}\\
=\frac{A\omega}{s^2+\omega^2}
$$
#### 5脉动函数

$$
f(t) =
\begin{cases}
\frac{A}{t_0} & 0<t<t_0 \\
0 &  0<t,t_0<t
\end{cases}
$$

$A,t_0$为常数

#### 6脉冲函数
脉冲函数是脉动函数的一种特殊极限情况
$$
f(t) =
\begin{cases}
\frac{A}{t_0} & 0<t<\Delta \\
0 &  0<t,\Delta<t
\end{cases}
$$

$A$为常数

$$
L[g(t)]=
\lim_{\Delta\rightarrow0}\bigg[\frac{A}{\Delta s}(a-e^{-s\Delta})\bigg]\\
=\lim_{\Delta\rightarrow0}\frac{\frac{\text{d}}{\text{d}\Delta}\big[A(1-e^{-s\Delta})\big]}
{\frac{\text{d}}{\text{d}\Delta}(\Delta s)}\\
=\frac{As}{s}=A
$$

$A=1,\varepsilon\rightarrow0$时，称为单位脉冲信号或狄拉克(Disac)函数$\delta(t)$
>$L[\delta(t)]=1$
#### 7加速度函数

$$
f(t) =
\begin{cases}
0 & t < 0 \\
At^2 &  t \geq 0
\end{cases}
$$

$A$为常数

$$
L[At^2]=\int_{0}^{+\infty }At^2e^{-st}\text{d}t\\
=\frac{A}{s}\bigg[t^2e^{-st}\Big|_0^{+\infty }-2\int_0^{+\infty }te^{-st}\text{d}t\bigg]\\
=2A\frac{1}{s^2}
$$

> 单位加速度信号 $A=1/2 \rightarrow a(t)$ 
> $$
f(t) =
\begin{cases}
0 & t < 0 \\
\frac{1}{2}t^2 &  t \geq 0
\end{cases} $$
> $$
L\bigg[\frac{1}{2}t^2\bigg]=\frac{1}{s^3}$$

## 拉普拉斯变换性质
### 1线性性质

$$
L[K_1f(t)\pm K_2f_2(t)]\\
=\int_{0}^{+\infty }[K_1f_1(t)\pm K_2f_2(t)]e^{-st}\text{d}t\\
=\int_{0}^{+\infty }K_1f_1(t)e^{-st}\text{d}t\pm\int_{0}^{+\infty }K_2f_2(t)e^{-st}\text{d}t\\
=K_1F_1(s)\pm K_2F_2(s)
$$

线性性质也称叠加性质，即各函数之和的拉氏变换等于各函数拉氏变换之和
### 2微分性质

$$
\text{若} L[f(t)]=F(s),\text{则} L \big[\frac{\text{d}f(t)}{\text{d}t}\big]=sF(s)-f(0)
$$

$f(0)$是$f(t)$在$t=0$时的初始值。

当 $f(t)$ 在 $t=0$ 处具有间断点时，$f(0_+)$ 和 $f(0_-)$ 之间的差别很重要，因为此时 $df(t)/dt$ 在 $t=0$ 处将包含一个脉冲函数 $\delta(t)$。即 $f(0_+) \neq f(0_-)$

则 $L_{+}\left[\frac{d f(t)}{d t}\right]=s F(s)-f(0_{+})$，$L_{-}\left[\frac{d f(t)}{d t}\right]=s F(s)-f(0_{-})$
### 3积分性质

### 4位移性质
### 5延迟性质
### 6尺度变换

### 7初值定理、终值定理

#### 初值定理
#### 终值定理


## 拉普拉斯逆变换
### 1只含不同极点的情况
### 2含共轭复极点的情况
### 3. 含多重极点的情况

## 卷积
### 1. 概念