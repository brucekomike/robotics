# 第二章 拉普拉斯变换

## 2.1 拉普拉斯变换概念

### 引言：复变量与复变函数

*   **复变量**: $s = \sigma + j\omega$
*   **复变函数**: $G(s) = U(\sigma, \omega) + jV(\sigma, \omega)$
*   **欧拉公式**:
    $$ \cos \theta = \frac{1}{2}(e^{j\theta} + e^{-j\theta}) $$
    $$ \sin \theta = \frac{1}{2j}(e^{j\theta} - e^{-j\theta}) $$

### 提出问题

*   将时域中的微分方程变换成复数域中的代数方程。
*   **傅里叶变换**：
    要求函数在区域内连续，有有限个第一类间断点，有限个极大极小值，且绝对可积。
*   **拉普拉斯变换**: $\varphi(t) \rightarrow \varphi(t)u(t)e^{-\beta t} (\beta > 0)$
    $$ \int_{-\infty}^{+\infty} \varphi(t)u(t)e^{-\beta t} e^{-j\omega t} dt $$

### 拉普拉斯（Laplace）变换 — 拉氏变换

*   $f(t) = \varphi(t)u(t)$ 为时间的函数，并且当 $t < 0$ 时 $f(t) = 0$。
*   $s = \beta + j\omega$ 为复变量。
*   $L$ 为运算符号，放在某个时间函数之前，表示该时间函数用拉氏积分进行变换。
*   $F(s)$ 为时间函数 $f(t)$ 的拉氏变换。
    $$ F(s) = L[f(t)] = \int_{0}^{\infty} e^{-st} dt [f(t)] = \int_{0}^{\infty} f(t)e^{-st} dt $$
    *   $F(s)$ 称为 **象函数**。
    *   $f(t)$ 称为 **原函数**。

### 拉氏逆变换

从拉氏变换 $F(s)$ 求时间函数 $f(t)$ 的过程。
$$ L^{-1}[F(s)] = f(t) = \frac{1}{2\pi j}\int_{\beta - j\infty}^{\beta + j\infty} F(s)e^{st} ds $$
$t \ge 0$

### 存在定理

1.  当 $t \ge 0$ 时，$f(t)$ 在任一有限区间上是分段连续的。
2.  当 $t \rightarrow \infty$ 时，$f(t)$ 的增长速度不超过某一指数函数，即存在常数 $M>0$ 及 $c \ge 0$，使得
    $$ |f(t)| \le Me^{ct}, \quad 0 \le t < \infty $$
    $$ \int_{0}^{+\infty} f(t)e^{-st} dt < \infty \quad f(t) = \varphi(t)u(t) $$

### 常用函数的拉氏变换

#### (1) 指数函数

$$ f(t) = \begin{cases} 0 & t < 0 \\ Ae^{-\alpha t} & t \ge 0 \end{cases} $$
$A$ 和 $\alpha$ 为常数。
$$ L[f(t)] = L[Ae^{-\alpha t}] = \int_{0}^{\infty} Ae^{-\alpha t} e^{-st} dt = A \int_{0}^{\infty} e^{-(a+s)t} dt = \frac{A}{s+\alpha} $$
其中，$\text{Re}(s+\alpha) > 0$，即 $\text{Re}(s) > -\alpha$。
指数函数在复平面内将产生一个极点。

#### (2) 阶跃函数

$$ f(t) = \begin{cases} 0 & t < 0 \\ A & t \ge 0 \end{cases} $$
$A$ 为常数。
$$ L[f(t)] = L[A] = \int_{0}^{\infty} Ae^{-st} dt = \frac{A}{s} $$
当 $A=1$ 时，单位阶跃信号 $u(t)$：
$$ u(t) = \begin{cases} 0 & t < 0 \\ 1 & t \ge 0 \end{cases} $$
$$ L[u(t)] = \frac{1}{s} $$

#### (3) 斜坡函数

$$ f(t) = \begin{cases} 0 & t < 0 \\ At & t \ge 0 \end{cases} $$
$A$ 为常数。
$$ L[At] = \int_{0}^{\infty} Ate^{-st} dt = A \int_{0}^{\infty} te^{-st} dt = A \left[ -\frac{te^{-st}}{s} \right]_0^{\infty} - A \int_{0}^{\infty} (-\frac{e^{-st}}{s}) dt = \frac{A}{s} \int_{0}^{\infty} e^{-st} dt = \frac{A}{s^2} $$
当 $A=1$ 时，单位斜坡信号 $r(t)$：
$$ f(t) = \begin{cases} 0 & t < 0 \\ t & t \ge 0 \end{cases} $$
$$ L[r(t)] = \frac{1}{s^2} $$

#### (4) 正弦函数

$$ f(t) = \begin{cases} 0 & t < 0 \\ A \sin \omega t & t \ge 0 \end{cases} $$
$A$ 和 $\omega$ 为常数。
欧拉公式： $\sin \omega t = \frac{1}{2j}(e^{j\omega t} - e^{-j\omega t})$
$$ L[A \sin \omega t] = \frac{A}{2j} \int_{0}^{\infty} (e^{j\omega t} - e^{-j\omega t})e^{-st} dt = \frac{A}{2j} \left[ \frac{1}{s-j\omega} - \frac{1}{s+j\omega} \right] = \frac{A}{2j} \frac{s+j\omega - (s-j\omega)}{s^2 + \omega^2} = \frac{A\omega}{s^2 + \omega^2} $$
$$ L[A \cos \omega t] = \frac{As}{s^2 + \omega^2} $$
其中，$\text{Re}(s) > 0$。

#### (5) 脉动函数

$$ f(t) = \begin{cases} \frac{A}{t_0} & 0 < t < t_0 \\ 0 & t < 0, t_0 < t \end{cases} $$
$A$ 和 $t_0$ 为常数。

#### (6) 脉冲函数

脉冲函数是脉动函数的一种特殊极限情况。
$$ g(t) = \begin{cases} \lim_{\Delta \to 0} \frac{A}{\Delta} & 0 < t < \Delta \\ 0 & t < 0, \Delta < t \end{cases} $$
$$ L[g(t)] = \lim_{\Delta \to 0} L\left[\frac{A}{\Delta}(1 - e^{-s\Delta})\right] = \lim_{\Delta \to 0} \frac{d}{d\Delta}\left[\frac{A(1-e^{-s\Delta})}{\Delta s}\right] = \lim_{\Delta \to 0} \frac{A}{s} = \frac{A}{s} $$
当 $A=1, \epsilon \to 0$ 时，称为单位脉冲信号或狄拉克（Dirac）函数 $\delta(t)$。
$$ L[\delta(t)] = 1 $$

#### (7) 加速度函数

$$ f(t) = \begin{cases} At^2 & t \ge 0 \\ 0 & t < 0 \end{cases} $$
$A$ 为常数。
$$ L[At^2] = \int_{0}^{\infty} At^2 e^{-st} dt = \frac{A}{s} \left[ t^2 e^{-st} \right]_0^{\infty} - 2 \int_{0}^{\infty} te^{-st} dt = 2A \frac{1}{s^3} $$
当 $A=1/2$ 时，单位加速度信号 $a(t)$：
$$ a(t) = \begin{cases} 0 & t < 0 \\ \frac{1}{2} t^2 & t \ge 0 \end{cases} $$
$$ L\left[\frac{1}{2}t^2\right] = \frac{1}{s^3} $$

## 2.2 拉普拉斯变换的性质

### 1. 线性性质

线性性质也称叠加性质，即各函数之和的拉氏变换等于各函数拉氏变换之和。当函数乘以 $K$ 时，其变换式也乘以相同的常数 $K$。

若 $L[f_1(t)] = F_1(s)$ 且 $L[f_2(t)] = F_2(s)$，$K_1, K_2$ 为常数，
则 $L[K_1 f_1(t) \pm K_2 f_2(t)] = K_1 F_1(s) \pm K_2 F_2(s)$。

**证明**:
$$ L[K_1 f_1(t) \pm K_2 f_2(t)] = \int_{0}^{\infty} [K_1 f_1(t) \pm K_2 f_2(t)] e^{-st} dt $$
$$ = K_1 \int_{0}^{\infty} f_1(t)e^{-st} dt \pm K_2 \int_{0}^{\infty} f_2(t)e^{-st} dt $$
$$ = K_1 F_1(s) \pm K_2 F_2(s) $$
这个性质表明了各函数线性组合的拉氏变换等于各函数拉氏变换的线性组合。

**例 2.1** 求 $f(t) = \sin \omega t$ 的拉氏变换。
**解**: 根据欧拉公式
$$ f(t) = \sin \omega t = \frac{1}{2j}(e^{j\omega t} - e^{-j\omega t}) $$
$$ L[e^{j\omega t}] = \frac{1}{s-j\omega}, \quad L[e^{-j\omega t}] = \frac{1}{s+j\omega} $$
由拉氏变换的线性性质可知
$$ L[\sin \omega t] = \frac{1}{2j} \left[ \frac{1}{s-j\omega} - \frac{1}{s+j\omega} \right] = \frac{1}{2j} \frac{s+j\omega - (s-j\omega)}{s^2 + \omega^2} = \frac{\omega}{s^2 + \omega^2} $$
用同样的方法可求得
$$ L[\cos \omega t] = \frac{s}{s^2 + \omega^2} $$

### 2. 微分性质

若 $L[f(t)] = F(s)$，则 $L\left[\frac{df(t)}{dt}\right] = sF(s) - f(0)$。
$f(0)$ 是 $f(t)$ 在 $t=0$ 时的初始值。
对于给定的时间函数，$f(0_+)$ 和 $f(0_-)$ 的值可能相同，也可能不同。
当 $f(t)$ 在 $t=0$ 处具有间断点时，$f(0_+)$ 和 $f(0_-)$ 之间的差别很重要，因为此时 $df(t)/dt$ 在 $t=0$ 处将包含一个脉冲函数 $\delta(t)$。即 $f(0_+) \ne f(0_-)$。
则 $L_+\left[\frac{df(t)}{dt}\right] = sF(s) - f(0_+)$， $L_-\left[\frac{df(t)}{dt}\right] = sF(s) - f(0_-)$。

**证明**:
$$ L\left[\frac{df(t)}{dt}\right] = \int_{0}^{\infty} \left[\frac{df(t)}{dt}\right]e^{-st} dt $$
$$ = \left[f(t)e^{-st}\right]_0^{\infty} - \int_{0}^{\infty} f(t)(-s)e^{-st} dt $$
$$ = s \int_{0}^{\infty} f(t)e^{-st} dt - f(0) = sL[f(t)] - f(0) = sF(s) - f(0) $$
这个性质表明了一个时间函数 $f(t)$ 求导后取拉氏变换等于这个函数的拉氏变换乘以 $s$，再减去这个函数的初始值 $f(0)$。
一阶导数的微分性质可以推广到高阶导数。
$$ L\left[\frac{d^2 f(t)}{dt^2}\right] = s^2F(s) - sf(0) - f'(0) $$
$f'(0)$ 是 $df(t)/dt$ 在 $t=0$ 时的值。

**证明**:
$$ L\left[\frac{d^2 f(t)}{dt^2}\right] = e^{-st} \frac{df(t)}{dt} \Big|_0^{\infty} + s \int_{0}^{\infty} \frac{df(t)}{dt} e^{-st} dt $$
$$ = -f'(0) + s[sF(s) - f(0)] $$
$$ = s^2F(s) - sf(0) - f'(0) $$

导出时间函数 $f(t)$ 的 $n$ 阶导数微分性质
$$ L\left[\frac{d^n f(t)}{dt^n}\right] = s^n F(s) - \sum_{r=0}^{n-1} s^{n-r-1} f^{(r)}(0) $$
$f^{(r)}(0)$ 是 $r$ 阶导数 $\frac{d^r f(t)}{dt^r}$ 在 $t=0$ 时的值。
为了保证 $f(t)$ 的各阶导数的拉氏变换存在，$\frac{d^n f(t)}{dt^n}$ ($n=1,2,3,\dots$) 必须能够进行拉氏变换。
如果 $f(t)$ 及其各阶导数的所有初始值全都为零，则
$$ L\left[\frac{d^n f(t)}{dt^n}\right] = s^n F(s) $$

**例 2.2** 求余弦函数 $g(t) = \begin{cases} 0 & t < 0 \\ \cos \omega t & t \ge 0 \end{cases}$ 的拉氏变换。
**解**: 由于正弦函数 $f(t) = \begin{cases} 0 & t < 0 \\ \sin \omega t & t \ge 0 \end{cases}$ 的拉氏变换 $F(s) = \frac{\omega}{s^2 + \omega^2}$。
根据拉氏变换的微分性质，可以求得
$$ L[\cos \omega t] = L\left[\frac{1}{\omega}\frac{d}{dt}(\sin \omega t)\right] = \frac{1}{\omega}[sF(s) - f(0)] $$
$$ = \frac{1}{\omega}\left[ s \frac{\omega}{s^2 + \omega^2} - 0 \right] = \frac{s}{s^2 + \omega^2} $$

**例 2.3** 求以下微分方程的拉氏变换，已知其各阶导数的初始值为零。
$$ \frac{d^3 x_0(t)}{dt^3} + 2\frac{d^2 x_0(t)}{dt^2} + 3\frac{dx_0(t)}{dt} + x_0(t) = 2\frac{dx_i(t)}{dt} + x_i(t) $$
**解**: 对上式两端取拉氏变换，得
$$ s^3 X_0(s) + 2s^2 X_0(s) + 3sX_0(s) + X_0(s) = 2sX_i(s) + X_i(s) $$
化简得 $(s^3 + 2s^2 + 3s + 1)X_0(s) = (2s+1)X_i(s)$
传递函数：$G(s) = \frac{X_0(s)}{X_i(s)} = \frac{2s+1}{s^3 + 2s^2 + 3s + 1}$

### 3. 积分性质

若 $L[f(t)] = F(s)$，则 $L\left[\int_{0}^{t} f(\tau) d\tau\right] = \frac{F(s)}{s} + \frac{f^{-1}(0)}{s}$。
$f^{-1}(0)$ 是 $\int f(t) dt$ 在 $t=0$ 的值。
如果 $f(t)$ 在 $t=0$ 处包含一个脉冲函数 $\delta(t)$，即 $f^{-1}(0_+) \ne f^{-1}(0_-)$。
则 $L_+ \left[\int_{0}^{t} f(\tau) d\tau\right] = \frac{F(s)}{s} + \frac{f^{-1}(0_+)}{s}$， $L_- \left[\int_{0}^{t} f(\tau) d\tau\right] = \frac{F(s)}{s} + \frac{f^{-1}(0_-)}{s}$。

**证明**: 借助分部积分法进行积分，得
$$ L\left[\int_{0}^{t} f(\tau) d\tau\right] = \int_{0}^{\infty} \left[\int_{0}^{t} f(\tau) d\tau\right] e^{-st} dt $$
$$ = \left[\int_{0}^{t} f(\tau) d\tau \cdot \frac{e^{-st}}{-s}\right]_0^{\infty} - \int_{0}^{\infty} f(t) \cdot \frac{e^{-st}}{-s} dt $$
$$ = \frac{1}{s} \left[\int_{0}^{t} f(\tau) d\tau\right]_{t=0} + \frac{1}{s} \int_{0}^{\infty} f(t)e^{-st} dt $$
$$ = \frac{f^{-1}(0)}{s} + \frac{F(s)}{s} $$

如果积分的初值为零，则 $L\left[\int_{0}^{t} f(\tau) d\tau\right] = \frac{F(s)}{s}$。
对于多重积分的拉氏变换，有
$$ L\left[\int_{0}^{t}\int_{0}^{\tau} f(\sigma) (d\sigma)^2\right] = \frac{F(s)}{s^2} + \frac{f^{-1}(0)}{s} + \frac{f^{-2}(0)}{s} $$
…
$$ L\left[\int_{0}^{t} \dots \int_{0}^{\tau_n} f(t_1) (dt_1)^n\right] = \frac{F(s)}{s^n} + \frac{f^{-1}(0)}{s^n} + \dots + \frac{f^{(-n)}(0)}{s} $$
$f^{-1}(0), f^{-2}(0), \dots, f^{(-n)}(0)$ 为 $f(t)$ 的各重积分在 $t=0$ 时的值。
如果 $f^{-1}(0) = f^{-2}(0) = \dots = f^{(-n)}(0) = 0$，则有
$$ L\left[\int_{0}^{t} \dots \int_{0}^{\tau_n} f(t_1) (dt_1)^n\right] = \frac{F(s)}{s^n} $$

### 4. 位移性质

若 $L[f(t)] = F(s)$，则 $L[f(t)e^{-at}] = F(s+a)$。
此性质表明：时间函数乘以 $e^{-at}$，相当于变换式在复频域内平移 $a$。

**证明**: 根据拉氏变换的定义，得
$$ L[e^{-at} f(t)] = \int_{0}^{\infty} e^{-at} f(t) e^{-st} dt = \int_{0}^{\infty} f(t) e^{-(s+a)t} dt $$
上式的右方只是在 $F(s)$ 中把 $s$ 换成 $s+a$，所以
$$ L[f(t)e^{-at}] = F(s+a) $$

**例**: 求 $e^{-\alpha t} \sin \omega t$ 和 $e^{-\alpha t} \cos \omega t$ 的拉氏变换。
**解**: 已知 $L[\sin \omega t] = \frac{\omega}{s^2 + \omega^2}$。
由拉氏变换的位移性质：
$$ L[e^{-\alpha t} \sin \omega t] = \frac{\omega}{(s+\alpha)^2 + \omega^2} $$
同理，因为 $L[\cos \omega t] = \frac{s}{s^2 + \omega^2}$，
故有
$$ L[e^{-\alpha t} \cos \omega t] = \frac{s+\alpha}{(s+\alpha)^2 + \omega^2} $$

### 5. 延迟性质

若 $L[f(t)] = F(s)$，则 $L[f(t-t_0)u(t-t_0)] = e^{-st_0} F(s)$。
此性质表明：如下图所示的时间函数 $f(t)u(t)$，若在时间轴上延迟 $t_0$ 得到时间函数 $f(t-t_0)u(t-t_0)$，则它的拉氏变换应乘以 $e^{-st_0}$。
[图像展示了原始函数 $f(t)u(t)$ 和延迟后的函数 $f(t-t_0)u(t-t_0)$]

**证明**:
$$ L[f(t-t_0)u(t-t_0)] = \int_{0}^{\infty} [f(t-t_0)u(t-t_0)] e^{-st} dt = \int_{t_0}^{\infty} f(t-t_0)e^{-st} dt $$
令 $\tau = t-t_0$，则 $t = \tau + t_0$，代入上式得
$$ L[f(t-t_0)u(t-t_0)] = \int_{0}^{\infty} f(\tau)e^{-s(\tau+t_0)} d\tau = e^{-st_0} \int_{0}^{\infty} f(\tau)e^{-s\tau} d\tau = e^{-st_0} F(s) $$

**例**: 求图(a)所示的时间函数的拉氏变换。
[包含三个图：(a) 梯形波形 $f(t)$，(b) $f'(t)$，(c) $f''(t)$]
**解**: 时间函数 $f(t)$ 的一阶、二阶导数，如图(b)、(c)所示。
$f''(t) = 2\delta(t) - 2\delta(t-1) - 2\delta(t-2) + 2\delta(t-3)$
由于 $L[\delta(t)] = 1$，
由延迟和线性性质得
$F_2(s) = L[f''(t)] = 2 - 2e^{-s} - 2e^{-2s} + 2e^{-3s}$
由积分性质得
$F(s) = L[f(t)] = \frac{F_2(s)}{s^2} = \frac{2(1-e^{-s} - e^{-2s} + e^{-3s})}{s^2}$

**例**: 已知 $f_1(t) = e^{-2(t-1)}u(t-1)$，$f_2(t) = e^{-2(t-1)}u(t)$，求 $f_1(t) + f_2(t)$ 的拉氏变换。
**解**: 因为 $L[e^{-2t}u(t)] = \frac{1}{s+2}$。
根据拉氏变换的延迟性质，得
$F_1(s) = L[e^{-2(t-1)}u(t-1)] = e^{-s} \cdot \frac{1}{s+2}$。
又因为 $f_2(t) = e^{-2(t-1)}u(t) = e^2 e^{-2t}u(t)$。
根据拉氏变换的线性性质，得 $F_2(s) = e^2 \cdot \frac{1}{s+2}$。
$L[f_1(t) + f_2(t)] = F_1(s) + F_2(s) = \frac{e^{-s} + e^2}{s+2}$。

### 6. 尺度变换

若 $L[f(t)] = F(s)$，则 $L[f(at)] = \frac{1}{a} F\left(\frac{s}{a}\right)$，$a > 0$。

**证明**:
$$ L[f(at)] = \int_{0}^{\infty} f(at) e^{-st} dt $$
令 $\tau = at$，则 $t = \frac{\tau}{a}$， $dt = \frac{1}{a} d\tau$。
$$ L[f(at)] = \int_{0}^{\infty} f(\tau) e^{-s(\tau/a)} \frac{1}{a} d\tau = \frac{1}{a} \int_{0}^{\infty} f(\tau) e^{-(s/a)\tau} d\tau = \frac{1}{a} F\left(\frac{s}{a}\right) $$

**例**: 已知 $L[f(t)] = F(s)$，若 $a>0, b>0$，求 $L[f(at-b)u(at-b)]$。
**解**: 此题既要用到尺度变换，也要用到延迟性质。

**解法一（不推荐）**:
由延迟性质得 $L[f(t-b)u(t-b)] = F(s)e^{-bs}$。
再由尺度变换，即可求得结果：
$L[f(at-b)u(at-b)] = \frac{1}{a} F\left(\frac{s}{a}\right) e^{-(\frac{s}{a})b} = \frac{1}{a} F\left(\frac{s}{a}\right) e^{-sb/a}$。

**解法二**:
先用尺度变换性质：$L[f(at)u(at)] = \frac{1}{a} F\left(\frac{s}{a}\right)$。
然后由延迟性质求出：
$L\left[f\left(at - \frac{b}{a}a\right)u\left(at - \frac{b}{a}a\right)\right] = L\left[f\left(a(t-\frac{b}{a})\right)u\left(a(t-\frac{b}{a})\right)\right]$
令 $t' = t - \frac{b}{a}$，则 $at' = a(t-\frac{b}{a}) = at - b$。
$L[f(at')u(at')] = \frac{1}{a} F\left(\frac{s}{a}\right)$。
将 $t'$ 延迟 $\frac{b}{a}$，即 $t \rightarrow t - \frac{b}{a}$，对应的 $s$ 变为 $s e^{-s(b/a)}$。
$L[f(a(t-\frac{b}{a}))u(a(t-\frac{b}{a}))] = \frac{1}{a} F\left(\frac{s}{a}\right) e^{-s(b/a)}$。
即 $L[f(at-b)u(at-b)] = \frac{1}{a} F\left(\frac{s}{a}\right) e^{-sb/a}$。
两种解法其结果一致。

### 7. 初值定理、终值定理

#### (1) 初值定理

**证明**: 根据拉氏变换的微分性质，有
$$ L\left[\frac{d}{dt}f(t)\right] = \int_{0}^{\infty} \frac{d}{dt}f(t)e^{-st} dt = sF(s) - f(0) $$
令 $s \rightarrow \infty$ 取极限得
$$ \lim_{s\to\infty} \int_{0}^{\infty} \frac{d}{dt}f(t)e^{-st} dt = \lim_{s\to\infty} [sF(s) - f(0)] $$
在时间区间 $[0_+, \infty)$ 内，$\lim_{s\to\infty} e^{-st} = 0$。
因此，
$$ \lim_{s\to\infty} \int_{0}^{\infty} \frac{d}{dt}f(t)e^{-st} dt = \int_{0+}^{\infty} \frac{d}{dt}f(t) \lim_{s\to\infty} e^{-st} dt = 0 $$
于是
$$ \lim_{s\to\infty} [sF(s) - f(0_+)] = 0 $$
即
$$ f(0_+) = \lim_{t\to 0_+} f(t) = \lim_{s\to\infty} sF(s) $$

#### (2) 终值定理

**证明**: 根据拉氏变换的微分性质，有
$$ L\left[\frac{d}{dt}f(t)\right] = \int_{0}^{\infty} \frac{d}{dt}f(t)e^{-st} dt = sF(s) - f(0) $$
令 $s \rightarrow 0$ 取极限得
$$ \lim_{s\to 0} \int_{0}^{\infty} \frac{d}{dt}f(t)e^{-st} dt = \lim_{s\to 0} [sF(s) - f(0)] $$
$$ \lim_{s\to 0} \int_{0}^{\infty} \frac{d}{dt}f(t)e^{-st} dt = \int_{0}^{\infty} \frac{d}{dt}f(t) \lim_{s\to 0} e^{-st} dt $$
$$ = \int_{0}^{\infty} df(t) = \lim_{t\to\infty} \int_0^t df(t) = \lim_{t\to\infty} [f(t) - f(0)] $$
于是
$$ \lim_{t\to\infty} f(t) = \lim_{s\to 0} sF(s) $$
终值定理表明：时间函数 $f(t)$ 的稳态值与复频域中 $s=0$ 附近的 $sF(s)$ 的值相同。因此，$f(t)$ 在 $t \rightarrow \infty$ 时的值可以直接从 $\lim_{s\to 0} sF(s)$ 得到。
利用该性质，可在复频域中得到控制系统在时间域中的稳态值，利用该性质还可以求得控制系统的稳态误差。

**例**: 已知 $f(t) = e^{-t} \cos t \cdot u(t)$，求 $f(0_+)$ 和 $f(\infty)$。
**解**: 由于 $L[\cos t \cdot u(t) e^{-t}] = \frac{s+1}{(s+1)^2 + 1}$。
由初值定理，得
$$ f(0_+) = \lim_{s\to\infty} sF(s) = \lim_{s\to\infty} \frac{s(s+1)}{(s+1)^2 + 1} = 1 $$
由终值定理，得
$$ f(\infty) = \lim_{s\to 0} sF(s) = \lim_{s\to 0} \frac{s(s+1)}{(s+1)^2 + 1} = 0 $$

## 2.3 拉普拉斯逆变换

由象函数 $F(s)$ 求原函数 $f(t)$，可根据公式
$$ f(t) = \frac{1}{2\pi j}\int_{\beta - j\infty}^{\beta + j\infty} F(s)e^{st} ds $$
$t \ge 0$
$$ f(t) = L^{-1}[F(s)] $$
对于简单的象函数，可直接应用拉氏变换对照表，查出相应的原函数。
对于有理分式这类复杂象函数，通常先用 **部分分式展开法**（也称海维赛德展开定理），将复杂函数展开成简单函数的和，再应用拉氏变换对照表，即可写出相应的原函数。

**试求** $F(s) = \frac{s}{s^2 + 2s + 5}$ **的拉氏逆变换**。
**解**
$$ L[F(s)] = L^{-1}\left[\frac{s}{s^2 + 2s + 5}\right] $$
$$ = L^{-1}\left[\frac{(s+1)-1}{(s+1)^2 + 2^2}\right] $$
$$ = L^{-1}\left[\frac{s+1}{(s+1)^2 + 2^2}\right] - L^{-1}\left[\frac{1}{2} \cdot \frac{2}{(s+1)^2 + 2^2}\right] $$
$$ = (e^{-t} \cos 2t - \frac{1}{2}e^{-t} \sin 2t) \cdot u(t) $$

一般系统，通常有如下形式的有理分式：
$$ F(s) = \frac{N(s)}{D(s)} = \frac{b_m s^m + b_{m-1}s^{m-1} + \dots + b_1 s + b_0}{a_n s^n + a_{n-1}s^{n-1} + \dots + a_1 s + a_0} $$
$a_1, a_2, \dots, a_n, b_1, b_2, \dots, b_m$ 都是实常数；$m, n$ 为正整数，通常 $m < n$。
$$ F(s) = \frac{N(s)}{D(s)} = K \frac{(s-z_1)(s-z_2)\dots(s-z_m)}{(s-p_1)(s-p_2)\dots(s-p_n)} $$
当分母 $D(s)=0$ 时 $s$ 的根，称为 $F(s)$ 的 **极点**；
当分子 $N(s)=0$ 时 $s$ 的根，称为 $F(s)$ 的 **零点**。

### 1. 只含不同极点的情况

$$ F(s) = \frac{N(s)}{D(s)} = \frac{a_1}{(s-p_1)} + \dots + \frac{a_k}{(s-p_k)} + \dots + \frac{a_n}{(s-p_n)} $$
$a_1, \dots, a_k, \dots, a_n$ 都是常数，叫做极点 $s=p_k$ 上的 **留数**。
$$ a_k = [(s-p_k)F(s)]_{s=p_k} $$
拉氏逆变换，求得：
$$ f(t) = L^{-1}[F(s)] = (a_1 e^{p_1 t} + \dots + a_k e^{p_k t} + \dots + a_n e^{p_n t}) \cdot u(t) $$

**试求** $F(s) = \frac{s+3}{s^2 + 3s + 2}$ **的拉氏逆变换**。
**解** $F(s)$ 的部分分式展开为
$$ F(s) = \frac{s+3}{(s+1)(s+2)} = \frac{a_1}{(s+1)} + \frac{a_2}{(s+2)} $$
又
$$ a_1 = \left[(s+1)\frac{s+3}{(s+1)(s+2)}\right]_{s=-1} = \left[\frac{s+3}{s+2}\right]_{s=-1} = \frac{2}{1} = 2 $$
$$ a_2 = \left[(s+2)\frac{s+3}{(s+1)(s+2)}\right]_{s=-2} = \left[\frac{s+3}{s+1}\right]_{s=-2} = \frac{1}{-1} = -1 $$
则 $f(t) = L^{-1}[F(s)] = L^{-1}\left[\frac{2}{s+1}\right] + L^{-1}\left[\frac{-1}{s+2}\right] = (2e^{-t} - e^{-2t})u(t)$。

**试求** $G(s) = \frac{s^3 + 5s^2 + 9s + 7}{s^2 + 3s + 2}$ **的拉氏逆变换**。
**解**: 因为分子多项式的阶次比分母多项式的阶次高，所以必须用分母去除分子。于是
$$ G(s) = s + 2 + \frac{s+3}{(s+1)(s+2)} $$
$$ g(t) = L^{-1}[G(s)] = \left(\frac{d}{dt}\delta(t) + 2\delta(t) + 2e^{-t} - e^{-2t}\right)u(t) $$

### 2. 含共轭复极点的情况

$$ F(s) = \frac{N(s)}{D(s)} = \frac{a_1 s + a_2}{s^2 + cs + d} = \frac{K_1}{s + \delta - j\beta} + \frac{K_2}{s + \delta + j\beta} $$
共轭复极点 $p_{1,2} = -\delta \pm j\beta$
分别求系数 $K_1, K_2$。 $K_1 = A + jB$， $K_2 = A - jB$。不难看出，$K_1$ 与 $K_2$ 成共轭关系。
$$ K_1 = [(s+\delta-j\beta)F(s)]_{s=-\delta+j\beta} = \frac{(a_2 - a_1\delta) + ja_1\beta}{2j\beta} = \frac{a_1\beta + (a_1\delta - a_2)j}{2\beta} $$
$$ K_2 = [(s+\delta+j\beta)F(s)]_{s=-\delta-j\beta} = \frac{(a_2 - a_1\delta) - ja_1\beta}{-2j\beta} = \frac{a_1\beta - (a_1\delta - a_2)j}{2\beta} $$
求出 $A, B$: $A = \frac{a_1}{2}$，$B = \frac{a_1\delta - a_2}{2\beta}$。
$$ L^{-1}[F(s)] = f(t) = 2e^{-\delta t} [A \cos(\beta t) - B \sin(\beta t)] \cdot u(t) $$

$$ F(s) = \frac{N(s)}{D(s)} = \frac{a_1 s + a_2}{s^2 + cs + d} + \frac{a_3}{s-p_3} + \dots + \frac{a_n}{s-p_n} $$
$$ = \frac{K_1}{s+\delta-j\beta} + \frac{K_2}{s+\delta+j\beta} + \frac{a_3}{s-p_3} + \dots + \frac{a_n}{s-p_n} $$
$$ f(t) = L^{-1}[F(s)] $$
$$ = 2e^{-\delta t} [A \cos(\beta t) - B \sin(\beta t)] \cdot u(t) + (a_3 e^{p_3 t} + a_4 e^{p_4 t} + \dots + a_n e^{p_n t}) \cdot u(t) $$
$$ = 2e^{-\delta t} [A \cos(\beta t) - B \sin(\beta t)] u(t) + \sum_{i=3}^{n} a_i e^{p_i t} u(t) $$

**例**: 试求 $F(s) = \frac{s^2 + 3}{(s^2 + 2s + 5)(s+2)}$ 的拉氏逆变换。
**解**: 将 $F(s)$ 分解因式为
$$ F(s) = \frac{s^2 + 3}{(s+1+j2)(s+1-j2)(s+2)} = \frac{K_1}{s+1-j2} + \frac{K_2}{s+1+j2} + \frac{a_3}{s+2} $$
则共轭复极点 $p_{1,2} = -1 \pm j2$， $\delta = 1$, $\beta = 2$。
分别求系数 $K_1, K_2$ 和 $a_3$:
$$ a_3 = [(s+2)F(s)]_{s=-2} = \left[(s+2)\frac{s^2+3}{(s^2+2s+5)(s+2)}\right]_{s=-2} = \left[\frac{s^2+3}{s^2+2s+5}\right]_{s=-2} = \frac{4+3}{4-4+5} = \frac{7}{5} $$
$$ K_1 = [(s+1-j2)F(s)]_{s=-1+j2} = \left[\frac{s^2+3}{(s+1+j2)(s+2)}\right]_{s=-1+j2} = \frac{(-1+j2)^2+3}{(-1+j2+1+j2)(-1+j2+2)} $$
$$ = \frac{1-4-2j+3}{(j4)(1+j2)} = \frac{-2j}{4j(1+j2)} = \frac{-1}{2(1+j2)} = \frac{-1}{2+4j} = \frac{-1(2-4j)}{4+16} = \frac{-2+4j}{20} = \frac{-1+2j}{10} $$
$K_2 = K_1^* = \frac{-1-2j}{10}$。
$A = \text{Re}(K_1) = -\frac{1}{10}$。
$B = \text{Im}(K_1) = \frac{2}{10} = \frac{1}{5}$。
$$ L^{-1}[F(s)] = f(t) = 2e^{-t} [A \cos(\beta t) - B \sin(\beta t)] + a_3 e^{-2t} $$
$$ = 2e^{-t} \left[-\frac{1}{10} \cos(2t) - \frac{1}{5} \sin(2t)\right] + \frac{7}{5} e^{-2t} $$
$$ = -\frac{1}{5}e^{-t} \cos(2t) - \frac{2}{5}e^{-t} \sin(2t) + \frac{7}{5} e^{-2t} $$

### 3. 含多重极点的情况

设 $A(s)=0$ 有 $r$ 个重极点 $p_1$ (即 $A(s)=0$ 有 $r$ 个重根 $s_1 = -p_1$)，则 $F(s)$ 可写成
$$ F(s) = \frac{B(s)}{A(s)} = \frac{B(s)}{(s+p_1)^r (s+p_{r+1})\cdots(s+p_n)} $$
$$ = \frac{a_r}{(s+p_1)^r} + \frac{a_{r-1}}{(s+p_1)^{r-1}} + \dots + \frac{a_1}{s+p_1} + \frac{a_{r+1}}{s+p_{r+1}} + \dots + \frac{a_n}{s+p_n} $$
$p_1$ 为 $F(s)$ 的重极点，$p_{r+1}, \dots, p_n$ 为 $F(s)$ 的 $(n-r)$ 个非重极点；$a_r, a_{r-1}, \dots, a_1$ 和 $a_{r+1}, a_{r+2}, \dots, a_n$ 为待定系数，按不同极点留数的方法求解。

系数 $a_r, a_{r-1}, \dots, a_1$ 可通过如下方法来计算：
令 $F_1(s) = \frac{a_r}{(s+p_1)^r} + \frac{a_{r-1}}{(s+p_1)^{r-1}} + \dots + \frac{a_1}{s+p_1}$。
用 $(s+p_1)^r$ 乘以上式的两边，得
$(s+p_1)^r F_1(s) = a_r + a_{r-1}(s+p_1) + \dots + a_1(s+p_1)^{r-1}$。
令 $s = -p_1$， $a_r = [(s+p_1)^r F_1(s)]_{s=-p_1}$。
将上式两边对 $s$ 进行微分，得
$\frac{d}{ds}[(s+p_1)^r F_1(s)] = a_{r-1} + 2a_{r-2}(s+p_1) + \dots + (r-1)a_1(s+p_1)^{r-2}$。
令 $s = -p_1$， $a_{r-1} = \frac{d}{ds}[(s+p_1)^r F_1(s)]_{s=-p_1}$。
再将上式两边对 $s$ 进行二阶微分，得
$\frac{d^2}{ds^2}[(s+p_1)^r F_1(s)] = 2a_{r-2} + 3 \cdot 2 a_{r-3}(s+p_1) + \dots + (r-1)(r-2)a_1(s+p_1)^{r-3}$。
再令 $s = -p_1$，得
$a_{r-2} = \frac{1}{2!} \frac{d^2}{ds^2}[(s+p_1)^r F_1(s)]_{s=-p_1}$。
类似地，对上式两边进行 $k$ 阶微分，并令 $s = -p_1$，
$a_{r-k} = \frac{1}{k!} \frac{d^{(k)}}{ds^{(k)}}[(s+p_1)^r F_1(s)]_{s=-p_1}$。

因此，得到递推公式：
$$ a_r = [(s+p_1)^r F_1(s)]_{s=-p_1} $$
$$ a_{r-1} = \frac{d}{ds}[(s+p_1)^r F_1(s)]_{s=-p_1} $$
$$ a_{r-2} = \frac{1}{2!} \frac{d^2}{ds^2}[(s+p_1)^r F_1(s)]_{s=-p_1} $$
$$ \dots $$
$$ a_{r-k} = \frac{1}{k!} \frac{d^{(k)}}{ds^{(k)}}[(s+p_1)^r F_1(s)]_{s=-p_1} $$
$$ \dots $$
$$ a_1 = \frac{1}{(r-1)!} \frac{d^{(r-1)}}{ds^{(r-1)}}[(s+p_1)^r F_1(s)]_{s=-p_1} $$

原函数 $f(t)$ 为
$$ f(t) = L^{-1}[F(s)] = L^{-1}[F_1(s)] + L^{-1}\left[\frac{a_{r+1}}{s+p_{r+1}} + \dots + \frac{a_n}{s+p_n}\right] $$
$$ = \left[\frac{a_r}{(r-1)!}t^{r-1} + \frac{a_{r-1}}{(r-2)!}t^{r-2} + \dots + a_2 t + a_1\right]e^{-p_1 t} \cdot u(t) + \sum_{i=r+1}^{n} a_i e^{-p_i t} \cdot u(t) $$

**【例 2.13】** 试求 $F(s) = \frac{s+2}{s(s+1)^2(s+3)}$ 的拉氏逆变换。
**解**: 当分母 $A(s)=0$， $s$ 有四个极点，即：二重极点 $p_1 = -1$ 和非重极点 $p_3 = 0, p_4 = -3$。
将 $F(s)$ 展开成部分分式
$$ F(s) = \frac{s+2}{s(s+1)^2(s+3)} = \frac{a_2}{(s+1)^2} + \frac{a_1}{s+1} + \frac{a_3}{s} + \frac{a_4}{s+3} $$
$$ a_1 = \frac{d}{ds}\left[(s+1)^2 \frac{s+2}{s(s+1)^2(s+3)}\right]_{s=-1} = \left[\frac{d}{ds}\frac{s+2}{s(s+3)}\right]_{s=-1} $$
$$ = \left[\frac{1 \cdot s(s+3) - (s+2)(2s+3)}{s^2(s+3)^2}\right]_{s=-1} = \left[\frac{s^2+3s - (2s^2+9s+6)}{s^2(s+3)^2}\right]_{s=-1} $$
$$ = \left[\frac{-s^2-6s-6}{s^2(s+3)^2}\right]_{s=-1} = \frac{-1+6-6}{1(2)^2} = -\frac{1}{4} $$
$$ a_2 = \left[(s+1)^2 \frac{s+2}{s(s+1)^2(s+3)}\right]_{s=-1} = \left[\frac{s+2}{s(s+3)}\right]_{s=-1} = \frac{-1+2}{(-1)(-1+3)} = \frac{1}{(-1)(2)} = -\frac{1}{2} $$
$$ a_3 = \left[s \frac{s+2}{s(s+1)^2(s+3)}\right]_{s=0} = \left[\frac{s+2}{(s+1)^2(s+3)}\right]_{s=0} = \frac{2}{1^2 \cdot 3} = \frac{2}{3} $$
$$ a_4 = \left[(s+3)\frac{s+2}{s(s+1)^2(s+3)}\right]_{s=-3} = \left[\frac{s+2}{s(s+1)^2}\right]_{s=-3} = \frac{-3+2}{(-3)(-3+1)^2} = \frac{-1}{(-3)(4)} = \frac{1}{12} $$
得到 $F(s)$ 的原函数 $f(t)$ 为
$$ f(t) = L^{-1}[F(s)] = -\frac{1}{2} e^{-t} t + (-\frac{1}{4}) e^{-t} + \frac{2}{3} + \frac{1}{12} e^{-3t}, \quad t \ge 0 $$

## 补充：卷积与卷积定理

### 1. 概念

由傅氏变换的卷积性质得知：两个函数 $f_1(t)$ 和 $f_2(t)$ 的卷积是指
$$ f_1(t) * f_2(t) = \int_{-\infty}^{\infty} f_1(\tau) f_2(t-\tau) d\tau $$
当 $t < 0$ 时，$f_1(t) = f_2(t) = 0$。
$$ f_1(t) * f_2(t) = \int_{-\infty}^{\infty} f_1(\tau) f_2(t-\tau) d\tau + \int_{0}^{t} f_1(\tau) f_2(t-\tau) d\tau + \int_{t}^{\infty} f_1(\tau) f_2(t-\tau) d\tau $$
$$ = \int_{0}^{t} f_1(\tau) f_2(t-\tau) d\tau $$

卷积定义：
而 $f_1(t) * f_2(t) = \int_{0}^{t} f_1(\tau) f_2(t-\tau) d\tau = \int_{0}^{t} f_2(\tau) f_1(t-\tau) d\tau = f_2(t) * f_1(t)$。
同时，卷积还满足结合律与对加法的分配律，即
$f_1(t) * [f_2(t) * f_3(t)] = [f_1(t) * f_2(t)] * f_3(t)$
$f_1(t) * [f_2(t) + f_3(t)] = f_1(t) * f_2(t) + f_1(t) * f_3(t)$

### 2. 卷积定理

若 $L[f_1(t)] = F_1(s)$， $L[f_2(t)] = F_2(s)$，则
$$ L[f_1(t) * f_2(t)] = F_1(s)F_2(s) $$

**证明**: 由拉氏变换和卷积的定义，可以写出
$$ L[f_1(t) * f_2(t)] = \int_{0}^{\infty} [f_1(t) * f_2(t)] e^{-st} dt $$
$$ = \int_{0}^{\infty} \left[ \int_{0}^{t} f_1(\tau) f_2(t-\tau) d\tau \right] e^{-st} dt $$
当 $t > \tau$ 时，有 $f_2(t-\tau)u(t-\tau) = f_2(t-\tau)$，即
$$ f_2(t - \tau)u(t - \tau) = \begin{cases} 0 & t < \tau \\ f_2(t-\tau) & t > \tau \end{cases} $$
因此，
$$ \int_{0}^{\infty} \left[ \int_{0}^{t} f_1(\tau) f_2(t-\tau) u(t-\tau) d\tau \right] e^{-st} dt = \int_{0}^{\infty} \int_{0}^{\infty} f_1(\tau) f_2(t-\tau) u(t-\tau) e^{-st} dt d\tau $$
$$ = \int_{0}^{\infty} f_1(\tau) d\tau \int_{0}^{\infty} f_2(t-\tau) u(t-\tau) e^{-st} dt $$
令 $t - \tau = \lambda$， $L[f_1(t) * f_2(t)] = \int_{0}^{\infty} f_1(\tau) d\tau \int_{0}^{\infty} f_2(\lambda) u(\lambda) e^{-s(\lambda+\tau)} d\lambda$
$$ = \int_{0}^{\infty} f_1(\tau)e^{-s\tau} d\tau \int_{0}^{\infty} f_2(\lambda)e^{-s\lambda} d\lambda $$
$$ = F_1(s)F_2(s) $$
同理可得复频域的卷积定理（也称时域相乘定理）
$$ L[f_1(t)f_2(t)] = \frac{1}{2\pi j} [F_1(s) * F_2(s)] = \frac{1}{2\pi j} \int_{\beta-j\infty}^{\beta+j\infty} F_1(p)F_2(s-p)dp $$

## 62. 拉氏变换表

1) 教材《自动控制原理》：
*   变换对照表序号：1-6, 14-15, 20, 22
*   第6版：P589-600
*   第7版：P630-633

2) 电子版材料：
*   拉氏变换的性质表2-1：P32
*   常用函数拉氏变换表2-2：P33

## 63. 本章作业 (请写题目)

《机械控制工程基础》玄兆燕编-第2版
页码：P36

**所有作业要求**:
*   图、公式用尺
*   正规作业纸
*   作业纸单面答题
*   字迹工整
*   包含推导过程

**√ 2.1 (奇数编号共4小题)**
**√ 2.2, 2.3**
**√ 2.5 (偶数编号共4小题)**

**习题 2**

2.1 试求下列函数的拉氏变换。(奇)(偶)
(1) $f(t) = (4t+5)\cdot \delta(t)+(t+2)\cdot u(t)$
(2) $f(t) = t^2 e^{at} u(t)$
(3) $f(t) = \begin{cases} \sin t & 0 < t \le \pi \\ 0 & t < 0, t > \pi \end{cases}$

## 64. 本章作业

**第2章 拉普拉斯变换 37**

**所有作业要求**:
*   图、公式用尺
*   正规作业纸
*   作业纸单面答题
*   字迹工整
*   包含推导过程

(4) $f(t) = \sin(5t+\frac{\pi}{3})u(t)$
(5) $f(t) = [4 \cos(2t)]\frac{t}{3} - u(t)$
(6) $f(t) = e^{-6t} (\cos 8t + \sin 8t)u(t)$
(7) $f(t) = (t e^{-3t} + e^{-t} \cos 2t + e^{-t} \sin 4t)u(t)$
(8) $f(t) = 5 - u(t-2) + (t-1)^2 e^{2t} u(t)$

2.2 已知 $F(s) = \frac{s}{s(s+1)}$，试求：
(1) 利用终值定理，求 $t \rightarrow \infty$ 时 $f(t)$ 的值。
(2) 通过 $F(s)$ 的拉氏逆变换，求 $t \rightarrow \infty$ 时 $f(t)$ 的值。

2.3 已知 $F(s) = \frac{1}{(s+2)^2}$，试求：
(1) 利用初值定理，求 $f(0_+)$ 和 $f'(0)$ 的值。
(2) 通过 $F(s)$ 的拉氏逆变换，求 $f(t)$ 和 $f'(t)$，然后求 $f(0_+)$ 和 $f'(0_+)$。

2.4 求如题 2.4 图所示的各种波形的拉氏变换。

2.5 试求下列函数的拉氏逆变换。(偶数)(奇)
(1) $F(s) = \frac{1}{s^2+4}$
(2) $F(s) = \frac{s}{s^2-2s+5} + \frac{s+1}{s^2+9}$
(3) $F(s) = e^{-t} u(t-1)$
(4) $F(s) = \frac{e^{-s}}{s-1}$
(5) $F(s) = \frac{4}{s^2+8s+4}$
(6) $F(s) = \frac{s}{s^2+9}$
(7) $F(s) = \frac{3}{(s+1)^2}$
(8) $F(s) = \frac{1}{(s+2)(s^2+2s+2)}$