# 第五章 线性系统的频域分析法

## 本章主要内容：

* 5.1 频率特性概述
* 5.2 频率特性的极坐标图（Nyquist图）
* 5.3 频率特性的对数坐标图（Bode图）
* 5.4 闭环的频率特性与指标
* 5.5 最小相位系统与非最小相位系统
* 5.6 根据频率特性曲线估计系统传递函数

## 时间响应 vs 频率响应

* **时间响应**：直观，但高阶系统困难，信息量少。瞬态响应 → 动态性能。
* **频率响应**：稳态输出响应 → 动态性能。将传递函数从复域引入频域来分析系统的特性。

## 5.1 频率特性

### 5.1.1 频率特性的概念

**1. 频率响应**
线性定常系统对正弦输入的稳态响应称为频率响应。

假设系统稳定，若输入为 $x_i(t) = X_i \sin \omega t$，则稳态输出信号为 $x_o(t) = X_o(\omega) \sin[\omega t + \phi(\omega)]$。

**输出信号特性：**
* 也是正弦信号，频率与输入信号相同。
* 幅值和相位发生了变化。

**证明（以传递函数 $G(s) = \frac{K}{Ts+1}$ 为例）：**
* 输入信号的拉普拉斯变换：$X_i(s) = \frac{X_i \omega}{s^2 + \omega^2}$
* 输出信号的拉普拉斯变换：$X_o(s) = G(s)X_i(s) = \frac{K}{Ts+1} \cdot \frac{X_i \omega}{s^2 + \omega^2}$
* 由于 $s = -1/T$ 是 $G(s)$ 的极点（或系统微分方程的特征根）且为负值，系统稳定，瞬态分量 $\frac{X_i KT \omega}{1+T^2 \omega^2} e^{-t/T}$ 随时间推移趋于零。
* 稳态分量为 $\frac{X_i K}{\sqrt{1+T^2 \omega^2}} \sin(\omega t - \arctan T \omega)$。

**此时系统只剩下稳态输出：**
稳态输出是一个与输入同频率的正弦信号。
* $X_o(\omega)$ - 输出幅值
* $\phi(\omega)$ - 相位差
* $X_o(\omega)$ 和 $\phi(\omega)$ 均为 $\omega$ 的非线性函数。

**2. 频率特性**
线性系统在正弦输入作用下，其稳态输出幅值和相位随频率 $\omega$ 的变化而变化，这反映了系统本身的特性。

将反映该特性的表达式 $\frac{X_o(\omega)}{X_i}$ (幅频特性 $A(\omega)$) 和 $\phi(\omega)$ (相频特性) 称为系统的频率特性。

* **幅频特性**：$A(\omega) = \frac{X_o(\omega)}{X_i} = |G(j\omega)|$
* **相频特性**：$\phi(\omega) = \angle G(j\omega)$ (例如：$-\arctan T\omega$)

**3. 频率特性与传递函数的关系**
设描述系统的微分方程为：
$a_n x_o^{(n)}(t) + a_{n-1} x_o^{(n-1)}(t) + \dots + a_1 \dot{x}_o(t) + a_o x_o(t) = b_m x_i^{(m)}(t) + b_{m-1} x_i^{(m-1)}(t) + \dots + b_1 \dot{x}_i(t) + b_o x_i(t)$
传递函数为：$G(s) = \frac{X_o(s)}{X_i(s)} = \frac{b_m s^m + b_{m-1} s^{m-1} + \dots + b_1 s + b_o}{a_n s^n + a_{n-1} s^{n-1} + \dots + a_1 s + a_o}$

当输入信号为正弦信号，即 $x_i(t) = X_i \sin \omega t$，其拉普拉斯变换为 $X_i(s) = \frac{X_i \omega}{s^2 + \omega^2}$。

若系统无重极点，则输出 $X_o(s)$ 可通过部分分式展开。对上式进行拉普拉斯逆变换可得系统输出 $x_o(t)$。

对于稳定系统，特征根 $s_i$ 均具有负实部，当 $t \to \infty$ 时，瞬态分量将衰减为零。

因此，上式只剩下稳态分量，系统稳态响应为：$x_o(t) = Be^{j\omega t} + B^*e^{-j\omega t}$。

待定系数 $B$ 和 $B^*$ 可通过 $G(s)$ 计算得到。
$B = \left. G(s) \frac{X_i \omega}{(s-j\omega)(s+j\omega)} (s-j\omega) \right|_{s=j\omega} = G(j\omega) \cdot \frac{X_i}{2j}$
$B^* = \left. G(s) \frac{X_i \omega}{(s-j\omega)(s+j\omega)} (s+j\omega) \right|_{s=-j\omega} = G(-j\omega) \cdot \frac{X_i}{-2j}$

则系统的稳态响应为：
$x_o(t) = |G(j\omega)| X_i \frac{e^{j[\omega t + \angle G(j\omega)]} - e^{-j[\omega t + \angle G(j\omega)]}}{2j} = |G(j\omega)| X_i \sin[\omega t + \angle G(j\omega)]$

根据频率特性定义可知，系统幅频特性和相频特性分别为：
$A(\omega) = \frac{X_o(\omega)}{X_i} = |G(j\omega)|$
$\phi(\omega) = \angle G(j\omega)$

故 $G(j\omega) = |G(j\omega)| e^{j\angle G(j\omega)}$ 就是系统的频率特性，它是将 $G(s)$ 中的 $s$ 用 $j\omega$ 取代后的结果，是 $\omega$ 的复变函数。

$G(j\omega)$ 是一个复变函数，写成实部和虚部之和，即 $G(j\omega) = \text{Re}[G(j\omega)] + \text{Im}[G(j\omega)] = u(\omega) + jv(\omega)$。

一个系统可以用微分方程或传递函数来描述，也可用频率特性来描述。

### 5.1.2 频率特性的特点和作用

1. **频率特性可通过频率响应试验求取**
   首先改变输入正弦信号的频率 $\omega$，并测出相应的输出幅值 $X_o(\omega)$ 与相移 $\phi(\omega)$。
   然后作出幅值比 $\frac{X_o(\omega)}{X_i}$ 对频率 $\omega$ 的函数曲线，即幅频特性曲线；作出相移 $\phi(\omega)$ 对频率 $\omega$ 的函数曲线，即相频特性曲线。

2. **频率特性是单位脉冲响应函数的频谱**
   设某系统的输出为 $X_o(s) = G(s)X_i(s)$。
   若输入为单位脉冲函数，则 $X_i(s) = 1$。
   故 $X_o(j\omega) = G(j\omega)X_i(j\omega) = G(j\omega)$。
   这表明系统的频率特性就是单位脉冲响应函数的 Fourier 变换或其频谱，所以对频率特性的分析就是对单位脉冲响应函数的频谱分析。

3. **在研究系统结构及参数的变化对系统性能的影响时，许多情况下（例如对于单输入、单输出系统），在频域中分析比在时域中分析要容易。**

4. **时间响应分析主要是通过分析线性系统过渡过程，以获得系统的动态特性，而频率特性分析则是通过分析不同的正弦输入时系统的稳态响应，来获得系统的动态特性。**

5. **若线性系统的阶次较高，求得系统的微分方程较困难时，用实验的方法获得频率特性会更方便。**

6. **若系统的输入信号中带有严重的噪声干扰，则对系统采用频率特性分析法可设计出合适的通频带，以抑制噪声的影响。**

## 5.2 频率特性的极坐标图 (Nyquist图)

实部和虚部分别表示为：
$U(\omega) = \text{Re}[G(j\omega)]$
$V(\omega) = \text{Im}[G(j\omega)]$

幅值和相角分别表示为：
$A(\omega) = \sqrt{U^2(\omega) + V^2(\omega)}$
$\phi(\omega) = \arctan \frac{V(\omega)}{U(\omega)}$

**逆时针方向旋转为正，顺时针方向旋转为负。**

当 $\omega$ 从 $0 \to \infty$ 时，$G(j\omega)$ 端点轨迹即为频率特性的极坐标图，或称 Nyquist 图。它不仅表示了幅频特性和相频特性，而且也表示了实频特性和虚频特性。图中的箭头方向为 $\omega$ 从小到大的方向。

### 5.2.1 典型环节的 Nyquist 图

**1. 比例环节**
   传递函数：$G(s) = \frac{X_o(s)}{X_i(s)} = K$
   频率特性：$G(j\omega) = K$
   实频特性 $u(\omega)$ 恒为 $K$，虚频特性 $v(\omega)$ 恒为 $0$。
   幅频特性 $|G(j\omega)| = K$，相频特性 $\angle G(j\omega) = 0^\circ$
   **比例环节频率特性的 Nyquist 图是实轴上的一个定点，其坐标为 $(K, j0)$**。

**2. 积分环节**
   传递函数：$G(s) = \frac{X_o(s)}{X_i(s)} = \frac{1}{Ts}$
   频率特性：$G(j\omega) = \frac{1}{jT\omega}$
   实频特性 $u(\omega)$ 为 $0$，虚频特性 $v(\omega) = -\frac{1}{T\omega}$。
   幅频特性 $|G(j\omega)| = \frac{1}{T\omega}$，相频特性 $\angle G(j\omega) = -90^\circ$
   **当 $\omega$ 从 $0 \to \infty$ 时，$|G(j\omega)|$ 从 $\infty \to 0$，相位总是 $-90^\circ$。积分环节频率特性的 Nyquist 图是虚轴下半轴，由无穷远点指向原点，具有恒定的相位滞后。**

**3. 微分环节**
   传递函数：$G(s) = \frac{X_o(s)}{X_i(s)} = Ts$
   频率特性：$G(j\omega) = jT\omega$
   实频特性 $u(\omega)$ 恒为 $0$，虚频特性 $v(\omega) = T\omega$。
   幅频特性 $|G(j\omega)| = T\omega$，相频特性 $\angle G(j\omega) = 90^\circ$
   **当 $\omega$ 从 $0 \to \infty$ 时，$|G(j\omega)|$ 从 $0 \to \infty$，相位总是 $90^\circ$。微分环节频率特性的 Nyquist 图是虚轴上半轴，由原点指向无穷远点，具有恒定的相位超前。**

**4. 惯性环节**
   传递函数：$G(s) = \frac{X_o(s)}{X_i(s)} = \frac{K}{Ts+1}$
   频率特性：$G(j\omega) = \frac{K}{jT\omega+1} = \frac{K}{1+T^2\omega^2} - j \frac{KT\omega}{1+T^2\omega^2}$
   实频特性：$u(\omega) = \frac{K}{1+T^2\omega^2}$
   虚频特性：$v(\omega) = -\frac{KT\omega}{1+T^2\omega^2}$
   幅频特性：$|G(j\omega)| = \frac{K}{\sqrt{1+T^2\omega^2}}$
   相频特性：$\angle G(j\omega) = -\arctan T\omega$

   **当 $\omega \to \infty$ 时，$|G(j\omega)| = 0$，$\angle G(j\omega) = -90^\circ$。**

   **当 $\omega$ 从 $0 \to \infty$ 时，惯性环节 Nyquist 图为一个半圆。**
   **惯性环节频率特性的 Nyquist 图是一个以 $\left(\frac{K}{2}, j0\right)$ 为圆心，以 $\frac{K}{2}$ 为半径的圆。**
   **惯性环节频率特性的幅值随着频率的增大而减小，因而具有低通滤波的性能。同时可以看出，它存在相位滞后，且滞后相位角随频率的增大而增大，最大相位滞后为 $90^\circ$。**

**5. 一阶微分环节**
   传递函数：$G(s) = Ts+1$
   频率特性：$G(j\omega) = Tj\omega+1 = 1 + j T\omega$
   实频特性 $u(\omega)$ 恒为 $1$，虚频特性 $v(\omega) = T\omega$。
   幅频特性 $|G(j\omega)| = \sqrt{1+T^2\omega^2}$。相频特性 $\angle G(j\omega) = \arctan T\omega$。
   **当 $\omega$ 从 $0 \to \infty$ 时，$|G(j\omega)|$ 从 $1 \to \infty$，相位由 $0^\circ \to 90^\circ$。**
   **一阶微分环节频率特性 Nyquist 图始于点 $(1, j0)$，平行于虚轴，是在第一象限的一条垂线。**

**6. 振荡环节**
   传递函数：$G(s) = \frac{\omega_n^2}{s^2 + 2\xi\omega_n s + \omega_n^2}$
   频率特性：$G(j\omega) = \frac{\omega_n^2}{-\omega^2 + j2\xi\omega_n\omega + \omega_n^2} = \frac{1}{(1-\frac{\omega^2}{\omega_n^2}) + j2\xi\frac{\omega}{\omega_n}}$
   令 $\frac{\omega}{\omega_n} = \lambda$。
   $G(j\omega) = \frac{1}{(1-\lambda^2) + j2\xi\lambda}$
   幅频特性：$|G(j\omega)| = \frac{1}{\sqrt{(1-\lambda^2)^2 + 4\xi^2\lambda^2}}$
   相频特性：$\angle G(j\omega) = \begin{cases} -\arctan \frac{2\xi\lambda}{1-\lambda^2}, & \omega \le \omega_n \text{ (第四象限)} \\ -180^\circ + \arctan \frac{2\xi\lambda}{\lambda^2-1}, & \omega > \omega_n \text{ (第三象限)} \end{cases}$
   **振荡环节频率特性的 Nyquist 图始于点 $(1, j0)$，终于点 $(0, j0)$。曲线和虚轴交点的频率是无阻尼固有频率 $\omega_r$，此时的幅值 $1/(2\xi)$。曲线在第三、四象限。**
   **当 $\omega \to \infty$ 时，$|G(j\omega)|=0$，$\angle G(j\omega)=-180^\circ$。**

   $\xi$ 取值不同，Nyquist 图的形状也不同。
   当阻尼比 $\xi < 0.707$ 时，幅频特性 $|G(j\omega)|$ 在频率为 $\omega_r$ (或频率比 $\lambda_r = \omega_r/\omega_n$) 处出现的峰值此峰值称为谐振峰值，对应的频率 $\omega_r$ 称为谐振频率。

**7. 延时环节**
   传递函数：$G(s) = e^{-s\tau}$
   频率特性：$G(j\omega) = e^{-j\tau\omega} = \cos \tau\omega - j \sin \tau\omega$
   幅频特性：$|G(j\omega)| = 1$，相频特性 $\angle G(j\omega) = -\tau\omega = -\frac{180}{\pi}\tau\omega$。
   **延时环节频率特性的 Nyquist 图是一单位圆。其幅值恒为 1，而相位 $\angle G(j\omega)$ 随 $\omega$ 顺时针方向的变化成正比变化，即端点在单位圆上无限循环。**

### 5.2.2 Nyquist 图的一般绘制方法

“描点法”— 绘制 Nyquist 概略图一般步骤如下：
(1) 由 $G(j\omega)$ 求出曲线上非特征点：
    * 实频特性 $u(\omega) = \text{Re}[G(j\omega)]$、虚频特性 $v(\omega) = \text{Im}[G(j\omega)]$
    * 幅频特性 $|G(j\omega)|$、相频特性 $\angle G(j\omega)$ 的表达式；
(2) 求出若干特征点：
    * 如起点 $(\omega=0)$、终点 $(\omega=\infty)$，
    * 与实轴的交点 $(\text{Im}[G(j\omega)]=0)$ 、
    * 与虚轴的交点 $(\text{Re}[G(j\omega)]=0)$，渐近线等。
   并在极坐标图上用“箭头”标注 $\omega$ 的变化方向。

(3) 补充必要的几点:
    * 根据 $\text{Re}[G(j\omega)]$、 $\text{Im}[G(j\omega)]$ 和 $|G(j\omega)|$、$\angle G(j\omega)$ 的变化趋势以及 $G(j\omega)$ 所处的象限，做出 Nyquist 的大致图形。

(4) $\omega$ 由 $0 \to -\infty$ 的 Nyquist 图与 $\omega$ 由 $0 \to +\infty$ 的 Nyquist 图关于实轴对称。
   （除比例环节外，各环节虚部是关于 $\omega$ 的奇函数）

## 5.3 频率特性的对数坐标图 (Bode图)

### 5.3.1 概述

频率特性的对数坐标图又称为 Bode 图。对数坐标图由对数幅频特性图和对数相频特性图组成，分别表示幅频特性和相频特性。

**1. 对数坐标图的横坐标**
   对数坐标图的横坐标表示频率 $\omega$，但按对数 ($\log \omega$) 分度，单位是弧度／秒或者 s$^{-1}$。

**2. 对数幅频、相频特性图的纵坐标**
   * 对数幅频特性图的纵坐标表示 $G(j\omega)$ 的幅值，用对数 $20 \lg |G(j\omega)|$ 表示，单位是分贝，记 dB，在坐标轴上采用线性刻度分度。
   * 对数相频特性纵坐标为相角 $\angle G$ 或者 $\phi(\omega)$，单位为度 ($^\circ$)。

### 5.3.2 典型环节的 Bode 图

**1. 比例环节**
   频率特性为 $G(j\omega) = K$
   对数幅频特性为 $20\lg|G(j\omega)| = 20\lg K$
   对数相频特性为 $\angle G(j\omega) = 0^\circ$
   * 对数幅频特性曲线：随 $\omega$ 的增加，在对数坐标中是一条高度为 $20\lg K$ 的水平直线；
   * 对数相频特性曲线：是与 $0^\circ$ 重合的一直线。

**2. 积分环节**
   频率特性 $G(j\omega) = \frac{1}{j\omega}$
   幅频特性 $|G(j\omega)| = \frac{1}{\omega}$，相频特性 $\angle G(j\omega) = -90^\circ$
   对数幅频特性为 $20\lg|G(j\omega)| = 20\lg \frac{1}{\omega} = -20\lg \omega$
   * **积分环节的对数幅频特性曲线：在整个频率范围内是一条斜率为 -20dB/dec 的直线。**
   * 当 $\omega = 1$ 时，$20\lg|G(j\omega)| = 0$，即在此频率时，积分环节的对数幅频特性曲线与 0 dB 线相交。
   * **积分环节的对数相频特性曲线在整个频率范围内为一条 $-90^\circ$ 的水平线。**

**3. 微分环节**
   频率特性 $G(j\omega) = j\omega$
   幅频特性 $|G(j\omega)| = \omega$，相频特性 $\angle G(j\omega) = 90^\circ$
   对数幅频特性为 $20\lg|G(j\omega)| = 20\lg \omega$
   **每当频率增加 10 倍时，对数幅频特性就增加 20dB。**
   * **微分环节的对数幅频特性曲线：在整个频率范围内是一条斜率为 20dB/dec 的直线。**
   * 当 $\omega = 1$ 时，$20\lg|G(j\omega)| = 0$，即在此频率时，微分环节的对数幅频特性曲线与 0dB 线相交。
   * **微分环节的对数相频特性曲线在整个频率范围内为一条 $90^\circ$ 的水平线。**

**4. 惯性环节**
   频率特性 $G(j\omega) = \frac{1}{Tj\omega+1}$
   若令 $\omega_T = \frac{1}{T}$，则 $G(j\omega) = \frac{1}{1 + j\frac{\omega}{\omega_T}}$
   幅频特性为 $|G(j\omega)| = \frac{1}{\sqrt{1 + (\frac{\omega}{\omega_T})^2}}$
   相频特性为 $\angle G(j\omega) = -\arctan \frac{\omega}{\omega_T}$
   对数幅频特性为 $20\lg|G(j\omega)| = 20\lg\omega_T - 20\lg\sqrt{\omega_T^2 + \omega^2}$
   * 当 $\omega \ll \omega_T$ 时，$20\lg|G(j\omega)| \approx 20\lg\omega_T - 20\lg\omega_T = 0 \text{ dB}$
   * 当 $\omega \gg \omega_T$ 时，$20\lg|G(j\omega)| \approx 20\lg\omega_T - 20\lg\omega$
   * 当 $\omega = \omega_T$ 时，$20\lg|G(j\omega_T)| = 0 \text{ dB}$

   * **对数幅频特性在 $\omega \ll \omega_T$ 的低频段近似为 0dB 水平线，称为低频渐近线。**
   * **在 $\omega \gg \omega_T$ 的高频段近似为一条斜率 -20dB/dec 直线，称为高频渐近线。**
   * $\omega_T$ 是低频渐近线和高频渐近线交点处的频率，称为转角频率，$\omega_T = 1/T$。
   * 由惯性环节的相频特性 $\angle G(j\omega) = -\arctan \frac{\omega}{\omega_T}$，有：
     * 当 $\omega = 0$ 时，$\angle G(j\omega) = 0^\circ$
     * 当 $\omega = \omega_T$ 时，$\angle G(j\omega) = -45^\circ$
     * 当 $\omega = \infty$ 时，$\angle G(j\omega) = -90^\circ$

**5. 一阶微分环节**
   频率特性 $G(j\omega) = T j\omega + 1 = \frac{T\omega_T + j\omega}{\omega_T}$ ($\omega_T = \frac{1}{T}$)
   对数幅频特性为 $20\lg|G(j\omega)| = 20\lg\sqrt{\omega_T^2 + \omega^2} - 20\lg\omega_T$
   相频特性 $\angle G(j\omega) = \arctan \frac{\omega}{\omega_T}$
   * **与惯性环节的对数幅频特性和相频特性比较，仅相差一个符号。**
   * **一阶微分环节的对数频率特性与惯性环节的对数频率特性呈镜像关系对称于横轴。**

**6. 振荡环节**
   幅频特性：$|G(j\omega)| = \frac{1}{\sqrt{(1-\lambda^2)^2 + 4\xi^2\lambda^2}}$
   相频特性：$\angle G(j\omega) = \begin{cases} -\arctan \frac{2\xi\lambda}{1-\lambda^2}, & \omega \le \omega_n \text{ (第四象限)} \\ -180^\circ + \arctan \frac{2\xi\lambda}{\lambda^2-1}, & \omega > \omega_n \text{ (第三象限)} \end{cases}$
   * **振荡环节的转角频率 $\omega_T = \omega_n$。**

**7. 延时环节**
   频率特性 $G(j\omega) = \cos \tau\omega - j \sin \tau\omega$
   对数幅频特性为 $20\lg|G(j\omega)| = 0 \text{ dB}$
   * **对数幅频特性为 0 dB 线。相频特性随 $\omega$ 增加而线性增加，在线性坐标中，应是一条直线；但对数相频特性是一曲线。**

### 5.3.3 Bode 图的一般绘制方法

绘制系统 Bode 图的一般步骤如下：
(1) 将系统传递函数转化为若干个标准形式的环节的传递函数；
(2) 由传递函数求出频率特性；
(3) 确定各典型环节的转角频率 (升序排列)；
(4) 做出各环节的对数幅频特性的渐近线；
(5) 根据误差修正曲线对渐近线进行修正 (不用)；
(6) 将各环节的对数幅频特性叠加 (不包括系统的增益 K)；
(7) 将叠加后的曲线垂直移动 $20\lg K$，得到系统的对数幅频特性；
(8) 作各环节的对数相频特性，叠加得到系统总的对数相频特性；
(9) 有延时环节时，对数幅频特性不变，对数相频特性则加上 $-\tau\omega$。

## 5.4 闭环系统的频率特性与指标

**1. 闭环系统的频率特性**
   单位反馈系统的闭环频率特性 $G_B(j\omega)$ 与开环频率特性 $G_K(j\omega)$ 的关系为：
   $G_B(j\omega) = \frac{X_o(j\omega)}{X_i(j\omega)} = \frac{G_K(j\omega)}{1+G_K(j\omega)}$

   $G_B(j\omega)$ 的幅值和相位可分别写为：
   $G_B(j\omega) = \frac{|G_K(j\omega)|}{|1+G_K(j\omega)|} = A_B(\omega)$
   $\angle G_B(j\omega) = \angle G_K(j\omega) - \angle [1+G_K(j\omega)] = \phi_B(j\omega)$
   因此，可以由系统的开环频率特性得到系统的闭环频率特性。

**2. 闭环系统的频域性能指标**
   * 零频值 $A(0)$
   * 复现频率 $\omega_M$ 与复现带宽 $0 \sim \omega_M$
   * 谐振频率 $\omega_r$ 与相对谐振峰值 $M_r$
   * 截止频率 $\omega_b$ 与截止带宽 $0 \sim \omega_b$
   通常，系统带宽越大，响应的快速性越好。

## 5.5 最小相位系统与非最小相位系统

* 在复平面 [s] 右半平面上没有极点和零点的开环传递函数称为最小相位传递函数；反之，在 [s] 右半平面上有极点和(或)零点的开环传递函数称为非最小相位传递函数。
* 具有最小相位传递函数的系统称为最小相位系统；反之，具有非最小相位传递函数的系统称为非最小相位系统。

**非最小相位系统：**
1. 比例环节 $G(s) = -K, (K > 0)$
2. 惯性环节 $G(s) = \frac{1}{-Ts+1}, (T>0)$
3. 一阶微分环节 $G(s) = -Ts+1, (T>0)$
4. 振荡环节 $G(s) = \frac{1}{T^2 s^2 - 2\xi Ts + 1}, (T>0, 0 < \xi < 1)$

* 最小相位系统，幅值特性和相角特性之间具有唯一的对应关系。
* 这意味着，如果系统的幅值曲线在从零到无穷大的全部频率范围上给定，则相角曲线被唯一确定。反之亦然。
* **但是，这个结论对于非最小相位系统不成立！**

## 5.6 根据频率特性曲线估计系统传递函数

### 5.6.1 确定放大倍数 K

**1. 0型系统**
   $v=0$，其低频段特性为 0 dB/dec 的直线，直线的
   高度为 $20 \lg K$，由此可以求得放大倍数 $K$。

**2. Ⅰ型系统**
   $v=1$，其低频段特性为 –20dB/dec。
   如果系统各转角频率均大于 $\omega=1$，Ⅰ型系统幅频
   特性 Bode 图在 $\omega=1$ 处的高度为 $20\lg K$。

   如果系统有转角频率小于 $\omega=1$，则首段 -20dB/dec 斜
   率线的延长线与 $\omega=1$ 线的交点高度为 $20\lg K$。

### 5.6.2 各环节传递函数确定

Bode 图在其转角频率处的变化量等于相应环节
的变化量，其变化原则为：
* 惯性环节： –20dB/dec
* 一阶微分环节： 20dB/dec
* 振荡环节： –40dB/dec