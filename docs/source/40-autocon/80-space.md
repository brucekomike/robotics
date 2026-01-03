# 第八章 线性系统的状态空间表示法

## 8.1 系统数学模型的基本描述方法

### 1. 系统动态过程的两类数学描述

#### (1) 系统的外部描述 (黑箱)

外部描述通常被称为输出-输入描述。对于单输入单输出 (SISO) 线性常系数系统，时间域的外部描述可以用高阶微分方程表示：

$ y^{(n)} + a_{n-1}y^{(n-1)} + \dots + a_1y^{(1)} + a_0y = b_{m}u^{(m)} + b_{m-1}u^{(m-1)} + \dots + b_1u^{(1)} + b_0u $

频率域的等价描述是传递函数：

$ g(s) = \frac{Y(s)}{U(s)} = \frac{b_m s^m + b_{m-1}s^{m-1} + \dots + b_1s + b_0}{s^n + a_{n-1}s^{n-1} + \dots + a_1s + a_0} $

其中 $m \le n$。

#### (2) 系统的内部描述 (白箱)

内部描述是系统内部变量之间的关系。状态空间描述是内部描述的基本形式，它由两个数学方程组表征：状态方程和输出方程。

#### (3) 外部描述和内部描述的比较

*   **外部描述**：仅是对系统的一种不完全描述，不能反映黑箱内部结构中不可控或不可观测的部分。
*   **内部描述**：是对系统的一种完全描述，能够反映系统的所有动力学特性。

#### 状态和状态空间的定义

*   **状态变量组**：一个动力学系统的状态变量组定义为能够完全表征其时间域行为的最小内部变量组。
*   **状态 (向量)**：一个动力学系统的状态定义为其状态变量组 $x_1(t), x_2(t), \dots, x_n(t)$ 所组成的一个列向量：
    $ \mathbf{x}(t) = \begin{bmatrix} x_1(t) \\ x_2(t) \\ \vdots \\ x_n(t) \end{bmatrix} $
    状态 $\mathbf{x}$ 的维数等于其组成状态变量的个数，即 $\dim \mathbf{x} = n$。
*   **状态空间**：定义为状态向量的集合。状态空间的维数等于状态的维数。

#### 几点解释 (七点)

*   状态变量组对系统行为的完全表征性：给定初始时刻 $t_0$ 的状态变量组和任意时刻 $t \ge t_0$ 的输入变量组 $u(t)$，系统的任意内部变量在 $t \ge t_0$ 时刻的运动行为都可以完全确定。
*   状态变量组最小性的物理特征：减少一个变量就会破坏系统行为表征的完整性。
*   状态变量组最小性的数学特征：它是系统所有内部变量中一个极大线性无关变量组。
*   状态变量组的不唯一性：系统任意两个状态变量组之间的关系是线性非奇异变换。
*   有穷维系统和无穷维系统。
*   状态空间的属性：建立在实数域 $\mathbb{R}$ 上的一个 $n$ 维向量空间 $\mathbb{R}^n$。

### 2. 线性系统的状态空间描述

*   系统的状态空间描述（动态方程或运动方程）是对系统输入、输出和状态变量之间关系的方程组。
*   它包括：
    *   **状态方程**（描述输入和状态变量之间的关系）：
        $ \dot{\mathbf{x}} = \mathbf{Ax} + \mathbf{Bu} $
    *   **输出方程**（描述输出和输入、状态变量之间的关系）：
        $ \mathbf{y} = \mathbf{Cx} + \mathbf{Du} $
    其中 $\mathbf{x}$ 是状态向量，$\mathbf{u}$ 是输入向量，$\mathbf{y}$ 是输出向量，$\mathbf{A}$, $\mathbf{B}$, $\mathbf{C}$, $\mathbf{D}$ 是矩阵。

---

**示例 1: 电路系统状态空间描述**

对于给定的 RLC 电路，确定其状态变量和状态方程。

解：
微分方程模型：
$ L\frac{di}{dt} + Ri + u_C = u, \quad u_C = \frac{1}{C}\int i dt $
选择 $i$ 和 $u_C$ 作为状态变量：
设 $x_1 = u_C$ 和 $x_2 = i$。
则状态方程为：
$ \begin{cases} \dot{x}_1 = \frac{1}{C}x_2 \\ \dot{x}_2 = -\frac{1}{L}x_1 - \frac{R}{L}x_2 + \frac{1}{L}u \end{cases} $
矩阵形式的状态方程：
$ \begin{bmatrix} \dot{x}_1 \\ \dot{x}_2 \end{bmatrix} = \begin{bmatrix} 0 & \frac{1}{C} \\ -\frac{1}{L} & -\frac{R}{L} \end{bmatrix} \begin{bmatrix} x_1 \\ x_2 \end{bmatrix} + \begin{bmatrix} 0 \\ \frac{1}{L} \end{bmatrix} u $
输出方程为 $y = u_C$，即 $y = x_1$。
矩阵形式的输出方程：
$ y = \begin{bmatrix} 1 & 0 \end{bmatrix} \begin{bmatrix} x_1 \\ x_2 \end{bmatrix} $

---

**状态空间描述的矩阵形式：**

$ \dot{\mathbf{x}} = \mathbf{Ax} + \mathbf{Bu} $
$ \mathbf{y} = \mathbf{Cx} + \mathbf{Du} $

其中：
*   $\mathbf{A}$ 是 **系统矩阵** ($n \times n$)
*   $\mathbf{B}$ 是 **输入矩阵** ($n \times p$)
*   $\mathbf{C}$ 是 **输出矩阵** ($q \times n$)
*   $\mathbf{D}$ 是 **传输矩阵** ($q \times p$)
*   $\mathbf{x}$ 是 $n$ 维 **状态向量**
*   $\mathbf{u}$ 是 $p$ 维 **输入向量**
*   $\mathbf{y}$ 是 $q$ 维 **输出向量**

这些矩阵由系统的结构和参数决定。

**线性时不变 (LTI) 系统的状态空间描述特点：**

1.  **线性属性和时不变属性**：描述形式是线性的，系数矩阵不随时间变化。
2.  **共性属性和个性属性**：不同系统的状态空间描述具有相同的形式，差别在于参数矩阵的不同。通常用 $(A, B, C, D)$ 来代表一个 LTI 系统。
3.  **简洁性**：用向量方程形式表示，当状态、输入、输出变量增加时，表达形式的复杂性不增加。

---

## 8.2 系统数学模型间的相互转换

### 1. 由系统输入输出描述导出状态空间描述

对于单输入单输出 (SISO) 线性时不变系统，其输入输出描述（高阶微分方程或传递函数）可以导出多种形式的状态空间描述（不唯一）。

**结论 1：由输入输出描述导出状态空间描述**

给定单输入单输出线性时不变系统的输入输出描述：
$ y^{(n)} + a_{n-1}y^{(n-1)} + \dots + a_1y^{(1)} + a_0y = b_m u^{(m)} + b_{m-1}u^{(m-1)} + \dots + b_1u^{(1)} + b_0u $
或传递函数：
$ g(s) = \frac{Y(s)}{U(s)} = \frac{b_m s^m + b_{m-1}s^{m-1} + \dots + b_1s + b_0}{s^n + a_{n-1}s^{n-1} + \dots + a_1s + a_0} $

#### (1) $m=0$ 情况

此时输入输出描述简化为：
$ y^{(n)} + a_{n-1}y^{(n-1)} + \dots + a_1y^{(1)} + a_0y = b_0u $
或传递函数：
$ g(s) = \frac{b_0}{s^n + a_{n-1}s^{n-1} + \dots + a_1s + a_0} $
对应的状态空间描述为：
$ \dot{\mathbf{x}} = \begin{bmatrix} 0 & 1 & & \\ & \ddots & \ddots & \\ & & 0 & 1 \\ -a_0 & -a_1 & \dots & -a_{n-1} \end{bmatrix} \mathbf{x} + \begin{bmatrix} 0 \\ \vdots \\ 0 \\ b_0 \end{bmatrix} u $
$ y = \begin{bmatrix} 1 & 0 & \dots & 0 \end{bmatrix} \mathbf{x} $

#### (2) $m \ne 0$ 情况

此时输入输出描述为：
$ y^{(n)} + a_{n-1}y^{(n-1)} + \dots + a_1y^{(1)} + a_0y = b_m u^{(m)} + b_{m-1}u^{(m-1)} + \dots + b_1u^{(1)} + b_0u $
或传递函数：
$ g(s) = \frac{Y(s)}{U(s)} = \frac{b_m s^m + b_{m-1}s^{m-1} + \dots + b_1s + b_0}{s^n + a_{n-1}s^{n-1} + \dots + a_1s + a_0} $
其中允许 $b_n=0$，包括 $m < n$ 和 $m=n$ 两种情形。
对应的状态空间描述为：
$ \dot{\mathbf{x}} = \begin{bmatrix} 0 & 1 & & \\ & \ddots & \ddots & \\ & & 0 & 1 \\ -a_0 & -a_1 & \dots & -a_{n-1} \end{bmatrix} \mathbf{x} + \begin{bmatrix} \beta_1 \\ \beta_2 \\ \vdots \\ \beta_n \end{bmatrix} u $
$ y = \begin{bmatrix} 1 & 0 & \dots & 0 \end{bmatrix} \mathbf{x} + \beta_0 u $
其中：
$ \beta_0 = b_n $
$ \beta_1 = b_{n-1} - a_{n-1}\beta_0 $
$ \beta_2 = b_{n-2} - a_{n-1}\beta_1 - a_{n-2}\beta_0 $
$ \dots $
$ \beta_n = b_0 - a_{n-1}\beta_{n-1} - a_{n-2}\beta_{n-2} - \dots - a_0\beta_0 $
(注意：当 $m < n$ 时，$b_n=0$。当 $m=n$ 时，$b_n \ne 0$。)

---

**例 1：** 给定一个单输入单输出线性时不变系统的输入输出描述 $3\ddot{y} + 6\dot{y} + 12y = 6\dot{u} + 3u$，导出其状态空间描述。

解：
$m \ne 0$，$n=3$。状态为 3 维。
系统方程归一化（使 $y^{(n)}$ 的系数为 1）：
$ \ddot{y} + 2\dot{y} + 4y = 2\dot{u} + u $
（注意：题目原方程为 $3\ddot{y} + 6\dot{y} + 12y = 6\dot{u} + 3u$。若按此计算，则 $a_3=1, a_2=2, a_1=4, a_0=3$；$b_2=0, b_1=6, b_0=3$。但题目中给出的 $\beta$ 系数公式中，$b$ 的下标似乎是与 $u$ 的阶数相关，并且与 $n$ 的大小有关。此处，按照标准形式，将 $u$ 的最高阶导数阶数记为 $m$，则 $m=1$。）
**若按照标准格式 $y^{(n)} + \dots = b_m u^{(m)} + \dots$ 来处理，且 $m=1$, $n=3$:**
$a_3 = 1, a_2 = 2, a_1 = 4, a_0 = 3$。
$b_1 = 2, b_0 = 1$。$m=1 < n=3$。
令 $b_2=0, b_3=0$。
计算 $\beta$ 系数：
$ \beta_0 = b_3 = 0 $
$ \beta_1 = b_2 - a_2\beta_0 = 0 - 2 \times 0 = 0 $
$ \beta_2 = b_1 - a_2\beta_1 - a_1\beta_0 = 2 - 2 \times 0 - 4 \times 0 = 2 $
$ \beta_3 = b_0 - a_2\beta_2 - a_1\beta_1 - a_0\beta_0 = 1 - 2 \times 2 - 4 \times 0 - 3 \times 0 = 1 - 4 = -3 $

状态空间描述：
$ \dot{\mathbf{x}} = \begin{bmatrix} 0 & 1 & 0 \\ 0 & 0 & 1 \\ -3 & -4 & -2 \end{bmatrix} \mathbf{x} + \begin{bmatrix} 0 \\ 2 \\ -3 \end{bmatrix} u $
$ y = \begin{bmatrix} 1 & 0 & 0 \end{bmatrix} \mathbf{x} $

**注意：** 如果按照题目中的 $a_3 = 1,a_2 = 2,a_1 = 4,a_0 = 3$ 和 $b_3 = 0,b_2 = 0,b_1 = 2,b_0 = 1$ 来计算，那么 $m=1$（因为 $b_1$ 和 $b_0$ 是 $u$ 的导数项），但 $b_3$ 和 $b_2$ 却是 0。通常 $b_m$ 是 $u$ 的最高阶导数系数。在这里，可能是 $b_1 = 6$, $b_0 = 3$ （原方程）。如果原方程是 $3\ddot{y} + 6\dot{y} + 12y = 6\dot{u} + 3u$，则归一化后为 $\ddot{y} + 2\dot{y} + 4y = 2\dot{u} + u$.
那么 $a_3=1, a_2=2, a_1=4, a_0=3$, $b_1=2, b_0=1$. $m=1$.
如果按照PPT中的示例计算 $a_3 = 1,a_2 = 2,a_1 = 4,a_0 = 3$ 和 $b_3 = 0,b_2 = 0,b_1 = 2,b_0 = 1$. 并且 $\beta_0 = b_3 = 0, \beta_1 = b_2 - a_2\beta_0 = 0, \beta_2 = b_1 - a_2\beta_1 - a_1\beta_0 = 2, \beta_3 = b_0 - a_2\beta_2 - a_1\beta_1 - a_0\beta_0 = 1 - 2(2) - 4(0) - 3(0) = -3$.
状态空间描述：
$ \dot{\mathbf{x}} = \begin{bmatrix} 0 & 1 & 0 \\ 0 & 0 & 1 \\ -3 & -4 & -2 \end{bmatrix} \mathbf{x} + \begin{bmatrix} 0 \\ 2 \\ -3 \end{bmatrix} u $
$ y = \begin{bmatrix} 1 & 0 & 0 \end{bmatrix} \mathbf{x} $
**PPT 的示例计算结果:**
$ \beta_0 = b_3 = 0 $
$ \beta_1 = b_2 - a_2\beta_0 = 0 - 2 \times 0 = 0 $
$ \beta_2 = b_1 - a_2\beta_1 - a_1\beta_0 = 2 - 2 \times 0 - 4 \times 0 = 2 $
$ \beta_3 = b_0 - a_2\beta_2 - a_1\beta_1 - a_0\beta_0 = 1 - 2 \times 2 - 4 \times 0 - 3 \times 0 = -3 $
PPT 的计算结果是正确的，但是 $b_3$ 的定义在这里可能有点误导，通常 $b_m$ 是 $u^{(m)}$ 的系数。
正确的理解应该是：
$y^{(3)} + 2\ddot{y} + 4\dot{y} + 3y = 2\dot{u} + u$.
$a_3=1, a_2=2, a_1=4, a_0=3$.
$b_2=2, b_1=1$. $m=2$.
所以 $m=2$.
**重新计算 $\beta$:**
$ \beta_0 = b_3 = 0 $ (因为 $m=2$, $b_3$ 应该为 0)
$ \beta_1 = b_2 - a_2\beta_0 = 2 - 2 \times 0 = 2 $
$ \beta_2 = b_1 - a_2\beta_1 - a_1\beta_0 = 1 - 2 \times 2 - 4 \times 0 = 1 - 4 = -3 $
$ \beta_3 = b_0 - a_2\beta_2 - a_1\beta_1 - a_0\beta_0 = 0 - 2 \times (-3) - 4 \times 2 - 3 \times 0 = 6 - 8 = -2 $

**看来，PPT的例子直接使用了 $a_i$ 和 $b_i$ 作为系数，而没有明确 $m$ 和 $n$。**
**最标准的做法是：**
$y^{(n)} + a_{n-1}y^{(n-1)} + \dots + a_0y = b_m u^{(m)} + \dots + b_0 u$.
这里 $n=3$. $a_2=2, a_1=4, a_0=3$.
$m=1$. $b_1=2, b_0=1$.
$\dot{y}^{(3)} \rightarrow y^{(3)}$
$3\ddot{y} \rightarrow \ddot{y}$
$6\dot{y} \rightarrow \dot{y}$
$12y \rightarrow y$
$6\dot{u} \rightarrow \dot{u}$
$3u \rightarrow u$
所以 $a_3=1, a_2=2, a_1=4, a_0=3$.
$b_1=2, b_0=1$. $m=1$.
**状态方程和输出方程**（标准形式，$m \le n$）：
$ \dot{x}_1 = x_2 $
$ \dot{x}_2 = x_3 $
$ \dot{x}_3 = -a_0 x_1 - a_1 x_2 - a_2 x_3 + b_0 u $
$ y = b_1 x_3 + b_0 u $
（这是另一种形式，由伯德图推导而来）

**回到 PPT 中的标准形式：**
$ y^{(n)} + a_{n-1}y^{(n-1)} + \dots + a_0y = b_m u^{(m)} + \dots + b_0u $
$g(s) = \frac{\hat{y}(s)}{\hat{u}(s)} = \frac{b_m s^m + \dots + b_0}{s^n + a_{n-1}s^{n-1} + \dots + a_0}$
**在例1中：**
$3\ddot{y} + 6\dot{y} + 12y = 6\dot{u} + 3u$
$\ddot{y} + 2\dot{y} + 4y = 2\dot{u} + u$
$n=3$. $a_2=2, a_1=4, a_0=3$.
$m=1$. $b_1=2, b_0=1$.
**结论 1 (2) 的形式：**
$ \dot{\mathbf{x}} = \begin{bmatrix} 0 & 1 & & \\ & \ddots & \ddots & \\ & & 0 & 1 \\ -a_0 & -a_1 & \dots & -a_{n-1} \end{bmatrix} \mathbf{x} + \begin{bmatrix} \beta_1 \\ \beta_2 \\ \vdots \\ \beta_n \end{bmatrix} u $
$ y = \begin{bmatrix} 1 & 0 & \dots & 0 \end{bmatrix} \mathbf{x} + \beta_0 u $
这里的 $\beta$ 系数计算方法是：
$\beta_0 = b_n$.  (这里 $b_3=0$)
$\beta_1 = b_{n-1} - a_{n-1}\beta_0 = b_2 - a_2\beta_0$.  (这里 $b_2=0$)
$\beta_2 = b_{n-2} - a_{n-1}\beta_1 - a_{n-2}\beta_0 = b_1 - a_2\beta_1 - a_1\beta_0$.  (这里 $b_1=2$)
$\beta_3 = b_{n-3} - a_{n-1}\beta_2 - a_{n-2}\beta_1 - a_{n-3}\beta_0 = b_0 - a_2\beta_2 - a_1\beta_1 - a_0\beta_0$.  (这里 $b_0=1$)

**按照 PPT 的计算：**
$a_3 = 1, a_2 = 2, a_1 = 4, a_0 = 3$.
$b_3 = 0, b_2 = 0, b_1 = 2, b_0 = 1$.
$\beta_0 = b_3 = 0$
$\beta_1 = b_2 - a_2\beta_0 = 0 - 2 \times 0 = 0$
$\beta_2 = b_1 - a_2\beta_1 - a_1\beta_0 = 2 - 2 \times 0 - 4 \times 0 = 2$
$\beta_3 = b_0 - a_2\beta_2 - a_1\beta_1 - a_0\beta_0 = 1 - 2 \times 2 - 4 \times 0 - 3 \times 0 = -3$

**状态空间描述 (PPT 结果):**
$ \dot{\mathbf{x}} = \begin{bmatrix} 0 & 1 & 0 \\ 0 & 0 & 1 \\ -3 & -4 & -2 \end{bmatrix} \mathbf{x} + \begin{bmatrix} 0 \\ 2 \\ -3 \end{bmatrix} u $
$ y = \begin{bmatrix} 1 & 0 & 0 \end{bmatrix} \mathbf{x} $

### 2. 由状态空间描述导出传递函数矩阵

对状态空间描述取拉普拉斯变换，并令初始状态 $x(0)=0$:
$ s\mathbf{X}(s) = \mathbf{AX}(s) + \mathbf{BU}(s) $
$ \mathbf{Y}(s) = \mathbf{CX}(s) + \mathbf{DU}(s) $
整理得到：
$ (s\mathbf{I} - \mathbf{A})\mathbf{X}(s) = \mathbf{BU}(s) $
$ \mathbf{X}(s) = (s\mathbf{I} - \mathbf{A})^{-1}\mathbf{BU}(s) $
代入输出方程：
$ \mathbf{Y}(s) = \mathbf{C}(s\mathbf{I} - \mathbf{A})^{-1}\mathbf{BU}(s) + \mathbf{DU}(s) $
$ \mathbf{Y}(s) = [\mathbf{C}(s\mathbf{I} - \mathbf{A})^{-1}\mathbf{B} + \mathbf{D}]\mathbf{U}(s) $

**结论 2：$G(s)$ 的基本关系式**

对于多输入多输出线性时不变系统，传递函数矩阵 $G(s)$ 与系数矩阵 $\{A, B, C, D\}$ 的基本关系式为：
$ G(s) = C(sI - A)^{-1}B + D $

*   **特征矩阵**：$\Delta(sI - A)$
*   **特征多项式**：$\det(sI - A) = s^n + a_{n-1}s^{n-1} + \dots + a_1s + a_0$
    (注意：特征矩阵 $sI-A$ 必为非奇异，其逆矩阵称为预解矩阵。)

**结论 3：$G(s)$ 的实用算式**

首先计算特征多项式 $a(s) = \det(sI - A) = s^n + a_{n-1}s^{n-1} + \dots + a_1s + a_0$。
然后计算一组系数矩阵 $E_{n-1}, E_{n-2}, \dots, E_0$:
$ E_{n-1} = CB $
$ E_{n-2} = CAB + a_{n-1}CB $
$ \dots $
$ E_1 = CA^{n-2}B + a_{n-1}CA^{n-3}B + \dots + a_2CB $
$ E_0 = CA^{n-1}B + a_{n-1}CA^{n-2}B + \dots + a_1CB $

则计算 $G(s)$ 的一个实用关系式为：
$ G(s) = \frac{1}{a(s)}[E_{n-1}s^{n-1} + E_{n-2}s^{n-2} + \dots + E_1s + E_0] + D $

---

**例 1：** 给定一个线性时不变系统的状态空间描述，计算系统的传递函数矩阵 $G(s)$。
$ \dot{\mathbf{x}} = \begin{bmatrix} 2 & 0 & 0 \\ 0 & 2 & 0 \\ 0 & 3 & 1 \end{bmatrix} \mathbf{x} + \begin{bmatrix} 1 & 2 \\ 1 & 0 \\ 3 & 1 \end{bmatrix} u $
$ y = \begin{bmatrix} 1 & 1 & 2 \end{bmatrix} \mathbf{x} $

解：
(1) 计算特征多项式：
$ A = \begin{bmatrix} 2 & 0 & 0 \\ 0 & 2 & 0 \\ 0 & 3 & 1 \end{bmatrix} $
$ sI - A = \begin{bmatrix} s-2 & 0 & 0 \\ 0 & s-2 & 0 \\ 0 & -3 & s-1 \end{bmatrix} $
$ \alpha(s) = \det(sI - A) = (s-2)(s-2)(s-1) = (s^2 - 4s + 4)(s-1) = s^3 - s^2 - 4s^2 + 4s + 4s - 4 = s^3 - 5s^2 + 8s - 4 $

(2) 计算系数矩阵 $n=3$:
$ B = \begin{bmatrix} 1 & 2 \\ 1 & 0 \\ 3 & 1 \end{bmatrix} $, $ C = \begin{bmatrix} 1 & 1 & 2 \end{bmatrix} $
$ E_2 = CB = \begin{bmatrix} 1 & 1 & 2 \end{bmatrix} \begin{bmatrix} 1 & 2 \\ 1 & 0 \\ 3 & 1 \end{bmatrix} = \begin{bmatrix} 1 \times 1 + 1 \times 1 + 2 \times 3 & 1 \times 2 + 1 \times 0 + 2 \times 1 \end{bmatrix} = \begin{bmatrix} 8 & 4 \end{bmatrix} $
$ E_1 = CAB + a_2CB $
$ a_2 = -5 $ (从 $s^3 - 5s^2 + 8s - 4$ 得到)
$ CA = \begin{bmatrix} 1 & 1 & 2 \end{bmatrix} \begin{bmatrix} 2 & 0 & 0 \\ 0 & 2 & 0 \\ 0 & 3 & 1 \end{bmatrix} = \begin{bmatrix} 2 & 8 & 2 \end{bmatrix} $
$ CAB = \begin{bmatrix} 2 & 8 & 2 \end{bmatrix} \begin{bmatrix} 1 & 2 \\ 1 & 0 \\ 3 & 1 \end{bmatrix} = \begin{bmatrix} 2 \times 1 + 8 \times 1 + 2 \times 3 & 2 \times 2 + 8 \times 0 + 2 \times 1 \end{bmatrix} = \begin{bmatrix} 16 & 6 \end{bmatrix} $
$ E_1 = \begin{bmatrix} 16 & 6 \end{bmatrix} + (-5)\begin{bmatrix} 8 & 4 \end{bmatrix} = \begin{bmatrix} 16 - 40 & 6 - 20 \end{bmatrix} = \begin{bmatrix} -24 & -14 \end{bmatrix} $
$ E_0 = CA^2B + a_2CAB + a_1CB $
$ a_1 = 8 $
$ A^2 = \begin{bmatrix} 2 & 0 & 0 \\ 0 & 2 & 0 \\ 0 & 3 & 1 \end{bmatrix} \begin{bmatrix} 2 & 0 & 0 \\ 0 & 2 & 0 \\ 0 & 3 & 1 \end{bmatrix} = \begin{bmatrix} 4 & 0 & 0 \\ 0 & 4 & 0 \\ 0 & 9 & 1 \end{bmatrix} $
$ CA^2 = \begin{bmatrix} 1 & 1 & 2 \end{bmatrix} \begin{bmatrix} 4 & 0 & 0 \\ 0 & 4 & 0 \\ 0 & 9 & 1 \end{bmatrix} = \begin{bmatrix} 4 & 13 & 2 \end{bmatrix} $
$ CA^2B = \begin{bmatrix} 4 & 13 & 2 \end{bmatrix} \begin{bmatrix} 1 & 2 \\ 1 & 0 \\ 3 & 1 \end{bmatrix} = \begin{bmatrix} 4 \times 1 + 13 \times 1 + 2 \times 3 & 4 \times 2 + 13 \times 0 + 2 \times 1 \end{bmatrix} = \begin{bmatrix} 17 & 10 \end{bmatrix} $
$ E_0 = \begin{bmatrix} 17 & 10 \end{bmatrix} + (-5)\begin{bmatrix} 16 & 6 \end{bmatrix} + 8\begin{bmatrix} 8 & 4 \end{bmatrix} = \begin{bmatrix} 17 - 80 + 64 & 10 - 30 + 32 \end{bmatrix} = \begin{bmatrix} 1 & 12 \end{bmatrix} $

**(3) 计算传递函数矩阵:**
$ G(s) = \frac{1}{a(s)}[E_2s^2 + E_1s + E_0] + D $
这里 $D = 0$ (因为输入和输出之间没有直接联系)。
$ G(s) = \frac{1}{s^3 - 5s^2 + 8s - 4} \left( \begin{bmatrix} 8 & 4 \end{bmatrix} s^2 + \begin{bmatrix} -24 & -14 \end{bmatrix} s + \begin{bmatrix} 1 & 12 \end{bmatrix} \right) $
$ G(s) = \frac{1}{s^3 - 5s^2 + 8s - 4} \begin{bmatrix} 8s^2 - 24s + 1 & 4s^2 - 14s + 12 \end{bmatrix} $
**注意：** PPT 中的 $E_0$ 计算结果是 $[16 \ 12]$。
$E_0 = CA^2B + a_2CAB + a_1CB$
$a_1=8$
$CA^2B = [17 \ 10]$
$a_2CAB = -5 \times [16 \ 6] = [-80 \ -30]$
$a_1CB = 8 \times [8 \ 4] = [64 \ 32]$
$E_0 = [17 - 80 + 64, 10 - 30 + 32] = [1, 12]$
**PPT 中的 $E_0$ 计算结果是 $[16 \ 12]$。** 存在计算错误。
**PPT 的 $E_0$ 计算:**
$E_0 = CA^2 B + a_2CAB + a_1CB = [16 \ 12]$
这里 $a_1=8$
$CA^2B$ 似乎是 $[16 \ 12]$。
$a_2CAB = -5 \times [16 \ 6] = [-80 \ -30]$
$a_1CB = 8 \times [8 \ 4] = [64 \ 32]$

**根据 PPT 中的结果：**
$ G(s) = \frac{1}{s^3 - 5s^2 + 8s - 4} \begin{bmatrix} 8s^2 - 24s + 16 & 4s^2 - 14s + 12 \end{bmatrix} $
（这里 PPT 示例中的 $E_1$ 是 $[-24 \ -14]$ 且 $E_0$ 是 $[16 \ 12]$）

### Matlab 命令

*   **传递函数矩阵 $\rightarrow$ 状态空间描述:**
    `[A, B, C, D] = tf2ss(num, den)`
*   **状态空间描述 $\rightarrow$ 传递函数矩阵:**
    `[num, den] = ss2tf(A, B, C, D, iu)`
    其中 `iu` 表示输入的代号。

---

## 2.3 系统的能控性与能观性分析

### 1. 线性连续系统的能控性与能观测性

*   **直观讨论**：
    *   **能控性**：研究系统内部状态是否可由输入影响。如果系统内部每个状态变量都可由输入完全影响，则系统状态完全能控。
    *   **能观测性**：研究系统内部状态是否可由输出反映。如果系统内部每个状态变量都可由输出完全反映，则系统状态完全能观测。
*   **定义**：
    *   **能控性**：在指定初始时刻 $t_0$，存在一个控制输入 $u(t)$ 和一个时刻 $t_1 > t_0$，使得系统状态能从任意初始状态 $x(t_0)$ 转移到零状态 $x(t_1) = 0$。
    *   **能达性**：在指定初始时刻 $t_0$，存在一个控制输入 $u(t)$ 和一个时刻 $t_1 > t_0$，使得系统状态能从零状态 $x(t_0) = 0$ 转移到任意目标状态 $x(t_1) = x_f$。
        *   对于连续时间线性时不变系统，能控性和能达性等价。
    *   **不能观测性**：如果存在一个时刻 $t_1 > t_0$，使得系统以 $x(t_0)$ 为初始状态的输出 $y(t)$ 恒为零（即 $\forall t \in [t_0, t_1], y(t) = 0$ 成立），则该初始状态 $x(t_0)$ 是不能观测的。
    *   **系统能观测性**：
        *   **完全能观测**：存在有限时刻 $t_1 > t_0$，使得输出 $y(t)$ 能唯一确定状态向量的初值 $x(t_0)$。
        *   **不完全能观测**：存在一个非零状态或一个非空状态集合在时刻 $t_0$ 为不能观测。

### 2. 连续时间线性时不变系统的能控性判据

#### (1) 格拉姆矩阵判据

系统完全能控的充分必要条件是存在时刻 $t_1 > 0$，使得格拉姆矩阵 $W_c[0, t_1] = \int_0^{t_1} e^{-A\tau}BB^T e^{-A^T\tau} d\tau$ 非奇异。

*   **意义**：理论分析和推导。
*   对于完全能控的 LTI 系统，格拉姆矩阵非奇异也意味着系统完全能达。

#### (2) 秩判据

构造能控性判别矩阵：
$ Q_c = [B, AB, A^2B, \dots, A^{n-1}B] $
系统完全能控的充分必要条件是 $\text{rank}(Q_c) = n$。

---

**例 1：** 判断系统是否完全能控。
$ \dot{\mathbf{x}} = \begin{bmatrix} 4 & 0 \\ 0 & -5 \end{bmatrix} \mathbf{x} + \begin{bmatrix} 1 \\ 2 \end{bmatrix} u $
$ Q_c = [B, AB] = \begin{bmatrix} 1 & 4 \\ 2 & -10 \end{bmatrix} $
$\text{rank}(Q_c) = 2 = n$。因此，系统完全能控。

**例 2：** 图示电路，判断系统能控性。
（需要根据电路图列出状态方程，然后计算 $Q_c$ 的秩。）

**例 3：** 判断系统是否完全能控。
$ \dot{\mathbf{x}} = \begin{bmatrix} -1 & -4 & -2 \\ 0 & 6 & -1 \\ 1 & 7 & -1 \end{bmatrix} \mathbf{x} + \begin{bmatrix} 2 & 0 \\ 0 & 1 \\ 1 & 1 \end{bmatrix} u $
$ Q_c = [B, AB, A^2B] $
计算 $A^2B$ 和 $Q_c$ 的秩。
$\text{rank}(Q_c) = 3 = n$。因此，系统完全能控。

### 3. 连续时间线性时不变系统的能观测性判据

能控性和能观测性在概念、特性和判据上都是对偶的。

#### (1) 格拉姆矩阵判据

系统完全能观测的充分必要条件是存在时刻 $t_1 > 0$，使得格拉姆矩阵 $W_o[0, t_1] = \int_0^{t_1} e^{A^T\tau}C^T C e^{A\tau} d\tau$ 非奇异。

#### (2) 秩判据

构造能观测性判别矩阵：
$ Q_o = \begin{bmatrix} C \\ CA \\ \vdots \\ CA^{n-1} \end{bmatrix} $
系统完全能观测的充分必要条件是 $\text{rank}(Q_o) = n$。

---

**例 4：** 判断系统是否完全能观测。
$ \dot{\mathbf{x}} = \begin{bmatrix} -1 & -4 & -2 \\ 0 & 6 & -1 \\ 1 & 7 & -1 \end{bmatrix} \mathbf{x} $
$ y = \begin{bmatrix} 0 & 2 & 1 \\ 1 & 1 & 0 \end{bmatrix} \mathbf{x} $
计算 $C, CA, CA^2$ 并构造 $Q_o$。
$\text{rank}(Q_o) = 3 = n$。因此，系统完全能观测。

---

### Matlab 命令

*   **能控性矩阵**：`Uc = ctrb(A, B)`
*   **能观测性矩阵**：`Vo = obsv(A, C)`
*   **能控性（秩）判据**：计算 `rank(ctrb(A, B))` 是否等于 $n$。
*   **能观测性（秩）判据**：计算 `rank(obsv(A, C))` 是否等于 $n$。
*   **格拉姆矩阵**：`Wc = gram(A, B)` (能控性)，`Wo = gram(A', C')` (能观测性)。

---

## 本章小结

1.  掌握线性连续控制系统的状态空间描述法。
2.  了解传递函数与状态空间描述之间的转换方法。
3.  理解系统的能控性与能观测性。