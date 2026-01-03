# 拉格朗日动力学方程（含约束）

图片顶部的第一个公式是包含约束力的拉格朗日方程的一种形式：

$$
\frac{d}{dt} \left( \frac{\partial T}{\partial \dot{q}} \right)^T - \left( \frac{\partial T}{\partial q} \right)^T + C_q^T \lambda + C \dot{q} + K q = Q_e
$$

**符号解释:**

*   $ T $: 系统的总**动能** (Kinetic Energy)。通常是广义坐标 $ q $ 和广义速度 $ \dot{q} $ 的函数，$ T = T(q, \dot{q}) $。
    *   LaTeX: `T`
*   $ q $: 系统的**广义坐标** (Generalized Coordinates) 列向量。它包含了描述系统位形所需的最少数量的独立或非独立的坐标。
    *   LaTeX: `q`
*   $ \dot{q} $: 系统的**广义速度** (Generalized Velocities) 列向量，即 $ \dot{q} = \frac{dq}{dt} $。
    *   LaTeX: `\dot{q}`
*   $ \frac{\partial T}{\partial \dot{q}} $: 动能 $ T $ 对广义速度 $ \dot{q} $ 的**偏导数**。这是一个与广义动量相关的向量。
    *   LaTeX: `\frac{\partial T}{\partial \dot{q}}`
*   $ \frac{\partial T}{\partial q} $: 动能 $ T $ 对广义坐标 $ q $ 的**偏导数**。当动能不仅依赖于速度，还显式地依赖于坐标时（例如在旋转坐标系或变质量系统中），此项不为零。
    *   LaTeX: `\frac{\partial T}{\partial q}`
*   $ \frac{d}{dt} (\cdot) $: 对时间 $ t $ 的**全导数**。
    *   LaTeX: `\frac{d}{dt}`
*   $ (\cdot)^T $: 向量或矩阵的**转置** (Transpose)。*注意：在此图中，转置符号 $ ^T $ 应用于偏导数项的方式可能与某些标准教材略有不同，这可能是由于符号约定或推导过程的差异。通常，拉格朗日方程最终会得到一个列向量形式的方程。一种可能的解释是 $ \frac{\partial T}{\partial \dot{q}} $ 和 $ \frac{\partial T}{\partial q} $ 在这里被视为行向量形式的梯度，然后转置为列向量。*
    *   LaTeX: `^T`
*   $ C(q, t) = 0 $: 系统的**约束方程** (Constraint Equations)。这是一个向量方程，描述了广义坐标之间必须满足的几何或运动学关系。这些通常是完整约束（Holonomic Constraints）。见图右侧的 $ C(q, t) = 0_{n_c \times 1} $。
    *   LaTeX: `C(q, t) = 0`
*   $ C_q $: 约束方程关于广义坐标 $ q $ 的**雅可比矩阵** (Jacobian Matrix)，即 $ C_q = \frac{\partial C}{\partial q} $。
    *   LaTeX: `C_q`
*   $ \lambda $: **拉格朗日乘子** (Lagrange Multipliers) 列向量。它与约束相关，$ C_q^T \lambda $ 代表了维持约束所需的广义约束力/力矩。
    *   LaTeX: `\lambda`
*   $ C $: **阻尼矩阵** (Damping Matrix)。$ C \dot{q} $ 代表系统中的粘性阻尼力或耗散力。假设是线性的。
    *   LaTeX: `C`
*   $ K $: **刚度矩阵** (Stiffness Matrix)。$ K q $ 代表系统中的弹性恢复力（类似弹簧力），假设是线性的，并且源于某个势能 $ V = \frac{1}{2} q^T K q $。
    *   LaTeX: `K`
*   $ Q_e $: **广义外力** (Generalized External Forces) 列向量。它包含了所有未被 $ K q $ (保守内力) 和 $ C \dot{q} $ (阻尼力) 以及 $ C_q^T \lambda $ (约束力) 描述的力/力矩，例如驱动力、重力（如果未包含在 $T$ 或 $K$ 中）等。
    *   LaTeX: `Q_e`

---

## 方程的展开与重组

接下来的几组公式是对第一个方程进行展开和重新整理：

**第二组公式:**

$$
\left\{
\begin{aligned}
M \ddot{q} + \dot{M} \dot{q} - \left( \frac{\partial T}{\partial q} \right)^T + C_q^T \lambda + C \dot{q} + K q &= Q_e \\
Q_v = -\dot{M} \dot{q} + \left( \frac{\partial T}{\partial q} \right)^T
\end{aligned}
\right.
$$

**新增符号解释:**

*   $ M $: 系统的**质量矩阵**或**惯性矩阵** (Mass Matrix / Inertia Matrix)。它通常是广义坐标 $ q $ 的函数，$ M = M(q) $，并且是对称正定的。
    *   LaTeX: `M`
*   $ \ddot{q} $: 系统的**广义加速度** (Generalized Accelerations) 列向量，即 $ \ddot{q} = \frac{d^2q}{dt^2} $。
    *   LaTeX: `\ddot{q}`
*   $ \dot{M} $: 质量矩阵 $ M(q) $ 对时间的**全导数**。由于 $ M $ 是 $ q $ 的函数，根据链式法则，$ \dot{M} = \frac{dM}{dt} = \frac{\partial M}{\partial q} \dot{q} $。
    *   LaTeX: `\dot{M}`

**推导关系:**
对于典型的机械系统，动能可以写成 $ T = \frac{1}{2} \dot{q}^T M(q) \dot{q} $。
那么，$ \frac{\partial T}{\partial \dot{q}} = M(q) \dot{q} $。
其对时间的全导数为 $ \frac{d}{dt} \left( \frac{\partial T}{\partial \dot{q}} \right) = \frac{d}{dt} (M \dot{q}) = \dot{M} \dot{q} + M \ddot{q} $。
将此代入第一个公式（并再次假设转置符号 $^T$ 的作用是得到列向量形式），得到第二组的第一个方程。

*   $ Q_v $: 定义为**科里奥利力/离心力**相关的广义力向量 (Generalized forces due to Coriolis and centrifugal effects)。它包含了 $ \dot{M} \dot{q} $ 项和 $ \frac{\partial T}{\partial q} $ 项。注意这里的符号定义 $ Q_v = -\dot{M} \dot{q} + (\frac{\partial T}{\partial q})^T $。*（在许多文献中，科里奥利/离心力项 $ C(q, \dot{q})\dot{q} $ 定义为 $ \dot{M}\dot{q} - \frac{1}{2} \dot{q}^T \frac{\partial M}{\partial q} \dot{q} $，或者通过克氏符号定义，这里的 $ Q_v $ 定义可能与具体推导有关，特别是 $ (\frac{\partial T}{\partial q})^T $ 项的处理）。*
    *   LaTeX: `Q_v`

**第三组公式:**

$$
\left\{
\begin{aligned}
M \ddot{q} + C_q^T \lambda &= Q_e + Q_v - C \dot{q} - K q \\
M \ddot{q} + C_q^T \lambda &= Q_{ev} \\
Q_{ev} &= Q_e + Q_v - C \dot{q} - K q
\end{aligned}
\right.
$$

**新增符号解释:**

*   $ Q_{ev} $: 定义的一个**等效广义力**向量。它将所有非惯性项（$ M \ddot{q} $）和非约束力项（$ C_q^T \lambda $）合并到方程右侧。它包括了外部力 $ Q_e $、科里奥利/离心力 $ Q_v $，以及阻尼力 $ -C \dot{q} $ 和弹性力 $ -K q $ 的贡献。
    *   LaTeX: `Q_{ev}`

这组公式将动力学方程整理成 $ M \ddot{q} + C_q^T \lambda = Q_{ev} $ 的形式，显式地分开了惯性项、约束力项和其他所有力的项。这对于数值求解或控制设计通常很有用。

---

## 约束与坐标变换

图片右侧的公式涉及约束和坐标的选择：

$$
\begin{cases}
C(q, t) = 0_{n_c \times 1} \\
q = [q^{1T} \ q^{2T} \ q^{3T}]^T_{1 \dots n} \\
q_n = [q_i^T \ q_d^T]^T_{1 \dots n}
\end{cases}
$$
$$
\begin{cases}
q_d = f(q_i) \\
\dot{q}_n = B \dot{q}_i \\
\ddot{q}_n = B \ddot{q}_i + \dot{B} \dot{q}_i
\end{cases}
$$

**符号解释:**

*   $ n_c $: **约束方程的数量**。
    *   LaTeX: `n_c`
*   $ n $: **广义坐标的总数**。
    *   LaTeX: `n`
*   $ q = [q^{1T} \ q^{2T} \ q^{3T}]^T $: 将广义坐标向量 $ q $ **分块**。具体含义 $ q^1, q^2, q^3 $ 需要上下文确定，可能代表不同子系统或不同类型的坐标。
    *   LaTeX: `q = [q^{1T} \ q^{2T} \ q^{3T}]^T`
*   $ q_n $: 可能是指**完整的 $ n $ 维广义坐标向量**（与之前的 $ q $ 相同，只是换了个符号强调维度）。
    *   LaTeX: `q_n`
*   $ q_i $: **独立广义坐标** (Independent Generalized Coordinates) 向量。其数量为 $ n - n_c $（假设约束是独立的）。
    *   LaTeX: `q_i`
*   $ q_d $: **相关（或从属）广义坐标** (Dependent Generalized Coordinates) 向量。其数量为 $ n_c $。
    *   LaTeX: `q_d`
*   $ q_n = [q_i^T \ q_d^T]^T $: 将完整坐标向量 $ q_n $ 划分为**独立坐标 $ q_i $ 和相关坐标 $ q_d $**。
    *   LaTeX: `q_n = [q_i^T \ q_d^T]^T`
*   $ q_d = f(q_i) $: 对于完整约束，相关坐标 $ q_d $ 可以表示为独立坐标 $ q_i $ 的**函数**。这个函数关系 $ f $ 由约束方程 $ C(q, t) = 0 $ 隐式定义。
    *   LaTeX: `q_d = f(q_i)`
*   $ \dot{q}_n = B \dot{q}_i $: 完整坐标的速度 $ \dot{q}_n $ 与独立坐标的速度 $ \dot{q}_i $ 之间的**线性关系**。$ B $ 是一个 $ n \times (n-n_c) $ 的**速度变换矩阵**（或雅可比矩阵），$ B = \begin{bmatrix} I \\ \frac{\partial f}{\partial q_i} \end{bmatrix} $。
    *   LaTeX: `\dot{q}_n = B \dot{q}_i`
*   $ B $: **速度/雅可比变换矩阵**。
    *   LaTeX: `B`
*   $ \ddot{q}_n = B \ddot{q}_i + \dot{B} \dot{q}_i $: 完整坐标的加速度 $ \ddot{q}_n $ 与独立坐标的加速度 $ \ddot{q}_i $ 之间的关系。注意这里出现了 $ \dot{B} \dot{q}_i $ 项，表示加速度层面的变换不是简单的线性关系。
    *   LaTeX: `\ddot{q}_n = B \ddot{q}_i + \dot{B} \dot{q}_i`
*   $ \dot{B} $: 矩阵 $ B $ 对时间的**全导数**。
    *   LaTeX: `\dot{B}`

这些关系允许将系统的动力学方程从包含 $ n $ 个（可能相关的）广义坐标和 $ n_c $ 个约束方程以及拉格朗日乘子 $ \lambda $ 的形式，转换为只包含 $ n - n_c $ 个独立坐标 $ q_i $ 的、去除了 $ \lambda $ 的最小坐标形式的动力学方程。

