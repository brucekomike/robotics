# scara2
为了计算 Atom ST620 (SCARA) 机械臂的自由度数目与类型，我们将采用基于旋量理论的方法。该机械臂的构型为 RRRP，这意味着它有三个转动关节（R）和一个移动关节（P）。

**1. 机械臂构型分析与假设**

*   **RRRP 构型：**
    *   R1：第一个转动关节。
    *   R2：第二个转动关节。
    *   R3：第三个转动关节。
    *   P4：第四个移动关节。
*   **SCARA 特点：** SCARA 机械臂通常具有平行于 Z 轴的转动关节，用于在 XY 平面内进行操作，以及一个沿 Z 轴的移动关节进行高度调整。因此，我们假设所有 R 关节的旋转轴都平行于 Z 轴，P 关节的移动轴也平行于 Z 轴。
*   **坐标系设定：** 建立一个基坐标系 {0}，其原点位于第一个转动关节 (R1) 处，Z 轴向上。
*   **初始位姿假设：** 为了在非奇异构型下分析，我们假设机械臂处于一个通用非奇异位姿。

**2. 列出每个关节的运动旋量 (Twists)**

每个关节都提供一个瞬时运动（扭曲旋量）。我们将这些旋量表示在基坐标系 {0} 中。一个旋量 $\$ = \begin{pmatrix} \mathbf{\omega} \\ \mathbf{v} \end{pmatrix}$，其中 $\mathbf{\omega}$ 是角速度向量，$\mathbf{v}$ 是原点处的线速度向量。

*   **R1 关节 (Revolute Joint 1):**
    *   旋转轴：平行于 Z 轴，方向 $\hat{\mathbf{k}} = \begin{pmatrix} 0 \\ 0 \\ 1 \end{pmatrix}$。
    *   轴上一点：$\mathbf{p}_1 = \begin{pmatrix} 0 \\ 0 \\ 0 \end{pmatrix}$。
    *   运动旋量：$\$_{R1} = \begin{pmatrix} \hat{\mathbf{k}} \\ \mathbf{p}_1 \times \hat{\mathbf{k}} \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \\ 1 \\ 0 \\ 0 \\ 0 \end{pmatrix}$

*   **R2 关节 (Revolute Joint 2):**
    *   旋转轴：平行于 Z 轴，方向 $\hat{\mathbf{k}} = \begin{pmatrix} 0 \\ 0 \\ 1 \end{pmatrix}$。
    *   轴上一点：$\mathbf{p}_2 = \begin{pmatrix} L_1 \cos \theta_1 \\ L_1 \sin \theta_1 \\ 0 \end{pmatrix}$ (假设第一连杆长度为 $L_1$，关节 R1 的角度为 $\theta_1$)。
    *   运动旋量：$\$_{R2} = \begin{pmatrix} \hat{\mathbf{k}} \\ \mathbf{p}_2 \times \hat{\mathbf{k}} \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \\ 1 \\ -L_1 \sin \theta_1 \\ L_1 \cos \theta_1 \\ 0 \end{pmatrix}$

*   **R3 关节 (Revolute Joint 3):**
    *   旋转轴：平行于 Z 轴，方向 $\hat{\mathbf{k}} = \begin{pmatrix} 0 \\ 0 \\ 1 \end{pmatrix}$。
    *   轴上一点：$\mathbf{p}_3 = \begin{pmatrix} L_1 \cos \theta_1 + L_2 \cos(\theta_1+\theta_2) \\ L_1 \sin \theta_1 + L_2 \sin(\theta_1+\theta_2) \\ 0 \end{pmatrix}$ (假设第二连杆长度为 $L_2$，关节 R2 的角度为 $\theta_2$)。
    *   运动旋量：$\$_{R3} = \begin{pmatrix} \hat{\mathbf{k}} \\ \mathbf{p}_3 \times \hat{\mathbf{k}} \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \\ 1 \\ -(L_1 \sin \theta_1 + L_2 \sin(\theta_1+\theta_2)) \\ L_1 \cos \theta_1 + L_2 \cos(\theta_1+\theta_2) \\ 0 \end{pmatrix}$

*   **P4 关节 (Prismatic Joint 4):**
    *   移动轴：平行于 Z 轴，方向 $\hat{\mathbf{k}} = \begin{pmatrix} 0 \\ 0 \\ 1 \end{pmatrix}$。
    *   运动旋量：$\$_{P4} = \begin{pmatrix} \mathbf{0} \\ \hat{\mathbf{k}} \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \\ 0 \\ 0 \\ 0 \\ 1 \end{pmatrix}$

**3. 确定独立运动旋量 (自由度)**

机械臂的运动旋量空间是所有关节运动旋量的线性组合（span）。为了找到自由度，我们需要找出这些旋量的线性独立基。我们假设机械臂处于非奇异构型（即 $L_1, L_2 \neq 0$ 且 $\theta_2 \neq 0, \pi$）。

1.  **$\$_{R1}$** 提供了绕 Z 轴的转动：
    $\$_{Rz} = \begin{pmatrix} 0 \\ 0 \\ 1 \\ 0 \\ 0 \\ 0 \end{pmatrix}$ (1个自由度)

2.  考虑 $\$_{R1}$ 和 $\$_{R2}$ 的组合：
    $\$_{R2} - \$_{R1} = \begin{pmatrix} 0 \\ 0 \\ 0 \\ -L_1 \sin \theta_1 \\ L_1 \cos \theta_1 \\ 0 \end{pmatrix}$
    这是一个纯平移旋量，方向垂直于连杆 1。在非奇异构型下，这个方向与 X 轴或 Y 轴均不平行。我们可以通过旋转 $\theta_1$ 得到沿 X 或 Y 的分量。为了得到标准基向量，我们知道两个平行转动关节可以产生沿垂直于它们连线方向的平移。
    若选择 $\theta_1 = 0$，则 $\$_{R2} - \$_{R1} = \begin{pmatrix} 0 \\ 0 \\ 0 \\ 0 \\ L_1 \\ 0 \end{pmatrix}$。这代表沿 Y 轴的平移。
    若选择 $\theta_1 = \pi/2$，则 $\$_{R2} - \$_{R1} = \begin{pmatrix} 0 \\ 0 \\ 0 \\ -L_1 \\ 0 \\ 0 \end{pmatrix}$。这代表沿 X 轴的平移。
    因此，$\$_{R1}$ 和 $\$_{R2}$ 共同提供了绕 Z 轴的转动和 XY 平面内的平移能力。

3.  更一般地，从 $\$_{R1}, \$_{R2}, \$_{R3}$ 的线性组合中，我们可以提取出：
    *   绕 Z 轴的转动：$\$_{Rz} = \begin{pmatrix} 0 \\ 0 \\ 1 \\ 0 \\ 0 \\ 0 \end{pmatrix}$
    *   沿 X 轴的平移：$\$_{Tx} = \begin{pmatrix} 0 \\ 0 \\ 0 \\ 1 \\ 0 \\ 0 \end{pmatrix}$
    *   沿 Y 轴的平移：$\$_{Ty} = \begin{pmatrix} 0 \\ 0 \\ 0 \\ 0 \\ 1 \\ 0 \end{pmatrix}$
    （这三者是线性独立的，并且在非奇异构型下，SCARA 机械臂的 RRR 部分可以实现这些运动。）

4.  **$\$_{P4}$** 提供了沿 Z 轴的平移：
    $\$_{Tz} = \begin{pmatrix} 0 \\ 0 \\ 0 \\ 0 \\ 0 \\ 1 \end{pmatrix}$ (1个自由度)

综合以上，我们得到 4 个线性独立的运动旋量：$\$_{Rz}, \$_{Tx}, \$_{Ty}, \$_{Tz}$。
这些旋量构成了机械臂末端执行器运动旋量空间的基。
因此，机械臂的运动旋量空间的秩为 4。

**自由度数目：** 4
**自由度类型：**
*   绕 Z 轴的转动 ($R_z$)
*   沿 X 轴的平移 ($T_x$)
*   沿 Y 轴的平移 ($T_y$)
*   沿 Z 轴的平移 ($T_z$)

**4. 列出所有约束 (Wrenches)**

约束是与所有独立运动旋量互易（reciprocal）的力/力矩旋量（wrench）。一个力/力矩旋量 $\$\text{w} = \begin{pmatrix} \mathbf{L} \\ \mathbf{F} \end{pmatrix}$ 与一个运动旋量 $\$ = \begin{pmatrix} \mathbf{\omega} \\ \mathbf{v} \end{pmatrix}$ 互易的条件是它们的互易积为零：
$\$\text{w} \cdot \$ = \mathbf{L} \cdot \mathbf{\omega} + \mathbf{F} \cdot \mathbf{v} = 0$

我们需要找到所有满足与上述 4 个独立运动旋量互易条件的力/力矩旋量 $\$\text{w} = \begin{pmatrix} L_x \\ L_y \\ L_z \\ F_x \\ F_y \\ F_z \end{pmatrix}$。

*   **与 $\$_{Rz}$ 互易：**
    $\begin{pmatrix} L_x \\ L_y \\ L_z \\ F_x \\ F_y \\ F_z \end{pmatrix} \cdot \begin{pmatrix} 0 \\ 0 \\ 1 \\ 0 \\ 0 \\ 0 \end{pmatrix} = L_z = 0$
    （这意味着不能有绕 Z 轴的力矩）

*   **与 $\$_{Tx}$ 互易：**
    $\begin{pmatrix} L_x \\ L_y \\ L_z \\ F_x \\ F_y \\ F_z \end{pmatrix} \cdot \begin{pmatrix} 0 \\ 0 \\ 0 \\ 1 \\ 0 \\ 0 \end{pmatrix} = F_x = 0$
    （这意味着不能有沿 X 轴的力）

*   **与 $\$_{Ty}$ 互易：**
    $\begin{pmatrix} L_x \\ L_y \\ L_z \\ F_x \\ F_y \\ F_z \end{pmatrix} \cdot \begin{pmatrix} 0 \\ 0 \\ 0 \\ 0 \\ 1 \\ 0 \end{pmatrix} = F_y = 0$
    （这意味着不能有沿 Y 轴的力）

*   **与 $\$_{Tz}$ 互易：**
    $\begin{pmatrix} L_x \\ L_y \\ L_z \\ F_x \\ F_y \\ F_z \end{pmatrix} \cdot \begin{pmatrix} 0 \\ 0 \\ 0 \\ 0 \\ 0 \\ 1 \end{pmatrix} = F_z = 0$
    （这意味着不能有沿 Z 轴的力）

从上述互易条件可知，任何约束旋量 $\$\text{w}$ 必须满足 $L_z=0$，$F_x=0$，$F_y=0$，$F_z=0$。
因此，约束旋量 $\$\text{w}$ 的形式为：
$\$\text{w} = \begin{pmatrix} L_x \\ L_y \\ 0 \\ 0 \\ 0 \\ 0 \end{pmatrix}$

其中 $L_x$ 和 $L_y$ 可以是任意值。我们可以识别出两个线性独立的约束旋量：

*   **绕 X 轴的力矩 ($M_x$)：**
    $\$\text{w}_1 = \begin{pmatrix} 1 \\ 0 \\ 0 \\ 0 \\ 0 \\ 0 \end{pmatrix}$

*   **绕 Y 轴的力矩 ($M_y$)：**
    $\$\text{w}_2 = \begin{pmatrix} 0 \\ 1 \\ 0 \\ 0 \\ 0 \\ 0 \end{pmatrix}$

这两个约束旋量是线性独立的。因此，约束旋量空间的秩为 2。

**5. 自由度数目验证**

机械臂的自由度 $M$ 可以通过以下公式计算：
$M = 6 - \text{rank}(\mathbf{W})$
其中 $\text{rank}(\mathbf{W})$ 是独立约束旋量矩阵的秩。

$M = 6 - 2 = 4$

这个结果与通过识别独立运动旋量直接得到的自由度数目一致。

**总结**

*   **自由度数目：** 4
*   **自由度类型：** 绕 Z 轴的转动 ($R_z$)，沿 X 轴的平移 ($T_x$)，沿 Y 轴的平移 ($T_y$)，沿 Z 轴的平移 ($T_z$)。
*   **约束类型：** 绕 X 轴的力矩 ($M_x$)，绕 Y 轴的力矩 ($M_y$)。

这些自由度使得 SCARA 机械臂能够在三维空间中进行位置定位，并绕其 Z 轴进行姿态调整，非常适合平面装配和拾取放置任务。