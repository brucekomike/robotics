# 四连杆机构运动学分析 Python 代码详解

## 简介

该 Python 脚本使用数值方法对一个平面四连杆机构进行运动学分析。给定输入杆（连杆1，AB）的初始角度和恒定角速度，脚本计算随动杆（连杆2 BC 和连杆3 CD）的角度、角速度和角加速度，以及各连杆质心的位置、速度和加速度。最后，它将结果可视化为时间序列图和机构运动动画。

分析的核心思想是利用**约束方程**来描述机构的几何关系，并通过求解这些方程及其时间导数来获得运动学量。

## 核心概念

1.  **广义坐标 (Generalized Coordinates):** 系统的构型由一组独立的坐标描述。对于这个四连杆机构，虽然有三个连杆角度 $ \theta_1, \theta_2, \theta_3 $，但它们不是独立的。通常选择输入杆的角度 $ \theta_1 $ 作为独立坐标（驱动坐标），而 $ \theta_2 $ 和 $ \theta_3 $ 作为相关坐标。整个系统的广义坐标向量可以表示为 $ \mathbf{q} = [\theta_1, \theta_2, \theta_3]^T $。

2.  **约束方程 (Constraint Equations):** 连杆的连接方式形成了闭环约束。这些约束可以用一组非线性代数方程表示，通常写作 $ \Phi(\mathbf{q}) = \mathbf{0} $ 。对于四连杆机构，矢量环路闭合方程为：

    $$ \vec{AB} + \vec{BC} = \vec{AD} + \vec{DC} $$

    将其投影到 x 和 y 轴，得到两个标量约束方程：
    
    $$
    \begin{cases}
    \Phi_1(\mathbf{q}) = L_1 \cos(\theta_1) + L_2 \cos(\theta_2) - L_3 \cos(\theta_3) - L_4 = 0 \\
    \Phi_2(\mathbf{q}) = L_1 \sin(\theta_1) + L_2 \sin(\theta_2) - L_3 \sin(\theta_3) = 0
    \end{cases}
    $$
    
    其中 $ L_1, L_2, L_3, L_4 $ 分别是连杆 AB, BC, CD, AD 的长度。

3.  **速度约束 (Velocity Constraints):** 对位置约束方程 $ \Phi(\mathbf{q}) = \mathbf{0} $ 关于时间求导，得到速度约束方程：
4.  
    $$ \frac{d\Phi}{dt} = \frac{\partial \Phi}{\partial \mathbf{q}} \frac{d\mathbf{q}}{dt} = \mathbf{\Phi_q} \dot{\mathbf{q}} = \mathbf{0} $$

    其中 $ \mathbf{\Phi_q} $ 是约束方程的**雅可比矩阵 (Jacobian Matrix)**，$ \dot{\mathbf{q}} = [\omega_1, \omega_2, \omega_3]^T $ 是广义速度（角速度）向量。

5.  **加速度约束 (Acceleration Constraints):** 对速度约束方程 $ \mathbf{\Phi_q} \dot{\mathbf{q}} = \mathbf{0} $ 再次关于时间求导，得到加速度约束方程：
6.  
    $$ \frac{d}{dt}(\mathbf{\Phi_q} \dot{\mathbf{q}}) = \dot{\mathbf{\Phi}_q} \dot{\mathbf{q}} + \mathbf{\Phi_q} \ddot{\mathbf{q}} = \mathbf{0} $$

    整理后得到：

    $$ \mathbf{\Phi_q} \ddot{\mathbf{q}} = - \dot{\mathbf{\Phi}_q} \dot{\mathbf{q}} = \gamma $$

    其中 $ \ddot{\mathbf{q}} = [\alpha_1, \alpha_2, \alpha_3]^T $ 是广义加速度（角加速度）向量，$ \gamma = - \dot{\mathbf{\Phi}_q} \dot{\mathbf{q}} $ 是**二次速度项 (Quadratic Velocity Term)**，它包含了科里奥利加速度和向心加速度相关的项。

7.  **奇异性 (Singularity):** 当机构运动到特定构型时（例如，两连杆共线），雅可比矩阵 $ \mathbf{\Phi_q} $ 的子矩阵（用于求解相关速度/加速度）可能变得奇异（行列式为零或接近零）。在这种构型下，机构的自由度会发生瞬时变化，速度和加速度的求解会失败或产生无穷大/不确定的结果。代码中通过检查相关子矩阵的行列式来检测奇异性。

## 代码结构详解

1.  **导入库 (Imports):**
    *   `numpy`: 用于高效的数值计算和数组操作。
    *   `matplotlib.pyplot`: 用于绘制图形。
    *   `scipy.optimize.root`: 用于求解非线性方程组（位置约束）。
    *   `matplotlib.animation.FuncAnimation`: 用于创建动画。

2.  **参数定义 (Parameter Definition):**
    *   `L1`, `L2`, `L3`, `L4`: 定义四个连杆的长度。
    *   `r_C*_local`: 定义每个运动连杆质心相对于其起点（A, B, D）的局部位置（假设在连杆中点）。

3.  **输入运动 (Input Motion):**
    *   `omega1`: 输入杆 AB 的恒定角速度。
    *   `theta1_initial`: 输入杆 AB 的初始角度（添加了一个小偏移量以避免初始奇异性）。
    *   `T_total`, `dt`: 总仿真时间和时间步长。
    *   `time_steps`, `t`: 计算时间步数和时间点数组。

4.  **结果存储 (Result Storage):**
    *   使用 `np.zeros` 初始化用于存储每个时间步的角度、角速度、角加速度以及各质心位置、速度、加速度的 NumPy 数组。预分配内存可以提高效率。

5.  **约束方程函数 (`constraint_equations`):**
    *   输入：待求解的角度 `vars` ($ [\theta_2, \theta_3] $) 和已知的输入角度 `theta1_val`。
    *   输出：包含两个位置约束方程 $ \Phi_1, \Phi_2 $ 计算结果的列表。`scipy.optimize.root` 会尝试找到使这两个方程等于零的 `vars`。

6.  **雅可比矩阵函数 (`jacobian`):**
    *   输入：当前所有连杆的角度 $ \theta_1, \theta_2, \theta_3 $。
    *   输出：约束方程 $ \Phi $ 对广义坐标 $ \mathbf{q} $ 的雅可比矩阵 $ \mathbf{\Phi_q} $。

        $$ \mathbf{\Phi_q} = \begin{bmatrix} \frac{\partial \Phi_1}{\partial \theta_1} & \frac{\partial \Phi_1}{\partial \theta_2} & \frac{\partial \Phi_1}{\partial \theta_3} \\ \frac{\partial \Phi_2}{\partial \theta_1} & \frac{\partial \Phi_2}{\partial \theta_2} & \frac{\partial \Phi_2}{\partial \theta_3} \end{bmatrix} = \begin{bmatrix} -L_1 \sin\theta_1 & -L_2 \sin\theta_2 & L_3 \sin\theta_3 \\ L_1 \cos\theta_1 & L_2 \cos\theta_2 & -L_3 \cos\theta_3 \end{bmatrix} $$

7.  **Gamma 向量函数 (`gamma_vector`):**
    *   输入：当前所有连杆的角度 $ \theta_1, \theta_2, \theta_3 $ 和角速度 $ \omega_1, \omega_2, \omega_3 $。
    *   输出：二次速度项 $ \gamma = - \dot{\mathbf{\Phi}_q} \dot{\mathbf{q}} $。代码中直接使用了 $ \gamma $ 的展开形式：
  
        $$ \gamma = \begin{bmatrix} L_1 \cos(\theta_1) \omega_1^2 + L_2 \cos(\theta_2) \omega_2^2 - L_3 \cos(\theta_3) \omega_3^2 \\ L_1 \sin(\theta_1) \omega_1^2 + L_2 \sin(\theta_2) \omega_2^2 - L_3 \sin(\theta_3) \omega_3^2 \end{bmatrix} $$

8.  **初始猜测 (Initial Guess):**
    *   提供一个 $ [\theta_2, \theta_3] $ 的初始猜测值。
    *   调用 `root` 函数求解初始时刻 $ t=0 $ 对应的 $ \theta_2, \theta_3 $，以获得更精确的起始构型。
    *   如果求解失败，会打印警告并检查机构是否能在该 $ \theta_1 $ 下组装。

9.  **数值仿真循环 (Numerical Simulation Loop):**
    *   遍历所有时间步 `i`。
    *   **步骤 1: 计算当前输入运动:**
        *   根据 $ \theta_1(t) = \theta_{1,initial} + \omega_1 t $ 计算当前 $ \theta_1 $。
        *   $ \dot{\theta}_1 = \omega_1 $ (恒定)。
        *   $ \ddot{\theta}_1 = \alpha_1 = 0 $ (因为 $ \omega_1 $ 恒定)。
        *   存储 $ \theta_1, \omega_1, \alpha_1 $。
    *   **步骤 2: 求解位置约束:**
        *   使用上一时间步的 $ [\theta_2, \theta_3] $ 作为当前步求解的初始猜测值。
        *   调用 `root(constraint_equations, ...)` 求解当前 $ \theta_1 $ 对应的 $ \theta_2, \theta_3 $。
        *   进行错误处理：如果求解失败，打印错误信息，检查是否可组装，设置 `simulation_successful = False`，用 `NaN` 填充剩余结果并退出循环。
        *   使用 `np.arctan2` 对求解出的 $ \theta_2, \theta_3 $ 进行归一化（保持在 $ [-\pi, \pi] $ 区间），以避免角度值无限增大。
        *   存储 $ \theta_2, \theta_3 $。
    *   **步骤 3: 求解速度约束:**
        *   从雅可比矩阵 $ \mathbf{\Phi_q} $ 和已知的 $ \omega_1 $ 建立线性方程组 $ \mathbf{A}_{vel} [\omega_2, \omega_3]^T = \mathbf{b}_{vel} $，其中 $ \mathbf{A}_{vel} = \mathbf{\Phi_q}_{(:, 1:3)} $ (雅可比矩阵的后两列)，$ \mathbf{b}_{vel} = - \mathbf{\Phi_q}_{(:, 0)} \omega_1 $ (雅可比矩阵的第一列乘以 $ -\omega_1 $)。
        *   **奇异性检查:** 计算 $ \mathbf{A}_{vel} $ 的行列式 `det_A_vel`。如果其绝对值小于阈值 `singularity_threshold`，则认为处于奇异或近奇异状态。
        *   如果检测到奇异性，打印警告，并将 $ \omega_2, \omega_3 $ 设为 `NaN`。
        *   如果非奇异，则使用 `np.linalg.solve` 求解 $ \omega_2, \omega_3 $。如果求解过程中发生 `LinAlgError` (例如矩阵病态)，则捕获异常，打印警告，并将 $ \omega_2, \omega_3 $ 设为 `NaN`。
        *   存储 $ \omega_2, \omega_3 $。
    *   **步骤 4: 求解加速度约束:**
        *   检查 $ \omega_2, \omega_3 $ 是否为 `NaN` 或者是否处于奇异状态。如果是，则无法计算 $ \gamma $ 或求解加速度，将 $ \alpha_2, \alpha_3 $ 设为 `NaN`。
        *   如果速度有效且非奇异：
            *   调用 `gamma_vector` 计算 $ \gamma $。
            *   建立线性方程组 $ \mathbf{A}_{acc} [\alpha_2, \alpha_3]^T = \mathbf{b}_{acc} $，其中 $ \mathbf{A}_{acc} = \mathbf{A}_{vel} $，$ \mathbf{b}_{acc} = \gamma - \mathbf{\Phi_q}_{(:, 0)} \alpha_1 $。
            *   使用 `np.linalg.solve` 求解 $ \alpha_2, \alpha_3 $。同样进行 `LinAlgError` 异常处理。
        *   存储 $ \alpha_2, \alpha_3 $。
    *   **步骤 5: 计算质心运动学:**
        *   检查当前步所需的所有角度和速度是否有效（非 `NaN`）。如果无效，则将该步所有质心运动学结果设为 `NaN`。
        *   **杆 1 质心 (C1):** 基于 $ \theta_1, \omega_1, \alpha_1 $ 计算 C1 的位置 $ \mathbf{r}_{C1} $、速度 $ \mathbf{v}_{C1} $ 和加速度 $ \mathbf{a}_{C1} $。加速度计算前检查 $ \alpha_1 $ 是否为 `NaN`。
        *   **杆 2 质心 (C2):** C2 的运动是基点 B 的运动与 C2 相对于 B 的运动的叠加。计算 B 点的运动学量（基于 $ \theta_1, \omega_1, \alpha_1 $），然后计算 C2 相对于 B 的运动学量（基于 $ \theta_2, \omega_2, \alpha_2 $），最后相加得到 C2 的绝对运动学量 $ \mathbf{r}_{C2}, \mathbf{v}_{C2}, \mathbf{a}_{C2} $。加速度计算前检查 $ \alpha_1, \alpha_2 $ 是否为 `NaN`。
        *   **杆 3 质心 (C3):** C3 的运动是相对于固定点 D 的转动。基于 $ \theta_3, \omega_3, \alpha_3 $ 计算 C3 相对于 D 的运动学量，由于 D 点固定，这就是 C3 的绝对运动学量 $ \mathbf{r}_{C3}, \mathbf{v}_{C3}, \mathbf{a}_{C3} $。加速度计算前检查 $ \alpha_3 $ 是否为 `NaN`。
        *   存储所有质心的 $ \mathbf{r}, \mathbf{v}, \mathbf{a} $。

10. **绘图 (Plotting):**
    *   如果仿真成功（或部分成功），则生成图表。
    *   **图 1:** 绘制三个连杆的角度 $ \theta_1, \theta_2, \theta_3 $ 随时间变化的曲线。
    *   **图 2:** 绘制三个连杆的角速度 $ \omega_1, \omega_2, \omega_3 $ 随时间变化的曲线。
    *   **图 3:** 绘制三个连杆的角加速度 $ \alpha_1, \alpha_2, \alpha_3 $ 随时间变化的曲线。
    *   **后续图:** 分别为每个运动连杆（L1, L2, L3）绘制其质心位置 (x, y)、速度 (vx, vy) 和加速度 (ax, ay) 随时间变化的曲线。`NaN` 值在绘图中通常表现为断点。

11. **可视化/动画 (Visualization/Animation):**
    *   如果仿真结果中存在有效数据（非 `NaN`），则创建动画。
    *   设置画布和坐标轴，使其比例相等 (`set_aspect('equal')`)。
    *   绘制固定的机架 AD。
    *   初始化表示连杆 AB, BC, CD 的线对象 (`link1`, `link2`, `link3`)、表示关节 B, C 的点对象 (`joint_B`, `joint_C`) 以及表示质心的标记对象 (`com1`, `com2`, `com3`)。
    *   定义 `animate` 函数，该函数在每一帧被调用：
        *   获取当前帧 `i_frame` 对应的角度 $ \theta_1, \theta_2, \theta_3 $。
        *   检查角度是否为 `NaN`。如果是，则隐藏所有运动部件，并在时间文本中注明。
        *   如果角度有效，则计算关节 B 和 C 的坐标。
        *   更新线对象 (`set_data`) 来绘制连杆。
        *   更新点对象 (`set_data`) 来绘制关节。
        *   更新质心标记对象 (`set_data`)，如果质心数据为 `NaN` 则隐藏。
        *   更新显示时间的文本。
        *   返回所有更新的绘图元素列表（为了 `blit=True` 优化）。
    *   使用 `FuncAnimation` 创建动画对象 `ani`。
        *   `frames`: 指定要渲染的帧（通过 `frame_skip` 控制帧率，避免动画过慢）。
        *   `interval`: 帧之间的延迟（毫秒）。
        *   `blit=True`: 优化渲染速度（只重绘改变的部分）。
        *   `repeat=False`: 动画不循环播放。
    *   (注释掉) 提供了保存动画为 MP4 文件的代码（需要 `ffmpeg`）。

12. **显示绘图和动画 (Display Plots and Animation):**
    *   调用 `plt.show()` 显示所有创建的图形窗口。程序会在此暂停，直到所有窗口被关闭。

## 总结

该脚本通过数值求解约束方程及其导数，对四连杆机构进行了详细的运动学分析。它考虑了奇异性问题，并通过将无效计算结果设为 `NaN` 来处理这些情况，保证了程序的鲁棒性。最终通过图形和动画直观地展示了分析结果。

## 脚本内容

```{literalinclude} ../../../../src/multi/dynamics.py
:caption: dynamics.py
:language: python
```
