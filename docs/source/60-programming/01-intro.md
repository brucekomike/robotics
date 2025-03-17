# 概论

## 数学基础
### 标记
- $O_{-xyz}$ 基坐标系
- $P_{-xyz}$ 刚体固连坐标系
- $a=(a_x,a_y,a_z)^T$ 加速度
- $a=(a_xi,a_yj,a_zk)^T$ 转换后的加速度

### 运算
#### 点积（dot/scalar product）
##### 行列式定义
${a} \times {b}$ 

$$
\mathbf{a} \times \mathbf{b} = 
\begin{vmatrix}
\mathbf{i} & \mathbf{j} & \mathbf{k} \\
a_x & a_y & a_z \\
b_x & b_y & b_z \\
\end{vmatrix}
$$(eqn:def1)

展开为

$$
\mathbf{a} \times \mathbf{b} = \left( a_y b_z - a_z b_y \right) \mathbf{i} - \left( a_x b_z - a_z b_x \right) \mathbf{j} + \left( a_x b_y - a_y b_x \right) \mathbf{k}
$$

##### 矩阵定义

$$
c = |a \times| b =
\begin{bmatrix}
0 & -a_z & a_y \\
a_z & 0 & -a_x \\
-a_y & a_x & 0 \\
\end{bmatrix}
\begin{bmatrix}
b_x\\ b_y \\ b_z
\end{bmatrix} \\
=\begin{bmatrix}
a_y b_z - a_z b_y \\
a_x b_z - a_z b_x \\ 
a_x b_y - a_y b_x
\end{bmatrix}
$$
