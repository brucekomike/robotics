# 正运动学

## DH 法
### Step 1: Define Joint Parameters

Define the joint parameters:
- $\theta_i$: Joint angle
- $d_i$: Link offset
- $a_i$: Link length
- $\alpha_i$: Link twist

### Step 2: Construct Transformation Matrix

For each link $i$, construct the transformation matrix $T_i$ using the Denavit-Hartenberg parameters:

$$
T_i = 
\begin{bmatrix}
\cos \theta_i & -\sin \theta_i \cos \alpha_i & \sin \theta_i \sin \alpha_i & a_i \cos \theta_i \\
\sin \theta_i & \cos \theta_i \cos \alpha_i & -\cos \theta_i \sin \alpha_i & a_i \sin \theta_i \\
0 & \sin \alpha_i & \cos \alpha_i & d_i \\
0 & 0 & 0 & 1 \\
\end{bmatrix}
$$

### Step 3: Compute Overall Transformation

Multiply the individual transformation matrices to get the overall transformation matrix $T$:

$$
T = T_1 \cdot T_2 \cdot \cdots \cdot T_n
$$

### Step 4: Extract Position and Orientation

Extract the position $(x, y, z)$ and orientation from the final transformation matrix $T$:

$$
\begin{bmatrix}
x \\
y \\
z \\
\end{bmatrix}
\quad \text{and} \quad
R = 
\begin{bmatrix}
r_{11} & r_{12} & r_{13} \\
r_{21} & r_{22} & r_{23} \\
r_{31} & r_{32} & r_{33} \\
\end{bmatrix}
$$

## 运动向量法

