import sympy as sp

def generalized_force(m, L, F, G1):
    # 定义符号变量
    a1, r11, r12, x1, da1, dr11, dr12, m1, g = sp.symbols('a1 r11 r12 x1 da1 dr11 dr12 m1 g')
    
    # 定义位置矢量和速度矢量
    R1 = sp.Matrix([r11, r12])
    dR1 = sp.Matrix([dr11, dr12])
    q1 = sp.Matrix([R1, a1])
    dq1 = sp.Matrix([dR1, da1])
    
    # 定义旋转矩阵 A1
    A1 = sp.Matrix([
        [sp.cos(a1), -sp.sin(a1)],
        [sp.sin(a1), sp.cos(a1)]
    ])
    
    # 计算旋转矩阵的导数
    dA1 = A1.diff(a1)
    
    # 定义位置偏移
    u_p = sp.Matrix([1, 0])
    
    # 计算新的位置矢量
    r1p = R1 + A1 * u_p
    
    # 计算新的速度矢量
    dr1p = dR1 + dA1 * u_p * da1
    
    # 计算广义力
    dr1p_dq = r1p.jacobian(q1)
    dr1p_ddq = dr1p.jacobian(dq1)
    
    Q = dr1p_dq.T * sp.Matrix(F) + dr1p_ddq.T * sp.Matrix(G1)
    
    # 打印结果
    print("广义力 Q:")
    sp.pprint(Q)

# 使用示例
F = [sp.symbols('f1'), sp.symbols('f2')]
G1 = sp.Matrix([0, -sp.symbols('m') * sp.symbols('g')])

generalized_force(m=1, L=1, F=F, G1=G1)