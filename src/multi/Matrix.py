import sympy as sp

def M(m,L):
  # 定义符号变量
  a1, r11, r12, x1, da1, dr11, dr12, m1 = sp.symbols('a1 r11 r12 x1 da1 dr11 dr12 m1')

  # 定义位置矢量和速度矢量
  R1 = sp.Matrix([r11, r12])
  dR1 = sp.Matrix([dr11, dr12])
  q1 = sp.Matrix([R1, a1])
  dq1 = sp.Matrix([dR1, da1])

  # 定义动能
  T1 = sp.Rational(1, 2) * m1 * dR1.dot(dR1)

  # 定义旋转矩阵 A1
  A1 = sp.Matrix([
      [sp.cos(a1), -sp.sin(a1)],
      [sp.sin(a1), sp.cos(a1)]
  ])

  # 计算旋转矩阵的导数
  dA1 = A1.diff(a1)

  # 定义位置偏移
  U_1 = sp.Matrix([x1, 0])

  # 计算新的位置矢量
  r1p = R1 + A1 * U_1

  # 计算新的速度矢量
  dr1p = dR1 + dA1 * U_1 * da1

  # 定义 G1 矩阵
  G1 = sp.Matrix.hstack(sp.eye(2), dA1 * U_1)

  # 打印结果
  print("动能 T1:")
  sp.pprint(T1)

  print("\n旋转矩阵 A1:")
  sp.pprint(A1)

  print("\n旋转矩阵的导数 dA1:")
  sp.pprint(dA1)

  print("\n新的位置矢量 r1p:")
  sp.pprint(r1p)

  print("\n新的速度矢量 dr1p:")
  sp.pprint(dr1p)

  print("\nG1 矩阵:")
  sp.pprint(G1)

  # 计算积分并简化
  M1_1 = sp.simplify(m1 / L * sp.integrate(G1.T * G1, (x1, -0.5 * L, 0.5 * L)))
  M1_2 = sp.simplify(m1 /L * sp.integrate(G1.T * G1, (x1, -0.2 * L, 0.8 * L)))

  # 打印结果
  print("M1{1}:")
  sp.pprint(M1_1)

  print("\nM1{2}:")
  sp.pprint(M1_2)

