import sympy as sp

theta1, theta2, theta3 = sp.symbols('theta1 theta2 theta3')
l1, l2, l3, l4 = sp.symbols('l1 l2 l3 l4')

def Rz(theta):
    return sp.Matrix([
        [sp.cos(theta), -sp.sin(theta), 0, 0],
        [sp.sin(theta), sp.cos(theta), 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1]])

def Ry(theta):
    return sp.Matrix([
        [sp.cos(theta), 0, sp.sin(theta), 0],
        [0, 1, 0, 0],
        [-sp.sin(theta), 0, sp.cos(theta), 0],
        [0, 0, 0, 1]])

def Rx(theta):
    return sp.Matrix([
        [1, 0, 0, 0],
        [0, sp.cos(theta), -sp.sin(theta), 0],
        [0, sp.sin(theta), sp.cos(theta), 0],
        [0, 0, 0, 1]])

def T(x,y,z):
    return sp.Matrix([
        [1,0,0,x],
        [0,1,0,y],
        [0,0,1,z],
        [0,0,0,1]])

T_total = Rz(theta1)*T(-l2,0,l1)*Rx(theta2)*T(0,l3,0)*Rx(theta3)*T(0,l4,0)
position = T_total[:3,3]
print("position")
sp.pprint(position)
