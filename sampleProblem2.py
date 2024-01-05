import sympy as sp
from continuumPy import *


x,y,z = sp.symbols(['x','y','z'])
t1,t2,t3 = sp.symbols(['t1', 't2', 't3'])
e1,e2,e3 = sp.symbols(['e1', 'e2', 'e3'])
s0 = sp.symbols('s0')
v = sp.symbols('v')

# PROBLEM 2
g_vectors = sp.Matrix([[1, 0, 0],[0, 1, 0],[0, 0, 1]])
cartesian = coordinateSystem(g_vectors, [e1,e2,e3], initVariance = 'co')
R = sp.Matrix([ t1*sp.cos(t2) , t1*sp.sin(t2) , t3 ])
trans = transformation(R, sp.Matrix([t1,t2,t3]))
cylindrical = coordinateSystem(trans.transformBaseVectors(cartesian, 'co'), [t1,t2,t3], initVariance='co')

v1 = vector(cartesian, [v,0,0], 'contra')
v1Prime = vector(cylindrical, trans.transformVector(v1, 'contra'), 'contra')
v2 = vector(cylindrical, [v,0,0], 'contra')

v1Prime_coDerivative = sp.trigsimp(sp.Matrix(v1Prime.coDerivative('co')))
v1pd = [[0,0,0],[0,0,0],[0,0,0]]
for i in range(3):
    for j in range(3):
        v1pd[i][j] = sp.diff(v1Prime.vec('co')[i],v1Prime.CS.thetaVariables_[j,0])
v1Prime_d = sp.trigsimp(sp.simplify(sp.Matrix(v1pd)))

v2_coDerivative = sp.trigsimp(sp.Matrix(v2.coDerivative('co')))
v2d = [[0,0,0],[0,0,0],[0,0,0]]
for i in range(3):
    for j in range(3):
        v2d[i][j] = sp.diff(v2.vec('co')[i],v2.CS.thetaVariables_[j,0])
v2_d = sp.trigsimp(sp.simplify(sp.Matrix(v2d)))

sp.pprint(v1Prime_coDerivative)
sp.pprint(v1Prime_d)
sp.pprint(v2_coDerivative)
sp.pprint(v2_d)