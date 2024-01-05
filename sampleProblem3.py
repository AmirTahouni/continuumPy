import sympy as sp
from continuumPy import *


x,y,z = sp.symbols(['x','y','z'])
t1,t2,t3 = sp.symbols(['t1', 't2', 't3'])
e1,e2,e3 = sp.symbols(['e1', 'e2', 'e3'])
s0 = sp.symbols('s0')
v = sp.symbols('v')

# PROBLEM 3
x1, x2, x3 = sp.symbols('x1 x2 x3')
g_vectors = sp.Matrix([[1, 0, 0],[0, 1, 0],[0, 0, 1]])
cartesian = coordinateSystem(g_vectors, [x1,x2,x3], initVariance = 'co')
R = sp.Matrix([ t1*sp.cos(t2) , t1*sp.sin(t2) , t3 ])
trans = transformation(R, sp.Matrix([t1,t2,t3]))
cylindrical = coordinateSystem(trans.transformBaseVectors(cartesian, 'co'), [t1,t2,t3], initVariance='co')

u1 = sp.Function('u1')(t1,t2,t3)
u2 = sp.Function('u2')(t1,t2,t3)
u3 = sp.Function('u3')(t1,t2,t3)
#u2/t1 and not u2 since we want physical strains and not tensorial strains
u = vector(cylindrical, [u1,u2/t1,u3], 'contra')
uCo = sp.Matrix(u.coDerivative('co'))
epsilon = (sp.Matrix(uCo)+sp.Matrix(uCo).T)/2
epsilon = sp.expand(sp.simplify(epsilon))
sp.pprint(epsilon)