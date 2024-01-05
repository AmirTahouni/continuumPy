import sympy as sp
import numpy as np
from continuumPy import *


x,y,z = sp.symbols(['x','y','z'])
t1,t2,t3 = sp.symbols(['t1', 't2', 't3'])
e1,e2,e3 = sp.symbols(['e1', 'e2', 'e3'])
s0 = sp.symbols('s0')
v = sp.symbols('v')

# PROBLEM 1
g_vectors = sp.Matrix([[1, 0, 0],[0, 1, 0],[0, 0, 1]])
cartesian = coordinateSystem(g_vectors, [e1,e2,e3], initVariance = 'co')
#R = sp.Matrix([ t1*sp.cos(t3)*sp.cos(t2) , t1*sp.cos(t3)*sp.sin(t2) , t1*sp.sin(t3) ])
R = sp.Matrix([ t1*sp.cos(t2) , t1*sp.sin(t2) , t3 ])
trans = transformation(R, sp.Matrix([t1,t2,t3]))
gPrime_vectors = trans.transformBaseVectors(cartesian, 'co')
cylindrical = coordinateSystem(gPrime_vectors, [t1,t2,t3], initVariance='co')

sigma = tensor2ndOrder(cartesian,[[-1.7538,-2.3896,-4.0515],[-2.3896,-3.1478,-5.4021],[-4.0515,-5.4021,-6.7526]], 'contra-contra')
sigmaPrime = tensor2ndOrder(cylindrical,trans.transformMatrix(sigma, 'contra-contra'),'contra-contra')

sigmaTransformed1 = sp.Matrix(sigmaPrime.matrix_['co-co'])
sigmaTransformed2 = sp.Matrix(sigmaPrime.matrix_['contra-contra'])

sp.pprint(sigmaTransformed1.subs({t1:3/0.6,t2:np.arccos(0.6),t3:5}).evalf())
sp.pprint(sigmaTransformed2.subs({t1:3/0.6,t2:np.arccos(0.6),t3:5}).evalf())