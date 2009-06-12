import numpy as np
import sympy as sp
#from x.py import coeff
x= np.arange(9)
print x
x.shape = (3,3)
print x

y = np.arange(24).reshape(4,6)
print y
print y[1:4,2:4]

z = np.arange(10,1,-1);print z
print z[np.array([3,3,1,8])]
print z[np.array([1,3,5,7])]

print y[np.array([0,1,2]),np.array([1,2,2])]

W,X,Y,Z = sp.symbols('WXYZ')
c = np.array([W,X,Y,Z])
print c.reshape(4,1)
print c.reshape(1,4).transpose()

x = sp.Symbol('x')
p = sp.Wild('p')
q = sp.Wild('q')
print (5*x**2 + 3*x).match(p*x**2 + q*x)

print (x**2).match(p*x**q)

print sp.Integral(x,(x,0,1)).evalf()

def printx():
    print 'x'