import sympy as sp
from sympy import exp,cos,I
x,y,z,a,b,c = sp.symbols('xyzabc')
alist = [a*b*exp(x*I - y),(a-c)*exp(x**2+y*I)]

n = sp.Symbol('n',integer=True)
sp.sum(alist[n],(n,0,len(alist)))