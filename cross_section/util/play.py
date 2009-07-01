import sympy as sp
from sympy import oo,I,cos,sin,exp
import numpy as np
from scipy.integrate import quad, inf

#alist = [1,2,3,4]
#blist = alist[:]
#blist.pop()
#print alist,blist
print quad(lambda x: sp.DiracDelta(x), -1, 1)
print quad(lambda x: sp.DiracDelta(x), -100, 100)
print quad(lambda x: sp.DiracDelta(x), -10000000, 10000000)
print quad(lambda x: sp.DiracDelta(x), -inf, inf)

x = sp.Symbol('x')
e = sp.DiracDelta(x)
print sp.integrate(e, (x,-oo,oo))
print sp.Integral(e, (x,-oo,oo))
print sp.Integral(e, (x,-oo,oo)).evalf()
print sp.Integral(e, (x,-oo,oo)).evalf(method='mpmath')
print sp.Integral(e, (x,-oo,oo)).evalf(method='scipy')
print sp.Integral(e, (x,-oo,oo)).evalf(method='numpy')
print sp.Integral(e, (x,-oo,oo)).evalf(method='')

#
#print 1e-6
#
#
#def permutations(l):
#    sz = len(l)
#    if sz <= 1:
#        return [l]
#    return [p[:i]+[l[0]]+p[i:] for i in xrange(sz) for p in permutations(l[1:])]
#def other_permutations(l):
#    all = permutations(l)
#    all.remove(l)
#    return all
#
#def perms(li):
#    sz = len(li)
#    if sz <= 1:
#        return li
#    perm_list = []
#    for i in range(sz):
#        if li[i].is_commutative:
#            el = li[i]
#            lis = li[:]; del lis[i]
#            for j in range(len(lis)):
#                lis.insert(j,el)
#                print lis
#                #if lis not in perm_list: 
#                perm_list.append(lis)
#                del lis[j]
##            lis.append(el)
##            if lis not in perm_list: perm_list.append(lis)
#        else: pass
#    return perm_list
#
#
#a,b,c,d = sp.symbols('abcd')
#l = [a,b,c]
#print perms(l)
#print ''
#for p in permutations([a,b,c]):
#    print p
#print
#for p in other_permutations([a,a,b]):
#    print p
#
#x,y,z = sp.symbols('xyz',commutative=False)
#
#print ''
#ex = (a*b*c).subs(a*b,d)
#print ex
#print''
#
#
#x = sp.Symbol("X",commutative=True,real=True)
#y = sp.Symbol("Y",commutative=False,real=True)
#
#
#print x*y
#
#zlist=[]
#
#for i in range(10):
#    zi = sp.Symbol("z%i"%(i,),commutative=False,real=True)
#    zlist.append(zi)
#    
#print zlist
#
#print x**(1/2)
#
#print exp(x)
#n=sp.Symbol("n",commutative=True,real=True)
##print sp.sum(6*n**2 + 2**n, (n, 0, 6))
#
#
##print sp.solve(x**2+x+1,x)
#
##Convolution
#tau = sp.Symbol('tau')
#t,T = sp.symbols('tT')
##print sp.integrate(tau**2*sp.DiracDelta(tau-(t-T)),(tau,-oo,oo))
