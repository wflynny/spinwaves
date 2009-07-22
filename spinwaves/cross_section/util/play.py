import matplotlib
matplotlib.use('WXAgg')
import pylab
from matplotlib._pylab_helpers import Gcf
import sympy as sp
from sympy import oo,I,cos,sin,exp
import numpy as np
import scipy as sci
from scipy.integrate import quad, inf
import sys
from scipy.optimize.slsqp import approx_jacobian
import periodictable
from periodictable import elements
from spinwaves.spinwavecalc.readfiles import atom
from numpy.random import uniform

if 1:
    #x=uniform(-2,2,200)
    #y=uniform(-2,2,200)
    #z=x**2+y**2
    #X,Y=np.meshgrid(x,y)
    xi=np.linspace(-2.1,2.1,100)
    yi=np.linspace(-2.1,2.1,100)
    zi=xi**2+yi**2
    Z=np.zeros((xi.shape[0],yi.shape[0]))
    #zi=matplotlib.mlab.griddata(x, y, z, xi, yi)
    #print zi.shape
    
    Z[range(len(xi)),range(len(yi))]=zi
    pylab.contourf(xi,yi,Z)
    pylab.colorbar()
    pylab.show()
    sys.exit()



if 1:
    a = np.array([1,1,1])
    b = 5
    print b/a
    print a/b

if 1:
    at1 = atom(pos=[0,0,0], atomicNum = 26, valence = 3)
    print at1.atomicNum
    el = elements[at1.atomicNum]
    print el
    Q=np.arange(0.,1.,1./100.)
    print Q
    Mq = el.magnetic_ff[at1.valence].M_Q(Q)
    print type(Mq)
    print type(Q)
    pylab.plot(Q,Mq)
    pylab.show()

if 1:
    x,y = sp.symbols('xy')
    print (x+1)*y

if 0:
    x,y,z = sp.symbols('xyz', real = False)
    a,b,c = sp.symbols('xyz')
    e = x+y
    f = a+b
    print (4*x+y-z)**(0.5) - (4*a+b-c)**(0.5) == 0
    print sp.cse((4*x+y-z)**(0.5) - (4*a+b-c)**(0.5))
    print sp.sqrt(4-4*sp.cos(x)**2).expand(mul = True, multinomial=True)

if 1:
    a,b,c = sp.symbols('abc')
    x,y,z = sp.symbols('xyz')
    e = x*y*sp.exp(a*x)*sp.exp(b*y)*sp.exp(c*z)
    print e
    print sp.powsimp(e, combine = 'exp')

if 0:
    x = sp.Symbol('x', real = True)
    kx = sp.Symbol('kx', real = True)
    ky = sp.Symbol('ky', real = True)
    kz = sp.Symbol('kz', real = True)
    J,S = sp.symbols('JS', real = True)
    a = (-8.0*J**2*S**2*sp.cos(kx) + 4*J**2*S**2 + 4.0*J**2*S**2*sp.cos(kx)**2)**(0.5)
    b = (-8.0*sp.cos(kx) + 4 + 4.0*sp.cos(kx)**2)**(0.5)
    c = (-8.0*sp.cos(kx) + 4 + 4.0*sp.cos(kx)**2)
    d = sp.sqrt(sp.cos(kx)**2+sp.sin(kx)**2)
    e = (1-sp.cos(kx))**2
    
    f = -8.0*J**2*S**2*cos(kx) + 4*J**2*S**2 + 4.0*J**2*S**2*cos(kx)**2
    g = -8.0*J**2*S**2*cos(kx) + 4*J**2*S**2 + 4.0*J**2*S**2*cos(kx)**2
    print f == g
    
#    print c
#    print sp.simplify(c)
#    print c.expand()
#    print sp.factor((c.subs(sp.cos(kx),x)))
#    print sp.simplify(sp.factor((a.args[0].subs(sp.cos(kx),x))).subs(x,sp.cos(kx))**a.args[1])

if 0:
    x = np.linspace(0.,1.,num=10)
    y = np.linspace(0.,1.,num=10)
    z = np.linspace(0.,1.,num=10)
    print [x,y,z]
    k = np.array([[x[i],y[i],z[i]] for i in range(len(x))])
    print k
    k = np.array([x,y,z])
    print k
    print np.indices((4,4))
    print k.shape[0] - 3

if 0:
    d = 3
    f = d != 0 
    print 'f',f,isinstance(f,bool)
    a = [1,2,3,4,5]
    c = [0,1,2,0,3]
    b = {'1':1,'2':2,'3':3}
    print isinstance(a,list)
    print isinstance(b,dict)
    print
    
    print (np.nonzero((np.array([True,False])!=0.) | (np.array([False,False])!=0.)))[0]
    print
    
    x = np.array([[True,False],[False,False]])
    y = (x != 0)
    print 'x',x
    print 'y',y
    print x[x>0]
    print x[y]
    
    y = np.array([False, False])
    print y
    z = np.array([True])
    print y[z]
    print np.where(y != True)[0]


if 0:
    f = sys._getframe(0)
    file_to_run = f.f_locals.get('__file__', None)
    assert file_to_run is not None
    print file_to_run

if 0:
    x,y,z = sp.symbols('xyz')
    e = x*y+x*z
    a = sp.Wild('a'); b = sp.Wild('b')
    print e.match(a*x)
    
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

if 0:
    print 1e-6
    
    
    def permutations(l):
        sz = len(l)
        if sz <= 1:
            return [l]
        return [p[:i]+[l[0]]+p[i:] for i in xrange(sz) for p in permutations(l[1:])]
    def other_permutations(l):
        all = permutations(l)
        all.remove(l)
        return all
    
    def perms(li):
        sz = len(li)
        if sz <= 1:
            return li
        perm_list = []
        for i in range(sz):
            if li[i].is_commutative:
                el = li[i]
                lis = li[:]; del lis[i]
                for j in range(len(lis)):
                    lis.insert(j,el)
                    print lis
                    #if lis not in perm_list: 
                    perm_list.append(lis)
                    del lis[j]
    #            lis.append(el)
    #            if lis not in perm_list: perm_list.append(lis)
            else: pass
        return perm_list
    
    
    a,b,c,d = sp.symbols('abcd')
    l = [a,b,c]
    print perms(l)
    print ''
    for p in permutations([a,b,c]):
        print p
    print
    for p in other_permutations([a,a,b]):
        print p
    
    x,y,z = sp.symbols('xyz',commutative=False)
    
    print ''
    ex = (a*b*c).subs(a*b,d)
    print ex
    print''
    
    
    x = sp.Symbol("X",commutative=True,real=True)
    y = sp.Symbol("Y",commutative=False,real=True)
    
    
    print x*y
    
    zlist=[]
    
    for i in range(10):
        zi = sp.Symbol("z%i"%(i,),commutative=False,real=True)
        zlist.append(zi)
        
    print zlist
    
    print x**(1/2)
    
    print exp(x)
    n=sp.Symbol("n",commutative=True,real=True)
    #print sp.sum(6*n**2 + 2**n, (n, 0, 6))
    
    
    #print sp.solve(x**2+x+1,x)
    
    #Convolution
    tau = sp.Symbol('tau')
    t,T = sp.symbols('tT')
    #print sp.integrate(tau**2*sp.DiracDelta(tau-(t-T)),(tau,-oo,oo))


print "End"