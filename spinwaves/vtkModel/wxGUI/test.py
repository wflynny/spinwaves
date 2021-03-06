import numpy as N
import scipy
import math
from math import sqrt
import scipy.optimize
from utilities.anneal import anneal
#from scikits.openopt import NLSP
#from scikits.openopt import NLP
from openopt.Kernel.NLSP import NLSP
from openopt.Kernel import NLP
from vtkModel.BondClass import JParam

class test:
    def __init__(self):
        print 'init'
        self.val = 1
    
    def __str__(self):
        return str(self.val)



def chisq_quaternion(p,a,b,c,a2,b2,c2,x0,y0,z0,x0p,y0p,z0p):
#For the quaternion case, the equations to be solved
    nq=N.sqrt(p[0]**2+p[1]**2+p[2]**2+p[3]**2)
    xq=p[0]
    yq=p[1]
    zq=p[2]
    wq=p[3]   
    cx=p[4]
    cy=p[5]
    cz=p[6]
    rmatrix=N.matrix([[1-2*yq**2-2*zq**2,2*xq*yq-2*wq*zq,2*zq*xq+2*wq*yq],\
                      [2*xq*yq+2*wq*zq,1-2*xq**2-2*zq**2,2*yq*zq-2*wq*xq],\
                      [2*zq*xq-2*wq*yq,2*yq*zq+2*wq*xq,1-2*xq**2-2*yq**2]],'d')
    translation=N.matrix([cx,cy,cz],'d').T
    xvec0=N.matrix([x0,y0,z0],'d').T
    xvec0p=transform_point(rmatrix,translation,xvec0)
    eqn4=normal_xq(xq,yq,zq,wq,cx,cy,cz,a,b,c,a2,b2,c2,x0,y0,z0)
    eqn5=normal_yq(xq,yq,zq,wq,cx,cy,cz,a,b,c,a2,b2,c2,x0,y0,z0)
    eqn6=normal_zq(xq,yq,zq,wq,cx,cy,cz,a,b,c,a2,b2,c2,x0,y0,z0)
    eqn1=N.ravel(xvec0p[0]-x0p)[0]
    eqn2=N.ravel(xvec0p[1]-y0p)[0]
    eqn3=N.ravel(xvec0p[2]-z0p)[0]
    normq=xq**2+yq**2+zq**2+wq**2-1  #require a unit quaternion
    fresult=N.array([eqn1,eqn2,eqn3,eqn4,eqn5,eqn6,normq],'d')
    
    return fresult


def chisq_anneal(p,a,b,c,a2,b2,c2,x0,y0,z0,x0p,y0p,z0p):
#similarly for the eulerian case    
    alpha=p[0]%(2*N.pi)
    beta=p[1]%(2*N.pi)
    gamma=p[2]%(2*N.pi)
    cx=p[3]
    cy=p[4]
    cz=p[5]
    rotx=N.matrix([[1,0,0],[0,N.cos(alpha),N.sin(alpha)],[0,-N.sin(alpha),N.cos(alpha)]],'d')
    roty=N.matrix([[N.cos(beta),0,-N.sin(beta)],[0,1,0],[N.sin(beta),0,N.cos(beta)]],'d')
    rotz=N.matrix([[N.cos(gamma),N.sin(gamma),0],[-N.sin(gamma),N.cos(gamma),0],[0,0,1]],'d')
    temp=rotz*rotx
    rmatrix=temp*roty
    translation=N.matrix([cx,cy,cz],'d').T
    xvec0=N.matrix([x0,y0,z0],'d').T
    xvec0p=transform_point(rmatrix,translation,xvec0)
    eqn4=normal_x(alpha,beta,gamma,cx,cy,cz,a,b,c,a2,b2,c2,x0,y0,z0)
    eqn5=normal_y(alpha,beta,gamma,cx,cy,cz,a,b,c,a2,b2,c2,x0,y0,z0)
    eqn6=normal_z(alpha,beta,gamma,cx,cy,cz,a,b,c,a2,b2,c2,x0,y0,z0)
    eqn7=normal_d(alpha,beta,gamma,cx,cy,cz,a,b,c,a2,b2,c2,x0,y0,z0)
    eqn1=N.ravel(xvec0p[0]-x0p)[0]
    eqn2=N.ravel(xvec0p[1]-y0p)[0]
    eqn3=N.ravel(xvec0p[2]-z0p)[0]
    fresult=N.array([eqn1,eqn2,eqn3,eqn4,eqn5,eqn6],'d')
    chisq=(fresult*fresult).sum()
    return chisq


def genmat(a,b,c,s):
    a11=1-2*b**2-2*c**2
    a12=2*a*b-2*s*c
    a13=2*a*c+2*s*b
    a21=2*a*b+2*s*c
    a22=1-2*a**2-2*c**2
    a23=2*b*c-2*s*a
    a31=2*a*c-2*s*b
    a32=2*b*c+2*s*a
    a33=1-2*a**2-2*b**2
    amat=N.matrix([[a11,a12,a13],[a21,a22,a23],[a31,a32,a33]],'float64')  
    return amat

def genmat2(x,y,z,s):
    t=1-c
    a11=t*x**2+c
    a12=t*x*y+s*z
    a13=t*x*z-s*y
    a21=t*x*y-s*z
    a22=t*y**2+c
    a23=t*y*z*s*x
    a31=t*x*z+s*y
    a32=t*y*z-s*x
    a33=t*z**2+c
    amat=N.matrix([[a11,a12,a13],[a21,a22,a23],[a31,a32,a33]],'float64')  
    return amat



def chisq(p,sx,sy,sz):
    a,b,c,s=p   
    eqn1=2*a*c+2*s*b-sx
    eqn2=2*b*c-2*s*a-sy
    eqn3=1-2*a**2-2*b**2-sz    
    eqn4=N.linalg.det(genmat(a,b,c,s))-1
    eqn5=1-a**2-b**2-c**2-s**2
    fresult=N.array([eqn1,eqn2,eqn3,eqn4,eqn5],'d')
    return fresult

def chisq_web(p,sx,sy,sz):
    x,y,z,s=p
    t=1-c   
    eqn1=t*x*z-s*y-sx
    eqn2=t*y*z+s*x-sy
    eqn3=t*z**2+c-sz    
    eqn4=N.linalg.det(genmat2(a,b,c,s))-1
    eqn5=1-x**2-y**2-z**2-s**2
    fresult=N.array([eqn1,eqn2,eqn3,eqn4,eqn5],'d')
    return fresult

def chisq_f(p,sx,sy,sz):
    a,b,c,s=p   
    eqn1=2*a*c+2*s*b-sx
    eqn2=2*b*c-2*s*a-sy
    eqn3=1-2*a**2-2*b**2-sz    
    eqn4=N.linalg.det(genmat(a,b,c,s))-1
    #eqn5=1-a**2-b**2-c**2-s**2
    fresult=N.array([eqn1,eqn2,eqn3,eqn4],'float64')
    return fresult



def chisq_an(p,sx,sy,sz):
    a,b,c,s=p   
    eqn1=2*a*c+2*s*b-sx
    eqn2=2*b*c-2*s*a-sy
    eqn3=1-2*a**2-2*b**2-sz
    eqn4=N.linalg.det(genmat(a,b,c,s))-1
    eqn5=1-a**2-b**2-c**2-s**2
    fresult=N.array([eqn1,eqn2,eqn3,eqn4,eqn5],'float64')
    chisq=(fresult*fresult).sum()
    return chisq


#def getmatrix(sx,sy,sz,mytol=1e-45,maxiter=2):
#    print "(sx, sy, sz) = ", (sx, sy, sz)
#    
#    p0=N.array([0,1,0,1],'d')
#    p0=N.array([0,0,0,1],'d')
#    #p0=N.array([0,1,0,1],'d')
#    lowerm=[-2,-2,-2,-2]
#    upperm=[2,2,2,2]
#    lowerm=[-1,-1,-1,-1]
#    upperm=[1,1,1,1]
#    flag=True
#    smat=N.empty((3,1),'float64')
#    smat=N.matrix(smat)
#    smat[0]=0
#    smat[1]=0
#    smat[2]=1
#    iter=0
#    pbest=p0
#    currbest=chisq_an(p0,sx,sy,sz)
#    while flag:
#        print 'iter',iter
#        p0,jmin=anneal(chisq_an,p0,args=(sx,sy,sz),\
#                      schedule='simple',lower=lowerm,upper=upperm,\
#                      maxeval=None, maxaccept=None,dwell=10,maxiter=600,T0=10000)
#        
#    
#        
#        print 'p0',p0
#        p = NLSP(chisq, p0,xtol=1e-37)
#        #p = NLP(chisq_an, p0,xtol=1e-17)
#        p.args.f=(sx,sy,sz)
#        p.iprint = 10
#        p.ftol=1e-30
#        p.maxFunEvals = 1e10
#        #p.checkdf()
#        r = p.solve('nssolve')
#        #r = p.solve('scipy_fsolve')
#        #r = p.solve('nlp:ralg')
#        #r = p.solve('ralg')
#        #r = p.solve('scipy_cobyla')
#        print 'solution:', r.xf
#        a,b,c,s=r.xf 
#        if 1:
#            p0=[a,b,c,s]
#            p=scipy.optimize.minpack.fsolve(chisq_f,p0,args=(sx,sy,sz),xtol=1e-37)
#            a,b,c,s=p
#            print '2nd soln',p
# 
#        amat=genmat(a,b,c,s)
#        if 1:
#            array_err=1e-2
#            for i in range(3):
#                for j in range(3):
#                    if N.absolute(amat[i,j])<array_err:
#                        amat[i,j]=0 
#        sout=amat*smat
#        curr_score=chisq_an(p0,sx,sy,sz)
#        if curr_score<currbest:
#            currbest=curr_score
#            pbest=p
#        if (N.abs(sout[0]-sx)<mytol and N.abs(sout[1]-sy)<mytol and N.abs(sout[2]-sz)<mytol) or iter>maxiter:
#            flag=False
#        iter=iter+1
#
#
#    a,b,c,s=pbest    
#    amat=genmat(a,b,c,s)
#    
#    if 1:
#            array_err=1e-5
#            for i in range(3):
#                for j in range(3):
#                    if N.absolute(amat[i,j])<array_err:
#                        amat[i,j]=0 
#    #print 'amat',amat
#    print "AMAT\n\n\nAMAT: \n", amat, "\nSMAT\n" , smat, "\nAMAT * SMAT\n\n", amat*smat, "\nsx, sy, sz\n", (sx,sy,sz)
#    return amat  


def getmatrix2(sx,sy,sz):
    """sx,sy,sz are components of a unit spin vector
    returns rotation matrix mat, such that mat*(0,0,1) = (sx,sy,sz), (both column matrices)"""

    #sx, sy, sz dot product with (0,0,1)=
    c = sz
    #also equal to cos(angle between vectors) becuase they are both unit vectors
    
    #cross product
    x = -sy
    y = sx
    #k = 0
    
    print "\ncross = \n\n", N.cross(N.array([0,0,1]), N.array([sx,sy,sz]))
    print "\n\ndot = \n", N.dot(N.array([0,0,1]), N.array([sx,sy,sz]))
    
    #sin(angle between them) becuase they are both unit vectors:
    s = sqrt(x**2 + y**2)
    
    #a11 = i**2 + (1-i**2)*c
    #a12 = i*j*(1-c)
    #a13 = j*s
    #a21 = i*j*(1-c)
    #a22 = j**2 + (1-j**2)*c
    #a23 = -i*s
    #a31 = -j*s
    #a32 = i*s
    #a33 = c
    
    a11 = 1+(1-c)*(x*x-1)
    a12 = (1-c)*x*y
    a13 = y*s
    a21 = (1-c)*x*y
    a22 = 1 + (1-c)*(y*y-1)
    a23 = -x*s
    a31 = -y*s
    a32 = x*s
    a33 = c
    
    mat = N.matrix([[a11, a12, a13],
                    [a21, a22, a23],
                    [a31, a32, a33]])
    
    return mat
    
    





def getmatrix(sx,sy,sz,mytol=1e-45,maxiter=2):
    print "(sx, sy, sz) = ", (sx, sy, sz)
    
    p0=N.array([0,1,0,1],'d')
    p0=N.array([0,0,0,1],'d')
    #p0=N.array([0,1,0,1],'d')
    lowerm=[-2,-2,-2,-2]
    upperm=[2,2,2,2]
    lowerm=[-1,-1,-1,-1]
    upperm=[1,1,1,1]
    flag=True
    smat=N.empty((3,1),'float64')
    smat=N.matrix(smat)
    smat[0]=0
    smat[1]=0
    smat[2]=1
    iter=0
    pbest=p0
    currbest=chisq_an(p0,sx,sy,sz)
    while flag:
        print 'iter',iter
        p0,jmin=anneal(chisq_an,p0,args=(sx,sy,sz),\
                      schedule='simple',lower=lowerm,upper=upperm,\
                      maxeval=None, maxaccept=None,dwell=10,maxiter=600,T0=10000)
        
    
        
        print 'p0',p0
        p = NLSP(chisq, p0,xtol=1e-37)
        #p = NLP(chisq_an, p0,xtol=1e-17)
        p.args.f=(sx,sy,sz)
        p.iprint = 10
        p.ftol=1e-30
        p.maxFunEvals = 1e10
        #p.checkdf()
        r = p.solve('nssolve')
        #r = p.solve('scipy_fsolve')
        #r = p.solve('nlp:ralg')
        #r = p.solve('ralg')
        #r = p.solve('scipy_cobyla')
        print 'solution:', r.xf
        a,b,c,s=r.xf 
        if 1:
            p0=[a,b,c,s]
            p=scipy.optimize.minpack.fsolve(chisq_f,p0,args=(sx,sy,sz),xtol=1e-37)
            a,b,c,s=p
            print '2nd soln',p
 
        amat=genmat(a,b,c,s)
        if 1:
            array_err=1e-2
            for i in range(3):
                for j in range(3):
                    if N.absolute(amat[i,j])<array_err:
                        amat[i,j]=0 
        sout=amat*smat
        curr_score=chisq_an(p0,sx,sy,sz)
        if curr_score<currbest:
            currbest=curr_score
            pbest=p
        if (N.abs(sout[0]-sx)<mytol and N.abs(sout[1]-sy)<mytol and N.abs(sout[2]-sz)<mytol) or iter>maxiter:
            flag=False
        iter=iter+1


    a,b,c,s=pbest    
    amat=genmat(a,b,c,s)
    
    if 1:
            array_err=1e-5
            for i in range(3):
                for j in range(3):
                    if N.absolute(amat[i,j])<array_err:
                        amat[i,j]=0 
    #print 'amat',amat
    print "AMAT\n\n\nAMAT: \n", amat, "\nSMAT\n" , smat, "\nAMAT * SMAT\n\n", amat*smat, "\nsx, sy, sz\n", (sx,sy,sz)
    return amat





mat1 = getmatrix(-1,0,0)
mat2 = getmatrix2(-1,0,0)

print "\n\nmat1\n", mat1, "\n\nmat2:\n", mat2
col = N.matrix([[0],
                [0],
                [1]])
print mat1*col, "\n\n"
print mat2 * col
print "\n\ndet:\n", N.linalg.det(mat1), "\n\n", N.linalg.det(mat2)
print "\n\n A*At= \n", mat1 * mat1.transpose()
print "\n\n A*At= \n", mat2 * mat2.transpose()

print "\nAt, A^-1\n\n", mat1.T, "\n", mat1.I
print "\nAt, A^-1\n\n", mat2.T, "\n", mat2.I

a = N.array([[JParam(), JParam()],
             [JParam(), JParam()]])
b = N.array([[JParam(), JParam()],
             [JParam(), JParam()]])
print a==b
a[0,1] = JParam(False, 5)
print a


