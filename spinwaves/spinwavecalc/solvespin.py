import numpy as N
import scipy
import math
import scipy.optimize
from utilities.anneal import anneal




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




def getmatrix(sx,sy,sz):
    """sx,sy,sz are components of a unit spin vector
    returns rotation matrix mat, such that mat*(0,0,1) = (sx,sy,sz), (both column matrices)"""

    
    
    
    #sx, sy, sz dot product with (0,0,1)=
    c = sz
    #also equal to cos(angle between vectors) because they are both unit vectors
    
    #cross product
    x = -sy
    y = sx
    #k = 0
    

    
    #print "\ncross = \n\n", N.cross(N.array([0,0,1]), N.array([sx,sy,sz]))
    #print "\n\ndot = \n", N.dot(N.array([0,0,1]), N.array([sx,sy,sz]))
    
    #sin(angle between them) because they are both unit vectors:
    s = math.sqrt(x**2 + y**2)
    
    #make it a unit vector
    #length = math.sqrt(x**2 + y**2)
    x = x/s
    y = y/s
    
    #if the vector (sx, sy, sz) points along the z axis, this method of finding the vector
    #to rotate about will not work because it will give a vector of (0,0,0)
    if (sx == 0 and sy == 0):
        #just rotate about x axis by pi
        return N.matrix([[1,0,0],
                         [0,c,s],
                         [0,-s,c]])
        #x = 1#so we will rotate about the x axis by an angle of pi if sz = -1 or 0 if sz = 1
#        if sz > 0: #sz = 1, theta = 0
#            s = 0
#            c = 1
#        if sz < 0: #sz = -1, theta = pi
#            s = 0
#            c = -1
        #sin(0) = 0 (sz)
        #cos(pi) = -1 (sz)
    
    
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
    
    print "\n\n\nmat: \n", mat
    
    return mat






if __name__=="__main__":
    p0=N.array([0,1,0,1],'d')  #rotation about x-axis
    p0=N.array([0,0,0,1],'d')
    lowerm=[-1,-1,-1,-1]
    upperm=[2,2,2,2]
    sx,sy,sz=[0,0,-1]
#use a simple nonlinear solver for finding A solution to our system of equations 
#    p=scipy.optimize.minpack.fsolve(chisq_quaternion,p0,args=(a,b,c,a2,b2,c2,x0,y0,z0,x0p,y0p,z0p))
    amat=getmatrix(sx,sy,sz)
    if 0:
        #p,jmin=anneal(chisq_an,p0,args=(sx,sy,sz),\
        #              schedule='simple',lower=lowerm,upper=upperm,\
        #              maxeval=None, maxaccept=None,dwell=500,maxiter=2000)
        
        p=scipy.optimize.minpack.fsolve(chisq,p0,args=(sx,sy,sz))
        print p
        a,b,c,s=p
        
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
    if 1:
        print amat
        smat=N.empty((3,1),'float32')
        smat=N.matrix(smat)
        smat[0]=0
        smat[1]=0
        smat[2]=1
        sout=amat*smat
        print sout
        