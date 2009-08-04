import numpy as N
import math
import unittest
eps=1e-3
pi=N.pi

def autovectorized(f):
     """Function decorator to do vectorization only as necessary.
     vectorized functions fail for scalar inputs."""
     def wrapper(input):
          if N.isscalar(input)==False:
               return N.vectorize(f)(input)
          return f(input)
     return wrapper



@autovectorized
def myradians(x):
     return math.radians(x)

@autovectorized
def mydegrees(x):
     return math.degrees(x)
#vecradians = N.vectorize(myradians, otypes=[double])

def sign(x):
     if x>0:
          ret=1
     if x<0:
          ret=-1
     if x==0:
          ret=0
     return ret

def blkdiag(g):
     "Returns a block diagonal matrix, given a list of square matrices, g"
     glen=len(g)
     n=0
     for i in range(glen):
          n=n+g[i].shape[0]


     gout=N.zeros((n,n))
     offset=0
     for i in range(glen):
          currblock=g[i]
          lenx,leny=currblock.shape
          for x in range(lenx):
               for y in range(leny):
                    gout[x+offset,y+offset]=currblock[x,y]
          offset=offset+lenx
     return gout


def similarity_transform(A,B):
     G=N.dot(B,A.transpose())
     G2=N.dot(A,G)
     #
     return G2

def CleanArgs(**args):
     """Takes as input a set of arrays and/or scalars.  It returns a set of the arrays that are the length as the longest entry
    For shorter arrays, it pads them with the last value in the array.  It returns a dictionary based on the calling args"""
     npts=[]
     for name,val in args.iteritems():
          if type(val) in [type(13),type(13.0)]:
               args[name]=N.array([val],'float64')
          if type(val)==type([1,2,3]):
               npts.append(len(val))
          else:
               npts.append(args[name].shape[0])
     maxpts=max(npts)
     for name,val in args.iteritems():
          if type(val)==type([1,2,3]):
               for i in range(maxpts-len(val)):
                    args[name].append(val[-1])
          else:
               if val.shape[0]<maxpts:
                    last_val=val[-1]
                    if len(val.shape)==1:
                         addendum=last_val*N.ones((maxpts-val.shape[0],),'float64')
                         args[name]=N.concatenate((val,addendum))
                    elif len(val.shape)==2:
                         addendum=N.tile(last_val,(maxpts-val.shape[0],1))
                         args[name]=N.concatenate((val,addendum))
               #print name, len(val.shape),args[name]
     return args


class Instrument(object):
     def __init__(self):
          self.tau_list={'pg(002)':1.87325, \
                         'pg(004)':3.74650, \
                         'ge(111)':1.92366, \
                         'ge(220)':3.14131, \
                         'ge(311)':3.68351, \
                         'be(002)':3.50702, \
                         'pg(110)':5.49806}
          return
     def get_tau(self,tau):
          return self.tau_list[tau]
     
class Orientation(object):
     def __init__(self, orient1,orient2):
          self.orient1=orient1
          self.orient2=orient2
     

class Lattice(object):
     
     def _get_a(self): return self._a
     def _get_b(self): return self._b
     def _get_c(self): return self._c
     def _get_alpha(self): return self._alpha
     def _get_beta(self): return self._beta
     def _get_gamma(self): return self._gamma
     #def _get_orient1(self): return self._orient1
     #def _get_orient2(self): return self._orient2
     def _get_orientation(self): return self._orientation
     def _set_a(self,x):
          self._a=x
          self.setvals()
     def _set_b(self,x):
          self._b=x
          self.setvals()
     def _set_c(self,x):
          self._c=x
          self.setvals()
     def _set_alpha(self,x):
          self._alpha=x
          self.setvals()
     def _set_beta(self,x):
          self._beta=x
          self.setvals()
     def _set_gamma(self,x):
          self._gamma=x
          self.setvals()
     #def _set_orient1(self,x):
     #     self._orient1=x
     #     self.setvals()
     #def _set_orient2(self,x):
     #     self._orient2=x
     #     self.setvals()
     def _set_orientation(self,x):
          self._orientation=x
          self._orient1=x.orient1
          self._orient2=x.orient2
          self.setvals()
     
     
     a=property(_get_a,_set_a)
     b=property(_get_b,_set_b)
     c=property(_get_c,_set_c)
     alpha=property(_get_alpha,_set_alpha)
     beta=property(_get_beta,_set_beta)
     gamma=property(_get_gamma,_set_gamma)
     #orient1=property(_get_orient1,_set_orient1)
     #orient2=property(_get_orient2,_set_orient2)
     orientation=property(_get_orientation, _set_orientation)
     
     def setvals(self):
          #if self._orient1.shape[0]==1:
          #     self.orient1=self.orient1.transpose()
          #if self._orient2.shape[0]==1:
          #     self.orient2=self.orient2.transpose()

          
          newinput=CleanArgs(a=self._a,
                             b=self._b,
                             c=self._c,
                             alpha=self._alpha,
                             beta=self._beta,
                             gamma=self._gamma,
                             orient1=self._orient1,
                             orient2=self._orient2)
          self._a=newinput['a']
          self._b=newinput['b']
          self._c=newinput['c']
          self._alpha=newinput['alpha']
          self._beta=newinput['beta']
          self._gamma=newinput['gamma']
          self.alphad=mydegrees(newinput['alpha'])
          self.betad=mydegrees(newinput['beta'])
          self.gammad=mydegrees(newinput['gamma'])
          self.star()
          self.gtensor('lattice')
          self.gtensor('latticestar')
          self.npts=N.size(self.a)
          self._orient1=newinput['orient1']
          self._orient2=newinput['orient2']
          #self._orientation.orient1=newinput['orient1']
          #self._orientation.orient2=newinput['orient2']
          
          #print self.npts
          
          #self._orient1=self._orient1.reshape((self.npts,1))
          #self._orient2=self._orient2.reshape((self.npts,1))
         
          #if newinput['orient1'].shape[0]==1:
          #     self.orient1=newinput['orient1'].transpose()
          #else:
          #     self.orient1=newinput['orient1']
          #if newinput['orient2'].shape[0]==1:
          #     self.orient2=newinput['orient2'].transpose()
          #else:
          #     self.orient2=newinput['orient2']

               
          self.StandardSystem()
          
     def __init__(self, a=None, \
                  b=None, \
                  c=None, \
                  alpha=None, \
                  beta=None, \
                  gamma=None, \
                  orientation=None, \
                  ):
          """a,b,c in Angstroms, alpha,beta, gamma in radians.  All are vectors """
          self._a=a
          self._b=b
          self._c=c
          self._alpha=alpha
          self._beta=beta
          self._gamma=gamma
          self._orient1=orientation.orient1
          self._orient2=orientation.orient2
          self.instrument=Instrument()
          self.setvals()
          return



     def scalar(self, x1, y1, z1, x2, y2, z2, lattice):
          "calculates scalar product of two vectors"
          if lattice=='lattice':
               a=self.a
               b=self.b
               c=self.c
               alpha=self.alpha
               beta=self.beta
               gamma=self.gamma
          if lattice=='latticestar':
               a=self.astar
               b=self.bstar
               c=self.cstar
               alpha=self.alphastar
               beta=self.betastar
               gamma=self.gammastar

          s=x1*x2*a**2+y1*y2*b**2+z1*z2*c**2+\
           (x1*y2+x2*y1)*a*b*N.cos(gamma)+\
           (x1*z2+x2*z1)*a*c*N.cos(beta)+\
           (z1*y2+z2*y1)*c*b*N.cos(alpha)
          return s

     def star(self):
          "Calculate unit cell volume, reciprocal cell volume, reciprocal lattice parameters"
          V=2*self.a*self.b*self.c*\
           N.sqrt(N.sin((self.alpha+self.beta+self.gamma)/2)*\
                  N.sin((-self.alpha+self.beta+self.gamma)/2)*\
                  N.sin((self.alpha-self.beta+self.gamma)/2)*\
                  N.sin((self.alpha+self.beta-self.gamma)/2))
          self.Vstar=(2*N.pi)**3/V;
          self.astar=2*N.pi*self.b*self.c*N.sin(self.alpha)/V
          self.bstar=2*N.pi*self.a*self.c*N.sin(self.beta)/V
          self.cstar=2*N.pi*self.b*self.a*N.sin(self.gamma)/V
          self.alphastar=N.arccos((N.cos(self.beta)*N.cos(self.gamma)-\
                                   N.cos(self.alpha))/ \
                                  (N.sin(self.beta)*N.sin(self.gamma)))
          self.betastar= N.arccos((N.cos(self.alpha)*N.cos(self.gamma)-\
                                   N.cos(self.beta))/ \
                                  (N.sin(self.alpha)*N.sin(self.gamma)))
          self.gammastar=N.arccos((N.cos(self.alpha)*N.cos(self.beta)-\
                                   N.cos(self.gamma))/ \
                                  (N.sin(self.alpha)*N.sin(self.beta)))
          self.V=V
          return

     def angle2(self, x, y, z, h, k, l):
          "Calculate the angle between vectors in real and reciprocal space"
          "x,y,z are the fractional cell coordinates of the first vector,"
          "h,k,l are Miller indices of the second vector"
          phi=N.arccos(2*pi*(h*x+k*y+l*z)/self.modvec(x, y, z, 'lattice')/self.modvec(h, k, l, 'latticestar'))
          return phi


     def angle(self, x1, y1, z1, x2, y2, z2, lattice):
          "Calculate the angle between vectors in real and reciprocal space"
          "xi,yi,zi are the fractional cell coordinates of the vectors"
          phi=N.arccos(self.scalar(x1, y1, z1, x2, y2, z2, lattice)/self.modvec(x1, y1, z1, lattice)/self.modvec(x1, y1, z1, lattice))
          return phi

     def modvec(self, x, y, z, lattice):
          "Calculates modulus of a vector defined by its fraction cell coordinates"
          "or Miller indexes"
          m=N.sqrt(self.scalar(x, y, z, x, y, z, lattice))
          return m

     def gtensor(self, lattice):
          "calculates the metric tensor of a lattice"
          g=N.zeros((3, 3, N.size(self.a)), 'Float64')
          #print 'shape ', g.shape
          if lattice=='lattice':
               a=self.a
               b=self.b
               c=self.c
               alpha=self.alpha
               beta=self.beta
               gamma=self.gamma
          if lattice=='latticestar':
               a=self.astar
               b=self.bstar
               c=self.cstar
               alpha=self.alphastar
               beta=self.betastar
               gamma=self.gammastar
          g[0, 0, :]=a**2;
          g[0, 1, :]=a*b*N.cos(gamma)
          g[0, 2, :]=a*c*N.cos(beta)

          g[1, 0, :]=g[0, 1, :]
          g[1, 1, :]=b**2
          g[1, 2, :]=c*b*N.cos(alpha)

          g[2, 0, :]=g[0, 2, :]
          g[2, 1, :]=g[1, 2, :]
          g[2, 2, :]=c**2
          if lattice=='lattice':
               self.g=g
          if lattice=='latticestar':
               self.gstar=g
          return

     def reciprocate(self, x, y, z, lattice):
          "calculate miller indexes of a vector defined by its fractional cell coords"
          if lattice=='lattice':
               g=self.g
          if lattice=='latticestar':
               g=self.gstar
          h=g[0, 0, :]*x+g[1, 0, :]*y+g[2, 0, :]*z;
          k=g[0, 1, :]*x+g[1, 1, :]*y+g[2, 1, :]*z;
          l=g[0, 2, :]*x+g[1, 2, :]*y+g[2, 2, :]*z;
          return h, k, l

     def vector(self, x1, y1, z1, x2, y2, z2, lattice):
          "calculates the fractional cell coordinates or Miller indexes of a vector"
          "product of two vectors, defined by their fractional cell coordinates or "
          "Miller idexes"
          if lattice=='lattice':
               g=self.gstar
               V=self.Vstar
          if lattice=='latticestar':
               g=self.g
               V=self.V
          g=g*V/(2*N.pi)**2
          x=y1*z2*g[0, 0, :]-z1*y2*g[0, 0, :]-x1*z2*g[1, 0, :]+z1*x2*g[1, 0, :]\
           +x1*y2*g[2, 0, :]-y1*x2*g[2, 0, :]
          y=y1*z2*g[0, 1, :]-z1*y2*g[0, 1, :]-x1*z2*g[1, 1, :]+z1*x2*g[1, 1, :]\
           +x1*y2*g[2, 1, :]-y1*x2*g[2, 1, :]
          z=y1*z2*g[0, 2, :]-z1*y2*g[0, 2, :]-x1*z2*g[1, 2, :]+z1*x2*g[1, 2, :]\
           +x1*y2*g[2, 2, :]-y1*x2*g[2, 2, :]
          return x,y,z

     def StandardSystem(self):
          orient1=self._orient1
          orient2=self._orient2
          try:
               modx=self.modvec(orient1[0, :], orient1[1, :], orient1[2, :], 'latticestar')
          except IndexError:
               orient1=orient1.transpose()
               orient2=orient2.transpose()
               modx=self.modvec(orient1[0, :], orient1[1, :], orient1[2, :], 'latticestar')
          x=N.copy(orient1)
          x[0, :]=x[0, :]/modx; # First unit basis vector
          x[1, :]=x[1, :]/modx;
          x[2, :]=x[2, :]/modx;

          proj=self.scalar(orient2[0, :], orient2[1, :], orient2[2, :], \
                           x[0, :], x[1, :], x[2, :], 'latticestar')

          y=N.copy(orient2)
          y[0, :]=y[0, :]-x[0, :]*proj;
          y[1, :]=y[1, :]-x[1, :]*proj;
          y[2, :]=y[2, :]-x[2, :]*proj;

          mody=self.modvec(y[0, :], y[1, :], y[2, :], 'latticestar');

#    check for collinearity of orienting vectors

          try:
               if N.where(mody<=eps)[0].size>0:
                    print 'ValueError'
                    raise ValueError
               y[0, :]=y[0, :]/mody; # Second unit basis vector
               y[1, :]=y[1, :]/mody;
               y[2, :]=y[2, :]/mody;

               z=N.copy(y);

               z[0, :]=x[1, :]*y[2, :]-y[1, :]*x[2, :];
               z[1, :]=x[2, :]*y[0, :]-y[2, :]*x[0, :];
               z[2, :]=-x[1, :]*y[0, :]+y[1, :]*x[0, :];

               proj=self.scalar(z[0, :], z[1, :], z[2, :], x[0, :], x[1, :], x[2, :], 'latticestar');

               z[0, :]=z[0, :]-x[0, :]*proj;
               z[1, :]=z[1, :]-x[1, :]*proj;
               z[2, :]=z[2, :]-x[2, :]*proj;

               proj=self.scalar(z[0, :], z[1, :], z[2, :], y[0, :], y[1, :], y[2, :], 'latticestar');

               z[0, :]=z[0, :]-y[0, :]*proj;
               z[1, :]=z[1, :]-y[1, :]*proj;
               z[2, :]=z[2, :]-y[2, :]*proj;

               modz=self.modvec(z[0, :], z[1, :], z[2, :], 'latticestar');

               z[0, :]=z[0, :]/modz; #% Third unit basis vector
               z[1, :]=z[1, :]/modz;
               z[2, :]=z[2, :]/modz;

               self.x=x
               self.y=y
               self.z=z
          except ValueError:
               print 'ORIENTATION VECTORS ARE COLLINEAR x,y,z not set'
          return

     def S2R(self, qx, qy, qz):
          "Given cartesian coordinates of a vector in the S System, calculate its Miller indexes."
          x=self.x
          y=self.y
          z=self.z
          H=qx*x[0, :]+qy*y[0, :]+qz*z[0, :];
          K=qx*x[1, :]+qy*y[1, :]+qz*z[1, :];
          L=qx*x[2, :]+qy*y[2, :]+qz*z[2, :];
          q=N.sqrt(qx**2+qy**2+qz**2);
          return H, K, L, q

     def R2S(self, H, K, L):
          "Given reciprocal-space coordinates of a vector, calculate its coordinates in the Cartesian space."
          x=self.x
          y=self.y
          z=self.z
          qx=self.scalar(H, K, L, x[0, :], x[1, :], x[2, :], 'latticestar');
          qy=self.scalar(H, K, L, y[0, :], y[1, :], y[2, :], 'latticestar');
          qz=self.scalar(H, K, L, z[0, :], z[1, :], z[2, :], 'latticestar');
          q=self.modvec(H, K, L, 'latticestar');
          return qx, qy, qz, q

     def SpecWhere(self,M2,S1,S2,A2,EXP):
          """ For given values of M3,S1,S2 and A2 spectrometer motors (AKA M2,M3,M4 and M6)
        and spectrometer and sample parameters specified in EXP calculates the wave vector
        transfer in the sample (H, K, L), Q=|(H,K,L)|, energy tranfer E, and incident
        and final neutron energies.  Angles are given in radians"""
          newinput=CleanArgs(a=self.a,b=self.b,c=self.c,alpha=self.alpha,beta=self.beta,\
                             gamma=self.gamma,orient1=self._orient1,orient2=self._orient2,M2=M2,S1=S1,S2=S2,A2=A2)
          orientation=Orientation(newinput['orient1'],newinput['orient2'])
          self.__init__(a=newinput['a'],b=newinput['b'],c=newinput['c'],alpha=newinput['alpha'],\
                        beta=newinput['beta'],gamma=newinput['gamma'],orientation=orientation\
                        )
          M2=newinput['M2']
          S1=newinput['S1']
          S2=newinput['S2']
          A2=newinput['A2']
          npts=len(EXP)
          taum=N.empty(npts,'Float64')
          taua=N.empty(npts,'Float64')
          for ind in range(npts):
               taum[ind]=self.instrument.get_tau(EXP[ind]['mono']['tau'])
          for ind in range(npts):
               taua[ind]=self.instrument.get_tau(EXP[ind]['ana']['tau'])
          ki=taum/N.sqrt(2.0-2*N.cos(M2))
          Ei=2.072142*ki**2
          kf=taua/N.sqrt(2.0-2*N.cos(A2))
          Ef=2.072142*kf**2
          E=Ei-Ef
          Q=N.sqrt(ki**2+kf**2-2*ki*kf*N.cos(S2))
          orienta=self.x
          orientb=self.y
          #phi=-atan2(-kf.*sin(S2), ki-kf.*cos(S2)); %Angle from ki to Q
          delta=N.absolute(N.arccos( (Q**2+ki**2-kf**2)/(2*ki*Q)))
          psi=S1+delta-pi/2 #Angle from first orienting vector to to Q
          qx=Q*N.cos(psi)
          qy=Q*N.sin(psi)
          H=qx*orienta[0]+qy*orientb[0]
          K=qx*orienta[1]+qy*orientb[1]
          L=qx*orienta[2]+qy*orientb[2]
          return H,K,L,E,Q,Ei,Ef


     def SpecGoTo(self,H,K,L,E,EXP):
          """Calculate shaft angles given momentum transfer H, K, L, energy transfer E, and
          spectrometer and smaple parameters in EXP.  The angles returned are in radians"""
          newinput=CleanArgs(a=self.a,b=self.b,c=self.c,alpha=self.alpha,beta=self.beta,\
                             gamma=self.gamma,orient1=self._orient1,orient2=self._orient2,H=H,K=K,L=L,E=E)
          orientation=Orientation(newinput['orient1'],newinput['orient2'])
          self.__init__(a=newinput['a'],b=newinput['b'],c=newinput['c'],alpha=newinput['alpha'],\
                        beta=newinput['beta'],gamma=newinput['gamma'],orientation=orientation\
                        )
          H=newinput['H']
          K=newinput['K']
          L=newinput['L']
          E=newinput['E']

          CONVERT2=2.072
          #mono=[EXP['mono']]
          npts=len(EXP)
          taum=N.empty((npts,),'float64')
          taua=N.empty((npts,),'float64')
          infin=-1*N.ones((npts,1),'float64')
          dir1=N.ones((npts,1),'float64')
          dir2=N.ones((npts,1),'float64')
          mondir=N.ones((npts,1),'float64')
          efixed=N.empty((npts,1),'float64')
          for ind in range(npts):
               taum[ind]=self.instrument.get_tau(EXP[ind]['mono']['tau'])
               taua[ind]=self.instrument.get_tau(EXP[ind]['ana']['tau'])
               if ('infin' in EXP[ind]):
                    infin[ind]=EXP[ind]['infin']
               if ('dir1' in EXP[ind]):
                    dir1[ind]=EXP[ind]['dir1']
               if ('dir2' in EXP):
                    dir2[ind]=EXP[ind]['dir2']
               if ('mondir' in EXP):
                    mondir[ind]=EXP[ind]['mondir']
               efixed[ind]=EXP[ind]['efixed']

          qx,qy,qz,Q=self.R2S(H,K,L)
          dir=N.zeros((3,npts),'float64')
          dir[0,:]=mondir
          dir[1,:]=-dir[0,:]*dir1
          dir[2,:]=-dir[1,:]*dir2

          ei=efixed+E
          ef=efixed+0*E
          change=N.where(infin>0)
          if N.size(change)!=0:
               ef[change]=efixed[change]-E[change]
               ei[change]=efixed[change]

          ki = N.sqrt(ei/CONVERT2)
          kf = N.sqrt(ef/CONVERT2)

          M1=N.arcsin(taum/(2*ki))#.*sign(dir(1,:));
          M2=2*M1
          A1=N.arcsin(taua/(2*kf))#.*sign(dir(3,:));
          A2=2*A1
          S2=N.arccos((ki**2+kf**2-Q**2)/(2*ki*kf))#.*sign(dir(2,:));
          delta=N.absolute(N.arccos( (Q**2+ki**2-kf**2)/(2*ki*Q)))
          #psi=S1+delta-pi/2 #Angle from first orienting vector to to Q
          psi=N.arctan2(qy,qx)
          S1=psi-delta+pi/2
          #TODO:  Add checks to make sure that the scattering triangle closed
##        bad=find(ei<0 | ef<0 | abs(taum./(2.*ki))>1 | abs(taua./(2.*kf))>1 | abs ( (ki.^2+kf.^2-Q.^2)./(2*ki.*kf))>1);
##        M1(bad)=NaN;
##        M2(bad)=NaN;
##        S1(bad)=NaN;
##        S2(bad)=NaN;
##        A1(bad)=NaN;
##        A2(bad)=NaN;
          return M1,M2,S1,S2,A1,A2

class TestLattice(unittest.TestCase):

     def setUp(self):
          a=N.array([2*pi],'Float64')
          b=N.array([2*pi],'Float64')
          c=N.array([2*pi],'Float64')
          alpha=myradians(N.array([90],'Float64'))
          beta=myradians(N.array([90],'Float64'))
          gamma=myradians(N.array([90],'Float64'))
          orient1=N.array([[1,0,0]],'Float64')
          orient2=N.array([[0,1,1]],'Float64')
          orientation=Orientation(orient1,orient2)
          self.fixture = Lattice(a=a,b=b,c=c,alpha=alpha,beta=beta,gamma=gamma,\
                                 orientation=orientation)

     def test_astar(self):
          self.assertAlmostEqual(self.fixture.astar[0],1.0,2,'astar Not equal to '+str(1.0))
     def test_bstar(self):
          self.assertAlmostEqual(self.fixture.bstar[0],1.0,2,'bstar Not equal to '+str(1.0))
     def test_cstar(self):
          self.assertAlmostEqual(self.fixture.cstar[0],1.0,2,'cstar '+str(self.fixture.cstar[0])+' Not equal to '+str(1.0))
     def test_alphastar(self):
          self.assertAlmostEqual(self.fixture.alphastar[0],pi/2,2,'alphastar Not equal to '+str(pi/2))
     def test_betastar(self):
          self.assertAlmostEqual(self.fixture.betastar[0],pi/2,2,'betastar Not equal to '+str(pi/2))
     def test_gammastar(self):
          self.assertAlmostEqual(self.fixture.gammastar[0],pi/2,2,'gammastar Not equal to '+str(pi/2))
     def test_V(self):
          self.assertAlmostEqual(self.fixture.V[0],248.0502,2,'V Not equal to '+str(248.0502))
     def test_Vstar(self):
          self.assertAlmostEqual(self.fixture.Vstar[0],1.0,2,'Vstar Not equal to '+str(1.0))
     def test_g(self):
          #print self.fixture.g
          self.assertAlmostEqual((self.fixture.g[:,:,0][0,0]),39.4784*(N.eye(3)[0,0]) ,2,'g Not equal to '+str(39.4784 ))
     def test_gstar(self):
          #print self.fixture.gstar
          self.assertAlmostEqual(self.fixture.gstar[:,:,0][0,0],1.0*N.eye(3)[0,0] ,2,'gstar Not equal to '+str(1.0 ))

     def test_StandardSystem_x(self):
     #       #print self.fixture.gstar
          self.assertAlmostEqual(self.fixture.x[0],1.0 ,2,'Standard System x Not equal to '+str(1.0 ))


class TestLatticeCubic(unittest.TestCase):

     def setUp(self):
          a=N.array([6.283],'Float64')
          b=N.array([6.283],'Float64')
          c=N.array([6.283],'Float64')
          alpha=myradians(N.array([90],'Float64'))
          beta=myradians(N.array([90],'Float64'))
          gamma=myradians(N.array([90],'Float64'))
          orient1=N.array([[1,0,0]],'Float64')
          orient2=N.array([[0,0,1]],'Float64')
          orientation=Orientation(orient1,orient2)
          self.fixture=Lattice(a=a,b=b,c=c,alpha=alpha,beta=beta,gamma=gamma,\
                                 orientation=orientation)
          EXP={}
          EXP['ana']={}
          EXP['ana']['tau']='pg(002)'
          EXP['mono']={}
          EXP['mono']['tau']='pg(002)';
          EXP['ana']['mosaic']=30
          EXP['mono']['mosaic']=30
          EXP['sample']={}
          EXP['sample']['mosaic']=10
          EXP['sample']['vmosaic']=10
          EXP['hcol']=N.array([40, 10, 20, 80],'Float64')
          EXP['vcol']=N.array([120, 120, 120, 120],'Float64')
          EXP['infix']=-1 #positive for fixed incident energy
          EXP['efixed']=14.7
          EXP['method']=0
          setup=[EXP]
          self.fixture.EXP=EXP

     def test_cubic1(self):
          
          #setup lattice
          orient1=N.array([[1,0,0]],'Float64')
          orient2=N.array([[0,0,1]],'Float64')
          orientation=Orientation(orient1,orient2)
          self.fixture.orientation=orientation
          
          #setup spectrometer
          self.fixture.EXP['infix']=1 # Fixed Ef
          self.fixture.EXP['efixed']=5.0
          
          #test the angles
          M2=myradians(N.array([74.169]))
          A2=myradians(N.array([74.169]))
          S1=myradians(N.array([97.958]))
          S2=myradians(N.array([89.131]))
          H,K,L,E,Q,Ei,Ef=self.fixture.SpecWhere(M2,S1,S2,A2,[self.fixture.EXP])
          print 'H ',H
          print 'K ',K
          print 'L ',L
          print 'E ',E
          print 'Q ',Q
          print 'Ei ',Ei
          print 'Ef ',Ef
          self.assertAlmostEqual(H[0],1.2999,4, 'H Not equal to '+ str(1.2999))
          self.assertAlmostEqual(K[0],0.0000,4, 'K Not equal to '+ str(0.0000))
          self.assertAlmostEqual(L[0],1.7499,4, 'L Not equal to '+ str(1.7499))
          self.assertAlmostEqual(E[0],0.0000,4, 'E Not equal to '+ str(0.0000))
          #self.assertAlmostEqual(Q[0],2.1799,4, 'Q Not equal to '+ str(2.1799))
          self.assertAlmostEqual(Ei[0],4.9995,4,'Ei Not equal to '+str(4.9995))
          self.assertAlmostEqual(Ef[0],4.9995,4,'Ef Not equal to '+str(4.9995))

          
     def test_cubic2(self):
          """test different Energy Transfer"""
          #setup lattice
          orient1=N.array([[1,0,0]],'Float64')
          orient2=N.array([[0,0,1]],'Float64')
          orientation=Orientation(orient1,orient2)
          self.fixture.orientation=orientation
          
          #setup spectrometer
          self.fixture.EXP['infix']=1 # Fixed Ef
          self.fixture.EXP['efixed']=5.0
          
          #test the angles
          M2=myradians(N.array([52.420]))
          A2=myradians(N.array([74.169]))
          S1=myradians(N.array([101.076]))
          S2=myradians(N.array([70.881]))
          H,K,L,E,Q,Ei,Ef=self.fixture.SpecWhere(M2,S1,S2,A2,[self.fixture.EXP])
          print 'H ',H
          print 'K ',K
          print 'L ',L
          print 'E ',E
          print 'Q ',Q
          print 'Ei ',Ei
          print 'Ef ',Ef
          self.assertAlmostEqual(H[0],1.2999,4, 'H Not equal to '+ str(1.2999))
          self.assertAlmostEqual(K[0],0.0000,4, 'K Not equal to '+ str(0.0000))
          self.assertAlmostEqual(L[0],1.7499,4, 'L Not equal to '+ str(1.7499))
          self.assertAlmostEqual(E[0],4.3195,4, 'E Not equal to '+ str(4.3195))
          #self.assertAlmostEqual(Q[0],2.1799,4, 'Q Not equal to '+ str(2.1799))
          self.assertAlmostEqual(Ei[0],9.3190,4,'Ei Not equal to '+str(9.3190))
          self.assertAlmostEqual(Ef[0],4.9995,4,'Ef Not equal to '+str(4.9995))          
          
     def test_cubic3(self):
          """test different Orientation"""
          #setup lattice
          orient1=N.array([[1,1,0]],'Float64')
          orient2=N.array([[0,0,1]],'Float64')
          orientation=Orientation(orient1,orient2)
          self.fixture.orientation=orientation
          
          #setup spectrometer
          self.fixture.EXP['infix']=1 # Fixed Ef
          self.fixture.EXP['efixed']=5.0
          
          #test the angles
          M2=myradians(N.array([74.169]))
          A2=myradians(N.array([74.169]))
          S1=myradians(N.array([98.375]))
          S2=myradians(N.array([109.575]))
          H,K,L,E,Q,Ei,Ef=self.fixture.SpecWhere(M2,S1,S2,A2,[self.fixture.EXP])
          print 'H ',H
          print 'K ',K
          print 'L ',L
          print 'E ',E
          print 'Q ',Q
          print 'Ei ',Ei
          print 'Ef ',Ef
          self.assertAlmostEqual(H[0],1.2999,4, 'H Not equal to '+ str(1.2999))
          self.assertAlmostEqual(K[0],1.2999,4, 'K Not equal to '+ str(1.2999))
          self.assertAlmostEqual(L[0],1.7499,4, 'L Not equal to '+ str(1.7499))
          self.assertAlmostEqual(E[0],0.0000,4, 'E Not equal to '+ str(0.0000))
          self.assertAlmostEqual(Ei[0],4.9995,4,'Ei Not equal to '+str(4.9995))
          self.assertAlmostEqual(Ef[0],4.9995,4,'Ef Not equal to '+str(4.9995))          
          
     def test_cubic4(self):
          """Switch order of orientations, compare with test 3"""
          #setup lattice
          orient1=N.array([[0,0,1]],'Float64')
          orient2=N.array([[1,1,0]],'Float64')
          orientation=Orientation(orient1,orient2)
          self.fixture.orientation=orientation
          
          #setup spectrometer
          self.fixture.EXP['infix']=1 # Fixed Ef
          self.fixture.EXP['efixed']=5.0
          
          #test the angles
          M2=myradians(N.array([74.169]))
          A2=myradians(N.array([74.169]))
          S1=myradians(N.array([101.200]))  #Note that this angle has changed
          #This is a consequence of the fact that ICP defines the first orientation vector 
          #to be at half of the detector two theta angle.  The second orientation vector is always
          #at higher a3
          S2=myradians(N.array([109.575]))
          H,K,L,E,Q,Ei,Ef=self.fixture.SpecWhere(M2,S1,S2,A2,[self.fixture.EXP])
          print 'H ',H
          print 'K ',K
          print 'L ',L
          print 'E ',E
          print 'Q ',Q
          print 'Ei ',Ei
          print 'Ef ',Ef
          self.assertAlmostEqual(H[0],1.2999,4, 'H Not equal to '+ str(1.2999))
          self.assertAlmostEqual(K[0],1.2999,4, 'K Not equal to '+ str(1.2999))
          self.assertAlmostEqual(L[0],1.7499,4, 'L Not equal to '+ str(1.7499))
          self.assertAlmostEqual(E[0],0.0000,4, 'E Not equal to '+ str(0.0000))
          self.assertAlmostEqual(Ei[0],4.9995,4,'Ei Not equal to '+str(4.9995))
          self.assertAlmostEqual(Ef[0],4.9995,4,'Ef Not equal to '+str(4.9995))          

     def test_cubic5(self):
          """Test another energy, compare with test 4"""
          #setup lattice
          orient1=N.array([[0,0,1]],'Float64')
          orient2=N.array([[1,1,0]],'Float64')
          orientation=Orientation(orient1,orient2)
          self.fixture.orientation=orientation
          
          #setup spectrometer
          self.fixture.EXP['infix']=1 # Fixed Ef
          self.fixture.EXP['efixed']=5.0
          
          #test the angles
          M2=myradians(N.array([48.661]))
          A2=myradians(N.array([74.169]))
          S1=myradians(N.array([99.257]))  
          S2=myradians(N.array([80.722]))
          H,K,L,E,Q,Ei,Ef=self.fixture.SpecWhere(M2,S1,S2,A2,[self.fixture.EXP])
          print 'H ',H
          print 'K ',K
          print 'L ',L
          print 'E ',E
          print 'Q ',Q
          print 'Ei ',Ei
          print 'Ef ',Ef
          self.assertAlmostEqual(H[0],1.2999,4, 'H Not equal to '+ str(1.2999))
          self.assertAlmostEqual(K[0],1.2999,4, 'K Not equal to '+ str(1.2999))
          self.assertAlmostEqual(L[0],1.7499,4, 'L Not equal to '+ str(1.7499))
          self.assertAlmostEqual(E[0],5.7097,4, 'E Not equal to '+ str(5.7097))
          self.assertAlmostEqual(Ei[0],10.7092,4,'Ei Not equal to '+str(10.7092))
          self.assertAlmostEqual(Ef[0],4.9995,4,'Ef Not equal to '+str(4.9995))    
          
          
     def test_cubic6(self):
          """test different energies, compare with test 3"""
          #setup lattice
          orient1=N.array([[1,1,0]],'Float64')
          orient2=N.array([[0,0,1]],'Float64')
          orientation=Orientation(orient1,orient2)
          self.fixture.orientation=orientation
          
          #setup spectrometer
          self.fixture.EXP['infix']=1 # Fixed Ef
          self.fixture.EXP['efixed']=5.0
          
          #test the angles
          M2=myradians(N.array([48.661],'Float64'))
          A2=myradians(N.array([74.169],'Float64'))
          S1=myradians(N.array([96.433],'Float64'))
          S2=myradians(N.array([80.722],'Float64'))
          H,K,L,E,Q,Ei,Ef=self.fixture.SpecWhere(M2,S1,S2,A2,[self.fixture.EXP])
          print 'H ',H
          print 'K ',K
          print 'L ',L
          print 'E ',E
          print 'Q ',Q
          print 'Ei ',Ei
          print 'Ef ',Ef
          self.assertAlmostEqual(H[0],1.2999,4, 'H Not equal to '+ str(1.2999))
          self.assertAlmostEqual(K[0],1.2999,4, 'K Not equal to '+ str(1.2999))
          self.assertAlmostEqual(L[0],1.7499,4, 'L Not equal to '+ str(1.7499))
          self.assertAlmostEqual(E[0],5.7097,4, 'E Not equal to '+ str(5.7097))
          self.assertAlmostEqual(Ei[0],10.7092,4,'Ei Not equal to '+str(10.7092))
          self.assertAlmostEqual(Ef[0],4.9995,4,'Ef Not equal to '+str(4.9995))   
          
     def test_tetragonal1(self):
          """test the tetragonal cell"""
          #setup lattice
          self.fixture.a=N.array([6.283],'Float64')
          self.fixture.b=N.array([6.283],'Float64')
          self.fixture.c=N.array([11.765],'Float64')
          orient1=N.array([[1,0,0]],'Float64')
          orient2=N.array([[0,0,1]],'Float64')
          orientation=Orientation(orient1,orient2)
          self.fixture.orientation=orientation
          
          #setup spectrometer
          self.fixture.EXP['infix']=1 # Fixed Ef
          self.fixture.EXP['efixed']=5.0
          
          #test the angles
          M2=myradians(N.array([74.169],'Float64'))
          A2=myradians(N.array([74.169],'Float64'))
          S1=myradians(N.array([76.720],'Float64'))
          S2=myradians(N.array([110.480],'Float64'))
          H,K,L,E,Q,Ei,Ef=self.fixture.SpecWhere(M2,S1,S2,A2,[self.fixture.EXP])
          print 'H ',H
          print 'K ',K
          print 'L ',L
          print 'E ',E
          print 'Q ',Q
          print 'Ei ',Ei
          print 'Ef ',Ef
          self.assertAlmostEqual(H[0],2.3749,4, 'H Not equal to '+ str(2.3749))
          self.assertAlmostEqual(K[0],0.0000,4, 'K Not equal to '+ str(0.0000))
          self.assertAlmostEqual(L[0],1.7499,4, 'L Not equal to '+ str(1.7499))
          self.assertAlmostEqual(E[0],0.0000,4, 'E Not equal to '+ str(0.0000))
          self.assertAlmostEqual(Ei[0],4.9995,4,'Ei Not equal to '+str(4.9995))
          self.assertAlmostEqual(Ef[0],4.9995,4,'Ef Not equal to '+str(4.9995))   
          

          
     def test_tetragonal2(self):
          """test the tetragonal cell, change E transfer"""
          #setup lattice
          self.fixture.a=N.array([6.283],'Float64')
          self.fixture.b=N.array([6.283],'Float64')
          self.fixture.c=N.array([11.765],'Float64')
          orient1=N.array([[1,0,0]],'Float64')
          orient2=N.array([[0,0,1]],'Float64')
          orientation=Orientation(orient1,orient2)
          self.fixture.orientation=orientation
          
          #setup spectrometer
          self.fixture.EXP['infix']=1 # Fixed Ef
          self.fixture.EXP['efixed']=5.0
          
          #test the angles
          M2=myradians(N.array([49.633],'Float64'))
          A2=myradians(N.array([74.169],'Float64'))
          S1=myradians(N.array([74.345],'Float64'))
          S2=myradians(N.array([82.717],'Float64'))
          H,K,L,E,Q,Ei,Ef=self.fixture.SpecWhere(M2,S1,S2,A2,[self.fixture.EXP])
          print 'H ',H
          print 'K ',K
          print 'L ',L
          print 'E ',E
          print 'Q ',Q
          print 'Ei ',Ei
          print 'Ef ',Ef
          self.assertAlmostEqual(H[0],2.3749,4, 'H Not equal to '+ str(2.3749))
          self.assertAlmostEqual(K[0],0.0000,4, 'K Not equal to '+ str(0.0000))
          self.assertAlmostEqual(L[0],1.7499,4, 'L Not equal to '+ str(1.7499))
          self.assertAlmostEqual(E[0],5.3197,4, 'E Not equal to '+ str(5.3197))
          self.assertAlmostEqual(Ei[0],10.3192,4,'Ei Not equal to '+str(10.3192))
          self.assertAlmostEqual(Ef[0],4.9995,4,'Ef Not equal to '+str(4.9995))   


     def test_tetragonal3(self):
          """test the tetragonal cell, change Ei, orientation vectors"""
          #setup lattice
          self.fixture.a=N.array([6.283],'Float64')
          self.fixture.b=N.array([6.283],'Float64')
          self.fixture.c=N.array([11.765],'Float64')
          orient1=N.array([[1,1,0]],'Float64')
          orient2=N.array([[0,0,1]],'Float64')
          orientation=Orientation(orient1,orient2)
          self.fixture.orientation=orientation
          
          #setup spectrometer
          self.fixture.EXP['infix']=-1 # Fixed Ei
          self.fixture.EXP['efixed']=3.7
          
          #test the angles
          M2=myradians(N.array([89.008],'Float64'))
          A2=myradians(N.array([89.008],'Float64'))
          S1=myradians(N.array([100.569],'Float64'))
          S2=myradians(N.array([148.389],'Float64'))
          H,K,L,E,Q,Ei,Ef=self.fixture.SpecWhere(M2,S1,S2,A2,[self.fixture.EXP])
          print 'H ',H
          print 'K ',K
          print 'L ',L
          print 'E ',E
          print 'Q ',Q
          print 'Ei ',Ei
          print 'Ef ',Ef
          self.assertAlmostEqual(H[0],1.6289,4, 'H Not equal to '+ str(1.6289))
          self.assertAlmostEqual(K[0],1.6289,4, 'K Not equal to '+ str(1.6289))
          self.assertAlmostEqual(L[0],2.1389,4, 'L Not equal to '+ str(2.1389))
          self.assertAlmostEqual(E[0],0.0000,4, 'E Not equal to '+ str(0.0000))
          self.assertAlmostEqual(Ei[0],3.6997,4,'Ei Not equal to '+str(3.6997))
          self.assertAlmostEqual(Ef[0],3.6997,4,'Ef Not equal to '+str(3.6997))   

     def test_tetragonal4(self):
          """test the tetragonal cell, change Energy transfer, compare with tetragonal3"""
          #setup lattice
          self.fixture.a=N.array([6.283],'Float64')
          self.fixture.b=N.array([6.283],'Float64')
          self.fixture.c=N.array([11.765],'Float64')
          orient1=N.array([[1,1,0]],'Float64')
          orient2=N.array([[0,0,1]],'Float64')
          orientation=Orientation(orient1,orient2)
          self.fixture.orientation=orientation
          
          #setup spectrometer
          self.fixture.EXP['infix']=-1 # Fixed Ei
          self.fixture.EXP['efixed']=3.7
          
          #test the angles
          M2=myradians(N.array([89.008],'Float64'))
          A2=myradians(N.array([98.663],'Float64'))
          S1=myradians(N.array([91.561],'Float64'))
          S2=myradians(N.array([117.979],'Float64'))
          H,K,L,E,Q,Ei,Ef=self.fixture.SpecWhere(M2,S1,S2,A2,[self.fixture.EXP])
          print 'H ',H
          print 'K ',K
          print 'L ',L
          print 'E ',E
          print 'Q ',Q
          print 'Ei ',Ei
          print 'Ef ',Ef
          self.assertAlmostEqual(H[0],1.3329,4, 'H Not equal to '+ str(1.3329))
          self.assertAlmostEqual(K[0],1.3329,4, 'K Not equal to '+ str(1.3329))
          self.assertAlmostEqual(L[0],2.1389,4, 'L Not equal to '+ str(2.1389))
          self.assertAlmostEqual(E[0],0.5400,4, 'E Not equal to '+ str(0.5400))
          self.assertAlmostEqual(Ei[0],3.6997,4,'Ei Not equal to '+str(3.6997))
          self.assertAlmostEqual(Ef[0],3.1597,4,'Ef Not equal to '+str(3.1597))   

     def test_tetragonal5(self):
          """test the tetragonal cell, swap orientation vectors"""
          #setup lattice
          self.fixture.a=N.array([6.283],'Float64')
          self.fixture.b=N.array([6.283],'Float64')
          self.fixture.c=N.array([11.765],'Float64')
          orient1=N.array([[0,0,1]],'Float64')
          orient2=N.array([[1,1,0]],'Float64')
          orientation=Orientation(orient1,orient2)
          self.fixture.orientation=orientation
          
          #setup spectrometer
          self.fixture.EXP['infix']=-1 # Fixed Ei
          self.fixture.EXP['efixed']=3.7
          
          #test the angles
          M2=myradians(N.array([89.008],'Float64'))
          A2=myradians(N.array([98.663],'Float64'))
          S1=myradians(N.array([119.133],'Float64'))  #recall how icp chooses the a3 angle based on first orientation vector
          S2=myradians(N.array([117.979],'Float64'))
          H,K,L,E,Q,Ei,Ef=self.fixture.SpecWhere(M2,S1,S2,A2,[self.fixture.EXP])
          print 'H ',H
          print 'K ',K
          print 'L ',L
          print 'E ',E
          print 'Q ',Q
          print 'Ei ',Ei
          print 'Ef ',Ef
          self.assertAlmostEqual(H[0],1.3329,4, 'H Not equal to '+ str(1.3329))
          self.assertAlmostEqual(K[0],1.3329,4, 'K Not equal to '+ str(1.3329))
          self.assertAlmostEqual(L[0],2.1389,4, 'L Not equal to '+ str(2.1389))
          self.assertAlmostEqual(E[0],0.5400,4, 'E Not equal to '+ str(0.5400))
          self.assertAlmostEqual(Ei[0],3.6997,4,'Ei Not equal to '+str(3.6997))
          self.assertAlmostEqual(Ef[0],3.1597,4,'Ef Not equal to '+str(3.1597))   

     def test_tetragonal6(self):
          """test the tetragonal cell, swap orientation vectors compared with tetragonal 4"""
          #setup lattice
          self.fixture.a=N.array([6.283],'Float64')
          self.fixture.b=N.array([6.283],'Float64')
          self.fixture.c=N.array([11.765],'Float64')
          orient1=N.array([[0,0,1]],'Float64')
          orient2=N.array([[1,1,0]],'Float64')
          orientation=Orientation(orient1,orient2)
          self.fixture.orientation=orientation
          
          #setup spectrometer
          self.fixture.EXP['infix']=-1 # Fixed Ei
          self.fixture.EXP['efixed']=3.7
          
          #test the angles
          M2=myradians(N.array([89.008],'Float64'))
          A2=myradians(N.array([89.008],'Float64'))
          S1=myradians(N.array([137.820],'Float64'))  #recall how icp chooses the a3 angle based on first orientation vector
          S2=myradians(N.array([148.389],'Float64'))
          H,K,L,E,Q,Ei,Ef=self.fixture.SpecWhere(M2,S1,S2,A2,[self.fixture.EXP])
          print 'H ',H
          print 'K ',K
          print 'L ',L
          print 'E ',E
          print 'Q ',Q
          print 'Ei ',Ei
          print 'Ef ',Ef
          self.assertAlmostEqual(H[0],1.6289,4, 'H Not equal to '+ str(1.6289))
          self.assertAlmostEqual(K[0],1.6289,4, 'K Not equal to '+ str(1.6289))
          self.assertAlmostEqual(L[0],2.1389,4, 'L Not equal to '+ str(2.1389))
          self.assertAlmostEqual(E[0],0.0000,4, 'E Not equal to '+ str(0.0000))
          self.assertAlmostEqual(Ei[0],3.6997,4,'Ei Not equal to '+str(3.6997))
          self.assertAlmostEqual(Ef[0],3.6997,4,'Ef Not equal to '+str(3.6997))   


     def test_orthorhombic1(self):
          """test the orthorhombic cell"""
          #setup lattice
          self.fixture.a=N.array([6.283],'Float64')
          self.fixture.b=N.array([5.7568],'Float64')
          self.fixture.c=N.array([11.765],'Float64')
          orient1=N.array([[1,0,0]],'Float64')
          orient2=N.array([[0,0,2]],'Float64')
          orientation=Orientation(orient1,orient2)
          self.fixture.orientation=orientation
          
          #setup spectrometer
          self.fixture.EXP['infix']=1 # Fixed Ef
          self.fixture.EXP['efixed']=4.5
          
          #test the angles
          M2=myradians(N.array([78.930],'Float64'))
          A2=myradians(N.array([78.930],'Float64'))
          S1=myradians(N.array([89.644],'Float64'))  
          S2=myradians(N.array([86.470],'Float64'))
          H,K,L,E,Q,Ei,Ef=self.fixture.SpecWhere(M2,S1,S2,A2,[self.fixture.EXP])
          print 'H ',H
          print 'K ',K
          print 'L ',L
          print 'E ',E
          print 'Q ',Q
          print 'Ei ',Ei
          print 'Ef ',Ef
          self.assertAlmostEqual(H[0],1.3919,4, 'H Not equal to '+ str(1.3919))
          self.assertAlmostEqual(K[0],0.0000,4, 'K Not equal to '+ str(0.0000))
          self.assertAlmostEqual(L[0],2.7379,4, 'L Not equal to '+ str(2.7379))
          self.assertAlmostEqual(E[0],0.0000,4, 'E Not equal to '+ str(0.0000))
          self.assertAlmostEqual(Ei[0],4.4996,4,'Ei Not equal to '+str(4.4996))
          self.assertAlmostEqual(Ef[0],4.4996,4,'Ef Not equal to '+str(4.4996))   
   
     def test_orthorhombic2(self):
          """test the orthorhombic cell, change E transfer"""
          #setup lattice
          self.fixture.a=N.array([6.283],'Float64')
          self.fixture.b=N.array([5.7568],'Float64')
          self.fixture.c=N.array([11.765],'Float64')
          orient1=N.array([[1,0,0]],'Float64')
          orient2=N.array([[0,0,2]],'Float64')
          orientation=Orientation(orient1,orient2)
          self.fixture.orientation=orientation
          
          #setup spectrometer
          self.fixture.EXP['infix']=1 # Fixed Ef
          self.fixture.EXP['efixed']=4.5
          
          #test the angles
          M2=myradians(N.array([46.305],'Float64'))
          A2=myradians(N.array([78.930],'Float64'))
          S1=myradians(N.array([98.405],'Float64'))  
          S2=myradians(N.array([57.515],'Float64'))
          H,K,L,E,Q,Ei,Ef=self.fixture.SpecWhere(M2,S1,S2,A2,[self.fixture.EXP])
          print 'H ',H
          print 'K ',K
          print 'L ',L
          print 'E ',E
          print 'Q ',Q
          print 'Ei ',Ei
          print 'Ef ',Ef
          self.assertAlmostEqual(H[0],1.3919,4, 'H Not equal to '+ str(1.3919))
          self.assertAlmostEqual(K[0],0.0000,4, 'K Not equal to '+ str(0.0000))
          self.assertAlmostEqual(L[0],2.7379,4, 'L Not equal to '+ str(2.7379))
          self.assertAlmostEqual(E[0],7.2594,4, 'E Not equal to '+ str(7.2594))
          self.assertAlmostEqual(Ei[0],11.7590,4,'Ei Not equal to '+str(11.7590))
          self.assertAlmostEqual(Ef[0],4.4996,4,'Ef Not equal to '+str(4.4996))   
          
     def test_orthorhombic3(self):
          """test the orthorhombic cell, change orientation vectors, compare with ortho2"""
          #setup lattice
          self.fixture.a=N.array([6.283],'Float64')
          self.fixture.b=N.array([5.7568],'Float64')
          self.fixture.c=N.array([11.765],'Float64')
          orient1=N.array([[1,1,0]],'Float64')
          orient2=N.array([[0,0,2]],'Float64')
          orientation=Orientation(orient1,orient2)
          self.fixture.orientation=orientation
          
          #setup spectrometer
          self.fixture.EXP['infix']=1 # Fixed Ef
          self.fixture.EXP['efixed']=4.5
          
          #test the angles
          M2=myradians(N.array([78.930],'Float64'))
          A2=myradians(N.array([78.930],'Float64'))
          S1=myradians(N.array([91.228],'Float64'))  
          S2=myradians(N.array([105.102],'Float64'))
          H,K,L,E,Q,Ei,Ef=self.fixture.SpecWhere(M2,S1,S2,A2,[self.fixture.EXP])
          print 'H ',H
          print 'K ',K
          print 'L ',L
          print 'E ',E
          print 'Q ',Q
          print 'Ei ',Ei
          print 'Ef ',Ef
          self.assertAlmostEqual(H[0],1.2339,4, 'H Not equal to '+ str(1.2339))
          self.assertAlmostEqual(K[0],1.2339,4, 'K Not equal to '+ str(1.2339))
          self.assertAlmostEqual(L[0],2.7379,4, 'L Not equal to '+ str(2.7379))
          self.assertAlmostEqual(E[0],0.0000,4, 'E Not equal to '+ str(0.0000))
          self.assertAlmostEqual(Ei[0],4.4996,4,'Ei Not equal to '+str(4.4966))
          self.assertAlmostEqual(Ef[0],4.4996,4,'Ef Not equal to '+str(4.4996))   

     def test_orthorhombic4(self):
          """test the orthorhombic cell, change E, compare with ortho3"""
          #setup lattice
          self.fixture.a=N.array([6.283],'Float64')
          self.fixture.b=N.array([5.7568],'Float64')
          self.fixture.c=N.array([11.765],'Float64')
          orient1=N.array([[1,1,0]],'Float64')
          orient2=N.array([[0,0,2]],'Float64')
          orientation=Orientation(orient1,orient2)
          self.fixture.orientation=orientation
          
          #setup spectrometer
          self.fixture.EXP['infix']=1 # Fixed Ef
          self.fixture.EXP['efixed']=4.5
          
          #test the angles
          M2=myradians(N.array([47.030],'Float64'))
          A2=myradians(N.array([78.930],'Float64'))
          S1=myradians(N.array([92.030],'Float64'))  
          S2=myradians(N.array([71.391],'Float64'))
          H,K,L,E,Q,Ei,Ef=self.fixture.SpecWhere(M2,S1,S2,A2,[self.fixture.EXP])
          print 'H ',H
          print 'K ',K
          print 'L ',L
          print 'E ',E
          print 'Q ',Q
          print 'Ei ',Ei
          print 'Ef ',Ef
          self.assertAlmostEqual(H[0],1.2339,4, 'H Not equal to '+ str(1.2339))
          self.assertAlmostEqual(K[0],1.2339,4, 'K Not equal to '+ str(1.2339))
          self.assertAlmostEqual(L[0],2.7379,4, 'L Not equal to '+ str(2.7379))
          self.assertAlmostEqual(E[0],6.9194,4, 'E Not equal to '+ str(6.9194))
          self.assertAlmostEqual(Ei[0],11.4190,4,'Ei Not equal to '+str(11.4190))
          self.assertAlmostEqual(Ef[0],4.4996,4,'Ef Not equal to '+str(4.4996))

          
     def test_orthorhombic5(self):
          """test the orthorhombic cell, swap orientation vectors, compare with ortho4"""
          #setup lattice
          self.fixture.a=N.array([6.283],'Float64')
          self.fixture.b=N.array([5.7568],'Float64')
          self.fixture.c=N.array([11.765],'Float64')
          orient1=N.array([[0,0,2]],'Float64')
          orient2=N.array([[1,1,0]],'Float64')
          orientation=Orientation(orient1,orient2)
          self.fixture.orientation=orientation
          
          #setup spectrometer
          self.fixture.EXP['infix']=1 # Fixed Ef
          self.fixture.EXP['efixed']=4.5
          
          #test the angles
          M2=myradians(N.array([47.030],'Float64'))
          A2=myradians(N.array([78.930],'Float64'))
          S1=myradians(N.array([104.676],'Float64'))  #recalll how icp defines a3 in terms of a4 of orient1  
          S2=myradians(N.array([71.391],'Float64'))
          H,K,L,E,Q,Ei,Ef=self.fixture.SpecWhere(M2,S1,S2,A2,[self.fixture.EXP])
          print 'H ',H
          print 'K ',K
          print 'L ',L
          print 'E ',E
          print 'Q ',Q
          print 'Ei ',Ei
          print 'Ef ',Ef
          self.assertAlmostEqual(H[0],1.2339,4, 'H Not equal to '+ str(1.2339))
          self.assertAlmostEqual(K[0],1.2339,4, 'K Not equal to '+ str(1.2339))
          self.assertAlmostEqual(L[0],2.7379,4, 'L Not equal to '+ str(2.7379))
          self.assertAlmostEqual(E[0],6.9194,4, 'E Not equal to '+ str(6.9194))
          self.assertAlmostEqual(Ei[0],11.4190,4,'Ei Not equal to '+str(11.4190))
          self.assertAlmostEqual(Ef[0],4.4996,4,'Ef Not equal to '+str(4.4996))

     def test_orthorhombic6(self):
          """test the orthorhombic cell, swap orient vectors, compare with ortho3"""
          #setup lattice
          self.fixture.a=N.array([6.283],'Float64')
          self.fixture.b=N.array([5.7568],'Float64')
          self.fixture.c=N.array([11.765],'Float64')
          orient1=N.array([[0,0,2]],'Float64')
          orient2=N.array([[1,1,0]],'Float64')
          orientation=Orientation(orient1,orient2)
          self.fixture.orientation=orientation
          
          #setup spectrometer
          self.fixture.EXP['infix']=1 # Fixed Ef
          self.fixture.EXP['efixed']=4.5
          
          #test the angles
          M2=myradians(N.array([78.930],'Float64'))
          A2=myradians(N.array([78.930],'Float64'))
          S1=myradians(N.array([103.874],'Float64'))  #recalll how icp defines a3 in terms of a4 of orient1  
          S2=myradians(N.array([105.102],'Float64'))
          H,K,L,E,Q,Ei,Ef=self.fixture.SpecWhere(M2,S1,S2,A2,[self.fixture.EXP])
          print 'H ',H
          print 'K ',K
          print 'L ',L
          print 'E ',E
          print 'Q ',Q
          print 'Ei ',Ei
          print 'Ef ',Ef
          self.assertAlmostEqual(H[0],1.2339,4, 'H Not equal to '+ str(1.2339))
          self.assertAlmostEqual(K[0],1.2339,4, 'K Not equal to '+ str(1.2339))
          self.assertAlmostEqual(L[0],2.7379,4, 'L Not equal to '+ str(2.7379))
          self.assertAlmostEqual(E[0],0.0000,4, 'E Not equal to '+ str(0.0000))
          self.assertAlmostEqual(Ei[0],4.4996,4,'Ei Not equal to '+str(4.4996))
          self.assertAlmostEqual(Ef[0],4.4996,4,'Ef Not equal to '+str(4.4996))     
          
          
     def test_monoclinic1(self):
          """test monoclinic 1"""
          #setup lattice
          self.fixture.a=N.array([6.283],'Float64')
          self.fixture.b=N.array([5.7568],'Float64')
          self.fixture.c=N.array([11.765],'Float64')
          self.fixture.beta=myradians(N.array([100.0],'Float64'))
          orient1=N.array([[1,0,0]],'Float64')
          orient2=N.array([[0,0,1]],'Float64')
          orientation=Orientation(orient1,orient2)
          self.fixture.orientation=orientation
          
          #setup spectrometer
          self.fixture.EXP['infix']=1 # Fixed Ef
          self.fixture.EXP['efixed']=4.5
          
          #test the angles
          M2=myradians(N.array([78.930],'Float64'))
          A2=myradians(N.array([78.930],'Float64'))
          S1=myradians(N.array([108.743],'Float64'))  #recalll how icp defines a3 in terms of a4 of orient1  
          S2=myradians(N.array([130.130],'Float64'))
          H,K,L,E,Q,Ei,Ef=self.fixture.SpecWhere(M2,S1,S2,A2,[self.fixture.EXP])
          print 'H ',H
          print 'K ',K
          print 'L ',L
          print 'E ',E
          print 'Q ',Q
          print 'Ei ',Ei
          print 'Ef ',Ef
          self.assertAlmostEqual(H[0],1.5829,4, 'H Not equal to '+ str(1.5829))
          self.assertAlmostEqual(K[0],0.0000,4, 'K Not equal to '+ str(0.0000))
          self.assertAlmostEqual(L[0],3.4558,4, 'L Not equal to '+ str(3.4558))
          self.assertAlmostEqual(E[0],0.0000,4, 'E Not equal to '+ str(0.0000))
          self.assertAlmostEqual(Ei[0],4.4996,4,'Ei Not equal to '+str(4.4996))
          self.assertAlmostEqual(Ef[0],4.4996,4,'Ef Not equal to '+str(4.4996))
 
     def test_monoclinic2(self):
          """test monoclinic 2, change E"""
          #setup lattice
          self.fixture.a=N.array([6.283],'Float64')
          self.fixture.b=N.array([5.7568],'Float64')
          self.fixture.c=N.array([11.765],'Float64')
          self.fixture.beta=myradians(N.array([100.0],'Float64'))
          orient1=N.array([[1,0,0]],'Float64')
          orient2=N.array([[0,0,1]],'Float64')
          orientation=Orientation(orient1,orient2)
          self.fixture.orientation=orientation
          
          #setup spectrometer
          self.fixture.EXP['infix']=1 # Fixed Ef
          self.fixture.EXP['efixed']=4.5
          
          #test the angles
          M2=myradians(N.array([74.186],'Float64'))
          A2=myradians(N.array([78.930],'Float64'))
          S1=myradians(N.array([106.473],'Float64'))  #recalll how icp defines a3 in terms of a4 of orient1  
          S2=myradians(N.array([123.991],'Float64'))
          H,K,L,E,Q,Ei,Ef=self.fixture.SpecWhere(M2,S1,S2,A2,[self.fixture.EXP])
          print 'H ',H
          print 'K ',K
          print 'L ',L
          print 'E ',E
          print 'Q ',Q
          print 'Ei ',Ei
          print 'Ef ',Ef
          self.assertAlmostEqual(H[0],1.5829,4, 'H Not equal to '+ str(1.5829))
          self.assertAlmostEqual(K[0],0.0000,4, 'K Not equal to '+ str(0.0000))
          self.assertAlmostEqual(L[0],3.4558,4, 'L Not equal to '+ str(3.4558))
          self.assertAlmostEqual(E[0],0.4979,4, 'E Not equal to '+ str(0.4979))
          self.assertAlmostEqual(Ei[0],4.9976,4,'Ei Not equal to '+str(4.9976))
          self.assertAlmostEqual(Ef[0],4.4996,4,'Ef Not equal to '+str(4.4996))

          
     def test_monoclinic3(self):
          """test monoclinic 2, change orient vectors"""
          #setup lattice
          self.fixture.a=N.array([6.283],'Float64')
          self.fixture.b=N.array([5.7568],'Float64')
          self.fixture.c=N.array([11.765],'Float64')
          self.fixture.beta=myradians(N.array([100.0],'Float64'))
          orient1=N.array([[1,2,0]],'Float64')
          orient2=N.array([[0,0,1]],'Float64')
          orientation=Orientation(orient1,orient2)
          self.fixture.orientation=orientation
          
          #setup spectrometer
          self.fixture.EXP['infix']=1 # Fixed Ef
          self.fixture.EXP['efixed']=4.5
          
          #test the angles
          M2=myradians(N.array([78.930],'Float64'))
          A2=myradians(N.array([78.930],'Float64'))
          S1=myradians(N.array([70.753],'Float64'))  #recalll how icp defines a3 in terms of a4 of orient1  
          S2=myradians(N.array([91.256],'Float64'))
          H,K,L,E,Q,Ei,Ef=self.fixture.SpecWhere(M2,S1,S2,A2,[self.fixture.EXP])
          print 'H ',H
          print 'K ',K
          print 'L ',L
          print 'E ',E
          print 'Q ',Q
          print 'Ei ',Ei
          print 'Ef ',Ef
          self.assertAlmostEqual(H[0],0.7650,4, 'H Not equal to '+ str(0.7650))
          self.assertAlmostEqual(K[0],1.5299,4, 'K Not equal to '+ str(1.5299))
          self.assertAlmostEqual(L[0],1.6539,4, 'L Not equal to '+ str(1.6539))
          self.assertAlmostEqual(E[0],0.0000,4, 'E Not equal to '+ str(0.0000))
          self.assertAlmostEqual(Ei[0],4.4996,4,'Ei Not equal to '+str(4.49996))
          self.assertAlmostEqual(Ef[0],4.4996,4,'Ef Not equal to '+str(4.4996))


     def test_monoclinic4(self):
          """test monoclinic 2, change E"""
          #setup lattice
          self.fixture.a=N.array([6.283],'Float64')
          self.fixture.b=N.array([5.7568],'Float64')
          self.fixture.c=N.array([11.765],'Float64')
          self.fixture.beta=myradians(N.array([100.0],'Float64'))
          orient1=N.array([[1,2,0]],'Float64')
          orient2=N.array([[0,0,1]],'Float64')
          orientation=Orientation(orient1,orient2)
          self.fixture.orientation=orientation
          
          #setup spectrometer
          self.fixture.EXP['infix']=1 # Fixed Ef
          self.fixture.EXP['efixed']=4.5
          
          #test the angles
          M2=myradians(N.array([56.239],'Float64'))
          A2=myradians(N.array([78.930],'Float64'))
          S1=myradians(N.array([73.059],'Float64'))  #recalll how icp defines a3 in terms of a4 of orient1  
          S2=myradians(N.array([73.305],'Float64'))
          H,K,L,E,Q,Ei,Ef=self.fixture.SpecWhere(M2,S1,S2,A2,[self.fixture.EXP])
          print 'H ',H
          print 'K ',K
          print 'L ',L
          print 'E ',E
          print 'Q ',Q
          print 'Ei ',Ei
          print 'Ef ',Ef
          self.assertAlmostEqual(H[0],0.7650,4, 'H Not equal to '+ str(0.7650))
          self.assertAlmostEqual(K[0],1.5299,4, 'K Not equal to '+ str(1.5299))
          self.assertAlmostEqual(L[0],1.6539,4, 'L Not equal to '+ str(1.6539))
          self.assertAlmostEqual(E[0],3.6838,4, 'E Not equal to '+ str(3.6838))
          self.assertAlmostEqual(Ei[0],8.1834,4,'Ei Not equal to '+str(8.1834))
          self.assertAlmostEqual(Ef[0],4.4996,4,'Ef Not equal to '+str(4.4996))

          
if __name__=="__main__":
     if 0:
          a=N.array([2*pi,2*pi],'Float64')
          b=N.array([8],'Float64')
          c=N.array([11],'Float64')
          alpha=myradians(N.array([87],'Float64'))
          beta=myradians(N.array([52],'Float64'))
          gamma=myradians(N.array([100],'Float64'))
          orient1=N.array([[0,1,0]],'Float64')
          orient2=N.array([[1,0,0]],'Float64')
          orientation=Orientation(orient1,orient2)
          self.fixture.orientation=orientation
          mylattice=Lattice(a=a,b=b,c=c,alpha=alpha,beta=beta,gamma=gamma,\
                            orientation=orientation)
          H=N.array([1],'Float64');K=N.array([0],'Float64');L=N.array([0],'Float64');W=N.array([0],'Float64')
          EXP={}
          EXP['ana']={}
          EXP['ana']['tau']='pg(002)'
          EXP['mono']={}
          EXP['mono']['tau']='pg(002)';
          EXP['ana']['mosaic']=30
          EXP['mono']['mosaic']=30
          EXP['sample']={}
          EXP['sample']['mosaic']=10
          EXP['sample']['vmosaic']=10
          EXP['hcol']=N.array([40, 10, 20, 80],'Float64')
          EXP['vcol']=N.array([120, 120, 120, 120],'Float64')
          EXP['infix']=-1 #positive for fixed incident energy
          EXP['efixed']=14.7
          EXP['method']=0
          setup=[EXP]
          M2=myradians(N.array([41.177]))
          A2=myradians(N.array([41.177]))
          S1=myradians(N.array([66.4363]))
          S2=myradians(N.array([37.6547]))
          H,K,L,E,Q,Ei,Ef=mylattice.SpecWhere(M2,S1,S2,A2,setup)
          print 'H ',H
          print 'K ',K
          print 'L ',L
          print 'E ',E
          print 'Q ',Q
          print 'Ei ',Ei
          print 'Ef ',Ef
          M1,M2,S1,S2,A1,A2=mylattice.SpecGoTo(H,K,L,E,setup)
          print 'M2 ',mydegrees(M2)
          print 'A2 ',mydegrees(A2)
          print 'M1 ',mydegrees(M1)
          print 'A1 ',mydegrees(A1)
          print 'S1 ',mydegrees(S1)
          print 'S2 ',mydegrees(S2)
     if 1:
          unittest.main()