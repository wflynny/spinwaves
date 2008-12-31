import sympy
import numpy as N
#I=sympy.I
#I=1.0j
#pi=sympy.pi
from sympy import exp,I,pi,sin,cos
import pylab
import readfiles
from sympy import pngview,latex
import scipy.linalg

#translations=[[0,0,0],
#              [0,0,1],[0,0,-1]
#              [0,1,0],[0,-1,0]
#              [0,1,1],[0,1,-1],[0,-1,1],[0,-1,-1]
#              [1,0,0],[-1,0,0],[]              
#              ]

def print_matplotlib(s):
    pylab.text(0,0,s)
    pylab.axis('off')
    pylab.show()
    return 


def coeff(expr, term):
   expr = sympy.collect(expr, term)
   #print 'expr',expr
   symbols = list(term.atoms(sympy.Symbol))
   #print 'symbols',symbols
   w = sympy.Wild("coeff", exclude=symbols)
   #print 'w',w
   m = expr.match(w*term+sympy.Wild("rest"))
   #print 'm',m
   m2=expr.match(w*term)
   #print 'm2',m2
   res=False
   if m2!=None:
       #print 'm2[w]',m2[w]
       res=m2[w]*term==expr
   if m and res!=True:
       return m[w]
   #added the next two lines
   elif m2:
       return m2[w]
       
   
class atom:
    def __init__(self,spin=[0,0,1],pos=[0,0,0],neighbors=None,interactions=None,label=0,Dx=0,Dy=0,Dz=0,cell=0,int_cell=[]):
        self.spin=spin
        if neighbors==None:
            neighbors=[]
        if interactions==None:
            interactions=[]
        self.pos=N.array(pos)
        self.neighbors=neighbors
        self.interactions=interactions
        self.label=label
        self.Dx=Dx
        self.Dy=Dy
        self.Dz=Dz
        self.cell=cell
        self.int_cell=[]


def generate_atoms():



    if 1:
        D=sympy.Symbol('D',real=True)
        spin0=N.matrix([[1,0,0],[0,1,0],[0,0,1]],'float64')
        pos0=[0,0,0]
        neighbors=[1]
        interactions=[0]
        cell=0
        int_cell=[5,21]
        atom0=atom(spin=spin0,pos=pos0,neighbors=neighbors,interactions=interactions,label=0,cell=cell,int_cell=int_cell,Dz=0)
        
        pos0=[1,0,0]
        spin0=N.matrix([[-1,0,0],[0,1,0],[0,0,1]],'float64')
        spin0=N.matrix([[1,0,0],[0,1,0],[0,0,1]],'float64')
        neighbors=[0]
        interactions=[0]
        cell=5
        int_cell=[0]
        atom1=atom(spin=spin0,pos=pos0,neighbors=neighbors,interactions=interactions,label=1,cell=cell,int_cell=int_cell,Dz=D)
        
        atomlist=[atom0,atom1]


 
       
    return atomlist


def generate_atoms_rot():



    if 1:
        D=sympy.Symbol('D',real=True)
        spin0=sympy.matrices.Matrix([[1,0,0],[0,1,0],[0,0,1]])
        #spin0=sympy.matrices.Matrix(spin0)
        pos0=[0,0,0]
        neighbors=[1]
        interactions=[0]
        cell=0
        int_cell=[5,21]
        atom0=atom(spin=spin0,pos=pos0,neighbors=neighbors,interactions=interactions,label=0,cell=cell,int_cell=int_cell,Dz=0)
        
        pos0=[1,0,0]
        #spin0=N.matrix([[-1,0,0],[0,1,0],[0,0,1]],'float64')
        #spin0=N.matrix([[-1,0,0],[0,1,0],[0,0,1]],'float64')
        P=sympy.Symbol('P',real=True,commutative=True)
        W=sympy.Symbol('W',real=True,commutative=True)
        X=sympy.Symbol('X',real=True,commutative=True)
        #spin0=sympy.matrices.Matrix([[W,-X,0],[W,X,0],[0,0,1]])
        #P=0
        spin0=sympy.matrices.Matrix([[cos(P),-sin(P),0],[sin(P),cos(P),0],[0,0,1]])
        spin0=sympy.matrices.Matrix([[1,0,0],[0,1,0],[0,0,1]])
        print 'spin0 analaytic', spin0
        #spin0=sympy.matrices.Matrix([[1,0,0],[0,1,0],[0,0,1]])
        #spin0=spin0.subs(P,sympy.pi)
        print 'spin0 converted',spin0
        neighbors=[0]
        interactions=[0]
        cell=5
        int_cell=[0]
        atom1=atom(spin=spin0,pos=pos0,neighbors=neighbors,interactions=interactions,label=1,cell=cell,int_cell=int_cell,Dz=0)
        
        atomlist=[atom0,atom1]


 
       
    return atomlist




def generate_sabn(N_atoms):
    "generate spins in local coordinate system, with Z as quantization axis"
    Sabn=[]
    S=sympy.Symbol("S",real=True)
    for i in range(N_atoms):
        c=sympy.Symbol('c%d'%(i,),commutative=False,real=True)
        cd=sympy.Symbol('cd%d'%(i,),commutative=False,real=True)
        curr=sympy.matrices.Matrix([sympy.sqrt(S/2.0)*(c+cd),sympy.sqrt(S/2.0)*(c-cd)/I,S-cd*c])
        Sabn.append(curr.reshape(3,1))
    return Sabn


def generate_sxyz(Sabn,atomlist):
    "transform spins from local coordinate system to global system"
    Sxyz=[]
    i=0
    for currS in Sabn:
        tempS=atomlist[i].spin*currS
        tempS=tempS.reshape(1,3)
        Sxyz.append(tempS)
        i=i+1
    return Sxyz




def generate_translations():
    translations=[[0,0,0]]
    for i in range(-1,2):
        for j in range(-1,2):
            for k in range(-1,2):
                if i==0 and j==0 and k==0:
                    continue
                else:
                    translations.append([i,j,k])
    return translations

def generate_sabnt(N_atoms,t=''):
    Sabn=[]
    S=sympy.Symbol("S",real=True)
    for i in range(N_atoms):
        ci=sympy.Symbol("c%d%s"%(i,t),commutative=False)
        cdi=sympy.Symbol("cd%d%s"%(i,t),commutative=False)
        curr=[sympy.sqrt(S/2)*(ci+cdi),sympy.sqrt(S/2)*(ci-cdi)/I,S-cdi*ci]
        Sabn.append(curr)
    return Sabn





def generate_hdef(atom_list,Jij,Sxyz,N_atoms_uc,N_atoms):
    "generate the hamiltonian for a set of interacting spins defined by Sxyz"
    N_atoms=len(atom_list)
    Hdef=0
    print 'Jij',Jij,len(Jij)

    
    
    #for i in range(N_atoms):
    for i in range(N_atoms_uc): #correct
        N_int=len(atom_list[i].interactions)
        for j in range(N_int):            
            print 'i',i,'j',j
            if 0:
                Hij=N.matrix(Sxyz[i])*atom_list[i].spin.T
                print 'making Ham'
                print 'spin i', atom_list[i].spin.T
                print 'Sxyz i', N.matrix(Sxyz[i]),N.matrix(Sxyz[i]).shape
                print 'Hijtemp',Hij
                Hij=Hij*Jij[atom_list[i].interactions[j]]
                print 'Jij', Jij[atom_list[i].interactions[j]]
                print 'Hij*Jij', Hij
                Hij=Hij*atom_list[atom_list[i].neighbors[j]].spin#
                print 'Hijtemp3', Hij.shape
            if 1:
                print 'Sxyz i', N.matrix(Sxyz[i]),N.matrix(Sxyz[i]).shape
                Hij=Sxyz[i]*Jij[atom_list[i].interactions[j]]
                print 'S*Jij', Hij,Hij.shape
            Sxyz_transpose=Sxyz[atom_list[i].neighbors[j]].reshape(3,1)
            #Sxyz_transpose=Sxyz_transpose.(3,1))
            print 'Sxyz.T',Sxyz_transpose.shape
            print 'Hij before multiply', Hij, Hij.shape
            Hij=Hij*Sxyz_transpose
            print 'Hij*Sxyz.T',Hij,Hij.shape
            Hij=Hij[0]
            #Hij=Hij+myterm
            Hij=-Hij-atom_list[i].Dx*Sxyz[i][0]**2-atom_list[i].Dy*Sxyz[i][1]**2-atom_list[i].Dz*Sxyz[i][2]**2
            Hdef=Hdef+Hij
    print 'generated hdef'
    print Hdef,Hdef.atoms(sympy.Symbol)
    return Hdef


def holstein(Hdef):
        S = sympy.Symbol('S',real=True)
        print 'holstein'
        print Hdef.atoms(sympy.Symbol)
        Hdef=Hdef.expand()
        #Hdef=Hdef.as_poly(S)
        p = sympy.Wild('p',exclude='S')
        q = sympy.Wild('q',exclude='S')
        r = sympy.Wild('r',exclude='S')
        l = sympy.Wild('l',exclude='S')
        #Hlin=Hdef.coeffs[0]*S**2+Hdef.coeffs[1]*S
        S2coeff=coeff(Hdef,S**2)
        Scoeff=coeff(Hdef,S)
        Hlin=None
        #Hlin=coeff(Hdef,S**2)*S**2+coeff(Hdef,S)*S
        #print 'S2Coeff', S2coeff
        #print 'Scoeff',Scoeff
        if Scoeff!=None and S2coeff!=None:
            Hlin=coeff(Hdef,S**2)*S**2+coeff(Hdef,S)*S
        elif Scoeff==None and S2coeff!=None:
            Hlin=coeff(Hdef,S**2)*S**2
        elif Scoeff!=None and S2coeff==None:
            #print 'S'
            Hlin=coeff(Hdef,S)*S
        return Hlin


def fouriertransform(atom_list,Jij,Hlin,k,N_atoms_uc,N_atoms):
    #N_atoms=len(atom_list)
    #N_atoms_uc=1
    #N_atoms_uc=N_atoms
    #Hdef=0
    #print 'atom_list',atom_list
    #print 'Sxyz',Sxyz
    #print 'Jij',Jij
    print 'fourier'
    print Hlin
    print Hlin.atoms(sympy.Symbol)
    print 'expand'
    #Hlin=Hlin.expand()
    print Hlin.atoms(sympy.Symbol)
    #print Hlin
    #for i in range(N_atoms):
    for i in range(N_atoms_uc): #correct
        N_int=len(atom_list[i].interactions)
        ci=sympy.Symbol('c%d'%(i,),commutative=False,real=True)
        cdi=sympy.Symbol('cd%d'%(i,),commutative=False,real=True)
        cki=sympy.Symbol('ck%d'%(i,),commutative=False,real=True)
        ckdi=sympy.Symbol('ckd%d'%(i,),commutative=False,real=True)
        cmki=sympy.Symbol('cmk%d'%(i,),commutative=False,real=True)
        cmkdi=sympy.Symbol('cmkd%d'%(i,),commutative=False,real=True)
        ri=atom_list[i].pos
        N_int=len(atom_list[i].interactions)
        for j in range(N_atoms):
            rj=atom_list[j].pos
            j2=i#atom_list[i].neighbors[j]
            cj=sympy.Symbol('c%d'%(j,),commutative=False,real=True)
            cdj=sympy.Symbol('cd%d'%(j,),commutative=False,real=True)
            ckj=sympy.Symbol('ck%d'%(j2,),commutative=False,real=True)
            ckdj=sympy.Symbol('ckd%d'%(j2,),commutative=False,real=True)
            cmkj=sympy.Symbol('cmk%d'%(j2,),commutative=False,real=True)
            cmkdj=sympy.Symbol('cmkd%d'%(j2,),commutative=False,real=True)
            diffr=ri-rj
            kmult=N.dot(k,diffr)
            t1=1.0/2*(ckdi*cmkdj*exp(-I*kmult)+cmkdi*ckdj*exp(I*kmult)               )
            t2=1.0/2*(cki*cmkj*exp(I*kmult)+cmki*ckj*exp(-I*kmult))
            t3=1.0/2*(ckdi*ckj*exp(-I*kmult)+cmkdi*cmkj*exp(I*kmult))
            t4=1.0/2*(cki*ckdj*exp(I*kmult)+cmki*cmkdj*exp(-I*kmult))
            t5=1.0/2*(ckdj*ckj+cmkdj*cmkj)
            print 'i',i,'j',j
            #print 'ci',ci,'cj',cj,'cdi',cdi,'cdj',cdj
            #print 't1',t1
            #print 't2',t2
            #print 't3',t3
            #print 't4',t4
            #print 't5',t5
            #a1=sympy.Symbol('a1')
            #a2=sympy.Symbol('a2')
            #a3=sympy.Symbol('a3')
            #a4=sympy.Symbol('a4')
            #a5=sympy.Symbol('a5')
            f1=cdi*cdj
            print 'f1',f1,f1.atoms(sympy.Symbol)
            Hlin=Hlin.subs(f1,t1)
            print 'H1',Hlin,Hlin.atoms(sympy.Symbol)
            f2=ci*cj
            print 'f2',f2,f2.atoms(sympy.Symbol)
            Hlin=Hlin.subs(f2,t2)
            print 'H2',Hlin,Hlin.atoms(sympy.Symbol)
            f3=cdi*cj
            print 'f3',f3,f3.atoms(sympy.Symbol)
            Hlin=Hlin.subs(f3,t3)
            print 'H3',Hlin,Hlin.atoms(sympy.Symbol)
            f4=ci*cdj
            print 'f4',f4,f4.atoms(sympy.Symbol)
            Hlin=Hlin.subs(f4,t4)
            print 'H4',Hlin,Hlin.atoms(sympy.Symbol)
            f5=cdj*cj
            print 'f5',f5,f5.atoms(sympy.Symbol)
            Hlin=Hlin.subs(f5,t5)
            print 'H5',Hlin,Hlin.atoms(sympy.Symbol)
    #print t1
    return Hlin#.expand()


def applycommutation(atom_list,Jij,Hfou,k,N_atoms_uc,N_atoms):
    """Operate commutation relations to put all the 2nd order term as ckd**ck, cmk**cmkd, cmk**ck and ckd**cmkd form"""

    for i in range(N_atoms_uc):
        N_int=len(atom_list[i].interactions)
        ci=sympy.Symbol("c%d"%(i,),commutative=False,real=True)
        cdi=sympy.Symbol("cd%d"%(i,),commutative=False,real=True)
        cki=sympy.Symbol("ck%d"%(i,),commutative=False,real=True)
        ckdi=sympy.Symbol("ckd%d"%(i,),commutative=False,real=True)
        cmki=sympy.Symbol("cmk%d"%(i,),commutative=False,real=True)
        cmkdi=sympy.Symbol("cmkd%d"%(i,),commutative=False,real=True)
        for j in range(N_atoms_uc):
            cj=sympy.Symbol("c%d"%(j,),commutative=False,real=True)
            cdj=sympy.Symbol("cd%d"%(j,),commutative=False,real=True)
            ckj=sympy.Symbol("ck%d"%(j,),commutative=False,real=True)
            ckdj=sympy.Symbol("ckd%d"%(j,),commutative=False,real=True)
            cmkj=sympy.Symbol("cmk%d"%(j,),commutative=False,real=True)
            cmkdj=sympy.Symbol("cmkd%d"%(j,),commutative=False,real=True)
            if i==j:
                Hfou=Hfou.subs(cki*ckdj,ckdj*cki+1)
                
                Hfou=Hfou.subs(cmkdi*cmkj,cmkj*cmkdi+1)
                
            else:
                Hfou=Hfou.subs(cki*ckdj,ckdj*cki)
                Hfou=Hfou.subs(cmkdi*cmkj,cmkj*cmkdi)
            Hfou=Hfou.subs(cki*cmkj,cmkj*cki)
                
            Hfou=Hfou.subs(cmkdi*ckdj,ckdj*cmkdi)
    
    
    return Hfou.expand()

def gen_operator_table(atom_list,N_atoms_uc):
    """Operate commutation relations to put all the 2nd order term as ckd**ck, cmk**cmkd, cmk**ck and ckd**cmkd form"""

    operator_table=[]
    operator_table_minus=[]
    
    for i in range(N_atoms_uc):
        N_int=len(atom_list[i].interactions)
        ci=sympy.Symbol("c%d"%(i,),commutative=False,real=True)
        cdi=sympy.Symbol("cd%d"%(i,),commutative=False,real=True)
        cki=sympy.Symbol("ck%d"%(i,),commutative=False,real=True)
        ckdi=sympy.Symbol("ckd%d"%(i,),commutative=False,real=True)
        cmki=sympy.Symbol("cmk%d"%(i,),commutative=False,real=True)
        cmkdi=sympy.Symbol("cmkd%d"%(i,),commutative=False,real=True)
        operator_table.append(cki)
        operator_table_minus.append(cmkdi)
    operator_table=[operator_table,operator_table_minus]
    operator_table=N.ravel(operator_table)
    return operator_table


def gen_operator_table_dagger(atom_list,N_atoms_uc):
    """Operate commutation relations to put all the 2nd order term as ckd**ck, cmk**cmkd, cmk**ck and ckd**cmkd form"""

    operator_table=[]
    operator_table_minus=[]
    
    for i in range(N_atoms_uc):
        N_int=len(atom_list[i].interactions)
        ci=sympy.Symbol("c%d"%(i,),commutative=False,real=True)
        cdi=sympy.Symbol("cd%d"%(i,),commutative=False,real=True)
        cki=sympy.Symbol("ck%d"%(i,),commutative=False,real=True)
        ckdi=sympy.Symbol("ckd%d"%(i,),commutative=False,real=True)
        cmki=sympy.Symbol("cmk%d"%(i,),commutative=False,real=True)
        cmkdi=sympy.Symbol("cmkd%d"%(i,),commutative=False,real=True)
        operator_table.append(ckdi)
        operator_table_minus.append(cmki)
    operator_table=[operator_table,operator_table_minus]
    operator_table=N.ravel(operator_table)
    return operator_table

def gen_XdX(atom_list,operator_table,operator_table_dagger,Hcomm,N_atoms_uc):
    """Operate commutation relations to put all the 2nd order term as ckd**ck, cmk**cmkd, cmk**ck and ckd**cmkd form"""

    exclude_list=[]
    coeff_list=[]
    Hcomm=Hcomm.expand()
    XdX=sympy.zeros(2*N_atoms_uc)
    g=sympy.zeros(2*N_atoms_uc)
    for i in range(2*N_atoms_uc):
        curr_row=[]
        for j in range(2*N_atoms_uc):
            mycoeff=operator_table_dagger[i]*operator_table[j]
            exclude_list.append(mycoeff)
            currcoeff=coeff(Hcomm,mycoeff)
            if currcoeff!=None:
                XdX[i,j]=currcoeff
                curr_row.append(currcoeff)
            if i!=j:
                g[i,j]=0
            else:
                if i>=N_atoms_uc:
                    g[i,j]=-1
                else:
                    g[i,j]=1
 
    return XdX,g



def calculate_dispersion(atom_list,N_atoms_uc,N_atoms,Jij,direction,steps,showEigs=False):
    Sabn=generate_sabn(N_atoms)       
    #print 'Sabn',Sabn 
    Sxyz=generate_sxyz(Sabn,atom_list)
    #print 'Sxyz', Sxyz
        
    if 1:
        #print len(translations)   
        J=sympy.Symbol('J',real=True)
        #Jij=[N.matrix([[J,0,0],[0,J,0],[0,0,J]])]
        #Hdef=generate_hdef(atom_list,Jij,Sabn,N_atoms_uc,N_atoms)
        Hdef=generate_hdef(atom_list,Jij,Sxyz,N_atoms_uc,N_atoms)
        print 'Hdef',Hdef
        print type(Hdef)
        #pngview(Hdef)
        #print_matplotlib(latex(Hdef)) 
    #if 0:
        Hlin=holstein(Hdef)
        print 'Hlin',Hlin
        kx=sympy.Symbol('kx',real=True)
        ky=sympy.Symbol('ky',real=True)
        kz=sympy.Symbol('kz',real=True)
        k=[kx,ky,kz]
        Hfou=fouriertransform(atom_list,Jij,Hlin,k,N_atoms_uc,N_atoms)
        print 'Hfou',Hfou
    if 1:
        Hcomm=applycommutation(atom_list,Jij,Hfou,k,N_atoms_uc,N_atoms)
        print 'Hcomm',Hcomm
        operator_table=gen_operator_table(atom_list,N_atoms_uc)
        print 'optable',operator_table
        operator_table_dagger=gen_operator_table_dagger(atom_list,N_atoms_uc)
        print 'optable_dagger',operator_table_dagger
        XdX,g=gen_XdX(atom_list,operator_table,operator_table_dagger,Hcomm,N_atoms_uc)
        print 'XdX',XdX
        print 'g',g
        TwogH2=2*g*XdX
        print 'TwogH2',TwogH2
        m,n=TwogH2.shape
        if 1:
                for i in range(m):
                    for j in range(n):
                        #print i,j
                        #print Ntwo[i,j]
                        #print 'matching'
                        #print 'kx',Ntwo[i,j].match(kx)
                        #print 'ky',Ntwo[i,j].match(ky)
                        #Ntwo[i,j]=sympy.re(Ntwo[i,j].evalf())
                        #Ntwo[i,j]=Ntwo[i,j].evalf()
                        TwogH2[i,j]=TwogH2[i,j].expand(complex=True)#.subs(I,1.0j)
        print 'trigified',TwogH2
        print TwogH2[0,1]
        
        if showEigs:
            #print 'calculating'
            x=sympy.Symbol('x')
            #eigspoly=TwogH2.berkowitz_charpoly(x)
            #print 'eigspoly'
            #print 'eigs poly',eigspoly
            print 'shape', TwogH2.shape
            print 'recalculating'
            eigs=TwogH2.eigenvals()
            x=sympy.Symbol('x')
            #eigs=TwogH2.berkowitz_charpoly(x)
            print 'eigs', eigs
            keys=eigs.keys()
            #print 'key',keys[0]
            #print keys[0].expand(complex=True)
            print 'charpoly',TwogH2.charpoly(x)
            #eigs=TwogH2.eigenvals()
            #print 'eigenvalues', sympy.simplify(eigs[1][0])        
        Hsave=TwogH2
        return Hsave
    
def calc_eigs(Hsave,direction,steps):
        kx=sympy.Symbol('kx',real=True)
        ky=sympy.Symbol('ky',real=True)
        kz=sympy.Symbol('kz',real=True)
        k=[kx,ky,kz]
        TwogH2=Hsave
        S=sympy.Symbol('S',real=True)
        D=sympy.Symbol('D',real=True)
        #TwogH2=TwogH2.subs(J,-1.0)
        TwogH2=TwogH2.subs(S,1.0)
        #eigs=TwogH2.eigenvals()
        #print 'subbed_eigs',eigs
        #TwogH2=TwogH2.subs(D,1.0)
        qrange=[]
        wrange0=[]
        wrange1=[]
        wrange=[]
        #wrangec=[]
        for q in N.arange(0,2*pi,2*pi/steps):
            
            currnum=q*direction['kx']
            #print 'currnum x', currnum
            TwogH3=TwogH2.subs(kx,currnum)
            currnum=q*direction['ky']
            #print 'currnum y',currnum
            TwogH3=TwogH3.subs(ky,currnum)
            currnum=q*direction['kz']
            TwogH3=TwogH3.subs(kz,currnum)
            #I=sympy.Symbol('I')
            Ntwo=TwogH3#.subs(I,1.0j)
            m,n=Ntwo.shape
            #print Ntwo.applyfunc(sympy.Basic.evalf)
            Nthree=N.empty([m,n],complex)
            if 1:
                for i in range(m):
                    for j in range(n):
                        #print i,j
                        #print Ntwo[i,j]
                        #print 'matching'
                        #print 'kx',Ntwo[i,j].match(kx)
                        #print 'ky',Ntwo[i,j].match(ky)
                        #Ntwo[i,j]=sympy.re(Ntwo[i,j].evalf())
                        #Ntwo[i,j]=Ntwo[i,j].evalf()
                        Nthree[i,j]=complex(Ntwo[i,j].expand(complex=True))#.subs(I,1.0j)
                        if 1:
                            if N.absolute(Nthree[i,j])<1e-5:
                                Nthree[i,j]=0
            #print 'Ntwo',Ntwo
            #print 'Nthree',Nthree
            if 1:
                
                l,v=scipy.linalg.eig(Nthree)
                #print l[1]
                qrange.append(q)
                #print 'num eigs', l.shape
                #wrange0.append(l[0])
                #wrange1.append(l[3])
                wrange.append(l)
            #print N.linalg.eigvals(Ntwo)
            #eigs=TwogH2.eigenvals()
            #print 'eigs', eigs
            #print 'eigenvalues', sympy.simplify(eigs[1][0])
        #wrange0=N.real(wrange0)
        print 'wrange',wrange
        wrange=N.real(wrange)
        #print qrange
        #print wrange
        wrange=N.array(wrange)
        
        wrange=N.real(wrange.T)
        for wrange1 in wrange:
            #pylab.plot(qrange,wrange0,'s')
            #print qrange
            #print wrange1
            #print len(qrange),len(wrange1)
            #print inum
            #inum=inum+1
            pylab.plot(qrange,wrange1,'s')
        



def multiply_ab(atom_list,Sxyz,a=0,b=0):
    N_atoms=len(atom_list)
    Sdef=0
    print 'atom_list',atom_list
    print 'Sxyz',Sxyz
    T=sympy.Symbol('T',commutative=False)
    Sij0=Sxyz[0][a]
    t=''
    c=sympy.Symbol("c%d%s"%(0,t),commutative=False)
    cd=sympy.Symbol("cd%d%s"%(0,t),commutative=False)
    t='t'
    ct=sympy.Symbol("c%d%s"%(0,t),commutative=False)
    cdt=sympy.Symbol("cd%d%s"%(0,t),commutative=False)
    Sij0=Sij0.subs(ct,c)
    Sij0=Sij0.subs(cdt,cd)  
    for i in range(1,N_atoms):
        Sdef=Sdef+Sij0*Sxyz[i][b]
 
    return Sdef




def Sfouriertransform(atom_list,Slin,k):
    N_atoms=len(atom_list)
    #Hdef=0
    #print 'atom_list',atom_list
    #print 'Sxyz',Sxyz
    #print 'Jij',Jij
    kxp=sympy.Symbol('kxp')
    kyp=sympy.Symbol('kyp')
    kzp=sympy.Symbol('kzp')
    kp=[kxp,kyp,kzp]
    wk=sympy.Symbol("wk")
    wkp=sympy.Symbol("wkp")
    t=sympy.Symbol("t")
    ri=atom_list[0].pos
    
    kmult=N.dot(k,ri)
    #kmultp=N.dot(kp,ri)
    
    ci=sympy.Symbol("c%d"%(0,),commutative=False)
    cdi=sympy.Symbol("cd%d"%(0,),commutative=False)
    cki=sympy.Symbol("ck%d"%(0,),commutative=False)
    ckdi=sympy.Symbol("ckd%d"%(0,),commutative=False)
    
    t1=sympy.exp(I*kmult)*cki
    t3=sympy.exp(-I*(kmult))*ckdi
    Slin=Slin.subs(ci,t1)
    Slin=Slin.subs(cdi,t3)


    for i in range(1,N_atoms):
        N_int=len(atom_list[i].interactions)
        #ci=sympy.Symbol("c%d"%(i,),commutative=False)
        #cdi=sympy.Symbol("cd%d"%(i,),commutative=False)
        cit=sympy.Symbol("c%dt"%(i,),commutative=False)
        cdit=sympy.Symbol("cd%dt"%(i,),commutative=False)

        cki=sympy.Symbol("ck%d"%(0,),commutative=False)
        ckdi=sympy.Symbol("ckd%d"%(0,),commutative=False)
        
        ri=atom_list[i].pos
        kmult=N.dot(k,ri)
        kmultp=N.dot(kp,ri)
        
        t2=sympy.exp(I*(kmultp-wk*t))*cki
                     
        
        t4=sympy.exp(-I*(kmultp-wkp*t))*ckdi
        
        Slin=Slin.subs(cit,t2)
        
        Slin=Slin.subs(cdit,t4)
    
    #Note that I have already assumed that k=kp because I didn't include cqp terms
    Slin=Slin.expand()
    Slin=Slin.subs(wkp,wk)
    Slin=Slin.subs(kxp,kx)
    Slin=Slin.subs(kyp,ky)
    Slin=Slin.subs(kzp,kz)
    return Slin

def Sapplycommutation(atom_list,Sfou,k):
    """Operate commutation relations to put all the 2nd order term as ckd**ck, cmk**cmkd, cmk**ck and ckd**cmkd form"""
    N_atoms=len(atom_list)
    #Hdef=0
    #print 'atom_list',atom_list
    #print 'Sxyz',Sxyz
    #print 'Jij',Jij
    for i in range(N_atoms):
        N_int=len(atom_list[i].interactions)
        ci=sympy.Symbol("c%d"%(0,),commutative=False)
        cdi=sympy.Symbol("cd%d"%(0,),commutative=False)
        cki=sympy.Symbol("ck%d"%(0,),commutative=False)
        ckdi=sympy.Symbol("ckd%d"%(0,),commutative=False)
        cmki=sympy.Symbol("cmk%d"%(0,),commutative=False)
        cmkdi=sympy.Symbol("cmkd%d"%(0,),commutative=False)
        for j in range(N_atoms):
            cj=sympy.Symbol("c%d"%(0,),commutative=False)
            cdj=sympy.Symbol("cd%d"%(0,),commutative=False)
            ckj=sympy.Symbol("ck%d"%(0,),commutative=False)
            ckdj=sympy.Symbol("ckd%d"%(0,),commutative=False)
            cmkj=sympy.Symbol("cmk%d"%(0,),commutative=False)
            cmkdj=sympy.Symbol("cmkd%d"%(0,),commutative=False)
            Sfou=Sfou.expand()
            if i==j:
                Sfou=Sfou.subs(cki*ckdj,ckdj*cki+1)
                #Sfou.subs(cmkdi*cmkj,cmkj*cmkdi)
                nkj=sympy.Symbol("nk%d"%(j,),commutative=True)
                
            else:
                Sfou=Sfou.subs(cki*ckdj,ckdj*cki) #just added
            
            Sfou=Sfou.expand()
            nkj=sympy.Symbol("nk%d"%(j,),commutative=True)
            Sfou=Sfou.subs(ckdj*cki,nkj)
            Sfou=Sfou.expand()
            Sfou=Sfou.subs(cki*ckj,0)
            Sfou=Sfou.expand()
            Sfou=Sfou.subs(ckdi*ckdj,0)

    
    
    return Sfou

def driver(spinfile,interactionfile,direction,steps):
    myfilestr=spinfile#r'c:\spins.txt'
    myspins=readfiles.read_spins(myfilestr)
    spins=readfiles.find_collinear(myspins)
    print 'driver spins',spins
    myfilestr=interactionfile#r'c:\montecarlo.txt'
    atom_list, jnums, jmats,N_atoms_uc=readfiles.read_interactions(myfilestr,spins)
    N_atoms=len(atom_list)
    #N_atoms_uc=1
    print 'N_atoms',N_atoms,'Natoms_uc',N_atoms_uc
    #atom_list=generate_atoms()
    #atom_list=generate_atoms_rot()
    Hsave=calculate_dispersion(atom_list,N_atoms_uc,N_atoms,jmats,direction,steps,showEigs=True)
    calc_eigs(Hsave,direction,steps)
    direction={}
    direction['kx']=0.
    direction['ky']=1.
    direction['kz']=0.
    #pylab.figure()
    #calc_eigs(Hsave,direction,steps)
    
    pylab.show()
    
    print jmats
    print direction
    print steps
    print spinfile
    print interactionfile
    
if __name__=='__main__':
    if 1:
        #translations=generate_translations()
        #atom_list=generate_atoms()
        #N_atoms_uc=1
        #N_atoms=2
         #print spins[0]
        #print spins[1]
        #print N.linalg.det(spins[0]), N.linalg.det(spins[1])
        spinfile=r'c:\spins.txt'
        #spinfile=r'c:\spinsp1.txt'
        #spins=readfiles.read_spins(myfilestr)
        interactionfile=r'c:\montecarlo.txt'
        #interactionfile=r'c:\montep11.txt'
        #interactionfile=r'c:\montecarlop1.txt'
        steps=100
        data={}
        data['kx']=1.
        data['ky']=0.
        data['kz']=0.
        direction=data
        driver(spinfile,interactionfile,direction,steps)
        #atom_list, jnums, jmats=readfiles.read_interactions(myfilestr,spins)
        #N_atoms=len(atom_list)
        #N_atoms_uc=1
        #print N_atoms
        #atom_list=generate_atoms()
        #calculate_dispersion(atom_list,N_atoms_uc,N_atoms,jmats)
        #print jmats



    if 0:
        print 'one magnon'
        print ''
        print ''
        Sabnt=generate_sabnt(N_atoms,t='t')
        SzSz=sympy.expand(multiply_ab(atom_list,Sabnt,a=2,b=2))
        print 'mult SzSz',SzSz
        SxSx=sympy.expand(multiply_ab(atom_list,Sabnt,a=0,b=0))
        SxSy=sympy.expand(multiply_ab(atom_list,Sabnt,a=0,b=1))
        print 'mult SxSx',SxSx
        SzSz_lin=holstein(sympy.expand(SzSz))
        print 'lin zz',SzSz_lin
        SxSx_lin=holstein(sympy.expand(SxSx))
        SxSy_lin=holstein(sympy.expand(SxSy))
        print 'lin xy',SzSz_lin
        print 'lin xx',SxSx_lin        
        if SzSz_lin!=None:
            SzSz_fou=Sfouriertransform(atom_list,SzSz_lin,k)
            print 'fourier', SzSz_fou
            Scomm=Sapplycommutation(atom_list,SzSz_fou,k)
            print 'Scomm',Scomm
            SxSx_fou=Sfouriertransform(atom_list,SxSx_lin,k)
            print 'fourier x',SxSx_fou
            Scommx=Sapplycommutation(atom_list,SxSx_fou,k)
            print 'Scommx',sympy.simplify(Scommx)
            SxSy_fou=Sfouriertransform(atom_list,SxSy_lin,k)
            print 'fourier xy',SxSy_fou
            Scommxy=Sapplycommutation(atom_list,SxSy_fou,k)
            print 'Scommxy',Scommxy            
