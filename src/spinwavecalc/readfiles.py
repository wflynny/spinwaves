import numpy as N
import solvespin
import sympy
import copy

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
        atom0=atom(spin=spin0,pos=pos0,neighbors=neighbors,interactions=interactions,label=0,cell=cell,int_cell=int_cell,Dz=D)
        
        pos0=[1,0,0]
        spin0=N.matrix([[-1,0,0],[0,1,0],[0,0,-1]],'float64')
        neighbors=[0]
        interactions=[0]
        cell=5
        int_cell=[0]
        atom1=atom(spin=spin0,pos=pos0,neighbors=neighbors,interactions=interactions,label=1,cell=cell,int_cell=int_cell,Dz=D)
        
        atomlist=[atom0,atom1]


 
       
    return atomlist

def get_tokenized_line(myfile,returnline=['']):
        lineStr=myfile.readline()
        returnline[0]=lineStr.rstrip()
        strippedLine=lineStr.lower().rstrip()
        tokenized=strippedLine.split()

        return tokenized


def read_interactions(myfilestr,spins):
    myfile = open(myfilestr, 'r')
    myFlag=True
        #self.metadata={}
    returnline=['']
    jmats=[]
    jnums=[]
    atomlist=[]
    numcell=0
    while myFlag:
        tokenized=get_tokenized_line(myfile,returnline=returnline)
        #print tokenized
        if not(tokenized):
            break
        #print tokenized
        if tokenized==[]:
            break
        if tokenized[0]=='#number':
            while 1:
                tokenized=get_tokenized_line(myfile,returnline=returnline)
                #print 'intoken ',tokenized
                if tokenized==[]:
                    break
                if tokenized[0]!='#atomnumber':
                    #print tokenized[0]
                    jnum=float(tokenized[0])
                    j11=float(tokenized[1])
                    j12=float(tokenized[2])
                    j13=float(tokenized[3])
                    j21=float(tokenized[4])
                    j22=float(tokenized[5])
                    j23=float(tokenized[6])
                    j31=float(tokenized[7])
                    j32=float(tokenized[8])
                    j33=float(tokenized[9]) 
                    #jij=N.matrix([[j11,j12,j13],[j21,j22,j23],[j31,j32,j33]],'Float64')
                    jij=sympy.matrices.Matrix([[j11,j12,j13],[j21,j22,j23],[j31,j32,j33]])
                    jnums.append(jnum)
                    jmats.append(jij)
                else:
                    currnum=0
                    while 1:
                        tokenized=get_tokenized_line(myfile,returnline=returnline)
                        if not(tokenized):
                            break
                        atom_num=tokenized[0]
                        x,y,z=float(tokenized[1]),float(tokenized[2]),float(tokenized[3])
                        Dx,Dy,Dz=float(tokenized[4]),float(tokenized[5]),float(tokenized[6])
                        #spin0=N.matrix([[1,0,0],[0,1,0],[0,0,1]],'float64')
                        pos0=[x,y,z]
                        if N.abs(x)<1.0 and N.abs(y)<1.0 and N.abs(z)<1.0:
                            numcell=numcell+1
                        atom0=atom(pos=pos0,Dx=Dx,Dy=Dy,Dz=Dz)
                        neighbors=[]
                        interactions=[]
                        #print 'range',range(7,len(tokenized),1)
                        for i in range(7,len(tokenized),2):
                            interacting_spin=int(tokenized[i])
                            #print interacting_spin
                            interaction_matrix=int(tokenized[i+1])
                            neighbors.append(interacting_spin)
                            interactions.append(interaction_matrix)
                        #print 'interactions', interactions
                        #print 'neighbors', neighbors
                        atom0.neighbors=neighbors
                        atom0.interactions=interactions
                        atom0.spin=spins[currnum]
                        currnum=currnum+1
                        atomlist.append(atom0)
    myfile.close()
    #for catom in atomlist:
    #    print 'pos', catom.pos
    #    print 'Dx,Dy,Dz',catom.Dx, catom.Dy,catom.Dz
    #    print 'interactions', catom.interactions
    #    print 'neighbors', catom.neighbors
    #print 'jnums', jnums
    #print 'jmats',jmats
    return atomlist, jnums, jmats,numcell
    

def read_spins(myfilestr):
    myfile = open(myfilestr, 'r')
    returnline=['']
    myFlag=True
        #self.metadata={}
    spins=[]
    while myFlag:
        tokenized=get_tokenized_line(myfile,returnline=returnline)
        #print tokenized
        if not(tokenized):
            break
        #print tokenized
        if tokenized[0]=='#atom_number':
            pass
        else:
            spin=N.array([float(tokenized[4]),float(tokenized[5]),float(tokenized[6])],'Float64')
            sx,sy,sz=spin
            spins.append(spin)
            #sm=N.sqrt(sx**2+sy**2+sz**2)
            #sx=sx/sm
            #sy=sy/sm
            #sz=sz/sm
            #print 'sx',sx,'sy',sy,'sz',sz
            #smat=solvespin.getmatrix(sx, sy, sz)
            #spins.append(smat)
    myfile.close()
    if 0:
        smat=N.empty((3,1),'float64')
        smat=N.matrix(smat,'float64')
        smat[0]=0
        smat[1]=0
        smat[2]=1
        #sout=spins[1]*smat
        #spins[1]=N.matrix([[1,0,0],[0,1,0],[0,0,1]],'Float64')
        #spins[0]=N.matrix([[1,0,0],[0,1,0],[0,0,1]],'Float64')
        #print sout
        icount=0
        for currspin in spins:
            print 'spin',icount,currspin*smat
            print 'mat', currspin
            print 'det', N.linalg.det(currspin)
            icount=icount+1
    return spins
    
def findmat(spin):
        sx,sy,sz=spin
        sm=N.sqrt(sx**2+sy**2+sz**2)
        sx=sx/sm
        sy=sy/sm
        sz=sz/sm
        #print 'sx',sx,'sy',sy,'sz',sz
        smat=solvespin.getmatrix(sx, sy, sz)
        return smat
        #spins.append(smat)
  
def iscollinear(a,b,eps=1e-6):
    cross=N.abs(N.cross(a,b))<=eps
    if cross[0]*cross[1]*cross[2]:
        return True
    else:
        return False

def isparallel(a,b,eps=1e-4):
    _a=N.array(a,'float64')/N.sqrt(N.dot(a,a))
    _b=N.array(b,'float64')/N.sqrt(N.dot(b,b))
    pr=N.dot(a,b)
    if pr>=1-eps:
        return 1
    if pr<=-1+eps:
        return -1
    else:
        return 0

class Collinear_group():
    def __init__(self,parallel=None,antiparallel=None):
        #print 'init'
        if parallel==None:
            parallel=[]
        if antiparallel==None:
            antiparallel=[]
        self.parallel=parallel
        self.antiparallel=antiparallel
        self.rmatrix=None
        self.armatrix=None
        #print 'parinit',self.parallel,parallel




def find_collinear(spins):
    #print 'find collinear'
    collinear_groups=[]
    rmatrices=[]
    flag=True
    spins=copy.deepcopy(spins)
    spins.reverse()
    #print 'spins',spins
    spin=N.array(spins.pop(),'float64')
    #print '1st spin',spin
    rmat=findmat(spin)
    myg=Collinear_group()
    myg.parallel.append(copy.deepcopy(spin))
    myg.rmatrix=rmat
    rmatrices.append(rmat)
    collinear_groups.append(myg)
    #print collinear_groups[0].parallel
    #print 'remaining', spins
    while flag:
        try:   
            spin=spins.pop()
            print 'spin',spin
            #print 'remain',spins
            ct=0
            for currgroup in collinear_groups:
                #print 'currgroup', currgroup.parallel
                test_spin=copy.deepcopy(currgroup.parallel[0])
                #print 'test spin',test_spin
                if iscollinear(spin,test_spin):
                    ispar=isparallel(spin,test_spin)
                    if ispar==1:
                        #print 'par'
                        currgroup.parallel.append(copy.deepcopy(spin))
                        ct=ct+1
                        rmatrices.append(currgroup.rmatrix)
                        break
                    elif ispar==-1:
                        #print 'antipar'
                        currgroup.antiparallel.append(copy.deepcopy(spin))
                        ct=ct+1
                        if currgroup.armatrix==None:
                            rmat=findmat(spin)
                            currgroup.armatrix=rmat
                        rmatrices.append(currgroup.armatrix)
                        break              
            if ct==0:
                #print 'new group'
                mycol=Collinear_group()
                mycol.parallel.append(copy.deepcopy(spin))
                rmat=findmat(spin)
                mycol.rmatrix=rmat
                rmatrices.append(rmat)
                collinear_groups.append(mycol)
                #for cgroup in collinear_groups:
                #    print 'collinear_groups par',cgroup.parallel
                #    print 'collinear_groups antipar',cgroup.antiparallel

        except IndexError:
            flag=False
    
    if 0:
        print 'final'
        for currgroup in collinear_groups:
            print 'groups'
            print 'parallel',currgroup.parallel
            print 'antiparallel',currgroup.antiparallel
    return rmatrices
        
    
    
    
if __name__=="__main__":
    myfilestr=r'c:\spins.txt'
    spins=read_spins(myfilestr)
    a=[1,0,0]
    b=[-2,0,0]
    print iscollinear(a,b)
    print isparallel(a,b)
    myspins=[]
    myspins.append(N.array([1,0,0],'float64'))
    myspins.append(N.array([-1,0,0],'float64'))
    myspins.append(N.array([0,1,0],'float64'))
    myspins.append(N.array([2,0,0],'float64'))
    myspins.append(N.array([0,-2,0],'float64'))
    myspins.append(N.array([1,1,0],'float64'))
    myspins=spins
    rmatrices=find_collinear(myspins)
    print rmatrices
    if 1:
        smat=N.empty((3,1),'float64')
        smat=N.matrix(smat,'float64')
        smat[0]=0
        smat[1]=0
        smat[2]=1
        #sout=spins[1]*smat
        #spins[1]=N.matrix([[1,0,0],[0,1,0],[0,0,1]],'Float64')
        #spins[0]=N.matrix([[1,0,0],[0,1,0],[0,0,1]],'Float64')
        #print sout
        icount=0
        for rmat in rmatrices:
            print 'spin',icount,rmat*smat,'realspin',myspins[icount]
            print 'mat', rmat
            print 'det', N.linalg.det(rmat)
            icount=icount+1

    #myfilestr=r'c:\montecarlo.txt'
    if 0:
        myfilestr=r'c:\montep11.txt'
        atomlist, jnums, jmats,numcell=read_interactions(myfilestr,spins)
        print 'numcell', numcell
        print jmats
        print atomlist[0].interactions
        print atomlist[1].neighbors
        #print N.linalg.det(atomlist[0].spin)