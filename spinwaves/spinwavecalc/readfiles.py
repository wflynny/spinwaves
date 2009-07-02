import numpy as N
import solvespin
import sympy
import copy

class atom:
    def __init__(self,spinRmatrix=sympy.Matrix([[1, 0, 0],
                                            [0, 1, 0],
                                            [0, 0, 1]]),
                    pos=[0,0,0],neighbors=None,interactions=None,label=0,Dx=0,Dy=0,Dz=0,cell=0,int_cell=[], orig_Index = None):
        self.spinRmatrix=spinRmatrix#found with findmat
        if neighbors==None:
            neighbors=[]
        if interactions==None:
            interactions=[]
        self.pos=N.array(pos)
        self.neighbors=neighbors#Interacting atoms
        self.interactions=interactions#the matrices describing the interactions^
        self.label=label
        self.Dx=Dx
        self.Dy=Dy
        self.Dz=Dz
        self.cell=cell
        self.int_cell=[]
        #Indices change when reading in only specific atoms from the file
        self.origIndex = orig_Index



def get_tokenized_line(myfile,returnline=['']):
        lineStr=myfile.readline()
        returnline[0]=lineStr.rstrip()
        strippedLine=lineStr.lower().rstrip()
        tokenized=strippedLine.split()

        return tokenized


def readFiles(interactionFileStr,spinFileStr):
    """modified from read_interactions.  Originally this(read_interactions) read in the
    atoms from the interaction file and matched the spin rotation matrices with the
    appropriate atoms based on indices.  Now it takes the interaction and spin file strings,
    reads in the atom information and matches it with spin rotation matrices based on coordinates"""
    #print interactionFileStr
    interactionFile = open(interactionFileStr, 'r')
    myFlag=True
        #self.metadata={}
    returnline=['']
    jmats=[]
    jnums=[]
    atomlist=[]
    numcell=0
    #print "here"
    while myFlag:
        tokenized=get_tokenized_line(interactionFile,returnline=returnline)
        #print tokenized
        if not(tokenized):
            break
        #print tokenized
        if tokenized==[]:
            break
        if tokenized[0]=='#number':
            while 1:
                tokenized=get_tokenized_line(interactionFile,returnline=returnline)
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
                        tokenized=get_tokenized_line(interactionFile,returnline=returnline)
                        if not(tokenized):
                            break
                        #print tokenized
                        atom_num=int(tokenized[0])
                        if tokenized[1] == 'x':  #If it is in the first interaction cell
                            print "atom in first interaction cell"
                            x,y,z=float(tokenized[2]),float(tokenized[3]),float(tokenized[4])
                            if x < 1 and y < 1 and z < 1:
                                numcell += 1
                            Dx,Dy,Dz=float(tokenized[5]),float(tokenized[6]),float(tokenized[7])
                            #spin0=N.matrix([[1,0,0],[0,1,0],[0,0,1]],'float64')
                            pos0=[x,y,z]
                            atom0=atom(pos=pos0,Dx=Dx,Dy=Dy,Dz=Dz, orig_Index = atom_num)
                            neighbors=[]
                            interactions=[]
                            for i in range(8,len(tokenized),2):
                                interacting_spin=int(tokenized[i])
                                #index number in export list not necessarily the same as index
                                #in list of atoms in first interacting 'cell'
                                interaction_matrix=int(tokenized[i+1])
                                neighbors.append(interacting_spin)
                                interactions.append(interaction_matrix)

                            atom0.neighbors=neighbors
                            atom0.interactions=interactions
                            currnum=currnum+1
                            atomlist.append(atom0)
    interactionFile.close()
    
    #interactions should contain the indices of interaction atoms in atomlist, but
    #the index they currently contain is the index from the file, which may have
    #changed.  Switch to new indices:
    #print "here2"
#    for atom1 in atomlist:
#        print atom1.neighbors
#        print atom1.interactions
#        print atom1.origIndex
    
    
    #change interaction indices to match indices in new list
    for a in atomlist:
        neighborList = a.neighbors
        i = 0
        while i < len(neighborList):
            neighbor_index = neighborList[i]
            for atom_index in range(len(atomlist)):
                atom1 = atomlist[atom_index]
                if atom1.origIndex == neighbor_index:
                    a.neighbors[i] = atom_index
                    i +=1
                    break
            else:#This interaction index is not in the interaction cell
                neighborList.pop(i)
                a.interactions.pop(i)
            

    #Match spin rotation matrices to atoms
    #print "here3"
   
    spinFile = open(spinFileStr, 'r')
    lines = spinFile.readlines()
    spinFile.close()
    for atom1 in atomlist:
        #search for the spin entry that matches this atoms coordinates
        for line in lines:
            tokenized_line = line.split()
            if not tokenized_line[0][0] == '#': #skip commented line
                xPos = float(tokenized_line[1])
                yPos = float(tokenized_line[2])
                zPos = float(tokenized_line[3])
                if (atom1.pos[0] == xPos and atom1.pos[1] == yPos and atom1.pos[2] == zPos):
                    spin=N.array([float(tokenized_line[4]),float(tokenized_line[5]),float(tokenized_line[6])],'Float64')
                    #spin=[float(tokenized_line[4]),float(tokenized_line[5]),float(tokenized_line[6])]
                    rmat = findmat(spin)
                    atom1.spinRmatrix = rmat
                    break
        else:
            raise Exception()
            
            
#    for atom1 in atomlist:
#        print "neighbors: ", atom1.neighbors
#        print "\ninteractions: ", atom1.interactions
#        print "\nroation matrix: ", atom1.spinRmatrix
#This method would use less memory, but it would be a little slower
#    while True:
#        line = get_tokenized_line(spinFile)
#        if line == []:
#            break
#        if line[0] == '#': #The line is commented out
#            pass
#        pos = (float(line[1]), float(line[2]), float(line[3]))
    
    #for test purposes
    #for spin in spins:
    #    print spin
    #print "numSpins: ", len(spins)
    return atomlist, jnums, jmats, numcell
    

def read_spins(myfilestr):
    """read spins from file and return as a list of numpy.array([sx,sy,sz])"""
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
    spins=copy.deepcopy(spins)#What's the point of this?
    spins.reverse()
    #print 'spins',spins
    spin=N.array(spins.pop(),'float64')#I believe it's already an array
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