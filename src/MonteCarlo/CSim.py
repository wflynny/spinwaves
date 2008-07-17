from simple import readFile, Timer, simpleAtom
import ctypes
from ctypes import c_float, c_int
import cstruct
import numpy as N
import sys

 
  
if __name__ == '__main__':      
    k = 10
    tMax = 10
    tMin = .01
    tFactor = .9
    timer = Timer()
    
    if sys.platform=='win32':
        print 'win32'
        monteCarloDll = N.ctypeslib.load_library('monteCarlo', 'C:/Dev-Cpp/workspace')
#    elif sys.platform=='mac':
#        monteCarloDll = N.ctypeslib.load_library('libpolarization2.so', '.')
#    else:
#        monteCarloDll = N.ctypeslib.load_library('libpolarization2.so', '.') #linux
        
    atoms, jMatrices = readFile("C:\Newexport1.txt")

    atomListPointer = monteCarloDll.new_atom_list(c_int(len(atoms)))
    print atomListPointer

    #to keep pointers
    matListList = []
    nbr_ListList = []

    #Create atom list in C
    for i in range(len(atoms)):
#    for i in range(5):
        atom = atoms[i]
        numInteractions = len(atom.interactions)
        matListClass = c_int * numInteractions
        nbrListClass = c_int * numInteractions
        s1Class = c_float * 3
        
        matList = matListClass()
        neighbors = nbrListClass()
        for j in range(numInteractions):
            neighbors[j] = c_int(atom.interactions[j][0])
            matList[j] = c_int(atom.interactions[j][1])
#            print "interaction: " , atom.interactions[j][0], atom.interactions[j][1],  " = " , neighbors[j], matList[j]
        
        s1 = s1Class()
        s1[0] = c_float(0)
        s1[1] = c_float(0)
        s1[2] = c_float(0)
        
        #to keep pointers
        matListList.append(matList)
        nbr_ListList.append(neighbors)
        
        monteCarloDll.set_atom(atomListPointer, c_int(i), matList, neighbors, c_int(numInteractions), s1)
    
    print "atoms added"
    timer.printTime()

    #Create JMatrix List in C
    matPointer = monteCarloDll.new_jMatrix_list(c_int(len(jMatrices)))
    
    print jMatrices
    print "matpointer: ", matPointer
    
    for i in range(len(jMatrices)):
        j = jMatrices[i]
        j11 = c_float(j[0][0])
        j12 = c_float(j[0][1])
        j13 = c_float(j[0][2])
        j21 = c_float(j[1][0])
        j22 = c_float(j[1][1])
        j23 = c_float(j[1][2])
        j31 = c_float(j[2][0])
        j32 = c_float(j[2][1])
        j33 = c_float(j[2][2])
        monteCarloDll.add_matrix(matPointer, c_int(i),j11, j12, j13, j21, j22, j23, j31, j32, j33)
    


    monteCarloDll.simulate(atomListPointer, c_int(len(atoms)), matPointer, c_int(k), c_float(tMax), c_float(tMin), c_float(tFactor))


    #for test purposes
#    monteCarloDll.atomTest(listPointer, c_int(len(atoms)))
    
#    print "atom1:"
#    for interaction in atoms[len(atoms)-3].interactions:
#        print interaction[0], " , " , interaction[1]
#    
#    print "atom2:"
#    for interaction in atoms[len(atoms)-2].interactions:
#        print interaction[0], " , " , interaction[1]
#    
#    print "atom3:"
#    for interaction in atoms[len(atoms)-1].interactions:
#        print interaction[0], " , " , interaction[1]
#        
#        
#    for j in jMatrices:
#        print j
#    monteCarloDll.matrixTest(matPointer, c_int(len(jMatrices)))
    
    monteCarloDll.del_jMat(matPointer, c_int(len(jMatrices)))
    monteCarloDll.del_atom(atomListPointer, c_int(len(atoms)))
    
    timer.printTime()