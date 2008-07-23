from simple import readFile, Timer, simpleAtom
import ctypes
from ctypes import c_float, c_int
import cstruct
import numpy as N
import sys
import time

 
  
if __name__ == '__main__':      
    k = 10
    tMax = 10
    tMin = .01
    tFactor = .90
    timer = Timer()
    inFilePath = "C:\Newexport1.txt"
    outFilePath = "C:\spins.txt"
    
    if sys.platform=='win32':
        print 'win32'
        monteCarloDll = N.ctypeslib.load_library('monteCarlo', 'C:/Dev-Cpp/workspace')
#    elif sys.platform=='mac':
#        monteCarloDll = N.ctypeslib.load_library('libpolarization2.so', '.')
#    else:
#        monteCarloDll = N.ctypeslib.load_library('libpolarization2.so', '.') #linux
        
    atoms, jMatrices = readFile(inFilePath)
    
#    print atoms
    print jMatrices
    

    successCode = c_int(0)
    atomListPointer = monteCarloDll.new_atom_list(c_int(len(atoms)), ctypes.byref(successCode))
    print "Success Code: ", successCode

    #to keep pointers
    matListList = []
    nbr_ListList = []

    #Create atom list in C
    s1Class = c_float * 3
    for i in range(len(atoms)):
#    for i in range(5):
        atom = atoms[i]
        numInteractions = len(atom.interactions)
        matListClass = c_int * numInteractions
        nbrListClass = c_int * numInteractions
        
        matList = matListClass()
        neighbors = nbrListClass()
        for j in range(numInteractions):
            neighbors[j] = c_int(atom.interactions[j][0])
            matList[j] = c_int(atom.interactions[j][1])
#            print "interaction: " , atom.interactions[j][0], atom.interactions[j][1],  " = " , neighbors[j], matList[j]
        
        s1 = s1Class()
        s1[0] = c_float(1)
        s1[1] = c_float(0)
        s1[2] = c_float(0)
        
        #to keep pointers
        matListList.append(matList)
        nbr_ListList.append(neighbors)
#        print "atom" + str(i) + ") numInteractions:", numInteractions, len(neighbors), len(matList)
        
        print str(i) + ")", atom.pos
        for interaction in atom.interactions:
            print atoms[interaction[0]].pos, " : ", interaction[1]
        monteCarloDll.set_atom(atomListPointer, c_int(i), matList, neighbors, c_int(numInteractions), s1)
    
    print "atoms added"
#    time.sleep(10)
    timer.printTime()

    #Create JMatrix List in C
    matPointer = monteCarloDll.new_jMatrix_list(c_int(len(jMatrices)), ctypes.byref(successCode))
    print "Matrix Success Code: ", successCode
#   time.sleep(5)

    
    
    print jMatrices
    print "matpointer python: ", matPointer
#    print 'mylen', len(jMatrices)
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
        print 'added'
        

#    print 'jmatrix python', jMatrices
#    monteCarloDll.del_jMat(matPointer)
#    print "jmat deleted"
#    time.sleep(5)
    
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

    
    spins = []
    for i in range(len(atoms)):
        spin = s1Class()
        monteCarloDll.getSpin(atomListPointer, c_int(i), ctypes.byref(spin))
        spins.append(spin)
#        print str(spin[0]) + " " + str(spin[1]) + " " + str(spin[2])
    
   
    monteCarloDll.del_jMat(matPointer)
    monteCarloDll.del_atom(atomListPointer)
    
    
    print "writing spins to file..."
    timer.printTime()
    
    #output the spins to a file
    outFile = open(outFilePath, 'w')
    outFile.write("#Atom_Number Position_X Position_Y Position_Z Spin_X Spin_Y Spin_Z\n")
    for i in range(len(atoms)):
        atom = atoms[i]
        Posx = str(atom.pos[0])
        Posy = str(atom.pos[1])
        Posz = str(atom.pos[2])
        spin = spins[i]
        Spinx = str(spin[0])
        Spiny = str(spin[1])
        Spinz = str(spin[2])
        atomStr = str(i) + " " + Posx + " " + Posy + " " + Posz + " " + Spinx + " " + Spiny + " " + Spinz + "\n"
        outFile.write(atomStr)
    
    outFile.close()
    
    
    
    timer.printTime()
    print "done"