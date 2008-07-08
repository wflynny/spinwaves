import numpy
import time
import random
import math
import pylab



class Timer():
    """for diagnostics"""
    def __init__(self):
        self.initialTime = time.clock()
    
    def printTime(self):
        seconds = time.clock() - self.initialTime
        minutes = int(seconds)/60
        seconds -= (minutes*60)
        hours = minutes/60
        minutes -= hours*60
        print "Time:", hours, "hours", minutes, "minutes", seconds, "seconds" 
 
    

class simpleAtom():
    def __init__(self, pos):
        self. pos = pos
        self.s = numpy.array([0,0,random.uniform(-1,1)])
        self.interactions = []

def readFile(filename):
    #read in the interactions file
    file = open(filename, 'r')
    
    
    #This will use a huge amount of memory
    atoms = []
    lines = file.readlines()
    #Read in J matrices
    index = 2
    line = lines[index]
    jMatrices = []
    while not line.startswith('#'):
        values = line.split()
    #    print values
        jMat = numpy.array([[float(values[1]), float(values[2]), float(values[3])],
                            [float(values[4]), float(values[5]), float(values[6])],
                            [float(values[7]), float(values[8]), float(values[9])]])
        jMatrices.append(jMat)
        if not (jMatrices[int(values[0])] == jMat).all():
            print "The Indeces are messed up!!!" #For now
        index +=1
        line = lines[index]
    
    for i in range(index, len(lines)):
        line = lines[i]
        if not line.startswith('#'):
            values = line.split()
#            print values
#            time.sleep(2)
            newAtom = simpleAtom((float(values[1]), float(values[2]), float(values[3])))
#            print "atom pos:", newAtom.pos
#            time.sleep(5)
            i = 4
            while i < len(values):
                otherAtomIndex = int(values[i])
                otherAtomPos = (float(values[i+1]), float(values[i+2]), float(values[i+3]))
#                print "other atom pos:", otherAtomPos
                jMatInt = int(values[i+4])
#               print "jamtrix", jMatInt
#               time.sleep(20)
                i += 5
                newAtom.interactions.append([otherAtomIndex, jMatInt])
            
            if int(values[0]) != len(atoms):
                print "problem, indeces don't match up!"
            atoms.append(newAtom)
            
    return atoms, jMatrices
  
def Energy(atom, S):
    E = 0
    for interaction in atom.interactions:
        s1 = S
        s2 = atoms[interaction[0]].s
        Jij = jMatrices[interaction[1]]
        E += numpy.dot(numpy.dot(s1,Jij),s2)
    return E
        

def flipSpins():
    randGen = random.Random()
    for atom in atoms:
        newS = numpy.array([0,0,randGen.uniform(-1,1)])
#        print newS
        oldS = atom.s
        newE = Energy(atom, newS)
        oldE = Energy(atom, oldS)
        if newE <= oldE:
            atom.s = newS
        else:
            deltE = newE - oldE
            probChange = math.exp(-deltE/T)
#            print probChange
            if randGen.random() < probChange:
                atom.s = newS
            #Otherwise it does not change
            
        
  
  
  
if __name__ == '__main__':      
    k = 1000
    tMax = 20
    tMin = .01
    tFactor = .95
    timer = Timer()
    atoms, jMatrices = readFile("C:\export.txt")
    print "File Read"
    timer.printTime()
    #double check
       #Check the simple atom list
#    checkBalanced(atoms)
    testFile = open("C:/testResults.txt", 'w')

    for i in range(10):
        T = tMax
        while T > tMin:  
            for i in range(k):
                flipSpins()
    #            pylab.plot([0], [atoms[0].s[2]], 'ro')
    #            pylab.axis([-1, 1, -1, 1])
    #            pylab.show()
            T = tFactor*T
            print "T:", T
            print "50:", atoms[50].s
            print "500:", atoms[500].s
            print "2000:", atoms[2000].s
            print "5000:", atoms[5000].s
            print "7500:", atoms[7500].s
            print "10000:", atoms[10000].s
            
            for atom in atoms:
                testFile.write(str(atom.s[2]) + "     ")
        
            #Just to see how much they change
        timer.printTime()
        testFile.write('\n')
        
    testFile.close()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    def checkBalanced(atoms):  
        def atomBalanced(atom):
            for otherAtom in atoms:
                if atomInteractsWithAtom(atom, otherAtom):
                    if atomInteractsWithAtom(otherAtom, atom):
                        return True
                    else:
                        return False
            return False
                  
        def atomInteractsWithAtom(atom, otherAtom):
            for interaction in otherAtom.interactions:
                if atom == atoms[interaction[0]]:
                    return True
            return False
        
        for atom in atoms:
            if not atomBalanced(atom):
                print "Not Balanced!!!"
                break
        else:
            print "Balanced!"