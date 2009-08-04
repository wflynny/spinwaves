import numpy
import time
import random
import math
#import pylab



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
    def __init__(self, pos, anisotropy, spinMag):
        self.pos = pos
        self.anisotropy = anisotropy
        self.spinMag = spinMag
        self.s = numpy.array([0,0,random.uniform(-1,1)])
        self.interactions = []



def readFile(filename):
    """This reads an export file from the main program and creates a list of
    simpleAtoms and jMatrices."""
    #read in the interactions file
    file = open(filename, 'r')
    
    #This will use a huge amount of memory
    atoms = []
    lines = file.readlines()#memory problem with huge files?
    
    #Read in J matrices
    index = 0
    while(lines[index] != "#Number J11 J12 J13 J21 J22 J23 J31 J32 J33\n"):
        index += 1
    index += 1 #now at line after "#Number J11 J12 J13 J21 J22 J23 J31 J32 J33"
    
    line = lines[index]
    jMatrices = []
    #This would stop reading if it hits a comment, rather than continuing past it
    #Should maybe include more after # as identifier, so comments can be added
    while not line.startswith('#'):
        values = line.split()
    #    print values
        jMat = numpy.array([[float(values[1]), float(values[2]), float(values[3])],
                            [float(values[4]), float(values[5]), float(values[6])],
                            [float(values[7]), float(values[8]), float(values[9])]])
        jMatrices.append(jMat)
        if not (jMatrices[int(values[0])] == jMat).all():
            print "The Indices are messed up!!!" #For now
        index +=1
        line = lines[index]
    
    for i in range(index, len(lines)):
        line = lines[i]
        if not line.startswith('#'):
            values = line.split()
            #skip the atom index and Whether it is n the first interaction cell or not
            #--also skip label, valence, cellnum, and atomic num
            newAtom = simpleAtom((float(values[6]), float(values[7]), float(values[8])), (float(values[9]), float(values[10]), float(values[11])), float(values[12]))
#            print "atom pos:", newAtom.pos
#            time.sleep(5)
            i = 13
            while i < len(values):
                otherAtomIndex = int(values[i])
   #format changed to no longer include position of other atom
#                otherAtomPos = (float(values[i+1]), float(values[i+2]), float(values[i+3]))
#                print "other atom pos:", otherAtomPos
#                jMatInt = int(values[i+4])
                jMatInt = int(values[i+1])
#                print "other atom index = ", otherAtomIndex
#                print "jamtrix", jMatInt
#                time.sleep(20)
#                i += 5
                i += 2
                newAtom.interactions.append([otherAtomIndex, jMatInt])
            
            if int(values[0]) != len(atoms):
                print "problem, indices don't match up!"
            atoms.append(newAtom)
#            print "newAtom:", newAtom.pos, newAtom.interactions
#            time.sleep(5)
    file.close()
    

    return atoms, jMatrices
  


#This has been moved to a Dll written in C (for speed)
def Energy(atom, S):
    E = 0
    for interaction in atom.interactions:
        s1 = S
        s2 = atoms[interaction[0]].s
        Jij = jMatrices[interaction[1]]
        E -= numpy.dot(numpy.dot(s1,Jij),s2) #E = -sum(S1*Jij*S2)
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
    k = 100
    tMax = 10
    tMin = .01
    tFactor = .9
    timer = Timer()
    atoms, jMatrices = readFile("C:\export.txt")
    print "File Read"
    timer.printTime()
    #double check
       #Check the simple atom list
#    checkBalanced(atoms)
    testFile = open("C:/testResults100.txt", 'w')

    magnetizations = []
    temperatures = []
    T = tMax
    while T > tMin:  
        for i in range(k):
            flipSpins()
#            pylab.plot([0], [atoms[0].s[2]], 'ro')
#            pylab.axis([-1, 1, -1, 1])
#            pylab.show()
        T = tFactor*T
#        print "T:", T
#        print "50:", atoms[50].s
#        print "500:", atoms[500].s
#        print "2000:", atoms[2000].s
#        print "5000:", atoms[5000].s
#        print "7500:", atoms[7500].s
#        print "10000:", atoms[10000].s
        
        #Add average magnetization
        avgMag = 0.0
        for atom in atoms:
            avgMag += atom.s[2]
        avgMag = avgMag/len(atoms)
        magnetizations.append(avgMag)
        temperatures.append(T)
    
        #Just to see how much they change
        timer.printTime()
        

    for i in range(len(atoms)):
        atom = atoms[i]
        testFile.write(str(i) + " " + str(atom.s[2]) + "   ")
    testFile.write('\n')
    testFile.close()
    #find max value
    max = 0
    min = 0
    for num in magnetizations:
        if num > max:
            max = num
        elif num < min:
            min = num

#    pylab.plot(range(len(magnetizations)), magnetizations, 'ro')
#    pylab.axis([0, len(magnetizations), min, max])
#    pylab.plot(temperatures, magnetizations, 'ro')
#    pylab.axis([temperatures[len(temperatures)-1], temperatures[0], min, max])
    #    savefig('secondfig.png')
#    pylab.show()
        
    
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
