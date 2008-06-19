   def generateMagnetCell(self, Na, Nb, Nc, renderer):
        """draws the magnet cell based on this unit cell and Na, Nb, Nc"""
        newAtoms = []
        newBonds = []
        
        for Nan in Na:
            for Nbn in Nb:
                for Ncn in Nc:
                    for atomn in Atoms:
                        renderer.AddActor()
                        
                        


       
        
        #copy bonds
        for bondn in self.Bonds:
            new_Bond = copy.deepcopy(bondn)
            new_Bond.setAtom1(new_cell.Atoms[self.Atoms.index(bondn.getAtom1())])
            new_Bond.setAtom2(new_cell.Atoms[self.Atoms.index(bondn.getAtom2())])
            
            
       
    def test(self, renderer, x):
        new_cell = Cell(x)
#        for atomn in self.Atoms:  #should preserve order of Atoms
        for i in range(0, len(self.Atoms)):
            atomn = self.Atoms[i]
            position = atomn.getPosition()
            color = atomn.getActor().GetProperty().GetColor()
            new_cell.addAtom(Atom(new_cell, position[0] + x, position[1], position[2], atomn.getDescription(), atomn.getSource().GetRadius(), color[0], color[1], color[2]))




print "here1"
cell = unitcell.translateCell(.1,0,0)

cell.drawCell(ren1)
print "here2"
print cell
print unitcell

for i in range(0, len(cell.Atoms)):
    print cell.Atoms[i]

for i in range(0, len(cell.Atoms)):
    print cell.Bonds[i]
    
print "here3"
for i in range(0, len(unitcell.Atoms)):
    print unitcell.Atoms[i]

for i in range(0, len(unitcell.Atoms)):
    print unitcell.Bonds[i]

print "here4"

cell.Atoms.pop()
cell.Atoms.pop()
cell.Atoms.pop()
cell.Atoms.pop()
for i in range(0, len(cell.Atoms)):
    print cell.Atoms[i]

for i in range(0, len(cell.Atoms)):
    print cell.Bonds[i]
    
print "here3"
for i in range(0, len(unitcell.Atoms)):
    print unitcell.Atoms[i]

for i in range(0, len(unitcell.Atoms)):
    print unitcell.Bonds[i]
    
    
    
    



def translateCell(cell, a, b, c):
    new_cell = Cell(a,b,c)
    for i in range(0, len(cell.Atoms)):
        atomn = cell.Atoms[i]
        position = atomn.getPosition()
        color = atomn.getActor().GetProperty().GetColor()
        new_cell.addAtom(Atom(new_cell, position[0], position[1], position[2], atomn.getDescription(), atomn.getSource().GetRadius(), color[0], color[1], color[2]))
    
    for i in range(0, len(cell.Bonds)):
        bondn = cell.Bonds[i]
        newAtom1 = new_cell.Atoms[ cell.Atoms.index(bondn.getAtom1()) ]
        newAtom2 = new_cell.Atoms[ cell.Atoms.index(bondn.getAtom2()) ]
        new_cell.addBond(Bond(new_cell, newAtom1, newAtom2, bondn.getRGBColor()))
    return new_cell





#Cell Class code when Atoms only dealt with their coordinates withthin their cell
   def drawCell(self, renderer):
        for atomn in self.Atoms:
#            AtomPos = atomn.getPosition()
 #           atomn.getActor().SetPosition(AtomPos[0] + self.PosX, AtomPos[1] + self.PosY, AtomPos[2] + self.PosZ)
            renderer.AddActor(atomn.getActor())
        for bondn in self.Bonds:
#           x,y,z = bondn.getActor().GetPosition()
#            bondn.getActor().SetPosition(x + self.PosX, y + self.PosY, z + self.PosZ)
            renderer.AddActor(bondn.getActor())
            
    def translateCell(self, a, b, c):
        new_cell = Cell(a,b,c)
        print len(self.Atoms)
        for atomn in self.Atoms:  #should preserve order of Atoms
            print "this"
            position = atomn.getPosition()
            color = atomn.getActor().GetProperty().GetColor()
            new_cell.addAtom(Atom(new_cell, position[0], position[1], position[2], atomn.getDescription(), atomn.getSource().GetRadius(), color[0], color[1], color[2]))
        
        print "here  in translate cell atoms done"
        for bondn in self.Bonds:
            newAtom1 = new_cell.Atoms[ self.Atoms.index(bondn.getAtom1()) ]
            newAtom2 = new_cell.Atoms[ self.Atoms.index(bondn.getAtom2()) ]
            new_cell.addBond(Bond(new_cell, newAtom1, newAtom2, bondn.getRGBColor()))
        
        return new_cell
    
    
    



class InterCellularBond(Bond):
    def __init__(self, Magnetic_Cell, Atom1, Atom2, r = 0,g = 0,b = 1):
        self.r = r#color
        self.g = g
        self.b = b
        
        self.Atom1 = Atom1
        self.Atom2 = Atom2
        self.Mag_cell = Magnetic_Cell
        self.cylinder_actor = self.makeCylinder(Atom1.getActor(), Atom2.getActor(), Atom1.getSource().GetRadius(), Atom2.getSource().GetRadius())
    
    def __str__(self):
        return "Bond between " + self.Atom1.__str__() + " and " + self.Atom2.__str__()
    
    
    
    
    


#        self.transformCellBondsToInterCell()
        
    
    #Can't just take random bond
    def transformCellBondsToInterCell(self):
        """transforms bond in the unit cell to bonds between cells using SymOps"""
        #pick a random bond from the first unit cell
        unitCellBond = self.unit_cell.getBonds()[0]
        if unitCellBond == None:
            return
        
        #Using most of the same code from addInterCellularBond
        #Only the atoms are in  the same cell and each translated pair needs
        #to be checked to see if they are in the same cell and discarded if they are
        newBonds = [] #bonds created by symmetry
        
        #Create Symmetry Bonds
        xyz = unitCellBond.getAtom1().getPosition()
        xyz2 = unitCellBond.getAtom2().getPosition()
        for symop in self.space_Group.iter_symops():
        # operate on coordinates in non-shifted spacegroup
            pos1 = symop(xyz)
            pos2 = symop(xyz2)
            mask1 = numpy.logical_or(pos1 < 0.0, pos1 >= 1.0)
            translation = numpy.floor(pos1[mask1])  #translates the first atom back to cell at (0,0,0)
            pos1[mask1] -= translation
            pos2[mask1] -= translation  #Uses same translation to translate other atom
            
            #if the second atom is not in the magnetic cell, do another translation so that it is
            if pos2[0] < 0 or pos2[1] < 0 or pos2[2] < 0:
                mask2 = pos2 < 0.0
                translation = numpy.floor(pos2[mask2])  #translates the first atom back to cell at (0,0,0)
                pos1[mask2] -= translation
                pos2[mask2] -= translation  #Uses same translation to translate other atom
                
                #If the atoms are in the same cell, discard
                if self.positionsInSameCell(pos1, pos2):
                    break
            
                atomAtPos1 = self.unit_cell.atomAtPosition(pos2) #Second Position was translated to the first cell (0,0,0)
                atomAtPos2 = self.atomAtPosition(pos1)
            else:
                if self.positionsInSameCell(pos1, pos2): #Ensuring bond is between cells
                    break
                
                atomAtPos1 = self.unit_cell.atomAtPosition(pos1) #first Position was translated to the first cell (0,0,0)
                atomAtPos2 = self.atomAtPosition(pos2)
            

#Both Atoms shouldn't be translated outside the bounds of the Magnetic Cell
            
            #Bond is between cells
            newBond = Bond(None, atomAtPos1, atomAtPos2)
            #check if the bond already exists
            for currentBond in newBonds:
                if newBond.sameBond(currentBond):
                    break
            else:  #if not, add the bond to the list of unique bonds
                newBonds.append(newBond)
                
                #translate new Bond to all cells
                originalAtom1 = newBond.getAtom1()
                originalAtom2 = newBond.getAtom2()
                origPos1 = originalAtom1.getPosition()
                origPos2 = originalAtom2.getPosition()
                #Using the fact that the translated Cells store the atoms in the same order:
                for i in range(0, self.Na - newBond.getAtom2().getUnitCell().getPosition()[0]): #translate in x direction (Na - Cell X position) times
                    for j in range(0, self.Nb - newBond.getAtom2().getUnitCell().getPosition()[1]): #translate in y direction (Nb - Cell Y position) times
                        for k in range(0, self.Nc - newBond.getAtom2().getUnitCell().getPosition()[2]): #translate in z direction (Nc - Cell Z position) times
                            translatedCell1 = self.cellAtPosition((i + origPos1[0],j + origPos1[1],k + origPos1[2]))
                            translatedCell2 = self.cellAtPosition((i + origPos2[0],j + origPos2[1],k + origPos2[2]))
                
                            translatedAtom1 = translatedCell1.atomAtIndex(originalAtom1.getIndexNumber())
                            translatedAtom2 = translatedCell2.atomAtIndex(originalAtom2.getIndexNumber())
                            self.IntercellularBonds.append(Bond(None,translatedAtom1, translatedAtom2))

        
        