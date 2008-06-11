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