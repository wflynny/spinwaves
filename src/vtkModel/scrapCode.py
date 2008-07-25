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

        


    def addBond(self, Atom1, Atom2):     
        newBonds = [] #this bond and bonds created by symmetry
        original = Bond(Atom1, Atom2)
        
        #Create Symmetry Bonds
        
        xyz = original.getAtom1().getPosition()
        xyz2 = original.getAtom2().getPosition()
        for symop in self.space_Group.iter_symops():
        # operate on coordinates in non-shifted spacegroup
            pos1 = symop(xyz)
            pos2 = symop(xyz2)
            mask1 = numpy.logical_or(pos1 < 0.0, pos1 >= 1.0)
            translation = numpy.floor(pos1[mask1])  #translates the first atom back to cell at (0,0,0)
            pos1[mask1] -= translation
            pos2[mask1] -= translation  #Uses same translation to translate other atom
            
            #if the second atom is not in the magntic cell, do another translation so that it is
            if pos2[0] < 0 or pos2[1] < 0 or pos2[2] < 0:
                mask2 = pos2 < 0.0
                translation = numpy.floor(pos2[mask2])  #translates the first atom back to cell at (0,0,0)
                print pos1
                print pos2
                pos1[mask2] -= translation
                pos2[mask2] -= translation  #Uses same translation to translate other atom
                atomAtPos1 = self.unit_cell.atomAtPosition(pos2) #Second Position was translated to the first cell (0,0,0)
                atomAtPos2 = self.atomAtPosition(pos1)
            else:
                atomAtPos1 = self.unit_cell.atomAtPosition(pos1) #first Position was translated to the first cell (0,0,0)
                atomAtPos2 = self.atomAtPosition(pos2)
            

#            translation = numpy.floor(pos1)
#            pos1 -= numpy.floor(pos1)
#            pos2 -= numpy.floor(pos1)
            
            
            #Right Now this allows Bonds to be created within a Unit Cell and stored as intercellular
            if atomAtPos2 != None and atomAtPos1 != None:  #Both Atoms could be translated outside the bounds of the Magnetic Cell
#            if atomAtPos2 != None:  #only atom2 could be out of bounds
                newBond = Bond(atomAtPos1, atomAtPos2)
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
                                self.Bonds.append(Bond(translatedAtom1, translatedAtom2))






myblock = CifFile.CifBlock()

cf['a_block'] = myblock

cf['a_block']['_data1'] = "data"
cf['a_block']['_data2'] = 'moredata'

listobj = []
listobj.append("this")
listobj.append("this2")
listobj.append("and this")
cf['a_block']['_listobj'] = listobj



outfile = open(filename, 'w')
outfile.write(cf.WriteOut())
outfile.close()

cf = CifFile.ReadCif(filename)















def __init__(self, parent, id):
        wx.Panel.__init__(self, parent, id)
        
        #Add Space Group
        spaceGroupLabel = wx.StaticText(self, -1, "Space Group:")
        self.spaceGroupSpinner = wx.SpinCtrl(self, -1, "")
        self.spaceGroupSpinner.SetRange(1,230)
        self.spaceGroupSpinner.SetValue(1)
        
        #Add Atom List
        self.atomList = atomListGrid(self, -1)
        
        #Add a, b, c, Alpha, gamma, beta, Na, Nb, Nc
        aLabel = wx.StaticText(self, -1, "a:")
        self.aText = wx.TextCtrl(self, -1, size = (60, -1))
        aSizer = wx.BoxSizer(wx.VERTICAL)
        aSizer.Add(aLabel, 0)
        aSizer.Add(self.aText, 0)
        
        bLabel = wx.StaticText(self, -1, "b:")
        self.bText = wx.TextCtrl(self, -1, size = (60, -1))
        bSizer = wx.BoxSizer(wx.VERTICAL)
        bSizer.Add(bLabel, 0)
        bSizer.Add(self.bText, 0)
        
        cLabel = wx.StaticText(self, -1, "c:")
        self.cText = wx.TextCtrl(self, -1, size = (60, -1))
        cSizer = wx.BoxSizer(wx.VERTICAL)
        cSizer.Add(cLabel, 0)
        cSizer.Add(self.cText, 0)
        
        alphaLabel = wx.StaticText(self, -1, "alpha:")
        self.alphaText = wx.TextCtrl(self, -1, size = (60, -1))
        alphaSizer = wx.BoxSizer(wx.VERTICAL)
        alphaSizer.Add(alphaLabel, 0)
        alphaSizer.Add(self.alphaText, 0)
        
        betaLabel = wx.StaticText(self, -1, "beta:")
        self.betaText = wx.TextCtrl(self, -1, size = (60, -1))
        betaSizer = wx.BoxSizer(wx.VERTICAL)
        betaSizer.Add(betaLabel, 0)
        betaSizer.Add(self.betaText, 0)
        
        gammaLabel = wx.StaticText(self, -1, "gamma:")
        self.gammaText = wx.TextCtrl(self, -1, size = (60, -1))
        gammaSizer = wx.BoxSizer(wx.VERTICAL)
        gammaSizer.Add(gammaLabel, 0)
        gammaSizer.Add(self.gammaText, 0)
        
        
        naLabel = wx.StaticText(self, -1, "Na:")
        self.naText = wx.TextCtrl(self, -1, size = (60, -1))
        naSizer = wx.BoxSizer(wx.VERTICAL)
        naSizer.Add(naLabel, 0)
        naSizer.Add(self.naText, 0)
        
        nbLabel = wx.StaticText(self, -1, "Nb:")
        self.nbText = wx.TextCtrl(self, -1, size = (60, -1))
        nbSizer = wx.BoxSizer(wx.VERTICAL)
        nbSizer.Add(nbLabel, 0)
        nbSizer.Add(self.nbText, 0)
        
        ncLabel = wx.StaticText(self, -1, "Nc:")
        self.ncText = wx.TextCtrl(self, -1, size = (60, -1))
        ncSizer = wx.BoxSizer(wx.VERTICAL)
        ncSizer.Add(ncLabel, 0)
        ncSizer.Add(self.ncText, 0)
        
        #Organize a,b,c and alpha, gamma, beta , Na, Nb, Nc into a grid
        dimSizer = wx.GridSizer(cols = 3, hgap = 15, vgap = 5)
        dimSizer.Add(aSizer)
        dimSizer.Add(bSizer)
        dimSizer.Add(cSizer)
        dimSizer.Add(alphaSizer)
        dimSizer.Add(betaSizer)
        dimSizer.Add(gammaSizer)
        dimSizer.Add(naSizer)
        dimSizer.Add(nbSizer)
        dimSizer.Add(ncSizer)
        
        sizer = wx.BoxSizer(wx.VERTICAL)
        
        spaceGroupSizer = wx.BoxSizer(wx.HORIZONTAL)
        spaceGroupSizer.Add(spaceGroupLabel, 0)
        spaceGroupSizer.Add(self.spaceGroupSpinner, 0)
        
        #Add a button on upper right to generate new image
        self.genButton = wx.Button(self, -1, "Generate")
        spaceGroupSizer.Add(self.genButton, flag = wx.ALIGN_RIGHT)
        
        sizer.Add(spaceGroupSizer, 0)
        sizer.Add(dimSizer, 0)
        sizer.Add(self.atomList, 1, wx.EXPAND)
        
        self.SetSizer(sizer)
        
        #For now, the generate button will create a new magnetic cell
        #which will then be drawn by the vtkPanel
        self.Bind(wx.EVT_BUTTON, self.OnGenerate, self.genButton)



















#-----------------------------------------------------------------------
    #Code from william being modified

#---------------------------------------------------------------------------



keyMap = {
    wx.WXK_BACK : "WXK_BACK",
    wx.WXK_TAB : "WXK_TAB",
    wx.WXK_RETURN : "WXK_RETURN",
    wx.WXK_ESCAPE : "WXK_ESCAPE",
    wx.WXK_SPACE : "WXK_SPACE",
    wx.WXK_DELETE : "WXK_DELETE",
    wx.WXK_START : "WXK_START",
    wx.WXK_LBUTTON : "WXK_LBUTTON",
    wx.WXK_RBUTTON : "WXK_RBUTTON",
    wx.WXK_CANCEL : "WXK_CANCEL",
    wx.WXK_MBUTTON : "WXK_MBUTTON",
    wx.WXK_CLEAR : "WXK_CLEAR",
    wx.WXK_SHIFT : "WXK_SHIFT",
    wx.WXK_ALT : "WXK_ALT",
    wx.WXK_CONTROL : "WXK_CONTROL",
    wx.WXK_MENU : "WXK_MENU",
    wx.WXK_PAUSE : "WXK_PAUSE",
    wx.WXK_CAPITAL : "WXK_CAPITAL",
    wx.WXK_PRIOR : "WXK_PRIOR",
    wx.WXK_NEXT : "WXK_NEXT",
    wx.WXK_END : "WXK_END",
    wx.WXK_HOME : "WXK_HOME",
    wx.WXK_LEFT : "WXK_LEFT",
    wx.WXK_UP : "WXK_UP",
    wx.WXK_RIGHT : "WXK_RIGHT",
    wx.WXK_DOWN : "WXK_DOWN",
    wx.WXK_SELECT : "WXK_SELECT",
    wx.WXK_PRINT : "WXK_PRINT",
    wx.WXK_EXECUTE : "WXK_EXECUTE",
    wx.WXK_SNAPSHOT : "WXK_SNAPSHOT",
    wx.WXK_INSERT : "WXK_INSERT",
    wx.WXK_HELP : "WXK_HELP",
    wx.WXK_NUMPAD0 : "WXK_NUMPAD0",
    wx.WXK_NUMPAD1 : "WXK_NUMPAD1",
    wx.WXK_NUMPAD2 : "WXK_NUMPAD2",
    wx.WXK_NUMPAD3 : "WXK_NUMPAD3",
    wx.WXK_NUMPAD4 : "WXK_NUMPAD4",
    wx.WXK_NUMPAD5 : "WXK_NUMPAD5",
    wx.WXK_NUMPAD6 : "WXK_NUMPAD6",
    wx.WXK_NUMPAD7 : "WXK_NUMPAD7",
    wx.WXK_NUMPAD8 : "WXK_NUMPAD8",
    wx.WXK_NUMPAD9 : "WXK_NUMPAD9",
    wx.WXK_MULTIPLY : "WXK_MULTIPLY",
    wx.WXK_ADD : "WXK_ADD",
    wx.WXK_SEPARATOR : "WXK_SEPARATOR",
    wx.WXK_SUBTRACT : "WXK_SUBTRACT",
    wx.WXK_DECIMAL : "WXK_DECIMAL",
    wx.WXK_DIVIDE : "WXK_DIVIDE",
    wx.WXK_F1 : "WXK_F1",
    wx.WXK_F2 : "WXK_F2",
    wx.WXK_F3 : "WXK_F3",
    wx.WXK_F4 : "WXK_F4",
    wx.WXK_F5 : "WXK_F5",
    wx.WXK_F6 : "WXK_F6",
    wx.WXK_F7 : "WXK_F7",
    wx.WXK_F8 : "WXK_F8",
    wx.WXK_F9 : "WXK_F9",
    wx.WXK_F10 : "WXK_F10",
    wx.WXK_F11 : "WXK_F11",
    wx.WXK_F12 : "WXK_F12",
    wx.WXK_F13 : "WXK_F13",
    wx.WXK_F14 : "WXK_F14",
    wx.WXK_F15 : "WXK_F15",
    wx.WXK_F16 : "WXK_F16",
    wx.WXK_F17 : "WXK_F17",
    wx.WXK_F18 : "WXK_F18",
    wx.WXK_F19 : "WXK_F19",
    wx.WXK_F20 : "WXK_F20",
    wx.WXK_F21 : "WXK_F21",
    wx.WXK_F22 : "WXK_F22",
    wx.WXK_F23 : "WXK_F23",
    wx.WXK_F24 : "WXK_F24",
    wx.WXK_NUMLOCK : "WXK_NUMLOCK",
    wx.WXK_SCROLL : "WXK_SCROLL",
    wx.WXK_PAGEUP : "WXK_PAGEUP",
    wx.WXK_PAGEDOWN : "WXK_PAGEDOWN",
    wx.WXK_NUMPAD_SPACE : "WXK_NUMPAD_SPACE",
    wx.WXK_NUMPAD_TAB : "WXK_NUMPAD_TAB",
    wx.WXK_NUMPAD_ENTER : "WXK_NUMPAD_ENTER",
    wx.WXK_NUMPAD_F1 : "WXK_NUMPAD_F1",
    wx.WXK_NUMPAD_F2 : "WXK_NUMPAD_F2",
    wx.WXK_NUMPAD_F3 : "WXK_NUMPAD_F3",
    wx.WXK_NUMPAD_F4 : "WXK_NUMPAD_F4",
    wx.WXK_NUMPAD_HOME : "WXK_NUMPAD_HOME",
    wx.WXK_NUMPAD_LEFT : "WXK_NUMPAD_LEFT",
    wx.WXK_NUMPAD_UP : "WXK_NUMPAD_UP",
    wx.WXK_NUMPAD_RIGHT : "WXK_NUMPAD_RIGHT",
    wx.WXK_NUMPAD_DOWN : "WXK_NUMPAD_DOWN",
    wx.WXK_NUMPAD_PRIOR : "WXK_NUMPAD_PRIOR",
    wx.WXK_NUMPAD_PAGEUP : "WXK_NUMPAD_PAGEUP",
    wx.WXK_NUMPAD_NEXT : "WXK_NUMPAD_NEXT",
    wx.WXK_NUMPAD_PAGEDOWN : "WXK_NUMPAD_PAGEDOWN",
    wx.WXK_NUMPAD_END : "WXK_NUMPAD_END",
    wx.WXK_NUMPAD_BEGIN : "WXK_NUMPAD_BEGIN",
    wx.WXK_NUMPAD_INSERT : "WXK_NUMPAD_INSERT",
    wx.WXK_NUMPAD_DELETE : "WXK_NUMPAD_DELETE",
    wx.WXK_NUMPAD_EQUAL : "WXK_NUMPAD_EQUAL",
    wx.WXK_NUMPAD_MULTIPLY : "WXK_NUMPAD_MULTIPLY",
    wx.WXK_NUMPAD_ADD : "WXK_NUMPAD_ADD",
    wx.WXK_NUMPAD_SEPARATOR : "WXK_NUMPAD_SEPARATOR",
    wx.WXK_NUMPAD_SUBTRACT : "WXK_NUMPAD_SUBTRACT",
    wx.WXK_NUMPAD_DECIMAL : "WXK_NUMPAD_DECIMAL",
    wx.WXK_NUMPAD_DIVIDE : "WXK_NUMPAD_DIVIDE"
    }



class AtomTable(wx.grid.PyGridTableBase):
    def __init__(self):
        wx.grid.PyGridTableBase.__init__(self)
        self.colLabels = ['Name', 'Atomic Number','x', 'y','z']
        self.rowLabels=['Atom 1']

        self.dataTypes = [wx.grid.GRID_VALUE_STRING, #Name
                          wx.grid.GRID_VALUE_FLOAT,#Atomic Number
                          wx.grid.GRID_VALUE_FLOAT,#x
                          wx.grid.GRID_VALUE_FLOAT, #y
                          wx.grid.GRID_VALUE_FLOAT, #z
                          ]
        self.data = []
        self.data.append(['', #Name
                        '', #Atomic #
                        '', #x
                        '', #y
                        '', #z
                        ])
        return

    #--------------------------------------------------
    # required methods for the wxPyGridTableBase interface
    def GetNumberRows(self):
        return len(self.data)
        #return len(self.data)
    def GetNumberCols(self):
        return len(self.colLabels)
    def IsEmptyCell(self, row, col):
        try:
            return not self.data[row][col]
        except IndexError:
            return True
    # Get/Set values in the table.  The Python version of these
    # methods can handle any data-type, (as long as the Editor and
    # Renderer understands the type too,) not just strings as in the
    # C++ version.
    def GetValue(self, row, col):
        try:
            return self.data[row][col]
        except IndexError:
            return ''   #ME - Should this be '' or None for ints?
    def SetValue(self, row, col, value):
        try:
            self.data[row][col] = value
            #print 'SetValue works',self.GetNumberRows(),self.data[row][1]
        except IndexError:
            # add a new row
            #print 'IndexError in SetValue',self.GetNumberRows()
            self.AppendRow()
            self.data[row][col]=value
            #print 'IndexError in SetValue after SetValue',self.GetNumberRows()
            #print 'setting row ',row,' col ',col, ' val ',value
            #print self.__dict__
            #self.SetValue(row, col, value)
        return

    def AppendRow(self):
            self.data.append([''] * self.GetNumberCols())
            #print 'After Append SetValue',self.GetNumberRows()
            #self.rowLabels[row]='File '+str(len(self.rowLabels))
            #self.rowLabels.append('File '+str(len(self.rowLabels)))

            # tell the grid we've added a row
            msg = wx.grid.GridTableMessage(self,            # The table
                    wx.grid.GRIDTABLE_NOTIFY_ROWS_APPENDED, # what we did to it
                    1                                       # how many
                    )
            #print 'size notified',self.GetNumberRows()
            self.GetView().ProcessTableMessage(msg)
            #print 'self.rowLabels', self.rowLabels
            #self.data[row][col] = value


    #--------------------------------------------------
    # Some optional methods
    # Called when the grid needs to display labels
    def GetColLabelValue(self, col):
        return self.colLabels[col]
    # Called when the grid needs to display labels
    def GetRowLabelValue(self, row):
        return 'Atom '+str(row)
        #return self.rowLabels[row]
    # Called to determine the kind of editor/renderer to use by
    # default, doesn't necessarily have to be the same type used
    # natively by the editor/renderer if they know how to convert.

    def GetTypeName(self, row, col):
        return self.dataTypes[col]
    # Called to determine how the data can be fetched and stored by the
    # editor and renderer.  This allows you to enforce some type-safety
    # in the grid.
    def CanGetValueAs(self, row, col, typeName):
        colType = self.dataTypes[col].split(':')[0]
        if typeName == colType:
            return True
        else:
            return False
    def CanSetValueAs(self, row, col, typeName):
        return self.CanGetValueAs(row, col, typeName)

    def DeleteRows(self,pos=0,numRows=1):
#        print 'Delete number',self.GetNumberRows()
#        print 'pos',pos
#        print 'numRows', numRows
        if numRows>=0 and numRows<=self.GetNumberRows():
#            print 'Delete',numRows
            #for i in range(numRows):
            #    self.data.pop()
            del self.data[pos:pos+numRows]
            msg = wx.grid.GridTableMessage(self,            # The table
            wx.grid.GRIDTABLE_NOTIFY_ROWS_DELETED, # what we did to it
            pos,numRows                                     # how many
            )
            #msg = wx.grid.GridTableMessage(self, 0, numRows)
            self.GetView().ProcessTableMessage(msg)
            
#            print 'Deleted'
            self.UpdateValues()
            return True
        else:
            return False
        
    def UpdateValues( self ):
            """Update all displayed values"""
            msg =wx.grid.GridTableMessage(self, wx.grid.GRIDTABLE_REQUEST_VIEW_GET_VALUES)
            self.GetView().ProcessTableMessage(msg)
#---------------------------------------------------------------------------
class AtomTableGrid(wx.grid.Grid):
    def __init__(self, parent):
        wx.grid.Grid.__init__(self, parent, -1)
        table = AtomTable()

#        gridbar.register(self) Don't know what this is



        # The second parameter means that the grid is to take ownership of the
        # table and will destroy it when done.  Otherwise you would need to keep
        # a reference to it and call it's Destroy method later.
        self.SetTable(table, True)
        #attr = wx.grid.GridCellAttr()
        #attr.SetReadOnly(True)
        #attr.SetRenderer(gridbar.GridCellBarRenderer())
        #self.SetColAttr(13, attr)
        #self.SetCellValue(1,13,'q')
        #self.SetCellRenderer(1,13,gridbar.GridCellBarRenderer)
        #self.SetRowLabelSize(0)
        self.SetMargins(0,0)
        self.AutoSize()
        #wx.grid.Grid.SetSelectionMode(self,wx.grid.Grid.SelectRows)
        wx.grid.Grid.EnableEditing(self,False)
        attr=wx.grid.GridCellAttr()
        attr.SetReadOnly(True)
        self.SetColAttr(0,attr)
        for col in range(1,14):
            attr=wx.grid.GridCellAttr()
            attr.SetReadOnly(True)
            #attr.SetBackgroundColour('grey' if col%2 else (139, 139, 122))
            #attr.SetTextColour((167,167,122) if col%2 else (139, 139, 122))
            self.SetColAttr(col,attr)
        #wx.grid.EVT_GRID_CELL_LEFT_DCLICK(self, self.OnLeftDClick)
        wx.grid.EVT_GRID_CELL_LEFT_CLICK(self,self.OnLeftClick)
        #wx.grid.EVT_GRID_CELL_CHANGE(self,self.OnCellChange)
        wx.grid.EVT_GRID_LABEL_LEFT_DCLICK(self,self.onLeftDClickRowCell)

    # I do this because I don't like the default behaviour of not starting the
    # cell editor on double clicks, but only a second click.
    def OnLeftClick(self, evt):
        print 'LeftClick'
        col=evt.GetCol()
        row=evt.GetRow()
        table=self.GetTable()
        if col<=0 and row >=0:
            currval=table.GetValue(row,0)
            if currval=='':
                table.SetValue(row,0,'x')
            else:
                table.SetValue(row,0,'')


        #if self.CanEnableCellControl():
        #    self.EnableCellEditControl()
        wx.grid.Grid.ForceRefresh(self)

    def OnCellChange(self, evt):
#        print 'Changed'
        if self.CanEnableCellControl():
            self.EnableCellEditControl()
        wx.grid.Grid.ForceRefresh(self)
        evt.Skip()



    def onLeftDClickRowCell(self,evt):
        col=evt.GetCol()
        table=self.GetTable()
        data=N.array(table.data)
#        print 'before ', data[:,0]
        col_to_sort=[(i,s) for i,s in enumerate(data[:,col])]
        col_to_sort.sort(lambda x,y: cmp(x[1],y[1]))
        g_col = [i for (i,s) in col_to_sort]
        #print col_to_sort
        if col >=0:
            if (N.diff(g_col)>0).all():
                g_col=g_col[::-1]

            #print 'col=',col
            #print 'sort '
            #print g
            for i in range(data.shape[1]):
                data[:,i]=data[g_col,i]
            table.data=data.tolist()
#            print 'after',data[:,0]
            wx.grid.Grid.AutoSize(self)
            wx.grid.Grid.ForceRefresh(self)
        #evt.Skip()


class AtomPanel(wx.Panel):
    ## Internal name for the AUI manager
    window_name = "AtomPanel"
    ## Title to appear on top of the window
    window_caption = "File Catalog Panel"
    CENTER_PANE = True

    def __init__(self,parent,id):
        wx.Panel.__init__(self,parent,id,style=0)
        cfstr='Spinwaves'
        #
        
        #print 'Global Config',wx.CONFIG_USE_GLOBAL_FILE
        wx.GetApp().SetAppName(cfstr)
        wx.GetApp().SetVendorName(cfstr)
        self.config=wx.ConfigBase.Get()
        #self.config=wx.Config(cfstr)
        #self.config=wx.Config(cfstr,style=wx.CONFIG_USE_GLOBAL_FILE)
        #self.config=wx.Config(cfstr,style=wx.CONFIG_USE_LOCAL_FILE)
        self.config.SetAppName(cfstr)
        self.config.SetVendorName(cfstr)
        self.config.Flush()
#        print 'mypath',self.config.GetPath().encode('ascii')
#        print 'myApp',self.config.GetAppName()
#        print 'myVendor',self.config.GetVendorName(),'internal',wx.GetApp().GetVendorName()
#        print 'mypath',self.config.Read('mypath')
        grid = AtomTableGrid(self)
        bs = wx.BoxSizer(wx.VERTICAL)
        bs.Add(grid, 1, wx.GROW|wx.ALL|wx.EXPAND, 5)
        self.SetSizer(bs)
        self.parent=parent
        self.grid=grid
        self.bs=bs
        self.tooltip = ''
        self.catalog=None
        self.log=sys.stdout
        self.grid.GetGridWindow().Bind(wx.EVT_MOTION, self.onMouseOver)
        self.grid.GetGridWindow().Bind(wx.EVT_CHAR,self.OnChar)
        #self.filetree_frame=parent.filetree_frame
        #self.filetree_panel=parent.filetree_frame.filetree_panel
        self.filetree_panel = None

        self.popupmenu = wx.Menu()
        item1 = self.popupmenu.Append(wx.ID_ANY, "Send Group for Reduction")
        self.popupmenu.AppendSeparator()
        item2 = self.popupmenu.Append(wx.ID_ANY, "Clear Selections")
        self.grid.Bind(wx.EVT_MENU, self.OnSendGroup, item1)
        self.grid.Bind(wx.EVT_MENU, self.OnClearGroups, item2)
        #self.grid.GetGridWindow().Bind(wx.EVT_CONTEXT_MENU,self.OnShowPopup)
        self.grid.GetGridWindow().Bind(wx.EVT_RIGHT_UP,self.OnShowPopup)

    def OnShowPopup(self,event):
#        print 'popping'
        pos=event.GetPosition()
        #pos=self.ScreenToClient(pos)
        self.grid.PopupMenu(self.popupmenu,pos)
        event.Skip()

    
    def OnSendGroup(self,event):
#        print 'SendingGroup'
        self.SendGroupToTree()
        
    def OnClearGroups(self,event):
#        print 'Clearing Group'
        table=self.grid.GetTable()
        nrows=table.GetNumberRows()
        #print 'old_nrows',old_nrows
        #for row in range(len(self.catalog.files)):
        for row in range(nrows):
            table.SetValue(row,0,'') # selected=False
        self.grid.ForceRefresh()
        

    def onMouseOver(self, event):
        '''
        Method to calculate where the mouse is pointing and
        then set the tooltip dynamically.
        '''

        # Use CalcUnscrolledPosition() to get the mouse position within the
        # entire grid including what's offscreen
        x, y =self.grid.CalcUnscrolledPosition(event.GetX(),event.GetY())

        coords = self.grid.XYToCell(x, y)
        #coords = grid.XYToCell(x, y)
        col = coords[1]
        table=self.grid.GetTable()
        # Example colum limit to apply the custom tooltip to
        if col>4:
            row = coords[0]
            bar = table.GetValue(row, col)
            try:
                low=bar.low
                high=bar.high
                event.GetEventObject().SetToolTipString('range=(%4.3f,%4.3f)'%(low,high))
                self.tooltip='range=(%4.3f,%4.3f)' %(low,high)
            except:
                event.GetEventObject().SetToolTipString('')
                self.tooltip = ''
        else:
            event.GetEventObject().SetToolTipString('')
            self.tooltip = ''

    def SendGroupToTree(self):
        if self.catalog!=None:
            table=self.grid.GetTable()
            nrows=table.GetNumberRows()
            ncols=table.GetNumberCols()
            sequence_selected=[]
            files_selected=[]
            pol_states=[]
            count_types=[]
            fileseq_orig=N.array(self.catalog.fileseq)
            treenode_data={}
            treedata=[]
            treedict={}
            for row in range(nrows):
                gridval=table.GetValue(row,0) #get selected rows
                if gridval=='x':
                    treenode_data={}
                    filename=table.GetValue(row,1)
                    sequence=table.GetValue(row,2)
                    polstate=table.GetValue(row,3)
                    #print 'file selected',filename
                    sequence_selected.append(sequence)
                    files_selected.append(filename)
                    pol_states.append(polstate)
                    loc=N.where(fileseq_orig==sequence)[0]
                    currdata=self.catalog.data[loc]['full_data']
                    count_type=currdata.metadata['count_info']['count_type']
                    count_types.append(count_type)
                    MonitorCorrect=0
                    PolMonitorCorrect=1
                    if count_type=='time':
                        PolMonitorCorrect=0
                        MonitorCorrect=0
                    treenode_data['data']=currdata
                    treenode_data['filename']=filename
                    treenode_data['absolute_filename']=os.path.join(wx.Config().GetPath(),filename)
                    treenode_data['sequence']=sequence
                    treenode_data['polstate']=polstate
                    #treenode_data['flags']=polcorrect.PBflags()
                    treedata.append(treenode_data)
                    treedict[polstate]=treenode_data
                    #print 'loc',loc,treenode_data['filename']

            #print 'selected files', files_selected
            #print 'sequences',sequence_selected
            if len(sequence_selected)>0:
                CurrentGroup=self.AddGroup(PolMonitorCorrect=PolMonitorCorrect,MonitorCorrect=MonitorCorrect)
                tree=self.filetree_panel.tree
                for curnode in treedata:
                    #print 'curnode',curnode['filename']
                    self.AddItem(CurrentGroup,curnode)
                tree.SelectItem(self.DataGroup)
                tree.Expand(self.DataGroup)
                tree.SelectItem(CurrentGroup)
                tree.Expand(CurrentGroup)



    def OnChar(self, event):

        keycode = event.GetKeyCode()
        keyname = keyMap.get(keycode, None)

        #if keycode == wx.WXK_BACK:
        #    self.log.write("OnKeyDown: HAHAHAHA! I Vetoed Your Backspace! HAHAHAHA\n")
        #    return

        if keyname is None:
            if "unicode" in wx.PlatformInfo:
                keycode = event.GetUnicodeKey()
                if keycode <= 127:
                    keycode = event.GetKeyCode()
                keyname = "\"" + unichr(event.GetUnicodeKey()) + "\""
                if keycode < 27:
                    keyname = "Ctrl-%s" % chr(ord('A') + keycode-1)

            elif keycode < 256:
                if keycode == 0:
                    keyname = "NUL"
                elif keycode < 27:
                    keyname = "Ctrl-%s" % chr(ord('A') + keycode-1)
                else:
                    keyname = "\"%s\"" % chr(keycode)
            else:
                keyname = "unknown (%s)" % keycode
        #keyname=self.keybuff+keyname
        if keyname=='Ctrl-R' and self.catalog!=None:
            self.SendGroupToTree()

        event.Skip()

    def AddGroup(self,MonitorCorrect=0,PolMonitorCorrect=1):
        tree=self.filetree_panel.tree
        root=tree.root
        DataGroup=tree.GetFirstChild(root)[0]
        self.DataGroup=DataGroup
        #print 'DataGroup',DataGroup
        if tree.HasChildren(DataGroup):
        #if 0:
            group_id=tree.GetChildrenCount(DataGroup, recursively=False)
        else:
            group_id=0
        child = tree.AppendItem(DataGroup, "Group"+str(group_id))
        tree.SetItemBold(child, True)
        pbflags=polcorrect.PBflags()
        pbflags.MonitorCorrect=MonitorCorrect
        pbflags.PolMonitorCorrect=PolMonitorCorrect
        pbflags.MonoSelect=1
        pbflags.Debug=0
        pbflags.SimFlux=0
        pbflags.SimDeviate=0
        pbflags.NoNegativeCS=0
        pbflags.HalfPolarized=0
        pbflags.CountsEnable=[0,0,0,0]
        pbflags.CountsAdd1=[0,0,0,0]
        pbflags.CountsAdd2=[0,0,0,0]
        pbflags.Sconstrain=[0,0,0,0]
        pbflags.Spp=[0,0,0,0]
        pbflags.Smm=[0,0,0,0]
        pbflags.Spm=[0,0,0,0]
        pbflags.Smp=[0,0,0,0]
        groupdata={}
        groupdata['pbflags']=pbflags
        groupdata['cellfile']=''
        #groupdata['absolute_files']=self.files
        tree.SetPyData(child, groupdata)
        tree.SetItemImage(child, 24, CT.TreeItemIcon_Normal)
        tree.SetItemImage(child, 13, CT.TreeItemIcon_Expanded)
        return child

    def AddItem(self,CurrentGroup,data):
        tree=self.filetree_panel.tree
        #CurrentGroup=tree.GetLastChild(self.DataGroup)[0]
        child = tree.AppendItem(CurrentGroup,data['filename'])#,ct_type=1) #ct_type=1 gives check box
        tree.SetItemBold(child, True)
        tree.SetPyData(child, data)
        tree.SetItemImage(child, 24, CT.TreeItemIcon_Normal)
        tree.SetItemImage(child, 13, CT.TreeItemIcon_Expanded)


    def UpdateCatalog(self):
            table=self.grid.GetTable()
            old_nrows=table.GetNumberRows()
            print 'old_nrows',old_nrows
            if old_nrows >0:
                print 'numRows before Deletion',old_nrows
                wx.grid.Grid.DeleteRows(self.grid,pos=0,numRows=old_nrows,updateLabels=True)
            for row in range(len(self.catalog.files)):
                #print 'row',row
                table.SetValue(row,0,'') # selected=False
                table.SetValue(row,1,self.catalog.files[row])
                table.SetValue(row,2,str(self.catalog.data[row]['fileseq_number']))
                table.SetValue(row,3,str(self.catalog.data[row]['polarization state']))
                table.SetValue(row,4,str(self.catalog.data[row]['hsample']))
                table.SetValue(row,5,str(self.catalog.data[row]['vsample']))
                i=6
                if self.catalog.data[row].has_key('h'):
                    range_column=(self.catalog.h_range.min,self.catalog.h_range.max)
                    range_cell=(self.catalog.data[row]['h']['min'],self.catalog.data[row]['h']['max'])
                    currbar=gridbar.Bar(self.catalog.data[row]['h']['min'],self.catalog.data[row]['h']['max'],range_column)
                    table.SetValue(row,i,currbar); #
                i=i+1
                if self.catalog.data[row].has_key('k'):
                    range_column=(self.catalog.k_range.min,self.catalog.k_range.max)
                    range_cell=(self.catalog.data[row]['k']['min'],self.catalog.data[row]['k']['max'])
                    currbar=gridbar.Bar(self.catalog.data[row]['k']['min'],self.catalog.data[row]['k']['max'],range_column)
                    table.SetValue(row,i,currbar); #
                i=i+1
                if self.catalog.data[row].has_key('l'):
                    range_column=(self.catalog.l_range.min,self.catalog.l_range.max)
                    range_cell=(self.catalog.data[row]['l']['min'],self.catalog.data[row]['l']['max'])
                    currbar=gridbar.Bar(self.catalog.data[row]['l']['min'],self.catalog.data[row]['l']['max'],range_column)
                    table.SetValue(row,i,currbar); #
                i=i+1
                if self.catalog.data[row].has_key('e'):
                    range_column=(self.catalog.e_range.min,self.catalog.e_range.max)
                    range_cell=(self.catalog.data[row]['e']['min'],self.catalog.data[row]['e']['max'])
                    currbar=gridbar.Bar(self.catalog.data[row]['e']['min'],self.catalog.data[row]['e']['max'],range_column)
                    table.SetValue(row,i,currbar);
                i=i+1
                if self.catalog.data[row].has_key('a3'):
                    range_column=(self.catalog.a3_range.min,self.catalog.a3_range.max)
                    range_cell=(self.catalog.data[row]['a3']['min'],self.catalog.data[row]['a3']['max'])
                    currbar=gridbar.Bar(self.catalog.data[row]['a3']['min'],self.catalog.data[row]['a3']['max'],range_column)
                    table.SetValue(row,i,currbar);
                i=i+1
                if self.catalog.data[row].has_key('a4'):
                    range_column=(self.catalog.a4_range.min,self.catalog.a4_range.max)
                    range_cell=(self.catalog.data[row]['a4']['min'],self.catalog.data[row]['a4']['max'])
                    currbar=gridbar.Bar(self.catalog.data[row]['a4']['min'],self.catalog.data[row]['a4']['max'],range_column)
                    table.SetValue(row,i,currbar);
                i=i+1
                if self.catalog.data[row].has_key('temp'):
                    #print 'temp'
                    range_column=(self.catalog.temp_range.min,self.catalog.temp_range.max)
                    range_cell=(self.catalog.data[row]['temp']['min'],self.catalog.data[row]['temp']['max'])
                    currbar=gridbar.Bar(self.catalog.data[row]['temp']['min'],self.catalog.data[row]['temp']['max'],range_column)
                    table.SetValue(row,i,currbar);
                i=i+1
                if self.catalog.data[row].has_key('magfield'):
                    range_column=(self.catalog.magfield_range.min,self.catalog.magfield_range.max)
                    range_cell=(self.catalog.data[row]['magfield']['min'],self.catalog.data[row]['magfield']['max'])
                    currbar=gridbar.Bar(self.catalog.data[row]['magfield']['min'],self.catalog.data[row]['magfield']['max'],range_column)
                    table.SetValue(row,i,currbar);
                i=i+1
            wx.grid.Grid.AutoSize(self.grid)
            wx.grid.Grid.ForceRefresh(self.grid)

    def OnOpen(self,event):
        # Create the dialog. In this case the current directory is forced as the starting
        # directory for the dialog, and no default file name is forced. This can easilly
        # be changed in your program. This is an 'open' dialog, and allows multitple
        # file selections as well.
        #
        # Finally, if the directory is changed in the process of getting files, this
        # dialog is set up to change the current working directory to the path chosen.

        #defaultDir=os.getcwd()
        #defaultDir=r'C:\polcorrecter\data'
        #defaultDir=wx.Config().GetPath()
        defaultDir=self.config.GetPath().encode('ascii')
        defaultDir=self.config.Read('mypath').encode('ascii')
        print 'defaultDir', defaultDir
        wildcard="bt7 files (*.bt7)|*.bt7|All files (*.*)|*.*"
        dlg = wx.FileDialog(
            self, message="Choose a file",
            defaultDir=defaultDir,
            defaultFile="",
            wildcard=wildcard,
            style=wx.OPEN | wx.MULTIPLE | wx.CHANGE_DIR
            )

        # Show the dialog and retrieve the user response. If it is the OK response,
        # process the data.
        if dlg.ShowModal() == wx.ID_OK:
            # This returns a Python list of files that were selected.
            paths = dlg.GetPaths()
            file0=paths[0]
#            print 'file0',file0
            cwd=os.path.dirname(file0)
#            print 'cwd',cwd.encode('ascii')
            self.config.SetPath(cwd.encode('utf8'))
#            print 'config style', self.config.GetStyle()
            self.config.Write('mypath',cwd.encode('utf8'))
            self.config.Flush()
            self.files=paths
#            print 'Opening', self.config.GetPath()
#            print 'mypath', self.config.Read('mypath')
            evt=myEVT_CLEAR_TREE(self.GetId())
            wx.PostEvent(self.filetree_panel.tree , evt)  #I'm not sure if this or the other is cleaner...  
            #self.filetree_panel.tree.GetEventHandler().ProcessEvent(evt)  
            self.catalog=classify_files.readfiles(self.files)
#            print 'event posted'
            #evt=ClearTreeEvent(myEVT_CLEAR_TREE,self.GetId())
            #self.GetEventHandler().ProcessEvent(evt)
            self.UpdateCatalog()      
            wx.grid.Grid.ForceRefresh(self.grid)  
#            print 'updated'
        # Destroy the dialog. Don't do this until you are done with it!
        # BAD things can happen otherwise!
        dlg.Destroy()


class AtomFrame(wx.Frame):
    def __init__(self,parent,id,log=None):
        wx.Frame.__init__(self,parent,id,'File Catalog',size=(640,200))
        self.Bind(wx.EVT_CLOSE,self.OnCloseWindow)
#        self.filetree_frame=FTC.FileTreeFrame(self,-1)
        self.atom_panel = AtomPanel(self,-1)
        self.atom_panel.log=log
#        self.filetree_frame.Show()


    def OnCloseWindow(self,event):
        self.Destroy()





        interCutoffBonds = []
        for eachBond in simpleBonds.list:
#            print "atom1: ", eachBond.pos1, "atom2: ", eachBond.pos2
            #Check if one bond is in the first cutoff cell
            if eachBond.pos1[0] < Na or eachBond.pos2[0] < Na: #x
                if eachBond.pos1[1] < Nb or eachBond.pos2[1] < Nb: #y
                    if eachBond.pos1[2] < Nc or eachBond.pos2[2] < Nc: #z
                        #check if the second bond is not in the first cutoff cell
                        if (not eachBond.pos1[0] < Na) or (not eachBond.pos2[0] < Na): #x
                            if (not eachBond.pos1[1] < Nb) or (not eachBond.pos2[1] < Nb): #y
                                if (not eachBond.pos1[2] < Nc) or (not eachBond.pos2[2] < Nc): #z
                                    interCutoffBonds.append(eachBond)
                                    
                                    

   
    def exportForMonteCarlo(self, filename):
        size = 5
        
        timer = Timer()
        
        file = open(filename, 'w')
   
        class SimpleBond():
            def __init__(self, pos1, pos2, jMatrix):
                self.pos1 = pos1
                self.pos2 = pos2
                self.jMatrix = jMatrix
                
            def sameBond(self, bond2):
                if self.pos1 == bond2.pos1 or self.pos1 == bond2.pos2:
                    if self.pos2 == bond2.pos2 or self.pos2 == bond2.pos1:
                        return True
                return False
        
        Na = self.getCutoffCell().getNa()
        Nb = self.getCutoffCell().getNb()
        Nc = self.getCutoffCell().getNc()
        
        class SimpleAtom():
            def __init__(self, pos, index):
                self.pos = pos
                self.index = index
                self.cellPos = []
                self.cellPosX = int(pos[0])/Na
                self.cellPosY = int(pos[1])/Nb
                self.cellPosZ = int(pos[2])/Nc
                self.interactions = []
                #self.interactions[position of other atom] = j number
                
            #might want to change this to position later when all atoms wont be created in same list
            def addInteraction(self, atom2, jMat):
#                self.interactions[atom2] = jMat
                self.interactions.append([atom2, jMat])
            
            #comparisons based on positions
            def __lt__(self, other):
                self.otherCellPosX = int(other.pos[0])/Na
                self.otherCellPosY = int(other.pos[1])/Nb
                self.otherCellPosZ = int(other.pos[2])/Nc
                
                if otherCellPosX > self.cellPosX:
                    return False
                if otherCellPosX < self.cellPosX:
                    return True
                #X's equal
                if otherCellPosY > self.cellPosY:
                    return False
                if otherCellPosY < self.cellPosY:
                    return True
                #Y's equal
                if otherCellPosZ > self.cellPosZ:
                    return False
                if otherCellPosZ < self.cellPosZ:
                    return True
                
                #They are in the same cell
                if other.pos[2] > self.pos[2]:
                    return False
                if other.pos[2] < self.pos[2]:
                    return True
                #Z's equal
                if other.pos[1] > self.pos[1]:
                    return False
                if other.pos[1] < self.pos[1]:
                    return True
                #Y's equal
                if other.pos[0] > self.pos[0]:
                    return False
                if other.pos[0] < self.pos[0]:
                    return True
                
                #Equal
                return False
            
            def __gt__(self, other):
                self.otherCellPosX = int(other.pos[0])/Na
                self.otherCellPosY = int(other.pos[1])/Nb
                self.otherCellPosZ = int(other.pos[2])/Nc
                
                if otherCellPosX > self.cellPosX:
                    return True
                if otherCellPosX < self.cellPosX:
                    return False
                #X's equal
                if otherCellPosY > self.cellPosY:
                    return True
                if otherCellPosY < self.cellPosY:
                    return False
                #Y's equal
                if otherCellPosZ > self.cellPosZ:
                    return True
                if otherCellPosZ < self.cellPosZ:
                    return False
                
                #All equal
                if other.pos[2] > self.pos[2]:
                    return True
                if other.pos[2] < self.pos[2]:
                    return False
                #Z's equal
                if other.pos[1] > self.pos[1]:
                    return True
                if other.pos[1] < self.pos[1]:
                    return False
                #Y's equal
                if other.pos[0] > self.pos[0]:
                    return True
                if other.pos[0] < self.pos[0]:
                    return False
                
                
                return False
            
            def __eq__(self, other):
                if self.pos[0] == other.pos[0]:
                    if self.pos[1] == other.pos[1]:
                        if self.pos[2] == other.pos[2]:
                            return True
                return False
            
            def  __le__(self, other):
                return (self.__eq__(other) or self.__lt__(other))
            
            def __ne__(self, other):
                return not self.__eq__(other)
            
            def __ge__(self, other):
                return (self.__eq__(other) or self.__gt__(other))
                
           
                
        class SimpleAtomList():
            def __init__(self):
                self.atoms = []
                
            def addBond(self, pos1, pos2, jMatrix):
                #Make sure each position is not yet represented the easy but inefficient way
                pos1Index = -1
                pos2Index = -1
                for i in range(len(self.atoms)):
                    currentPos = self.atoms[i].pos
                    if currentPos == pos1:
                        pos1Index = i
                    if currentPos == pos2:
                        pos2Index = i
                
                if pos1Index < 0:
                    atom1 = SimpleAtom(pos1, len(self.atoms))
                    self.atoms.append(atom1)
#                    print self.atoms[atom1.index] == atom1
                else:
                    atom1 = self.atoms[pos1Index]
                    
                if pos2Index < 0:
                    atom2 = SimpleAtom(pos2, len(self.atoms))
                    self.atoms.append(atom2)
#                    print self.atoms[atom2.index] == atom2
                else:
                    atom2 = self.atoms[pos2Index]
                    
                atom1.addInteraction(atom2, jMatrix)
                atom2.addInteraction(atom1, jMatrix)
                
        
        
        
        class sortedAtomList():
            def __init__(self):
                self.atoms = []
                
            def append(self, item):
                self.atoms.append(item)
            
            def extend(self, otherList):
                self.atoms.extend(otherList)
                
            def translate(self, x, y, z):
                newList = sortedAtomList()
                for atom in self.atoms:
                    pos = atom.pos
                    newAtom = simpleAtom((pos[0] + x, pos[1] + y, pos[2] + z))
                    newInteractions = {}
                    for position in atom.interactions:
                        newInteractions[(position[0] + x, position[1] + y, position[2] + z)] = atom.interactions[position]
                    newAtom.interactions = newInteractions
                    newList.append(newAtom)
                
            def add(self, item):
                index = addAux(0, len(self.atoms), item)
                
#            def addAux(self, start, end):
                
            
            def find(self, item):
                index = findAux(0, len(self.atoms), item)
                if index >= 0:
                    return self.atoms[index]
                return None
                
            def findAux(self, start, end):
                if end < start:
                    return -1
                mid = (end + start)/2
                if item < self.atoms[mid]:
                    return findAux(start, mid - 1, item)
                elif item > self.atoms[mid]:
                    return findAux(mid + 1, end, item)
                
                return mid
                
        
        
        class SimpleBondList():
            def __init__(self):
                self.list = []
                
            def sort(self):
                self.list.sort()
                
            def addBond(self, bond):
                if not self.containsBond(bond):
                    self.list.append(bond)
                    
            def containsBond(self, bond):
                for eachBond in self.list:
                    if eachBond.sameBond(bond):
                        return True
                return False
        
        
        def contains(list, element):
            for item in list:
                if (item == element).all():
                    return True
            return False
        
        def indexOf(list, item):
            for i in range(len(list)):
                if (item == list[i]).all():
                    return i
            return -1
        
        
        #Create list of matrices
        matrices = []
        for bond in self.getCutoffCell().getBonds():
#           pos1 = bond.getAtom1().getPosition()
#           pos2 = bond.getAtom2().getPosition()
            jMat = bond.getJMatrix()
#            count = matrices.count(jMat)
            if not contains(matrices, jMat):
                matrices.append(jMat)
        
        
        #create simple bonds within cutoff cell
        simpleCellBonds = []
        for bond in self.getCutoffCell().getBonds():
            pos1 = bond.getAtom1().getPosition()
            pos2 = bond.getAtom2().getPosition()
            jMat = bond.getJMatrix()
            newBond = SimpleBond(pos1, pos2, indexOf(matrices,jMat))
            simpleCellBonds.append(newBond)
        
        
        simpleBonds = SimpleBondList()
        for bond in simpleCellBonds:
            pos1 = bond.pos1
            pos2 = bond.pos2
            jMatInt = bond.jMatrix
            for i in range(2):
                for j in range(2):
                    for k in range(2):
                        for a in range(Na):
                            for b in range(Nb):
                                for c in range(Nc):
                                    x1 = pos1[0] + a + (Na * i)
                                    y1 = pos1[1] + b + (Nb * j)
                                    z1 = pos1[2] + c + (Nc * k)
                                    
                                    x2 = pos2[0] + a + (Na * i)
                                    y2 = pos2[1] + b + (Nb * j)
                                    z2 = pos2[2] + c + (Nc * k)    
                                    bond = SimpleBond( (x1,y1,z1), (x2,y2,z2), jMatInt )
                                    simpleBonds.addBond(bond)                              

        def PosInFirstCutoff(pos):
            return (pos[0] < Na and pos[1] < Nb and pos[2] < Nc)
        
        
        #Pick out bonds that link first cutoff cell and another
        interCutoffBonds = []
        for eachBond in simpleBonds.list:
#            print "atom1: ", eachBond.pos1, "atom2: ", eachBond.pos2
            #Check if bond1 is in the first cutoff cell, but not bond 2
            if PosInFirstCutoff(eachBond.pos1) and not PosInFirstCutoff(eachBond.pos2):
                interCutoffBonds.append(eachBond)
            #Check if bond2 is in the first cutoff cell, but not bond 1
            if PosInFirstCutoff(eachBond.pos2) and not PosInFirstCutoff(eachBond.pos1):
                interCutoffBonds.append(eachBond)
                                    
     
#        print "intercutoffBonds: " , len(interCutoffBonds)
#        for nBond in interCutoffBonds:
#            print bond.pos1, bond.pos2
#        return
     
     
     #The easy but slow way
         
        file.write("#J Matrices\n#Number J11 J12 J13 J21 J22 J23 J31 J32 J33\n")
        for i in range(len(matrices)):
            jMat = matrices[i]
            jStr = str(i) + " " + str(jMat[0][0]) + " " + str(jMat[0][1]) + " " + str(jMat[0][2]) + " " + str(jMat[1][0]) + " " + str(jMat[1][1]) + " " + str(jMat[1][2]) + " " + str(jMat[2][0]) + " " + str(jMat[2][1]) + " " + str(jMat[2][2])
            file.write(jStr + "\n")
        
        
        finalBondList = []
        for bond in simpleCellBonds:   
            pos1 = bond.pos1
            pos2 = bond.pos2
#            jStr += " " 
#            jStr += str(jMat[0][1]) 
#            jStr += " " + str(jMat[0][2]) 
#            jStr += " " + str(jMat[1][0]) + " " + str(jMat[1][1]) + " " + str(jMat[1][2]) + " " + str(jMat[2][0]) + " " + str(jMat[2][1]) + " " + str(jMat[2][2])
            for i in range(size):
                for j in range(size):
                    for k in range(size):
                        x1 = pos1[0] + (Na * i)
                        y1 = pos1[1] + (Nb * j)
                        z1 = pos1[2] + (Nc * k)
                        
                        x2 = pos2[0] + (Na * i)
                        y2 = pos2[1] + (Nb * j)
                        z2 = pos2[2] + (Nc * k)    
#                        smplBond = SimpleBond( (x1,y1,z1), (x2,y2,z2), jMat )
 #                       finalBondList.append(smplBond)
                        pos1Str = str(x1) + " " + str(y1) + " " + str(z1)
                        pos2Str = str(x2) + " " + str(y2) + " " + str(z2)
                        finalBondList.append(SimpleBond((x1,y1,z1),(x2,y2,z2),bond.jMatrix))

        print "Added all bond within cutoff Cells"
        
        for smplBond in interCutoffBonds:
            pos1 = smplBond.pos1
            pos2 = smplBond.pos2
            jMat = smplBond.jMatrix
            jStr = str(jMat)
#            jStr = str(jMat[0][0]) + " " + str(jMat[0][1]) + " " + str(jMat[0][2]) + " " + str(jMat[1][0]) + " " + str(jMat[1][1]) + " " + str(jMat[1][2]) + " " + str(jMat[2][0]) + " " + str(jMat[2][1]) + " " + str(jMat[2][2])
            aDisp = abs(pos1[0] - pos2[0])
            bDisp = abs(pos1[1] - pos2[1])
            cDisp = abs(pos1[2] - pos2[2])
            for i in range(size - aDisp):
                for j in range(size - bDisp):
                    for k in range(size - cDisp):
                        x1 = pos1[0] + (Na * i)
                        y1 = pos1[1] + (Nb * j)
                        z1 = pos1[2] + (Nc * k)
                        
                        x2 = pos2[0] + (Na * i)
                        y2 = pos2[1] + (Nb * j)
                        z2 = pos2[2] + (Nc * k)    
#                        smplBond = SimpleBond( (x1,y1,z1), (x2,y2,z2), jMat )
 #                       finalBondList.append(smplBond)
                        pos1Str = str(x1) + " " + str(y1) + " " + str(z1)
                        pos2Str = str(x2) + " " + str(y2) + " " + str(z2)
                        finalBondList.append(SimpleBond((x1,y1,z1),(x2,y2,z2),smplBond.jMatrix))
        
        print "added all bonds, checking for reapeats"
        timer.printTime()
        
        #Check for reapeats in finalBond List just for testing
        def isRepeat(finalBondList):
            for i in range(0, len(finalBondList)):
                for j in range(i + 1, len(finalBondList)):
                    if finalBondList[i].sameBond(finalBondList[j]):
                        return True
            return False
        
        if isRepeat(finalBondList):
            print "There is a repeat!"
        else:
            print "NO repeats!"
            
        timer.printTime()
        
        atomList = SimpleAtomList()
        for bond in finalBondList:
            atomList.addBond(bond.pos1, bond.pos2, bond.jMatrix)
            
        print "bonds added to associative atom list, checking if it is balanced"
        timer.printTime()
        
        #Check the simple atom list
        def atomBalanced(atom):
            for otherAtom in atomList.atoms:
                if atomInteractsWithAtom(atom, otherAtom):
                    if atomInteractsWithAtom(otherAtom, atom):
                        return True
                    else:
                        return False
            return False
                
                
        def atomInteractsWithAtom(atom, otherAtom):
            for interaction in otherAtom.interactions:
                if atom == interaction[0]:
                    return True
            return False
        
        
        for atom in atomList.atoms:
            if not atomBalanced(atom):
                print "Not Balanced!!!"
                break
        else:
            print "Balanced!"
            
        
        print "number of atoms: ", len(atomList.atoms)
            
        
        #print out the simple atom list
        file.write("#AtomNumber AtomPosition(X Y Z) OtherIndex OtherPos Jmatrix OtherIndex OtherPos Jmatrix...\n")
        for atom in atomList.atoms:
            atomStr = str(atom.index) + " " + str(atom.pos[0]) + " " + str(atom.pos[1]) + " " + str(atom.pos[2])
            for interactionKey in atom.interactions:
                otherAtom = interactionKey[0]
                jMat = interactionKey[1]
                atomStr += " " + str(otherAtom.index) + " " + str(otherAtom.pos[0]) + " " + str(otherAtom.pos[1]) + " " + str(otherAtom.pos[2])
                atomStr += " " + str(jMat)
            file.write(atomStr + "\n")
        
        file.close()
            
        timer.printTime()
            
        
        #More efficient method with a sorted list and eventually in layers
#        atoms = sortedAtomList()
#        for bond in simpleCellBonds:
#            pos1 = bond.pos1
#            pos2 = bond.pos2
#            jMat = bond.jMatrix
#            atom1 = atoms.find(SimpleAtom(pos1))
#            atom2 = atoms.find(SimpleAtom(pos2))
#            if not atom1:
#                atom1 = SimpleAtom(pos1)
#                atoms.append(atom1)
#            if not atom2:
#                atom2 = SimpleAtom(pos2)
#                atoms.append(atom2)
#            atom1.add(atom2, jMat)
#            atom2.add(atom1, jMat)
 #       atoms.sort()
 #       
 #       
 #       #This can all be done in layers of the Z plane later to save memory
 #       #rather than dealing with the entire cube at once
 #       
 #       #translate the cutoff cell and the bonds contained within it
 #       allAtoms = sortedAtomList()
 #       for i in range(size):
 #           for j in range(size):
 #               for k in range(size):
 #                   atoms2 = atoms.translate(i, j, k)
 #                   allAtoms.extend(atoms2.atoms)
 #                   
 #       
 # 
 #       for smplBond in interCutoffBonds:
 #           pos1 = smplBond.pos1
 #           pos2 = smplBond.pos2
 #           jMat = smplBond.jMatrix
##            jStr = str(jMat[0][0]) + " " + str(jMat[0][1]) + " " + str(jMat[0][2]) + " " + str(jMat[1][0]) + " " + str(jMat[1][1]) + " " + str(jMat[1][2]) + " " + str(jMat[2][0]) + " " + str(jMat[2][1]) + " " + str(jMat[2][2])
 #           aDisp = abs(pos1[0] - pos2[0])
 #           bDisp = abs(pos1[1] - pos2[1])
 #           cDisp = abs(pos1[2] - pos2[2])
 #           for i in range(size - aDisp):
 #               for j in range(size - bDisp):
 #                   for k in range(size - cDisp):
                        #The one that was translated from the original cell can be simply appended
  
            
            #check if it is sorted






    def exportForMonteCarlo(self, filename):
        size = 2
        
        timer = Timer()
        
        file = open(filename, 'w')
   
        class SimpleBond():
            def __init__(self, pos1, pos2, jMatrix):
                self.pos1 = pos1
                self.pos2 = pos2
                self.jMatrix = jMatrix
                
            def sameBond(self, bond2):
                if self.pos1 == bond2.pos1 or self.pos1 == bond2.pos2:
                    if self.pos2 == bond2.pos2 or self.pos2 == bond2.pos1:
                        return True
                return False
        
        Na = self.getCutoffCell().getNa()
        Nb = self.getCutoffCell().getNb()
        Nc = self.getCutoffCell().getNc()
        
        class SimpleAtom():
            def __init__(self, pos):
                self.pos = pos
#                self.cellPos = []
#                self.cellPosX = int(pos[0])/Na
#                self.cellPosY = int(pos[1])/Nb
#                self.cellPosZ = int(pos[2])/Nc
                self.interactions = []
                #self.interactions[position of other atom] = j number
                
            #might want to change this to position later when all atoms wont be created in same list
            def addInteraction(self, atom2, jMat):
#                self.interactions[atom2] = jMat
                self.interactions.append([atom2, jMat])
            
            #comparisons based on positions
            def __lt__(self, other):
                self.otherCellPosX = int(other.pos[0])/Na
                self.otherCellPosY = int(other.pos[1])/Nb
                self.otherCellPosZ = int(other.pos[2])/Nc
                
                if otherCellPosX > self.cellPosX:
                    return False
                if otherCellPosX < self.cellPosX:
                    return True
                #X's equal
                if otherCellPosY > self.cellPosY:
                    return False
                if otherCellPosY < self.cellPosY:
                    return True
                #Y's equal
                if otherCellPosZ > self.cellPosZ:
                    return False
                if otherCellPosZ < self.cellPosZ:
                    return True
                
                #They are in the same cell
                if other.pos[2] > self.pos[2]:
                    return False
                if other.pos[2] < self.pos[2]:
                    return True
                #Z's equal
                if other.pos[1] > self.pos[1]:
                    return False
                if other.pos[1] < self.pos[1]:
                    return True
                #Y's equal
                if other.pos[0] > self.pos[0]:
                    return False
                if other.pos[0] < self.pos[0]:
                    return True
                
                #Equal
                return False
            
            def __gt__(self, other):
                self.otherCellPosX = int(other.pos[0])/Na
                self.otherCellPosY = int(other.pos[1])/Nb
                self.otherCellPosZ = int(other.pos[2])/Nc
                
                if otherCellPosX > self.cellPosX:
                    return True
                if otherCellPosX < self.cellPosX:
                    return False
                #X's equal
                if otherCellPosY > self.cellPosY:
                    return True
                if otherCellPosY < self.cellPosY:
                    return False
                #Y's equal
                if otherCellPosZ > self.cellPosZ:
                    return True
                if otherCellPosZ < self.cellPosZ:
                    return False
                
                #All equal
                if other.pos[2] > self.pos[2]:
                    return True
                if other.pos[2] < self.pos[2]:
                    return False
                #Z's equal
                if other.pos[1] > self.pos[1]:
                    return True
                if other.pos[1] < self.pos[1]:
                    return False
                #Y's equal
                if other.pos[0] > self.pos[0]:
                    return True
                if other.pos[0] < self.pos[0]:
                    return False
                
                
                return False
            
            def __eq__(self, other):
                if self.pos[0] == other.pos[0]:
                    if self.pos[1] == other.pos[1]:
                        if self.pos[2] == other.pos[2]:
                            return True
                return False
            
            def  __le__(self, other):
                return (self.__eq__(other) or self.__lt__(other))
            
            def __ne__(self, other):
                return not self.__eq__(other)
            
            def __ge__(self, other):
                return (self.__eq__(other) or self.__gt__(other))
                
           
                
        class SimpleAtomList():
            def __init__(self):
                self.atoms = []
                
            def addBond(self, pos1, pos2, jMatrix):
                #Make sure each position is not yet represented the easy but inefficient way
                pos1Index = -1
                pos2Index = -1
                for i in range(len(self.atoms)):
                    currentPos = self.atoms[i].pos
                    if currentPos == pos1:
                        pos1Index = i
                    if currentPos == pos2:
                        pos2Index = i
                
                if pos1Index < 0:
                    atom1 = SimpleAtom(pos1)
                    self.atoms.append(atom1)
#                    print self.atoms[atom1.index] == atom1
#                else:
#                    atom1 = self.atoms[pos1Index]
                    
                if pos2Index < 0:
                    atom2 = SimpleAtom(pos2)
                    self.atoms.append(atom2)
#                    print self.atoms[atom2.index] == atom2
#                else:
#                    atom2 = self.atoms[pos2Index]
                    
#Unnecesary the way it is currently beig used
#                atom1.addInteraction(atom2, jMatrix)
#                atom2.addInteraction(atom1, jMatrix)
                  
        
        
        class sortedAtomList():
            def __init__(self):
                self.atoms = []
                
            def append(self, item):
                self.atoms.append(item)
            
            def extend(self, otherList):
                self.atoms.extend(otherList)
                
            def translate(self, x, y, z):
                newList = sortedAtomList()
                for atom in self.atoms:
                    pos = atom.pos
                    newAtom = simpleAtom((pos[0] + x, pos[1] + y, pos[2] + z))
                    newInteractions = {}
                    for position in atom.interactions:
                        newInteractions[(position[0] + x, position[1] + y, position[2] + z)] = atom.interactions[position]
                    newAtom.interactions = newInteractions
                    newList.append(newAtom)
                
            def add(self, item):
                index = addAux(0, len(self.atoms), item)
                
#            def addAux(self, start, end):
                
            
            def find(self, item):
                index = findAux(0, len(self.atoms), item)
                if index >= 0:
                    return self.atoms[index]
                return None
                
            def findAux(self, start, end):
                if end < start:
                    return -1
                mid = (end + start)/2
                if item < self.atoms[mid]:
                    return findAux(start, mid - 1, item)
                elif item > self.atoms[mid]:
                    return findAux(mid + 1, end, item)
                
                return mid
                
        
        class SimpleBondList():
            def __init__(self):
                self.list = []
                
            def sort(self):
                self.list.sort()
                
            def addBond(self, bond):
                if not self.containsBond(bond):
                    self.list.append(bond)
#                else:
#                    print "Duplicate Bonds!" #should not get here
                    
            def containsBond(self, bond):
                for eachBond in self.list:
                    if eachBond.sameBond(bond):
                        return True
                return False
        
        
        def contains(list, element):
            for item in list:
                if (item == element).all():
                    return True
            return False
        
        def atomListContains(list, element):
            for item in list:
                if item == element:
                    return True
            return False
        
        def indexOf(list, item):
            for i in range(len(list)):
                if (item == list[i]).all():
                    return i
            return -1
        
        
        #Create list of matrices
        matrices = []
        for bond in self.getCutoffCell().getBonds():
#           pos1 = bond.getAtom1().getPosition()
#           pos2 = bond.getAtom2().getPosition()
            jMat = bond.getJMatrix()
#            count = matrices.count(jMat)
            if not contains(matrices, jMat):
                matrices.append(jMat)
        
        
        #create simple bonds within cutoff cell
        simpleCellBonds = []
        for bond in self.getCutoffCell().getBonds():
            pos1 = bond.getAtom1().getPosition()
            pos2 = bond.getAtom2().getPosition()
            jMat = bond.getJMatrix()
            newBond = SimpleBond(pos1, pos2, indexOf(matrices,jMat))
            simpleCellBonds.append(newBond)
        
        
        def PosInFirstCutoff(pos):
            return (pos[0] < Na and pos[1] < Nb and pos[2] < Nc)
        
        
        cellAtoms = []
        cellBonds = SimpleBondList()
        for bond in simpleCellBonds:
            pos1 = bond.pos1
            pos2 = bond.pos2
            jMatInt = bond.jMatrix
            for i in range(2):
                for j in range(2):
                    for k in range(2):
                        for a in range(Na):
                            for b in range(Nb):
                                for c in range(Nc):
                                    x1 = pos1[0] + a + (Na * i)
                                    y1 = pos1[1] + b + (Nb * j)
                                    z1 = pos1[2] + c + (Nc * k)
                                    
                                    x2 = pos2[0] + a + (Na * i)
                                    y2 = pos2[1] + b + (Nb * j)
                                    z2 = pos2[2] + c + (Nc * k)  
                                    newPos1 = (x1,y1,z1)
                                    newPos2 = (x2,y2,z2)
                                    if PosInFirstCutoff(newPos1):
                                        newAtom = SimpleAtom(newPos1)
                                        bond = SimpleBond( (x1,y1,z1), (x2,y2,z2), jMatInt )
                                        if not atomListContains(cellAtoms, newAtom):
                                            cellAtoms.append(newAtom)
                                        cellBonds.addBond(bond)
                                    if PosInFirstCutoff(newPos2):
                                        newAtom = SimpleAtom(newPos2)
                                        bond = SimpleBond( (x1,y1,z1), (x2,y2,z2), jMatInt )
                                        if not atomListContains(cellAtoms, newAtom):
                                            cellAtoms.append(newAtom)
                                        cellBonds.addBond(bond)
                
        

         
        file.write("#J Matrices\n#Number J11 J12 J13 J21 J22 J23 J31 J32 J33\n")
        for i in range(len(matrices)):
            jMat = matrices[i]
            jStr = str(i) + " " + str(jMat[0][0]) + " " + str(jMat[0][1]) + " " + str(jMat[0][2]) + " " + str(jMat[1][0]) + " " + str(jMat[1][1]) + " " + str(jMat[1][2]) + " " + str(jMat[2][0]) + " " + str(jMat[2][1]) + " " + str(jMat[2][2])
            file.write(jStr + "\n")
        
        
        
        allAtoms = []
        numAtomsPerCell = len(cellAtoms)
        
        print "atoms in cell: ", numAtomsPerCell
        
        for i in range(size):
            for j in range(size):
                for k in range(size):
                    for index in range(len(cellAtoms)):
                        pos = cellAtoms[index].pos
                        x = pos[0] + (Na * i)
                        y = pos[1] + (Nb * j)
                        z = pos[2] + (Nc * k)
                        newAtom = SimpleAtom((x,y,z))
#                        print (len(allAtoms)%numAtomsPerCell == index)#just a check, should always be true
                        allAtoms.append(newAtom)

 
        #Add bonds cellBonds to allAtoms (currently not most efficient way, but its a short list
        for bond in cellBonds.list:
            #Make sure each position is not yet represented the easy but inefficient way
            pos1 = bond.pos1
            pos2 = bond.pos2
            pos1Index = -1
            pos2Index = -1
            for i in range(len(allAtoms)):
                currentPos = allAtoms[i].pos
                if currentPos == pos1:
                    pos1Index = i
                    break
                
            for i in range(len(allAtoms)):
                currentPos = allAtoms[i].pos  
                if currentPos == pos2:
                    pos2Index = i
                    break
            
            if pos1Index < 0 or pos2Index < 0:
                print "Atom list does not contain all atoms!" 
            else:
                allAtoms[pos1Index].addInteraction(pos2Index, bond.jMatrix)
                allAtoms[pos2Index].addInteraction(pos1Index, bond.jMatrix)


        timer.printTime()
        print"translating..."
        
        #translate bonds
        for i in range(len(allAtoms)- numAtomsPerCell):
            for interaction in allAtoms[i].interactions:
                newInteraction = interaction[0] + numAtomsPerCell
                if newInteraction < len(allAtoms):
                    allAtoms[i+numAtomsPerCell].addInteraction(newInteraction, interaction[1])
        
        
        
#        print "done translating, checking list"
#        timer.printTime()
        
        #Check for reapeats in finalBond List just for testing
#        def isRepeat(finalBondList):
#            for i in range(0, len(finalBondList)):
#                for j in range(i + 1, len(finalBondList)):
#                    if finalBondList[i].sameBond(finalBondList[j]):
#                        return True
#            return False
#        
#        if isRepeat(finalBondList):
#            print "There is a repeat!"
#        else:
#            print "NO repeats!"
#            
#        timer.printTime()
        
        
        #Check the simple atom list
        def atomBalanced(atomIndex):
            atom = allAtoms[atomIndex]
            for otherAtomIndex in range(len(allAtoms)):
                otherAtom = allAtoms[otherAtomIndex]
                if atomInteractsWithAtom(atomIndex, otherAtom):
                    if atomInteractsWithAtom(otherAtomIndex, atom):
                        return True
                    else:
                        return False
            return False
        
        def atomInteractsWithAtom(atomIndex, otherAtom):
            for interaction in otherAtom.interactions:
                if atomIndex == interaction[0]:
                    return True
            return False
        
        
#        for atomIndex in range(len(allAtoms)):
#            if not atomBalanced(atomIndex):
#                print "Not Balanced!!!"
#                break
#        else:
#            print "Balanced!"
            
        
        timer.printTime()
        print "number of atoms: ", len(allAtoms), "\n writing to disk..."
            
        
        #print out the simple atom list
        file.write("#AtomNumber AtomPosition(X Y Z) OtherIndex Jmatrix OtherIndex Jmatrix...\n")
        for atomIndex in range(len(allAtoms)):
            atom = allAtoms[atomIndex]
            atomStr = str(atomIndex) + " " + str(atom.pos[0]) + " " + str(atom.pos[1]) + " " + str(atom.pos[2])
            for interaction in atom.interactions:
                otherAtom = interaction[0]
                jMat = interaction[1]
                atomStr += " " + str(otherAtom)
                atomStr += " " + str(jMat)
            file.write(atomStr + "\n")
        
        file.close()
        print "done"
        timer.printTime()

