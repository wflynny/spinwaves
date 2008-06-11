import copy

class MagneticCell():
    
    def __init__(self, Unit_Cell, Na, Nb, Nc):
        self.unit_cell = Unit_Cell
        
        self.AllUnitCells = [Unit_Cell]
        
        self.IntercellularBonds = []
        
        #generate the magnetic Cell by translating the unit cell
        for i in range(0, Na):
            for j in range(0, Nb):
                for k in range(0, Nc):
                    self.AllUnitCells.append(self.unit_cell.translateCell(i,j,k))
    
    
    
    
    def addInterCellularBond(self, Atom1, Atom2):
        self.IntercellularBonds.append(Bond())
    
    def drawCell(self, Renderer):
        #draw all Cells
        for cell in self.AllUnitCells:
            cell.drawCell(Renderer)
    
    def getAllUnitCells(self):
        return self.AllUnitCells
    
