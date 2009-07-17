import sympy as sp
import sympy.matrices as spm
import numpy as N
from spinwaves.spinwavecalc.spinwave_calc_file import calculate_dispersion
import spinwaves.spinwavecalc.readfiles as rf
import unittest

#------------------------------------------------------------------------------------------------------
#-----------------T E S T S----------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------

class TestFM(unittest.TestCase):   
    def setUp(self):
        self.S = sp.Symbol('S', real = True)
        self.J = sp.Symbol('J', real = True)
        self.kx = sp.Symbol('kx', real = True)
        self.ky = sp.Symbol('ky', real = True)
        self.kz = sp.Symbol('kz', real = True)

    def testChain(self, num = 3, num_uc = 1):
        ats = [rf.atom(pos = [i,0,0]) for i in range(num)]
        spin = rf.findmat(N.array([0,0,1]))
        
        ints = [[0],[0,1],[1]]
        nbrs = [[1],[0,2],[1]]
        for i in range(num): 
            ats[i].interactions = ints[i]
            ats[i].neighbors = nbrs[i]
            ats[i].spinRmatrix = spin

        inters = [item for atom in ats for item in atom.interactions]
        max_inter = max(inters)
        
        J = sp.Symbol('J', real = True)
        Jij = [N.matrix([[J,0,0],[0,J,0],[0,0,J]]) for i in range(max_inter+1)]
        (Hsave, charpoly, eigs) = calculate_dispersion(ats, num_uc, num, Jij, showEigs = True)

        eigs = eigs.keys()
        for item in eigs:
            item = item.expand(mul = True, multinomial=True)
        
        testeig = (-8.0*self.J**2*self.S**2*sp.cos(self.kx) + 4.0*self.J**2*self.S**2 + 4.0*self.J**2*self.S**2*sp.cos(self.kx)**2)**(0.5)
        
        self.assertAlmostEqual(Hsave.shape[0],2,4, 'H wrong shape')
        self.assertAlmostEqual(Hsave.shape[1],2,4, 'H wrong shape')
        self.assert_(len(eigs[0].args)==len(testeig.args), 'Eig.arg lengths are different sizes')
        f = sp.simplify(eigs[0].args[0].expand())
        f = f.subs([(self.J,1),(self.S,2),(self.kx,3),(self.ky,4),(self.kz,5)])
        g = sp.simplify(testeig.args[0].expand())
        g = g.subs([(self.J,1),(self.S,2),(self.kx,3),(self.ky,4),(self.kz,5)])
        myFlag = f-g
        self.assertAlmostEqual(myFlag,0,7, 'Eigs are wrong')

    def testSquare(self, num = 5, num_uc = 1):
        ats = [rf.atom() for i in range(num)]
        spin = rf.findmat(N.array([0,0,1]))
        # Positions
        #        0        1        2        3        4
        pstn = [[0,0,0],[0,1,0],[1,0,0],[0,-1,0],[-1,0,0]]
        # Neighbors
        #        0        1    2   3   4
        nbrs = [[1,2,3,4],[0],[0],[0],[0]]
        # Interactions
        #        0        1    2   3   4
        ints = [[0,1,2,3],[0],[1],[2],[3]]
    
        for i in range(len(ats)):
            ats[i].pos = N.array(pstn[i])
            ats[i].neighbors = nbrs[i]
            ats[i].interactions = ints[i]
            ats[i].spinRmatrix = spin
            ats[i].spinMagnitude = 1

        inters = [item for atom in ats for item in atom.interactions]
        max_inter = max(inters)
        
        J = sp.Symbol('J', real = True)
        Jij = [N.matrix([[J,0,0],[0,J,0],[0,0,J]]) for i in range(max_inter+1)]
        (Hsave, charpoly, eigs) = calculate_dispersion(ats, num_uc, num, Jij, showEigs = True)

        eigs = eigs.keys()
        for item in eigs:
            item = item.expand(mul = True, multinomial=True)

        testeig = (-64.0*self.J**2*self.S**2*sp.cos(self.kx) - 64.0*self.J**2*self.S**2*sp.cos(self.ky) + 32.0*self.J**2*self.S**2*sp.cos(self.kx)*sp.cos(self.ky) + 64.0*self.J**2*self.S**2 + 16.0*self.J**2*self.S**2*sp.cos(self.kx)**2 + 16.0*self.J**2*self.S**2*sp.cos(self.ky)**2)**(0.5)

        self.assertAlmostEqual(Hsave.shape[0],2,4, 'H wrong shape')
        self.assertAlmostEqual(Hsave.shape[1],2,4, 'H wrong shape')
        self.assert_(len(eigs[1].args)==len(testeig.args), 'Eig.arg lengths are different sizes')
        f = sp.simplify(eigs[1].args[0].expand())
        f = f.subs([(self.J,1),(self.S,2),(self.kx,3),(self.ky,4),(self.kz,5)])
        g = sp.simplify(testeig.args[0].expand())
        g = g.subs([(self.J,1),(self.S,2),(self.kx,3),(self.ky,4),(self.kz,5)])
        myFlag = f-g
        self.assertAlmostEqual(myFlag,0,7, 'Eigs are wrong')

    def testCube(self, num = 7, num_uc = 1):
        ats = [rf.atom() for i in range(num)]
        spin = rf.findmat(N.array([0,0,1]))
        # Positions
        #          0        1         2        3        4         5        6
        pstn = [[0,0,0], [1,0,0], [-1,0,0], [0,1,0], [0,-1,0], [0,0,1], [0,0,-1]]
        # Neighbors
        #            0          1    2    3    4    5    6
        nbrs = [[1,2,3,4,5,6], [0], [0], [0], [0], [0], [0]]
        # Interactions
        #            0          1    2    3    4    5    6
        ints = [[0,0,0,0,0,0], [1], [2], [3], [4], [5], [6]]
    
        for i in range(len(ats)):
            ats[i].pos = N.array(pstn[i])
            ats[i].neighbors = nbrs[i]
            ats[i].interactions = ints[i]
            ats[i].spinRmatrix = spin
            ats[i].spinMagnitude = 1

        inters = [item for atom in ats for item in atom.interactions]
        max_inter = max(inters)
        
        J = sp.Symbol('J', real = True)
        Jij = [N.matrix([[J,0,0],[0,J,0],[0,0,J]]) for i in range(max_inter+1)]
        (Hsave, charpoly, eigs) = calculate_dispersion(ats, num_uc, num, Jij, showEigs = True)

        eigs = eigs.keys()
        for item in eigs:
            item = item.expand(mul = True, multinomial=True)

        testeig = (-96.0*self.J**2*self.S**2*sp.cos(self.kx) - 96.0*self.J**2*self.S**2*sp.cos(self.ky) - 96.0*self.J**2*self.S**2*sp.cos(self.kz) + 32.0*self.J**2*self.S**2*sp.cos(self.kx)*sp.cos(self.ky) + 32.0*self.J**2*self.S**2*sp.cos(self.kx)*sp.cos(self.kz) + 32.0*self.J**2*self.S**2*sp.cos(self.ky)*sp.cos(self.kz) + 144.0*self.J**2*self.S**2 + 16.0*self.J**2*self.S**2*sp.cos(self.kx)**2 + 16.0*self.J**2*self.S**2*sp.cos(self.ky)**2 + 16.0*self.J**2*self.S**2*sp.cos(self.kz)**2)**(0.5)

        self.assertAlmostEqual(Hsave.shape[0],2,4, 'H wrong shape')
        self.assertAlmostEqual(Hsave.shape[1],2,4, 'H wrong shape')
        self.assert_(len(eigs[1].args)==len(testeig.args), 'Eig.arg lengths are different sizes')
        f = sp.simplify(eigs[1].args[0].expand())
        f = f.subs([(self.J,1),(self.S,2),(self.kx,3),(self.ky,4),(self.kz,5)])
        g = sp.simplify(testeig.args[0].expand())
        g = g.subs([(self.J,1),(self.S,2),(self.kx,3),(self.ky,4),(self.kz,5)])
        myFlag = f-g
        self.assertAlmostEqual(myFlag,0,7, 'Eigs are wrong')
#
    def testBCC(self, num = 9, num_uc = 1):
        ats = [rf.atom() for i in range(num)]
        spin = rf.findmat(N.array([0,0,1]))
        # Positions
        #          0            1               2               3               4              
        pstn = [[0,0,0], [-0.5,0.5,-0.5], [0.5,0.5,-0.5], [0.5,0.5,0.5], [-0.5,0.5,0.5],
        #               5                6                7               8               
                [-0.5,-0.5,-0.5], [0.5,-0.5,-0.5], [0.5,-0.5,0.5], [-0.5,-0.5,0.5]]
        # Neighbors
        #          0               1   2   3   4              
        nbrs = [[1,2,3,4,5,6,7,8], [], [], [], [],
        #        5   6   7   8               
                [], [], [], []]
        # Interactions
        #          0               1   2   3   4              
        ints = [[0,1,2,3,4,5,6,7], [], [], [], [],
        #        5   6   7   8               
                [], [], [], []]

        for i in range(len(ats)):
            ats[i].pos = N.array(pstn[i])
            ats[i].neighbors = nbrs[i]
            ats[i].interactions = ints[i]
            ats[i].spinRmatrix = spin
            ats[i].spinMagnitude = 1

        inters = [item for atom in ats for item in atom.interactions]
        max_inter = max(inters)
        
        J = sp.Symbol('J', real = True)
        Jij = [N.matrix([[J,0,0],[0,J,0],[0,0,J]]) for i in range(max_inter+1)]
        (Hsave, charpoly, eigs) = calculate_dispersion(ats, num_uc, num, Jij, showEigs = True)

        eigs = eigs.keys()
        for item in eigs:
            item = item.expand(mul = True, multinomial=True)

        testeig = (-2048.0*self.J**2*self.S**2*sp.cos(0.5*self.kx)*sp.cos(0.5*self.ky)*sp.cos(0.5*self.kz) + 1024.0*self.J**2*self.S**2 + 1024.0*self.J**2*self.S**2*sp.cos(0.5*self.kx)**2*sp.cos(0.5*self.ky)**2*sp.cos(0.5*self.kz)**2)**(0.5)*0.5

        self.assertAlmostEqual(Hsave.shape[0],2,4, 'H wrong shape')
        self.assertAlmostEqual(Hsave.shape[1],2,4, 'H wrong shape')
        f = eigs[0].args[2].args[1].args[0].expand()
        f = f.subs([(self.J,1),(self.S,2),(self.kx,3),(self.ky,4),(self.kz,5),(sp.I,0)])
        g = sp.simplify(testeig.args[1].args[0].expand())
        g = g.subs([(self.J,1),(self.S,2),(self.kx,3),(self.ky,4),(self.kz,5)])
        myFlag = f-g
        self.assertAlmostEqual(myFlag,0,7, 'Eigs are wrong')

    def testFCC(self, num = 13, num_uc = 1):
        ats = [rf.atom() for i in range(num)]
        spin = rf.findmat(N.array([0,0,1]))
        # Positions
        #          0          1             2            3              4              
        pstn = [[0,0,0], [0.5,0,0.5], [-0.5,0,0.5], [0.5,0,-0.5], [-0.5,0,-0.5],
        #           5             6             7             8               9
                [0,0.5,0.5], [0,-0.5,0.5], [0,0.5,-0.5], [0,-0.5,-0.5], [0.5,0.5,0], 
        #             10            11            12
                [-0.5,0.5,0], [0.5,-0.5,0], [-0.5,-0.5,0]]
        # Neighbors
        #          0                               1             2            3              4              
        nbrs = [[1,2,3,4,5,6,7,8,9,10,11,12], [], [], [], [],
        #           5             6             7             8               9
                [], [], [], [], [], 
        #             10            11            12
                [], [], []]
        # Interactions
        #          0                               1             2            3              4              
        ints = [[0,1,2,3,4,5,6,7,8,9,10,11], [], [], [], [],
        #           5             6             7             8               9
                [], [], [], [], [], 
        #             10            11            12
                [], [], []]
        
        for i in range(len(ats)):
            ats[i].pos = N.array(pstn[i])
            ats[i].neighbors = nbrs[i]
            ats[i].interactions = ints[i]
            ats[i].spinRmatrix = spin
            ats[i].spinMagnitude = 1

        inters = [item for atom in ats for item in atom.interactions]
        max_inter = max(inters)
        
        J = sp.Symbol('J', real = True)
        Jij = [N.matrix([[J,0,0],[0,J,0],[0,0,J]]) for i in range(max_inter+1)]
        (Hsave, charpoly, eigs) = calculate_dispersion(ats, num_uc, num, Jij, showEigs = True)

        eigs = eigs.keys()
        for item in eigs:
            item = item.expand(mul = True, multinomial=True)

        testeig = (-384.0*self.J**2*self.S**2*sp.cos(0.5*self.kx)*sp.cos(0.5*self.ky) - 384.0*self.J**2*self.S**2*sp.cos(0.5*self.kx)*sp.cos(0.5*self.kz) - 384.0*self.J**2*self.S**2*sp.cos(0.5*self.ky)*sp.cos(0.5*self.kz) + 128.0*self.J**2*self.S**2*sp.cos(0.5*self.kx)**2*sp.cos(0.5*self.ky)*sp.cos(0.5*self.kz) + 128.0*self.J**2*self.S**2*sp.cos(0.5*self.ky)**2*sp.cos(0.5*self.kx)*sp.cos(0.5*self.kz) + 128.0*self.J**2*self.S**2*sp.cos(0.5*self.kz)**2*sp.cos(0.5*self.kx)*sp.cos(0.5*self.ky) + 576.0*self.J**2*self.S**2 + 64.0*self.J**2*self.S**2*sp.cos(0.5*self.kx)**2*sp.cos(0.5*self.ky)**2 + 64.0*self.J**2*self.S**2*sp.cos(0.5*self.kx)**2*sp.cos(0.5*self.kz)**2 + 64.0*self.J**2*self.S**2*sp.cos(0.5*self.ky)**2*sp.cos(0.5*self.kz)**2)**(0.5)

        self.assertAlmostEqual(Hsave.shape[0],2,4, 'H wrong shape')
        self.assertAlmostEqual(Hsave.shape[1],2,4, 'H wrong shape')
        self.assert_(len(eigs[0].args)==len(testeig.args), 'Eig.arg lengths are different sizes')
        f = sp.simplify(eigs[0].args[0].expand())
        f = f.subs([(self.J,1),(self.S,2),(self.kx,3),(self.ky,4),(self.kz,5)])
        g = sp.simplify(testeig.args[0].expand())
        g = g.subs([(self.J,1),(self.S,2),(self.kx,3),(self.ky,4),(self.kz,5)])
        myFlag = f-g
        self.assertAlmostEqual(myFlag,0,7, 'Eigs are wrong')

class TestAFM(unittest.TestCase):
    def setUp(self):
        self.S = sp.Symbol('S', real = True)
        self.J = sp.Symbol('J', real = True)
        self.kx = sp.Symbol('kx', real = True)
        self.ky = sp.Symbol('ky', real = True)
        self.kz = sp.Symbol('kz', real = True)

    def testChain(self, num = 3, num_uc = 1):
        ats = [rf.atom(pos = [i,0,0]) for i in range(num)]
        spin = rf.findmat(N.array([0,0,1]))
        aspin = rf.findmat(N.array([0,0,-1]))
        
        ints = [[0],[0,1],[1]]
        nbrs = [[1],[0,2],[1]]
        for i in range(num): 
            ats[i].interactions = ints[i]
            ats[i].neighbors = nbrs[i]
            ats[i].spinRmatrix = spin
            ats[i].spinMagnitude = 1
        for i in range(1,num,2): 
            ats[i].spinRmatrix = aspin
            ats[i].spinMagnitude = 1

        inters = [item for atom in ats for item in atom.interactions]
        max_inter = max(inters)
        
        J = sp.Symbol('J', real = True)
        Jij = [N.matrix([[J,0,0],[0,J,0],[0,0,J]]) for i in range(max_inter+1)]
        
        (Hsave, charpoly, eigs) = calculate_dispersion(ats, num_uc, num, Jij, showEigs = True)
        eigs = eigs.keys()
        for item in eigs:
            item = item.expand(mul = True, multinomial=True)
    
        testeig = (4*self.J**2*self.S**2 - 4.0*self.J**2*self.S**2*sp.cos(self.kx)**2)**(0.5)

        self.assertAlmostEqual(Hsave.shape[0],2,4, 'H wrong shape')
        self.assertAlmostEqual(Hsave.shape[1],2,4, 'H wrong shape')
        self.assert_(len(eigs[0].args)==len(testeig.args), 'Eig.arg lengths are different sizes')
        f = sp.simplify(eigs[0].args[0].expand())
        f = f.subs([(self.J,1),(self.S,2),(self.kx,3),(self.ky,4),(self.kz,5)])
        g = sp.simplify(testeig.args[0].expand())
        g = g.subs([(self.J,1),(self.S,2),(self.kx,3),(self.ky,4),(self.kz,5)])
        myFlag = f-g
        self.assertAlmostEqual(myFlag,0,7, 'Eigs are wrong')

    def testSquare(self, num = 5, num_uc = 1):
        ats = [rf.atom() for i in range(num)]
        spin = rf.findmat(N.array([0,0,1]))
        aspin = rf.findmat(N.array([0,0,-1]))
        # Positions
        #        0        1        2        3        4
        pstn = [[0,0,0],[0,1,0],[1,0,0],[0,-1,0],[-1,0,0]]
        # Neighbors
        #        0        1    2   3   4
        nbrs = [[1,2,3,4],[0],[0],[0],[0]]
        # Interactions
        #        0        1    2   3   4
        ints = [[0,1,2,3],[0],[1],[2],[3]]
    
        for i in range(len(ats)):
            ats[i].pos = N.array(pstn[i])
            ats[i].neighbors = nbrs[i]
            ats[i].interactions = ints[i]
            ats[i].spinRmatrix = spin
            ats[i].spinMagnitude = 1
        
        ats[0].spinMagnitude = 1
        ats[0].spinRmatrix = aspin

        inters = [item for atom in ats for item in atom.interactions]
        max_inter = max(inters)
        
        J = sp.Symbol('J', real = True)
        Jij = [N.matrix([[J,0,0],[0,J,0],[0,0,J]]) for i in range(max_inter+1)]
        (Hsave, charpoly, eigs) = calculate_dispersion(ats, num_uc, num, Jij, showEigs = True)

        eigs = eigs.keys()
        for item in eigs:
            item = item.expand(mul = True, multinomial=True)

        testeig = (-32.0*self.J**2*self.S**2*sp.cos(self.kx)*sp.cos(self.ky) + 64.0*self.J**2*self.S**2 - 16.0*self.J**2*self.S**2*sp.cos(self.kx)**2 - 16.0*self.J**2*self.S**2*sp.cos(self.ky)**2)**(0.5)

        self.assertAlmostEqual(Hsave.shape[0],2,4, 'H wrong shape')
        self.assertAlmostEqual(Hsave.shape[1],2,4, 'H wrong shape')
        self.assert_(len(eigs[1].args)==len(testeig.args), 'Eig.arg lengths are different sizes')
        f = sp.simplify(eigs[1].args[0].expand())
        f = f.subs([(self.J,1),(self.S,2),(self.kx,3),(self.ky,4),(self.kz,5)])
        g = sp.simplify(testeig.args[0].expand())
        g = g.subs([(self.J,1),(self.S,2),(self.kx,3),(self.ky,4),(self.kz,5)])
        myFlag = f-g
        self.assertAlmostEqual(myFlag,0,7, 'Eigs are wrong')

    def testCube(self, num = 7, num_uc = 1):
        ats = [rf.atom() for i in range(num)]
        spin = rf.findmat(N.array([0,0,1]))
        aspin = rf.findmat(N.array([0,0,-1]))
        # Positions
        #          0        1         2        3        4         5        6
        pstn = [[0,0,0], [1,0,0], [-1,0,0], [0,1,0], [0,-1,0], [0,0,1], [0,0,-1]]
        # Neighbors
        #            0          1    2    3    4    5    6
        nbrs = [[1,2,3,4,5,6], [0], [0], [0], [0], [0], [0]]
        # Interactions
        #            0          1    2    3    4    5    6
        ints = [[0,0,0,0,0,0], [1], [2], [3], [4], [5], [6]]
    
        for i in range(len(ats)):
            ats[i].pos = N.array(pstn[i])
            ats[i].neighbors = nbrs[i]
            ats[i].interactions = ints[i]
            ats[i].spinRmatrix = spin
            ats[i].spinMagnitude = 1

        ats[0].spinMagnitude = 1
        ats[0].spinRmatrix = aspin

        inters = [item for atom in ats for item in atom.interactions]
        max_inter = max(inters)
        
        J = sp.Symbol('J', real = True)
        Jij = [N.matrix([[J,0,0],[0,J,0],[0,0,J]]) for i in range(max_inter+1)]
        (Hsave, charpoly, eigs) = calculate_dispersion(ats, num_uc, num, Jij, showEigs = True)

        eigs = eigs.keys()
        for item in eigs:
            item = item.expand(mul = True, multinomial=True)

        testeig = (-32.0*self.J**2*self.S**2*sp.cos(self.kx)*sp.cos(self.ky) - 32.0*self.J**2*self.S**2*sp.cos(self.kx)*sp.cos(self.kz) - 32.0*self.J**2*self.S**2*sp.cos(self.ky)*sp.cos(self.kz) + 144.0*self.J**2*self.S**2 - 16.0*self.J**2*self.S**2*sp.cos(self.kx)**2 - 16.0*self.J**2*self.S**2*sp.cos(self.ky)**2 - 16.0*self.J**2*self.S**2*sp.cos(self.kz)**2)**(0.5)

        self.assertAlmostEqual(Hsave.shape[0],2,4, 'H wrong shape')
        self.assertAlmostEqual(Hsave.shape[1],2,4, 'H wrong shape')
        self.assert_(len(eigs[0].args)==len(testeig.args), 'Eig.arg lengths are different sizes')
        f = sp.simplify(eigs[0].args[0].expand())
        f = f.subs([(self.J,1),(self.S,2),(self.kx,3),(self.ky,4),(self.kz,5)])
        g = sp.simplify(testeig.args[0].expand())
        g = g.subs([(self.J,1),(self.S,2),(self.kx,3),(self.ky,4),(self.kz,5)])
        myFlag = f-g
        self.assertAlmostEqual(myFlag,0,7, 'Eigs are wrong')

    def testBCC(self, num = 9, num_uc = 1):
        ats = [rf.atom() for i in range(num)]
        spin = rf.findmat(N.array([0,0,1]))
        aspin = rf.findmat(N.array([0,0,-1]))
        # Positions
        #          0            1               2               3               4              
        pstn = [[0,0,0], [-0.5,0.5,-0.5], [0.5,0.5,-0.5], [0.5,0.5,0.5], [-0.5,0.5,0.5],
        #               5                6                7               8               
                [-0.5,-0.5,-0.5], [0.5,-0.5,-0.5], [0.5,-0.5,0.5], [-0.5,-0.5,0.5]]
        # Neighbors
        #          0               1   2   3   4              
        nbrs = [[1,2,3,4,5,6,7,8], [], [], [], [],
        #        5   6   7   8               
                [], [], [], []]
        # Interactions
        #          0               1   2   3   4              
        ints = [[0,1,2,3,4,5,6,7], [], [], [], [],
        #        5   6   7   8               
                [], [], [], []]

        for i in range(len(ats)):
            ats[i].pos = N.array(pstn[i])
            ats[i].neighbors = nbrs[i]
            ats[i].interactions = ints[i]
            ats[i].spinRmatrix = spin
            ats[i].spinMagnitude = 1
            
        ats[0].spinMagnitude = 1
        ats[0].spinRmatrix = aspin

        inters = [item for atom in ats for item in atom.interactions]
        max_inter = max(inters)
        
        J = sp.Symbol('J', real = True)
        Jij = [N.matrix([[J,0,0],[0,J,0],[0,0,J]]) for i in range(max_inter+1)]
        (Hsave, charpoly, eigs) = calculate_dispersion(ats, num_uc, num, Jij, showEigs = True)

        eigs = eigs.keys()
        for item in eigs:
            item = item.expand(mul = True, multinomial=True)

        testeig = (256.0*self.J**2*self.S**2 - 256.0*self.J**2*self.S**2*sp.cos(0.5*self.kx)**2*sp.cos(0.5*self.ky)**2*sp.cos(0.5*self.kz)**2)**(0.5)

        self.assertAlmostEqual(Hsave.shape[0],2,4, 'H wrong shape')
        self.assertAlmostEqual(Hsave.shape[1],2,4, 'H wrong shape')
        self.assert_(len(eigs[0].args)==len(testeig.args), 'Eig.arg lengths are different sizes')
        f = eigs[0].args[0].expand()
        f = f.subs([(self.J,1),(self.S,2),(self.kx,3),(self.ky,4),(self.kz,5),(sp.I,0)])
        g = testeig.args[0].expand()
        g = g.subs([(self.J,1),(self.S,2),(self.kx,3),(self.ky,4),(self.kz,5)])
        myFlag = f-g
        self.assertAlmostEqual(myFlag,0,7, 'Eigs are wrong')

    def testFCC(self, num = 13, num_uc = 1):
        ats = [rf.atom() for i in range(num)]
        spin = rf.findmat(N.array([0,0,1]))
        aspin = rf.findmat(N.array([0,0,-1]))
        # Positions
        #          0          1             2            3              4              
        pstn = [[0,0,0], [0.5,0,0.5], [-0.5,0,0.5], [0.5,0,-0.5], [-0.5,0,-0.5],
        #           5             6             7             8               9
                [0,0.5,0.5], [0,-0.5,0.5], [0,0.5,-0.5], [0,-0.5,-0.5], [0.5,0.5,0], 
        #             10            11            12
                [-0.5,0.5,0], [0.5,-0.5,0], [-0.5,-0.5,0]]
        # Neighbors
        #          0                               1             2            3              4              
        nbrs = [[1,2,3,4,5,6,7,8,9,10,11,12], [], [], [], [],
        #           5             6             7             8               9
                [], [], [], [], [], 
        #             10            11            12
                [], [], []]
        # Interactions
        #          0                               1             2            3              4              
        ints = [[0,1,2,3,4,5,6,7,8,9,10,11], [], [], [], [],
        #           5             6             7             8               9
                [], [], [], [], [], 
        #             10            11            12
                [], [], []]
        
        for i in range(len(ats)):
            ats[i].pos = N.array(pstn[i])
            ats[i].neighbors = nbrs[i]
            ats[i].interactions = ints[i]
            ats[i].spinRmatrix = spin
            ats[i].spinMagnitude = 1
        
        ats[0].spinMagnitude = 1
        ats[0].spinRmatrix = aspin

        inters = [item for atom in ats for item in atom.interactions]
        max_inter = max(inters)
        
        J = sp.Symbol('J', real = True)
        Jij = [N.matrix([[J,0,0],[0,J,0],[0,0,J]]) for i in range(max_inter+1)]
        (Hsave, charpoly, eigs) = calculate_dispersion(ats, num_uc, num, Jij, showEigs = True)

        eigs = eigs.keys()
        for item in eigs:
            item = item.expand(mul = True, multinomial=True)

        testeig = (-128.0*self.J**2*self.S**2*sp.cos(0.5*self.kx)**2*sp.cos(0.5*self.ky)*sp.cos(0.5*self.kz) - 128.0*self.J**2*self.S**2*sp.cos(0.5*self.ky)**2*sp.cos(0.5*self.kx)*sp.cos(0.5*self.kz) - 128.0*self.J**2*self.S**2*sp.cos(0.5*self.kz)**2*sp.cos(0.5*self.kx)*sp.cos(0.5*self.ky) + 576.0*self.J**2*self.S**2 - 64.0*self.J**2*self.S**2*sp.cos(0.5*self.kx)**2*sp.cos(0.5*self.ky)**2 - 64.0*self.J**2*self.S**2*sp.cos(0.5*self.kx)**2*sp.cos(0.5*self.kz)**2 - 64.0*self.J**2*self.S**2*sp.cos(0.5*self.ky)**2*sp.cos(0.5*self.kz)**2)**(0.5)

        self.assertAlmostEqual(Hsave.shape[0],2,4, 'H wrong shape')
        self.assertAlmostEqual(Hsave.shape[1],2,4, 'H wrong shape')
        self.assert_(len(eigs[1].args)==len(testeig.args), 'Eig.arg lengths are different sizes')
        f = sp.simplify(eigs[1].args[0].expand())
        f = f.subs([(self.J,1),(self.S,2),(self.kx,3),(self.ky,4),(self.kz,5)])
        g = sp.simplify(testeig.args[0].expand())
        g = g.subs([(self.J,1),(self.S,2),(self.kx,3),(self.ky,4),(self.kz,5)])
        myFlag = f-g
        self.assertAlmostEqual(myFlag,0,7, 'Eigs are wrong')

#------------------------------------------------------------------------------------------------------
#-----------------Non-testing Methods------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------

def gen_chain(num, anti = False):
    ats = [rf.atom(pos = [i,0,0]) for i in range(num)]
    spin = rf.findmat(N.array([0,0,1]))
    
    if num == 2: 
        ints = [[0],[0]]
        nbrs = [[1],[0]]
        for i in range(num): 
            ats[i].interactions = ints[i]
            ats[i].neighbors = nbrs[i]
            ats[i].spinRmatrix = spin
            ats[i].spinMagnitude = 1
    elif num == 3:
        ints = [[0],[0,1],[1]]
        nbrs = [[1],[0,2],[1]]
        for i in range(num): 
            ats[i].interactions = ints[i]
            ats[i].neighbors = nbrs[i]
            ats[i].spinRmatrix = spin
            ats[i].spinMagnitude = 1
    if anti == True:
        spin = rf.findmat(N.array([0,0,-1]))
        for i in range(1,num,2): 
            ats[i].spinRmatrix = spin
            ats[i].spinMagnitude = 1
    for i in range(len(ats)):
        print ats[i].spinRmatrix
    return ats

def square(anti = False):
    ats = [rf.atom() for i in range(5)]
    spin = rf.findmat(N.array([0,0,1]))
    # Positions
    #        0        1        2        3        4
    pstn = [[0,0,0],[0,1,0],[1,0,0],[0,-1,0],[-1,0,0]]
    # Neighbors
    #        0        1    2   3   4
    nbrs = [[1,2,3,4],[0],[0],[0],[0]]
    # Interactions
    #        0        1    2   3   4
    ints = [[0,1,2,3],[0],[1],[2],[3]]

    for i in range(len(ats)):
        ats[i].pos = N.array(pstn[i])
        ats[i].neighbors = nbrs[i]
        ats[i].interactions = ints[i]
        ats[i].spinRmatrix = spin
        ats[i].spinMagnitude = 1
        
    if anti == True:
        spin = rf.findmat(N.array([0,0,-1]))
        ats[0].spinMagnitude = 1
        ats[0].spinRmatrix = spin
    return ats

def gen_cube(body = False, face = False, anti = False):
    ats = []
    spin = rf.findmat(N.array([0,0,1]))
    # Cubic
    if not body and not face:

        ats = [rf.atom() for i in range(7)]
#        # Positions
#        #          0        1         2        3        4         5        6
#        pstn = [[0,0,0], [1,0,0], [-1,0,0], [0,1,0], [0,-1,0], [0,0,1], [0,0,-1]]
#
#        # Neighbors
#        #            0          1    2    3    4    5    6
#        nbrs = [[1,2,3,4,5,6], [0], [0], [0], [0], [0], [0]]
#
#        # Interactions
#        #            0          1    2    3    4    5    6
#        ints = [[0,1,2,3,4,5], [0], [1], [2], [3], [4], [5]]

        # Positions
        #          0        1         2        3        4         5        6
        pstn = [[0,0,0], [1,0,0], [-1,0,0], [0,1,0], [0,-1,0], [0,0,1], [0,0,-1]]

        # Neighbors
        #            0          1    2    3    4    5    6
        nbrs = [[1,2,3,4,5,6], [0], [0], [0], [0], [0], [0]]

        # Interactions
        #            0          1    2    3    4    5    6
        ints = [[0,0,0,0,0,0], [1], [2], [3], [4], [5], [6]]

        # Plug in values for positions and neighbors
        for i in range(len(ats)):
            ats[i].pos = N.array(pstn[i])
            ats[i].neighbors = nbrs[i]
            ats[i].interactions = ints[i]
            ats[i].spinRmatrix = spin
            ats[i].spinMagnitude = 1

    # BCC - Body Centered Cubic
    elif body and not face:   

        # Positions
        #          0            1               2               3               4              
        pstn = [[0,0,0], [-0.5,0.5,-0.5], [0.5,0.5,-0.5], [0.5,0.5,0.5], [-0.5,0.5,0.5],
        #               5                6                7               8               
                [-0.5,-0.5,-0.5], [0.5,-0.5,-0.5], [0.5,-0.5,0.5], [-0.5,-0.5,0.5]]

        # Neighbors
        #          0               1   2   3   4              
        nbrs = [[1,2,3,4,5,6,7,8], [], [], [], [],
        #        5   6   7   8               
                [], [], [], []]
        
        # Interactions
        #          0               1   2   3   4              
        ints = [[0,1,2,3,4,5,6,7], [], [], [], [],
        #        5   6   7   8               
                [], [], [], []]

        ats = [rf.atom() for i in range(len(pstn))]
        for i in range(len(ats)):
            ats[i].pos = N.array(pstn[i])
            ats[i].neighbors = nbrs[i]
            ats[i].interactions = ints[i]
            ats[i].spinRmatrix = spin
            ats[i].spinMagnitude = 1

    # FCC - Face Centered Cubic
    elif face and not body:
        # Positions
        #          0          1             2            3              4              
        pstn = [[0,0,0], [0.5,0,0.5], [-0.5,0,0.5], [0.5,0,-0.5], [-0.5,0,-0.5],
        #           5             6             7             8               9
                [0,0.5,0.5], [0,-0.5,0.5], [0,0.5,-0.5], [0,-0.5,-0.5], [0.5,0.5,0], 
        #             10            11            12
                [-0.5,0.5,0], [0.5,-0.5,0], [-0.5,-0.5,0]]
        # Neighbors
        #          0                               1             2            3              4              
        nbrs = [[1,2,3,4,5,6,7,8,9,10,11,12], [], [], [], [],
        #           5             6             7             8               9
                [], [], [], [], [], 
        #             10            11            12
                [], [], []]

        # Interactions
        #          0                               1             2            3              4              
        ints = [[0,1,2,3,4,5,6,7,8,9,10,11], [], [], [], [],
        #           5             6             7             8               9
                [], [], [], [], [], 
        #             10            11            12
                [], [], []]
        
        ats = [rf.atom() for i in range(len(pstn))]
        for i in range(len(ats)):
            ats[i].pos = N.array(pstn[i])
            ats[i].neighbors = nbrs[i]
            ats[i].interactions = ints[i]
            ats[i].spinRmatrix = spin
            ats[i].spinMagnitude = 1

    if anti == True:
        spin = rf.findmat(N.array([0,0,-1]))
        ats[0].spinMagnitude = 1
        ats[0].spinRmatrix = spin

    return ats

#calculate_dispersion(atom_list,N_atoms_uc,N_atoms,Jij,showEigs=False):
def run_dispersion(lattice, num_uc = 1, eigs = False):
    num = len(lattice)

    inters = [item for atom in lattice for item in atom.interactions]
    max_inter = max(inters)
    J = sp.Symbol('J', real = True)
    Jij = [N.matrix([[J,0,0],[0,J,0],[0,0,J]]) for i in range(max_inter+1)]
    if eigs == True:
        (Hsave, charpoly, eigs) = calculate_dispersion(lattice, num_uc, num, Jij, showEigs = eigs)
        print eigs
        eigs = eigs.keys()
        for eig in eigs:
            print eig.args
            eig = sp.simplify(eig.expand(deep = True, mul = True, multinomial=True))
        return (Hsave, charpoly, eigs)
    else: 
        Hsave = calculate_dispersion(lattice, num_uc, num, Jij , showEigs =eigs)
        return Hsave

def run_dispersion_from_file(interactionfile,spinfile, eigs=False):
    atom_list, jnums, jmats,N_atoms_uc=rf.readFiles(interactionfile,spinfile)
    N_atoms=len(atom_list)
    Hsave = calculate_dispersion(atom_list,N_atoms_uc,N_atoms,jmats,showEigs=eigs)
    if eigs == True:
        temp = 1
        x = sp.Symbol('x')
        for i in range(Hsave.cols): temp *= x-Hsave[i,i]
        poly = sp.simplify(temp.expand())
        print poly
    return sp.simplify(Hsave[0])

def print_it(it):
    if not isinstance(it, tuple): print
    if len(it) != 3: print 
    print '\nH    = ', it[0]
    print 'poly = ', it[1]
    print 'eigs = ', it[2]

#------------------------------------------------------------------------------------------------------
#-----------------M A I N------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------

if __name__ == "__main__":
    # Run Tests
    if 0:
        unittest.main()    
    
    # Run Dispersion from files
    if 1:   
        interfile = 'C:/montecarlo_cube.txt'
        spinfile = 'C:/Spins_cube.txt'    
        
        run_dispersion_from_file(interfile,spinfile)
    
    # Ferromagnetic Cases
    if 0:
        chain3 = gen_chain(3)
        sq = square()
        cube = gen_cube()
        bcc = gen_cube(body = True)
        fcc = gen_cube(face = True)
        
        Hchain3 = run_dispersion(chain3, num_uc = 1, eigs=True)
        Hsq = run_dispersion(sq, num_uc = 1, eigs=True)
        Hcube = run_dispersion(cube, num_uc = 1, eigs=True)
        Hbcc = run_dispersion(bcc, num_uc = 1, eigs=True)
        Hfcc = run_dispersion(fcc, num_uc = 1, eigs=True)
    
        print '\nFerro'
        print '\nchain3', 
        print_it(Hchain3)
        print '\nsq    ', 
        print_it(Hsq)
        print '\ncube  ', 
        print_it(Hcube)
        print '\nbcc   ', 
        print_it(Hbcc)
        print '\nfcc   ', 
        print_it(Hfcc)
        raw_input("\nPress <enter> to Continue\n")
    
    # Antiferromagnetic Cases
    if 0:
        achain3 = gen_chain(3, anti=True)
        asq = square(anti=True)
        acube = gen_cube(anti=True)
        abcc = gen_cube(body = True, anti=True)
        afcc = gen_cube(face = True, anti=True)
        
        aHchain3 = run_dispersion(achain3, num_uc = 1, eigs=True)
        aHsq = run_dispersion(asq, num_uc = 1, eigs=True)
        aHcube = run_dispersion(acube, num_uc = 1, eigs=True)
        aHbcc = run_dispersion(abcc, num_uc = 1, eigs=True)
        aHfcc = run_dispersion(afcc, num_uc = 1, eigs=True)
    
        print '\nAntiFerro'
        print '\nachain3', 
        print_it(aHchain3)
        print '\nasq    ', 
        print_it(aHsq)
        print '\nacube  ', 
        print_it(aHcube)
        print '\nabcc   ', 
        print_it(aHbcc)
        print '\nafcc   ', 
        print_it(aHfcc)
