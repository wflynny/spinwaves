import os
from multiprocessing import Process, Queue
from threading import Thread
from copy import deepcopy
import wx
#from wx.py.dispatcher import send
import numpy as np
import matplotlib.pyplot as plt
import spinwaves.spinwavecalc.spinwave_calc_file as spinwave_calc_file
import spinwaves.cross_section.util.printing as printing
from spinwaves.MonteCarlo.CSim import ShowSimulationFrame
from spinwaves.vtkModel.wxGUI.Session import Session


    
class ProcessManager():
    def __init__(self, parentWindow):
        """All results will be displayed in windows which are direct children of parentWindow.  This should probably be the main GUI Frame."""
        self.parent = parentWindow
        #Threads to poll the process Queues
        self._analyticDispThread = None
        self._numericDispThread = None
        self._fitThread = None
        self._analyticCrossSecThread = None
        
        self._analyticDispQueue = Queue()
        self._numericDispQueue = Queue()
        self._fitQueue = Queue()
        self._analyticCrossSecQueue = Queue()
        
        self._analyticDispProcesses= []
        self._numericDispProcesses= []
        self._fitProcesses= []
        self._analyticCrossSecProcesses= []
        
    def startAnalyticDispersion(self, interaction_file, spin_file):
        p = Process(target=AnalyticDispFunc, args=(self._analyticDispQueue, interaction_file, spin_file))
        self._analyticDispProcesses.append(p)
        p.start()
        #AnalyticDispFunc(self._analyticDispQueue, interaction_file, spin_file)
        if self._analyticDispThread == None:
            self._analyticDispThread = AnalyticDispersionThread(self.parent, self._analyticDispQueue)
            self._analyticDispThread.start()
       
    
    def startNumericDispersion(self, interaction_file, spin_file, direction, k_min, k_max, steps): 
        p = Process(target = NumericDispFunc, args = (self._numericDispQueue, interaction_file, spin_file, direction, k_min, k_max, steps))
        self._numericDispProcesses.append(p)
        p.start()
        #NumericDispFunc(self._numericDispQueue, interaction_file, spin_file, direction, k_min, k_max, steps)
        if self._numericDispThread == None:
            self._numericDispThread = NumericDispersionThread(self.parent, self._numericDispQueue)
            self._numericDispThread.start()

    
    def startFit(self, session, fileName, k, tmin, tmax, size, tfactor, useMC):
        #bondTable = deepcopy(session.bondTable)
        #It's easier just to copy the session although it is less memory efficient
        #sess = deepcopy(session)
        #wx objects are not picklable, so only the data list from the bondata will be sent
        #bondTable = session.bondTable.__deepcopy__()
        sessXML = session.createXMLStr()
        #FitFunc(self._fitQueue, sessXML, fileName, k, tmin, tmax, size, tfactor, useMC)
        p = Process(target = FitFunc, args = (self._fitQueue, sessXML, fileName, k, tmin, tmax, size, tfactor, useMC))
        self._fitProcesses.append(p)
        p.start()
        if self._fitThread == None:
            self._fitThread = FitThread(self.parent, self._fitQueue)
            self._fitThread.start()
        
        
    def startAnalyticCrossSection(self):
        print "method stub"

        
class AnalyticDispersionThread(Thread):
   def __init__ (self, parentWindow, queue):
      Thread.__init__(self)
      self.parent = parentWindow
      self.queue = queue
    
   def run(self):
       while(True):
           ans = self.queue.get()
           wx.CallAfter(showAnalyticEigs, ans)
           #send(signal = "Analytic Dispersion Complete", answer = ans)
           #evt = AnalyticDispCompleteEvent(myAnalyticDispEvt, None)
           #evt.SetAns(ans)
           
           #eig_frame = printing.LaTeXDisplayFrame(self.parent, ans, 'Dispersion Eigenvalues')
           #eig_frame.Show()
           
def showAnalyticEigs(ans):
    eig_frame = printing.LaTeXDisplayFrame(None, ans, 'Dispersion Eigenvalues')
    eig_frame.Show()
    

      
class NumericDispersionThread(Thread):
    def __init__ (self, parentWindow, queue):
      Thread.__init__(self)
      self.parent = parentWindow
      self.queue = queue
      
    def run(self):
        while(True):
            ans = self.queue.get()
            qrange = ans[0]
            wrange = ans[1]
            fig = plt.figure()
            ax = fig.add_subplot(111)
            for wrange1 in wrange:
                ax.plot(qrange, wrange1)
                plt.hold(True)
            plt.show()

from spinwaves.utilities.fitting import showFitResultFrame, fitFromFile
class FitThread(Thread):
    def __init__ (self, parentWindow, queue):
       Thread.__init__(self)
       self.parent = parentWindow
       self.queue = queue
    
    def run(self):
        #poll the fit queue while there are still fitting processes running
        #while(len(self._fitProcesses)>0):#This could also be infinite
        while(True):
            result = self.queue.get()
            data = result[1]
            pid = result[0]
            wx.CallAfter(showFitResultFrame,data, pid)




class AnalyticCrossSectionThread(Thread):
    def __init__ (self):
        Thread.__init__(self)
      
    def run(self):
        print "method stub"
        
def AnalyticDispFunc(queue, int_file, spin_file):
    #Since calculating Hsave is most of what the numeric process does, it might be better not to do this twice.
    Hsave = spinwave_calc_file.driver1(spin_file, int_file)
    myeigs=printing.eig_process(deepcopy(Hsave))
    queue.put(printing.create_latex(myeigs, "eigs"))
        
def NumericDispFunc(queue, int_file, spin_file, direction, k_min, k_max, steps):
    Hsave = spinwave_calc_file.driver1(spin_file, int_file)
    qrange, wranges = spinwave_calc_file.driver2(Hsave, direction, steps, k_min, k_max)
    queue.put((qrange, wranges))

def FitFunc(queue, sessionXML, fileName, k, tmin, tmax, size, tfactor, useMC):
    sess = Session()
    sess.loadXMLStr(sessionXML)
    ans = fitFromFile(fileName, sess, k = k, tMax = tmax, tMin = tmin, size = size, tFactor = tfactor, MCeveryTime = useMC)
    #The session which contains the fitted parameters is more useful
    pid = os.getpid()
    queue.put((pid, sess.bondTable.data))