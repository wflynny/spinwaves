import wx
import wxaddons.sized_controls as sc
import  wx.lib.intctrl
import  wx.grid as  gridlib
import numpy as N
import sys,os
import spinwave_calc_file as spinwave_calc_file
import wx.richtext
from sympy import pi
import spinwaves.cross_section.util.printing as printing
from multiprocessing import Process, Pipe
import copy
from spinwaves.utilities.Processes import ProcessManager

class MyApp(wx.App):
    def __init__(self, redirect=False, filename=None, useBestVisual=False, clearSigInt=True):
        wx.App.__init__(self,redirect,filename,clearSigInt)


    def OnInit(self):
        return True
    
    
class mFormValidator(wx.PyValidator):
    def __init__(self,data,key):
        wx.PyValidator.__init__(self)
        #print 'FormValidator init', key
        #self.TransferToWindow()
        self.data=data
        self.key=key
    def Clone(self):
        return mFormValidator(self.data,self.key)

    def Validate(self,win):
        print 'validating'
        textCtrl=self.GetWindow()
        text=textCtrl.GetValue()
        
        try:
            value=float(text)
            textCtrl.SetBackgroundColour(wx.SystemSettings_GetColour(wx.SYS_COLOUR_WINDOW))
            textCtrl.Refresh()
            return True
            
        except ValueError:
            wx.MessageBox("This field must be a number","error")
            textCtrl.SetBackgroundColour("pink")
            textCtrl.SetFocus()
            textCtrl.Refresh()
            return False
        

    def TransferToWindow(self):
        print 'Form TransferToWindow',self.key
        
        textCtrl=self.GetWindow()
        ##print 'checkctrl',checkctrl
        print self.__dict__

        textCtrl.SetValue(self.data.get(self.key,""))
        return True

    def TransferFromWindow(self):
        print 'TransferFromWindow'
        textCtrl=self.GetWindow()
        self.data[self.key]=float(textCtrl.GetValue())
        #self.qfloat=float(textCtrl.GetValue())
        return True    
    
    
def WalkTree(parent):
        print 'walking', parent
        parent.SetExtraStyle(wx.WS_EX_VALIDATE_RECURSIVELY)
        for child in parent.GetChildren():
            if child==None:
                print 'child',child.size
            else:
                WalkTree(child)

    
class RichTextFrame(wx.Frame):
    def __init__(self, *args, **kw):
        wx.Frame.__init__(self, *args, **kw)

     #   self.MakeMenuBar()
     #   self.MakeToolBar()
     #   self.CreateStatusBar()
     #   self.SetStatusText("Welcome to wx.richtext.RichTextCtrl!")
        sizer = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(sizer)
        self.spinsRtc = wx.richtext.RichTextCtrl(self, style=wx.VSCROLL|wx.HSCROLL|wx.NO_BORDER, size = (620,120));
        self.interactionsRtc = wx.richtext.RichTextCtrl(self, style=wx.VSCROLL|wx.HSCROLL|wx.NO_BORDER, size = (620,120));
     #   wx.CallAfter(self.rtc.SetFocus)
        interactionsTitle = wx.StaticText(self, -1, " Interaction File ")
        spinsTitle = wx.StaticText(self, -1, " Spins File ")
        
        
        self.interactionSave = wx.Button(self, -1, "Save", size = (50,20))
        interactionSizer = wx.BoxSizer(wx.HORIZONTAL)
        interactionSizer.Add(interactionsTitle)
        interactionSizer.Add(self.interactionSave)
        
        self.spinSave = wx.Button(self, -1, "Save", size = (50,20))
        spinSizer = wx.BoxSizer(wx.HORIZONTAL)
        spinSizer.Add(spinsTitle)
        spinSizer.Add(self.spinSave)
        
        
        sizer.Add(interactionSizer)
        sizer.Add(self.interactionsRtc)
        sizer.Add(spinSizer)
        sizer.Add(self.spinsRtc)
        
        self.Bind(wx.EVT_BUTTON, self.OnSaveInteractions, self.interactionSave)
        self.Bind(wx.EVT_BUTTON, self.OnSaveSpins, self.spinSave)
        
        self.Fit()
        self.SetMinSize(self.GetSize())
        
    
    def OnSaveInteractions(self, evt):
        self.interactionsRtc.SaveFile()
    
    def OnSaveSpins(self, evt):
        self.spinsRtc.SaveFile()
    
    def loadSpins(self, file):
        self.spinsRtc.LoadFile(file)
    
    def loadInteractions(self, file):
        self.interactionsRtc.LoadFile(file)



class FormDialog(sc.SizedPanel):
    def __init__(self, parent, id, procManager):
        self.parent = parent
        #self.process_list = []
        self.processManager = procManager
        
        #valstyle=wx.WS_EX_VALIDATE_RECURSIVELY
        sc.SizedPanel.__init__(self, parent, -1,
                        style= wx.RESIZE_BORDER)#| wx.WS_EX_VALIDATE_RECURSIVELY)
        self.SetExtraStyle(wx.WS_EX_VALIDATE_RECURSIVELY)
        pane = self
        pane.SetSizerType("vertical")
               
        print 'FormDialog called'
        FilePane = sc.SizedPanel(pane, -1)
        FilePane.SetSizerType("vertical")
        FilePane.SetSizerProps(expand=True)
   
        self.interactionfile=''
        self.interactionfilectrl=wx.StaticText(FilePane, -1, "Interaction File:%s"%(self.interactionfile,))
        b = wx.Button(FilePane, -1, "Browse")
        self.Bind(wx.EVT_BUTTON, self.OnOpenInt, b)
        
        self.spinfile=''
        self.spinfilectrl=wx.StaticText(FilePane, -1, "Spin File:%s"%(self.spinfile,))
        b = wx.Button(FilePane, -1, "Browse")
        self.Bind(wx.EVT_BUTTON, self.OnOpen, b)
        



        self.data={}
        self.data['kx']=str(0.0)
        self.data['ky']=str(0.0)
        self.data['kz']=str(1.0)
        DirectionPane = sc.SizedPanel(pane, -1)
        DirectionPane.SetSizerType("vertical")
        DirectionPane.SetSizerProps(expand=True)
        wx.StaticText(DirectionPane, -1, "Scan Direction")
        
        DirectionsubPane = sc.SizedPanel(pane, -1)
        DirectionsubPane.SetSizerType("horizontal")
        DirectionsubPane.SetSizerProps(expand=True)

        

        DirectionsubPane.SetExtraStyle(wx.WS_EX_VALIDATE_RECURSIVELY)
        WalkTree(pane)
        
        wx.StaticText(DirectionsubPane, -1, "qx")
        qx=wx.TextCtrl(DirectionsubPane, -1, self.data['kx'],validator=mFormValidator(self.data,'kx'))
        #print 'qx', qx
        #self.Bind(wx.EVT_TEXT, self.Evtqx, qx)
        wx.StaticText(DirectionsubPane, -1, "qy")
        qy=wx.TextCtrl(DirectionsubPane, -1, self.data['ky'],validator=mFormValidator(self.data,'ky'))
        #self.Bind(wx.EVT_TEXT, self.Evtqx, qy)
        wx.StaticText(DirectionsubPane, -1, "qz")
        qz=wx.TextCtrl(DirectionsubPane, -1, self.data['kz'],validator=mFormValidator(self.data,'kz'))
        #self.Bind(wx.EVT_TEXT, self.Evtqx, qz)
        #print 'Directed'


        wx.StaticText(DirectionsubPane, -1, "Number of divisions")
        self.data['step']=8
        spinctrl = wx.SpinCtrl(DirectionsubPane, -1, "", (30, 50))
        spinctrl.SetRange(1,100)
        spinctrl.SetValue(self.data['step'])
        self.Bind(wx.EVT_SPINCTRL,self.EvtSpinCtrl)
        self.spinctrl=spinctrl
        
        
        #Add Range
        self.kRange = {}
        self.kRange['kMin'] = str(0)
        self.kRange['kMax'] = str(2)
        
        wx.StaticText(pane, -1, "k Range")
        kRangePane = sc.SizedPanel(pane, -1)
        kRangePane.SetSizerType("horizontal")
        
        wx.StaticText(kRangePane, -1, "Min =")
        self.kMinCtrl = wx.TextCtrl(kRangePane, -1, self.kRange['kMin'],validator=mFormValidator(self.kRange,'kMin'))
        wx.StaticText(kRangePane, -1, "*pi")
        
        wx.StaticText(kRangePane, -1, "    Max =")
        self.kMaxCtrl = wx.TextCtrl(kRangePane, -1, self.kRange['kMax'],validator=mFormValidator(self.kRange,'kMax'))
        wx.StaticText(kRangePane, -1, "*pi")
        
        
        
        # add dialog buttons
        
        #Getting rid of standard button function
        #self.SetButtonSizer(self.CreateStdDialogButtonSizer(wx.OK | wx.CANCEL))
        
        btnPane = sc.SizedPanel(pane, -1)
        btnPane.SetSizerType("horizontal")
        btnPane.SetSizerProps(expand=True)
        
        self.btnOk = wx.Button(btnPane,-1, "Ok")
        #self.btnOk.SetDefault()
        #btnsizer.Add(self.btnOk)

        btnCancel = wx.Button(btnPane, -1, "Cancel")
        #btnsizer.Add(btn)
        
        #self.GetSizer().Add(btnsizer)
        
        #intercept OK button click event
        self.Bind(wx.EVT_BUTTON, self.OnOk, self.btnOk)
        
        #intercept Cancel button click event
        self.Bind(wx.EVT_BUTTON, self.OnCancel, btnCancel)
        

        
        #self.TransferDataToWindow()
        # a little trick to make sure that you can't resize the dialog to
        # less screen space than the controls need
        self.Fit()
        size = self.GetSize()
        size = (size[0]+5, size[1]+35)#I think borders may not be included because size is consisitantly too small
        self.parent.SetSize(size)
        self.parent.SetMinSize(size)
        #This trick is not working because the self.GetSize() is too small
        #self.parent.SetMinSize((645, 245))
        #self.parent.SetMinSize((400,400))
        
           
                
        #Text editor
        self.editorWin = RichTextFrame(self, -1, "Editor",
                            size=(620, 250),
                            style = wx.DEFAULT_FRAME_STYLE)
        self.editorWin.Show(True)

        
    def OnCancel(self, evt):
        """Closes window"""
        for pro in self.process_list:
            pro.terminate()
        self.parent.Close()
        
    
    
    def OnOk(self, event):
        #Formerly, when this was modal, this would be done in calling function
        self.Fit()
        self.Validate()
        print "OK"
        self.TransferDataFromWindow()
        print 'data',self.data
        print self.data['step']
        print self.interactionfile
        print self.spinfile

        self.processManager.startAnalyticDispersion(self.interactionfile, self.spinfile)
        self.processManager.startNumericDispersion(self.interactionfile, self.spinfile, self.data, float(self.kRange['kMin'])*pi, float(self.kRange['kMax'])*pi, self.data['step'])
        

        
        #Hsave = spinwave_calc_file.driver1(self.spinfile,self.interactionfile)
        #myeigs=printing.eig_process(copy.deepcopy(Hsave))
        
        
    def EvtSpinCtrl(self,evt):
        print self.__dict__
        print 'event',evt.__dict__
        #spinCtrl=evt.control
        self.data['step']=self.spinctrl.GetValue()
        #sc.SizedDialog.GetW
        #self.G
        #self.kx=text
 


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
        confBase = wx.ConfigBase.Create()
        confBase.SetStyle(wx.CONFIG_USE_LOCAL_FILE)
        defaultDir=confBase.Get().GetPath()
        wildcard="files (*.txt)|*.txt|All files (*.*)|*.*"
        dlg = wx.FileDialog(
            self, message="Choose a spin configuration file",
            defaultDir=defaultDir,
            defaultFile="",
            wildcard=wildcard,
            style=wx.OPEN | wx.CHANGE_DIR
            )

        # Show the dialog and retrieve the user response. If it is the OK response,
        # process the data.
        if dlg.ShowModal() == wx.ID_OK:
            # This returns a Python list of files that were selected.
            paths = dlg.GetPaths()
            #self.log.WriteText('You selected %d files:' % len(paths))
            self.spinfile=paths[0].encode('ascii')
            self.spinfilectrl.SetLabel("Spin File:%s"%(self.spinfile,))
            #display in richtextcontrol
            #self.spinsRtc.LoadFile(paths[0])
            self.editorWin.loadSpins(paths[0])
            
            #wx.StaticText(FilePane, -1, "CellFile:%s"%(self.groupdata['cellfile'],))
        # Destroy the dialog. Don't do this until you are done with it!
        # BAD things can happen otherwise!
        
        dlg.Destroy()


    def OnOpenInt(self,event):
        # Create the dialog. In this case the current directory is forced as the starting
        # directory for the dialog, and no default file name is forced. This can easilly
        # be changed in your program. This is an 'open' dialog, and allows multitple
        # file selections as well.
        #
        # Finally, if the directory is changed in the process of getting files, this
        # dialog is set up to change the current working directory to the path chosen.

        #defaultDir=os.getcwd()
        #defaultDir=r'C:\polcorrecter\data'
        confBase = wx.ConfigBase.Create()
        confBase.SetStyle(wx.CONFIG_USE_LOCAL_FILE)
        defaultDir=confBase.Get().GetPath()
        #defaultDir=wx.ConfigBase.Get().GetPath()
        wildcard="files (*.txt)|*.txt|All files (*.*)|*.*"
        dlg = wx.FileDialog(
            self, message="Choose an Interaction file",
            defaultDir=defaultDir,
            defaultFile="",
            wildcard=wildcard,
            style=wx.OPEN | wx.CHANGE_DIR
            )

        # Show the dialog and retrieve the user response. If it is the OK response,
        # process the data.
        if dlg.ShowModal() == wx.ID_OK:
            # This returns a Python list of files that were selected.
            paths = dlg.GetPaths()
            #self.log.WriteText('You selected %d files:' % len(paths))
            self.interactionfile=paths[0].encode('ascii')
            self.interactionfilectrl.SetLabel("Interaction File:%s"%(self.interactionfile,))
            #wx.StaticText(FilePane, -1, "CellFile:%s"%(self.groupdata['cellfile'],))
            
            #display in richtextcontrol
            #self.interactionsRtc.LoadFile(paths[0])
            self.editorWin.loadInteractions(paths[0])
            
        # Destroy the dialog. Don't do this until you are done with it!
        # BAD things can happen otherwise!
        dlg.Destroy()


if __name__=='__main__':
    app=MyApp()
    frame1 = wx.Frame(None, -1, "Spinwaves")
    dlg=FormDialog(parent=frame1,id=-1, procManager = ProcessManager(frame1))
    frame1.Show()
    if 0:
        frame1 = wx.Frame(None, -1, "Spinwaves")
        dlg=FormDialog(parent=frame1,id=-1)
        frame1.Show()
        result=dlg.ShowModal()
        if result==wx.ID_OK:
            dlg.Validate()
            print "OK"
            dlg.TransferDataFromWindow()
            print dlg.data
        else:
            print "Cancel"
    
        dlg.Destroy()

    app.MainLoop()
