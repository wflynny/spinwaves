import wx
#import wxaddons.sized_controls as sc
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

    
class RichTextPanel(wx.Panel):
    def __init__(self, allowEdit, *args, **kw):
        wx.Panel.__init__(self, *args, **kw)

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
        
        if not allowEdit:
            self.spinSave.Enable(False)
            self.interactionSave.Enable(False)
        
    
    def OnSaveInteractions(self, evt):
        self.interactionsRtc.SaveFile()
    
    def OnSaveSpins(self, evt):
        self.spinsRtc.SaveFile()
    
    def loadSpins(self, file):
        print "\n\n\n\nfilename: " , file
        self.spinsRtc.LoadFile(file)
    
    def loadInteractions(self, file):
        print "\n\n\n\nfilename: " , file
        self.interactionsRtc.LoadFile(file)

def showEditorWindow(parent, title, allowEditting = True):
    """Creates and displays a simple frame containing the RichTextPanel."""
    frame = wx.Frame(parent, -1, title, size=(630, 320), style = wx.DEFAULT_FRAME_STYLE)
    panel = RichTextPanel(allowEditting, frame, -1)
    #frame.Fit()
    #frame.SetMinSize(frame.GetSize())
    frame.Show()
    return panel



class SpinwavePanel(wx.Panel):
    def __init__(self, procManager, *args, **kwds):
        # begin wxGlade: SpinwavePanel.__init__
        kwds["style"] = wx.TAB_TRAVERSAL
        wx.Panel.__init__(self, *args, **kwds)
        self.sizer_7_staticbox = wx.StaticBox(self, -1, "Scan Direction")
        self.int_file_label = wx.StaticText(self, -1, " Interaction File: ")
        self.int_file_txtCtrl = wx.TextCtrl(self, -1, "")
        self.int_browse_btn = wx.Button(self, -1, "Browse")
        self.spin_file_label = wx.StaticText(self, -1, " Spin File: ")
        self.spin_file_txtCtrl = wx.TextCtrl(self, -1, "")
        self.spin_browse_btn = wx.Button(self, -1, "Browse")
        self.static_line_1 = wx.StaticLine(self, -1)
        self.kMin_label = wx.StaticText(self, -1, "k Min:")
        self.kMin_txtCtrl = wx.TextCtrl(self, -1, "")
        self.label_9 = wx.StaticText(self, -1, "*pi")
        self.kMax_label = wx.StaticText(self, -1, "k Max:")
        self.kMax_txtCtrl = wx.TextCtrl(self, -1, "")
        self.label_10 = wx.StaticText(self, -1, "*pi")
        self.label_7 = wx.StaticText(self, -1, "Steps:")
        self.steps_spin_ctrl = wx.SpinCtrl(self, -1, "", min=0, max=100)
        self.kx_label = wx.StaticText(self, -1, "kx:")
        self.kx_txtCtrl = wx.TextCtrl(self, -1, "")
        self.ky_label = wx.StaticText(self, -1, "ky:")
        self.ky_txtCtrl = wx.TextCtrl(self, -1, "")
        self.kz_label = wx.StaticText(self, -1, "kz:")
        self.kz_txtCtrl = wx.TextCtrl(self, -1, "")
        self.ok_btn = wx.Button(self, -1, "Ok")
        self.cancel_btn = wx.Button(self, -1, "Cancel")

        self.__set_properties()
        self.__do_layout()

        self.Bind(wx.EVT_BUTTON, self.OnIntFileBrowse, self.int_browse_btn)
        self.Bind(wx.EVT_BUTTON, self.OnSpinFileBrowse, self.spin_browse_btn)
        self.Bind(wx.EVT_BUTTON, self.OnOk, self.ok_btn)
        self.Bind(wx.EVT_BUTTON, self.OnCancel, self.cancel_btn)
        # end wxGlade
        self.editorWin = showEditorWindow(self, "Spinwave File Editor")
        self.processManager = procManager

    def __set_properties(self):
        # begin wxGlade: SpinwavePanel.__set_properties
        self.int_file_txtCtrl.SetMinSize((180, 27))
        self.spin_file_txtCtrl.SetMinSize((180, 27))
        # end wxGlade
        
        self.kx_txtCtrl.SetValue(str(0))
        self.ky_txtCtrl.SetValue(str(0))
        self.kz_txtCtrl.SetValue(str(1))
        self.kMin_txtCtrl.SetValue(str(0))
        self.kMax_txtCtrl.SetValue(str(2))

    def __do_layout(self):
        # begin wxGlade: SpinwavePanel.__do_layout
        sizer_4 = wx.BoxSizer(wx.VERTICAL)
        sizer_12 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_7 = wx.StaticBoxSizer(self.sizer_7_staticbox, wx.HORIZONTAL)
        sizer_10 = wx.BoxSizer(wx.VERTICAL)
        sizer_9 = wx.BoxSizer(wx.VERTICAL)
        sizer_8 = wx.BoxSizer(wx.VERTICAL)
        sizer_13 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_11 = wx.BoxSizer(wx.VERTICAL)
        sizer_21 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_20 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_15 = wx.BoxSizer(wx.VERTICAL)
        sizer_17 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_19 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_14 = wx.BoxSizer(wx.VERTICAL)
        sizer_16 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_18 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_6 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_5 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_5.Add(self.int_file_label, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        sizer_5.Add(self.int_file_txtCtrl, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        sizer_5.Add(self.int_browse_btn, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        sizer_4.Add(sizer_5, 2, wx.EXPAND, 0)
        sizer_6.Add(self.spin_file_label, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        sizer_6.Add(self.spin_file_txtCtrl, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        sizer_6.Add(self.spin_browse_btn, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        sizer_4.Add(sizer_6, 2, wx.EXPAND, 0)
        sizer_4.Add(self.static_line_1, 0, wx.EXPAND, 0)
        sizer_4.Add((15, 15), 0, 0, 0)
        sizer_13.Add((15, 15), 0, 0, 0)
        sizer_18.Add(self.kMin_label, 0, wx.ALIGN_BOTTOM|wx.ALIGN_CENTER_HORIZONTAL, 0)
        sizer_14.Add(sizer_18, 1, wx.EXPAND, 0)
        sizer_16.Add(self.kMin_txtCtrl, 0, wx.ALIGN_CENTER_HORIZONTAL|wx.ALIGN_CENTER_VERTICAL, 0)
        sizer_16.Add(self.label_9, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        sizer_14.Add(sizer_16, 1, wx.EXPAND, 0)
        sizer_13.Add(sizer_14, 1, wx.EXPAND, 0)
        sizer_19.Add(self.kMax_label, 0, wx.ALIGN_BOTTOM|wx.ALIGN_CENTER_HORIZONTAL, 0)
        sizer_15.Add(sizer_19, 1, wx.EXPAND, 0)
        sizer_17.Add(self.kMax_txtCtrl, 0, wx.ALIGN_CENTER_HORIZONTAL|wx.ALIGN_CENTER_VERTICAL, 0)
        sizer_17.Add(self.label_10, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        sizer_15.Add(sizer_17, 1, wx.EXPAND, 0)
        sizer_13.Add(sizer_15, 1, wx.EXPAND, 0)
        sizer_20.Add(self.label_7, 0, wx.ALIGN_BOTTOM|wx.ALIGN_CENTER_HORIZONTAL, 0)
        sizer_11.Add(sizer_20, 1, wx.EXPAND, 0)
        sizer_21.Add(self.steps_spin_ctrl, 0, wx.ALIGN_CENTER_HORIZONTAL|wx.ALIGN_CENTER_VERTICAL, 0)
        sizer_11.Add(sizer_21, 1, wx.EXPAND, 0)
        sizer_13.Add(sizer_11, 1, wx.EXPAND, 0)
        sizer_4.Add(sizer_13, 3, wx.EXPAND, 0)
        sizer_4.Add((15, 15), 0, 0, 0)
        sizer_8.Add(self.kx_label, 0, wx.ALIGN_CENTER_HORIZONTAL, 0)
        sizer_8.Add(self.kx_txtCtrl, 0, wx.ALIGN_CENTER_HORIZONTAL, 0)
        sizer_7.Add(sizer_8, 1, wx.EXPAND, 0)
        sizer_9.Add(self.ky_label, 0, wx.ALIGN_CENTER_HORIZONTAL, 0)
        sizer_9.Add(self.ky_txtCtrl, 0, wx.ALIGN_CENTER_HORIZONTAL, 0)
        sizer_7.Add(sizer_9, 1, wx.EXPAND, 0)
        sizer_10.Add(self.kz_label, 0, wx.ALIGN_CENTER_HORIZONTAL, 0)
        sizer_10.Add(self.kz_txtCtrl, 0, wx.ALIGN_CENTER_HORIZONTAL, 0)
        sizer_7.Add(sizer_10, 1, wx.EXPAND, 0)
        sizer_4.Add(sizer_7, 3, wx.EXPAND, 0)
        sizer_12.Add(self.ok_btn, 0, wx.ALIGN_BOTTOM|wx.ALIGN_CENTER_HORIZONTAL, 0)
        sizer_12.Add(self.cancel_btn, 0, wx.ALIGN_BOTTOM|wx.ALIGN_CENTER_HORIZONTAL, 0)
        sizer_4.Add(sizer_12, 2, wx.ALIGN_CENTER_HORIZONTAL, 0)
        self.SetSizer(sizer_4)
        sizer_4.Fit(self)
        # end wxGlade

    def OnIntFileBrowse(self, event): # wxGlade: SpinwavePanel.<event_handler>
        #defaultDir=os.getcwd()
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
            interactionfile=paths[0].encode('ascii')
            self.int_file_txtCtrl.SetValue(interactionfile)
            #display in richtextcontrol
            self.editorWin.loadInteractions(interactionfile)

        dlg.Destroy()

    def OnSpinFileBrowse(self, event): # wxGlade: SpinwavePanel.<event_handler>
        #defaultDir=os.getcwd()
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
            spinfile=paths[0].encode('ascii')
            self.spin_file_txtCtrl.SetValue(spinfile)
            #display in richtextcontrol
            self.editorWin.loadSpins(spinfile)
        
        dlg.Destroy()

    def OnCancel(self, evt):
        """Closes window"""
        self.GetParent().Close()
        
    def OnOk(self, event):
        #Formerly, when this was modal, this would be done in calling function
        failed, data, kMin, kMax = self.Validate()
        
        if not failed:
            int_file = self.int_file_txtCtrl.GetValue()
            spin_file = self.spin_file_txtCtrl.GetValue()
            self.processManager.startAnalyticDispersion(int_file, spin_file)
            self.processManager.startNumericDispersion(int_file, spin_file, data, kMin*pi, kMax*pi, data['step'])
    
    def Validate(self):
        """Checks that all values are the right type. Any field that is not of the right
        type will be turned pink.
        
        Returns failed, data, kMin, kMax
        failed is True if validation fails and false otherwise."""
        kx = self.kx_txtCtrl.GetValue()
        ky = self.ky_txtCtrl.GetValue()
        kz = self.kz_txtCtrl.GetValue()
        kMin = self.kMin_txtCtrl.GetValue()
        kMax = self.kMax_txtCtrl.GetValue()
        
         
        bgColor = "pink"
        failed = False
        #Validate kx(must be a float)
        numKx = None
        try:
            numKx = float(kx)
            self.kx_txtCtrl.SetBackgroundColour("white")
        except:
            self.kx_txtCtrl.SetBackgroundColour(bgColor)
            failed = True
            
        #Validate ky(must be a float)
        numKy = None
        try:
            numKy = float(ky)
            self.ky_txtCtrl.SetBackgroundColour("white")
        except:
            self.ky_txtCtrl.SetBackgroundColour(bgColor)
            failed = True
            
        #Validate kz(must be a float)
        numKz = None
        try:
            numKz = float(kz)
            self.kz_txtCtrl.SetBackgroundColour("white")
        except:
            self.kz_txtCtrl.SetBackgroundColour(bgColor)
            failed = True
        
        #Validate kMin(must be a float)
        numKMin = None
        try:
            numKMin = float(kMin)
            self.kMin_txtCtrl.SetBackgroundColour("white")
        except:
            self.kMin_txtCtrl.SetBackgroundColour(bgColor)
            failed = True
            
        #Validate kMax(must be a float)
        numKMax = None
        try:
            numKMax = float(kMax)
            self.kMax_txtCtrl.SetBackgroundColour("white")
        except:
            self.kMax_txtCtrl.SetBackgroundColour(bgColor)
            failed = True
            
        
        data={}
        data['kx']=numKx
        data['ky']=numKy
        data['kz']=numKz
        data['step']=int(self.steps_spin_ctrl.GetValue())
        
        return failed, data, numKMin, numKMax
        
# end of class SpinwavePanel

#class FormDialog(sc.SizedPanel):
class FormDialog(wx.Panel):
    def __init__(self, parent, id, procManager):
        self.parent = parent
        #self.process_list = []
        self.processManager = procManager
        
        #valstyle=wx.WS_EX_VALIDATE_RECURSIVELY
        #sc.SizedPanel.__init__(self, parent, -1,
        #                style= wx.RESIZE_BORDER)#| wx.WS_EX_VALIDATE_RECURSIVELY)
        wx.Panel(self,parent, -1, style = wx.RESIZE_BORDER)
        self.SetExtraStyle(wx.WS_EX_VALIDATE_RECURSIVELY)
        pane = self
        pane.SetSizerType("vertical")
               
        print 'FormDialog called'
        #FilePane = sc.SizedPanel(pane, -1)
        FilePane = wx.Panel(pane, -1)
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
        #DirectionPane = sc.SizedPanel(pane, -1)
        DirectionPane = wx.Panel(pane, -1)
        DirectionPane.SetSizerType("vertical")
        DirectionPane.SetSizerProps(expand=True)
        wx.StaticText(DirectionPane, -1, "Scan Direction")
        
        #DirectionsubPane = sc.SizedPanel(pane, -1
        DirectionsubPane = wx.Panel(pane, -1)
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
        #kRangePane = sc.SizedPanel(pane, -1
        kRangePane = wx.Panel(pane, -1)
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
        
        #btnPane = sc.SizedPanel(pane, -1)
        btnPane = wx.Panel(pane, -1)
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
        #self.editorWin = RichTextFrame(self, -1, "Editor",
        #                    size=(620, 250),
        #                    style = wx.DEFAULT_FRAME_STYLE)
        #self.editorWin.Show(True)
        self.editorWin = showEditorWindow(self, "Spinwave File Editor")

        
        
    def EvtSpinCtrl(self,evt):
        print self.__dict__
        print 'event',evt.__dict__
        #spinCtrl=evt.control
        self.data['step']=self.spinctrl.GetValue()
        #sc.SizedDialog.GetW
        #self.G
        #self.kx=text
 
       


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


class MyApp(wx.App):
    def __init__(self, redirect=False, filename=None, useBestVisual=False, clearSigInt=True):
        wx.App.__init__(self,redirect,filename,clearSigInt)


    def OnInit(self):
        return True
    

if __name__=='__main__':
    from spinwaves.utilities.Processes import ProcessManager
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
