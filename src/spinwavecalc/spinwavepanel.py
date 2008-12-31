import wx
import wxaddons.sized_controls as sc
import  wx.lib.intctrl
import  wx.grid as  gridlib
import numpy as N
import sys,os
import spinwave_calc_file
import wx.richtext

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
                print 'child',child
            else:
                WalkTree(child)

    
class FormDialog(sc.SizedDialog):
    def __init__(self, parent, id):
        valstyle=wx.WS_EX_VALIDATE_RECURSIVELY
        sc.SizedDialog.__init__(self, None, -1, "Calculate Dispersion",
                        style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER)#| wx.WS_EX_VALIDATE_RECURSIVELY)
        self.SetExtraStyle(wx.WS_EX_VALIDATE_RECURSIVELY)
        pane = self.GetContentsPane()
        pane.SetSizerType("vertical")
        
        #Text editor
        interactionsEditorPanel = wx.Panel(pane, -1)
        sizer = wx.BoxSizer(wx.VERTICAL)
        interactionsEditorPanel.SetSizer(sizer)
        title = wx.StaticText(interactionsEditorPanel, -1, "Interaction File")
        self.interactionsRtc = wx.richtext.RichTextCtrl(interactionsEditorPanel, style=wx.VSCROLL | wx.TE_PROCESS_ENTER, size = (620,120));
        
        sizer.Add(title)
        sizer.Add(self.interactionsRtc, wx.EXPAND)
        
        spinsEditorPanel = wx.Panel(pane, -1)
        spinSizer = wx.BoxSizer(wx.VERTICAL)
        spinsEditorPanel.SetSizer(spinSizer)
        spinsTitle = wx.StaticText(spinsEditorPanel, -1, "Spins File")
        self.spinsRtc = wx.richtext.RichTextCtrl(spinsEditorPanel, style=wx.VSCROLL, size = (620,120));
        
        spinSizer.Add(spinsTitle)
        spinSizer.Add(self.spinsRtc, wx.EXPAND)
        
        
        
        FilePane = sc.SizedPanel(pane, -1)
        FilePane.SetSizerType("vertical")
        FilePane.SetSizerProps(expand=True)

#Never used?
#        InteractionFilePane = sc.SizedPanel(pane, -1)
#        InteractionFilePane.SetSizerType("vertical")
#        InteractionFilePane.SetSizerProps(expand=True)

        
        self.interactionfile=''
        self.interactionfilectrl=wx.StaticText(FilePane, -1, "Interaction File:%s"%(self.interactionfile,))
        b = wx.Button(FilePane, -1, "Browse", (50,50))
        self.Bind(wx.EVT_BUTTON, self.OnOpenInt, b)
        
        self.spinfile=''
        self.spinfilectrl=wx.StaticText(FilePane, -1, "Spin File:%s"%(self.spinfile,))
        b = wx.Button(FilePane, -1, "Browse", (50,50))
        self.Bind(wx.EVT_BUTTON, self.OnOpen, b)
        



        self.data={}
        self.data['kx']=str(1.0)
        self.data['ky']=str(0.0)
        self.data['kz']=str(0.0)
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
        # add dialog buttons
        
        #Getting rid of standard button function
        #self.SetButtonSizer(self.CreateStdDialogButtonSizer(wx.OK | wx.CANCEL))
        
        btnsizer = wx.StdDialogButtonSizer()
        
        self.btnOk = wx.Button(self, wx.ID_OK)
        self.btnOk.SetDefault()
        btnsizer.AddButton(self.btnOk)

        btn = wx.Button(self, wx.ID_CANCEL)
        btnsizer.AddButton(btn)
        btnsizer.Realize()
        
        self.SetButtonSizer(btnsizer)
        
        #intercept OK button click event
        self.Bind(wx.EVT_BUTTON, self.OnOk, self.btnOk)
        
        
        #self.TransferDataToWindow()
        # a little trick to make sure that you can't resize the dialog to
        # less screen space than the controls need
        self.Fit()
        self.SetMinSize(self.GetSize())

        
    def OnOk(self, event):
        #Formerly, when this was modal, this would be done in calling function
        self.Validate()
        print "OK"
        self.TransferDataFromWindow()
        print 'data',self.data
        print self.data['step']
        print self.interactionfile
        print self.spinfile
        spinwave_calc_file.driver(self.spinfile,self.interactionfile,self.data,self.data['step'])
        #button click event is not passed on, so the window does not close when OK is clicked

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
        defaultDir=wx.ConfigBase.Get().GetPath()
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
            self.spinsRtc.LoadFile(paths[0])
            
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
        defaultDir=wx.ConfigBase.Get().GetPath()
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
            self.interactionsRtc.LoadFile(paths[0])
            
        # Destroy the dialog. Don't do this until you are done with it!
        # BAD things can happen otherwise!
        dlg.Destroy()



if __name__=='__main__':
    app=MyApp()
    dlg=FormDialog(parent=None,id=-1)
    result=dlg.ShowModal()
    if result==wx.ID_OK:
        dlg.Validate()
        print "OK"
        dlg.TransferDataFromWindow()
        print dlg.data
    else:
        print "Cancel"

    dlg.Destroy()


