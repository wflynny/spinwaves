#!/usr/bin/env python

"""This is a test of wxVTKRenderWindow"""

import wx
import wx.grid
from picker import Picker
from wxVTKRenderWindowInteractor import *
from vtkModel.SpaceGroups import *
from vtkModel.VTKDrawer import *
from vtkModel.MagneticCellClass import *
from vtkModel.CellClass import *
#import random









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










































#________________________________________________
#My Code






class atomPanel(wx.Panel):
    def __init__(self, parent, id):
        wx.Panel.__init__(self, parent, id)
        
        #Add Space Group
        spaceGroupLabel = wx.StaticText(self, -1, "Space Group:")
        self.spaceGroupSpinner = wx.SpinCtrl(self, -1, "")
        self.spaceGroupSpinner.SetRange(1,230)
        self.spaceGroupSpinner.SetValue(1)
        
        #Add Atom List
        self.atomList = atomListGrid(self, -1)
         
        #Add a button on upper right to generate new image
        self.genButton = wx.Button(self, -1, "Generate")
        
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
        
        #Magnetic Cell
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
        
        #Cutoff
        cutoffNaLabel = wx.StaticText(self, -1, "Na:")
        self.cutoffNaText = wx.TextCtrl(self, -1, size = (60, -1))
        cutoffNaSizer = wx.BoxSizer(wx.VERTICAL)
        cutoffNaSizer.Add(cutoffNaLabel, 0)
        cutoffNaSizer.Add(self.cutoffNaText, 0)
        
        cutoffNbLabel = wx.StaticText(self, -1, "Nb:")
        self.cutoffNbText = wx.TextCtrl(self, -1, size = (60, -1))
        cutoffNbSizer = wx.BoxSizer(wx.VERTICAL)
        cutoffNbSizer.Add(cutoffNbLabel, 0)
        cutoffNbSizer.Add(self.cutoffNbText, 0)
        
        cutoffNcLabel = wx.StaticText(self, -1, "Nc:")
        self.cutoffNcText = wx.TextCtrl(self, -1, size = (60, -1))
        cutoffNcSizer = wx.BoxSizer(wx.VERTICAL)
        cutoffNcSizer.Add(cutoffNcLabel, 0)
        cutoffNcSizer.Add(self.cutoffNcText, 0)
        
        
        #Organize a,b,c and alpha, gamma, beta , Na, Nb, Nc into a grid
        dimSizer = wx.GridSizer(cols = 3, hgap = 15, vgap = 5)
        dimSizer.Add(aSizer)
        dimSizer.Add(bSizer)
        dimSizer.Add(cSizer)
        dimSizer.Add(alphaSizer)
        dimSizer.Add(betaSizer)
        dimSizer.Add(gammaSizer)

        
        unitCellBox = wx.StaticBox(self, -1, "Unit Cell")
        unitCellSizer = wx.StaticBoxSizer(unitCellBox, wx.VERTICAL)
        unitCellSizer.Add(dimSizer)
        
        leftTopSizer = wx.GridBagSizer(2,2)
        leftTopSizer.Add(unitCellSizer, (0,0), (1,2))
        leftTopSizer.Add(wx.StaticText(self, -1, "Atoms:"), (1,0))
        self.atomSpinner = wx.SpinCtrl(self, -1, "")
        self.atomSpinner.SetRange(1,100)
        self.atomSpinner.SetValue(1)
        self.atomSpinner.Bind(wx.EVT_TEXT, self.OnGridResize, self.atomSpinner)
        leftTopSizer.Add(self.atomSpinner, (1,1))
        
        magCellBox = wx.StaticBox(self, -1, "Magnetic Cell")
        magCellSizer = wx.StaticBoxSizer(magCellBox, wx.HORIZONTAL)
        magCellSizer.Add(naSizer)
        magCellSizer.Add(nbSizer)
        magCellSizer.Add(ncSizer)
        
        cutoffBox = wx.StaticBox(self, -1, "Cutoff")
        cutoffSizer = wx.StaticBoxSizer(cutoffBox, wx.HORIZONTAL)
        cutoffSizer.Add(cutoffNaSizer)
        cutoffSizer.Add(cutoffNbSizer)
        cutoffSizer.Add(cutoffNcSizer)
        
        
        spaceGroupSizer = wx.BoxSizer(wx.HORIZONTAL)
        spaceGroupSizer.Add(spaceGroupLabel, 0)
        spaceGroupSizer.Add(self.spaceGroupSpinner, 0)
        
        MagCutoffSizer = wx.BoxSizer(wx.VERTICAL)
        MagCutoffSizer.Add(magCellSizer)
        MagCutoffSizer.Add(cutoffSizer)
        
        topSizer = wx.FlexGridSizer(cols = 2, hgap = 5, vgap = 5)
        topSizer.Add(spaceGroupSizer)
        topSizer.Add(self.genButton, 0, wx.ALIGN_RIGHT)
        topSizer.Add(leftTopSizer)
        topSizer.Add(MagCutoffSizer)
        

        
        
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(topSizer, 0)
        sizer.Add(self.atomList, 0, wx.EXPAND)
        
        self.SetSizer(sizer)
        
        #For now, the generate button will create a new magnetic cell
        #which will then be drawn by the vtkPanel
        self.Bind(wx.EVT_BUTTON, self.OnGenerate, self.genButton)

    def OnGenerate(self, event):
         self.validate()
         
         print self.atomList.GetTable().data
        
    def validate(self):
         a = self.aText.GetLabel() #str
         b = self.bText.GetLabel() #str
         c = self.cText.GetLabel() #str
         alpha = self.alphaText.GetLabel() #str
         beta = self.betaText.GetLabel() #str
         gamma = self.gammaText.GetLabel() #str
         spaceGroup = self.spaceGroupSpinner.GetValue() #int
         
    def OnGridResize(self, event):
        rows = self.atomSpinner.GetValue()
        self.atomList.SetNumberRows(rows)

#        self.atomList.GetTable().SetNumberRows(rows)
        event.Skip()
        
    
        
    
class atomTable(wx.grid.PyGridTableBase):
    def __init__(self):
        wx.grid.PyGridTableBase.__init__(self)
        self.colLabels = ['Name', 'Atomic Number','x', 'y','z']
        self.rowLabels=['Atom 1']
        
        self.data = [
                     ['','','','','']#Row 1
                     ]
    
    def GetNumberRows(self):
        return len(self.data)
    def AppendRows(self, num):
        self.AppendRow()
        print "here"
        return True
    def GetNumberCols(self):
        return len(self.colLabels)
    def GetColLabelValue(self, col):
        return self.colLabels[col]
    def GetRowLabelValue(self, row):
        return self.rowLabels[row]
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
            return ''
    def SetValue(self, row, col, value):
        try:
            self.data[row][col] = value
        except IndexError:
            # add a new row
            self.AppendRow()
            self.data[row][col]=value
        return
    
    def AppendRow(self):
            self.data.append([''] * self.GetNumberCols())

            # tell the grid we've added a row
            msg = wx.grid.GridTableMessage(self,            # The table
                    wx.grid.GRIDTABLE_NOTIFY_ROWS_APPENDED, # what we did to it
                    1                                       # how many
                    )
            self.GetView().ProcessTableMessage(msg)
            return True

    def DeleteRows(self,pos=0,numRows=1):
        if numRows>=0 and numRows<=self.GetNumberRows():
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


class atomListGrid(wx.grid.Grid):
    def __init__(self, parent, id):
        wx.grid.Grid.__init__(self, parent, id)
        self.table = atomTable()
        self.SetTable(self.table)
#        self.AutoSize()
    def SetNumberRows(self, num):
        diff = num - self.table.GetNumberRows()
        if diff > 0:
            for i in range(diff):
                self.table.AppendRow() #Add a blank row
        elif diff < 0:
            self.table.DeleteRows(self.table.GetNumberRows(), num)
        return diff

class vtkPanel(wx.Panel):
    def __init__(self, parent, id):
        wx.Panel.__init__(self, parent)
        self.initVTKWindow()

    
    def initVTKWindow(self):
        #Code from wxRenderWindowInterActor Sample
       
        self.window = wxVTKRenderWindowInteractor(self, -1)
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.window, 1, wx.EXPAND)
        self.SetSizer(sizer)
        self.Layout()
        # It would be more correct (API-wise) to call widget.Initialize() and
        # widget.Start() here, but Initialize() calls RenderWindow.Render().
        # That Render() call will get through before we can setup the
        # RenderWindow() to render via the wxWidgets-created context; this
        # causes flashing on some platforms and downright breaks things on
        # other platforms.  Instead, we call widget.Enable().  This means
        # that the RWI::Initialized ivar is not set, but in THIS SPECIFIC CASE,
        # that doesn't matter.
 #       self.window.Enable(1)
    
        self.window.AddObserver("ExitEvent", lambda o,e,f=self: f.Close())
        self.window.Render()
        
#        self.initializeVTKData()
#        self.draw()
        self.window.Bind(wx.EVT_KEY_DOWN, self.OnKeyEvent)
 
    
    
    def OnKeyEvent(self, event):
        event.Skip()
        
        #Delete if key pressed is Del key  Del = 127
        if event.GetKeyCode() == 127:
            self.OnDelete(event)
        
    
    
    def OnDelete(self,event):
        event.Skip()
        selectedObj = self.drawer.getObjFromActor(self.picker.getPicked())
        if isinstance(selectedObj, Bond):
            self.MagCell.deleteBond(selectedObj)
            #handle atoms next
        
        self.draw()
    
#    def initializeVTKData(self):
                #My Code
#        Space_Group = sg141
#        unitcell = Cell(Space_Group)
#        atomPos = [0, 0, 0]
    
        #Create the unit cell
#        unitcell.generateAtoms(atomPos, "atom1")

        
        #Create the Magnetic Cell
#        self.MagCell = MagneticCell(unitcell, 1,1,1, Space_Group)
        
#        AllAtoms = self.MagCell.getAllAtoms()
#        for i in range(0, len(AllAtoms)):
#            print i, AllAtoms[i]
#        self.MagCell.addBond(AllAtoms[0], AllAtoms[1])
        
    def setMagCell(self, MagCell):
        """Right now, other frames can change the magcell and update here"""
        self.MagCell = Magcell    
        
    def draw(self):
         
        # a renderer for the data
        ren1 = vtkRenderer()
        ren1.SetBackground(1,1,1)
        
        #Add the renderer to the window
        self.window.GetRenderWindow().AddRenderer(ren1)
            
        #Create vtkDrawer
        self.drawer = vtkDrawer(ren1)
        
        #Add my picker
        self.picker = Picker(self.drawer, self.window._Iren, ren1)

        #Draw the Magnetic Cell
        self.drawer.drawMagneticCell(self.MagCell)
        
        self.window.setUpRender()    

        self.drawer.addAxes()
        self.drawer.labelAtoms(self.MagCell)
        self.window.Render()
        
    
    def openCif(self, filename):
        self.MagCell = magneticCellFromCif(filename)
        self.draw()
    
    
    def getStatusText(self):
        return self.picker.getPicked()




class Frame(wx.Frame):
    def __init__(self, parent, id):
        wx.Frame.__init__(self, parent, id, 'Magnetic Cell', size= (900,600))

        self.vtkPanel = vtkPanel(self, -1)
        
        #Add Menus
        self.AddMenus()
        
        #Add Tool Bar
#        self.AddToolBar()
        
        #Add Status Bar                     
#        self.AddStatusBar()

   
   
#    def AddToolBar(self):
#        toolbar = self.CreateToolBar()
#        toolbar.AddSimpleTool()
        
#    def AddStatusBar(self):
#        self.statusBar = self.CreateStatusBar()
#        self.Bind(wx.EVT_LEFT_UP, self.OnClick)
#        
#    
#    def OnClick(self, event):
#        self.statusBar.SetStatusText("This")
#        #self.statusBar.SetStatusText(self.vtkPanel.getStatusText())
#        event.Skip()
        
    def AddMenus(self):
                #Add Menus
        menuBar = wx.MenuBar()
        
        #Add File Menu
        fileMenu = wx.Menu()
        newMenuItem = fileMenu.Append(wx.NewId(), "&New Magnetic Cell")
        openMenuItem = fileMenu.Append(wx.NewId(), "&Open")
        saveMenuItem = fileMenu.Append(wx.NewId(), "&Save")
        quitMenuItem = fileMenu.Append(wx.NewId(), "&Quit")
        menuBar.Append(fileMenu, "&File")
        
        #Add Model Menu
        modelMenu = wx.Menu()
        addCellMenuItem = modelMenu.Append(wx.NewId(), "Add Atom")
        addBondMenuItem = modelMenu.Append(wx.NewId(), "Add Bond")
        deleteMenuItem = modelMenu.Append(wx.NewId(), "Delete")
        menuBar.Append(modelMenu, "Model")
        self.SetMenuBar(menuBar)
        
        #Bind Events
        self.Bind(wx.EVT_MENU, self.OnCloseMe, quitMenuItem)
        self.Bind(wx.EVT_MENU, self.OnSave, saveMenuItem)
        self.Bind(wx.EVT_MENU, self.OnOpenFile, openMenuItem)
        self.Bind(wx.EVT_MENU, self.vtkPanel.OnDelete, deleteMenuItem)
#        self.Bind(wx.EVT_MENU, self.OnNew, newMenuItem)
    
 #   def OnNew(self, event):
        #Add drawing panel
#        self.vtkPanel.draw()
#        self.GetEventHandler().ProcessEvent(wx.SizeEvent())
    
    def OnCloseMe(self, event):
        self.Close(True)
    
    def OnSave(self, event):
        saveDialog = wx.FileDialog(self, "Save File", style = wx.SAVE)
        if saveDialog.ShowModal() == wx.ID_OK:
            print saveDialog.GetPath()
        saveDialog.Destroy()
        
    def OnOpenFile(self, event):
        saveDialog = wx.FileDialog(self, "Open File", style = wx.OPEN, wildcard = "*.cif")
        if saveDialog.ShowModal() == wx.ID_OK:
            self.vtkPanel.openCif(saveDialog.GetPath())
        saveDialog.Destroy()
        



class App(wx.App):
    def __init__(self, redirect = False, filename = None):
        wx.App.__init__(self, redirect, filename)
    
    def OnInit(self):
        self.frame = Frame(None, -1)
        self.frame.Show()
        self.SetTopWindow(self.frame)
        frame1 = wx.Frame(self.frame, -1)
        atomPanel(frame1, -1)
        frame1.Show()
        
        
        #Williams
#       frame2 = AtomFrame(self.frame, -1)
#        frame2.Show()
        
        return True
    


if __name__ == '__main__':
    app = App()
    app.MainLoop()
    
    
    


