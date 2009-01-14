#!/usr/bin/env python

"""A Hello program in wxPython"""


import wx

class Frame(wx.Frame):
    """Frame class that displays image"""
    
    def __init__(self, image, parent = None, id = -1,
                 pos = wx.DefaultPosition,
                 title = 'Hello!'):
        """Create a frame instance and display an image"""
        temp = image.ConvertToBitmap()
        size = temp.GetWidth(), temp.GetHeight()
        wx.Frame.__init__(self, parent, id, title, pos, size)
        self.bmp = wx.StaticBitmap(parent = self, bitmap = temp)

class App(wx.App):
    """Application Class"""
    
    def OnInit(self):
        image = wx.Image('C:\\menu.jpg', wx.BITMAP_TYPE_JPEG)
        self.frame = Frame(image)
        sself.frame.Show()
        self.SetTopWindow(self.frame)
        return True

def main():
    app = App()
    app.MainLoop()
    
if __name__ == '__main__':
    main()