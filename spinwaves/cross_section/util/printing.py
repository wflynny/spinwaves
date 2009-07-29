import sys
import os
import sympy as sp
import numpy as np
import subprocess as sub
from subprocess import *
import sympy.galgebra.latex_ex as tex
import sympy.galgebra.GA as GA
import wx
from multiprocessing import Process

class LaTeXDisplayFrame(wx.Frame):
    def __init__(self, parent, ID, input_text, title):
        self.parent = parent
        self.PID = ID
        self.input = input_text

        wx.Frame.__init__(self, parent, -1, title, size=(300,250))
        panel = wx.Panel(self,-1)
        multiLabel = wx.StaticText(panel, -1, "Analytic")
        multiText = wx.TextCtrl(panel, -1, self.input, size=(225,200), style=wx.TE_MULTILINE)
        multiText.SetInsertionPoint(0)
        
        sizer = wx.FlexGridSizer(cols=2, hgap=6, vgap=6)
        sizer.AddMany([multiLabel, multiText])
        panel.SetSizer(sizer)

def eig_process(mat):
    process_info('function eig_process')
    if isinstance(mat, np.ndarray):
        print mat.eig
        return mat.eig
    elif isinstance(mat, sp.matrices.Matrix):
        #print mat.eigenvals().keys()
        return mat.eigenvals().keys()

def process_info(title):
    print title
    print 'module name:', __name__
    print 'process id:', os.getpid()        
        
def create_latex(conn, input, name = None):
    pieces = []
    print type(input)
    print input
    if not isinstance(input, str):
        if isinstance(input, list):
            try:
                print 'list'
                for i in range(len(input)):
                    pieces.append(sp.latex(input[i]))
            except: e
        else:
            try:
                output = sp.latex(input)+"\n"
            except: 
                e
                print e
    else: output = output+"\n\n"
    if pieces != []:
        output = '\n\n'.join(pieces)
    else:
        if output.find(',') > 0:
            output = '\n'.join(output.split(','))
    output = output.replace('operatorname','mbox')
    if name != None:
        output = '\section{'+ name + '}' +'\n\n'+ output
    conn.send(output)
    conn.close()
    
    print "Output OUT"


def generate_output_3(output):
    
    m = Mathtex(output, u)
    m.save('testnew.png', 'png')


GA.set_main(sys.modules[__name__])

def generate_output_2(output):
    tex.Format()
    tex.sym_format(1)
    x = tex.print_LaTeX(output)
    print x
    #tex.xdvi(debug=True)

def generate_output(output, out = ".dvi"):
    """ Prints the given output in LaTeX """
    
    # Standard tex inputs required for compiling
    filename = os.path.join("c:","output")
    tex = ".tex"; pdf = ".pdf"; dvi = ".dvi"; ps = ".ps"
    begin = ["\documentclass[12pt]{article}\n",
                "\usepackage{amsmath,url}\n",
                "\\begin{document}\n",
                "\section{Cross-Section}\n\n"]
    end = ["\end{document}"]
    
    if not isinstance(output, str): 
        try:
            output = sp.latex(output)+"\n"
        except: e
    else: output = output+"\n\n"
    output = '\n'.join(output.split(','))
    #split = output.split('+')
    #reformatted = '$ \n\n$+'.join(split)
    #reformatted = reformatted.replace('$','')
    #print reformatted
    #output = reformatted

    # Create file and write to it
    FILE = open(filename+tex, "w")
    FILE.writelines(begin)
    FILE.writelines(output)
    FILE.writelines(end)
    FILE.close()

    if 1:
        # Create commands
        compile = ["latex",filename+tex]
        disdvi = ["yap", filename+dvi]
    
        # Process commands
        a = sub.Popen(compile,stdin=PIPE,stdout=PIPE,stderr=STDOUT)
        a.communicate()
        a.wait()
    
        # BROKEN
        if out == "pdf":
            tops = ["dvips", filename+dvi]
            topdf = ["ps2pdf", filename+ps]
            dispdf = ["C:/Program Files/Adobe/Reader 9.0/Reader/AcroRd32", filename+pdf]
            c = sub.check_call(tops)
    #        c = sub.Popen(tops,stdin=PIPE,stdout=PIPE,stderr=STDOUT)
    #        c.communicate
    #        c.wait()
            d = sub.Popen(topdf,stdin=PIPE,stdout=PIPE,stderr=STDOUT)
            d.communicate
            d.wait()
            e = sub.Popen(dispdf,stdin=PIPE,stdout=PIPE,stderr=STDOUT)
            e.communicate
        else:
            b = sub.Popen(disdvi,stdin=PIPE,stdout=PIPE,stderr=STDOUT)
            b.communicate()

if __name__ == "__main__":
    if 1:
        x,y = sp.symbols('xy')
        create_latex(x+2*y)