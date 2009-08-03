"""
This file contains 4 methods to create LaTeX output. 
1.) The first method creates the latex popup wxpython window
    - uses the create_latex, eig_process, process_info methods and LaTeXDisplayFrame class
2.) The second method creates a .tex file, then compiles it using latex.exe. Finally it is displayed using Yap
    - uses the generate_output method
3.) The third method creates LaTeX output using sympy's GA module. It can generate an .tex file and pdf in the method
    - uses the generate_output_2 method
4.) The fourth method creates LaTeX output using Freddie's MathTex package.
    - uses the generate_output_3 method
    - not functional 

"""

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

#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------

def create_latex(conn, input, name = None):
    """ This method creates LaTeXified output given:
             - a sympy expression
             - a list of sympy expressions
             - a string
        It is currently compatible with the multiprocessing module, hence the 'conn' argument it takes.
        Arguments:
             - conn - connection from a multiprocessing pipe. probably needs to be reworked to accept a queue instead of a pipe
             - input - a sympy expression, a list of sympy expressions or a string
             - name - name for the LaTeX \section{} command
    """
    
    pieces = []
    # Crappy method to find out the type of the input, and then LaTeXify it
    if not isinstance(input, str):
        
        # Input is a list. Break it up and try to LaTeXify each piece
        if isinstance(input, list):
            try:
                print 'list'
                for i in range(len(input)):
                    pieces.append(sp.latex(input[i]))
            except: e
        # Input is probably just a sympy expression
        else:
            try:
                output = sp.latex(input)+"\n"
            except: 
                e
                print e
    
    # Input is a string
    else: output = output+"\n\n"

    
    # If the input was a list, join all the pieces into one string with 2 spaces between them. 
    if pieces != []:
        output = '\n\n'.join(pieces)
    # If the LaTeXifed input has any commas in it, split the expression at those commas and put some blank lines in between
    else:
        if output.find(',') > 0:
            output = '\n'.join(output.split(','))
    
    # Replace operatorname with mbox incase you don't want to use amsmath LaTeX package
    output = output.replace('operatorname','mbox')
    
    # Create a section if a name is supplied
    if name != None:
        output = '\section{'+ name + '}' +'\n\n'+ output
        
    # Send the output through the pipe
    conn.send(output)
    conn.close()
    print "Output OUT"

def eig_process(mat):
    """ This method takes a matrix and returns the eigenvalues from the matrix""" 
    process_info('function eig_process')
    if isinstance(mat, np.ndarray):
        print mat.eig
        return mat.eig
    elif isinstance(mat, sp.matrices.Matrix):
        # Currently, sympy's quartic/cubic polynomial solvers suck so this is currently out of commission. 
        #print mat.eigenvals().keys()
        #return mat.eigenvals().keys()
        return 1

def process_info(title):
    """ This method just prints process names and IDs """
    print title
    print 'module name:', __name__
    print 'process id:', os.getpid()         
    
class LaTeXDisplayFrame(wx.Frame):
    """ Basic text box frame class """
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

#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------        
#---------------------------------------------------------------------------------
        
def generate_output(output, out = ".dvi"):
    """ Generates a .tex file by LaTeXifying the input.
        Then, using latex.exe, it compiles the .tex file
        and displays it using Yap
    """
    
    # Standard tex inputs required for compiling .tex file
    filename = os.path.join("c:","output")
    tex = ".tex"; pdf = ".pdf"; dvi = ".dvi"; ps = ".ps"
    begin = ["\documentclass[12pt]{article}\n",
                "\usepackage{amsmath,url}\n",
                "\\begin{document}\n",
                "\section{Cross-Section}\n\n"]
    end = ["\end{document}"]
    
    pieces = []
    # Crappy method to find out the type of the input, and then LaTeXify it
    if not isinstance(input, str):
        
        # Input is a list. Break it up and try to LaTeXify each piece
        if isinstance(input, list):
            try:
                print 'list'
                for i in range(len(input)):
                    pieces.append(sp.latex(input[i]))
            except: e
        # Input is probably just a sympy expression
        else:
            try:
                output = sp.latex(input)+"\n"
            except: 
                e
                print e
    
    # Input is a string
    else: output = output+"\n\n"

    # If the input was a list, join all the pieces into one string with 2 spaces between them. 
    if pieces != []:
        output = '\n\n'.join(pieces)
    # If the LaTeXifed input has any commas in it, split the expression at those commas and put some blank lines in between
    else:
        if output.find(',') > 0:
            output = '\n'.join(output.split(','))

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


#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------        
#---------------------------------------------------------------------------------
GA.set_main(sys.modules[__name__])

def generate_output_2(output):
    """ See sympy's GA module online """ 
    tex.Format()
    tex.sym_format(1)
    x = tex.print_LaTeX(output)
    print x
    #tex.xdvi(debug=True)
    
#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------        
#---------------------------------------------------------------------------------
    
def generate_output_3(output):
    """ Uses the MathTex summer of code project """ 
    m = Mathtex(output, u)
    m.save('testnew.png', 'png')

#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------        
#---------------------------------------------------------------------------------    
    
if __name__ == "__main__":
    if 1:
        x,y = sp.symbols('xy')
        generate_output(x+2*y)