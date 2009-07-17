import sys
import os
import sympy as sp
import numpy as np
import subprocess as sub
from subprocess import *
import sympy.galgebra.latex_ex as tex
import sympy.galgebra.GA as GA

def generate_output_3(output):
    
    m = Mathtex(output, u)
    m.save('testnew.png', 'png')


GA.set_main(sys.modules[__name__])

def generate_output_2(output):
    tex.Format()
    tex.sym_format(1)
    tex.print_LaTeX(output)
    tex.xdvi(debug=True)

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
        generate_output_2(x+2*y)