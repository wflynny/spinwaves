import sys
import os
import sympy as sp
import numpy as np
import subprocess as sub
from subprocess import *


def generate_output(output, out = ".dvi"):
    """ Prints the given output in LaTeX """
    
    # Standard tex inputs required for compiling
    filename = os.path.join("c:","output")
    tex = ".tex"; pdf = ".pdf"; dvi = ".dvi"; ps = ".ps"
    begin = ["\documentclass[12pt]{article}\n",
                "\usepackage{palatino,url}\n",
                "\\begin{document}\n",
                "\section*{Cross-Section}\n\n"]
    end = ["\end{document}"]
    
    if not isinstance(output, str): 
        try:
            output = sp.latex(output)+"\n"
        except: e
    else: output = output+"\n\n"

    # Create file and write to it
    FILE = open(filename+tex, "w")
    FILE.writelines(begin)
    FILE.writelines(output)
    FILE.writelines(end)
    FILE.close()

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
