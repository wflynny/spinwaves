import sympy as sp
import numpy as np
import subprocess as sub
from subprocess import *

p=sub.Popen(["C;/Program Files/MiKTeX 2.7/miktex/bin/latex.exe",],stdin=PIPE,stdout=PIPE,stderr=STDOUT)
p.stdin.write("new\r\n")
p.stdin.write("exit\r\n")