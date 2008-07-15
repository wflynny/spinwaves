import sys
import ctypes
from ctypes import c_int, c_float
import numpy as N

if sys.platform=='win32':
    print 'win32'
    mypolcorrect = N.ctypeslib.load_library('monteCarlo', 'C:\Dev-Cpp\workspace')
elif sys.platform=='mac':
    mypolcorrect = N.ctypeslib.load_library('libpolarization2.so', '.')
else:
    mypolcorrect = N.ctypeslib.load_library('libpolarization2.so', '.') #linux



neighbors1Class = c_int * 3
neighbors1 = neighbors1Class()
neighbors1[0] = 1
neighbors1[1] = 2
neighbors1[2] = 3

neighbors2Class = c_int * 2
neighbors2 = neighbors2Class()
neighbors2[0] = 4
neighbors2[1] = 5

sClass = c_float * 3
s1 = sClass()
s1[0] = 1.1
s1[1] = 2.2
s1[2] = 3.3
s2 = sClass()
s1[0] = 4.4
s1[1] = 5.5
s1[2] = 6.6

mat1Class = c_int * 3
mat1 = mat1Class()
mat1[0] = 1
mat1[1] = 2
mat1[2] = 3

mat2Class = c_int * 2
mat2 = mat2Class()
mat2[0] = 4
mat2[1] = 5


listPointer = mypolcorrect.new_atom_list(c_int(2))
mypolcorrect.set_atom(listPointer, 0, mat1, neighbors1, s1)
mypolcorrect.set_atom(listPointer, 1, mat2, neighbors2, s2)
mypolcorrect.atomTest(listPointer, c_int(2))