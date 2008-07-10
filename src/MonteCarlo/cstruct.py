import ctypes, numpy, os

# Factory for ctypes/numpy struct classes
def cstruct(*fields):
  class SpecificStruct(object):
    """
    C structure usable from python.

    Use a._pointer to pass the data as a pointer to structure in
    ctypes.
    """
    _dtype = numpy.dtype(list(fields))
    def __init__(self,**kw):
        self.__dict__['_array'] = numpy.zeros((),dtype=self._dtype)
        for k,v in kw.iteritems(): self._array[k] = v

    def _getdata(self): return self._array.ctypes.data
    _pointer = property(_getdata,doc='ctypes data pointer to struct')
    def __getattr__(self,field):
        return self._array[field]
    def __setattr__(self,field,value):
        self._array[field] = value
  return SpecificStruct

if __name__ == "__main__":
    # Create a structure for:
    #    struct {
    #       int n, k;
    #       double A[4][4],B[4][4];
    #    }
    Embed = cstruct(('n',numpy.intc),('k',numpy.intc),
                    ('A',float,(4,4)),('B',float,(4,4)))

    # Initialize the structure and invoke call() in embedarray.dll
    a = Embed(k=6)
    a.A = numpy.array([[1,2,3,4],[5,6,7,8],[9,10,11,12],[13,14,15,16]])
    a.A[0,0] = 12
    lib = numpy.ctypeslib.load_library('embedarray',os.path.dirname(__file__))
    lib.call(a._pointer)
    print a.A,a.k,a.n

