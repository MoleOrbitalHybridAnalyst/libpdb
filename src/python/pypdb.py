from pypdb_core import *

def make_vector_size_t(arr):
    """
    convert python iterable into c++ vector of size_t
    """
    v = std_vector_size_t()
    v.extend(arr)
    return v

def make_Vector(arr):
    """
    convert python iterable into c++ PDB_NS::Vector
    """
    v = Vector()
    v[0] = arr[0]
    v[1] = arr[1]
    v[2] = arr[2]
    return v

def __repr__(v):
   arr = [__ for __ in v]
   return arr.__repr__()

setattr(std_vector_size_t, "__repr__", __repr__)
setattr(std_vector_int, "__repr__", __repr__)
setattr(std_vector_float, "__repr__", __repr__)
setattr(std_vector_string, "__repr__", __repr__)
setattr(Vector, "__repr__", __repr__)
