from pypdb_core import *

def make_vector_size_t(arr):
    v = std_vector_size_t()
    v.extend(arr)
    return v

def __repr__(v):
    arr = [__ for __ in v]
    return arr.__repr__()

setattr(std_vector_size_t, "__repr__", __repr__)
