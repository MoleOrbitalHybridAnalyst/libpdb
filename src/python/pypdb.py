from pypdb_core import *

# enable syntax completion
# https://stackoverflow.com/questions/246725/how-do-i-add-tab-completion-to-the-python-shell
try:
    import readline
except ImportError:
    print("Module readline not available.")
else:
    import rlcompleter
    readline.parse_and_bind("tab: complete")

class HelpingFunctions:
    def __init__(self):
        self.repr_iterable = lambda v : [__ for __ in v].__repr__()

    def make_vector_size_t(self, arr):
        """
        convert python iterable into c++ vector of size_t
        """
        v = std_vector_size_t()
        v.extend(arr)
        return v
    
    def make_Vector(self, arr):
        """
        convert python iterable into c++ PDB_NS::Vector
        """
        v = Vector()
        v[0] = arr[0]
        v[1] = arr[1]
        v[2] = arr[2]
        return v
    

hfs = HelpingFunctions()
setattr(std_vector_size_t, "__repr__", hfs.repr_iterable)
setattr(std_vector_int, "__repr__", hfs.repr_iterable)
setattr(std_vector_float, "__repr__", hfs.repr_iterable)
setattr(std_vector_string, "__repr__", hfs.repr_iterable)
setattr(Vector, "__repr__", hfs.repr_iterable)

class PDB:
    def __init__(self, fname):
        self.pdb = pdb_obj(fname)
        self.natoms = self.pdb.natom
    def select_atoms(self, atom_def):
        if type(atom_def) is str:
            return self.pdb.select_atoms(pdb_def(atom_def))
        else:
            return self.pdb.select_atoms(atom_def)
    def write2file(self, *args):
        """
        write2file(fname)
        write2file(fname, selection_string)
        write2file(fname, pdbdef)
        write2file(fname, vector_of_indexes)
        write2file(fname, iterable_of_indexes)
        """
        if len(args) == 1:
            self.pdb.write2file(args[0])
        elif len(args) == 2:
            try:
                self.pdb.write2file(args[0], args[1])
            except:
                if type(args[1]) is str:
                    self.pdb.write2file(args[0], pdb_def(args[1]))
                else:
                    self.pdb.write2file(\
                            args[0], hfs.make_vector_size_t(args[1]))
        else:
            raise Exception("write2file accepts 1 or 2 args")
    def show(self, *args):
        """
        show()
        show(pdbdef)
        """
        if len(args) == 0:
            self.pdb.show()
        elif len(args) == 1:
            self.pdb.show(args[0])
        else:
            raise Exception("show accepts 0 or 1 arg")
    def print_atoms(self, arg):
        """
        print_atoms()
        print_atoms(index)
        print_atoms(vector_of_indexes)
        print_atoms(pdbdef)
        print_atoms(selection_string)
        """
        if type(arg) == str:
            self.pdb.print_atom(pdb_def(arg))
        else:
            self.pdb.print_atom(arg)
    def geo_center(self, arg):
        """
        geo_center(pdbdef)
        geo_center(selection_string)
        geo_center(vector_of_indexes)
        geo_center(iterable_of_indexes)
        """
        pass
