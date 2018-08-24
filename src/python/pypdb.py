from pypdb_core import *
from functools import partial

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

class PDB(object):

    pdb_attrs = ["atomnames", "resnames", "segnames", "atomtypes", \
                 "chainids", "resids", "residues", "resnames", \
                 "xs", "ys", "zs", "occs", "tempfs"]
#    pdb_fields = ["atomname", "resname", "segname", "atomtype", \
#                 "chainid", "resid", "residue", "resname", \
#                 "x", "y", "z", "occ", "tempf"]

    
    def __init__(self, fname):
        self.core_data = pdb_obj(fname)
        self.natoms = self.core_data.natom
    def __setattr__(self, name, value):
        if name in self.pdb_attrs:
            object.__setattr__(self.core_data, name, value)
        else:
            object.__setattr__(self, name, value)
    def __getattr__(self, name):
        if name in self.pdb_attrs:
            return self.core_data.__getattribute__(name)
        else:
            object.__getattribute__(name)
    def get_atomname(self, arg):
        """
        get_atomname(index)
        get_atomname(selection_string)
        get_atomname(pdbdef)
        get_atomname(iterable_of_indexes)
        """
        if hasattr(arg, "__iter__"):
            return [self.core_data.atomnames[_] for _ in arg]
        else:
            try:
                return self.core_data.atomnames[arg]
            except IndexError as err: raise err
            except:
                return self.get_atomname(self.select_atoms(arg))
    def get_resid(self, arg):
        """
        get_resid(index)
        get_resid(selection_string)
        get_resid(pdbdef)
        get_resid(iterable_of_indexes)
        """
        if hasattr(arg, "__iter__"):
            return [self.core_data.resids[_] for _ in arg]
        else:
            try:
                return self.core_data.resids[arg]
            except IndexError as err: raise err
            except:
                return self.get_resid(self.select_atoms(arg))
    def get_chainid(self, arg):
        """
        get_chainid(index)
        get_chainid(selection_string)
        get_chainid(pdbdef)
        get_chainid(iterable_of_indexes)
        """
        if hasattr(arg, "__iter__"):
            return [self.core_data.chainids[_] for _ in arg]
        else:
            try:
                return self.core_data.chainids[arg]
            except IndexError as err: raise err
            except:
                return self.get_chainid(self.select_atoms(arg))
    def get_atomtype(self, arg):
        """
        get_atomtype(index)
        get_atomtype(selection_string)
        get_atomtype(pdbdef)
        get_atomtype(iterable_of_indexes)
        """
        if hasattr(arg, "__iter__"):
            return [self.core_data.atomtypes[_] for _ in arg]
        else:
            try:
                return self.core_data.atomtypes[arg]
            except IndexError as err: raise err
            except:
                return self.get_atomtype(self.select_atoms(arg))
    def get_chainid(self, arg):
        """
        get_chainid(index)
        get_chainid(selection_string)
        get_chainid(pdbdef)
        get_chainid(iterable_of_indexes)
        """
        if hasattr(arg, "__iter__"):
            return [self.core_data.chainids[_] for _ in arg]
        else:
            try:
                return self.core_data.chainids[arg]
            except IndexError as err: raise err
            except:
                return self.get_chainid(self.select_atoms(arg))
    def get_segname(self, arg):
        """
        get_segname(index)
        get_segname(selection_string)
        get_segname(pdbdef)
        get_segname(iterable_of_indexes)
        """
        if hasattr(arg, "__iter__"):
            return [self.core_data.segnames[_] for _ in arg]
        else:
            try:
                return self.core_data.segnames[arg]
            except IndexError as err: raise err
            except:
                return self.get_segname(self.select_atoms(arg))
    def get_resname(self, arg):
        """
        get_resname(index)
        get_resname(selection_string)
        get_resname(pdbdef)
        get_resname(iterable_of_indexes)
        """
        if hasattr(arg, "__iter__"):
            return [self.core_data.resnames[_] for _ in arg]
        else:
            try:
                return self.core_data.resnames[arg]
            except IndexError as err: raise err
            except:
                return self.get_resname(self.select_atoms(arg))
    def get_coordinate(self, arg):
        """
        get_coordinate(index)
        get_coordinate(selection_string)
        get_coordinate(pdbdef)
        get_coordinate(iterable_of_indexes)
        """
        if hasattr(arg, "__iter__"):
            return [self.core_data.get_coordinates(_) for _ in arg]
        else:
            try:
                return self.core_data.get_coordinates(arg)
            except IndexError as err: raise err
            except:
                return self.get_coordinate(self.select_atoms(arg))
    def select_atoms(self, atom_def):
        if type(atom_def) is str:
            return self.core_data.select_atoms(pdb_def(atom_def))
        else:
            return self.core_data.select_atoms(atom_def)
    def write2file(self, *args):
        """
        write2file(fname)
        write2file(fname, selection_string)
        write2file(fname, pdbdef)
        write2file(fname, vector_of_indexes)
        write2file(fname, iterable_of_indexes)
        """
        if len(args) == 1:
            self.core_data.write2file(args[0])
        elif len(args) == 2:
            try:
                self.core_data.write2file(args[0], args[1])
            except Exception as err:
                if not str(err).startswith("Python argument types in"): 
                    raise err
                if type(args[1]) is str:
                    self.core_data.write2file(args[0], pdb_def(args[1]))
                else:
                    self.core_data.write2file(\
                            args[0], hfs.make_vector_size_t(args[1]))
        else:
            raise TypeError("write2file accepts 1 or 2 args")
    def show(self, *args):
        """
        show()
        show(pdbdef)
        """
        if len(args) == 0:
            self.core_data.show()
        elif len(args) == 1:
            self.core_data.show(args[0])
        else:
            raise TypeError("show accepts 0 or 1 arg")
    def print_atoms(self, arg):
        """
        print_atoms()
        print_atoms(index)
        print_atoms(vector_of_indexes)
        print_atoms(pdbdef)
        print_atoms(selection_string)
        """
        if type(arg) == str:
            self.core_data.print_atom(pdb_def(arg))
        else:
            self.core_data.print_atom(arg)
    def geo_center(self, arg):
        """
        geo_center(pdbdef)
        geo_center(selection_string)
        geo_center(vector_of_indexes)
        geo_center(iterable_of_indexes)
        """
        try:
            # try to call geo_center(pdbdef)
            #         and geo_center(vector_of_indexes)
            return self.core_data.geo_center(arg)
        except Exception as err:
            if not str(err).startswith("Python argument types in"): 
                raise err
        if hasattr(arg, "__iter__"):
            # convert iterable to vector
            return self.core_data.geo_center(hfs.make_vector_size_t(arg))
        else:
            try:
                # try to convert string to pdb_def
                return self.core_data.geo_center(pdb_def(arg))
            except Exception as err:
                if not str(err).startswith("Python argument types in"):
                    raise
                else:
                    raise ValueError("illegal argument received by geo_center")

#for fs in PDB.pdb_fields:
#    setattr(PDB, "get_" + fs, \
#      lambda *args: \
#      fs_ = fs + "s"; \
#      args[0].__getattr__(fs_+"s") if len(args) == 1 \
#      else ( \
#        args[0].__getattr__(fs_+"s")[args[1]] if not hasattr(args[1],"__iter__")\
#        else \
#          [args[0].__getattr__(fs_+"s")[_] for _ in args[1]] \
#        ) \
#    )
