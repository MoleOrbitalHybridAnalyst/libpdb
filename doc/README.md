# Quick Start
Read a PDB file by creating a `pdb_obj`
```
import pypdb as pp
pdb = pp.pdb_obj("input.pdb")
```
Save a `pdb_obj` to file
```
pdb.write2file("output.pdb")
```
Show atom information with an index
```
# print the information of atom 0
pdb.print_atom(0)
```
Acess the PDB data by the name of fields
```
# access the atom names
pdb.atomnames
# get the atom name of atom 0
pdb.atomnames[0]
# modify the atom name of atom 0
pdb.atomnames[0] = "XXX"
# loop in all the atom names
for _ in pdb.atomnames:
  print(_)
# the same thing works for atomtypes, chainids, resids, resnames
# segnames, xs, ys, zs, tempfs, occs
```
# Use `pdb_def` Object
`pdb_def` provides a way of selecting atoms
```
# create a pdb_def by giving a constraint
defhyd = pp.pdb_def("resname H3O")
# get a vector of atom indexes of the hydronium
idxhyd = pdb.select_atoms(defhyd)
# convert the vector to python list
[_ for _ in idxhyd]
# print the information of selected atoms
pdb.print_atom(defhyd)
# print the information of all atoms in the pdb
pdb.print_atom(pp.defhyd("all"))
# output the select atoms to a pdb file
pdb.write2file("hyd.pdb", defhyd)
# or equivalently
pdb.write2file(defhyd, "hyd.pdb")
```
Manipulate a `pdb_def`
```
# add other constraints
defhyd.push_back("atomtype O")
defhyd.push_back("resname TIP3")
# now defhyd selects atoms whose resname is H3O or TIP3 and atomtype is O
# this can be checked by print out what is in the pdb_def
defhyd.show()
```
