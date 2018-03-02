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
# the same thing works for atomtypes, chainids, resids, resnames
# segnames, xs, ys, zs, tempfs, occs
```
# Use `pdb_def` Object
