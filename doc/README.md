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
## `pdb_def` provides a way of selecting atoms
Create a `pdb_def` by giving a constraint or read from a file
```
# defo selects atoms whose atomtype is O and chainid is W
defo = pp.pdb_def("atomtype O and chainid W")
# defhyd selects atoms whose resname is H3O
defhyd = pp.pdb_def("resname H3O")
```
Get a vector object containing atom indexes of the hydronium
```
idxhyd = pdb.select_atoms(defhyd)
```
Then the vector can be converted into a python list
```
[_ for _ in idxhyd]
```
`print_atom` method accepts a `pdb_def` to print the selected atoms
```
pdb.print_atom(defhyd)
# print the atoms whose type is O or H
pdb.print_atom(pp.pdb_def("atomtype O H"))
```
If one wants to do something to all atoms in the pdb, create a `pdb_def` that selects all the atoms
```
pdb.print_atom(pp.pdb_def("all"))
```
Output the selected atoms to a pdb file
```
pdb.write2file("hyd.pdb", defhyd)
# or equivalently
pdb.write2file(defhyd, "hyd.pdb")
```
## Manipulate a `pdb_def`
Other constraints can be added to an existing `pdb_obj`
```
defhyd.push_back("atomtype O")
defhyd.push_back("resname TIP3")
```
Now defhyd selects atoms whose resname is H3O or TIP3 and atomtype is O
this can be checked by print out what is in the pdb_def
````
defhyd.show()
````
A constraint can be cleared by using `clear`
````
# let defhyd also select resid
defhyd.push_back("resid 1 2 3 4")
# recover the original defhyd
defhyd.clear("resid")
````
# Some Useful Methods
## `assemble_water`
Some people like to make a pdb of proton in water in the order of "OOOO...HHHHHHHH...", making people hard to see which H is bonded to which O. This kind of pdb can be converted into an order of "OHHOHHOHH...OHHH" by using `assemble_water`
````
# define what is water oxygen
defo = pp.pdb_def("chainid W and atomtype O")
# define what is water(including hydronium) hydrogen
defh = pp.pdb_def("chainid W and atomtype H")
pdb.assemble_water(defo, defh)
````
## `reorder_water`
Reactive MD or MC can result in messy pdbs where an O can be followed by two non-bonded H's. `reorder_water` can adjust the position of H in pdb to ensure the best bonding topology
````
pdb.reorder_water(defo, defh, defhyd)
````
Now, two H's just below one O are bonded to it. Consecutive OHHH is a hydronium
## `get_solvation_shells`
Get `n` solvation shells of waters of a hydronium based on a distance criteria
````
sol = pdb.get_solvation_shells(3,2.8,defo,defhyd)
````
`sol` is a vector containing the atom indexes of 3 solvation shells of water of the hydronium. Water oxygens which are within 2.8 A of the bonded H of the central species are regarded in the first solvation shell of it.
