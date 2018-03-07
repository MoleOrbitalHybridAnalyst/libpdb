# PyPDB
A c++ library for reading and modifying [Protein Data Bank (PDB)](https://www.wwpdb.org/) files with a Python wrapper
## Compile
```
mkdir build
cd build
```
### Compile as a pure c++ library
```
camke ..
make
```
### Compile with the Python wrapper
[Boost.Python](http://www.boost.org/doc/libs/1_66_0/libs/python/doc/html/index.html) is needed
```
cmake -DMAKE_PYTHON=ON ..
make
```
## Install
Ensure the library binary can be searched as `pypdb_core.so` in `$PYTHONPATH` (softlink or mv whatever). Then the module can be imported by `import pypdb_core`. A possibly useful wrapper `src/python/pypdb.py` is also provided.
