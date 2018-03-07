# PyPDB
PyPDB is a Python library for reading and modifying [Protein Data Bank (PDB)](https://www.wwpdb.org/) files. The core is written in C++, which makes it also work as an independent C++ PDB library.
## Compile
```
mkdir build
cd build
```
### Compile as a pure C++ library
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
### Install as a Python module
Ensure the library binary can be searched as `pypdb_core.so` in `$PYTHONPATH` (softlink or mv whatever). Then the module can be imported by `import pypdb_core`. A possibly useful wrapper `src/python/pypdb.py` is also provided enabling one to import the module by `import pypdb`
