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

