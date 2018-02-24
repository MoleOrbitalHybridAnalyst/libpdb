cmake_minimum_required(VERSION 3.8)
project(libpdb)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES
        general.h
        pdb.cpp
        pdb.h
        pdbdef.cpp
        pdbdef.h
        utili.h
        utili.cpp)

#add_executable(libpdb ${SOURCE_FILES})
add_library(pdb SHARED ${SOURCE_FILES})
