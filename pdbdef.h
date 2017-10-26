#ifndef __LIBPDB_PDBDEF_H
#define __LIBPDB_PDBDEF_H

//#include <vector>
#include <unordered_map>
#include <string>
#include "general.h"

namespace PDB_NS {

class PDBDef {
   //std::vector<std::pair<PDBField,std::string>> defstr;
   //std::vector<std::pair<PDBField,float>> defflt;
   //std::vector<std::pair<PDBField,int>> defint;
   //std::vector<std::pair<PDBField,char>> defc;
   std::unordered_map<PDBField,std::string> defstr;
   std::unordered_map<PDBField,float> defflt;
   std::unordered_map<PDBField,int> defint;
   std::unordered_map<PDBField,char> defc;
}

#endif
