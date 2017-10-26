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
   std::unordered_multimap<PDBField,std::string> defstr;
   std::unordered_multimap<PDBField,float> defflt;
   std::unordered_multimap<PDBField,int> defint;

public:
   std::unordered_multimap<PDBField,std::string>& getDefstr() const;
};

inline
std::unordered_multimap<PDBField,std::string>& getDefstr() const
{   
   return defstr;
}

}

#endif
