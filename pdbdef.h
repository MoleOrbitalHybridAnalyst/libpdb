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
   std::unordered_multimap<PDBField,char> defchr;

public:
   PDBDef() = default;
   PDBDef(const std::string& deffn);
   void pushBack(PDBField,const std::string&);
   void pushBack(PDBField,const float);
   void pushBack(PDBField,const int);
   void pushBack(PDBField,const char);
   const std::unordered_multimap<PDBField,std::string>& getDefstr() const;
   const std::unordered_multimap<PDBField,float>& getDefflt() const;
   const std::unordered_multimap<PDBField,int>& getDefint() const;
   const std::unordered_multimap<PDBField,char>& getDefchr() const;
};

inline
const std::unordered_multimap<PDBField,std::string>& PDBDef::getDefstr() const
{   
   return defstr;
}

inline
const std::unordered_multimap<PDBField,float>& PDBDef::getDefflt() const
{   
   return defflt;
}

inline
const std::unordered_multimap<PDBField,int>& PDBDef::getDefint() const
{   
   return defint;
}

inline
const std::unordered_multimap<PDBField,char>& PDBDef::getDefchr() const
{   
   return defchr;
}

}

#endif
