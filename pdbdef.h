#ifndef __LIBPDB_PDBDEF_H
#define __LIBPDB_PDBDEF_H

//#include <vector>
#include <map>
#include <string>
#include "general.h"

namespace PDB_NS {

class PDBDef {
   //std::vector<std::pair<PDBField,std::string>> defstr;
   //std::vector<std::pair<PDBField,float>> defflt;
   //std::vector<std::pair<PDBField,int>> defint;
   //std::unordered_multimap<PDBField,std::string> _defstr;
   //std::unordered_multimap<PDBField,float> _defflt;
   //std::unordered_multimap<PDBField,int> _defint;
   //std::unordered_multimap<PDBField,char> _defchr;
   std::multimap<PDBField,std::string> _defstr;
   std::multimap<PDBField,float> _defflt;
   std::multimap<PDBField,int> _defint;
   std::multimap<PDBField,char> _defchr;

public:
   PDBDef() = default;
   PDBDef(const std::string& deffn);
   void pushBack(PDBField,const std::string&);
   void pushBack(PDBField,const float);
   void pushBack(PDBField,const int);
   void pushBack(PDBField,const char);
   //const std::unordered_multimap<PDBField,std::string>& getDefstr() const;
   //const std::unordered_multimap<PDBField,float>& getDefflt() const;
   //const std::unordered_multimap<PDBField,int>& getDefint() const;
   //const std::unordered_multimap<PDBField,char>& getDefchr() const;
   const std::multimap<PDBField,std::string>& getDefstr() const;
   const std::multimap<PDBField,float>& getDefflt() const;
   const std::multimap<PDBField,int>& getDefint() const;
   const std::multimap<PDBField,char>& getDefchr() const;
};

inline
const std::multimap<PDBField,std::string>& PDBDef::getDefstr() const
{   
   return _defstr;
}

inline
const std::multimap<PDBField,float>& PDBDef::getDefflt() const
{   
   return _defflt;
}

inline
const std::multimap<PDBField,int>& PDBDef::getDefint() const
{   
   return _defint;
}

inline
const std::multimap<PDBField,char>& PDBDef::getDefchr() const
{   
   return _defchr;
}

}

#endif
