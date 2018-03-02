#ifndef __LIBPDB_PDBDEF_H
#define __LIBPDB_PDBDEF_H

//#include <vector>
#include <map>
#include <string>
#include "general.h"

namespace PDB_NS {

class PDBDef {
   std::multimap<PDBField,std::string> _defstr;
   std::multimap<PDBField,float> _defflt;
   std::multimap<PDBField,int> _defint;
   //std::multimap<PDBField,char> _defchr;
   bool _all;

public:
   PDBDef() = default;
   PDBDef(const std::string& deffn);
   void pushBack(PDBField,const std::string&);
   void pushBack(PDBField,const float);
   void pushBack(PDBField,const int);
/// push back by giving one line
   void pushBack(const std::string&);
   void popBack(PDBField);
   void clear(PDBField f) {_defstr.erase(f);_defint.erase(f);_defflt.erase(f);}
   //void pushBack(PDBField,const char);
   const std::multimap<PDBField,std::string>& getDefstr() const;
   const std::multimap<PDBField,float>& getDefflt() const;
   const std::multimap<PDBField,int>& getDefint() const;
   //const std::multimap<PDBField,char>& getDefchr() const;
   bool empty() const;
   void print() const;
   bool all() const {return _all; }
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

//inline
//const std::multimap<PDBField,char>& PDBDef::getDefchr() const
//{   
//   return _defchr;
//}

inline
bool PDBDef::empty() const
{
   return 
      (!_all) && _defstr.empty() && _defflt.empty() && _defint.empty();
      //_defstr.empty() && _defflt.empty() && _defint.empty() && _defchr.empty();
}

}

#endif
