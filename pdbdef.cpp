#include "pdbdef.h"
#include <fstream>
#include <sstream>
#include <iostream>

using namespace std;

namespace PDB_NS {

PDBDef::PDBDef(const std::string& deffn)
{
   ifstream fs(deffn);
   if(!fs.is_open()) 
      throw runtime_error("cannot open "+deffn);
   string line;
   while(getline(fs,line)) {
      stringstream ss(line);
      string substring;
      ss >> substring;
      PDBField field = transString(substring);
      if(isString(field)) {
         while(ss >> substring) {
            _defstr.emplace(field, substring);
         }
      } else if(isChar(field)) {
         while(ss >> substring) {
            _defchr.emplace(field, substring[0]);
         }
      } else if(isFloat(field)) {
         float value;
         while(ss >> value) {
            _defflt.emplace(field, value);
         }
      } else if(isInt(field)) {
         int value;
         while(ss >> value) {
            _defint.emplace(field, value);
         }
      }
   }
}

void PDBDef::pushBack(PDBField f, const std::string& s)
{
   if(isString(f)) _defstr.emplace(f,s);
}

void PDBDef::pushBack(PDBField f, const float x)
{
   if(isFloat(f)) _defflt.emplace(f,x);
}

void PDBDef::pushBack(PDBField f, const int n)
{
   if(isInt(f)) _defint.emplace(f,n);
}

void PDBDef::pushBack(PDBField f, const char c)
{
   if(isChar(f)) _defchr.emplace(f,c);
}

void PDBDef::print() const
{
   {
      vector<PDBField> strfs;
      strfs.push_back(PDBField::atomname);
      strfs.push_back(PDBField::resname);
      strfs.push_back(PDBField::segname);
      strfs.push_back(PDBField::atomtype);
      for(auto f : strfs)
      {
         cout << transField(f) << endl;
         auto range = _defstr.equal_range(f);
         for(auto it = range.first; it != range.second; ++it) {
            cout << it->second << ' ';
         }
         cout << endl;
      }
   }

   {
      auto f = PDBField::chainid;
      cout << transField(f) << endl;
      auto range = _defchr.equal_range(f);
      for(auto it = range.first; it != range.second; ++it) {
         cout << it->second << ' ';
      }
      cout << endl;
   }

   {
      auto f = PDBField::resid;
      cout << transField(f) << endl;
      auto range = _defint.equal_range(f);
      for(auto it = range.first; it != range.second; ++it) {
         cout << it->second << ' ';
      }
      cout << endl;
   }
}

}
