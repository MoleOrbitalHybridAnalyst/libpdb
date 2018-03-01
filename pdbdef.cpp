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
      //throw runtime_error("cannot open "+deffn);
      pushBack(deffn);
   string line;
   while(getline(fs,line)) 
      pushBack(line);
}

void PDBDef::pushBack(const std::string& s)
{
   stringstream ss(s);
   string substring;
   ss >> substring;
   PDBField field = transString(substring);
   if(isString(field)) {
      while(ss >> substring) {
         _defstr.emplace(field, substring);
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

void PDBDef::pushBack(PDBField f, const std::string& s)
{
   if(isString(f)) _defstr.emplace(f,s);
   else throw invalid_argument("adding a string to non-string field");
}

void PDBDef::pushBack(PDBField f, const float x)
{
   if(isFloat(f)) _defflt.emplace(f,x);
   else throw invalid_argument("adding a float to non-float field");
}

void PDBDef::pushBack(PDBField f, const int n)
{
   if(isInt(f)) _defint.emplace(f,n);
   else throw invalid_argument("adding an int to non-int field");
}

void PDBDef::popBack(PDBField f)
{
   if(isString(f)) {
      auto range = _defstr.equal_range(f);
      if(range.first != range.second)
         _defstr.erase(--range.second);
   }
   else if(isInt(f)) {
      auto range = _defint.equal_range(f);
      if(range.first != range.second)
         _defint.erase(--range.second);
   }
   else if(isFloat(f)) {
      auto range = _defflt.equal_range(f);
      if(range.first != range.second)
         _defflt.erase(--range.second);
   }
}

void PDBDef::print() const
{
   {
      vector<PDBField> strfs;
      strfs.push_back(PDBField::atomname);
      strfs.push_back(PDBField::resname);
      strfs.push_back(PDBField::segname);
      strfs.push_back(PDBField::atomtype);
      strfs.push_back(PDBField::chainid);
      for(auto f : strfs)
      {
         auto range = _defstr.equal_range(f);
         if(range.first == range.second) continue;
         cout << transField(f) << ' ';
         for(auto it = range.first; it != range.second; ++it) {
            cout << it->second << ' ';
         }
         cout << endl;
      }
   }

   {
      auto f = PDBField::resid;
      auto range = _defint.equal_range(f);
      if(range.first == range.second) return;
      cout << transField(f) << ' ';
      for(auto it = range.first; it != range.second; ++it) {
         cout << it->second << ' ';
      }
      cout << endl;
   }
}

}
