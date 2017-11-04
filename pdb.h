#ifndef __LIBPDB_H
#define __LIBPDB_H
#include <string>
#include <vector>
#include <array>
#include <cmath>
#include <unordered_map>
#include "pdbdef.h"
#include "general.h"

namespace PDB_NS {

class PDB {
//pdb file format:
//ATOM serial atomname resname chainid resid x y z occ tempf segname atomtytpe
   std::vector<std::string> atomnames, resnames, segnames, atomtypes;
   std::vector<size_t> linenumbers;
   std::vector<char> chainids;
   std::vector<int> resids;
   std::vector<float> xs, ys ,zs, occs, tempfs;
   std::vector<std::pair<size_t,std::string>> nonatomlines;
   std::vector<std::vector<bool>> defineds;
   float boxlens[3]; size_t nAtoms;
   void centerAlignedPrint4(FILE *fp, const std::string& s) const;
/// calculate pbc distance pbc(x2-x1)
   float pbcDiff(double x1, double x2, int dim) const;
/// call moveTo and handle with two index lists //TODO using ...
   bool moveToWithIndexes(const size_t i1, const size_t i2,
                        std::vector<size_t>& list1, std::vector<size_t>& list2);
   //std::string transField(const PDBField& pdbfield) const;

   //void eraseSpace(std::string& str);

public:
   PDB() = default;
   explicit PDB(const std::string& fname);
   void write2file(const std::string& fname) const;
/// get all the undefined pairs(atom_index, str(pdbfiled))
   //std::vector<std::pair<size_t,std::string>> checkUndefined() const;
   std::vector<std::pair<size_t,PDBField>> checkUndefined() const;
/// calculate pbc distance pbc(x2-x1)
   std::array<float,3> pbcDistance(const std::array<float,3>& x1,
         const std::array<float,3>& x2) const;
/// calculate pbc distance of two atoms
   std::array<float,3> pbcDistance(size_t i1, size_t i2) const;
   float pbcDistance2(size_t i1, size_t i2) const;
/// get things
   size_t getNatoms() const;
   const std::vector<std::string>& getAtomnames() const;
   const std::vector<std::string>& getResnames() const;
   const std::vector<std::string>& getSegnames() const;
   const std::vector<std::string>& getAtomtypes() const;
   const std::vector<size_t>& getLinenumbers() const;
   const std::vector<char>& getChainids() const;
   const std::vector<float>& getXs() const;
   const std::vector<float>& getYs() const;
   const std::vector<float>& getZs() const;
   const std::vector<float>& getOccs() const;
   const std::vector<float>& getTempfs() const;
   const std::string& getAtomname(size_t index) const;
   const std::string& getResname(size_t index) const;
   const std::string& getSegname(size_t index) const;
   const std::string& getAtomtype(size_t index) const;
   size_t getLinenumber(size_t index) const;
   char getChainid(size_t index) const;
   float getX(size_t index) const;
   float getY(size_t index) const;
   float getZ(size_t index) const;
   float getOcc(size_t index) const;
   float getTempf(size_t index) const;
   std::array<float,3> getCoordinates(size_t index) const;
/// set things
   bool setAtomname(size_t index, const std::string & s);
   bool setResid(size_t index, const int i);
   bool setResname(size_t index, const std::string & s);
   bool setSegname(size_t index, const std::string & s);
   bool setAtomtype(size_t index, const std::string & s);
   bool setChainid(size_t index, char c);
   bool setX(size_t index, float f);
   bool setY(size_t index, float f);
   bool setZ(size_t index, float f);
   bool setOcc(size_t index, float f);
   bool setTempf(size_t index, float f);
/// guess one chainid according to segname; return false if segname undefined
   bool guessOneChainid(const size_t index);
/// guess all the chainids; return false if anyone fails
   bool guessAllChainids();
/// guess one segname according to chainid; return false if chainid undefined
   bool guessOneSegname(const size_t index);
/// guess all the segnames; return false if anyone fails
   bool guessAllSegnames();
/// guess one atomtype according to atomname; return false if atomname undefined
   bool guessOneAtomtype(const size_t index);
/// guess all the atomtypes; return false if anyone fails
   bool guessAllAtomtypes();
/// swap specific fields for two atoms
   void swapFields(const size_t i1, const size_t i2,
         const std::vector<PDBField>& fields); //swap multi
   void swapFields(const size_t i1, const size_t i2, const PDBField& field);
   void swapFields(const size_t i1, const size_t i2); //swap all
/// swap coordinates of two atoms
   void swapCoordinates(const size_t i1, const size_t i2);
/// check whether an atom matches the def in atomname resname segname atomtype
/// chainid resid
   bool isMatched(size_t index, const PDBDef& def) const;
   bool isMatched(const std::string& s, const PDBDef& def, PDBField f) const;
   bool isMatched(char c, const PDBDef& def, PDBField f) const;
   bool isMatched(int n, const PDBDef& def, PDBField f) const;
   bool isMatched(float x, const PDBDef& def, PDBField f) const;
/// check whehter the given string field could match
   //bool isStrMatched(size_t index, const PDBDef& def, PDBField f) const;
   //bool isStrMatched(size_t index, const PDBDef& def) const;

   //bool isAtomnameMatched(size_t index, const PDBDef& def) const;
/// reorder water atoms
/// useful with a poorly made PDB
   size_t reorderWater(bool guess, bool check, bool reorder,
         const PDBDef& defo, const PDBDef& defh, const PDBDef& defhyd);
/// useful with a very strict PDB (noguess, nocheck, noreorder)
   size_t reorderWater(const PDBDef& defo, 
         const PDBDef& defh, const PDBDef& defhyd); 
/// move an atom (i1) to the given position (i2) keeping the order of others
   bool moveTo(const size_t i1, const size_t i2);
/// assemble water such that they are in the order OHHOHHOHH...OHHH
   bool assembleWater( bool guess,
         bool check, const PDBDef& defo, const PDBDef& defh);
};

inline
const std::vector<std::string>& PDB::getAtomnames() const 
{
   return atomnames;
}

inline
const std::vector<std::string>& PDB::getResnames() const
{
   return resnames;
}

inline
const std::vector<std::string>& PDB::getSegnames() const 
{
   return segnames;
}

inline
const std::vector<std::string>& PDB::getAtomtypes() const
{ 
   return atomtypes;
}

inline
const std::vector<size_t>& PDB::getLinenumbers() const
{
   return linenumbers;
}

inline
const std::vector<char>& PDB::getChainids() const
{
   return chainids;
}

inline
const std::vector<float>& PDB::getXs() const
{
   return xs;
}

inline
const std::vector<float>& PDB::getYs() const
{
   return ys;
}

inline
const std::vector<float>& PDB::getZs() const
{
   return zs;
}

inline
const std::vector<float>& PDB::getOccs() const
{
   return occs;
}

inline
const std::vector<float>& PDB::getTempfs() const
{
   return tempfs;
}

inline
const std::string& PDB::getAtomname(size_t index) const
{
   return atomnames[index];
}

inline
const std::string& PDB::getResname(size_t index) const
{
   return resnames[index];
}

inline
const std::string& PDB::getSegname(size_t index) const 
{
   return segnames[index];
}

inline
const std::string& PDB::getAtomtype(size_t index) const
{ 
   return atomtypes[index];
}

inline
size_t PDB::getLinenumber(size_t index) const
{
   return linenumbers[index];
}

inline
char PDB::getChainid(size_t index) const
{
   return chainids[index];
}

inline
float PDB::getX(size_t index) const
{
   return xs[index];
}

inline
float PDB::getY(size_t index) const
{
   return ys[index];
}

inline
float PDB::getZ(size_t index) const
{
   return zs[index];
}

inline
float PDB::getOcc(size_t index) const
{
   return occs[index];
}

inline
float PDB::getTempf(size_t index) const
{
   return tempfs[index];
}

inline
bool PDB::setAtomname(size_t index, const std::string& s)
{
   if(index >= nAtoms)
      return false;
   else {
      atomnames[index] = s;
      defineds[index][static_cast<size_t>(PDBField::atomname)] = true;
      return true;
   }
}

inline
bool PDB::setResid(size_t index, const int i)
{
   if(index >= nAtoms)
      return false;
   else {
      resids[index] = i;
      defineds[index][static_cast<size_t>(PDBField::resid)] = true;
      return true;
   }
}

inline
bool PDB::setResname(size_t index, const std::string& s)
{
   if(index >= nAtoms)
      return false;
   else {
      resnames[index] = s;
      defineds[index][static_cast<size_t>(PDBField::resname)] = true;
      return true;
   }
}

inline
bool PDB::setSegname(size_t index, const std::string& s)
{
   if(index >= nAtoms)
      return false;
   else {
      segnames[index] = s;
      defineds[index][static_cast<size_t>(PDBField::segname)] = true;
      return true;
   }
}

inline
bool PDB::setAtomtype(size_t index, const std::string& s)
{
   if(index >= nAtoms)
      return false;
   else {
      atomtypes[index] = s;
      defineds[index][static_cast<size_t>(PDBField::atomtype)] = true;
      return true;
   }
}

inline
bool PDB::setChainid(size_t index, char c)
{
   if(index >= nAtoms)
      return false;
   else {
      chainids[index] = c;
      defineds[index][static_cast<size_t>(PDBField::chainid)] = true;
      return true;
   }
}

inline
bool PDB::setX(size_t index, float f)
{
   if(index >= nAtoms)
      return false;
   else {
      xs[index] = f;
      defineds[index][static_cast<size_t>(PDBField::x)] = true;
      return true;
   }
}

inline
bool PDB::setY(size_t index, float f)
{
   if(index >= nAtoms)
      return false;
   else {
      ys[index] = f;
      defineds[index][static_cast<size_t>(PDBField::y)] = true;
      return true;
   }
}

inline
bool PDB::setZ(size_t index, float f)
{
   if(index >= nAtoms)
      return false;
   else {
      zs[index] = f;
      defineds[index][static_cast<size_t>(PDBField::z)] = true;
      return true;
   }
}

inline
bool PDB::setOcc(size_t index, float f)
{
   if(index >= nAtoms)
      return false;
   else {
      occs[index] = f;
      defineds[index][static_cast<size_t>(PDBField::occ)] = true;
      return true;
   }
}

inline
bool PDB::setTempf(size_t index, float f)
{
   if(index >= nAtoms)
      return false;
   else {
      tempfs[index] = f;
      defineds[index][static_cast<size_t>(PDBField::tempf)] = true;
      return true;
   }
}

inline
float PDB::pbcDiff(double x1, double x2, int dim) const
{
   float diff = x2 - x1;
   int n = static_cast<int>(std::floor(diff / boxlens[dim] + 0.5));
   //int n = 0;
   return diff - n * boxlens[dim];
}

inline bool PDB::isMatched(const std::string& s, const PDBDef& def, PDBField f) const {
   //TODO is there a way to get corresponding vector name by f?
   auto range = def.getDefstr().equal_range(f);
   if(range.first == range.second) return true;
   for(auto iter = range.first; iter != range.second; ++iter) {
      if(s == iter->second) return true;
   }
   return false;
}

inline bool PDB::isMatched(char c, const PDBDef& def, PDBField f) const {
   auto range = def.getDefchr().equal_range(f);
   if(range.first == range.second) return true;
   for(auto iter = range.first; iter != range.second; ++iter) {
      if(c == iter->second) return true;
   }
   return false;
}

inline bool PDB::isMatched(int n, const PDBDef& def, PDBField f) const {
   auto range = def.getDefint().equal_range(f);
   if(range.first == range.second) return true;
   for(auto iter = range.first; iter != range.second; ++iter) {
      if(n == iter->second) return true;
   }
   return false;
}

inline bool PDB::isMatched(float x, const PDBDef& def, PDBField f) const {
   auto range = def.getDefflt().equal_range(f);
   if(range.first == range.second) return true;
   for(auto iter = range.first; iter != range.second; ++iter) {
      if(x == iter->second) return true;
   }
   return false;
}

inline float PDB::pbcDistance2(size_t i1, size_t i2) const
{
   auto x = pbcDistance(i1, i2);
   return x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
}

inline size_t PDB::getNatoms() const
{
   return nAtoms;
}

inline 
std::array<float,3> PDB::getCoordinates(size_t index) const
{
   return std::array<float,3>({xs[index], ys[index], zs[index]});
}

}
#endif
