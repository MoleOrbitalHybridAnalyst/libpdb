#ifndef __LIBPDB_H
#define __LIBPDB_H
#include <string>
#include <vector>
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
   std::string transField(const PDBField& pdbfield) const;

   //void eraseSpace(std::string& str);

public:
   PDB() = default;
   explicit PDB(const std::string& fname);
   void write2file(const std::string& fname) const;
/// get all the undefined pairs(atom_index, pdbfiled)
   std::vector<std::pair<size_t,std::string>> checkUndefined() const;
/// get things
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
   const std::string& getResnames(size_t index) const;
   const std::string& getSegnames(size_t index) const;
   const std::string& getAtomtypes(size_t index) const;
   size_t getLinenumbers(size_t index) const;
   char getChainid(size_t index) const;
   float getX(size_t index) const;
   float getY(size_t index) const;
   float getZ(size_t index) const;
   float getOcc(size_t index) const;
   float getTempf(size_t index) const;
/// set things
   bool setAtomname(size_t index, const std::string & s);
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
         const std::vector<PDBField>& fields);
   void swapFields(const size_t i1, const size_t i2, const PDBField& field);
/// swap coordinates of two atoms
   void swapCoordinates(const size_t i1, const size_t i2);
/// reorder water atoms
//TODO PDBDef class
typedef PDBField PDBDef;
/// useful with a poorly made PDB
   size_t reorderWater(bool guess, bool check, bool reorder,
         const PDBDef& defo, const PDBDef& defh);
/// useful with a very strict PDB (noguess, nocheck, noreorder)
   size_t reorderWater(const PDBDef& defo, const PDBDef& defh); 
/// useful with a very strict PDB but waters are not in the order of OHHOHH...
/// (noguess nocheck reorder)
   size_t reorderWater(const PDBDef& defo, const PDBDef& defh); 

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
   else
      atomnames[index] = s;
}

inline
bool PDB::setResname(size_t index, const std::string& s)
{
   if(index >= nAtoms)
      return false;
   else
      resnames[index] = s;
}

inline
bool PDB::setSegname(size_t index, const std::string& s)
{
   if(index >= nAtoms)
      return false;
   else
      segnames[index] = s;
}

inline
bool PDB::setAtomtype(size_t index, const std::string& s)
{
   if(index >= nAtoms)
      return false;
   else
      atomtypes[index] = s;
}

inline
bool PDB::setChainid(size_t index, char c)
{
   if(index >= nAtoms)
      return false;
   else
      chainids[index] = c;
}

inline
bool PDB::setX(size_t index, float f)
{
   if(index >= nAtoms)
      return false;
   else
      xs[index] = f;
}

inline
bool PDB::setY(size_t index, float f)
{
   if(index >= nAtoms)
      return false;
   else
      ys[index] = f;
}

inline
bool PDB::setZ(size_t index, float f)
{
   if(index >= nAtoms)
      return false;
   else
      zs[index] = f;
}

inline
bool PDB::setOcc(size_t index, float f)
{
   if(index >= nAtoms)
      return false;
   else
      occs[index] = f;
}

inline
bool PDB::setTempf(size_t index, float f)
{
   if(index >= nAtoms)
      return false;
   else
      tempfs[index] = f;
}

}
#endif
