#ifndef __LIBPDB_H
#define __LIBPDB_H
#include <string>
#include <vector>
#include <array>
#include <cmath>
#include <unordered_map>
#include <utility>
#include <exception>
#include "pdbdef.h"
#include "global.h"
#include "utili.h"

#ifdef _HAS_BOOST_PYTHON
#include <boost/python.hpp>
#endif

namespace PDB_NS {

class PDB {
//pdb file format:
//ATOM serial atomname resname chainid resid x y z occ tempf segname atomtytpe
   std::vector<std::string> atomnames, resnames, segnames, atomtypes;
   std::vector<size_t> linenumbers;
   //std::vector<char> chainids;
   std::vector<std::string> chainids;
   std::vector<int> resids;
/// Global version of resids
   std::vector<int> residues;
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
/// check whether all atoms have defined fields required by def
   void checkDefined(const PDBDef& def) const;
   //std::string transField(const PDBField& pdbfield) const;

   //void eraseSpace(std::string& str);

public:
   PDB() = default;
   explicit PDB(const std::string& fname);
   void write2file(const std::string& fname) const;
   void write2file(const std::string& fname, const PDBDef& def) const;
   void write2file(const PDBDef& def, 
         const std::string& fname) const {write2file(fname, def); }
   void write2file(FILE* fp) const;
   void write2file(FILE* fp, const PDBDef& def) const;
   void write2file(FILE* fp, const std::vector<size_t>& indexes) const;
   void write2file(const std::string& fname, 
         const std::vector<size_t>& indexes) const;
   void write2file(const std::vector<size_t>& indexes,
         const std::string& fname) const {write2file(fname, indexes); }
   void show() const { write2file(stdout); }
   void show(const PDBDef& def) const { write2file(stdout, def); }
/// get all the undefined pairs(atom_index, str(pdbfiled))
   //std::vector<std::pair<size_t,std::string>> checkUndefined() const;
   std::vector<std::pair<size_t,PDBField>> checkUndefined() const;
/// calculate pbc distance pbc(x2-x1)
   Vector pbcDistance(const Vector& x1, const Vector& x2) const;
/// calculate pbc distance of two atoms
   Vector pbcDistance(size_t i1, size_t i2) const;
   float pbcDistance2(size_t i1, size_t i2) const;
   float pbcDistance2(const Vector& x1, const Vector& x2) const;
/// cos(Angle(i1-i2-i3))
   float cosAngle(size_t i1, size_t i2, size_t i3) const;
/// calculate pbc distance of one atom and a group atoms
   std::pair<float,size_t> 
      pbcDistance2(size_t i, std::vector<size_t> group) const;
/// get things
   size_t getNatoms() const;
   const std::vector<std::string>& getAtomnames() const;
   const std::vector<std::string>& getResnames() const;
   const std::vector<std::string>& getSegnames() const;
   const std::vector<std::string>& getAtomtypes() const;
   const std::vector<int>& getResids() const;
   const std::vector<int>& getResidues() const;
   const std::vector<size_t>& getLinenumbers() const;
   //const std::vector<char>& getChainids() const;
   const std::vector<std::string>& getChainids() const;
   const std::vector<float>& getXs() const;
   const std::vector<float>& getYs() const;
   const std::vector<float>& getZs() const;
   const std::vector<float>& getOccs() const;
   const std::vector<float>& getTempfs() const;
   int getResid(size_t index) const;
   const std::string& getAtomname(size_t index) const;
   const std::string& getResname(size_t index) const;
   const std::string& getSegname(size_t index) const;
   const std::string& getAtomtype(size_t index) const;
   size_t getLinenumber(size_t index) const;
   //char getChainid(size_t index) const;
   const std::string& getChainid(size_t index) const;
   float getX(size_t index) const;
   float getY(size_t index) const;
   float getZ(size_t index) const;
   float getOcc(size_t index) const;
   float getTempf(size_t index) const;
   float getBox(int i) const;
   Vector getCoordinates(size_t index) const;
   std::pair<Vector,Vector> getBoundary() const;
/// set things
   void setAtomname(size_t index, const std::string & s);
   void setResid(size_t index, const int i);
   void setResname(size_t index, const std::string & s);
   void setSegname(size_t index, const std::string & s);
   void setAtomtype(size_t index, const std::string & s);
   //bool setChainid(size_t index, char c);
   void setChainid(size_t index, const std::string & s);
   void setX(size_t index, float f);
   void setY(size_t index, float f);
   void setZ(size_t index, float f);
   void setOcc(size_t index, float f);
   void setTempf(size_t index, float f);
   void setCoordinates(size_t index, const Vector& v) {
      setX(index, v[0]); setY(index, v[1]); setZ(index, v[2]);
   }
   void setSegname(const PDBDef& def, const std::string & s);
   //bool setChainid(const PDBDef& def, char c);
   void setChainid(const PDBDef& def, const std::string & c);
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
   //bool isMatched(char c, const PDBDef& def, PDBField f) const;
   bool isMatched(int n, const PDBDef& def, PDBField f) const;
   bool isMatched(float x, const PDBDef& def, PDBField f) const;
/// check whehter the given string field could match
   //bool isStrMatched(size_t index, const PDBDef& def, PDBField f) const;
   //bool isStrMatched(size_t index, const PDBDef& def) const;

   //bool isAtomnameMatched(size_t index, const PDBDef& def) const;
/// reorder water atoms
/// useful with a poorly made PDB
   size_t reorderWater(bool guess, bool check, bool assemble,
         const PDBDef& defo, const PDBDef& defh, const PDBDef& defhyd);
/// useful with a very strict PDB (noguess, nocheck, noreorder)
   size_t reorderWater(const PDBDef& defo, 
         const PDBDef& defh, const PDBDef& defhyd); 
   size_t reorderWaterFast(bool guess, bool check, bool assemble,
         const PDBDef& defo, const PDBDef& defh, const PDBDef& defhyd);
/// move an atom (i1) to the given position (i2) keeping the order of others
   bool moveTo(const size_t i1, const size_t i2);
/// assemble water such that they are in the order OHHOHHOHH...OHHH
   void assembleWater( bool guess,
         bool check, const PDBDef& defo, const PDBDef& defh);
   void assembleWater(const PDBDef& defo, const PDBDef& defh)
       {assembleWater(false, false, defo, defh); }
/// select atoms according to PDBDef; return vector of indexes
   std::vector<size_t> selectAtoms(const PDBDef& def) const;
/// select atoms within certain distance
   std::vector<size_t> atomsWithin(const Vector& cen, float r) const;
#ifdef _HAS_BOOST_PYTHON
   std::vector<size_t> atomsWithin(const boost::python::list& cen, float r) 
      const {return atomsWithin(Vector(cen), r); }
#endif
/// compute geometric center of a group of atoms according to indexes
   Vector geoCenter(const std::vector<size_t>&) const;
/// compute geometric center of a group of atoms according to PDBDef
   Vector geoCenter(const PDBDef& def) const;
/// shift all atoms by a vector keeping the box boundaries unchanged
   void shiftBy(const Vector& offset);
/// shift all atoms by a vector so that the center of selected group atom
/// will go to the middle of the box while keeping the original box
/// boundaries unchanged
   Vector shiftToMiddle(const std::vector<size_t>&);
   Vector shiftToMiddle(const PDBDef& def) {
                            return shiftToMiddle(selectAtoms(def)); }
/// write groups to a gmx index file
   void writeIndexFile(const std::string& fname, 
                 const std::vector<Group>& grps) const;
/// write all atoms to a gmx index file
   void writeIndexFile(const std::string& fname, 
                     const std::string& grpname) const;
/// wrap a pdb by given lx hx ly hy lz hz
   void pbcWrap(float lx, float hx, float ly, float hy, float lz, float hz);
/// write groups to a xyz file
   void writeXYZ(const std::string& fname, const Group& grp) const;
/// write all atoms to a xyz file
   void writeXYZ(const std::string& fname, const std::string& grpname) const;
   void printOneAtom(FILE* fp, size_t index) const;
   void printAtoms(FILE* fp, const PDBDef& def) const;
   void printAtoms(FILE* fp, const std::vector<size_t>& indexes) const;
   void printOneAtom(size_t index) const {printOneAtom(stdout, index); }
   void printAtoms(const PDBDef& def) const {printAtoms(stdout, def); }
   void printAtoms(const std::vector<size_t>& indexes) 
                   const {printAtoms(stdout, indexes); }
/// Get n solvation shells of water oxygen of a hydronium (no reoderWater
/// will be done). Returning a vector containing indexes of both oxygen and 
/// hydrogen instead of returning a vector of WaterNode is for the convenience
/// in Python interface. WE NEED OVERLOAD BY RETURN TYPE!!!
   std::vector<size_t> getSolvationShells(int n, float cutoff, 
         const std::vector<size_t>& oindexes, size_t hydindex, int, bool) ;
   std::vector<size_t> getSolvationShells(int n, float cutoff, 
         const PDBDef& defo, const PDBDef& defhyd, int, bool make_whole) ;
/// Get the HB distance of two water nodes
   std::pair<double,double> adjacencyWaterNode(
                WaterNode n1, WaterNode n2) const;
/// Check if there is a HB from n1 to n2 as VMD does
   bool isWaterNodeHBonded (WaterNode n1, WaterNode n2, 
         float roo, float theta) const;
/// Get HB acceptors or donors according to EVB state searching
   std::vector<size_t> getHBAcceptors(int n, float cutoff, 
         const PDBDef& defo, const PDBDef& defhyd)  {
      return getSolvationShells(n, cutoff, defo, defhyd, 0, true); }
   std::vector<size_t> getHBAcceptors(int n, float cutoff, 
         const std::vector<size_t>& oindexes, size_t hydindex)  {
      return getSolvationShells(n, cutoff, oindexes, hydindex, 0, true); }
   std::vector<size_t> getHBDonors(int n, float cutoff, 
         const PDBDef& defo, const PDBDef& defhyd)  {
      return getSolvationShells(n, cutoff, defo, defhyd, 1, true); }
   std::vector<size_t> getHBDonors(int n, float cutoff, 
         const std::vector<size_t>& oindexes, size_t hydindex)  {
      return getSolvationShells(n, cutoff, oindexes, hydindex, 1, true); }
/// Get HB network using DA distance and DHA angle criteria
   std::vector<size_t> getHBNetwork(int n, float roo, float theta,
      size_t root,
      const std::vector<size_t>& ocindexes, const std::vector<size_t>& owcindexes,
      size_t hydindex, Direction d, bool make_whole);
   std::vector<size_t> getHBNetwork(int n, float roo, float theta,
      const PDBDef& defroot,
      const PDBDef& defoc, const PDBDef& defow, const PDBDef& defhyd, 
      Direction d, bool make_whole);
/// make atom 2 in the same pbc images of atom1
   void make_connect(size_t i1, size_t i2) {
      setCoordinates(i2, getCoordinates(i1) + pbcDistance(i1, i2));
   }
};

inline
const std::vector<std::string>& PDB::getAtomnames() const 
{
   return atomnames;
}

inline
const std::vector<int>& PDB::getResids() const
{
   return resids;
}

inline
const std::vector<int>& PDB::getResidues() const
{
   return residues;
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
const std::vector<std::string>& PDB::getChainids() const
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
int PDB::getResid(size_t index) const
{
   return resids[index];
}

inline
const std::string& PDB::getChainid(size_t index) const
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
void PDB::setAtomname(size_t index, const std::string& s)
{
   if(index >= nAtoms)
      throw std::out_of_range("index out of range in setAtomname");
   else {
      atomnames[index] = s;
      defineds[index][static_cast<size_t>(PDBField::atomname)] = true;
   }
}

inline
void PDB::setResid(size_t index, const int i)
{
   if(index >= nAtoms)
      throw std::out_of_range("index out of range in setResid");
   else {
      resids[index] = i;
      defineds[index][static_cast<size_t>(PDBField::resid)] = true;
   }
}

inline
void PDB::setResname(size_t index, const std::string& s)
{
   if(index >= nAtoms)
      throw std::out_of_range("index out of range in setResname");
   else {
      resnames[index] = s;
      defineds[index][static_cast<size_t>(PDBField::resname)] = true;
   }
}

inline
void PDB::setSegname(size_t index, const std::string& s)
{
   if(index >= nAtoms)
      throw std::out_of_range("index out of range in setSegname");
   else {
      segnames[index] = s;
      defineds[index][static_cast<size_t>(PDBField::segname)] = true;
   }
}

inline
void PDB::setAtomtype(size_t index, const std::string& s)
{
   if(index >= nAtoms)
      throw std::out_of_range("index out of range in setAtomtype");
   else {
      atomtypes[index] = s;
      defineds[index][static_cast<size_t>(PDBField::atomtype)] = true;
   }
}

inline
void PDB::setChainid(size_t index, const std::string& c)
{
   if(c.length() != 1)
      throw std::invalid_argument("a string of length 1 is expected");
   if(index >= nAtoms)
      throw std::out_of_range("index out of range in setChainid");
   else {
      chainids[index] = c;
      defineds[index][static_cast<size_t>(PDBField::chainid)] = true;
   }
}

inline
void PDB::setX(size_t index, float f)
{
   if(index >= nAtoms)
      throw std::out_of_range("index out of range in setX");
   else {
      xs[index] = f;
      defineds[index][static_cast<size_t>(PDBField::x)] = true;
   }
}

inline
void PDB::setY(size_t index, float f)
{
   if(index >= nAtoms)
      throw std::out_of_range("index out of range in setY");
   else {
      ys[index] = f;
      defineds[index][static_cast<size_t>(PDBField::y)] = true;
   }
}

inline
void PDB::setZ(size_t index, float f)
{
   if(index >= nAtoms)
      throw std::out_of_range("index out of range in setZ");
   else {
      zs[index] = f;
      defineds[index][static_cast<size_t>(PDBField::z)] = true;
   }
}

inline
void PDB::setOcc(size_t index, float f)
{
   if(index >= nAtoms)
      throw std::out_of_range("index out of range in setOcc");
   else {
      occs[index] = f;
      defineds[index][static_cast<size_t>(PDBField::occ)] = true;
   }
}

inline
void PDB::setTempf(size_t index, float f)
{
   if(index >= nAtoms)
      throw std::out_of_range("index out of range in setTempf");
   else {
      tempfs[index] = f;
      defineds[index][static_cast<size_t>(PDBField::tempf)] = true;
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

inline bool PDB::isMatched(
      const std::string& s, const PDBDef& def, PDBField f) const {
   //TODO is there a way to get corresponding vector name by f?
   auto range = def.getDefstr().equal_range(f);
   if(range.first == range.second) return true;
   for(auto iter = range.first; iter != range.second; ++iter) {
      if(s == iter->second) return true;
   }
   return false;
}

//inline bool PDB::isMatched(char c, const PDBDef& def, PDBField f) const {
//   auto range = def.getDefchr().equal_range(f);
//   if(range.first == range.second) return true;
//   for(auto iter = range.first; iter != range.second; ++iter) {
//      if(c == iter->second) return true;
//   }
//   return false;
//}

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
   //return x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
   return x.modulo2();
}

inline float PDB::pbcDistance2(const Vector& x1, const Vector& x2) const
{
   auto x = pbcDistance(x1, x2);
   //return x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
   return x.modulo2();
}

inline size_t PDB::getNatoms() const
{
   return nAtoms;
}

inline 
Vector PDB::getCoordinates(size_t index) const
{
   return Vector(xs[index], ys[index], zs[index]);
}

inline
float PDB::getBox(int i) const
{
   return boxlens[i];
}

inline
Vector PDB::geoCenter(const PDBDef& def) const
{
   return geoCenter(selectAtoms(def));
}

inline
void PDB::write2file(FILE* fp, const PDBDef& def) const 
{
   write2file(fp, selectAtoms(def));
}

inline
void PDB::write2file(const std::string& fname) const
{
   FILE *fp = fopen(fname.c_str(),"w");
   if(!fp) {
      throw std::runtime_error("cannot open " + fname);
   }
   write2file(fp);
}

inline
void PDB::write2file(FILE* fp) const
{
   std::vector<size_t> indexes;
   for(size_t i = 0; i < nAtoms; ++i)
      indexes.push_back(i);
   write2file(fp, indexes);
}

inline
void PDB::write2file(const std::string& fname, const PDBDef& def) const 
{
   FILE *fp = fopen(fname.c_str(),"w");
   if(!fp) {
      throw std::runtime_error("cannot open " + fname);
   }
   write2file(fp, def);
}

inline
void PDB::write2file(const std::string& fname, 
      const std::vector<size_t>& indexes) const
{
   FILE *fp = fopen(fname.c_str(),"w");
   if(!fp) {
      throw std::runtime_error("cannot open " + fname);
   }
   write2file(fp, indexes);
}

} // end of namespace PDB_NS

#endif
