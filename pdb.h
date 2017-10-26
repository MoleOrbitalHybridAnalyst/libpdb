#ifndef __LIBPDB_H
#define __LIBPDB_H
#include<string>
#include<vector>

namespace PDB_NS {

/// line size for complete atomtype
const size_t minLineSize = 78;
const size_t pdbLineSize = 78;
enum PDBField { atomname = 0, resname = 1, segname = 9, atomtype = 10, 
      chainid = 2, resid = 3, x = 4, y = 5, z = 6, occ = 7, tempf = 8};
/// line size for complete z coordinate
//const size_t errLineSize = 54;
//const bool autoComplete = true;

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
/// reorder 
};

}
#endif
