#include<string>
#include<vector>

namespace PDB_NS {

/// line size for complete atomtype
const size_t minLineSize = 78;
const size_t pdbLineSize = 78;
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
   float boxlens[3];

   //void eraseSpace(std::string& str);

public:
   PDB() = default;
   explicit PDB(const std::string& fname);
   void write2file(const std::string& fname) const;
};

}
