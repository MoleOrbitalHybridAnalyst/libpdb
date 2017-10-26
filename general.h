#ifndef __LIBPDB_GENERAL_H
#define __LIBPDB_GENERAL_H

namespace PDB_NS {
/// line size for complete atomtype
const size_t minLineSize = 78;
const size_t pdbLineSize = 78;
//enum class PDBField : size_t { atomname = 0, resname = 1, segname = 9, atomtype = 10, 
enum  PDBField { atomname = 0, resname = 1, segname = 9, atomtype = 10, 
      chainid = 2, resid = 3, x = 4, y = 5, z = 6, occ = 7, tempf = 8};
/// line size for complete z coordinate
//const size_t errLineSize = 54;
//const bool autoComplete = true;

}

namespace std {

template<>
struct hash<PDB_NS::PDBField>
{
   typedef PDB_NS::PDBField argument_type;
   typedef size_t result_type;

   result_type operator () (const argument_type& x) const
   {
      using type = typename std::underlying_type<argument_type>::type;
      return std::hash<type>()(static_cast<type>(x));
   }
};

}

#endif
