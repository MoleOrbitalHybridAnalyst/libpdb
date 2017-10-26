#ifndef __LIBPDB_GENERAL_H
#define __LIBPDB_GENERAL_H

/// line size for complete atomtype
const size_t minLineSize = 78;
const size_t pdbLineSize = 78;
enum PDBField { atomname = 0, resname = 1, segname = 9, atomtype = 10, 
      chainid = 2, resid = 3, x = 4, y = 5, z = 6, occ = 7, tempf = 8};
/// line size for complete z coordinate
//const size_t errLineSize = 54;
//const bool autoComplete = true;

#endif
