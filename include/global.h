#ifndef __LIBPDB_GENERAL_H
#define __LIBPDB_GENERAL_H

#include <string>
#include <vector>
#include <limits>
#include <exception>
#include <iomanip>
#include <set>

namespace PDB_NS {
/// line size for complete atomtype
const size_t minLineSize = 78;
const size_t pdbLineSize = 78;
const size_t size_tMax = std::numeric_limits<size_t>::max();
//enum class PDBField : size_t { atomname = 0, resname = 1, segname = 9, atomtype = 10, 
enum class PDBField { atomname = 0, resname = 1, segname = 9, atomtype = 10, 
      chainid = 2, resid = 3, x = 4, y = 5, z = 6, occ = 7, tempf = 8, unknown = 11};
enum class Direction { backward = -1, both = 0, forward = 1, either = 2};
typedef std::pair<std::string,std::vector<size_t>> Group;
typedef std::pair<size_t,int> WaterNode;
/// line size for complete z coordinate
//const size_t errLineSize = 54;
//const bool autoComplete = true;

inline
bool isString(PDBField f) {
   return (f==PDBField::atomname or f==PDBField::resname or f==PDBField::segname
         or f==PDBField::atomtype or f==PDBField::chainid);
}

inline
bool isString(const std::string& s) {
   return (s=="atomname" or s=="resname" or s=="segname" or s=="atomtype"
         or s=="chainid");
}

//inline bool isChar(PDBField f) {
//   return f==PDBField::chainid;
//}

//inline bool isChar(const std::string& s) {
//   return s=="chainid";
//}

inline bool isInt(PDBField f) {
   return f==PDBField::resid;
}

inline bool isInt(const std::string& s) {
   return s=="resid";
}

inline bool isFloat(PDBField f) {
   return (f==PDBField::x or f==PDBField::y or f==PDBField::z or 
         f==PDBField::occ or f==PDBField::tempf);
}

inline bool isFloat(const std::string& s) {
   return (s=="x" or s=="y" or s=="z" or s=="occ" or s=="tempf");
}

inline const std::string transField(PDBField pdbfield) {
   switch(pdbfield) {
      case PDBField::atomname: return "atomname";
      case PDBField::resname: return "resname";
      case PDBField::segname: return "segname";
      case PDBField::atomtype: return "atomtype";
      case PDBField::chainid: return "chainid";
      case PDBField::resid: return "resid";
      case PDBField::x: return "x";
      case PDBField::y: return "y";
      case PDBField::z: return "z";
      case PDBField::occ: return "occ";
      case PDBField::tempf: return "tempf";
      //default: return "\0";
      default: throw std::invalid_argument("unkown field");
   }
}

inline PDBField transString(const std::string& s) {
   if(s=="atomname") return PDBField::atomname; 
   else if(s=="resname") return PDBField::resname;
   else if(s=="segname") return PDBField::segname;
   else if(s=="atomtype") return PDBField::atomtype;
   else if(s=="chainid") return PDBField::chainid;
   else if(s=="resid") return PDBField::resid;
   else if(s=="x") return PDBField::x;
   else if(s=="y") return PDBField::y;
   else if(s=="z") return PDBField::z;
   else if(s=="occ") return PDBField::occ;
   else if(s=="tempf") return PDBField::tempf;
   //else return PDBField::unknown; 
   else throw std::invalid_argument("unknown field " + s);
}

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
