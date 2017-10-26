#ifndef __PDB_UTILI_H
#define __PDB_UTILI_H
#include <string>
#include <sstream>
#include <algorithm>
#include <locale>
namespace PDB_NS {

template <class T>
inline bool readField(const std::string& s, T& t) {
   if(s.find_first_not_of(" \t\r\n\f\v") != std::string::npos) {
      std::stringstream ss(s);
      ss >> t; 
      return true;
   } else {
      return false;
   }
}

template <>
inline bool readField(const std::string& s, int& t) {
   if(s.find_first_not_of(" \t\r\n\f\v") != std::string::npos) {
      t = std::stoi(s); 
      return true;
   } else {
      return false;
   }
}

template <>
inline bool readField(const std::string& s, float& t) {
   if(s.find_first_not_of(" \t\r\n\f\v") != std::string::npos) {
      t = std::stof(s); 
      return true;
   } else {
      return false;
   }
}

template <>
inline bool readField(const std::string& s, std::string& t) {
   t = s;
   t.erase(std::remove_if(t.begin(),t.end(),::isspace),t.end());
   if(t=="\0") {
      t = " "; return false;
   } else {
      return true;
   }
}


}
#endif
