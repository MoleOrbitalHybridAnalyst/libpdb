#include <string>
#include <sstream>
namespace PDB_NS {

template <class T>
bool readField(const std::string& s, T& t) {
   if(s.find_first_not_of(" \t\r\n\f\v") != std::string::npos) {
      std::stringstream ss(s);
      ss >> t; 
      return true;
   } else {
      return false;
   }
}

}
