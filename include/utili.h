#ifndef __PDB_UTILI_H
#define __PDB_UTILI_H
#include <string>
#include <sstream>
#include <algorithm>
#include <locale>
#include <array>
#include <vector>
#include <set>
#include <exception>

#ifdef _HAS_BOOST_PYTHON
#include <boost/python.hpp>
#endif

namespace PDB_NS {

template <class T>
inline bool readField(const std::string& s, T& t) 
{
   if(s.find_first_not_of(" \t\r\n\f\v") != std::string::npos) {
      std::stringstream ss(s);
      ss >> t; 
      return true;
   } else {
      return false;
   }
}

class Vector {
   std::array<float,3> d;

public:
   Vector() = default;
/// create a vector with same value
   Vector(float d0);
/// create a vector by giving three values
   Vector(float d0, float d1, float d2);
#ifdef _HAS_BOOST_PYTHON
/// create a vector from python::list
   Vector(const boost::python::list& bplist);
#endif
/// braket access
   float & operator[] (unsigned i);
   const float& operator[] (unsigned i) const;
/// increment
   Vector& operator+= (const Vector& v);
/// decrement
   Vector& operator-= (const Vector& v);
/// scale
   Vector& operator*= (double s);
   Vector& operator/= (double s);
/// sign -
   Vector operator- () const;
/// v1 + v2
   friend Vector operator+ (const Vector&, const Vector&);
/// v1 - v2
   friend Vector operator- (const Vector&, const Vector&);
/// s * v
   friend Vector operator* (double, const Vector&);
/// v * s
   friend Vector operator* (const Vector&, double);
/// v / s
   friend Vector operator/ (const Vector&, double);
/// v1 dot v2
   friend float dotProduct(const Vector&, const Vector&);
/// compute modulo square
   float modulo2() const;
/// compute modulo 
   float modulo() const;
   int size() const {return 3; }
};

template <class Node, typename F>
void DFS(Node root, int depth, F&& isAdj, 
      const std::vector<Node>& list, std::set<Node>& results)
{
   if(depth <= 0) return;
   for(const Node& i : list) 
      if(isAdj(root, i)) {
         results.insert(i);
         DFS(i, depth - 1, isAdj, list, results);
      }
}

//template <>
//inline bool readField(const std::string& s, int& t) 
//{
//   if(s.find_first_not_of(" \t\r\n\f\v") != std::string::npos) {
//      t = std::stoi(s); 
//      return true;
//   } else {
//      return false;
//   }
//}
//
//template <>
//inline bool readField(const std::string& s, float& t) 
//{
//   if(s.find_first_not_of(" \t\r\n\f\v") != std::string::npos) {
//      t = std::stof(s); 
//      return true;
//   } else {
//      return false;
//   }
//}
//
//template <>
//inline bool readField(const std::string& s, std::string& t) 
//{
//   t = s;
//   t.erase(std::remove_if(t.begin(),t.end(),::isspace),t.end());
//   if(t=="\0") {
//      t = " "; return false;
//   } else {
//      return true;
//   }
//}


}
#endif
