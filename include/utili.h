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
/// Things needed for boost::vector_indexing_suite
   typedef float value_type;
   typedef size_t size_type;
   typedef std::ptrdiff_t difference_type;
   typedef std::array<float,3>::iterator iterator;
   typedef std::array<float,3>::const_iterator const_iterator;
   iterator begin() {return d.begin();}
   iterator end() {return d.end();}
   iterator erase(iterator it) {
      (void) it;
      std::runtime_error("Vector::erase should not be called");
      return d.end();
   }
   iterator erase(iterator f, iterator l) {
      (void) f; (void) l;
      std::runtime_error("Vector::erase should not be called");
      return d.end();
   }
   iterator insert(iterator pos, const float& value) {
      (void) pos; (void) value;
      std::runtime_error("Vector::insert should not be called");
      return d.end();
   }
   template< class InputIt >
   iterator insert(const_iterator pos, InputIt f, InputIt l) {
      (void) pos; (void) f; (void) l;
      std::runtime_error("Vector::insert should not be called");
      return d.end();
   }
   void push_back(const float& value) {
      (void) value;
      std::runtime_error("Vector::push_back should not be called");
   }
   template< class InputIt >
   Vector(InputIt f, InputIt l) {
      (void) f; (void) l;
      std::runtime_error("Vector::(InputIt f, InputIt l) should not be called");
   }
/// default constructor
   Vector() = default;
/// create a vector with same value
   Vector(float d0);
/// create a vector by giving three values
   Vector(float d0, float d1, float d2);
#ifdef _HAS_BOOST_PYTHON
/// create a vector from python::list
   Vector(const boost::python::list& bplist)
   {
      int n = boost::python::len(bplist);
      if(n != 3)
         throw std::invalid_argument("input list is not 3 dimensional");
      for(int i = 0; i < n; ++i)
         d[i] = boost::python::extract<float>(bplist[i]);
   }
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
