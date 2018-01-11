#include "utili.h"

namespace PDB_NS {

Vector::Vector(float d0) {
   for(int i = 0 ; i < 3; ++i) {
      d[i] = d0;
   }
}

Vector::Vector(float d0, float d1, float d2)
{
   d[0] = d0;
   d[1] = d1;
   d[2] = d2;
}

float& Vector::operator[] (unsigned i)
{
   return d[i];
}

Vector& Vector::operator+= (const Vector& v)
{
   for(int i = 0; i < 3; ++i) {
      d[i] += v.d[i];
   }
   return *this;
}

Vector& Vector::operator-= (const Vector& v)
{
   for(int i = 0; i < 3; ++i) {
      d[i] -= v.d[i];
   }
   return *this;
}

Vector& Vector::operator*= (double s)
{
   for(int i = 0; i < 3; ++i) {
      d[i] *= s;
   }
   return *this;
}

Vector& Vector::operator/= (double s)
{
   for(int i = 0; i < 3; ++i) {
      d[i] /= s;
   }
   return *this;
}

Vector Vector::operator- () const
{
   return Vector(-d[0],-d[1],-d[2]);
}

Vector operator+ (const Vector& v1, const Vector& v2)
{
   Vector v(v1);
   return v += v2;
}

Vector operator- (const Vector& v1, const Vector& v2)
{
   Vector v(v1);
   return v -= v2;
}

Vector operator* (double s, const Vector& v)
{
   Vector _v(v);
   return _v *= s;
}

Vector operator* (const Vector& v, double s)
{
   Vector _v(v);
   return _v *= s;
}

Vector operator/ (const Vector& v, double s)
{
   Vector _v(v);
   return _v /= s;
}

float dotProduct(const Vector& v1, const Vector& v2)
{
   float x = 0.0f;
   for(int i = 0; i < 3; ++i)
      x += (v1.d[i] * v2.d[i]);
   return x;
}

float Vector::modulo2() const
{
   return dotProduct(*this, *this);
}

float Vector::modulo() const
{
   return std::sqrt(modulo2());
}

} 
