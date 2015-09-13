/*
 *
 *   PURPOSE: To implement class Ten4Ten4
 *
 */

#include "ten4ten4.hpp"
#include "mvector.hpp"

Ten4Ten4::Ten4Ten4(int m, int n, int p, int q) : w(m), x(n), y(p), z(q)
{
  arr.resize(m*n*p*q);
}

Ten4Ten4::~Ten4Ten4()
{
 
}

void Ten4Ten4::clean()
{
  arr.erase(arr.begin(), arr.end());
  w = x = y = z = 0;
}

void Ten4Ten4::set(int i, int j, int k, int l, Tensor4 tens)
{
  arr.at(i*x*y*z+j*y*z+k*z+l) = tens;
}

Tensor4& Ten4Ten4::operator()(int i, int j, int k, int l)
{
  return arr.at(i*x*y*z+j*y*z+k*z+l);
}

Tensor4 Ten4Ten4::operator()(int i, int j, int k, int l) const
{
  return arr.at(i*x*y*z+j*y*z+k*z+l);
}
