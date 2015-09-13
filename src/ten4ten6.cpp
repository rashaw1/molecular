/*
 *
 *   PURPOSE: To implement class Ten4Ten6
 *
 */

#include "ten4ten6.hpp"

Ten4Ten6::Ten4Ten6(int m, int n, int p, int q) : w(m), x(n), y(p), z(q)
{
  arr.resize(m*n*p*q);
}

Ten4Ten6::~Ten4Ten6()
{
  
}

void Ten4Ten6::clean()
{
  arr.erase(arr.begin(), arr.end());
  w = x = y = z = 0; 
}

void Ten4Ten6::set(int i, int j, int k, int l, Tensor6 tens)
{
  arr.at(i*x*y*z+j*y*z+k*z+l) = tens;
}

Tensor6& Ten4Ten6::operator()(int i, int j, int k, int l)
{
  return arr.at(i*x*y*z+j*y*z+k*z+l);
}

Tensor6 Ten4Ten6::operator()(int i, int j, int k, int l) const
{
  return arr.at(i*x*y*z+j*y*z+k*z+l);
}
