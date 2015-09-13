/*
 *
 *   PURPOSE: To declare a class for a tensor of rank 6, mostly for storage purposes.
 *            It's a class instead of a struct, as I'm currently undecided about
 *            how superficial/deep it will be.
 *
 */

#ifndef TENSOR6HEADERDEF
#define TENSOR6HEADERDEF

#include "matrix.hpp"

class Tensor6
{
private:
  Matrix data;
  int u, v, w, x, y, z;
public:
  Tensor6() { } 
  Tensor6(int a, int b, int c, int d, int e, int f);
  Tensor6(int a, int b, int c, int d, int e, int f, double val);
  Tensor6(const Tensor6& other);
  void resize(int a, int b, int c, int d, int e, int f);
  void assign(int a, int b, int c, int d, int e, int f, double val);
  void print() const;
  double& operator()(int i, int j, int k, int l, int m, int n);
  double operator()(int i, int j, int k, int l, int m, int n) const;
  Tensor6& operator=(const Tensor6& other);
  Tensor6 operator+(const Tensor6& other) const;
  Tensor6& operator*=(double scalar) { return *this; }
  void multiply(double scalar);
  Matrix& getData() { return data; }
  int getu() const { return u; }
  int getv() const { return v; }
  int getw() const { return w; }
  int getx() const { return x; }
  int gety() const { return y; }
  int getz() const { return z; }
};

// Scalar multiplication
inline Tensor6 operator*(const Tensor6& tens, double scalar)
{
  Tensor6 rtens;
  rtens = tens;
  rtens.multiply(scalar);
  return rtens;
}

inline Tensor6 operator*(double scalar, const Tensor6& tens)
{
  Tensor6 rtens;
  rtens= tens;
  rtens.multiply(scalar);
  return rtens;
}


#endif
