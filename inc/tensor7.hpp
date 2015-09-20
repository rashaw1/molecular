/*
 *
 *   PURPOSE: To declare a class for a tensor of rank 7, mostly for storage purposes.
 *            It's a class instead of a struct, as I'm currently undecided about
 *            how superficial/deep it will be.
 *
 */

#ifndef TENSOR7HEADERDEF
#define TENSOR7HEADERDEF

#include "matrix.hpp"

class Tensor7
{
private:
  Matrix data;
  int t, u, v, w, x, y, z;
public:
  Tensor7() { } 
  Tensor7(int a, int b, int c, int d, int e, int f, int g);
  Tensor7(int a, int b, int c, int d, int e, int f, int g, double val);
  Tensor7(const Tensor7& other);
  void resize(int a, int b, int c, int d, int e, int f, int g);
  void assign(int a, int b, int c, int d, int e, int f, int g, double val);
  void print() const;
  double& operator()(int i, int j, int k, int l, int m, int n, int p);
  double operator()(int i, int j, int k, int l, int m, int n, int p) const;
  Tensor7& operator=(const Tensor7& other);
  Tensor7 operator+(const Tensor7& other) const;
  Tensor7& operator*=(double scalar) { return *this; }
  void multiply(double scalar);
  Matrix& getData() { return data; }
  int gett() const { return t; }
  int getu() const { return u; }
  int getv() const { return v; }
  int getw() const { return w; }
  int getx() const { return x; }
  int gety() const { return y; }
  int getz() const { return z; }
};

// Scalar multiplication
inline Tensor7 operator*(const Tensor7& tens, double scalar)
{
  Tensor7 rtens;
  rtens = tens;
  rtens.multiply(scalar);
  return rtens;
}

inline Tensor7 operator*(double scalar, const Tensor7& tens)
{
  Tensor7 rtens;
  rtens= tens;
  rtens.multiply(scalar);
  return rtens;
}


#endif
