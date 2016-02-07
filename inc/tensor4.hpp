/*
 *
 *   PURPOSE: To declare a class for a tensor of rank 4, mostly for storage purposes.
 *            It's a class instead of a struct, as I'm currently undecided about
 *            how superficial/deep it will be.
 *
 */

#ifndef TENSOR4HEADERDEF
#define TENSOR4HEADERDEF

#include "matrix.hpp"

class Tensor4
{
private:
  Matrix data;
  int w, x, y, z;
public:
  Tensor4() { } 
  Tensor4(int a, int b, int c, int d);
  Tensor4(int a, int b, int c, int d, double val);
  Tensor4(const Tensor4& other);
	int getW() const { return w; }
	int getX() const { return x; }
	int getY() const { return y; }
	int getZ() const { return z; }
	void resize(int a, int b, int c, int d);
  void assign(int a, int b, int c, int d, double val);
  void print() const;
  double& operator()(int i, int j, int k, int l);
  double operator()(int i, int j, int k, int l) const;
  Tensor4& operator=(const Tensor4& other);
  Tensor4 operator+(const Tensor4& other) const;
};

#endif
