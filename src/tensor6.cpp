/* 
 *
 *   PURPOSE: To implement class Tensor6
 *
 */

#include "tensor6.hpp"
#include <iostream>

Tensor6::Tensor6(int a, int b, int c, int d, int e, int f) : u(a), v(b), w(c), x(d), y(e), z(f)
{
  // Resize the data matrix
  data.resize(a*b*c, d*e*f);
}

Tensor6::Tensor6(int a, int b, int c, int d, int e, int f, double val) : u(a), v(b), w(c), x(d), y(e), z(f)
{
  // Resize and assign the data matrix
  data.assign(a*b*c, d*e*f, val);
}

Tensor6::Tensor6(const Tensor6& other)
{
  u = other.u; v = other.v; w = other.w;
  x = other.x; y = other.y; z = other.z;
  data.resize(u*v*w, x*y*z);
  for (int i = 0; i < u; i++){
    for (int j = 0; j < v; j++){
      for (int k = 0; k < w; k++){
        for (int l = 0; l < x; l++){
          for (int m = 0; m < y; m++){
            for (int n = 0; n < z; n++){
              data(i*v*w+j*w+k, l*y*z+m*z+n) = other(i, j, k, l, m, n);
            }
          }
        }
      }
    }
  }

}

void Tensor6::resize(int a, int b, int c, int d, int e, int f)
{
  u = a; v= b; w = c; x = d; y = e; z = f;
  data.resize(a*b*c, d*e*f);
}

void Tensor6::assign(int a, int b, int c, int d, int e, int f, double val)
{
  u = a; v= b; w = c; x = d; y = e; z =f;
  data.assign(a*b*c, d*e*f, val);
}

void Tensor6::print() const
{
  data.print(); std::cout << "\n\n";
}

double& Tensor6::operator()(int i, int j, int k, int l, int m, int n)
{
  return data(i*v*w+j*w+k, l*y*z+m*z+n);
}

double Tensor6::operator()(int i, int j, int k, int l, int m, int n) const
{
  return data(i*v*w+j*w+k, l*y*z+m*z+n);
}

Tensor6& Tensor6::operator=(const Tensor6& other)
{
  resize(other.u, other.v, other.w, other.x, other.y, other.z);
  for (int i = 0; i < u; i++){
    for (int j = 0; j < v; j++){
      for (int k = 0; k < w; k++){
        for (int l = 0; l < x; l++){
          for (int m = 0; m < y; m++){
            for (int n = 0; n < z; n++){
              data(i*v*w+j*w+k, l*y*z+m*z+n) = other(i, j, k, l, m, n);
            }
          }
        }
      }
    }
  }
  
  return *this;
}

Tensor6 Tensor6::operator+(const Tensor6& other) const
{
  Tensor6 rtens(u, v, w, x, y, z);

  for (int i = 0; i < u; i++){
    for (int j = 0; j < v; j++){
      for (int k = 0; k < w; k++){
	for (int l = 0; l < x; l++){
	  for (int m = 0; m < y; m++){
	    for (int n = 0; n < z; n++){
	      rtens(i, j, k, l, m, n) = data(i*v*w + j*w + k, l*y*z + m*z + n) + other(i, j, k, l, m, n);
	    }
	  }
	}
      }
    }
  }
  
  return rtens;
}

void Tensor6::multiply(double scalar)
{
  data = scalar*data;
}
