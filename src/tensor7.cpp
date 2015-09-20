/* 
 *
 *   PURPOSE: To implement class Tensor7
 *
 */

#include "tensor7.hpp"
#include <iostream>

Tensor7::Tensor7(int a, int b, int c, int d, int e, int f, int g) : t(g), u(a), v(b), w(c), x(d), y(e), z(f)
{
  // Resize the data matrix
  data.resize(a*b*c, d*e*f*g);
}

Tensor7::Tensor7(int a, int b, int c, int d, int e, int f, int g, double val) : t(g), u(a), v(b), w(c), x(d), y(e), z(f)
{
  // Resize and assign the data matrix
  data.assign(a*b*c, d*e*f*g, val);
}

Tensor7::Tensor7(const Tensor7& other)
{
  u = other.u; v = other.v; w = other.w;
  x = other.x; y = other.y; z = other.z; t = other.t;
  data.resize(u*v*w, x*y*z);
  for (int i = 0; i < u; i++){
    for (int j = 0; j < v; j++){
      for (int k = 0; k < w; k++){
        for (int l = 0; l < x; l++){
          for (int m = 0; m < y; m++){
            for (int n = 0; n < z; n++){
	      for (int p = 0; p < t; p++){
		data(i*v*w+j*w+k, l*y*z*t+m*z*t+n*t +p) = other(i, j, k, l, m, n, p);
	      }
	    }
          }
        }
      }
    }
  }

}

void Tensor7::resize(int a, int b, int c, int d, int e, int f, int g)
{
  u = a; v= b; w = c; x = d; y = e; z = f; t = g;
  data.resize(a*b*c, d*e*f*g);
}

void Tensor7::assign(int a, int b, int c, int d, int e, int f, int g, double val)
{
  u = a; v= b; w = c; x = d; y = e; z =f; t= g;
  data.assign(a*b*c, d*e*f*g, val);
}

void Tensor7::print() const
{
  data.print(); std::cout << "\n\n";
}

double& Tensor7::operator()(int i, int j, int k, int l, int m, int n, int p)
{
  return data(i*v*w+j*w+k, l*y*z*t+m*z*t+n*t+p);
}

double Tensor7::operator()(int i, int j, int k, int l, int m, int n, int p) const
{
  return data(i*v*w+j*w+k, l*y*z*t+m*z*t+n*t+p);
}

Tensor7& Tensor7::operator=(const Tensor7& other)
{
  resize(other.u, other.v, other.w, other.x, other.y, other.z, other.t);
  for (int i = 0; i < u; i++){
    for (int j = 0; j < v; j++){
      for (int k = 0; k < w; k++){
        for (int l = 0; l < x; l++){
          for (int m = 0; m < y; m++){
            for (int n = 0; n < z; n++){
	      for (int p = 0; p < t; p++){
		data(i*v*w+j*w+k, l*y*z*t+m*z*t+n*t+p) = other(i, j, k, l, m, n, p);
	      }
	    }
          }
        }
      }
    }
  }
  
  return *this;
}

Tensor7 Tensor7::operator+(const Tensor7& other) const
{
  Tensor7 rtens(u, v, w, x, y, z, t);

  for (int i = 0; i < u; i++){
    for (int j = 0; j < v; j++){
      for (int k = 0; k < w; k++){
	for (int l = 0; l < x; l++){
	  for (int m = 0; m < y; m++){
	    for (int n = 0; n < z; n++){
	      for (int p = 0; p < t; p++){
		rtens(i, j, k, l, m, n, p) = data(i*v*w + j*w + k, l*y*z*t + m*z*t + n*t + p) + other(i, j, k, l, m, n, p);
	      }
	    }
	  }
	}
      }
    }
  }
  
  return rtens;
}

void Tensor7::multiply(double scalar)
{
  data = scalar*data;
}
