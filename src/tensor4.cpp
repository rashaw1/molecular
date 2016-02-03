/* 
 *
 *   PURPOSE: To implement class Tensor4
 *
 */

#include "tensor4.hpp"
#include <iostream>

Tensor4::Tensor4(int a, int b, int c, int d) : w(a), x(b), y(c), z(d)
{
  // Resize the data matrix
  data.resize(a*b, c*d);
}

Tensor4::Tensor4(int a, int b, int c, int d, double val) : w(a), x(b), y(c), z(d)
{
  // Resize and assign the data matrix
  data.assign(a*b, c*d, val);
}

Tensor4::Tensor4(const Tensor4& other)
{
  w = other.w; x = other.x; y = other.y; z = other.z;
  data.resize(w*x, y*z);
  for (int i = 0; i < w; i++){
    for (int j = 0; j < x; j++){
      for (int k = 0; k < y; k++){
        for (int l = 0; l < z; l++){
          data(i*x+j, k*z+l) = other(i, j, k, l);
        }
      } 
    }
  }
}  

void Tensor4::resize(int a, int b, int c, int d)
{
  w = a; x = b; y = c; z = d;
  data.resize(a*b, c*d);
}

void Tensor4::assign(int a, int b, int c, int d, double val)
{
  w = a; x = b; y = c; z = d;
  data.assign(a*b, c*d, val);
}

void Tensor4::print() const 
{
  for (int i = 0; i < w; i++){
    for (int j = 0; j < x; j++){
      for (int k = 0; k < y; k++){
	for (int l = 0; l < z; l++){
	  std::cout << i << " " << j << " " << k << " " << l << "   " << data(i*x+j, k*z+l) << "\n";
	}
      }
    }
  }
  std::cout << "\n\n";
}

double& Tensor4::operator()(int i, int j, int k, int l)
{
  return data(i*x+j, k*z+l);
}

double Tensor4::operator()(int i, int j, int k, int l) const
{
  return data(i*x+j, k*z+l);
}

Tensor4& Tensor4::operator=(const Tensor4& other)
{
  resize(other.w, other.x, other.y, other.z);
  for (int i = 0; i < w; i++){
    for (int j = 0; j < x; j++){
      for (int k = 0; k < y; k++){
	for (int l = 0; l < z; l++){
	  data(i*x+j, k*z+l) = other(i, j, k, l);
	}
      }
    }
  }
  return *this;
}

Tensor4 Tensor4::operator+(const Tensor4& other) const
{
	Tensor4 retVal(w, x, y, z);
	for (int i = 0; i < w; i++){
		for (int j = 0; j < x; j++){
			for (int k = 0; k < y; k++){
				for (int l = 0; l < z; l++){
					retVal(i, j, k, l) = data(i*x+j, k*z+l) + other(i, j, k, l);
				}
			}
		}
	}
	return retVal;
}
