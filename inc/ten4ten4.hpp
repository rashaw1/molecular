/*
 *
 *   PURPOSE:  A container for a 4-tensor of 4-tensors
 *
 */

#ifndef TEN4TEN4HEADERDEF
#define TEN4TEN4HEADERDEF

#include "tensor4.hpp"
#include <vector>

class Ten4Ten4
{
private:
  std::vector<Tensor4> arr;
  int w, x, y, z;
public:
  Ten4Ten4() : w(0), x(0), y(0), z(0) { }
  Ten4Ten4(int m, int n, int p, int q);
  ~Ten4Ten4();
  void clean();
  void set(int i, int j, int k, int l, Tensor4 tens);
  Tensor4& operator()(int i, int j, int k, int l);
  Tensor4 operator()(int i, int j, int k, int l) const;
};

#endif
