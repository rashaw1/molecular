/*
 *
 *   PURPOSE:  A container for a 4-tensor of 6-tensors
 *
 */

#ifndef TEN4TEN6HEADERDEF
#define TEN4TEN6HEADERDEF

#include "tensor6.hpp"
#include <vector>

class Ten4Ten6
{
private:
  std::vector<Tensor6> arr;
  int w, x, y, z;
public:
  Ten4Ten6() : w(0), x(0), y(0), z(0) { }
  Ten4Ten6(int m, int n, int p, int q);
  ~Ten4Ten6();
  void clean();
  void set(int i, int j, int k, int l, Tensor6 tens);
  Tensor6& operator()(int i, int j, int k, int l);
  Tensor6 operator()(int i, int j, int k, int l) const;
};

#endif
