#include <iostream>
#include "atom.hpp"
#include "vector.hpp"
#include "ioutil.hpp"
#include "matrix.hpp"

int main (int argc, char* argv[])
{
  Vector coord(3); Vector c; 
  coord[0] = 1.50; coord[1] = 0.79; coord[2] = 3.43;
  Atom a(coord, 16, getAtomMass(16));
  std::cout << "Atom a has :\n";
  std::cout << "Mass " << a.getMass() << "\n";
  std::cout << "Charge " << a.getCharge() << "\n";
  c = a.getCoords();
  std::cout << "Coordinates "; c.print(); std::cout  << "\n";
  Atom b(coord, 12, getAtomMass(12));
  a.translate(-1.5, 0, 0.57);
  c = a.getCoords();
  std::cout << "Moved to "; c.print(); std::cout << "\n";
  std::cout << "b has mass, charge, coords: \n";
  std::cout << b.getMass() << " " << b.getCharge();
  c = b.getCoords();
  std::cout << " "; c.print(); std::cout << "\n";
  b = a;
  std::cout << "b has mass, charge, coords: \n";
  std::cout << b.getMass() << " " << b.getCharge();
  c = b.getCoords();
  std::cout << " " ; c.print(); std::cout << "\n";

  Matrix U(3, 3);
  U(0, 0) = 1.0; U(0, 1) = U(0, 2) = U(1, 0) = U(2, 0) = 0.0;
  U(1, 1) = U(2, 2) = 0.707; U(2, 1) = -0.707; U(1, 2) = 0.707;
  a.rotate(U);
  c = a.getCoords();
  std::cout << "a has coordinates "; c.print(); std::cout << "\n";
  c = b.getCoords();
  std::cout << "b has coordinates "; c.print(); std::cout << "\n";
  return 0;
}
