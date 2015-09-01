/*
 *
 *    PURPOSE: Implements atom.hpp, defining class Atom
 * 
 *    DATE         AUTHOR            CHANGES
 *    =======================================================================
 *    27/08/15     Robert Shaw       Original code.
 *
 */

// Includes
#include "atom.hpp"
#include "matrix.hpp"
#include <iostream>

//Constructors
Atom::Atom(const Vector& coords, int q, double m)
{
  x = coords(0);
  y = coords(1);
  z = coords(2);
  charge = q;
  mass = m;
  nbfs = 0;
  nshells = 0;
}

// Copy constructor
Atom::Atom(const Atom& other)
{
  x = other.x;
  y = other.y;
  z = other.z;
  mass = other.mass;
  charge = other.charge;
  nbfs = other.nbfs;
  nshells = other.nshells;
  
  // Copy in basis functions if they exist
  if (nbfs > 0){
    bfs = new BF[nbfs];
    for (int i = 0; i < nbfs; i++){
      bfs[i] = other.bfs[i];
    }
    shells = other.shells;
    lnums = other.lnums;
  }
}

// Destructor
Atom::~Atom()
{
  // Delete bfs if it was initialised
  if (nbfs > 0){
    delete[] bfs;
  }
}

// Coordinate accessor
Vector Atom::getCoords() const
{
  Vector v(3);
  v[0] = x; v[1] = y; v[2] = z;
  return v;
}

// Set the basis functions
void Atom::setBasis(Basis& bs)
{
  // Look up how many contracted basis functions an atom of this type has
  // in the given basis set, and how many shells
  nbfs = bs.getSize(charge);
  nshells = bs.getShellSize(charge);

  // Allocate memory, resize vectors
  bfs = new BF[nbfs];
  
  // Fill the array and the vectors
  for (int i = 0; i < nbfs; i++){
    bfs[i] = bs.getBF(charge, i); // Returns the ith BF for this atom type
  }
  shells = bs.getShells(charge);
  lnums = bs.getLnums(charge);
}

// Routines

// Translate(dx, dy, dz) translates the coordinates by [dx, dy, dz]
// while rotate(Matrix U) applies the unitary transformation U to the coords
void Atom::rotate(const Matrix& U)
{
  // Brute force is quicker than matrix-vector multiplication for
  // a 3x3 situation such as this
  x = U(0, 0)*x + U(0, 1)*y + U(0, 2)*z;
  y = U(1, 0)*x + U(1, 1)*y + U(1, 2)*z;
  z = U(2, 0)*x + U(2, 1)*y + U(2, 2)*z;
}

void Atom::translate(double dx, double dy, double dz)
{
  x += dx; y += dy; z += dz;
}

// Overloaded operators

// Overload the assignment operator, =
Atom& Atom::operator=(const Atom& other)
{
  // Check to see if this atom already has bfs
  // deallocate memory if it does

  // Assign attributes
  charge = other.charge;
  nbfs = other.nbfs;
  nshells = other.nshells;
  x = other.x; y = other.y; z = other.z;
  mass = other.mass;
  shells = other.shells;
  lnums = other.lnums;
 
  // Copy over basis functions
  if (nbfs > 0) {
    bfs = new BF[nbfs];
    for (int i = 0; i < nbfs; i++){
      bfs[i] = other.bfs[i];
    }
  }
  return *this;
}

