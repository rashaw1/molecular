/*
 *
 *    PURPOSE: Implements atom.hpp, defining class Atom
 * 
 *    DATE         AUTHOR            CHANGES
 *    =======================================================================
 *    27/08/15     Robert Shaw       Original code.
 *    05/09/15     Robert Shaw       Added function to count number of cgbfs
 *                                   in spherical basis.
 */

// Includes
#include "atom.hpp"
#include "pbf.hpp"
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

// Get the number of distinct primitives in a given shell
int Atom::getNShellPrims(int shell) const
{
  // Calculate the max. possible number of prims in this shell
  int nP = 0;
  for (int i = 0; i < nbfs; i++){
    nP += bfs[i].getNPrims();
  }

  // Work out starting and ending bf for this shell
  int start = 0;
  for (int i = 0; i < shell; i++)
    start += shells(i);
  int end = start + shells(shell);

  // Make a list of the ids of all prims on all bfs
  // in this shell
  Vector ids(nP);
  Vector temp;
  int k = 0;

  for (int i = start; i < end; i++){
    temp = bfs[i].getPrimList();
    for (int j = 0; j < temp.size(); j++){
      ids[k] = temp(j);
      k++;
    }
  }

  // Resize and sort
  ids.resizeCopy(k);
  ids.sort();

  // Count unique entries
  int count = 1;
  for (int i = 1; i < k; i++){
    if(ids(i)!=ids(i-1)){
      count++;
    }
  }

  return count;
}

// Get the ith BF in a given shell
BF& Atom::getShellBF(int shell, int i)
{
  // Work out the position in the bfs array
  // First find start position of this shell
  int start = 0;
  for (int i = 0; i < shell; i++){
    start += shells(i);
  }
  return bfs[i+start];
}

// Get the ith unique primitive bf of a given shell
PBF& Atom::getShellPrim(int shell, int i)
{
  bool found = false;
  int bf = -1; // Counter for cgbfs
  // Find the start of the shell
  for (int j = 0; j < shell; j++){
    bf += shells(j);
  }
  int pbf;
  Vector pList;
  // Loop until the prim with id i is found
  while(!found){
    bf++;
    pbf = 0;
    pList = bfs[bf].getPrimList();

    while(pbf < pList.size() && !found){
      found = (pList(pbf) == i ? true : false);
      if(!found){ pbf++; }
    }
  }
  return bfs[bf].getPBF(pbf);
}

int Atom::getNSpherical() const
{
  int scount=0, pcount=0, dcount=0, fcount=0, gcount=0; // Currently only cope with up to g-type bfs
  // Loop over all bfs
  for (int i = 0; i < nbfs; i++){
    switch(bfs[i].getLnum()){
    case 0: { // s type 
      scount++;
      break;
    }
    case 1: { // p type
      pcount++;
      break;
    }
    case 2: { // d type
      dcount++;
      break;
    }
    case 3: { // f type
      fcount++;
      break;
    }
    case 4: { // g type
      gcount++;
      break;
    }
    default: scount++; // Assume s type
    }
  }

  return scount+pcount+(5*dcount/6)+(7*fcount/10)+(9*gcount/15);
}

int Atom::getNSpherShellBF(int shell) const
{
	// Get angular momentum of this shell
	int L = lnums(shell);
	
	// Get no of cart. bfs in this shell
	int nc = shells(shell);
	
	// Calculate the corresponding number of spherical bfs
	int ns = 0;
	switch(L){
	case 2: { ns = 5*(nc/6); break; }
	case 3: { ns = 7*(nc/10); break; }
	case 4: { ns = 9*(nc/15); break; } 
	default: ns = nc;
	}
	
	return ns;
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

