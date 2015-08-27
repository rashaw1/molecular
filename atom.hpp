/*
 *
 *       PURPOSE: To define a class Atom, which will contain the basic info
 *                and methods needed by the Molecule class, for each atom in
 *                a molecular calculation. 
 *    
 *       class Atom:
 *                owns: bfs - a set of contracted basis functions
 *                data: x, y, z - the cartesian coordinates of the atom
 *                      nbfs - the number of basis functions
 *                      nshells - the number of shells of basis functions
 *                      shells - a list corresponding to which bfs are in 
 *                               which shell
 *                      lnums - a list of what angular momentum type each
 *                              shell in shells has
 *                      charge - the atomic number (i.e. the charge in a.u.)
 *                      mass - the atomic mass (i.e. the mass in a.m.u)
 *                accessors: all of the above have get... routines
 *                           note that getCoords() returns an [x, y, z] vector
 *                           bfs has a setBasis(Basis) routine
 *                routines:
 *                      rotate(Matrix U) - rotate coords according to unitary
 *                                         transformation matrix, U
 *                      translate(x, y, z) - translate coords [+x, +y, +z]
 *                      
 *       DATE           AUTHOR             CHANGES 
 *       =====================================================================
 *       27/08/15       Robert Shaw        Original code.
 *
 */

#ifndef ATOMHEADERDEF
#define ATOMHEADERDEF

// Includes
#include "vector.hpp"

// Declare forward dependencies
class Matrix;
class BF;

// Begin class definition
class Atom
{
private:
  BF* bfs;
  int charge, nbfs, nshells;
  Vector shells, lnums;
  double x, y, z, mass;
public:
  // Constructors
  Atom(const Vector& coords, int q, double m); // q = charge, m = mass
  Atom(const Atom& other); // Copy constructor
  ~Atom(); // Destructor - gets rid of array bfs
  // Accessors
  int getCharge() const { return charge; }
  double getMass() const { return mass; }
  int getNbfs() const { return nbfs; }
  int getNshells() const { return nshells; }
  Vector getCoords() const;
  Vector& getShells const { return shells; }
  Vector& getLnums const { return lnums; }
  BF& getBF(int i) const { return bfs[i]; } // Return bf i - no bounds check
  void setBasis(const Basis& bs); // Set the basis functions using basis set bs
  // Routines
  void rotate(const Matrix& U); 
  void translate(double dx, double dy, double dz);
  // Overloaded operators
  Atom& operator=(const Atom& other);
};  
  
#endif
