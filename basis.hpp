/*
 *
 *     PURPOSE: To define the classes Basis, BF, and PBF, representing a
 *              basis set, a contracted gaussian basis function, and a
 *              primitive gaussian basis function, respectively.
 *
 *     class Basis:
 *              owns: bfs - a set of BFs
 *              data: name - the name of the basis set, needed for file io
 *                    charges - a list of atomic numbers corresponding to
 *                            each BF in bfs.
 *                    shells - a list representing which bfs are in which
 *                             shells - i.e. if the first 3 are in one shell,
 *                             the next two in another, etc, it would be
 *                             3, 2, ...
 *                    lnums - a list with the angular quantum number 
 *                            of each shell in shells                 
 *              routines:
 *                    findPosition(q) - find the first position in the bfs 
 *                                      array of an atom of atomic number q
 *                    findShellPosition(q) - find the first position in the 
 *                                           shells array
 *                    getBF(charge, i) - return the ith BF corresponding to
 *                                       an atom of atomic number charge   
 *                    getSize(charge) - return how many basis functions an atom
 *                                      of atomic number charge has
 *                    getShellSize(charge) - return how many shells it has
 *                    getShells(charge) - return a vector of the correct subset
 *                                        of shells
 *                    getLnums(charge) - return vector of subset of lnums
 *                  
 *     class BF: 
 *              owns: pbfs - a set of primitive basis functions
 *              data: coeffs - a list of contraction coefficients
 *                    norm - the normalisation constant of the function
 *                    lx, ly, lz - the angular momentum quantum number 
 *                                 components in the cartesian directions
 *                    id - a unique identifier for integral indexing
 *              accessors: all of the above have get routines
 *                          in addition, getCoeff(i) will return just the ith
 *                          coefficient
 *                          setID(id) will set the id
 *              routines:
 *                       none
 *                  
 *     class PBF:
 *              data: exponent - the gaussian exponent of the basis function
 *                    norm - the normalisation constant
 *                    lx, ly, lz - the angular momentum quantum numbers in the 
 *                                 cartesian directions
 *              accessors: all have get routines
 *              routines:
 *                    normalise - calculate the normalisation constant
 *
 *     DATE        AUTHOR            CHANGES
 *     ====================================================================
 *     27/08/15    Robert Shaw       Original code.
 *
 */

#ifndef BASISHEADERDEF
#define BASISHEADERDEF

// Includes
#include "vector.hpp"
#include <string>

// Begin class definitions

class Basis
{
private:
  BF* bfs;
  std::string name;
  Vector charges;
  Vector shells;
  Vector lnums;
public:
  // Constructors and destructor
  // Note - no copy constructor, as it doesn't really seem necessary
  // Need to specify the name of the basis, n, and a list of the 
  // distinct atoms that are needed (as a vector of atomic numbers)
  Basis(std::string n, Vector& atoms);
  ~Basis(); // Destructor
  // Accessors
  int findPosition(int q) const;
  int findShellPosition(int q) const;
  BF& getBF(int q, int i) const;
  int getSize(int q) const;
  int getShellSize(int q) const;
  Vector& getShells(int q) const;
  Vector& getLnums(int q) const;
  // Overloaded operators
  Basis& operator=(const Basis& other);
};
  
class BF
{
private:
  PBF* pbfs;
  Vector coeffs;
  double norm;
  int lx, ly, lz, id;
public:
  // Constructors and destructor
  // Need to specify a vector of contraction coefficients, c, the angular
  // momentum quantum numbers, lx = l1, ly = l2, lz = l3,
  // and a vector of primitive exponents, exps
  BF(Vector& c, int l1, int l2, int l3, Vector& exps);
  BF(const BF& other); // Copy constructor
  ~BF(); // Delete the pbfs
  // Accessors
  PBF& getPBF(int i) const { return pbfs[i]; }
  Vector& getCoeffs() const { return coeffs; }
  double getCoeff(int i) const { return coeffs[i]; }
  double getNorm() const { return norm; }
  int getLnum() const { return lx+ly+lz; }
  int getLx() const { return lx; }
  int getLy() const { return ly; }
  int getLz() const { return lz; }
  int getID() const { return id; }
  void setID(int i);
  // Overloaded operators
  BF& operator=(const BF& other);
};

class PBF
{
private:
  double exponent, norm;
  int lx, ly, lz;
public:
  // Constructors
  // Need to specify an exponent, e, and the angular momentum quantum numbers
  PBF(double e, int l1, int l2, int l3);
  PBF(const PBF& other); // Copy constructor
  // Accessors
  double getExponent() const { return exponent; }
  double getNorm() const { return norm; }
  int getLnum() const { return lx+ly+lz; }
  int getLx() const { return lx; }
  int getLy() const { return ly; }
  int getLz() const { return lz; }
  // Routines
  void normalise();
  // Overloaded operators
  PBF& operator=(const PBF& other);
};

#endif
