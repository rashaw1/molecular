/*
 *
 *     PURPOSE: To define the class Basis representing a
 *              basis set.
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
 *                    getNBFs() - returns the total number of basis functions
 *                    getName() - returns the name
 *                    getCharges() - returns the vector of charges
 *                    getBF(charge, i) - return the ith BF corresponding to
 *                                       an atom of atomic number charge   
 *                    getSize(charge) - return how many basis functions an atom
 *                                      of atomic number charge has
 *                    getShellSize(charge) - return how many shells it has
 *                    getShells(charge) - return a vector of the correct subset
 *                                        of shells
 *                    getLnums(charge) - return vector of subset of lnums
 *                 
 *
 *     DATE        AUTHOR            CHANGES
 *     ====================================================================
 *     27/08/15    Robert Shaw       Original code.
 *
 */

#ifndef BASISHEADERDEF
#define BASISHEADERDEF

// Includes
#include "mvector.hpp"
#include <string>

// Forward declarations
class BF;

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
  Basis() : name("Undefined") { } // Default constructor
  Basis(std::string n, Vector& atoms);
  ~Basis(); // Destructor
  // Accessors
  int getNBFs() const { return charges.size(); }
  std::string getName() const { return name; }
  Vector getCharges() const { return charges; }
  int findPosition(int q) const;
  int findShellPosition(int q) const;
  BF& getBF(int i);
  BF& getBF(int q, int i);
  int getSize(int q) const;
  int getShellSize(int q) const;
  Vector getShells(int q) const;
  Vector getLnums(int q) const;
  // Overloaded operators
  Basis& operator=(const Basis& other);
};

#endif
