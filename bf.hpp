/*
 *
 *     PURPOSE: To define the class BF representing a
 *              contracted cartesian gaussian basis function.
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
 *                          getExps() will return a vector of the exponents
 *                          getNPrims() will return the size of coeffs, i.e
 *                          the number of primitives
 *                          setID(id) will set the id
 *              routines:
 *                       none
 *
 *     DATE        AUTHOR            CHANGES
 *     ====================================================================
 *     27/08/15    Robert Shaw       Original code.
 *
 */
 
#ifndef BFHEADERDEF
#define BFHEADERDEF

// Includes
#include "vector.hpp"
#include "pbf.hpp" 
 
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
  BF() : id(-1) { } // Default constructor
  BF(Vector& c, int l1, int l2, int l3, Vector& exps);
  BF(const BF& other); // Copy constructor
  ~BF(); // Delete the pbfs
  // Accessors
  PBF& getPBF(int i) { return pbfs[i]; }
  int getNPrims() const { return coeffs.size(); }
  Vector getCoeffs() { return coeffs; }
  Vector getExps() const; 
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

#endif
