/* 
 * 
 *   PURPOSE: To declare a class IntegralEngine, which will store and calculate
 *            the molecular integrals necessary for ab initio calculations.
 * 
 *   class IntegralEngine:
 *            owns: molecule - a reference to a molecule on which the calculations
 *                             need to be carried out.
 *                  sints - a matrix of overlap integrals
 *                  tints - a matrix of kinetic integrals.
 *            data: sizes - a vector of the number of integrals needed for 
 *                          [1e cartesian, 2e cartesian, 1e spherical, 2e spherical]
 *                          assuming none can be neglected
 *            routines: 
 *                  getEstimates - returns a vector of estimates of the memory needed
 *                           to store each of the integral types in sizes
 *                  getVals(exponents, centres) - return the centre-of-charge 
 *                           coordinates, total exponents, reduced exponents,
 *                           and pre-exponential factors between two cgbfs
 *                  getN(l, m) - returns the spherical normalisation for a GTO with
 *                               angular and magnetic quantum numbers l, m
 *                  getC(l, m, t, u, v) - returns the Clebsch-Gordon coefficient
 *                  getS(l1, l2, Si0, AB) - calculates the overlap of two primitives via recursion
 *                                           on the set Si0.
 *                  makeContracted(coeffs1, coeffs2, ints) - contracts the given set of integrals
 *                           with the given sets of coefficients (1e- integrals)
 *                  makeContracted(coeffs1, coeffs2, coeffs2, coeffs4, ints) - same, but for
 *                           2e- integrals
 *                  makeSpherical(l1, m1, l2, m2, ints) - convert to spherical gaussians, 1e- ints
 *                  makeSpherical(l1, m1, l2, m2, l3, m3, l4, m4, ints) - same, 2e- ints
 *                  formOverlap(), formKinetic() - forms the matrices sints, tints
 *                  makeMultipole(pole) - returns a matrix of multipole integrals of order pole
 *                  makeNucAttract() - returns a matrix of nuclear attraction integrals
 *                  makeERI() - returns a matrix of electron repulsion integrals
 *
 *   DATE          AUTHOR            CHANGES 
 *   =============================================================================
 *   02/09/15      Robert Shaw       Original code.
 * 
 */ 

#ifndef INTEGRALSHEADERDEF
#define INTEGRALSHEADERDEF

// Includes
#include "matrix.hpp"
#include "vector.hpp"
#include "molecule.hpp"

// Declare forward dependencies

//Begin class declaration
class IntegralEngine
{
private:
  Molecule& molecule;
  Matrix sints;
  Matrix tints;
  Vector sizes;
public:
  IntegralEngine(Molecule& m); //Constructor

  // Accessors
  Vector getEstimates() const;
  double getOverlap(int i, int j) const { return sints(i, j); }
  double getKinetic(int i, int j) const { return tints(i, j); }
  Matrix makeMultipole(int pole) const;
  Matrix makeNucAttract() const;
  Matrix makeERI() const;

  // Intrinsic routines
  Vector getVals(double a, double b, const Vector& A, const Vector& B) const;
  double getN(int l, int m) const;
  double getC(int l, int m, int t, int u, double v) const;
  double getS(int l1, int l2, Vector& Si0, double AB) const;
  double makeContracted(Vector& c1, Vector& c2, Vector& ints) const;
  double makeContracted(Vector& c1, Vector& c2, Vector& c3, 
			Vector& c4, Matrix& ints) const;
  double makeSpherical(int l1, int m1, int l2, int m2, Matrix& ints) const;
  double makeSpherical(int l1, int m1, int l2, int m2, int l3, int m3,
		       int l4, int m4, Matrix& ints) const;
  void formOverlap();
  void formKinetic();
};

#endif  
    
