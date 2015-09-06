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
 *                  naints - a matrix of nuclear attraction integrals.
 *            data: sizes - a vector of the number of integrals needed for 
 *                          [1e cartesian, 2e cartesian, 1e spherical, 2e spherical]
 *                          assuming none can be neglected
 *            routines: 
 *                  getEstimates - returns a vector of estimates of the memory needed
 *                           to store each of the integral types in sizes
 *                  getVals(exponents, centres) - return the centre-of-charge 
 *                           coordinates, total exponents, reduced exponents,
 *                           and pre-exponential factors between two cgbfs
 *                  overlapKinetic(u, v, ucoords, vcoords) - calculates the overlap and kinetic 
 *                                         integrals between two primitives, u, v, given the 
 *                                         coordinates of their atomic centres.
 *                  makeContracted(coeffs1, coeffs2, ints) - contracts the given set of integrals
 *                           with the given sets of coefficients (1e- integrals)
 *                  makeContracted(coeffs1, coeffs2, coeffs2, coeffs4, ints) - same, but for
 *                           2e- integrals
 *                  makeSpherical(l1, m1, l2, m2, ints) - convert to spherical gaussians, 1e- ints
 *                  makeSpherical(l1, m1, l2, m2, l3, m3, l4, m4, ints) - same, 2e- ints
 *                  formOverlapKinetic() - forms the matrices sints, tints
 *                  multipoleComponent(a, b, acoord, bcoord, ccoord, powers) - calculates the multipole
 *                                     integral about c-coordinates to the power powers in each coordinate
 *                                     for the basis functions a, b
 *                  formNucAttract() - forms the matrix of nuclear attraction integrals, naints
 *                  makeERI() - returns a matrix of electron repulsion integrals
 *
 *   DATE          AUTHOR            CHANGES 
 *   =============================================================================
 *   02/09/15      Robert Shaw       Original code.
 *   03/09/15      Robert Shaw       Kinetic integrals merged with overlap.
 *   04/09/15      Robert Shaw       Multipole integrals added.
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
  Matrix naints;
  Vector sizes;
public:
  IntegralEngine(Molecule& m); //Constructor

  // Accessors
  Vector getEstimates() const;
  double getOverlap(int i, int j) const { return sints(i, j); }
  double getKinetic(int i, int j) const { return tints(i, j); }
  double getNucAttract(int i, int j) const { return naints(i, j); }
  Matrix makeERI() const;

  // Intrinsic routines
  Vector getVals(double a, double b, const Vector& A, const Vector& B) const;
  Vector overlapKinetic(const PBF& u, const PBF& v, const Vector& ucoords,
			const Vector& vcoords) const;
  double makeContracted(Vector& c1, Vector& c2, Vector& ints) const;
  double makeContracted(Vector& c1, Vector& c2, Vector& c3, 
			Vector& c4, Matrix& ints) const;
  Matrix makeSpherical(int l1, int l2, Matrix& ints) const;
  double makeSpherical(int l1, int m1, int l2, int m2, int l3, int m3,
		       int l4, int m4, Matrix& ints) const;
  void formOverlapKinetic();
  void formNucAttract();
  double multipole(BF& a,  BF& b, const Vector& acoords,
		   const Vector& bcoords, const Vector& ccoords, 
		   const Vector& powers) const;
  double multipole(PBF& u, PBF& v, const Vector& ucoords,
		   const Vector& vcoords, const Vector& ccoords,
		   const Vector& powers) const;
};

#endif  
    
