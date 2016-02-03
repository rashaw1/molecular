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
 *                  nucAttract(u, v, ucoords, vcoords) - same but calculates nuclear attraction ints.
 *                  makeContracted(coeffs1, coeffs2, ints) - contracts the given set of integrals
 *                           with the given sets of coefficients (1e- integrals)
 *                  makeSpherical(ints, lnums) - transform a matrix of 1e cartesian integrals to a 
 *                                               spherical harmonic basis
 *                  formOverlapKinetic() - forms the matrices sints, tints
 *                  multipoleComponent(a, b, acoord, bcoord, ccoord, powers) - calculates the multipole
 *                                     integral about c-coordinates to the power powers in each coordinate
 *                                     for the basis functions a, b
 *                  formNucAttract() - forms the matrix of nuclear attraction integrals, naints
 *                  printERI(output) - prints a sorted list of ERIs to the ostream output
 *                  twoe(A, B, C, D, shellA, shellB, shellC, shellD) - calculate the (ab|cd) two electron
 *                                     contracted spherical integrals over a shell quartet on atoms A,B,C,D
 *                  twoe(u, v, w, x, ucoords, vcoords, wcoords, xcoords) - calculate the [u0|w0]
 *                                     2e- primitive cartesian integrals
 *
 *   DATE          AUTHOR            CHANGES 
 *   =============================================================================
 *   02/09/15      Robert Shaw       Original code.
 *   03/09/15      Robert Shaw       Kinetic integrals merged with overlap.
 *   04/09/15      Robert Shaw       Multipole integrals added.
 *   06/09/15      Robert Shaw       Nuclear attraction integrals.
 *   08/09/15      Robert Shaw       Auxiliary two elec. ints. added.
 *   09/09/15      Robert Shaw       Shell quartet 2e- ints added.
 *   10/09/15      Robert Shaw       Removed 2e makeContracted/Spherical.
 */ 

#ifndef INTEGRALSHEADERDEF
#define INTEGRALSHEADERDEF

// Includes
#include "matrix.hpp"
#include "mvector.hpp"
#include "molecule.hpp"
#include <iostream>
#include "tensor4.hpp"

// Declare forward dependencies
class Atom;
class Tensor6;

//Begin class declaration
class IntegralEngine
{
private:
  Molecule& molecule;
  Matrix sints;
  Matrix tints;
  Matrix naints;
  Matrix prescreen;
  Vector sizes;
  Tensor4 twoints;
public:
  IntegralEngine(Molecule& m); //Constructor

  // Accessors
  Vector getEstimates() const;
  double getOverlap(int i, int j) const { return sints(i, j); }
  Matrix getOverlap() const { return sints; }
  double getKinetic(int i, int j) const { return tints(i, j); }
  Matrix getKinetic() const { return tints; }
  double getNucAttract(int i, int j) const { return naints(i, j); }
  Matrix getNucAttract() const { return naints; }
  double getERI(int i, int j, int k, int l) const;
  Tensor4 getERI() const { return twoints; }

  // Intrinsic routines
  void printERI(std::ostream& output, int NSpher) const;
  void formERI(bool tofile);
  void diagERIThread(int start, int end, int NS, Vector &atoms, Vector &shells,
				Vector &bfs, Tensor4 &twints, Matrix &pscreen);
  void offDiagERIThread(int start, int end, int NS, Vector &atoms, Vector &shells,
			  	Vector &bfs, Tensor4 &twints);
  Vector getVals(double a, double b, const Vector& A, const Vector& B) const;
  Vector overlapKinetic(const PBF& u, const PBF& v, const Vector& ucoords,
			const Vector& vcoords) const;
  double nucAttract(const PBF& u, const PBF& v, const Vector& ucoords, 
		    const Vector& vcoords, const Vector& ccoords) const;
  double mmNucAttract(const PBF& u, const PBF& v, const Vector& ucoords,
  			const Vector& vcoords, const Vector& ccoords) const;
  Tensor4 makeE(int u, int v, double K, double p, double PA, double PB) const;
  Tensor4 twoe(Atom& A, Atom& B, Atom& C, Atom& D, int shellA, int shellB,
	      int shellC, int shellD) const;
  Tensor6 twoe(const PBF& u, const PBF& v, const PBF& w, const PBF& x, 
	      const Vector& ucoords, const Vector& vcoords, const Vector& wcoords,
	      const Vector& xcoords) const;
  double makeContracted(Vector& c1, Vector& c2, Vector& ints) const;
  Matrix makeSpherical(const Matrix& ints, const Vector& lnums) const;
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
    
