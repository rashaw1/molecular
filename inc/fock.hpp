/*
 *
 *   PURPOSE:  To define a class Fock, which contains the data and routines
 *             needed for methods based in Fock space. 
 *
 *       class Fock
 *              owns: hcore,  JK, orthog - the core hamiltonian, coulomb/ exchange, orthogonalising matrices
 *                          matrices.
 *                    dens, dens_1, dens_2 - the current and two previous density
 *                          matrices (previous needed for DIIS).
 *                    integrals - the integral engine
 *              data:
 * 
 *              routines:
 *   
 * 
 *     DATE        AUTHOR              CHANGES
 *    ==========================================================================
 *    14/09/15     Robert Shaw         Original code.
 * 
 */

#ifndef FOCKHEADERDEF
#define FOCKHEADERDEF

// Includes
#include "matrix.hpp"
#include "integrals.hpp"
#include "molecule.hpp"
#include "mvector.hpp"
#include <vector>

// Forward declarations

// Begin class definition
class Fock
{
private:
  Matrix hcore;
  Matrix jkints;
  Matrix orthog;
  Matrix fockm;
  Matrix CP;
  Vector eps;
  std::vector<Matrix> focks;
  std::vector<Vector> errs;
  Matrix dens;
  IntegralEngine& integrals;
  Molecule& molecule;
  bool direct, twoints, fromfile, diis;
  int nbfs, iter, MAX;
public:
  Fock(IntegralEngine& ints, Molecule& m);
  Matrix& getHCore() { return hcore; }
  Matrix& getFock() { return fockm; }
  Matrix& getCP() { return CP; }
  Vector& getEps() { return eps; }
  Matrix& getJK() { return jkints; }
  Matrix& getDens() { return dens; }
  void formHCore();
  void formOrthog();
  void addErr(Vector e);
  void guessDens();
  void makeJK();
  void formJK();
  void formJKdirect();
  void formJKfile();
  void makeFock();
  void makeDens(int nocc, bool first = false);
  void DIIS();


};
#endif
