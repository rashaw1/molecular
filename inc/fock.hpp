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
  Matrix jints;
  Matrix kints;
  Matrix orthog;
  Matrix fockm;
  Matrix focka;
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
  IntegralEngine& getIntegrals() { return integrals; }
  Molecule& getMolecule() { return molecule; }
  Matrix& getHCore() { return hcore; }
  Matrix& getFockAO() { return focka; }
  Matrix& getFockMO() { return fockm; }
  Matrix& getOrthog() { return orthog; }
  Matrix& getCP() { return CP; }
  Vector& getEps() { return eps; }
  Matrix& getJK() { return jkints; }
  Matrix& getJ() { return jints; } 
  Matrix& getK() { return kints; }
  Matrix& getDens() { return dens; }
  void setDIIS(bool d) { diis = d; } 
  void formHCore();
  void formOrthog();
  void addErr(Vector e);
  void transform(bool first = false);
  void diagonalise();
  void makeJK();
  void formJK();
  void formJK(Matrix& jbints);
  void formJKdirect();
  void formJKfile();
  void makeFock();
  void makeFock(Matrix& jbints);
  void makeDens(int nocc);
  void DIIS();
  void simpleAverage(Matrix& D0, double weight = 0.5);
};
#endif
