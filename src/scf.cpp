/*
 *
 *   PURPOSE: To implement class SCF, which carries out HF self-consistent field calculations.
 *
 *   DATE        AUTHOR         CHANGES
 *   ===============================================================
 *   15/09/15    Robert Shaw    Original code.
 *
 */

#include "scf.hpp"
#include "logger.hpp"
#include "integrals.hpp"
#include <cmath>
#include "mvector.hpp"

// Constructor
SCF::SCF(Molecule& m, Fock& f) : molecule(m), focker(f), energy(0.0), last_energy(0.0), error(0.0), last_err(0.0)
{
}

// Routines

// Calculate the scf energy
void SCF::calcE()
{
  // Get the necessary matrices
  Matrix& dens = focker.getDens();
  Matrix& hcore = focker.getHCore();
  Matrix& fock = focker.getFockAO();
  
  // Calculate the energy
  last_energy = energy;
  energy = calcE(hcore, dens, fock) + molecule.getEnuc();
}

// Do the same but as an external function, where matrices are given as arguments
double SCF::calcE(const Matrix& hcore, const Matrix& dens, const Matrix& fock) const
{
  double e1 = (dens*hcore).trace();
  double e2 = (dens*fock).trace();
  return 0.5*(e1+e2);
}

// Calculate the error vector from the difference between
// the diagonalised MO fock matrix and the previous one
Vector SCF::calcErr()
{
  Matrix& F = focker.getFockAO();
  Matrix& D = focker.getDens();
  Matrix S = focker.getIntegrals().getOverlap();
  
  S = (F*(D*S) - S*(D*F));
  error = fnorm(S);
  Vector err(S.nrows()*S.nrows(), 0.0);
  for (int u = 0; u < S.nrows(); u++){
    for (int v = 0; v < S.nrows(); v++){
      err[u*S.nrows() + v] = S(u, v);
    }
  } 
  return err;
}

// Determine the distance between the current and previous density
// matrices for convergence testing
bool SCF::testConvergence()
{
  bool result = (fabs(error-last_err) < molecule.getLog().converge() ? true : false);
  last_err = error;
  return result;
}

// Do an rhf calculation
// Algorithm:
//    - Fock has formed orthog, hcore
//    - Make initial guess density
//    Until convergence:
//       - Form JK matrix
//       - Form fock matrix F = H + JK
//       - Calculate the electronic energy
//       - Make new density matrix - also gives orbitals
//       - Test for convergence
void SCF::rhf()
{
  // Check the multiplicity and number of electrons
  int nel = molecule.getNel();
  int mult = molecule.getMultiplicity();
  
  if ( (nel%2 != 0) ) {
    Error e1("RHF", "Molecule has an odd number of electrons.");
    molecule.getLog().error(e1);
  } else if (mult != 1) {
    Error e2("RHF", "Molecule is not a singlet state.");
    molecule.getLog().error(e2);
  } else { // All is fine
    molecule.getLog().title("RHF SCF Calculation");
    molecule.getLog().initIteration();
    bool converged = false;
    // Get initial guess
    focker.transform(true); // Get guess of fock from hcore
    focker.diagonalise();
    focker.makeDens(nel/2);
    focker.makeJK();
    focker.makeFock();
    focker.addErr(calcErr());
    calcE();
    molecule.getLog().iteration(0, energy, 0.0);
    focker.transform(false);
    int iter = 1;
    double delta;
    while (!converged && iter < molecule.getLog().maxiter()) {
      // Recalculate
      focker.diagonalise();
      focker.makeDens(nel/2);
      focker.makeJK();
      focker.makeFock();
      focker.addErr(calcErr());
      calcE();
      delta = fabs(energy-last_energy);
      molecule.getLog().iteration(iter, energy, delta);
      focker.transform(false);
      converged = testConvergence();
      if ( delta > molecule.getLog().converge()/100.0 ) { converged = false; }
      iter++;
    }
    
    if (!converged) { 
      molecule.getLog().result("SCF failed to converge.");
    } else {
      molecule.getLog().result("RHF Energy = " + std::to_string(energy) + " Hartree");
      focker.getEps().print();
   }
  }
}

// UHF
void SCF::uhf()
{
  // Make a second focker instance
  Fock focker2(focker.getIntegrals(), molecule);
  
  // Get number of alpha/beta electrons
  int nalpha = molecule.nalpha();
  int nbeta = molecule.nbeta();

  // Start logging
  molecule.getLog().title("UHF SCF Calculation");
  molecule.getLog().print("# alpha = " + std::to_string(nalpha));
  molecule.getLog().print("# beta = " + std::to_string(nbeta));
  molecule.getLog().print("\n");
  molecule.getLog().initIteration();
  bool converged = false;
  
  // Get initial guess                                                                                                                                                                      
  focker.transform(true); focker2.transform(true);
  focker.diagonalise(); focker2.diagonalise();
  int iter = 1;
  double delta, ea, eb, dist;
  Matrix DA; Matrix DB;
  bool average = molecule.getLog().diis();
  while (!converged && iter < molecule.getLog().maxiter()) {
    focker.makeDens(nalpha); focker2.makeDens(nbeta);
    if (iter != 1 && average ) { 
      focker.simpleAverage(DA, 0.5); 
      focker2.simpleAverage(DB, 0.5);
    }
    DA = focker.getDens(); DB = focker2.getDens();
    
    focker.makeJK(); focker2.makeJK();
    focker.makeFock(focker2.getJ()); focker2.makeFock(focker.getJ());
    ea = calcE(focker.getHCore(), focker.getDens(), focker.getFockAO());
    eb = calcE(focker2.getHCore(), focker2.getDens(), focker2.getFockAO());

    focker.transform(); focker2.transform();
    focker.diagonalise(); focker2.diagonalise();
    
    last_energy = energy;
    energy = (ea + eb)/2.0 + molecule.getEnuc();
    delta = fabs(energy - last_energy);
    molecule.getLog().iteration(iter, energy, delta);

    dist = fnorm((focker.getDens() + focker2.getDens()) - (DA + DB));
    if (delta < molecule.getLog().converge()/100.0 && dist < molecule.getLog().converge()) { converged = true; }
    iter++;
  }
  if (converged) {
    molecule.getLog().result("UHF Energy = " + std::to_string(energy));
  } else {
    molecule.getLog().result("UHF failed to converge");
  }
}
