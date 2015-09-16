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
SCF::SCF(Molecule& m, Fock& f) : molecule(m), focker(f), energy(0.0), last_energy(0.0)
{
}

// Routines

// Calculate the scf energy
void SCF::calcE()
{
  // Get the necessary matrices
  Matrix& fock = focker.getFock();
  Matrix& dens = focker.getDens();
  Matrix& hcore = focker.getHCore();
  Matrix& jkints = focker.getJK();

  // Calculate the energy
  last_energy = energy;
  energy = 0.0;
  double e1 = 0.0;
  double e2 = 0.0;
  for (int u = 0; u < dens.nrows(); u++){
    for (int v = 0; v < dens.ncols(); v++){
      e1 += 2*dens(v, u)*hcore(u,v);
      e2 += dens(v, u)*jkints(u, v);
    }
  }  
  std::cout << "E1/E2: " << e1 << " " << e2 << "\n";
  energy = 0.5*(e1+e2) + molecule.getEnuc();
}

// Calculate the error vector from the difference between
// the diagonalised MO fock matrix and the previous one
Vector SCF::calcErr() const
{
  Matrix& CP = focker.getCP();

  Vector err(CP.nrows());
  for (int i =0; i< CP.nrows(); i++)
    err[i] = CP(i, i) - last_CP(i, i);
  
  return err;
}

// Determine the distance between the current and previous density
// matrices for convergence testing
bool SCF::testConvergence() const
{
  Matrix& dens = focker.getDens();
  
  // Calculate the Frobenius norm of the difference
  Matrix test = dens - last_dens;
  double dist = fnorm(test);
  bool result = (dist < molecule.getLog().converge() ? true : false);
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
    focker.makeDens(nel/2, true);
    int iter = 1;
    while (!converged && iter < molecule.getLog().maxiter()) {
      // Set the previous density and coeff matrices
      last_dens = focker.getDens();
      last_CP = focker.getCP();

      // Recalculate
      focker.makeJK();
      focker.makeFock();
      calcE();
      double delta = fabs(energy - last_energy);
      molecule.getLog().iteration(iter, energy, delta);
      focker.makeDens(nel/2);
      focker.addErr(calcErr());
      converged = testConvergence();
      if ( delta > molecule.getLog().converge()/10.0 ) { converged = false; }
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
}
