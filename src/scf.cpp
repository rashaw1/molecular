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
SCF::SCF(Molecule& m, Fock& f) : molecule(m), focker(f), energy(0.0), last_energy(0.0), one_E(0.0), two_E(0.0), error(0.0), last_err(0.0)
{
	diis.init(8, molecule.getLog().diis());
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
double SCF::calcE(const Matrix& hcore, const Matrix& dens, const Matrix& fock) 
{
  one_E = 0.5*(dens*hcore).trace();
  two_E = 0.5*(dens*fock).trace();
  return (one_E+two_E);
}

// Calculate the error vector from the difference between
// the diagonalised MO fock matrix and the previous one
Vector SCF::calcErr(const Matrix& F, const Matrix& D, Matrix S, const Matrix& orthog)
{
  Matrix temp = (F*(D*S) - S*(D*F));
  temp = orthog.transpose() * temp * orthog;
  error = fnorm(temp);
  Vector err(temp.nrows()*temp.nrows(), 0.0);
  for (int u = 0; u < temp.nrows(); u++){
    for (int v = 0; v < temp.nrows(); v++){
      err[u*temp.nrows() + v] = temp(u, v);
    }
  } 
  return err;
}

Vector SCF::calcErr()
{
  Matrix& F = focker.getFockAO();
  Matrix& D = focker.getDens();
  Matrix S = focker.getIntegrals().getOverlap();
  Matrix& orthog = focker.getOrthog();

  Vector err = calcErr(F, D, S, orthog);
  return err;
}

// Determine the distance between the current and previous density
// matrices for convergence testing
bool SCF::testConvergence(double val)
{
  bool result = (val < molecule.getLog().converge() ? true : false);
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
    Matrix old_dens = focker.getDens();
    focker.makeJK();
    focker.makeFock();
	
	std::vector<Vector> errs;
	errs.push_back(calcErr());
   	Vector weights = diis.compute(errs);
	errs.clear();
	
    calcE();
    molecule.getLog().iteration(0, energy, 0.0, 0.0);
	focker.average(weights);
    focker.transform(false);
    int iter = 1;
    double delta, dd;
    
    while (!converged && iter < molecule.getLog().maxiter()) {
      // Recalculate
      focker.diagonalise();
      focker.makeDens(nel/2);
      dd = fnorm(focker.getDens() - old_dens);
      old_dens = focker.getDens();
      focker.makeJK();
      focker.makeFock();
      
	  errs.push_back(calcErr());
	  weights = diis.compute(errs);
	  errs.clear();
	  
      calcE();
      delta = fabs(energy-last_energy);
      molecule.getLog().iteration(iter, energy, delta, dd);
	  focker.average(weights);
      focker.transform(false);
      converged = testConvergence(dd);
      if ( delta > molecule.getLog().converge()/100.0 ) { converged = false; }
      iter++; 
    }
	focker.diagonalise();
	
    if (!converged) { 
      molecule.getLog().result("SCF failed to converge.");
    } else {
      molecule.getLog().print("\nOne electron energy (Hartree) = " + std::to_string(one_E));
      molecule.getLog().print("\nTwo electron energy (Hartree) = " + std::to_string(two_E));
      molecule.getLog().print("\n");
      molecule.getLog().orbitals(focker.getEps(), nel, false);
      molecule.getLog().result("RHF Energy = " + std::to_string(energy) + " Hartree");
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
  Matrix DA(focker.getFockMO().nrows(), focker.getFockMO().nrows(), 0.0); 
  Matrix DB(focker2.getFockMO().nrows(), focker2.getFockMO().nrows(), 0.0);
  //bool average = molecule.getLog().diis();
  double err1 = 0.0, err2 = 0.0, err1_last = 0.0, err2_last = 0.0;
  std::vector<Vector> errs;
  while (!converged && iter < molecule.getLog().maxiter()) {
    if (iter!= 1) {
      DA = focker.getDens(); DB = focker2.getDens();
    }
    focker.makeDens(nalpha); focker2.makeDens(nbeta);
    /*if (iter != 1 && average ) { 
      focker.simpleAverage(DA, 0.5); 
      focker2.simpleAverage(DB, 0.5);
    }*/
    focker.makeJK(); focker2.makeJK();
    focker.makeFock(focker2.getJ()); focker2.makeFock(focker.getJ());    

    errs.push_back(calcErr(focker.getFockAO(), focker.getDens(), focker.getIntegrals().getOverlap(), focker.getOrthog()));
    err1_last = err1;
    err1 = error;

    errs.push_back(calcErr(focker2.getFockAO(), focker2.getDens(), focker2.getIntegrals().getOverlap(), focker.getOrthog()));
    err2_last = err2;
    err2 = error;
    
	Vector weights = diis.compute(errs);
	errs.clear();
	focker.average(weights);
	focker2.average(weights);
	
    ea = calcE(focker.getHCore(), focker.getDens(), focker.getFockAO());
    eb = calcE(focker2.getHCore(), focker2.getDens(), focker2.getFockAO());

    focker.transform(); focker2.transform();
    focker.diagonalise(); focker2.diagonalise();
    
    last_energy = energy;
    energy = (ea + eb)/2.0 + molecule.getEnuc();
    delta = fabs(energy - last_energy);
    
    dist = fnorm((focker.getDens() + focker2.getDens()) - (DA + DB));
    //dist = 0.5*(err1+err2-err1_last-err2_last);
    
    molecule.getLog().iteration(iter, energy, delta, dist);

    if (delta < molecule.getLog().converge()/100.0 && dist < molecule.getLog().converge()) { converged = true; }
    iter++;
  }

  focker.diagonalise();
  focker2.diagonalise();
  if (converged) {
    // Construct the orbital energies
    molecule.getLog().print("\nALPHA ORBITALS");
    molecule.getLog().orbitals(focker.getEps(), nalpha, true);
    molecule.getLog().print("\nBETA ORBITALS");
    molecule.getLog().orbitals(focker2.getEps(), nbeta, true);
    molecule.getLog().result("UHF Energy = " + std::to_string(energy) + " Hartree");
  } else {
    molecule.getLog().result("UHF failed to converge");
  }
}
