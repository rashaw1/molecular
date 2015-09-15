/*
 *
 *   PURPOSE: To implement class Fock, a class containing the data and routines
 *            needed for Fock space based methods.
 *
 *   DATE          AUTHOR              CHANGES
 *  ===========================================================================
 *  14/09/15       Robert Shaw         Original code.
 * 
 */

#include "fock.hpp"
#include "error.hpp"
#include "factors.hpp"
#include "solvers.hpp"
#include "tensor4.hpp"
#include <iostream>
#include <cmath>
#include "logger.hpp"

// Constructor
Fock::Fock(IntegralEngine& ints, Molecule& m) : integrals(ints), molecule(m)
{
  // Make the core hamiltonian matrix
  formHCore();

  // Form the orthogonalising matrix and get initial density guess
  try {
    formOrthog();
  } catch (Error e) {
    molecule.getLog().error(e);
  }

  // Retrieve whether this calculation should be done direct, or whether
  // the twoints matrix has been formed, or if the 2e integrals need to be
  // read from file.
  direct = molecule.getLog().direct();
  diis = molecule.getLog().diis();
  iter = 0;
  MAX = 5;
  twoints = false;
  if (!direct){
    Vector ests = integrals.getEstimates();
    if (ests[3] < molecule.getLog().getMemory())
      twoints = true;
  }

  fromfile = false;
  if (!twoints && !direct)
    fromfile = true;

}

// Form the core hamiltonian matrix
void Fock::formHCore()
{
  // The core Hamiltonian matrix is defined to be
  // the sum of the kinetic and nuclear attraction
  // matrices
  hcore = integrals.getKinetic() + integrals.getNucAttract();
  nbfs = hcore.nrows();
}

void Fock::formOrthog()
{
  Matrix U; Vector lambda;
  // Diagonalise the overlap matrix into lambda and U,
  // so that U(T)SU = lambda
  if (symqr(integrals.getOverlap(), lambda, U, molecule.getLog().precision())) {
    // We can now form S^(-1/2) - the orthogonalising matrix
    
    orthog.assign(nbfs, nbfs, 0.0);
    for (int i = 0; i < nbfs; i++)
      lambda[i] = 1.0/(std::sqrt(lambda(i)));

    // S^-1/2  = U(lambda^-1/2)U(T)
    for (int i = 0; i < nbfs; i++){
      for (int j = 0; j < nbfs; j++){
	for (int k = 0; k < nbfs; k++){
	  orthog(i, j) += U(i, k)*lambda(k)*U(j, k);
	}
      }
    }
  } else { // Throw error if didn't work
    Error e("ORTHOG", "Failed to diagonalise overlap matrix.");
    throw(e);
  }
  
}

// Construct the density matrix - first => initial guess from hcore
void Fock::makeDens(int nocc, bool first)
{
  if (diis) DIIS(); // Do diis averaging if required
  
  if (first) { 
  	// Form the core Fock matrix as (S^-1/2)(T)H(S^-1/2)
  	fockm = (orthog.transpose()) * ( hcore * orthog);
  } else {
	// Form the orthogonalised fock matrix
  	fockm = orthog.transpose() * (fockm * orthog);
  }
  
  // Diagonalise the initial fock matrix
  if (symqr(fockm, eps, CP, molecule.getLog().precision())) {
    // Sort the eigenvalues and eigenvectors
    // using a selection sort
    int k;
    for (int i = 0; i < nbfs; i++){
      k=i;
      // Find the smallest element
      for (int j = i+1; j < nbfs; j++)
	if (eps(j) < eps(k)) { k=j; }
      
      // Swap rows of eps and columns of CP
      eps.swap(i, k);
      CP.swapCols(i, k);
    }

    // Form the initial SCF eigenvector matrix
    CP = orthog*CP;
    
    // Form the density matrix
    dens.assign(nbfs, nbfs, 0.0);
    for (int u = 0; u < nbfs; u++){
      for (int v = 0; v < nbfs; v++){
	for (int t = 0; t < nocc; t++){
	  dens(u, v) += CP(u, t)*CP(v, t);
	}
      }
    }
    dens = 2.0*dens;
    
  } else { // Throw error
    Error e("DENS", "Unable to diagonalise initial Fock matrix.");
  }
}

// Make the JK matrix, depending on how two electron integrals are stored/needed
void Fock::makeJK()
{
	if (twoints){
    formJK(); 
  } else if (direct) {
    formJKdirect();
  } else {
    try {
      formJKfile();
    } catch (Error e) {
      molecule.getLog().error(e);
    }
  }
}

// Form the 2J-K matrix, given that twoints is stored in memory
void Fock::formJK()
{
  jkints.assign(nbfs, nbfs, 0.0);
  for (int u = 0; u < nbfs; u++){
    for (int v = 0; v < nbfs; v++){
      for (int p = 0; p < nbfs; p++){
	for (int s = 0; s < nbfs; s++){
	  jkints(u, v) += dens(p, s)*(integrals.getERI(u, v, p, s) 
	  	- 0.5*integrals.getERI(u, s, p, v));
	}
      }
    }
  }
  		
}

// Form JK using integral direct methods
void Fock::formJKdirect()
{
}

// Form the JK matrix from two electron integrals stored on file
void Fock::formJKfile()
{
}

// Add an error vector
void Fock::addErr(Vector e)
{
	if (iter+1 > MAX) {
		errs.erase(errs.begin()); // Remove first element
	}
	errs.push_back(e); // Push e onto the end of errs
}	
		

void Fock::makeFock()
{
	fockm = hcore + jkints;
	if (diis) { // Archive for averaging
		if (iter+1 > MAX) {
			focks.erase(focks.begin());
		}
		focks.push_back(fockm);
	}		
	iter++;
}

// Perform DIIS averaging
// Direct inversion of the iterative subspace
// Greatly improves convergence and numerical behaviour
// of the scf iterations.
void Fock::DIIS()
{
	if (iter > 1) {
		int lim = (iter+1 < MAX ? iter+1 : MAX);
		Matrix B(lim+1, lim+1, -1.0); // Error norm matrix
		B(lim, lim) = 0.0;
		
		// The elements of B are <e_i | e_j >
		for (int i = 0; i < lim; i++){
			for (int j = i; j < lim; j++){
				B(i, j) = inner(errs[i], errs[j]);
				B(j, i) = B(i, j);
			}
		}
		
		// Solve the linear system of equations for the weights
		Vector w(lim+1, 0.0); w[lim] = -1.0;
		w = choleskysolve(B, w);
		
		// Average the fock matrices according to the weights
		fockm.assign(nbfs, nbfs, 0.0);
		for (int i = 0; i < lim; i++)
			fockm = fockm + w(i)*focks[i];
			
	}
}

