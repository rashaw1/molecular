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
#include "tensor4.hpp"
#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
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
  MAX = 10;
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
  Matrix S = integrals.getOverlap();
  Eigen::MatrixXd temp(S.nrows(), S.nrows());
  for (int i = 0; i < S.nrows(); i++){
    for (int j = 0; j < S.nrows(); j++){
      temp(i, j) = S(i, j);
    }
  }
  Matrix U(S.nrows(), S.nrows(), 0.0); Vector lambda(S.nrows(), 0.0);
  // Diagonalise the overlap matrix into lambda and U,
  // so that U(T)SU = lambda
  // We can now form S^(-1/2) - the orthogonalising matrix
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(temp);
  temp = es.eigenvectors();
  for (int i = 0; i < S.nrows(); i++){
    for (int j = 0; j < S.nrows(); j++){
      U(i, j) = temp(i, j);
    }
  }
  temp = es.eigenvalues().asDiagonal();
  for (int i = 0; i < S.nrows(); i++)
    lambda[i] = temp(i, i);
  
  orthog.assign(nbfs, nbfs, 0.0);
  for (int i = 0; i < nbfs; i++) {
    orthog(i, i) = 1.0/(std::sqrt(lambda(i)));
  }
  
  // S^-1/2  = U(lambda^-1/2)U(T)
  orthog = U * orthog * U.transpose();
}
  


// Transform the AO fock matrix to the MO basis 
void Fock::transform(bool first)
{
  if (first) { 
    // Form the core Fock matrix as (S^-1/2)(T)H(S^-1/2)
    fockm = (orthog.transpose()) * ( hcore * orthog);
  } else {
    // Form the orthogonalised fock matrix
    if (diis) DIIS();
   fockm = orthog.transpose() * (focka * orthog);
  }
}

// Diagonalise the MO fock matrix to get CP and eps
void Fock::diagonalise() 
{
  Eigen::MatrixXd temp(fockm.nrows(), fockm.nrows());
  for (int i = 0; i < fockm.nrows(); i++){
    for (int j = 0; j < fockm.nrows(); j++){
      temp(i, j) = fockm(i, j);
    }
  }

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(temp);
  temp = es.eigenvectors();
  CP.assign(fockm.nrows(), fockm.nrows(), 0.0); eps.assign(fockm.nrows(), 0.0);
  for(int i = 0; i < fockm.nrows(); i++){
    for (int j = 0; j < fockm.nrows(); j++){
      CP(i, j) = temp(i, j);
    }
  }
  
  temp = es.eigenvalues().asDiagonal();
  for (int i = 0; i < fockm.nrows(); i++)
    eps[i] = temp(i, i);

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
}

// Construct the density matrix from CP, 
// for nocc number of occupied orbitals
void Fock::makeDens(int nocc)
{
  // Form the density matrix
  dens.assign(nbfs, nbfs, 0.0);
  for (int u = 0; u < nbfs; u++){
    for (int v = 0; v < nbfs; v++){
      for (int t = 0; t < nocc; t++){
	dens(u, v) += CP(u, t)*CP(v,t);
      }
    }
  }
  dens = 2.0*dens;
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
  jints.assign(nbfs, nbfs, 0.0);
  kints.assign(nbfs, nbfs, 0.0);
  for (int u = 0; u < nbfs; u++){
    for (int v = 0; v < nbfs ; v++){
      for (int s = 0; s < nbfs; s++){
	for (int l = 0; l < nbfs; l++){
	  jints(u, v) += dens(s, l)*integrals.getERI(u, v, l, s);
	  kints(u, v) += dens(s, l)*integrals.getERI(u, s, l, v);
	}
      }
    }
  }
  jkints = jints - 0.5*kints;
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
  if (iter > MAX) {
    errs.erase(errs.begin()); // Remove first element
  }
  errs.push_back(e); // Push e onto the end of errs
  iter++;
}	
		

void Fock::makeFock()
{
  focka = hcore + jkints;
  if (diis) { // Archive for averaging
    if (iter > MAX) {
      focks.erase(focks.begin());
    }
    focks.push_back(focka);
  }		
}

void Fock::makeFock(Matrix& jbints)
{
  focka = hcore + 0.5*(jints + jbints - kints);
  if (diis) { // Archive for averaging
    if (iter > MAX) {
      focks.erase(focks.begin());
    }
    focks.push_back(focka);
  }
}

// Perform DIIS averaging
// Direct inversion of the iterative subspace
// Greatly improves convergence and numerical behaviour
// of the scf iterations.
void Fock::DIIS()
{
  if (iter > 1) {
    int lim = (iter < MAX ? iter : MAX);
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
    
    Eigen::MatrixXd temp(lim+1, lim+1);
    for (int i = 0; i < lim+1; i++){
      for (int j = 0; j < lim+1; j++){
	temp(i, j) = B(i, j);
      }
    }

    Eigen::JacobiSVD<Eigen::MatrixXd> svd(temp, Eigen::ComputeFullU | Eigen::ComputeFullV);
    temp = svd.singularValues().asDiagonal();
    for (int i = 0; i < lim+1; i++){
      if (fabs(temp(i, i)) > molecule.getLog().precision())
	temp(i, i) = 1.0/temp(i, i);
    }
    temp = svd.matrixV() * temp * svd.matrixU().transpose();
    
    for (int i = 0; i < lim+1; i++){
      for (int j = 0; j < lim+1; j++){
	B(i, j) = temp(i, j);
      }
    }

    w = B*w;
    
    // Average the fock matrices according to the weights
    focka.assign(nbfs, nbfs, 0.0);
    for (int i = 0; i < lim; i++) {
      focka = focka + w(i)*focks[i];   
    } 
  }
}

void Fock::simpleAverage(Matrix& D0, double weight)
{
  dens = weight*dens + (1.0-weight)*D0;
}
