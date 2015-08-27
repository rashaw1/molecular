/*
 *
 *     Implements molecule.hpp, defining class Molecule
 *
 *     DATE         AUTHOR           CHANGES  
 *     ==============================================================
 *     26/08/15     Robert Shaw      Original code.
 *
 */

#include "atom.hpp"
#include "logger.hpp"
#include "molecule.hpp"
#include "matrix.hpp"
#include "vector.hpp"
#include "error.hpp"
#include "solvers.hpp"
#include <cmath>

// An initialisation function, for code reuse purposes
void Molecule::init()
{
  // Set the data variables
  charge = log.charge();
  multiplicity = log.multiplicity();

  // Get the basis set     
  bfset = log.getBasis();
  
  // Now declare the array of atoms
  int natoms = log.natoms();
  if (natoms != 0){

    atoms = new Atom[natoms];
    nel = 0; // Initialise to no electrons

    // Populate the array, and give the atoms
    // their basis functions. At the same time,
    // calculate the number of electrons
    for (int i = 0; i < natoms; i++) {
      atoms[i] = log.getAtom();
      atoms[i].setBasis(bfset);
      nel += atoms[i].getCharge();
    }
    
    // Account for overall charge
    nel -= charge;
    
    // Calculate the nuclear energy
    calcEnuc();

  } else { // Nothing to do
    log.error(Error("INIT", "There are no atoms!"));
    atoms = NULL;
  }
}

// Constructor
Molecule::Molecule(Logger& logger)
{
  // Set the log and then initialise
  log = logger;
  init();
}

// Copy constructor
Molecule::Molecule(const Molecule& other)
{
  log = other.log;
  init();
}


// Destructor
Molecule::~Molecule()
{
  if(log.natoms()!=0){
    delete [] atoms;
  }
}

// Routines

// rotate(Matrix U) rotates the coordinate system according
// to the transformation specified by the unitary matrix U
// U needs to be 3x3
void Molecule::rotate(Matrix& U)
{
  // Check if 3 x 3
  if(U.nrows() == 3 && U.ncols() == 3){
    for (int i = 0; i < log.natoms(); i++){
      atoms[i].rotate(U); // Do the rotation
    }
  } else {
    log.error(Error("ROTATE", "Unsuitable rotation matrix."));
  }
}

// translate(x, y, z) translates the coordinate system by
// adding the vector (x, y, z) to the coordinate vector
// of each atom.
void Molecule::translate(double x, double y, double z)
{
  for (int i = 0; i < log.natoms(); i++){
    atoms[i].translate(x, y, z);
  }
}

// nalpha and nbeta return the number of alpha and beta
// (spin up/spin down) electrons
int Molecule::nalpha() const
{
  // multiplicity = 2s+1; each alpha electron contributes
  // s = 1/2, therefore s = nalpha*(1/2), so
  return multiplicity - 1;
}

int Molecule::nbeta() const
{
  // The number of beta electrons is (nel - nalpha)/2
  return (nel - multiplicity + 1)/2;
}

// Calculate the nuclear energy, i.e. the sum of terms
// (Zi * Zj)/Rij for each distinct pair of nucleii, i, j
void Molecule::calcEnuc() const
{
  double zi;
  enuc = 0.0;

  // Outer loop over all atoms
  for (int i = 0; i < log.natoms(); i++) {
    zi = atoms[i].getCharge();
    // Inner loop over all atoms with index
    // greater than i, to avoid double counting
    for (int j = i+1; j < log.natoms(); j++){
      enuc += (zi*atoms[j].getCharge())/dist(i, j);
    }
  }
}
  
// Compute the centre of mass of the molecule
Vector Molecule::com() const 
{
  // Coordinate vector of centre of mass, set to all
  // zeroes for summations
  Vector c(3, 0.0); 
  
  // Loop over all atoms
  double sum = 0.0;
  double m;
  for (int i = 0; i < log.natoms(); i++){
    m = atoms[i].getMass();
    c = c + m*atoms[i].getCoords();
    sum += m;
  }

  c = (1.0/sum)*c;
  return c;
}

// Calculate the principal moments of inertia;
// optionally, shift to the inertial coordinate system
// First the elements of the inertia tensor are
// calculated, and then this is diagonalised to 
// give the eigenvalues - the principal moments - 
// which are returned as a 3-vector
Vector Molecule::getInertia(bool shift) 
{
  Vector v(3); // To return the moments in
  Matrix inertia(3, 3, 0.0); // To contain the elements

  // We must translate to centre of mass coordinates
  Vector c(3);
  c = -1.0*com();
  translate(c(0), c(1), c(2));

  // Loop over all atoms
  double m;
  for (int i = 0; i < log.natoms(); i++){
    m = atoms[i].getMass();
    c = atoms[i].getCoords();

    // The diagonal elements are the mass times the sum
    // of the squares of the other coordinates
    // e.g. I(0, 0) = m(y^2 + z^2)
    inertia(0, 0) += m*(c(1)*c(1) + c(2)*c(2));
    inertia(1, 1) += m*(c(0)*c(0) + c(2)*c(2));
    inertia(2, 2) += m*(c(0)*c(0) + c(1)*c(1));

    // The off diagonal elements are minus the 
    // mass times the two coordinates
    // e.g. I(1, 0) = -m*x*y = I(0 , 1)
    inertia(0, 1) -= m*c(0)*c(1);
    inertia(0, 2) -= m*c(0)*c(2);
    inertia(1, 2) -= m*c(1)*c(2);
  }
  // Finally fill out the symmetric elements
  inertia(1, 0) = inertia(0, 1); 
  inertia(2, 0) = inertia(0, 2); 
  inertia(2, 1) = inertia(1, 2);
  
  // Now we diagonalise, only need eigenvalues
  bool success;
  if(shift){
    Matrix U; 
    // U is the unitary transformation matrix
    // i.e. the eigenvectors of I
    success = symqr(I, v, U, log.precision());
    if (success) {
      rotate(U);
    }
  } else {
    // only calculate the eigenvalues
    success = symqr(I, v, log.precision());
  }

  if(!success){
    v[0] = v[1] = v[2] = 0.0;
    log.error(Error("INERTIA", "Inertia matrix failed to diagonalise."));
  }
  
  return v;
}


// Calculate the (Euclidean) distance between the atoms i and j
// dist2 gives the square of this distance
double Molecule::dist2(int i, int j) const
{
  Vector dvec(3);
  dvec = atoms[j].getCoords() - atoms[i].getCoords(); 
  return inner(dvec, dvec);
}

double Molecule::dist2(int i, int j) const
  return sqrt(dist2(i, j));
}

// Calculate bond, out of plane, and torsional angles
double Molecule::bondAngle(int i, int j, int k) const
{
  // Get the differences of the position vectors
  // and calculate the angle between them
  Vector rij(3); Vector rjk(3);
  rij = atoms[j].getCoords();
  rjk = rij;
  rij = rij - atoms[i].getCoords();
  rjk = atoms[k].getCoords() - rjk;
  return angle(rij, rjk); // Vector class friend function
}

double Molecule::oopAngle(int i, int j, int k, int l) const
{
  // The out of plane angle is given by the formula:
  // [(elj x elk).eli]/sin(anglejlk)
  // and sin(anglejlk) = ||ejl x elk||
  // therefore start by forming the unit vectors
  Vector elk(3); Vector elj(3); Vector eli(3);
  Vector l(3);
  l = atoms[l].getCoords();
  elk = atoms[k].getCoords() - l;
  elj = atoms[j].getCoords() - l;
  eli = atoms[i].getCoords() - l;
  
  // Normalise
  elk = (1.0/pnorm(elk))*elk;
  elj = (1.0/pnorm(elj))*elj;
  eli = (1.0/pnorm(eli))*eli;

  // Get the sine
  double sinangle = pnorm(cross(elj, elk));
  
  // Do the triple product and divide by the sine
  return triple(elj, elk, eli)/sinangle;
}

double Molecule::torsionAngle(int i, int j, int k, int l) const
{
  // The torsional angle is given by the formula:
  // [(eij x ejk).(ejk x ekl)]/[sin(ijk)*sin(jkl)]

  // Form the identity vectors
  Vector eij(3); Vector ejk(3); Vector ekl(3);
  Vector j(3); Vector k(3);
  j = atoms[j].getCoords(); k = atoms[k].getCoords();
  eij = j - atoms[i].getCoords();
  ejk = k - j;
  ekl = atoms[l].getCoords() - k;
  
  // Normalise
  eij = (1.0/pnorm(eij))*eij;
  ejk = (1.0/pnorm(ejk))*ejk;
  ekl = (1.0/pnorm(ekl))*ekl;

  // Get angles and cross products
  double sinijk = pnorm(cross(eij, ejk));
  double sinjkl = pnorm(cross(ejk, ekl));
  j = cross(eij, ejk);
  k = cross(ejk, ekl);
 
  // Return the inner product over the product of angles
  return inner(j, k)/(sinijk*sinjkl);
}

// Compute the rotational type and constants of the molecule
// using the getInertia method
std::string Molecule::rType() const
{
  std::string rstring = "asymmetric"; // Asymm. is default case
  // Diatomic case
  if(log.natoms() == 2){
    rstring = "diatomic";
  } else {
    // Calculate the principal moments of inertia
    Vector I(3);
    I = getInertia();
    
    // All moments are >= 0, so no need to use fabs
    // Linear if Ic = Ib > Ia = 0
    // Prolate symm. if Ic = Ib > Ia /= 0
    // Oblate symm. if Ic > Ib = Ia /= 0
    // Spherical if Ic = Ib = Ia /= 0
    double CUTOFF = log.precision();
    if ( (I(2)-I(1)) < CUTOFF && I(0) < CUTOFF ){
      rstring = "linear";
    }   else if ( (I(2) - I(1)) < CUTOFF && I(0) >= CUTOFF){
      rstring = "prolate";
    } else if ( (I(1) - I(0)) < CUTOFF && I(0) >= CUTOFF){
      rstring = "oblate";
    } else if ( (I(1) - I(0)) < CUTOFF && (I(2) - I(1)) < CUTOFF){
      rstring = "spherical";
    }
  }
  return rstring;
}

Vector Molecule::rConsts(int units)
{
  // First get the principal moments of inertia
  Vector I(3);
  I = getInertia();
  
  double K;
  // Choose units
  if (units == 0){ // reciprocal centimetres
    K = Logger::RTOCM; // ~27.99319901 cm-1 . A^2 . amu
  } else { // MHz
    K = Logger::RTOMHZ; // ~8.39214994 x 10^5 MHz . A^2 . amu
  }
  // Constant Bi = K/Ii 
  for (int i = 0; i < 3; i++){
    I[i] = K/I(i);
  }
  return I;
}
  


  
  
