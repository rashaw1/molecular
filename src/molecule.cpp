/*
 *
 *     Implements molecule.hpp, defining class Molecule
 *
 *     DATE         AUTHOR           CHANGES  
 *     ==============================================================
 *     26/08/15     Robert Shaw      Original code.
 *
 */

#include <iostream>
#include "logger.hpp"
#include "molecule.hpp"
#include "matrix.hpp"
#include "mvector.hpp"
#include "error.hpp"
#include "solvers.hpp"
#include <cmath>

// An initialisation function, for code reuse purposes
void Molecule::init()
{
  // Set the data variables
  charge = log.getCharge();
  multiplicity = log.getMultiplicity();

  // Get the basis set     
  bfset = log.getBasis();

  // Now declare the array of atoms
  natoms = log.getNatoms();
  if (natoms != 0){
    
    atoms = new Atom[natoms];
    nel = 0; // Initialise to no electrons

    // Populate the array, and give the atoms
    // their basis functions. At the same time,
    // calculate the number of electrons
    for (int i = 0; i < natoms; i++) {
      atoms[i] = log.getAtom(i);
      atoms[i].setBasis(bfset);
      nel += atoms[i].getCharge();
    }

    // Account for overall charge
    nel -= charge;
    
    // Calculate the nuclear energy
    calcEnuc();

  } else { // Nothing to do
    Error e("INIT", "There are no atoms!");
    log.error(e);
    atoms = NULL;
  }
}

// Constructor
Molecule::Molecule(Logger& logger, int q) : log(logger)
{
  // Set the log and then initialise
  init();
}

// Copy constructor
Molecule::Molecule(const Molecule& other) : log(other.log)
{
  init();
}


// Destructor
Molecule::~Molecule()
{
  if(log.getNatoms()!=0){
    delete [] atoms;
  }
}

// Routines

// rotate(Matrix U) rotates the coordinate system according
// to the transformation specified by the unitary matrix U
// U needs to be 3x3
void Molecule::rotate(const Matrix& U)
{
  // Check if 3 x 3
  if(U.nrows() == 3 && U.ncols() == 3){
    for (int i = 0; i < log.getNatoms(); i++){
      atoms[i].rotate(U); // Do the rotation
    }
  } else {
    Error e("ROTATE", "Unsuitable rotation matrix.");
    log.error(e);
  }
}

// translate(x, y, z) translates the coordinate system by
// adding the vector (x, y, z) to the coordinate vector
// of each atom.
void Molecule::translate(double x, double y, double z)
{
  for (int i = 0; i < log.getNatoms(); i++){
    atoms[i].translate(x, y, z);
  }
}

// nalpha and nbeta return the number of alpha and beta
// (spin up/spin down) electrons
int Molecule::nalpha() const
{
  // multiplicity = 2s+1; each unpaired alpha electron contributes
  // s = 1/2, therefore s = nunpairedalpha*(1/2), so
  
  return (nel-multiplicity+1)/2 + multiplicity -1;
}

int Molecule::nbeta() const
{
  // The number of beta electrons is (nel - nunpairedalpha)/2
  return (nel - multiplicity + 1)/2;
}

// Calculate the nuclear energy, i.e. the sum of terms
// (Zi * Zj)/Rij for each distinct pair of nucleii, i, j
void Molecule::calcEnuc()
{
  double zi;
  enuc = 0.0;

  // Outer loop over all atoms
  for (int i = 0; i < log.getNatoms(); i++) {
    zi = atoms[i].getCharge();
    // Inner loop over all atoms with index
    // greater than i, to avoid double counting
    for (int j = i+1; j < log.getNatoms(); j++){
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
  for (int i = 0; i < log.getNatoms(); i++){
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
  for (int i = 0; i < log.getNatoms(); i++){
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
    success = symqr(inertia, v, U, log.precision());
    if (success) {
      rotate(U);
    }
  } else {
    // only calculate the eigenvalues
    success = symqr(inertia, v, log.precision());
  }

  if(!success){
    v[0] = v[1] = v[2] = 0.0;
    Error e("INERTIA", "Inertia matrix failed to diagonalise.");
    log.error(e);
  }

  v.sort();
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

double Molecule::dist(int i, int j) const {
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
  Vector co(3);
  co = atoms[l].getCoords();
  elk = atoms[k].getCoords() - co;
  elj = atoms[j].getCoords() - co;
  eli = atoms[i].getCoords() - co;
  
  // Normalise
  elk = (1.0/pnorm(elk, 2))*elk;
  elj = (1.0/pnorm(elj, 2))*elj;
  eli = (1.0/pnorm(eli, 2))*eli;

  // Get the sine
  double sinangle = pnorm(cross(elj, elk), 2);
  
  // Do the triple product and divide by the sine
  return triple(elj, elk, eli)/sinangle;
}

double Molecule::torsionAngle(int i, int j, int k, int l) const
{
  // The torsional angle is given by the formula:
  // [(eij x ejk).(ejk x ekl)]/[sin(ijk)*sin(jkl)]

  // Form the identity vectors
  Vector eij(3); Vector ejk(3); Vector ekl(3);
  Vector jv(3); Vector kv(3);
  jv = atoms[j].getCoords(); kv = atoms[k].getCoords();
  eij = jv - atoms[i].getCoords();
  ejk = kv - jv;
  ekl = atoms[l].getCoords() - kv;
  
  // Normalise
  eij = (1.0/pnorm(eij, 2))*eij;
  ejk = (1.0/pnorm(ejk, 2))*ejk;
  ekl = (1.0/pnorm(ekl, 2))*ekl;

  // Get angles and cross products
  double sinijk = pnorm(cross(eij, ejk), 2);
  double sinjkl = pnorm(cross(ejk, ekl), 2);
  jv = cross(eij, ejk);
  kv = cross(ejk, ekl);
 
  // Return the inner product over the product of angles
  return inner(jv, kv)/(sinijk*sinjkl);
}

// Compute the rotational type and constants of the molecule
// using the getInertia method
std::string Molecule::rType()
{
  std::string rstring = "asymmetric"; // Asymm. is default case
  // Diatomic case
  if(log.getNatoms() == 2){
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
    double CUTOFF = 0.01;

    if ( (I(2)-I(1)) < CUTOFF && I(0) < CUTOFF ){
      rstring = "linear";
    }   else if ( (I(2) - I(1)) < CUTOFF && (I(1)-I(0)) > CUTOFF &&  I(0) >= CUTOFF){
      rstring = "prolate";
    } else if ( (I(1) - I(0)) < CUTOFF && (I(2) - I(1)) > CUTOFF &&  I(0) >= CUTOFF){
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
    K = Logger::RTOCM; 
  } else { // GHz
    K = Logger::RTOGHZ; // ~8.39214994 x 10^2 GHz . A^2 . amu
  }
  // Constant Bi = K/Ii 
  for (int i = 0; i < 3; i++){
    if (I(i) < log.precision()) { I[i] = 0.0; }  
    else { I[i] = K/I(i); }
  }
  return I;
}
  


  
  
