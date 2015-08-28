/*
 *    PURPOSE: To implement basis.hpp, defining classes Basis, BF, and PBF
 *   
 *    DATE            AUTHOR              CHANGES
 *    =====================================================================
 *    27/08/15        Robert Shaw         Original code.
 *
 */

// Includes
#include "ioutil.hpp"
#include <cmath>
#include "mathutil.hpp"

// Constructors and destructors

Basis::Basis(std::string n, Vector& atoms)
{
  name = n;

  int natoms = atoms.size(); // Get how many different atoms there are
  // Now determine how many basis functions are needed
  int nbfs = 0;
  Vector qnbfs(natoms); // Store the nbfs for each atom
  for (int i = 0; i < natoms; i++){
    qnbfs[i] = readNbfs(name, atoms(i));
    nbfs += qnbfs(i);
  }
  
  // Allocate memory
  if(nbfs > 0){
    bfs = new BF[nbfs];
  } else {
    bfs = NULL;
  }
  charges.resize(nbfs);  

  // Fill the array, and charges
  int k = 0; // Index counter for bfs array
  for (int i = 0; i < natoms; i++){
    for (int j = 0; j < qnbfs(i); j++){
      bfs[k] = readBF(atoms(i), j);
      bfs[k].setID(k); // Index the basis functions
      charges[k] = atoms(i);
      k++;
    }
  }

  // Read in the shell data
  shells = readShells(name, charges);
  lnums = readLnums(name, charges);
}

Basis::~Basis()
{
  if (charges.size() > 0){
    delete[] bfs;
  }
}

BF::BF(Vector& c, int l1, int l2, int l3, Vector& exps)
{
  coeffs = c;
  lx = l1; ly = l2; lz = l3;
  
  // Construct the set of pbfs
  int npbfs = exps.size();
  if (npbfs > 0){
    pbfs = new PBF[npbfs];
    for (int i = 0; i < npbfs; i++){
      PBF temp(exps(i), l1, l2, l3);
      pbfs[i] = temp;
    }
  } else {
    pbfs = NULL;
  }
  normalise();
}

// Copy constructor
BF::BF(const BF& other)
{
  coeffs = other.coeffs;
  norm = other.norm;
  lx = other.lx; ly = other.ly; lz = other.lz;
  
  int npbfs = coeffs.size();
  if(npbfs > 0){
    pbfs = new PBF[npbfs];
    for (int i = 0; i < npbfs; i++){
      pbfs[i] = other.pbfs[i];
    }
  } else {
    pbfs = NULL;
  }
}

// Destructor
BF::~BF()
{
  if (coeffs.size() > 0){
    delete[] pbfs;
  }
}

PBF::PBF(double e, int l1, int l2, int l3) : exponent(e), lx(l1), ly(l2), lz(l3)
{
  normalise();
}

// Copy constructor
PBF::PBF(const PBF& other)
{
  exponent = other.exponent;
  lx = other.lx; ly = other.ly; lz = other.lz;
  norm = other.norm;
}


// Accessors

// Find the first position in bfs at which atom of atomic number q occurs
int Basis::findPosition(int q) const
{
  // Loop to find q
  int position = 0;
  bool found = false;
  while (!found && position < charges.size()) {
    found = (charges(position) == q ? true : false);
    position++;
  }
  return position;
}

// Same as above but for shells instead of bfs
int Basis::findShellPosition(int q) const
{
  int position = findPosition(q);
 
  int sum = 0;
  int i = 0;
  // Locate the entries in shells corresponding to q
  while(sum < position && i < shells.size()){
    sum+= shells(i);
    i++;
  }
  // i now points to the first element in shells corresponding to q
  return i;
}


// Return the ith basis function corresponding to an atom
// with atomic number q
BF& Basis::getBF(int q, int i) const
{
  int position = findPosition(q);
  // No bounds checking
  return bfs[position+i];
}

// Return the number of basis functions that an atom of atomic number q has
int Basis::getSize(int q) const
{
  int position = findPosition(q);
  // Now find the last position of q
  int size = 1;
  bool found = false;
  while(!found && (size+position) < charges.size()){
    found = (charges(position+size) != q ? true : false);
    size++;
  }
  return size;
}

// Find the number of shells that atom of atomic number q has
int Basis::getShellSize(int q) const
{
  int nbfs = getSize(q);
  int i = findShellPosition(q);
  
  // Starting 
  int sum = 0;
  int size = 0;
  while(sum < nbfs && i < shells.size()){
      sum += shells(i);
      i++; size++;
  }
  return size;
}

// Return the subset of shells corresponding to q
Vector& Basis::getShells(int q) const
{
  int position = findShellPosition(q);
  int size = getShellSize(q);
  Vector s(size);
  // Fill in the values
  for (int i = 0; i < size; i++){
    s[i] = shells[position+i];
  }
  return s;
}

// Find th subset of lnums corresponding to atom with atomic number q
Vector& Basis::getLnums(int q) const
{
  int position = findShellPosition(q);
  int size = getShellSize(q);
  Vector l(size);
  // Fill in the values
  for (int i = 0; i < size; i++){
    l[i] = lnums[position+i];
  }
  return l;
}

void PBF::setID(int i)
{
  id = i;
}

// Routines

// Calculate the normalisation constant for a primitive cartesian gaussian
void PBF::normalise()
{
  // The formula can be found in Taketa, Huzinaga, and O-ohata, Journal of
  // the Physical Society of Japan, Vol. 21, No. 11, Nov 1966:
  // Gaussian-Expansion Methods for Molecular Integrals
  norm = pow(2, 2*(lx+ly+lz) + 1.5)*pow(exponent, lx+ly+lz+1.5);
  // Calculate double factorials
  norm = norm / ( (double) (fact2(2*lx-1) * fact2(2*ly-1) * fact2(2*lz-1)) );
  norm = norm / pow(M_PI, 1.5);
  norm = sqrt(norm);
}
  
// Overloaded operators
Basis& Basis::operator=(const Basis& other)
{
  // If basis functions already exist, deallocate memory
  if(charges.size() > 0){
    delete[] bfs;
  }

  // Assign attributes
  name = other.name;
  charges = other.charges;
  shells = other.shells;
  lnums = other.lnums;
  
  // Copy across bfs
  int nbfs = charges.size();
  if (nbfs > 0){
    bfs = new BF[nbfs];
    for (int i = 0; i < nbfs; i++){
      bfs[i] = other.bfs[i];
    }
  }
  return *this;
}

BF& BF::operator=(const BF& other)
{
  // If PBFs already exist, deallocate memory
  if(coeffs.size() > 0){
    delete[] pbfs;
  }
  
  // Assign attributes
  coeffs = other.coeffs;
  norm = other.norm;
  lx = other.lx; ly = other.ly; lz = other.lz;

  // Copy across PBFs
  int npbfs = coeffs.size();
  if (npbfs > 0){
    pbfs = new PBF[npbfs];
    for (int i = 0; i < npbfs; i++){
      pbfs[i] = other.pbfs[i];
    }
  }
  
  return *this;
}

PBF& PBF::operator=(const PBF& other)
{
  // Assign attributes
  exponent = other.exponent;
  norm = other.norm;
  lx = other.lx; ly = other.ly; lz = other.lz;
  return *this;
}
  

  


