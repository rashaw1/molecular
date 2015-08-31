/*
 *    PURPOSE: To implement basis.hpp, defining class Basis, representing a basis set.
 *   
 *    DATE            AUTHOR              CHANGES
 *    =====================================================================
 *    27/08/15        Robert Shaw         Original code.
 *
 */

// Includes
#include "basis.hpp"
#include "bf.hpp"
#include "basisreader.hpp"
#include "ioutil.hpp"

// Constructors and destructors
Basis::Basis(std::string n, Vector& atoms)
{
  BasisReader input(n); // Make a basis reading object
  name = n;

  int natoms = atoms.size(); // Get how many different atoms there are
  // Now determine how many basis functions are needed
  int nbfs = 0;
  Vector qnbfs(natoms); // Store the nbfs for each atom
  for (int i = 0; i < natoms; i++){
    qnbfs[i] = input.readNbfs(atoms(i));
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
      bfs[k] = input.readBF(atoms(i), j);
      bfs[k].setID(k); // Index the basis functions
      charges[k] = atoms(i);
      k++;
    }
  }

  // Read in the shell data
  //shells = input.readShells(charges);
  //lnums = input.readLnums(charges);
}

Basis::~Basis()
{
  if (charges.size() > 0){
    delete[] bfs;
  }
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
BF& Basis::getBF(int q, int i)
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
Vector Basis::getShells(int q) const
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
Vector Basis::getLnums(int q) const
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


  


