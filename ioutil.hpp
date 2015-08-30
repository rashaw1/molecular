/*
 *
 *   PURPOSE: To declare a selection of functions and classes that are needed for
 *            the input/output interface of the molecular suite of programs.
 * 
 *   CONTAINS:
 *            class FileReader - a class that reads and stores information from
 *                               the user input file.
 *                    data: parameters (charge, multiplicity, basis, precision,
 *                          maxiter, natoms), geometry string
 *                          file positions: geomstart, geomend
 *                    routines: get for all parameters, getGeomLine
 *                            readParameters - reads in all parameters
 *                            readGeometry - reads in the geometry
 *
 *            class BasisReader - reads in basis specifications:
 *                    data: name, ifstream input
 *                    routines:
 *                       readNbfs(q) - read the number of basis functions an atom
 *                                     with charge q has in this basis set
 *                       readBF(q, i) - reads the ith basis function corresponding to
 *                                      atom q.
 *                       readShells(qs) - create a vector of the #shells that the atoms
 *                                        specified in qs would have
 *                       readLnums(qs) - same, but gives the ang. momentum. numbers
 *           
 *            General functions: getAtomMass, getAtomCharge, getAtomName
 *                               getShellName
 *                    
 *                                          
 *            
 *   DATE         AUTHOR         CHANGES 
 *   ==========================================================================
 *   30/08/15     Robert Shaw    Original code.
 *
 */

#ifndef IOUTILHEADERDEF
#define IOUTILHEADERDEF

// Includes
#include <string>
#include <fstream>

// Declare forward dependencies
class Vector;
class BF;

class FileReader
{
private:
  ifstream& input;
  int charge, multiplicity, maxiter, natoms;
  int geomstart, geomend;
  double precision;
  std::string basis;
  std::string* geometry;
public:
  FileReader(ifstream& in) : input(in) {} // Constructor
  void readParameters();
  void readGeometry();
  int getCharge() const { return charge; }
  int getMultiplicity() const { return multiplicity; }
  int getMaxIter() const { return maxiter; }
  int getNAtoms() const { return natoms; }
  int getBasis() const { return basis;}
  int getPrecision() const { return precision; }
  std::string getGeomLine() const;
};

class BasisReader
{
private:
  ifstream input;
  std::string name;
public:
  BasisReader(std::string n) : name(n) {} // Constructor
  int readNbfs(int q);
  BF& readBF(int q, int i);
  Vector readShells(const Vector& qs);
  Vector readLnums(const Vector& qs);
};

// General routines

// Get the atomic mass/name of an atom with atomic number q
double getAtomMass(int q);
std::string getAtomName(int q);

// Get the atomic number of an atom with name n
int getAtomCharge(const std::string& n);

// Get the text name for a shell of angular momentum l
std::string getShellName(int l);
  
#endif
