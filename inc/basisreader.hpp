/*
 *
 *   PURPOSE: To declare class BasisReader, for reading a basis set in from files.
 * 
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
 *
 *   DATE         AUTHOR         CHANGES 
 *   ==========================================================================
 *   30/08/15     Robert Shaw    Original code.
 *
 */

#ifndef BASISREADERHEADERDEF
#define BASISREADERHEADERDEF

// Includes
#include <string>
#include <fstream>

// Declare forward dependencies
class Vector;
class BF;

class BasisReader
{
private:
  std::string name;
  std::ifstream input;
  void openFile(int q);
  void closeFile();
public:
  BasisReader(std::string n) : name(n) {} // Constructor
  int readNbfs(int q);
  BF readBF(int q, int i);
  Vector readShells(int q);
  Vector readLnums(int q);
};

#endif