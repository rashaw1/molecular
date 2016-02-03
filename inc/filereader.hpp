/*
 *
 *   PURPOSE: To declare a class FileReader, for reading the input file.
 * 
 *                    data: parameters (charge, multiplicity, basis, precision,
 *                          maxiter, natoms), geometry string
 *                          file positions: geomstart, geomend
 *                    routines: get for all parameters, getGeomLine(i) return ith line
 *                              of geometry.
 *                            readParameters - reads in all parameters
 *                            readGeometry - reads in the geometry
 *            
 *   DATE         AUTHOR         CHANGES 
 *   ==========================================================================
 *   30/08/15     Robert Shaw    Original code.
 *
 */
 #ifndef FILEREADERHEADERDEF
 #define FILEREADERHEADERDEF
 
// includes
#include <fstream>
#include <string>
#include <vector>
 
 // Forward dependencies
 
 class FileReader
{
private:
  std::ifstream& input;
  int charge, multiplicity, maxiter, natoms, nthreads;
  int geomstart, geomend;
  double precision, thrint, memory, converge;
  bool direct, twoprint, diis, bprint, angstrom;
  std::string basis, intfile;
  std::vector<std::string> geometry;
  std::vector<std::string> commands; 
  int findToken(std::string t); // Find the command being issued
public:
  FileReader(std::ifstream& in) : input(in), natoms(0) {} // Constructor
  ~FileReader(); // Destructor
  void readParameters();
  void readGeometry();
  int getCharge() const { return charge; }
  int getNThreads() const { return nthreads; }
  int getMultiplicity() const { return multiplicity; }
  int getMaxIter() const { return maxiter; }
  int getNAtoms() const { return natoms; }
  std::string getBasis() const { return basis;}
  std::string getIntFile() const { return intfile; }
  std::vector<std::string> getCmds() const { return commands; }
  bool getDirect() const { return direct; }
  bool getTwoPrint() const { return twoprint; }
  bool getDIIS() const { return diis; }
  bool getBPrint() const { return bprint; }
  bool getAngstrom() const { return angstrom; }
  double getMemory() const { return memory; }
  double getPrecision() const { return precision; }
  double getThrint() const { return thrint; }
  double getConverge() const { return converge; }
  std::string& getGeomLine(int i) { return geometry[i]; }
};

#endif
