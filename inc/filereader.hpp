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
 
 // Forward dependencies
 
 class FileReader
{
private:
  std::ifstream& input;
  int charge, multiplicity, maxiter, natoms;
  int geomstart, geomend;
  double precision, thrint, memory;
  bool direct, twoprint;
  std::string basis, intfile;
  std::string* geometry;
  int findToken(std::string t); // Find the command being issued
public:
  FileReader(std::ifstream& in) : input(in), natoms(0) {} // Constructor
  ~FileReader(); // Destructor
  void readParameters();
  void readGeometry();
  int getCharge() const { return charge; }
  int getMultiplicity() const { return multiplicity; }
  int getMaxIter() const { return maxiter; }
  int getNAtoms() const { return natoms; }
  std::string getBasis() const { return basis;}
  std::string getIntFile() const { return intfile; }
  bool getDirect() const { return direct; }
  bool getTwoPrint() const { return twoprint; }
  double getMemory() const { return memory; }
  double getPrecision() const { return precision; }
  double getThrint() const { return thrint; }
  std::string& getGeomLine(int i) { return geometry[i]; }
};

#endif
