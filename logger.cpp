/*
 *
 *   PURPOSE: To implement the file logger.hpp, defining the class Logger, a communication and
 *            storage device used throughout the molecular suite.
 * 
 *   DATE          AUTHOR           CHANGES
 *   =========================================================================================
 *   29/08/15      Robert Shaw      Original code.
 *
 */

// Includes
#include <iostream>
#include <fstream>
#include "basis.hpp"
#include "matrix.hpp"
#include "vector.hpp"
#include "atom.hpp"
#include "error.hpp"
#include "ioutil.hpp"

// Define static constants
const double Logger::RTOCM = 16.857630400054006;
const double Logger::RTOMHZ = 0.5052274767280815;
const double Logger::TOKCAL = 627.509469;
const double Logger::TOKJ = 2625.49962;
const double Logger::TOBOHR = 0.52917721092;
const double Logger::TOANG = 1.889726124565;

// Constructor
Logger::Logger(ifstream& in, ofstream& out, ostream& e) : infile(in), outfile(out), errstream(e)
{
  // Start the timer!
  

  // Initialise error array - fixed maximum of 20 errors, because any more than that and the 
  // whole idea is pointless really! Dynamic memory seems an unecessary overhead for just storing
  // errors.
  errs = new Error[20];
  nerr = 0; // No errors as of yet (hopefully)!
  
  FileReader input(infile); // Declare a file reader object
  input.readParameters(); 

  // Do a whole heap of input file reading
  // The single variables are pretty easy
  charge = input.getCharge();
  multiplicity = input.getMultiplicity();
  PRECISION = input.getPrecision();
  MAXITER = input.getMaxIter();

  // Now we deal with the arrays
  natoms = input.getNAtoms(); // Get how many atoms there are
  if (natoms > 0){
    atoms = new Atom[natoms]; // Allocate memory for atoms array
 
    // Now loop to get all the atoms 
    std::string temp, token;
    input.readGeometry();
    std::string delimiter = ","; // Define the separation delimiter
    int position, q; Vector coords(3); double m;
    Vector qs(natoms); // Hold the qs of all the atoms, for finding unique qs laterv
    for (int i = 0; i < natoms; i++){
      temp = input.getGeomLine(i); // Get a line from the geometry
      position = temp.find(delimiter); // Find first delimiter
      token = temp.substr(0, position); // Tokenise

      // Get rid of whitespace, if any
      token.erase(std::remove(token.begin(), token.end(), ' '), token.end());

      // Get the atom type and mass
      q = getAtomCharge(token);
      m = getAtomMass(token);
      qs[i] = q;

      // Move on to next token
      temp.erase(0, position+delimiter.length());
      position = temp.find(delimiter);
      token = temp.substr(0, position); 

      // This should now be the x coord
      token.erase(std::remove(token.begin(), token.end(), ' '), token.end());
      coords[0] = std::stod(token);  // Convert it to double

      // Repeat for y and z
      temp.erase(0, position+delimiter.length());
      position = temp.find(delimiter);
      token = temp.substr(0, position);
      token.erase(std::remove(token.begin(), token.end(), ' '), token.end());
      coords[1] = std::stod(token);
      temp.erase(0, position+delimiter.length());
      temp.erase(std::remove(temp.begin(), temp.end(), ' '), temp.end());
      coords[2] = std::stod(temp);

      // We can now initialise the atom and add to array
      Atom a(coords, q, m);
      atoms[i] = a;
    }
    
    // Next, the basis set
    // First find all the unique qs
    Vector tempqs(natoms);
    tempqs = qs.sorted();
    qs[0] = tempqs(0);
    int k = 1;
    for (int i = 1; i < natoms; i++){
      if (tempqs(i) != tempqs(i-1)){
	qs[k] = tempqs[k];
	k++;
      }
    }
    // k is now the number of unique qs, and all these unique qs are stored in qs
    // resize to get rid of extra weight
    qs.resizeCopy(k);
    // Find the basis name and initialise the basis set
    std::string bname = input.getBasis();
    Basis b(bname, qs);
    basisset = b;

  } else { // Our first error :( 
    error(Error("NOATOMS", "Nothing to see here."));
  }
}

// Destructor - get rid of the atoms and errors arrays
Logger::~Logger()
{
  if (natoms > 0){
    delete[] atoms;
  }
  delete[] errs;
}
