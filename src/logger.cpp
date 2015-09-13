/*
 *
 *   PURPOSE: To implement the file logger.hpp, defining the class Logger, a communication and
 *            storage device used throughout the molecular suite.
 * 
 *   DATE          AUTHOR           CHANGES
 *   =========================================================================================
 *   29/08/15      Robert Shaw      Original code.
 *   30/08/15      Robert Shaw      Implemented the many, many print functions
 */

// Includes
#include <algorithm>
#include <iomanip>
#include "logger.hpp"
#include "molecule.hpp"
#include "bf.hpp"
#include "pbf.hpp"
#include "matrix.hpp"
#include "mvector.hpp"
#include "error.hpp"
#include "ioutil.hpp"
#include "filereader.hpp"

// Define static constants
const double Logger::RTOCM = 16.857630400054006;
const double Logger::RTOMHZ = 0.5052274767280815;
const double Logger::TOKCAL = 627.509469;
const double Logger::TOKJ = 2625.49962;
const double Logger::TOBOHR = 0.52917721092;
const double Logger::TOANG = 1.889726124565;

// Constructor
Logger::Logger(std::ifstream& in, std::ofstream& out, std::ostream& e) : infile(in), outfile(out), errstream(e)
{
  // Timer is started on initialisation of logger.
  last_time = 0;

  // Initialise error array - fixed maximum of 20 errors, because any more than that and the 
  // whole idea is pointless really! Dynamic memory seems an unecessary overhead for just storing
  // errors.
  errs = new Error[20];
  nerr = 0; // No errors as of yet (hopefully)!
  
  FileReader input(infile); // Declare a file reader object
  try {
    input.readParameters(); 
  } catch (Error e) {
    error(e);
  }
  // Do a whole heap of input file reading
  // The single variables are pretty easy
  charge = input.getCharge();
  multiplicity = input.getMultiplicity();
  PRECISION = input.getPrecision();
  MAXITER = input.getMaxIter();
  THRINT = input.getThrint();

  // Now we deal with the arrays
  natoms = input.getNAtoms(); // Get how many atoms there are
  if (natoms > 0){
    atoms = new Atom[natoms]; // Allocate memory for atoms array
 
    // Now loop to get all the atoms 
    std::string temp, token;

    try {
    input.readGeometry();
    } catch (Error e) {
      error(e);
    }

    std::string delimiter = ","; // Define the separation delimiter
    std::size_t position; int q; Vector coords(3); double m;
    Vector qs(natoms); // Hold the qs of all the atoms, for finding unique qs laterv
    for (int i = 0; i < natoms; i++){
      temp = input.getGeomLine(i); // Get a line from the geometry
      position = temp.find(delimiter); // Find first delimiter
      token = temp.substr(0, position); // Tokenise

      // Get rid of whitespace, if any
      token.erase(std::remove(token.begin(), token.end(), ' '), token.end());

      // Get the atom type and mass
      q = getAtomCharge(token);
      m = getAtomMass(q);
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
      if (tempqs(i) != qs(k-1)){
	qs[k] = tempqs[i];
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
      Error e("NOATOMS", "Nothing to see here.");
    error(e);
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

// Overloaded print functions

// Print out a basic string message
void Logger::print(const std::string& msg) const
{
  outfile << msg << "\n";
}

// Print out a vector to a given precision, either horizontally or vertically
void Logger::print(const Vector& v, int digits, bool vertical) const
{
  std::string ender = "";
  if(vertical){
    ender = "\n";
  }
  for (int i = 0; i < v.size(); i++){
    outfile << std::fixed << std::setprecision(digits) << std::setw(digits+5) << v[i] << ender;
  }
  outfile << "\n";
}

// Print out a matrix, row by row, to a given precision
void Logger::print(const Matrix& m, int digits) const
{
  // Print out row by row
  Vector temp(m.ncols());
  for (int i = 0; i < m.nrows(); i++){
    temp = m.rowAsVector(i);
    print (temp, digits, false);
  }
}

// Print out the basis set specification in the format:
// (example is C atom in 3-21G basis set)
// BASIS: 6-31G
// Total no. of cgbfs: 9
// Total no. of primitives: 9
// ===============
//  Specification
// ===============
// Atom      Shell     #CGBFs    #Prims
// ......................................
// C         s          3         6
//           p          6         3
// (if full = true, then continue with this)
// ===============
// Basis Functions
// ===============
// Atom      Shell     BF    Coeff        Exponent
// .................................................
// C         s         1     0.0617669    172.2560
//                           0.358794     25.91090
//                           0.700713     5.533350
// etc...
void Logger::print(Basis& b, bool full) const
{
  // Collect the data needed for printing
  int nbfs = b.getNBFs(); // Store number of cgbfs and prims
  int nprims = 0;

  Vector qs = b.getCharges();
  // Sort the qs and get rid of duplicates
  qs.sort();
  Vector qtemp(qs.size());
  qtemp[0] = qs(0);
  int k = 1;
  for (int i = 1; i < qs.size(); i++){
    if (qs(i) != qtemp[k-1]){
      qtemp[k] = qs(i);
      k++;
    }
  }
  qs = qtemp;
  qs.resizeCopy(k);

  // Now sum over all basis functions to get the number of prims
  Vector c(3); c[0] = 0.0; c[1] = 0.0; c[2] = 0.0;
  BF bftemp(c, 0, 0, 0, c, c);

  for (int i = 0; i < nbfs; i++){
    bftemp = b.getBF(i);
    nprims += bftemp.getNPrims();
  }

  // Start printing
  title("Basis Set");
  outfile << "BASIS: " << b.getName() << "\n";
  outfile << "Total no. of cgbfs: " << nbfs << "\n";
  outfile << "Total no. of prims: " << nprims << "\n";
  title("Specification");
  outfile << std::setw(8) << "Atom";
  outfile << std::setw(8) << "Shell";
  outfile << std::setw(8) << "#CGBFs";
  outfile << std::setw(8) << "#Prims\n";
  outfile << std::string(35, '.') << "\n";
  
  // loop over the atom types
  outfile << std::setprecision(2);
  Vector subshells; Vector sublnums; 
  for (int i = 0; i < k; i++){
    int nc = 0; int np = 0;
    outfile << std::setw(8) << getAtomName(qs(i));
    subshells = b.getShells(qs[i]);

    sublnums = b.getLnums(qs[i]);
    outfile << std::setw(8) << getShellName(sublnums[0]);
    outfile << std::setw(8) << subshells[0];

    for (int j = 0; j < subshells[0]; j++){
      np += b.getBF(qs[i], j).getNPrims();
    }

    nc += subshells[0];
    outfile << std::setw(8) << np << "\n";
    for (int j = 1; j < subshells.size(); j++){
      outfile << std::setw(8) << "";
      outfile << std::setw(8) << getShellName(sublnums[j]);
      outfile << std::setw(8) << subshells[j];
      np = 0;
      for (int l = 0; l < subshells[j]; l++){
	np += b.getBF(qs[i], nc + l).getNPrims();	
      }
      nc += subshells[j];
      outfile << std::setw(8) << np << "\n";
    }

  }
  outfile << std::setprecision(8);
  // Now print out basis functions if required
  if (full) {
    title("Basis Functions");
    outfile << std::setw(8) << "Atom";
    outfile << std::setw(8) << "Shell";
    outfile << std::setw(5) << "BF";
    outfile << std::setw(18) << "Coeff";
    outfile << std::setw(18) << "Exponent\n";
    outfile << std::string(58, '.') << "\n";
    // Loop over all the basis functions
    Vector subshell; Vector sublnums; 
    Vector coeffs; Vector exps;
    std::string filler = "";
    for (int i = 0; i < k; i++){
      subshell = b.getShells(qs(i));
      sublnums = b.getLnums(qs(i));

      // Loop over shells
      int sum = 0;
      for (int r = 0; r < subshell.size(); r++){ 
	// Loop over bfs
	for (int s = 0; s < subshell[r]; s++){
	  bftemp = b.getBF(qs(i), s+sum);
	  coeffs = bftemp.getCoeffs();
	  exps = bftemp.getExps();

	  // Loop over coeffs/exps
	  for (int t = 0; t < coeffs.size(); t++){
	    filler = ((r == 0 && s==0 && t==0) ? getAtomName(qs[i]) : "");
	    outfile << std::setw(8) << filler;
	    filler = ((s == 0 && t == 0) ? getShellName(sublnums[r]) : "");
	    outfile << std::setw(8) << filler;
	    filler = (t == 0 ? std::to_string(s+1) : "");
	    outfile << std::setw(5) << filler;
	    outfile << std::setw(18) << std::setprecision(8) << coeffs(t);
	    outfile << std::setw(18) << std::setprecision(8) << exps(t) << "\n";
	  }
	}
	sum += subshell[r];
      }
    }
  }
}      

// Print out the details of an atom, taking the form:
// Atom Type   Atomic Number    Atomic Mass(amu)  #CGBFS   x    y    z
void Logger::print(const Atom& a) const
{
  int q = a.getCharge();
  Vector c(3);
  c = a.getCoords();
  outfile << std::fixed << std::setprecision(6);
  outfile << std::setw(10) << getAtomName(q);
  outfile << std::setw(10) << q;
  outfile << std::setw(10) << a.getMass();
  outfile << std::setw(10) << a.getNbfs();
  outfile << std::setw(4) << "(" << std::setw(8) << c(0);
  outfile << ", " << std::setw(8) << c(1);
  outfile << ", " << std::setw(8) << c(2) << ")\n";
}

// Print out details of Molecule, taking the form:
// ========
// MOLECULE
// ========
// No. of e- = nel,  charge = q,  singlet/doublet/triplet etc.
// (if inertia = true then also prints the section:
// ............................
// Principal Moments of Inertia
// ............................
// Ia = ...,  Ib = ...,  Ic = ...
// Rotational type: ...
// ....................
// Rotational Constants
// ....................
// (the coordinate system is also transformed to inertial coords)
// =====
// ATOMS
// =====
// Atom     z     Mass (a.m.u)     Coordinates
// ...............................................................
// C        6     12.014           (0.0, 3.53, -1.24)
// etc...
void Logger::print(Molecule& mol, bool inertia) const
{
  title("Molecule");
  // Print out basic details
  outfile << "# electrons = " << mol.getNel() << ",  ";
  outfile << "charge = " << mol.getCharge() << ",  ";
 
  std::string temp;
  // Get the state type
  switch (mol.getMultiplicity()) {
  case 1: { temp = "Singlet"; break; }
  case 2: { temp = "Doublet"; break; }
  case 3: { temp = "Triplet"; break; }
  case 4: { temp = "Quartet"; break; }
  case 5: { temp = "Quintet"; break; }
  default: { temp = "Spintacular"; break; } // Are hextets even a thing?!
  }
  outfile << temp << "\n";
  
  // Print inertial details if needed, and rotate
  // into the inertial coordinate system
  if (inertia) {
    // Get it
    temp = mol.rType();
    Vector rconsts(3);
    rconsts = mol.rConsts(0); // cm-1
    Vector inert(3);
    inert = mol.getInertia(true);
    // Print it out
    outfile << std::string(30, '.') << "\n";
    outfile << "Principal Moments of Inertia\n";
    outfile << std::string(30, '.') << "\n";
    outfile << "Ia = " << std::setw(12) << inert(0);
    outfile << ",  Ib = " << std::setw(12) << inert(1);
    outfile << ",  Ic = " << std::setw(12) << inert(2) << "\n";
    outfile << "Rotational type: " << temp << "\n";
    outfile << std::string(29, '.') << "\n";
    outfile << "Rotational Constants / cm-1\n";
    outfile << std::string(29, '.') << "\n";
    outfile << "A = " << std::setw(12) << rconsts(0);
    outfile << ",  B = " << std::setw(12) << rconsts(1);
    outfile << ",  C = " << std::setw(12) << rconsts(2) << "\n";
    
  }
  
  // Finally, print out all the atoms
  title("Atoms");
  outfile << std::setw(10) << "Atom" << std::setw(10) << "z";
  outfile << std::setw(10) << "Mass" << std::setw(10) << "#CGBFs"; 
  outfile << std::setw(30) << "Coordinates" << "\n";
  outfile << std::string(70, '.') << "\n";
  for (int i = 0; i < mol.getNAtoms(); i++){
    print(mol.getAtom(i));
  }
}

// Print out a contracted gaussian basis function, in the form:
// #Prims = ..., s/px/py/pz/dxy/... type orbital
// Coefficients: ......
// Exponents: ......
void Logger::print(BF& bf) const
{
  // Get information
  Vector coef;
  coef = bf.getCoeffs();
  Vector exp;
  exp = bf.getExps();
  int lx = bf.getLx(); int ly = bf.getLy(); int lz = bf.getLz();
  
  // Work out the orbital type
  std::string temp;
  switch(lx + ly + lz){
  case 0: { temp = "s"; break; }
  case 1: { 
    if (lx == 1){
      temp = "px";
    } else if (ly == 1){
      temp = "py";
    } else {
      temp = "pz";
    }
    break;
  }
  case 2: { 
    if (lx == 1 && ly == 1){
      temp = "dxy";
    } else if (lx == 1 && lz == 1){
      temp = "dxz";
    } else if (ly == 1 && lz == 1){
      temp = "dyz";
    } else if (lz == 2) {
      temp = "dz2";
    } else {
      temp = "d(x2 - y2)";
    }
      break;
  }
  case 3: { temp = "f"; break; } 
  case 4: { temp = "g"; break; }
  case 5: { temp = "h"; break; }
  case 6: { temp = "Very angular"; break; }
  }

  // Print it out
  outfile << "#Primitives = " << coef.size();
  outfile << ",  " << temp << " type basis function\n";
  outfile << std::setw(15) << "Coefficients: ";
  for (int i = 0; i < coef.size(); i++){
    outfile << std::setw(12) << std::setprecision(9) << coef(i);
    if ((i % 4) == 0 && i != 0){ // Print 4 per line
      outfile << "\n" << std::setw(15) << "";
    }
  }
  outfile << "\n";
  outfile << std::setw(15) <<  "Exponents: ";
  for (int i = 0; i < exp.size(); i++){
    outfile << std::setw(12) << std::setprecision(9) << exp(i);
    if ( (i%4) == 0 && i!=0) {
      outfile << "\n" << std::setw(15) << "";
    }
  }
}

// Print out a primitive gaussian basis function in the form:
// Lx = ..., Ly = ..., Lz = ...
// norm = ..., exponent = ...
void Logger::print(const PBF& pbf) const
{
  outfile << "Lx = " << pbf.getLx();
  outfile << ",  Ly = " << pbf.getLy();
  outfile << ",  Lz = " << pbf.getLz() << "\n";
  outfile << "Norm = " << pbf.getNorm();
  outfile << ",  Exponent = " << pbf.getExponent() << "\n";
}
 

// Specialised printing routines

// Print a title like so:
//
// 
// =====
// TITLE
// =====
//
void Logger::title(const std::string& msg) const 
{
  std::string temp = msg;
  // Make upper case
  std::transform(temp.begin(), temp.end(), temp.begin(), ::toupper);
  
  // Print
  outfile << "\n\n";
  outfile << std::string(temp.length(), '=') << "\n";
  outfile << temp << "\n";
  outfile << std::string(temp.length(), '=') << "\n";
  outfile << "\n";
}

// Print a result like so:
//
// **********************
// result goes here
// **********************
// 
void Logger::result(const std::string& msg) const
{
  outfile << "\n";
  outfile << std::string(msg.length(), '*') << "\n";
  outfile << msg << "\n";
  outfile << std::string(msg.length(), '*') << "\n\n";
}

// Log an error
void Logger::error(Error& e)
{
  nerr++;
  errs[nerr%20] = e;
  errstream << "ERROR: " << e.getCode() << "\n";
  errstream << "Message: " << e.getMsg() << "\n";
  errTime();
}

// Timing functions

// Print time elapsed since last call
void Logger::localTime()
{
  boost::timer::nanosecond_type temp = timer.elapsed().wall;
  outfile << "Time taken: " <<  ((double)(temp - last_time))/(1e9)
	  << " seconds\n";
  last_time = temp;
}
 
// Print total time taken
void Logger::globalTime()
{
  outfile << "Total time: " << ((double)(timer.elapsed().wall))/(1e9) 
	  << " seconds\n";
}

// Print time at which error occured
void Logger::errTime()
{
  errstream << "Error after " << ((double)(timer.elapsed().wall))/(1e9) 
	    << " seconds\n";
}

// Return the localTime/globalTime, instead of printing
double Logger::getLocalTime() 
{
  boost::timer::nanosecond_type temp = timer.elapsed().wall;
  boost::timer::nanosecond_type rval = temp - last_time;
  last_time = temp;
  return ((double)(rval))/(1e9);
}

double Logger::getGlobalTime()
{
  return ((double)(timer.elapsed().wall))/(1e9);
}

// Finalise the output
// Currently just prints the time and the number of errors
// that occurred
void Logger::finalise()
{
  outfile << std::string(30, '-') << "\n";
  globalTime(); // Print time elapsed
  outfile << "Number of errors: " << nerr << "\n";
}
  
