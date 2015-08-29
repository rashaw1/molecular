/* 
 *   
 *     PURPOSE: To define a class Logger, which will perform the majority of IO and data storage
 *              functions in the molecular suite. It is esentially acts as the communicator 
 *              between all the different units, 'logging' everything that happens. 
 *              For convenience, fundamental constants are stored here as well. 
 *
 *     class Logger:
 *              owns: outfile - an ofstream, for all primary logging functions to permanent
 *                              record.
 *                    infile - an ifstream for the input file
 *                    errstream - an ostream for logging all error messages - could be a file,
 *                                the console, any ostream
 *                    errs - an array of Error messages that have been thrown
 *              data: nerr, natoms - the number of errors accumulated, and the num. of atoms
 *                    clock - a boost::timer::cpu_timer for keeping track of time elapsed, and the time
 *                            that the log was instantiated at
 *              input storage: charge, multiplicity, atoms, basisset
 *              user defined constants: 
 *                    PRECISION - the numerical precision to be used throughout the program
 *                    MAXITER - the maximum number of iterations that will be performed
 *              fundamental constants:
 *                    M_PI - a definition of pi, in case one is not available for some reason
 *              conversion factors:
 *                    RTOCM, RTOMHZ - convert the rotation constants calculated by molecule into
 *                                    units of cm-1 or MHz, respectively.
 *                    TOKCAL, TOKJ - convert energies in Hartrees to kcal/mol or kj/mol
 *                    TOBOHR, TOANG - convert distances either from angstrom to bohr, or bohr 
 *                                    to angstrom, respectively
 *              routines:
 *                    INPUT:
 *                         getBasis(), charge(), multiplicity(), getAtom(charge), precision(),
 *                         maxiter() - all return the data/constant/object to which they refer
 *                         natoms() - returns the number of atoms
 *                    OUTPUT:
 *                         print() - overloaded function that will just print a general message
 *                                   or object in a generic way
 *                         result() - will print a message specifically formatted to stand out
 *                                    as a result
 *                         error() - will log an Error passed to it, in errs, and print it to
 *                                   the errstream.
 *                         getTime() - will call the timer functions, and return how much time
 *                                     passed since the last call.
 *                         finalise() - will close the output streams, flushing the buffers, and
 *                                      adding the customary final parts of the output.
 *
 *        DATE            AUTHOR              CHANGES    
 *        ===========================================================================
 *        27/08/15        Robert Shaw         Original code.
 *
 */

#ifndef LOGGERHEADERDEF
#define LOGGERHEADERDEF

// Exceptional use of preprocessor here for fundamental constants
// purely because they aren't always defined by every compiler
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Includes
#include <boost/timer.hpp>
#include <string>

// Declare forward dependencies
class ostream;
class ofstream;
class ifstream;
class Basis;
class Matrix;
class Vector;
class Atom;
class Error;

// Begin class declaration
class Logger
{
private:
  ifstream& infile;
  ofstream& outfile;
  ostream& errstream;
  Error* errs;
  Atom* atoms;
  int nerr, charge, multiplicity, natoms;
  boost::timer::cpu_timer clock;
  Basis basisset;
  // User defined constants
  double PRECISION;
  int MAXITER;
public:
  // Conversion factors
  static const double RTOCM;
  static const double RTOMHZ;
  static const double TOKCAL;
  static const double TOKJ;
  static const double TOBOHR;
  static const double TOANG;
  // Constructor/destructor
  Logger(ifstream& in, ofstream& out, ostream& e);
  ~Logger(); // Delete the various arrays
  // Accessors
  Basis getBasis() const;
  int charge() const { return charge; }
  int multiplicity() const { return multiplicity; }
  Atom& getAtom() const;
  double precision() const { return PRECISION; }
  int maxiter() const { return MAXITER; }
  int natoms() const { return natoms; }
  // Overloaded print functions
  void print(const std::string& msg) const; // Print a string message
  void print(const Vector& v) const; // Print out a formatted vector
  void print(const Matrix& m) const; // Matrix
  void print(const Basis& b) const; // Basis set - spec., no. of bfs, etc.
  void print(const Atom& a) const; // Atom - i.e coords, etc.
  void print(const Molecule& mol) const; // Molecule - geometry, charge, multiplicity, etc.
  void print(const BF& bf) const; // Basis function - coeffs and each pbf
  void print(const PBF& pbf) const; // Primitive gaussian - exponent, norm, ang. momenta
  // Specific logging formats
  void result(const std::string& msg) const;
  void error(const Error& e);
  double getTime(); 
  void finalise();
};
 
#endif

