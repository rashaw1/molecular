/*
 *      PURPOSE: To define a class Molecule, which will be the work-horse of
 *               the entire suite of programs. It is defined as follows:
 *        
 *       class Molecule:
 *           owns:
 *                atoms: an array of the atoms contained in the molecule
 *                logger: an object for logging all and any information,
 *                        errors, global data, etc.
 *           data:
 *                charge: the overall charge of the molecule
 *                nel: the total number of electrons in the molecule
 *                multiplicity: the spin multiplicity of the molecule
 *                enuc: the nuclear energy of the molecule
 *           accessors: - all data has a get... accessor, e.g. getNel
 *           routines:
 *                rotate(Matrix U): rotate the coordinate system, according
 *                                  to the transformation in the unitary
 *                                  matrix, U
 *                translate(x, y, z): translate the coordinate system 
 *                nalpha: returns the number of alpha electrons
 *                nbeta: returns the number of beta electrons
 *                calcEnuc: calculates the nuclear energy, and stores in enuc
 *                com: calculates the centre of mass of the molecule
 *                getInertia(shift): calculates the principal moments 
 *                                  of inertia, optionally shifting to the
 *                                  inertial frame of reference if shift=true
 *                dist(i, j): calculates the distance between atoms i and j
 *                dist2(i, j): same as above, but squared
 *                bondAngle(i, j, k): calc. bond angle between atoms i,j,k
 *                oopAngle(i, j, k, l): calc out of plane angle
 *                torsionAngle(i, j, k, l): torsional angle
 *                rType: returns the molecules rotational type, giving one of
 *                       ['diatomic', 'linear', 'asymmetric', 'oblate',
 *                        'prolate', 'spherical']
 *                rConsts(units): returns rotational constants, A,B,C
 *                                units = 0 -> cm-1; 1 -> MHz
 *
 *     DATE            AUTHOR                CHANGES
 *     =======================================================================
 *     26/08/15        Robert Shaw           Original code.
 * 
 */

#ifndef MOLECULEHEADERDEF
#define MOLECULEHEADERDEF

// Includes
#include "basis.hpp"
#include <string>

// Declare forward dependcies
class Atom;
class Logger;
class Matrix;
class Vector;

// Begin class definition
class Molecule
{
private:
  Basis bfset;
  Logger& log;
  Atom* atoms;
  int charge, nel, multiplicity;
  double enuc;
public:
  // Constructors and destructor
  void init(); // An initialisation function
  Molecule(Logger& logger, int q = 0); // Need the log for input, q is charge
  Molecule(const Molecule& other); // Copy constructor
  ~Molecule(); // Deletes the atom array
  // Accessors
  int getCharge() const { return charge; }
  int getNel() const { return nel; }
  int getMultiplicity() const { return multiplicity; }
  double getEnuc() const { return enuc; }
  Atom& getAtom(int i) const { return atoms[i]; } // Return atom i
  BF& getBF(int i) const { return bfset.getBF(i); } // Return basis func. i
  // Routines
  void rotate(const Matrix& U);
  void translate(double x, double y, double z);
  int nalpha() const;
  int nbeta() const;
  void calcEnuc(); 
  Vector com() const;
  Vector getInertia(bool shift = false);
  double dist(int i, int j) const;
  double dist2(int i, int j) const;
  double bondAngle(int i, int j, int k) const;
  double oopAngle(int i, int j, int k, int l) const;
  double torsionAngle(int i, int j, int k, int l) const;
  std::string rType() const;
  Vector rConsts(int units) const;
};

#endif
