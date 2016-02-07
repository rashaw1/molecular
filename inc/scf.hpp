/*
 *
 *   PURPOSE: To declare a class SCF which is used for doing HF-SCF calculations
 *            of both the rhf and uhf kind.
 *
 *       class SCF
 *             owns: focker - a Fock class instance, for doing the bulk of the legwork
 *                   molecule - a reference to the molecule in question
 *             data: last_dens - the previous density matrix, for convergence checking
 *                   last_CP - the previous coefficient matrix, for DIIS
 *                   energy - the SCF energy
 *             routines: calcE - calculates the energy
 *                       rhf - does a restricted HF calculation
 *                       uhf - does an unrestricted HF calculation
 *
 *   DATE             AUTHOR               CHANGES
 *   =================================================================================
 *   15/09/15         Robert Shaw          Original code.
 *
 */

#ifndef SCFHEADERDEF
#define SCFHEADERDEF

#include "fock.hpp"
#include "matrix.hpp"
#include "molecule.hpp"

// Declare forward dependencies
class IntegralEngine;
class Vector;

// Begin class
class SCF
{
private:
  Molecule& molecule;
  Fock& focker;
  double energy, last_energy, one_E, two_E, error, last_err;
public:
  // Constructor
  SCF(Molecule& m, Fock& f);
  // Routines
  void calcE();
	double getEnergy() const { return energy; } 
  double calcE(const Matrix& hcore, const Matrix& dens, const Matrix& fock); 
  Vector calcErr(const Matrix& F, const Matrix& D, Matrix S);
  Vector calcErr();
  bool testConvergence();
  void rhf();
  void uhf();
};

#endif
