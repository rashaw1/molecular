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
  Matrix last_dens;
  Matrix last_CP;
  double energy, last_energy;
public:
  // Constructor
  SCF(Molecule& m, Fock& f);
  // Routines
  void calcE();
  Vector calcErr() const;
  bool testConvergence() const;
  void rhf();
  void uhf();
};

#endif
