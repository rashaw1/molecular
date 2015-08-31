/*
 *
 *   PURPOSE: To declare a selection of functions that are needed for
 *            the input/output interface of the molecular suite of programs.
 * 
 *           
 *            Routines: getAtomMass(q) - gets the mass (in amu) of an atom with
 *										with atomic number q.
 *                      getAtomName(q) - get the name, e.g. 'H' or 'Ca', of that atom.
 *                      getAtomCharge(n) - does the reverse; gives the atomic number
 *                                         of an atom with name n.
 *                      getShellName(l) - returns the name of a shell with angular mom.
 *                                        l - e.g. 's', or 'd'.
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

// General routines
#include <string>

// Get the atomic mass/name of an atom with atomic number q
double getAtomMass(int q);
std::string getAtomName(int q);

// Get the atomic number of an atom with name n
int getAtomCharge(const std::string& n);

// Get the text name for a shell of angular momentum l
std::string getShellName(int l);
  
#endif
