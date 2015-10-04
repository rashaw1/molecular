/*
 *
 *   PURPOSE: The main program loop for the molecular suite of ab initio quantum
 *            chemistry routines.
 * 
 *   DATE              AUTHOR                CHANGES
 *   ===========================================================================
 *   23/09/15          Robert Shaw           Original code.
 *
 */

#include <iostream>
#include <fstream>
#include <string>
#include "logger.hpp"
#include "molecule.hpp"
#include "error.hpp"
#include "integrals.hpp"
#include "fock.hpp"
#include "scf.hpp"

int main (int argc, char* argv[])
{
  if (argc == 1) { 
    std::cerr << "You must provide an input file.\nUSAGE: ./molecular inputfile.mol\n";
  } else {
    // Get the input file name, and create the error and output file names
    std::string ifname = argv[1];
    std::string ofname = ifname;
    std::size_t pos = ofname.find('.');
    if (pos != std::string::npos) { 
      ofname.erase(pos, ofname.length());
    }    
    std::string efname = ofname + ".log";
    ofname += ".out";

    // Open the file streams
    std::ifstream input(ifname);
    if (!input.is_open()) {
      std::cerr << "Could not open input file.\n";
    } else {
      std::ofstream output(ofname);
      std::ofstream err(efname);
      
      // Create the logger
      Logger log(input, output, err);
      log.init();

      try{
	// Create the molecule and print initial output
	Molecule mol(log, 1);
	log.print(mol, true);
	log.print(log.getBasis(), log.bprint());
	log.print("\nPRELIMINARIES FINISHED");
	log.localTime();
	output.flush();
	// Make an integral engine
	IntegralEngine integral(mol);

	// All calculations will need some form of HF
	Fock focker(integral, mol);
	SCF hf(mol, focker);

	// Run the commands given
	int cmd = log.nextCmd();
	while (cmd!=0) {
	  switch(cmd) {
	  case 1: { // HF
	    if(mol.getMultiplicity() > 1){
	      hf.uhf();
	    } else if ((mol.getNel()%2) != 0) {
	      hf.uhf();
	    } else {
	      hf.rhf();
	    }
	    break;
	  }
	  case 2: { // RHF
	    hf.rhf();
	    break;
	  }
	  case 3: { // UHF
	    hf.uhf();
	    break;
	  }
	  default: { }
	  }
	  cmd = log.nextCmd();
	}

	// Finalise the run
	log.finalise();
      } catch (Error e){
	log.error(e);
      }
      
      // Close file streams
      input.close();
      output.close();
      err.close();
    }
  }
  return 0;
}
