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
#include "mp2.hpp"
#include "gaussquad.hpp"
#include <functional>
#include <cmath>

double quadratic(double x) {
  double value = 1.5 + 2.3*x + 0.8*x*x;
  return value;
}

double exponential(double x) {
  return exp(x);
}

double logar(double x) {
  return log((x+3.0)/2.0);
}

double hard_test(double x) {
  double x2 = x * x;
  return x2 * x2 * exp(-2.1*x2);
}

double gaussiantest(double x) {
  return exp(-x*x);
}

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
	  case 4: { // MP2
		  log.title("MP2 CALCULATION");
		  MP2 mp2obj(focker);
		  mp2obj.transformIntegrals();
		  log.print("Integral transformation complete.\n");
		  log.localTime();
		  mp2obj.calculateEnergy();
		  log.result("MP2 Energy Correction = " + std::to_string(mp2obj.getEnergy()) + " Hartree");
		  log.result("Total Energy = " + std::to_string(hf.getEnergy() + mp2obj.getEnergy()) + " Hartree");
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

      // GCQuadrature test
      log.title("GC QUAD TEST");
      GCQuadrature gc;
      gc.initGrid(1024, TWOPOINT);
      log.localTime();
      std::function<double(double)> integrand = exponential;
      gc.integrate(integrand, 1e-12);
      log.print(std::to_string(gc.getI()));
      log.localTime();
      integrand = logar;
      gc.integrate(integrand, 1e-12);
      log.print(std::to_string(gc.getI()));
      log.localTime();
      integrand = hard_test;
      gc.integrate(integrand, 1e-12);
      log.print(std::to_string(gc.getI()));
      log.localTime();
      integrand = exponential;
      gc.transformRMinMax(16, 3.75);
      gc.integrate(integrand, 1e-12);
      log.print(std::to_string(gc.getI()));
      log.localTime();

      // Close file streams
      input.close();
      output.close();
      err.close();
    }
  }

  return 0;
}

