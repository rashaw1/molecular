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
#include "gshell.hpp"
#include "ecp.hpp"
#include <functional>
#include <cmath>
#include "ecpint.hpp"

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
	  
	  // Radial integral test
	  log.localTime();
	  double C[3] = {0.8, 0.8, 0.8};
	  ECP U1(C);
	  U1.addPrimitive(0, 0, 70.03441, 49.953835, false);
	  U1.addPrimitive(0, 0, 31.895971, 369.978238, false);
	  U1.addPrimitive(0, 0, 7.353728, 10.054544, false);
	  U1.addPrimitive(0, 1, 47.265183, 297.3055, false);
	  U1.addPrimitive(0, 1, 23.008133, 84.614009, false);
	  U1.addPrimitive(0, 1, 46.664268, 297.348786, false);
	  U1.addPrimitive(0, 1, 22.100979, 84.775187, false);
	  U1.addPrimitive(0, 2, 50.768165, -46.472192, false);
	  U1.addPrimitive(0, 2, 15.479497, -0.733528, false);
	  U1.addPrimitive(0, 2, 2.877154, 0.172203, false);
	  U1.addPrimitive(0, 2, 50.782228, -46.458997, false);
	  U1.addPrimitive(0, 2, 15.559383, -0.782332, false);
	  U1.addPrimitive(0, 2, 1.991196, 0.220450, false);
	  U1.addPrimitive(0, 3, 15.437567, -2.854575, false);
	  U1.addPrimitive(0, 3, 22.055742, -5.234658);
	  
	  double centerA[3] = { 1.5, 0.1, 0.9 };
	  double centerB[3] = { 0.8, 0.8, 0.8 };
	  //double centerA[3] = {0.0, 0.0, 0.0};
	  //double centerB[3] = {0.0, 0.0, 0.0};
	  GaussianShell shellA(centerA, 0);
	  shellA.addPrim(2808.6, 0.001606);
	  shellA.addPrim(421.18, 0.008393);
	  shellA.addPrim(50.3457, 0.069578);
	  shellA.addPrim(17.9133, -0.389908);
	  shellA.addPrim(3.80531, 0.694497);
	  shellA.addPrim(1.74968, 0.491354);
	  shellA.addPrim(0.448555, 0.022637);
	  shellA.addPrim(0.164498, -0.003723);
	  GaussianShell shellB(centerB, 0);
	  shellB.addPrim(2808.6, 0.001606);
	  shellB.addPrim(421.18, 0.008393);
	  shellB.addPrim(50.3457, 0.069578);
	  shellB.addPrim(17.9133, -0.389908);
	  shellB.addPrim(3.80531, 0.694497);
	  shellB.addPrim(1.74968, 0.491354);
	  shellB.addPrim(0.448555, 0.022637);
	  shellB.addPrim(0.164498, -0.003723);
	 
	  ECPIntegral ecpint;
	  Matrix values;
	  /*log.title("Type 1 test");
	  ecpint.type1(U1, shellA, shellB, centerA, centerB, values);
	  log.localTime();
	  log.print(values);
	  
	  log.title("Type 2 test");
	  int l = 0;
	  ThreeIndex values2(shellA.ncartesian(), shellB.ncartesian(), 2*l + 1);
	  log.localTime();
	  ecpint.type2(l, U1, shellA, shellB, centerA, centerB, values2);
	  log.localTime();
	  for (int na = 0; na < shellA.ncartesian(); na++) {
		  for (int nb = 0; nb < shellB.ncartesian(); nb++) {
			  for (int mu = -l; mu<=l; mu++) {
				  std::cout << na << " " << nb << " " << mu << " " << values2(na, nb, l + mu) <<"\n";
			  }
		  }
	  }*/
	  
	  log.title("ECP Test");
	  log.localTime();
	  ecpint.compute_shell_pair(U1, shellA, shellA, values);
	  log.localTime();
	  log.print(values);
	 
      // Close file streams
      input.close();
      output.close();
      err.close();
    }
  }

  return 0;
}

