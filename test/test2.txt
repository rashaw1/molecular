#include <iostream>
#include <fstream>
#include <iomanip>
#include "mathutil.hpp"
#include "mvector.hpp"
#include "logger.hpp"
#include "tensor6.hpp"
#include "tensor4.hpp"
#include "molecule.hpp"
#include "error.hpp"
#include "integrals.hpp"
#include <stdexcept>
#include "fock.hpp"
#include "scf.hpp"

int main (int argc, char* argv[])
{
  std::ifstream input("inputfile.mol");
  if (input.is_open()) { std::cout << "input opened\n"; }
  std::ofstream output("out.out");
  if (output.is_open()) { std::cout << "output opened\n"; }
  Logger log(input, output, std::cout);
  std::cout << "log made\n";
  try {
  Molecule mol(log, 1);
  std::cout << "molecule made\n";

  log.print(mol, true);
  log.print(log.getBasis(), true);
  //log.print(mol.getBF(7, 6));
  //log.print(mol.getBF(7, 6).getPBF(0));
  log.result("This is a result.");
  /*  std::cout << mol.getAtom(0).getNShellPrims(1);
  std::cout << "\n";
  std::cout << mol.getAtom(0).getBF(0).getPBF(0).getNorm() << " ";
  std::cout << mol.getAtom(0).getBF(0).getPBF(1).getNorm() << " ";
  std::cout << mol.getAtom(0).getBF(0).getPBF(2).getNorm() << "\n";

  //Vector c = mol.getAtom(0).getCoords();  
  //Vector cd = mol.getAtom(1).getCoords();
  //Vector d = mol.getAtom(2).getCoords();
  */
  log.print("PRELIMINARIES FINISHED\n");
  log.localTime();
  IntegralEngine integral(mol);
  Fock focker(integral, mol);
  SCF hf(mol, focker);
  hf.rhf();
  //Vector tempi;
  //tempi = integral.twoe(mol.getAtom(0), mol.getAtom(0), mol.getAtom(0), mol.getAtom(0), 0, 0, 0, 1);
  //tempi.print(); std::cout << "\n\n";
  //Matrix twoprims;
  //twoprims = integral.twoe(mol.getAtom(0).getShellPrim(0, 0), mol.getAtom(0).getShellPrim(0, 1), mol.getAtom(1).getShellPrim(0, 0), mol.getAtom(1).getShellPrim(0, 1),
  //			   mol.getAtom(0).getCoords(), mol.getAtom(0).getCoords(), mol.getAtom(1).getCoords(), mol.getAtom(1).getCoords());
  //std::cout << "Prim test:\n"; twoprims.print(); std::cout <<"\n\n";
  /*//std::cout << "Nuc Attract = " << integral.nucAttract(mol.getAtom(0).getBF(0).getPBF(2), mol.getAtom(2).getBF(0).getPBF(2), 
  //						       c, d, cd) << "\n";
  std::cout << integral.getOverlap(6, 6) << "  " << integral.getKinetic(6, 6) << "\n";
    log.print("\n");
  //log.print(mol.getAtom(2).getBF(28));
  //log.print("\n");
  //log.print(mol.getAtom(2).getBF(24));
  //log.print("\n");
  */
  log.finalise();
  } catch (Error e) {
    log.error(e);
  }
  input.close();
  std::cout << "input closed\n";
  output.close();
  std::cout << "output closed\n";

  return 0;
}
